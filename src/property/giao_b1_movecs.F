      subroutine giao_b1_movecs(rtdb,basis,geom, ncomp,  ! IN
     &                          g_vecB1, g_sket1) ! OUT
c $Id$

c     This routine is a modification of hnd_giaox.F. Instead of
c     computing NMR shielding tensors we simply calculate the MO vectors
c     perturbed by a STATIC external magnetic field. I left most of the
c     comments from hnd_giaox and the variable names unmodified.

c     g_vecB1 should have at least dimensions (nbf, nocc, 3)
c     g_sket1 should have dimensions (nbf, nbf, 3)

c     Note: integrals have to be initialized by the calling routine.
c     Note: the CPHF call below terminates the integrals
c
c  Written by J. Autschbach, SUNY Buffalo
c
c  Extension to spin-unrestricted case 
c          by F. Aquino,     Northwestern University 
c          03-15-12

      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "stdio.fh"
#include "apiP.fh"
#include "prop.fh"
#include "bgj.fh"
#include "case.fh"
      integer rtdb       ! [in] rtdb handle
      integer basis      ! [in] basis handle
      integer geom       ! [in] geometry handle
      integer ncomp      ! [in] number of spin components
      integer g_vecB1(ncomp) ! [out] B-field perturbed MO coeffs GIAO
      integer g_sket1 ! [out] GIAO right hand side overlap derivative
      integer nclosed(2), nopen(2), nvirt(2), ndens, nbf, nmo
      integer ixy, ix, iy, iatom, iocc, ifld, ioff
      integer alo(3), ahi(3), blo(3), bhi(3), clo(3), chi(3)
      integer dlo(3), dhi(3)
      integer l_occ, k_occ, l_eval, k_eval
      integer g_dens(3), g_s10, g_d1, g_rhs, g_fock, g_u(2)
      integer vectors(2), geomnew, i, j, ij, g_xc(3)
      integer vectors1(2)
      integer ndata,m1,m2
      double precision atn, tol2e, val
      double precision jfac(12),kfac(12),a(6),xfac
      character*3 scftyp
      character*16 tag
      character*32 element
      character*256 cphf_rhs, cphf_sol
      character*2 symbol
c ===== added for unrestricted calc ===== START
      integer ndir,    ! nr directions (x,y,z)
     &        ntot,    ! sum_{i=1,npol} nocc(i)*nvirt(i)
     &        ispin,disp,shift,nind_jk,
     &        nocc(2), ! store nr occupations 
     &        npol     ! nr of polarizations =1 (  restricted) 
                       !                     =2 (unrestricted) calc    
      external get_d1_giao_b1,update_rhs_shfock,
     &         get_fock2e,update_rhs_threeAOints,get_vecB1,
     &         get_nocc,update_rhs_eS10,get_vecB1_opt2
c ===== added for unrestricted calc ===== END
      integer nbq, nextbq, ncosbq
      integer g_rhs_im,read_grhs_giaob1,n_data
      character*255 aorespfilename
      logical dft_CPHF1_read,dft_CPHF1_write
      character*(*) lbl_cphfaoresp
      parameter(lbl_cphfaoresp='aoresp_giao_b1')
      logical  cphf2, file_write_ga, file_read_ga, cphf
      external cphf2, file_write_ga, file_read_ga, cphf

      logical  oskel, status, debug
      logical  xc_gotxc
      external xc_gotxc
      double precision ppm
      data tol2e   /1.0d-10/
c     ==================================================================
      ndir=3 ! nr directions (x,y,z)
      if (ga_nodeid().eq.0) write(luout,9999)
      debug = .false. .and. (ga_nodeid().eq.0) ! special debugging
      
c      debug=.true.

      if (debug) then
        write (luout,*) 'giao_b1_movecs: xc_gotxc =',xc_gotxc()
        write (luout,*) 'giao_b1_movecs: use_theory =',use_theory
      end if

c     there is a possibility that this routine is called from tddft
c     in which case use_theory is not set. We set it to 'dft' in that case,
c     assuming that we are indeed calling from some DFT response code

      if ((use_theory.ne.'scf').and.(use_theory.ne.'dft')) then
        if (ga_nodeid().eq.0) write (luout,*)
     &     'giao_b1_movecs: assuming DFT/TDDFT'
        use_theory = 'dft'
      end if
c
c     Current CPHF does not handle symmetry 
c     Making C1 geometry and store it on rtdb
c     
      oskel = .false.
c
c     If DFT get part of the exact exchange defined
      xfac = 1.0d0
      if (use_theory.eq.'dft') xfac = bgj_kfac()
      nind_jk=12
      do ifld = 1,nind_jk
        jfac(ifld) =  0.0d0       ! used in update_rhs_shfock()
        kfac(ifld) = -1.0d0*xfac  ! used in update_rhs_shfock()
c        if (ga_nodeid().eq.0) then
c         write(*,144) ifld,jfac(ifld),kfac(ifld)
c  144    format('(j,k)(',i3,')=(',f15.8,',',f15.8,')')
c        endif
      enddo
c
c     Integral initialization
      call int_init(rtdb,1,basis)
      call schwarz_init(geom,basis)
      call hnd_giao_init(basis,1)
      call scf_get_fock_param(rtdb,tol2e)
      status = rtdb_parallel(.true.)
      
c     Get Unperturbed MO vectors and eigenvalues
c     First allocate some memory for occupation numbers and eigenvalues

      if (.not. bas_numbf(basis,nbf)) call
     &    errquit('giao_b1: could not get nbf',0, BASIS_ERR)
      if (.not. ma_push_get(mt_dbl,2*nbf,'occ num',l_occ,k_occ)) call
     &    errquit('giao_b1: ma_push_get failed k_occ',0,MA_ERR)
      if (.not. ma_push_get(mt_dbl,2*nbf,'eigenval',l_eval,k_eval)) call
     &    errquit('giao_b1: ma_push_get failed k_eval',0,MA_ERR)

      call hnd_prp_vec_read(rtdb,geom,basis,     ! in : handles
     &                      nbf,                 ! out: nr basis fcts.
     &                      nclosed,nopen,nvirt, ! out: occ/unocc
     &                      scftyp,              ! out: type calc
     &                      vectors,             ! out: MO vectors
     &                      dbl_mb(k_occ),       ! out: occupations
     &                      dbl_mb(k_eval),      ! out: DFT energies
     &                      nmo)                 ! out: nr MOs

      call get_nocc(rtdb,   ! in : rtdb handle
     &              nocc,   ! out: nr occupations
     &              npol,   ! out: nr of polarization
     &              nclosed,! in : nr closed shells
     &              nopen,  ! in : nr open shells
     &              nvirt,  ! in : nr virtual MOs
     &              scftyp, ! in : string = UHF or RHF
     &              ntot)   ! out: sum_{i,npol} nocc(i)*nvirt(i)

      if (debug) then
        write(*,10) nocc(1)   ,nocc(2),
     &              nopen(1)  ,nopen(2),
     &              nclosed(1),nclosed(2),
     &              nvirt(1)  ,nvirt(2),scftyp,ntot
 10    format('giao_b1_mov: nocc =(',i3,',',i3,') ',
     &        'nopen=(',i3,',',i3,') ',
     &        'nclos=(',i3,',',i3,') ',
     &        'nvirt=(',i3,',',i3,') ',
     &        'scftyp=',a,' ntot=',i3)
      endif

      if (npol.lt.1 .or. npol.gt.2) call errquit(
     &   'giao_b1_mov: npol out of range',0,CALC_ERR)
      if (npol.ne.ncomp) call errquit(
     &   'giao_b1_mov: npol not in sysnc with calling routine',
     &   0,CALC_ERR)

c
c     Get Unperturbed Density Matrix

      call hnd_prp_get_dens(rtdb,geom,basis,     ! in : handles
     &                      g_dens,ndens,        ! out: electron density
     &                      scftyp,              ! in : type calc
     &                      nclosed,nopen,nvirt) ! in : occupations numbers

c      if (ga_nodeid().eq.0) then
c        write(*,21) npol,nocc(1)   ,nocc(2),
c     &              nopen(1)  ,nopen(2),
c     &              nclosed(1),nclosed(2),
c     &              nvirt(1)  ,nvirt(2),scftyp,ntot
c 21    format('npol=',i3,' nocc =(',i3,',',i3,') ',
c     &        'nopen=(',i3,',',i3,') ',
c     &        'nclos=(',i3,',',i3,') ',
c     &        'nvirt=(',i3,',',i3,') ',
c     &        'scftyp=',a,' ntot=',i3)
c      endif

      if (debug) write (luout,*) 'unpertubed MOs and Pmat assembled'
c
c     Create U matrix of dimension (nbf,nmo,3) and zero
c     Use ahi for dimension and ahi array for chunking/blocking
      do ispin=1,npol
       alo(1) = nbf
       alo(2) = -1
       alo(3) = -1
       ahi(1) = nbf
       ahi(2) = nocc(ispin)
       ahi(3) = ndir
       if (.not.nga_create(MT_DBL,3,ahi,'U matrix',alo,g_u(ispin))) call 
     &    errquit('giao_b1: nga_create failed g_u',0,GA_ERR)
       call ga_zero(g_u(ispin))
      enddo ! end-loop-ispin
c 
c     Construction of right-hand side CPHF
c     Create CPHF array of proper dimension : (nocc*nvirt,3)
      ndata=2 ! 1st subspace corresponds to g_b, 
c             ! 2nd subspace corresponds to sol (if exists)
        if (.not. rtdb_put(rtdb,'cphf2-aores:ndata', 
     &          mt_int, 1,ndata)) call
     $     errquit('fiao_b1: failed to write skew ', 0, RTDB_ERR)
      if(.not.ga_create(MT_DBL,ntot,ndata*ndir,
     &   'RHS',-1,-1,g_rhs))
     &   call errquit('giao_b1: ga_create failed g_rhs',0,GA_ERR)
      call ga_zero(g_rhs)

c      if (ga_nodeid().eq.0)
c     &  write(*,*) 'FA-BEF update_rhs_eS10'

c      if (ga_nodeid().eq.0) then
c        write(*,70) nocc(1)   ,nocc(2),
c     &              nvirt(1)  ,nvirt(2),npol,nbf,nmo
c 70    format('BEF-update_rhs_eS10 nocc =(',i5,',',i5,') ',
c     &        'nvirt=(',i5,',',i5,') ',
c     &        '(npol,nbf,nmo)=(',i3,',',i3,',',i3,')')
c      endif

      call update_rhs_eS10(
     &            g_rhs,         !in/out:
     &            g_u,           !out:
     &            g_sket1,       !out:
     &            dbl_mb(k_eval),!in : energy values
     &            vectors,       !in : MO vectors
     &            nocc,          !in : nr.   occupied MOs
     &            nvirt,         !in : nr. unoccupied MOs
     &            npol,          !in : nr. polarizations
     &            nbf,           !in : nr. basis functions
     &            nmo,           !in : nr. MOs
     &            basis,         !in : basis handle
     &            debug)         !in : logical var for debugging

c
c     We also need the occupied-occupied contribution of g_u
c     contributing to the first order density matrix. As this block does
c     not change during the CPHF we can calculate it once and subtract
c     it from the RHS. We will reuse g_s10 as scratch space.

      call get_d1_giao_b1(g_d1,   ! out:
     &                    g_u,    ! in :
     &                    vectors,! in : MO vectors
     &                    nocc,   ! in : nr. occ shells
     &                    npol,   ! in : nr. polarizations
     &                    nbf,    ! in : nr. basis functions
     &                    nmo,    ! in : nr. MOs
     &                    debug)  ! in : =.true. for debugging


      call update_rhs_shfock(
     &                    g_rhs,  ! in/out: RHS used for cphf2/3
     &                    g_d1,   ! in    :
     &                    vectors,! in    : MO vectors
     &                    rtdb,   ! in    : rtdb  handle
     &                    geom,   ! in    : geom  handle
     &                    basis,  ! in    : basis handle 
     &                    jfac,   ! in    : exch factors
     &                    kfac,   ! in    : exch factors
     &                    tol2e,  ! in    : tolerance coeff
     &                    nocc,   ! in    : nr. occ  shells
     &                    nvirt,  ! in    : nr. virt shells
     &                    npol,   ! in    : nr. polarizations
     &                    nbf,    ! in    : nr. basis functions
     &                    nmo,    ! in    : nr. MOs
     &                    debug)  ! in    : =.true. for debugging

c     add the one-electron integrals for the GIAO
c     Fock operator derivative to the RHS:
      
      call update_rhs_threeAOints(
     &                    g_rhs,  ! in/out: RHS used for cphf2/3
     &                    vectors,! in    : MO vectors
     &                    rtdb,   ! in    : rtdb  handle 
     &                    basis,  ! in    : basis handle 
     &                    nocc,   ! in    : nr occ  shells
     &                    nvirt,  ! in    : nr virt shells
     &                    npol,   ! in    : nr. polarizations
     &                    nbf,    ! in    : nr. basis functions
     &                    nmo,    ! in    : nr. MOs
     &                    debug)  ! in    : logical for debugging

c     add the two-electron part of the GIAO Fock matrix
c     derivative to the RHS:

      call update_rhs_fock2e(
     &                    g_rhs,  ! in/out: 
     &                    vectors,! in : MO vectors
     &                    rtdb,   ! in : rtdb  handle
     &                    basis,  ! in : basis handle
     &                    geom,   ! in : geom  handle
     &                    g_dens, ! in : e-density
     &                    nocc,   ! in : nr. occ  shells
     &                    nvirt,  ! in : nr. virt shells
     &                    npol,   ! in : nr. polarizations
     &                    nbf,    ! in : nr. basis functions
     &                    nmo,    ! in : nr. MOs   
     &                    xfac,   ! in : exchange factor
     &                    tol2e,  ! in : tolerance coeff.
     &                    debug)  ! in : logical for debugging

      call schwarz_tidy()
      call int_terminate()

      blo(1) = 1
      bhi(1) = ntot
      blo(2) = 1
      bhi(2) = ndir  
      if      (npol.eq.1) then
        call nga_scale_patch(g_rhs,blo,bhi,-4.0d0)
      else if (npol.eq.2) then
        call nga_scale_patch(g_rhs,blo,bhi,-2.0d0)
      endif  

      if (debug) then
       if (ga_nodeid().eq.0)
     &  write(*,*) 'FA-AFT update_rhs_fock2e'
       if (ga_nodeid().eq.0)
     &  write(*,*) '---- g_rhs-AFT-fock2e-------- START'
        call ga_print(g_rhs)
       if (ga_nodeid().eq.0)
     &  write(*,*) '---- g_rhs-AFT-fock2e--------  END'
      endif ! end-if-debug

       call util_file_name(lbl_cphfaoresp,
     &                     .false.,.false.,aorespfilename)
      read_grhs_giaob1=0 
      if (.not. dft_CPHF1_read( ! file exists and read g_rhs guess
     &           aorespfilename,! in: filename
     &           npol,          ! in: nr polarization
     &           nocc,          ! in: nr occupied MOs
     &           nvirt,         ! in: nr virtual  MOs
     &           1,             ! in: nr. components
     &           g_rhs,         ! in: (ntot,3)       GA matrix
     &           g_rhs_im,      ! in: dummy
     &           .false.))      ! in: =T if (RE,IM) =F if RE
     & then
         read_grhs_giaob1=1
       else
        if (.not. rtdb_put(rtdb,'cphf2-aores:guess', 
     &          mt_log, 1,.true.)) call
     $     errquit('giao_b1: failed to write skew ', 0, RTDB_ERR)
      endif

      if (debug) then
       if (ga_nodeid().eq.0)
     &  write(*,*) '---- g_rhs-AFT-readfile-------- START'
        call ga_print(g_rhs)
       if (ga_nodeid().eq.0)
     &  write(*,*) '---- g_rhs-AFT-readfile--------  END'
      endif ! end-if-debug
  
       if (ga_nodeid().eq.0)
     &  write(*,*) 'COMPUTE cphf giao_b1 data ...'   
c
c     Write ga_rhs to disk 
      call cphf_fname('cphf_rhs',cphf_rhs)
      call cphf_fname('cphf_sol',cphf_sol)
      if(.not.file_write_ga(cphf_rhs,g_rhs)) call errquit
     $  ('giao_b1: could not write cphf_rhs',0, DISK_ERR)
c
c     Call the CPHF routine
c     
c     We do need to tell the CPHF that the density is skew symmetric.
c     Done via rtdb, put cphf:skew .false. on rtdb and later remove it.
      if (.not. rtdb_put(rtdb, 'cphf:skew', mt_log, 1,.false.)) call
     $   errquit('giao_b1: failed to write skew ', 0, RTDB_ERR)
      if (debug) write (luout,*) 'calling cphf'
      if (.not.cphf2(rtdb)) call errquit
     $  ('giao_b1: failure in cphf ',0, RTDB_ERR)
      if (.not. rtdb_delete(rtdb, 'cphf:skew')) call
     $   errquit('giao_b1: rtdb_delete failed ', 0, RTDB_ERR)
      if (debug) write (luout,*) 'cphf done'
c
c     Occ-virt blocks are the solution pieces of the CPHF
c     Read solution vector from disk and put solutions in U matrices
      call ga_zero(g_rhs)
      if(.not.file_read_ga(cphf_sol,g_rhs)) call errquit
     $  ('giao_b1: could not read cphf_rhs',0, DISK_ERR)  

      if (debug) then
c       if (ga_nodeid().eq.0)
c     &  write(*,*) '---- g_rhs-BEF-write2file------- START'
        call ga_print(g_rhs)
c       if (ga_nodeid().eq.0)
c     &  write(*,*) '---- g_rhs-BEF-write2file--------  END'
      endif ! end-if-debug

       call util_file_name(lbl_cphfaoresp,
     &                     .false.,.false.,aorespfilename)

       status=dft_CPHF1_write(
     &           aorespfilename,! in: filename
     &           npol,          ! in: nr polarization
     &           nocc,          ! in: nr occupied MOs
     &           nvirt,         ! in: nr virtual  MOs
     &           1,             ! in: nr. components
     &           g_rhs,         ! in: (ntot,3)       GA matrix
     &           g_rhs_im,      ! in: dummy
     &           .false.)       ! in: =T if (RE,IM) =F if RE

c 000000000000 move 2nd subspace to 1st 00000 START
          shift=ndir
          m1=shift+1
          m2=shift+ndir
          call ga_copy_patch('n',g_rhs,1,ntot,m1,m2, 
     $                           g_rhs,1,ntot,1 ,ndir)     
c 000000000000 move 2nd subspace to 1st 00000 END
c      if (ga_nodeid().eq.0)
c     & write(*,*) 'FA-BEF get_vecB1_opt2'

      call get_vecB1_opt2(
     &                  g_vecB1,    ! out:
     &                  g_rhs,      ! in : g_rhs vector (occ-virt of g_u)
     &                  g_u,        ! in : occ-occ of g_u
     &                  vectors,    ! in : MO vectors
     &                  nbf,        ! in : nr. basis functions
     &                  nmo,        ! in : nr. MOs
     &                  npol,       ! in : nr polarizations
     &                  nocc,       ! in : nr. occupied MOs
     &                  nvirt,      ! in : nr. virtual  MOs
     &                  debug)      ! in : = .true. allow debugging

c      if (ga_nodeid().eq.0)
c     & write(*,*) 'FA-AFT get_vecB1_opt2'

      if (debug) then
       if (ga_nodeid().eq.0)
     &  write(*,*) '------- g_vecB1-gb1-nw ---- START'
       do ispin=1,npol
        call ga_print(g_vecB1(ispin))
       enddo
       if (ga_nodeid().eq.0)
     &  write(*,*) '------- g_vecB1-gb1-nw ---- END'
      endif ! end-if-debug

      do ispin=1,npol
       if (.not.ga_destroy(g_u(ispin))) call 
     &   errquit('giao_b1: ga_destroy failed vectors',0,GA_ERR)
       if (.not.ga_destroy(vectors(ispin))) call 
     &   errquit('giao_b1: ga_destroy failed vectors',0,GA_ERR)
      enddo
      do ispin=1,ndens
       if (.not.ga_destroy(g_dens(ispin))) call 
     &    errquit('giao_b1: ga_destroy failed g_dens',0,GA_ERR)
      enddo ! end-loop-ispin
c      RHS arrays are no longer needed
        if (.not.ga_destroy(g_rhs)) call 
     &     errquit('fiao_f1: ga_destroy failed g_rhs',0,GA_ERR)
      call ga_sync()

      if (.not.ma_pop_stack(l_eval)) call
     &    errquit('giao_b1: ma_pop_stack failed k_eval',0,MA_ERR)
      if (.not.ma_pop_stack(l_occ)) call
     &    errquit('giao_b1: ma_pop_stack failed k_occ',0,MA_ERR)

      return

 7000 format(/,10x,'B-field perturbed MOs cannot be calculated for',
     1      ' UHF or ROHF wave functions')
 9983 format(6x,i1,3x,3f12.4)
 9985 format(10x,3f12.4,/)
 9999 format(
     1 /,10x,54(1h-),/,
     2 10x,'Calculating magnetic field perturbed MO vectors (GIAO)',/,
     3 10x,54(1h-),/)
      end
