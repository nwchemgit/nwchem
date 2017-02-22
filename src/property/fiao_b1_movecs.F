      subroutine fiao_b1_movecs(rtdb,basis,geom, ncomp, ! IN
     &                          g_vecB1)         ! OUT

c $Id$

c     This routine is a modification of fiao_b1_movecs.F. We calculate
c     the MO vectors perturbed by a STATIC external magnetic field but
c     not with GIAOs. This is a heavily stripped-down version of
c     hnd_giaox.

c     g_vecB1 should have at least dimensions (nbf, nocc, 3)

c     Note: integrals have to be initialized by the calling routine.
c     Note: the CPHF call below terminates the integrals
c
c  Written by J. Autschbach, SUNY Buffalo
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

      integer rtdb      ! [in] rtdb handle
      integer basis     ! [in] basis handle
      integer geom      ! [in] geometry handle
      integer ncomp     ! [in] number of spin components
      integer g_vecB1(ncomp)! [out] B-field perturbed MO coeffs FIAO
      integer nclosed(2), nopen(2), nvirt(2), ndens, nbf, nmo
      integer sh_atom, ixy, ix, iy, iatom, iocc, ifld, ioff
      integer alo(3), ahi(3), blo(3), bhi(3), clo(3), chi(3)
      integer dlo(3), dhi(3)
      integer l_occ, k_occ, l_eval, k_eval
      integer g_dens(3), g_s10, g_d1, g_rhs, g_fock, g_u(2)
      integer vectors(2), geomnew, i, j, ij, g_xc(3)
      integer vectors1(2)
      integer ndata,m1,m2
      double precision atn, tol2e, val
      character*3 scftyp
      character*16 tag
      character*32 element
      character*256 cphf_rhs, cphf_sol
      character*2 symbol
      double precision origin(3)
      data origin/0d0,0d0,0d0/
      integer g_rhs_im,read_grhs_fiaob1
      character*255 aorespfilename
      logical dft_CPHF1_read,dft_CPHF1_write
      character*(*) lbl_cphfaoresp
      parameter(lbl_cphfaoresp='aoresp_fiao_b1')
      integer nat
      parameter (nat=1)
c ------- Added for unrestricted calc ----- START
      integer ndir,    ! nr directions (x,y,z)
     &        ntot,    ! sum_{i=1,npol} nocc(i)*nvirt(i)
     &        ispin,disp,shift,
     &        nocc(2), ! store nr occupations 
     &        npol     ! nr of polarizations =1 (  restricted) 
                       !                     =2 (unrestricted) calc
c ------- Added for unrestricted calc ----- END
      logical  cphf2, file_write_ga, file_read_ga, cphf
      external cphf2, file_write_ga, file_read_ga, cphf
      external get_nocc,update_rhs_angmomAO,
     &         get_vecB1_opt1

      logical  status, debug
      double precision ppm
      data tol2e   /1.0d-10/
      ndir=3 ! nr directions (x,y,z)
c     ==================================================================
c
      if (ga_nodeid().eq.0) write(luout,9999)
      debug = .false. .and. (ga_nodeid().eq.0) ! special debugging
c     there is a possibility that this routine is called from tddft
c     in which case use_theory is not set. We set it to 'dft' in that case,
c     assuming that we are indeed calling from some DFT response code

      if ((use_theory.ne.'scf').and.(use_theory.ne.'dft')) then
        if (ga_nodeid().eq.0) write (luout,*)
     &     'fiao_b1_movecs: assuming DFT/TDDFT'
        use_theory = 'dft'
      end if
c
c     Current CPHF does not handle symmetry 
c     Making C1 geometry and store it on rtdb
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
     &    errquit('fiao_b1: could not get nbf',0, BASIS_ERR)
      if (.not. ma_push_get(mt_dbl,2*nbf,'occ num',l_occ,k_occ)) call
     &    errquit('fiao_b1: ma_push_get failed k_occ',0,MA_ERR)
      if (.not. ma_push_get(mt_dbl,2*nbf,'eigenval',l_eval,k_eval)) call
     &    errquit('fiao_b1: ma_push_get failed k_eval',0,MA_ERR)

      call hnd_prp_vec_read(rtdb,geom,basis,     ! in : handles
     &                      nbf,                 ! out: nr basis functions
     &                      nclosed,nopen,nvirt, ! out: occupation numbers
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

      if (ga_nodeid().eq.0) then
        write(*,10) nocc(1)   ,nocc(2),
     &              nopen(1)  ,nopen(2),
     &              nclosed(1),nclosed(2),
     &              nvirt(1)  ,nvirt(2),scftyp,ntot
 10    format('nocc =(',i3,',',i3,') ',
     &        'nopen=(',i3,',',i3,') ',
     &        'nclos=(',i3,',',i3,') ',
     &        'nvirt=(',i3,',',i3,') ',
     &        'scftyp=',a,' ntot=',i3)
      endif
      if (.not.ma_pop_stack(l_eval)) call
     &    errquit('fiao_b1: ma_pop_stack failed k_eval',0,MA_ERR)
      if (.not.ma_pop_stack(l_occ)) call
     &    errquit('fiao_b1: ma_pop_stack failed k_occ',0,MA_ERR)

      if (npol.lt.1 .or. npol.gt.2) call errquit(
     &   'giao_b1_mov: npol out of range',0,CALC_ERR)
      if (npol.ne.ncomp) call errquit(
     &   'giao_b1_mov: npol not in sysnc with calling routine',
     &   0,CALC_ERR)

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
     &   call errquit('fiao_b1: ga_create failed g_rhs',0,GA_ERR)
      call ga_zero(g_rhs)  

      call update_rhs_angmomAO(
     &               g_rhs,  !in/out:
     &               vectors,!in : MO vectors
     &               nocc,   !in : nr. occup   MOs
     &               nvirt,  !in : nr. virtual MOs   
     &               npol,   !in : nr. polarizations
     &               nbf,    !in : nr. basis functions
     &               nmo,    !in : nr. MOs
     &               basis,  !in : basis handle
     &               debug)  !in : logical var for debugging

      call schwarz_tidy()
      call int_terminate()
      if (debug) then
       if (ga_nodeid().eq.0)
     &  write(*,*) '---- g_rhs-AFT-angmomAO-------- START'
        call ga_print(g_rhs)
       if (ga_nodeid().eq.0)
     &  write(*,*) '---- g_rhs-AFT-angmomAO--------  END'
      endif ! end-if-debug
c     there is a division by -4 somewhere in the
c     CPKS solver so we scale the rhs by -4 in order to
c     get the correct magnitude of the result back
      blo(1) = 1
      bhi(1) = ntot
      blo(2) = 1
      bhi(2) = ndir  
      if      (npol.eq.1) then
        call nga_scale_patch(g_rhs,blo,bhi,-4.0d0)
      else if (npol.eq.2) then
        call nga_scale_patch(g_rhs,blo,bhi,-2.0d0)
      endif    
c 0000000000000 read (g_rhs,g_rhs_im) from file 000000 START
       call util_file_name(lbl_cphfaoresp,
     &                     .false.,.false.,aorespfilename)
      read_grhs_fiaob1=0 

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
         read_grhs_fiaob1=1
       else
        if (.not. rtdb_put(rtdb,'cphf2-aores:guess', 
     &          mt_log, 1,.true.)) call
     $     errquit('giao_b1: failed to write skew ', 0, RTDB_ERR)
      endif
c 0000000000000 read (g_rhs,g_rhs_im) from file 000000 END
c       if (ga_nodeid().eq.0)
c     &  write(*,*) 'COMPUTE cphf fiao_b1 data ...'    
c     Write ga_rhs to disk 
       call cphf_fname('cphf_rhs',cphf_rhs)
       call cphf_fname('cphf_sol',cphf_sol)
       if(.not.file_write_ga(cphf_rhs,g_rhs)) call errquit
     $  ('fiao_b1: could not write cphf_rhs',0, DISK_ERR)

c     Call the CPHF routine 
c     We do need to tell the CPHF that the density is skew symmetric.
c     Done via rtdb, put cphf:skew .false. on rtdb and later remove it.
       if (.not. rtdb_put(rtdb, 'cphf:skew', mt_log, 1,.false.)) call
     $   errquit('fiao_b1: failed to write skew ', 0, RTDB_ERR)
       if (.not.cphf2(rtdb)) call errquit
     $  ('fiao_b1: failure in cphf ',0, RTDB_ERR)
       if (.not. rtdb_delete(rtdb, 'cphf:skew')) call
     $   errquit('fiao_b1: rtdb_delete failed ', 0, RTDB_ERR)
c
c     Occ-virt blocks are the solution pieces of the CPHF
c     Read solution vector from disk and put solutions in U matrices
       call ga_zero(g_rhs)
       if(.not.file_read_ga(cphf_sol,g_rhs)) call errquit
     $  ('fiao_b1: could not read cphf_rhs',0, DISK_ERR) 

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
c     & write(*,*) 'FA-BEF get_vecB1-fb1'

      call get_vecB1_opt1(
     &               g_vecB1,    ! out:
     &               g_rhs,      ! in : g_rhs vector
     &               vectors,    ! in : MO vectors
     &               nbf,        ! in : nr. basis functions
     &               nmo,        ! in : nr. MOs
     &               npol,       ! in : nr polarizations
     &               nocc,       ! in : nr. occupied MOs
     &               nvirt,      ! in : nr. virtual  MOs
     &               debug)      ! in : = .true. allow debugging

c      if (ga_nodeid().eq.0)
c     & write(*,*) 'FA-AFT get_vecB1-fb1'
      if (debug) then
       if (ga_nodeid().eq.0)
     &  write(*,*) '------- g_vecB1-fb1-nw ---- START'
        do ispin=1,npol
         call ga_print(g_vecB1(ispin))
        enddo
       if (ga_nodeid().eq.0)
     &  write(*,*) '------- g_vecB1-fb1-nw ---- END'
      endif ! end-if-debug
      if (.not.ga_destroy(g_rhs)) call 
     &    errquit('fiao_b1: ga_destroy failed g_rhs',0,GA_ERR)
      do ispin=1,npol
       if (.not.ga_destroy(vectors(ispin))) call 
     &   errquit('giao_b1: ga_destroy failed vectors',0,GA_ERR)
      enddo ! end-loop-ispin
      call ga_sync()
      return

 7000 format(/,10x,'B-field perturbed MOs cannot be calculated for',
     1      ' UHF or ROHF wave functions at present')
 9983 format(6x,i1,3x,3f12.4)
 9985 format(10x,3f12.4,/)
 9999 format(
     1 /,10x,54(1h-),/,
     2 10x,'Calculating magnetic field perturbed MO vectors (FIAO)',/,
     3 10x,54(1h-),/)

      end
