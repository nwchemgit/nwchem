      subroutine giao_b1_movecs(rtdb,basis,geom, g_vecB1, g_sket1)

c $Id$

c     This routine is a modification of hnd_giaox.F. Instead of
c     computing NMR shielding tensors we simply calculate the MO vectors
c     perturbed by a STATIC external magnetic field. I left most of the
c     comments from hnd_giaox and the variable names unmodified.

c     g_vecB1 should have at least dimensions (nbf, nocc, 3)
c     g_sket1 should have dimensions (nbf, nbf, 3)

c     Note: integrals have to be initialized by the calling routine.
c     Note: the CPHF call below terminates the integrals

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
c
      integer rtdb    ! [in] rtdb handle
      integer basis   ! [in] basis handle
      integer geom    ! [in] geometry handle
      integer g_vecB1 ! [out] B-field perturbed MO coefficients GIAO
      integer g_sket1 ! [out] GIAO right hand side overlap derivative

      integer nclosed(2), nopen(2), nvirt(2), ndens, nbf, nmo
      integer ixy, ix, iy, iatom, iocc, ifld, ioff
      integer alo(3), ahi(3), blo(3), bhi(3), clo(3), chi(3)
      integer dlo(3), dhi(3)
      integer l_occ, k_occ, l_eval, k_eval

      integer g_dens(3), g_s10, g_d1, g_rhs, g_fock, g_u
      integer vectors(2), geomnew, i, j, ij, g_xc(3)
      double precision atn, tol2e, val
      double precision jfac(3),kfac(3),a(6),xfac
      character*3 scftyp
      character*16 tag
      character*32 element
      character*256 cphf_rhs, cphf_sol
      character*2 symbol

      double precision origin(3)
      data origin/0d0,0d0,0d0/

      integer nat
      parameter (nat=1)

      integer nbq, nextbq, ncosbq

      logical  cphf2, file_write_ga, file_read_ga, cphf
      external cphf2, file_write_ga, file_read_ga, cphf

      logical     oskel, status, debug
      logical  xc_gotxc
      external xc_gotxc
      double precision ppm
      data tol2e   /1.0d-10/

c     ==================================================================
c
      if (ga_nodeid().eq.0) write(luout,9999)

      call ga_zero(g_vecB1) ! initialize array for main result

      debug = .false. .and. (ga_nodeid().eq.0) ! special debugging

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
c
      xfac = 1.0d0
      if (use_theory.eq.'dft') xfac = bgj_kfac()
c
c     Integral initialization
c
c      call int_init(rtdb,1,basis)
c      call schwarz_init(geom,basis)
      call hnd_giao_init(basis,1)
      call scf_get_fock_param(rtdb,tol2e)
c
      status = rtdb_parallel(.true.)

c     Get Unperturbed MO vectors and eigenvalues
c     First allocate some memory for occupation numbers and eigenvalues
c
      if (.not. bas_numbf(basis,nbf)) call
     &    errquit('giao_b1: could not get nbf',0, BASIS_ERR)
      if (.not. ma_push_get(mt_dbl,2*nbf,'occ num',l_occ,k_occ)) call
     &    errquit('giao_b1: ma_push_get failed k_occ',0,MA_ERR)
      if (.not. ma_push_get(mt_dbl,2*nbf,'eigenval',l_eval,k_eval)) call
     &    errquit('giao_b1: ma_push_get failed k_eval',0,MA_ERR)
      call hnd_prp_vec_read(rtdb,geom,basis,nbf,nclosed,nopen,
     &                      nvirt,scftyp,vectors,dbl_mb(k_occ),
     &                      dbl_mb(k_eval),nmo)
c
c     Get Unperturbed Density Matrix
c
      call hnd_prp_get_dens(rtdb,geom,basis,g_dens,ndens,scftyp,
     &                      nclosed,nopen,nvirt)

      if (debug) write (luout,*) 'unpertubed MOs and Pmat assembled'
c
c     Error exit if scftyp equals UHF (= ROHF)
c
      if (scftyp.eq.'UHF') then
          if (ga_nodeid().eq.0) write(luout,7000)
          call errquit('giao_b1: incompatible SCF type for Response',
     &       0,INPUT_ERR)
      endif
c
c     Create U matrix of dimension (nbf,nmo,3) and zero
c     Use ahi for dimension and ahi array for chunking/blocking
c
      alo(1) = nbf
      alo(2) = -1
      alo(3) = -1
      ahi(1) = nbf
      ahi(2) = nclosed(1)
      ahi(3) = 3
      if (.not.nga_create(MT_DBL,3,ahi,'U matrix',alo,g_u)) call 
     &    errquit('giao_b1: nga_create failed g_u',0,GA_ERR)
      call ga_zero(g_u)
c
c     Construction of right-hand side CPHF
c     Create CPHF array of proper dimension : (nocc*nvirt,3)
c
      if(.not.ga_create(MT_DBL,nclosed(1)*nvirt(1),3,'RHS',-1,-1,g_rhs))
     &   call errquit('giao_b1: ga_create failed g_rhs',0,GA_ERR)
      call ga_zero(g_rhs)
c
c     NGA dimension arrays for copying will be the same every time
c     Also third NGA dimension for any of the three dimensional
c     arrays will be the same everytime (running from 1 to 3)
c     So, lets define them once and for all in blo and bhi
c
      blo(1) = 1
      bhi(1) = nclosed(1)*nvirt(1)
      blo(2) = 1
      bhi(2) = 3
c    
c     Get S10 in GA and transform to MO set (virt,occ)
c
      alo(1) = nbf
      alo(2) = -1
      alo(3) = -1
      ahi(1) = nbf
      ahi(2) = nbf
      ahi(3) = 3
      if (.not.nga_create(MT_DBL,3,ahi,'s10 matrix',alo,g_s10)) call 
     &    errquit('giao_b1: nga_create failed g_s01',0,GA_ERR)
      call ga_zero(g_s10)
      call int_giao_1ega(basis,basis,g_s10,'s10',origin,
     &                   nat,oskel)
      call giao_aotomo(g_s10,vectors(1),nclosed,nvirt,1,3,nbf)

      if (debug) write (luout,*) 'S10 done'

c     while we are calculating integrals, let's also determine
c     the 'half' overlap derivative, used later in the calling
c     routine

      call ga_zero(g_sket1)
      call int_giao_1ega(basis,basis,g_sket1,'srxRb',origin,
     &   nat,oskel)

      if (debug) write (luout,*) 'S1-ket done'

c     g_sket1 will not be used further here. It is one of the 
c     output results of this routine.
c     Broceed with the computation of the B-field perturbed
c     MO coefficients.

c
c     ga_rhs(a,i) = ga_rhs(a,i) - e(i) * S10(a,i)
c     Scale (occ,virt) block g_s10 with - (minus) eigenvalues 
c
      alo(1) = nclosed(1)+1
      ahi(1) = nmo
      alo(3) = 1
      ahi(3) = 3
      do iocc = 1, nclosed(1)
         alo(2) = iocc
         ahi(2) = iocc
         call nga_scale_patch(g_s10,alo,ahi,-dbl_mb(k_eval+iocc-1)) 
      enddo
      if (.not.ma_pop_stack(l_eval)) call
     &    errquit('giao_b1: ma_pop_stack failed k_eval',0,MA_ERR)
      if (.not.ma_pop_stack(l_occ)) call
     &    errquit('giao_b1: ma_pop_stack failed k_occ',0,MA_ERR)
c
c     Copy to ga_rhs 
c     alo(1) and ahi(1) the same as before
c
      alo(2) = 1
      ahi(2) = nclosed(1)
      call nga_copy_patch('n',g_s10,alo,ahi,g_rhs,blo,bhi)

      if (debug) write (luout,*) 'S10 to Umat done'
c
c     Construct occ-occ part of the three U matrices
c     Occ-occ blocks for each field direction are defined as -1/2 S10
c     Scale (occ,occ) block g_s10 with -1/2 and add to g_u
c
c     alo(2) and ahi(2) will stay as 1 and nclosed(1) for a while
c
      alo(1) = 1
      ahi(1) = nclosed(1)
      call nga_scale_patch(g_s10,alo,ahi,-0.5d0)
      call nga_copy_patch('n',g_s10,alo,ahi,g_u,alo,ahi)

      if (debug) write (luout,*) 'S10 in occ-occ done'
c
c     We also need the occupied-occupied contribution of g_u
c     contributing to the first order density matrix. As this block does
c     not change during the CPHF we can calculate it once and subtract
c     it from the RHS. We will reuse g_s10 as scratch space.
c
      call ga_zero(g_s10)
      clo(1) = 3
      clo(2) = nbf
      clo(3) = nbf
      chi(1) = 1  
      chi(2) = -1 
      chi(3) = -1
      if (.not.nga_create(MT_DBL,3,clo,'Fock matrix',chi,g_fock)) call 
     &    errquit('giao_b1: nga_create failed g_fock',0,GA_ERR)
      if (.not.nga_create(MT_DBL,3,clo,'D10 matrix',chi,g_d1)) call 
     &    errquit('giao_b1: nga_create failed g_d1',0,GA_ERR)
      call ga_zero(g_fock)
      call ga_zero(g_d1)
      alo(1) = 1
      alo(2) = 1
      blo(1) = 1
      blo(2) = 1
      bhi(1) = nbf
      clo(2) = 1
      clo(3) = 1
      chi(2) = nbf
      chi(3) = nbf
      dlo(1) = 1
      dlo(2) = 1
      dhi(1) = nbf
      dhi(2) = nclosed(1)
c
c     Create "perturbed density matrix" for closed-closed g_u block
c
      do ifld = 1, 3
         jfac(ifld) = 0.0d0
         kfac(ifld) = -1.0d0*xfac
         alo(3) = ifld
         ahi(3) = ifld
         clo(1) = ifld
         chi(1) = ifld
         dlo(3) = ifld
         dhi(3) = ifld
         ahi(1) = nmo
         ahi(2) = nclosed(1)
         bhi(2) = nmo 
         call nga_matmul_patch('n','n',1.0d0,0.0d0,vectors(1),blo,bhi,  
     &                         g_u,alo,ahi,g_s10,dlo,dhi)  
         ahi(2) = nbf
         ahi(1) = nclosed(1)
         bhi(2) = nclosed(1)
c
c     Minus sign as we subtract it from the RHS as we do not include 
c     it in the LHS
c
         call nga_matmul_patch('n','t',-1.0d0,0.0d0,vectors(1),blo,bhi,
     &                         g_s10,alo,ahi,g_d1,clo,chi)  
      enddo
c
c     Build "first order fock matrix"
c
      if (use_theory.eq.'dft') then
         if(.not. rtdb_put(rtdb,'bgj:xc_active', MT_LOG, 1, .true.))
     $     call errquit('hess_cphf: rtdb_put of xc_active failed',0,
     &       RTDB_ERR)
         if(.not. rtdb_put(rtdb,'fock_xc:calc_type', MT_INT, 1, 2))
     $     call errquit('hess_cphf: rtdb_put of calc_type failed',0,
     &       RTDB_ERR)
         if(.not. rtdb_put(rtdb,'fock_j:derfit', MT_LOG, 1, .false.))
     $     call errquit('hess_cphf: rtdb_put of j_derfit failed',0,
     &       RTDB_ERR)
      endif
c      call shell_fock_build(geom, basis, 0, 3,
c     $     jfac, kfac,tol2e, g_d1, g_fock,.false.)
c      if(use_theory.eq.'dft') then
c         if (.not. rtdb_put(rtdb, 'fock_xc:calc_type', mt_int, 1, 0))
c     $      call errquit('giaox: rtdb_put failed',0,RTDB_ERR)
c      endif

c from hnd_giaox.F]
c
c     Note: Just the exchange: jfac = 0.d0 (see above)
c
      if (.not.cam_exch) then
         call shell_fock_build(geom, basis, 0, 3,
     $     jfac, kfac, tol2e, g_d1, g_fock, .false.)
      else
         call shell_fock_build_cam(geom, basis, 0, 3,
     $     jfac, kfac, tol2e, g_d1, g_fock, .false.)
      end if
c
      if(use_theory.eq.'dft') then
         if (.not. rtdb_put(rtdb, 'fock_xc:calc_type', mt_int, 1, 0))
     $      call errquit('giaox: rtdb_put failed',0,RTDB_ERR)
      endif


c
c     Transform to the occ-virt MO basis and add to RHS
c
      call ga_zero(g_s10)
      alo(1) = 1
      ahi(1) = nbf
      alo(2) = 1
      ahi(2) = nclosed(1)
      do ifld = 1, 3
         alo(3) = ifld
         ahi(3) = ifld
         clo(1) = ifld
         chi(1) = ifld
         call nga_matmul_patch('n','n',2.0d0,0.0d0,
     $        g_fock,clo,chi,
     $        vectors(1),alo,ahi,
     $        g_s10,alo,ahi)
      enddo
      call ga_zero(g_fock)
      clo(2) = nclosed(1)+1
      clo(3) = 1
      chi(2) = nmo
      chi(3) = nclosed(1)
      do ifld = 1, 3
         blo(1) = nclosed(1)+1
         blo(2) = 1
         bhi(1) = nmo
         bhi(2) = nbf
         alo(3) = ifld
         ahi(3) = ifld
         clo(1) = ifld
         chi(1) = ifld
         call nga_matmul_patch('t','n',1.0d0,0.0d0,
     $        vectors(1), blo,bhi,
     $        g_s10, alo,ahi,
     $        g_fock, clo,chi )
         blo(1) = 1
         bhi(1) = nclosed(1)*nvirt(1)
         blo(2) = ifld
         bhi(2) = ifld
         call nga_add_patch(1.0d0,g_rhs,blo,bhi,1.0d0,g_fock,clo,chi,
     &                      g_rhs,blo,bhi)

!      call nga_print_patch(g_fock,clo,chi,1)

      enddo
c
c     Cleanup of g_d1 and g_fock, not needed for now
c
      if (.not.ga_destroy(g_d1)) call 
     &    errquit('giao_b1: ga_destroy failed g_d1',0,GA_ERR)
      if (.not.ga_destroy(g_fock)) call 
     &    errquit('giao_b1: ga_destroy failed g_fock',0,GA_ERR)

      if (debug) write (luout,*) 'S10 Fockop done'
c
c     Get H10 in GA, reusing g_s10 array
c
      call ga_zero(g_s10)
      call int_giao_1ega(basis,basis,g_s10,'l10',origin,
     &                   nat,oskel)
      call int_giao_1ega(basis,basis,g_s10,'tv10',origin,
     &                   nat,oskel)
c
c     Get external and cosmo bq contribution
c
      nbq = 0
      nextbq = 0
      ncosbq = 0
      if(geom_extbq_on()) nextbq = geom_extbq_ncenter()
      nbq = nextbq ! external bq's
c
      if (rtdb_get(rtdb,'cosmo:nefc',mt_int,1,ncosbq))
     &    nbq = ncosbq ! cosmo bq's
c
      if (nextbq.gt.0.and.ncosbq.gt.0)
     &    nbq = nextbq + ncosbq  ! tally up cosmo and external bqs
c
c     if (ga_nodeid().eq.0) write(6,*) "nbq: ", nbq
      if (nbq.gt.0) then
        call int_giao_1ega(basis,basis,g_s10,'bq10',origin,
     &                   nat,oskel)
      end if
c
c     ga_rhs(a,i) = ga_rhs(a,i) + H10(a,i)
c     Transform H10 to MO and add to g_rhs
c
      call giao_aotomo(g_s10,vectors(1),nclosed,nvirt,1,3,nbf)
      alo(1) = nclosed(1)+1
      ahi(1) = nmo
      alo(2) = 1
      ahi(2) = nclosed(1)
      alo(3) = 1
      ahi(3) = 3
      blo(1) = 1
      bhi(1) = nclosed(1)*nvirt(1)
      blo(2) = 1
      bhi(2) = 3
      call nga_add_patch(1.0d0,g_rhs,blo,bhi,1.0d0,g_s10,alo,ahi,
     &                   g_rhs,blo,bhi)
c
c     Cleanup g_s10 as we do not need it right now
c
      if (.not.ga_destroy(g_s10)) call 
     &    errquit('giao_b1: ga_destroy failed g_s10',0,GA_ERR)

      if (debug) write (luout,*) 'H10 to Umat done'
c
c     Remaining term is Perturbed (GIAO) two-electron term times
c     Unperturbed density Calculate Sum(r,s) D0(r,s) * G10(m,n,r,s) in
c     AO basis
c
      alo(1) = -1 
      alo(2) = -1
      alo(3) = 1
      ahi(1) = nbf
      ahi(2) = nbf
      ahi(3) = 3
      if (.not.nga_create(MT_DBL,3,ahi,'Fock matrix',alo,g_fock)) call 
     &    errquit('giao_b1: nga_create failed g_fock',0,GA_ERR)
      call ga_zero(g_fock)
      if(use_theory.eq.'dft') then
         ifld = 4
         if (.not. rtdb_put(rtdb, 'fock_xc:calc_type', mt_int, 1, ifld))
     $      call errquit('giao_b1: rtdb_put failed',0,RTDB_ERR)
      endif
      if (debug) write (luout,*) 'calling new_giao_2e'
c     call new_giao_2e(geom, basis, nbf, tol2e, g_dens, g_fock, xfac)
      call new_giao_2e(geom, basis, nbf, tol2e, g_dens, g_fock, xfac,1) ! FA-restricted calc
      if (debug) write (luout,*) 'done new_giao_2e'
      if(use_theory.eq.'dft') then
         ifld = 0
         if (.not. rtdb_put(rtdb, 'fock_xc:calc_type', mt_int, 1, ifld))
     $      call errquit('giao_b1: rtdb_put failed',0,RTDB_ERR)
         if(.not. rtdb_put(rtdb,'bgj:xc_active', MT_LOG, 1, .false.))
     $       call errquit('giao_b1: rtdb_put of xc_active failed',0,
     &       RTDB_ERR)
      endif
c
c     Transform to MO basis and add to right-hand-side
c
      call giao_aotomo(g_fock,vectors(1),nclosed,nvirt,1,3,nbf)
      alo(1) = nclosed(1)+1
      ahi(1) = nmo
      alo(2) = 1
      ahi(2) = nclosed(1)
      alo(3) = 1
      ahi(3) = 3
      call nga_add_patch(1.0d0,g_rhs,blo,bhi,1.0d0,g_fock,alo,ahi,
     &                   g_rhs,blo,bhi)
      if (.not.ga_destroy(g_fock)) call 
     &    errquit('giao_b1: ga_destroy failed g_fock',0,GA_ERR)
      call nga_scale_patch(g_rhs,blo,bhi,-4.0d0)

      if (debug) write (luout,*) 'D00 GIAO Fock terms done'
c
c     Write ga_rhs to disk 
c
      call util_file_name('cphf_rhs',.true.,.true.,cphf_rhs)
      call util_file_name('cphf_sol',.true.,.true.,cphf_sol)
      if(.not.file_write_ga(cphf_rhs,g_rhs)) call errquit
     $  ('giao_b1: could not write cphf_rhs',0, DISK_ERR)
c
      call schwarz_tidy()
      call int_terminate()
c
c     Call the CPHF routine
c     
c     We do need to tell the CPHF that the density is skew symmetric.
c     Done via rtdb, put cphf:skew .false. on rtdb and later remove it.
c
      if (.not. rtdb_put(rtdb, 'cphf:skew', mt_log, 1,.false.)) call
     $   errquit('giao_b1: failed to write skew ', 0, RTDB_ERR)

      if (debug) write (luout,*) 'calling cphf'
c
      if (.not.cphf2(rtdb)) call errquit
     $  ('giao_b1: failure in cphf ',0, RTDB_ERR)
c
      if (.not. rtdb_delete(rtdb, 'cphf:skew')) call
     $   errquit('giao_b1: rtdb_delete failed ', 0, RTDB_ERR)

      if (debug) write (luout,*) 'cphf done'
c
c     Occ-virt blocks are the solution pieces of the CPHF
c     Read solution vector from disk and put solutions in U matrices
c
      call ga_zero(g_rhs)
      if(.not.file_read_ga(cphf_sol,g_rhs)) call errquit
     $  ('giao_b1: could not read cphf_rhs',0, DISK_ERR)      
      call nga_copy_patch('n',g_rhs,blo,bhi,g_u,alo,ahi)
c
      if (.not.ga_destroy(g_rhs)) call 
     &    errquit('giao_b1: ga_destroy failed g_rhs',0,GA_ERR)
c
c     From U matrices, generate the perturbed density matrices D1x,y,z
c     C1 = C0 * U10
c     D1 = 2[(C1*C0+) - (C0*C1+)]
c
      alo(1) = nbf
      alo(2) = -1
      alo(3) = -1
      ahi(1) = nbf
      ahi(2) = nbf
      ahi(3) = 3
      if (.not.nga_create(MT_DBL,3,ahi,'D10 matrix',alo,g_d1)) call 
     &    errquit('giao_b1: nga_create failed g_d1',0,GA_ERR)
c
      alo(1) = 1
      alo(2) = 1
      blo(1) = 1
      blo(2) = 1
      clo(1) = 1
      chi(1) = nbf
      clo(2) = 1
      chi(2) = nbf
      dlo(1) = 1
      dlo(2) = 1
      dhi(1) = nbf
      dhi(2) = nclosed(1)
      do ifld = 1, 3
         alo(3) = ifld
         ahi(3) = ifld
         blo(3) = ifld
         bhi(3) = ifld
         clo(3) = ifld
         chi(3) = ifld
         dlo(3) = ifld
         dhi(3) = ifld
         bhi(1) = nbf
         bhi(2) = nmo 
         ahi(1) = nmo
         ahi(2) = nclosed(1)
c
c     Make C1
c
         call nga_matmul_patch('n','n',1.0d0,0.0d0,vectors(1),blo,bhi,  
     &                         g_u,alo,ahi,g_d1,dlo,dhi)  
         call nga_copy_patch('n',g_d1,dlo,dhi,g_vecB1,dlo,dhi)

c        This patch of g_vecB1 now has the perturbed MO
c        coefficients. let's print them for debug purposes:

         if (.false.) then ! change .false. to debug to print them
           write (luout,*) 'giao_b1: perturbed C, direction ',ifld
           call nga_print_patch(g_vecB1,dlo,dhi,1)
         end if

       enddo ! ifld

      if (.not.ga_destroy(g_u)) call 
     &    errquit('giao_b1: ga_destroy failed g_u',0,GA_ERR)
      if (.not.ga_destroy(g_d1)) call 
     &    errquit('giao_b1: ga_destroy failed g_d1',0,GA_ERR)
      if (.not.ga_destroy(vectors(1))) call 
     &   errquit('giao_b1: ga_destroy failed vectors',0,GA_ERR)
      if (.not.ga_destroy(g_dens(1))) call 
     &    errquit('giao_b1: ga_destroy failed g_dens',0,GA_ERR)

      if (debug) write (luout,*) 'Perturbed C and Pmat done'

c     Alll done.
c     At this point, we don't need to terminate the integrals.
c     They were terminated by the cphf.

      call ga_sync()

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

