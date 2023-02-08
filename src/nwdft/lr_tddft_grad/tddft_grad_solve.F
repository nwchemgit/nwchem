c
c     Solution of the equation (A+B)Z=-R (Eq. 18, Furche & Alrichs)
c     HvD 11/2007, NG 11/2012
c
      subroutine tddft_grad_solve(rtdb,ihdl_geom,ihdl_bfao,ipol,nroots,
     +           nao,nocc,naoc,nav,nfc,nfv,g_mo,iptr_mo_e,tda,oskel,
     +           kfac,lhashf,otriplet,solve_thresh,solve_maxsub,
     +           solve_maxit,fn_solve,solve_when,tol2e,g_r,g_z,lhascd)
c
c     At this point Ria has just been computed. The quantities
c     are stored in 3D global arrays as R(1)(1:nroot,1:nocca,1:nvirta),
c     R(2)(1:nroot,1:noccb,1:nvirtb). 
c     The NWChem CPHF solver GA_LKAIN assumes that the quantities are 
c     stored in the form R(1:nocca*nvirta+noccb*nvirtb,1:nroot). 
c
c     Furthermore GA_LKAIN takes a subroutine for the matrix vector
c     products as argument as well as a subroutine for the 
c     preconditioner.
c
c     The matrix vector routine has the interface:
c
c        subroutine product(accuracy,g_x,g_Ax)
c
c     And the preconditioner:
c
c        subroutine precond(g_x,shift)
c
c     See the subroutine GA_LKAIN (nwchem/src/util/ga_it2.F) for more
c     details. In the product routine the quantities will need to be 
c     transformed to AO basis and back to be able to use the subroutine 
c     tddft_nga_cont. The most efficient way to do this is to use the 
c     transformation tddft_grad_trans_mo2ao and tddft_grad_trans_ao2mo.
c     These however expect 3-dimensional global arrays again. The
c     "transformation" between 2D and 3D representations can be handled
c     through the routines ga_mat_to_vec and ga_vec_to_mat.
c
c     So the overall scheme becomes:
c
c       1. Transform R 3D->2D
c       2. Call ga_lkain
c       2.1  Call product
c       2.1.1  Transform g_x 2D->3D
c       2.1.2  Transform g_x 3D MO->AO
c       2.1.3  Call tddft_nga_cont -> g_Ax 3D AO
c       2.1.4  Transform g_Ax 3D AO->MO
c       2.1.5  Transform g_Ax 3D->2D
c       3. Transform Z 2D->3D
c
c     This is a bit tedious but should work...
c
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "global.fh"
#include "tddft_grad_util.fh"
#include "stdio.fh"
#include "tddft_grad_solve_common.fh"
#include "rtdb.fh"
c
      logical odiff,oprint
      parameter (odiff=.true.)
c
c     Input:
c
      integer rtdb        ! the runtime database handle
      integer ipol        ! =1 (restricted), =2 (unrestricted)
      logical tda         ! True if Tamm-Dancoff approximation
      logical oskel       ! Use of symmetry in Fock build
      logical lhashf      ! has the functional got HF exchange
      logical otriplet    ! are we doing triplet excitations
      logical lhascd
      integer nroots      ! the number of roots
      integer nao         ! the number of AOs
      integer nfc(2)      ! the number of frozen core orbitals
      integer nfv(2)      ! the number of frozen virtual orbitals
      integer nocc(2)     ! the number of occupied orbitals
      integer naoc(2)     ! the number of active occupied orbitals
      integer nav(2)      ! the number of active virtual orbitals
      integer iptr_mo_e   ! mo energy memory index
      integer ihdl_bfao   ! the AO basis set handle
      integer ihdl_geom   ! the geometry handle
c
      double precision kfac  ! Prefactor for HF exchange
      double precision tol2e ! 2-electron integral tolerance
c
      double precision solve_thresh ! convergence threshold for solver
c
      integer solve_maxit  ! the maximum number of iterations in solver
      integer solve_maxsub ! the maximum number expansion vectors
                           ! per root
c
      character*(*) fn_solve   ! filename for linear system solutions
      character*(*) solve_when ! when to store solutions
c
      integer g_mo(2)      ! global arrays holding MOs
      integer g_r(2)       ! linear system RHS
c
c     Output:
c
      integer g_z(2)       ! linear system solutions
c
c     Subroutines:
c
      external tddft_grad_solve_product
      external tddft_grad_solve_precond
c
      logical  tddft_grad_solve_load_x
      external tddft_grad_solve_load_x
c
c     Local:
c
      integer idim,jdim   ! GA dimensions
      integer ichk,jchk   ! GA chunk sizes
      integer g_b         ! 2D GA for right-hand-sides
      integer g_x         ! 2D GA for solutions
      integer klo         ! location in vectors in 2D GA
      integer col         ! number of columns in 2D GA
      integer alo(3)      ! block lower limits in 3D GA
      integer ahi(3)      ! block upper limits in 3D GA
      integer ip          ! counter over spin components
      integer ir          ! counter over roots
      logical oroot
c
      character*32 pname
      logical tdaloc, tddft_grad_tdaloc
      external tddft_grad_tdaloc
      logical oprint_debug
c
      pname= "tddft_grad_solve: "
      oprint=util_print('tddft_grad_solve',print_medium)
      oprint_debug=util_print('tddft_grad_solvecphf',print_debug)
c
c     Initialize the common block for tddft_grad_solve_product and
c     tddft_grad_solve_precond
c
      oroot = ga_nodeid().eq.0
c
      tdaloc=tddft_grad_tdaloc(rtdb,tda,lhashf)
c End test
      tdg_tol2e   = tol2e
      tdg_kfac    = kfac
      tdg_rtdb    = rtdb
      tdg_geom    = ihdl_geom
      tdg_bfao    = ihdl_bfao
      tdg_ipol    = ipol
      tdg_nroots  = nroots
      tdg_nao     = nao
      tdg_nocc(1) = nocc(1)
      tdg_nocc(2) = nocc(2)
      tdg_naoc(1) = naoc(1)
      tdg_naoc(2) = naoc(2)
      tdg_nav(1)  = nav(1)
      tdg_nav(2)  = nav(2)
      tdg_nfc(1)  = nfc(1)
      tdg_nfc(2)  = nfc(2)
      tdg_nfv(1)  = nfc(1)
      tdg_nfv(2)  = nfc(2)
      tdg_nmo(1)  = nfc(1)+naoc(1)+nav(1)+nfv(1)
      tdg_nmo(2)  = nfc(2)+naoc(2)+nav(2)+nfv(2)
      tdg_g_mo(1) = g_mo(1)
      tdg_g_mo(2) = g_mo(2)
      tdg_e_mo    = iptr_mo_e
      tdg_tda     = tdaloc
      tdg_oskel   = oskel
      tdg_lhashf  = lhashf
cold  tdg_triplet = otriplet
      tdg_triplet = .false.
      tdg_lhascd  = lhascd
c Daniel (2-26-13): It was not obvious that we need to unset
c fock_xc:triplet here for restricted triplet calculations to work.
      if (otriplet) then
        if (.not.rtdb_put(rtdb,'fock_xc:triplet',mt_log,1,.false.))
     1    call errquit(pname//'failed to set triplet',0,RTDB_ERR)
      endif    
c
c     g_r 3D MO -> g_b 2D MO
c
      idim = naoc(1)*nav(1)+(ipol-1)*naoc(2)*nav(2)
      jdim = nroots
      ichk = -1
      jchk = -1
      if (.not.ga_create(mt_dbl,idim,jdim,"g_b",ichk,jchk,g_b)) 
     +  call errquit(pname//'failed to create g_b',0, GA_ERR)
      if (.not.ga_create(mt_dbl,idim,jdim,"g_x",ichk,jchk,g_x)) 
     +  call errquit(pname//'failed to create g_x',0, GA_ERR)
c
      do ip = 1, tdg_ipol
        do ir = 1, tdg_nroots
          alo(1) = ir
          ahi(1) = ir
          alo(2) = 1
          ahi(2) = tdg_naoc(ip)
          alo(3) = 1
          ahi(3) = tdg_nav(ip)
          klo = (ip-1)*(tdg_naoc(1)*tdg_nav(1))+1
          col = ir
          call tddft_grad_mat_to_vec(g_r(ip),alo,ahi,g_b,klo,col,
     +                               -1.0d0,"=")
        enddo
      enddo
c
c     Solve the linear system of equations
c
      if (.not.tddft_grad_solve_load_x(g_x,ipol,nao,nfc,naoc,nocc,nav,
     +         nfv,jdim,g_mo,fn_solve)) then
        call ga_copy(g_b,g_x)
        call tddft_grad_solve_precond(g_x,0.0d0)
      endif
      if(oprint_debug) then
         if(oroot) write(luout,*) ' LHS'
         call util_flush(luout)
         call ga_print(g_x)
         if(oroot) write(luout,*) ' RHS'
         call util_flush(luout)
         call ga_print(g_b)
      endif
c
      call ga_lkain(rtdb,g_x,g_b,tddft_grad_solve_product,
     +              tddft_grad_solve_precond,solve_thresh,
     +              solve_maxsub*nroots,solve_maxit,odiff,oprint)
c
      call tddft_grad_solve_save_x(g_x,ipol,nao,nfc,naoc,nocc,nav,nfv,
     +              jdim,g_mo,fn_solve,solve_when,ihdl_geom,ihdl_bfao,
     +              oskel)
c
c     g_x 2D MO -> g_z 3D MO
c
      do ip = 1, tdg_ipol
        do ir = 1, tdg_nroots
          alo(1) = ir
          ahi(1) = ir
          alo(2) = 1
          ahi(2) = tdg_naoc(ip)
          alo(3) = 1
          ahi(3) = tdg_nav(ip)
          klo = (ip-1)*(tdg_naoc(1)*tdg_nav(1))+1
          col = ir
          call tddft_grad_vec_to_mat(g_z(ip),alo,ahi,g_x,klo,col)
        enddo
      enddo
c
      if (.not.ga_destroy(g_b)) 
     +   call errquit(pname//'failed to destroy g_b',0, GA_ERR)
      if (.not.ga_destroy(g_x)) 
     +   call errquit(pname//'failed to destroy g_x',0, GA_ERR)
c
      if (util_print('tddft grad z',print_debug)) then
        if (oroot) write(LuOut,*)'DEBUG: '//pname//'Z'
        call tddft_grad_print_array(ipol,nroots,g_z,dble(ipol))
      endif
c
c     Reset fock_xc:triplet here for restricted triplet calculations
      if (otriplet) then
        if (.not.rtdb_put(rtdb,'fock_xc:triplet',mt_log,1,.true.))
     1    call errquit(pname//'failed to set triplet',0,RTDB_ERR)
      endif
c
      end
c  We manipulate the code here because the
c R vector has the same number of terms for RPA and CIS.  This is a
c consequence of (X-Y) = X.  It might be a good idea to avoid doing this
c part for the Y vector since Y = 0.  Note that the coupling matrix
c expressions H^+[V] and H^-[V] can both be nonzero for CIS, so it isn't
c okay to skip the anti-symmetric part in the tddft_nga_cont routine.  
c What CIS does is makes the Y vector contribution zero in the following
c routines.
c Set the local TDA variable so that we don't change
c the global one.
      logical function tddft_grad_tdaloc(rtdb,tda,lhashf)
      implicit none
#include "rtdb.fh"
#include "mafdecls.fh"
      integer rtdb
      logical tda
      logical lhashf
c
      logical tdaloc_rtdb
c      
      if (lhashf) then
        if (tda) then
          tddft_grad_tdaloc = .false. ! For CIS calculations
        else
          tddft_grad_tdaloc = .false. ! For RPA calculations
        endif
      else 
        if (tda) then
cold          tddft_grad_tdaloc = .true.  ! For TDDFT/TDA calculations
          tddft_grad_tdaloc = .false.  ! For TDDFT/TDA calculations
        else
          tddft_grad_tdaloc = .false. ! For TDDFT calculations
        endif
      endif
      if (rtdb_get(rtdb,'tddft:tdaloc',mt_log,1,tdaloc_rtdb))
     +     tddft_grad_tdaloc=tdaloc_rtdb
      return
      end
