      subroutine wfn1_input(rtdb)
      implicit none
!>
!> \brief Read the input for the RDMFT calculations
!>
      integer rtdb !> [Input] The runtime database
      end subroutine wfn1_input
!>
!>----------------------------------------------------------------------
!>
      logical function wfn1_task_energy(rtdb)
      implicit none
#include "wfn1_param.fh"
#include "wfn1_wfn.fh"
#include "wfn1_drv.fh"
#include "wfn1_matcnst.fh"
#include "wfn1_size.fh"
#include "wfn1_trace.fh"
#include "wfn1_en.fh"
#include "bas.fh"
#include "errquit.fh"
#include "geom.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "util.fh"
!>
!> \brief Minimize the free energy
!>
!> The wavefunction we are using is outlined in [1]. The energy
!> expression is quite different though. The reason is that attempts to
!> reconstruct a valid 2-electron density matrix from an arbitrary
!> 1-electron density matrix never work. Instead when starting from 
!> any wavefunction the 2-electron density matrix is N-representable
!> by construction (this includes the wavefunction from [1]). Hence
!> it would seem that the correlation energy can be more easily
!> formulated in other ways than in terms of 1-electron density
!> matrices. For this reason RDMFT was abandoned and we explore a
!> rather different approach here.
!>
      integer rtdb !> [Input] The runtime database
!>
!> \return `true`: if the energy minimizartion converged, and `false`
!> !otherwise.
!>
!> ### References ###
!>
!> [1] Hubertus van Dam,
!>     "Generalization of the Kohn-Sham system enabling representing
!>     arbitary one electron density matrices",
!>     arXiv: <a href="http://arxiv.org/abs/1512.06120">
!>     http://arxiv.org/abs/1512.06120</a>
!>
c     Local variables
c
      logical lstatus ! the status of the runtime database
      logical oprint  ! should something be printed
      integer iproc   ! the processor rank
      integer nproc   ! the number of processors
      integer nelec   ! the total number of electrons
      integer nopen   ! the number of unpaired electrons
      integer nclosed ! the number of closed shell orbitals
      integer maxit   ! the maximum number of iterations
      integer curit   ! the current iteration number
c
      logical oconverged ! is this calculation converged
c
      double precision charge     ! the molecular charge
      double precision charge_nuc ! the total nuclear charge
c
      character*12  pname ! the subroutine name
      character*256 title ! the job title
c
      integer g_bfbf     ! nbf x nbf scratch global array
      integer g_bfmo     ! nbf x nmo scratch global array
      integer g_momo1    ! nmo x nmo scratch global array
      integer g_momo2    ! nmo x nmo scratch global array
c
      integer ndep_nmo          ! the number of linearly independent
                                ! vectors set as a user parameter
      double precision ndep_tol ! the tolerance on the eigenvalues of
                                ! the overlap matrix set as a user
                                ! parameter
      double precision ndep_tol_default ! default setting for ndep_tol
      parameter (ndep_tol_default = 1.0d-5)
c
      double precision conv_tol ! the convergence tolerance
c
      double precision, allocatable :: eig(:)
      double precision              :: error
c
      type(wfn1_prminst) :: wfn1_inst    ! paramater instances
      type(wfn1_prm)     :: wfn1_param   ! calculation parameters
      type(wfn1_matcnst) :: wfn1_mat     ! matrix constants
      type(wfn1_wfn)     :: wfn1_wave    ! the current wave function
      type(wfn1_drv)     :: wfn1_deriv   ! the derivative of the current
                                         ! wave function
      type(wfn1_en)      :: wfn1_energy  ! the energy terms
c
      type(wfn1_wfn)     :: wfn1_trial_wv  ! the trial wave function
      type(wfn1_drv)     :: wfn1_trial_drv ! the derivative wrt the trial
                                           ! wave function
      type(wfn1_en)      :: wfn1_trial_en  ! the trial energy
c
      type(wfn1_wfn) :: wfn1_t_wv  ! temporary variable
      type(wfn1_drv) :: wfn1_t_drv ! temporary variable
c
      iproc = ga_nodeid()
      nproc = ga_nnodes()
c
      wfn1_inst%rtdb = rtdb
      pname = "wfn1_energy:"
      wfn1_task_energy = .false.
      lstatus = rtdb_parallel(.true.)
      call util_print_push
      call util_print_rtdb_load(wfn1_inst%rtdb, 'wfn1') 
c
c     == Print Header ==
c
      oprint = util_print('information',print_low)
      if (oprint) then
        if (.not. rtdb_cget(wfn1_inst%rtdb, 'title', 1, title))
     &       title = ' '
        if(iproc.eq.0) then
          write(LuOut,*)
          call util_print_centered(
     &         LuOut, 'NWChem WFN1 Module', 40, .true.)
          write(LuOut,*)
          write(LuOut,*)
          if (title .ne. ' ') then
            call util_print_centered(LuOut, title, 40, .false.)
            write(LuOut,*)
            write(LuOut,*)
          endif
        endif
      endif
c
c     == Load geometry and basis set ==
c
      if (.not.geom_create(wfn1_inst%geom,'geometry'))
     &  call errquit(pname//' geom_create failed',0,GEOM_ERR)
      if (.not.geom_rtdb_load(wfn1_inst%rtdb,wfn1_inst%geom,'geometry'))
     &  call errquit(pname//' no geometry found',0,RTDB_ERR)
      if (.not.bas_create(wfn1_inst%basis,'ao basis'))
     &  call errquit(pname//' bas_create failed',0, BASIS_ERR)
      if (.not.bas_rtdb_load(wfn1_inst%rtdb,wfn1_inst%geom,
     &                       wfn1_inst%basis,'ao basis'))
     &  call errquit(pname//' no ao basis found',0, RTDB_ERR)
      if (.not.bas_numbf(wfn1_inst%basis,wfn1_param%nbf))
     &  call errquit(pname//" no basis functions",0,BASIS_ERR)
c
c     == Load linearly independence settings ==
c
      if (rtdb_get(rtdb, 'lindep:n_dep', MT_INT, 1, ndep_nmo)) then
        if (oprint) then
          write(6, 2) ndep_nmo
 2        format(/' Forcing the number of linearly dependent ',
     &            'vectors to be', i6)
          call util_flush(6)
        endif
      else
        ndep_nmo = -1
      endif
      if (rtdb_get(rtdb, 'lindep:tol', MT_DBL, 1, ndep_tol)) then
        if (oprint) then
          write(6, 3) ndep_tol
 3        format(/' Setting tolerance for linear dependence to',
     &          1p,d9.2)
          call util_flush(6)
        endif
      else
         ndep_tol = -1.0d0
      endif
c
      if ((ndep_tol.ne.-1.0d0) .and. (ndep_nmo.ne.-1)) call errquit
     &  (pname//' only one of lindep:n_dep and lindep:tol permitted',0,
     &  UNKNOWN_ERR)
c
      if (ndep_nmo .eq. -1) then
        if (ndep_tol .le. -1.0d0) ndep_tol = ndep_tol_default
      endif
c
c     == Print the geometry ==
c
      if(iproc.eq.0.and.util_print('geometry', print_high)) then
        if(.not.geom_print(wfn1_inst%geom)) call errquit(
     &       pname//'geom print failed',0, GEOM_ERR)
      endif
c
c     == Evaluate the nuclear repulsion energy ==
c
      lstatus = geom_nuc_rep_energy(wfn1_inst%geom, wfn1_param%Enuc)
      if (.not.lstatus)
     &  call errquit(pname//'Error in nuclear repulsion',1, GEOM_ERR)
c
c     == Work out numbers of electrons ==
c
      if (.not.geom_nuc_charge(wfn1_inst%geom,charge_nuc))
     &  call errquit(pname//" geometry nuclear charge failed",
     &               0,GEOM_ERR)
      if (.not.rtdb_get(wfn1_inst%rtdb,"charge",MT_DBL,1,charge))
     &  charge = 0.0d0
      if (.not. rtdb_get(wfn1_inst%rtdb,'wfn1:nopen',MT_INT,1,nopen))
     &  nopen = 0
      nelec = int(charge_nuc - charge)
      if (nelec.lt.0) then
        call errquit(pname//" negative number of electrons",
     &               nelec,INPUT_ERR)
      endif
      nclosed =  (nelec-nopen)/2
      wfn1_param%nea = nclosed + nopen
      wfn1_param%neb = nclosed
c
c     == Get the convergence criterion
c
      if (.not.rtdb_get(wfn1_inst%rtdb,"conv_tol",MT_DBL,1,conv_tol))
     &  conv_tol = 1.0d-4
c
c     == Deal with 1-electron integrals ==
c
      call wfn1_allocate_mat(wfn1_inst,wfn1_mat)
c
c     == Initialize the integral code ==
c
      call int_init(wfn1_inst%rtdb,1,wfn1_inst%basis)
c     call schwarz_init(geom,basis)
c     call int_1e_cache_ga(wfn1_inst%basis,oskel)
      call wfn1_gen_mat(wfn1_inst,wfn1_mat)
c
c     == create scratch global arrays that depend on nbf only ==
c
      if (.not.ga_create(MT_DBL,wfn1_param%nbf,wfn1_param%nbf,"g_bfbf",
     &                   -1,-1,g_bfbf)) then
        call errquit(pname//" allocate g_bfbf failed",
     &               wfn1_ibd_proc(wfn1_param%nbf**2,nproc),GA_ERR)
      endif
      if (.not.ga_create(MT_DBL,wfn1_param%nbf,wfn1_param%nmo,"g_bfmo",
     &                   -1,-1,g_bfmo)) then
        call errquit(pname//" allocate g_bfmo failed",
     &               wfn1_ibd_proc(wfn1_param%nbf*wfn1_param%nmo,nproc),
     &               GA_ERR)
      endif
      if (.not.ga_create(MT_DBL,wfn1_param%nmo,wfn1_param%nmo,"g_momo1",
     &                   -1,-1,g_momo1)) then
        call errquit(pname//" allocate g_momo1 failed",
     &               wfn1_ibd_proc(wfn1_param%nmo**2,nproc),GA_ERR)
      endif
      if (.not.ga_create(MT_DBL,wfn1_param%nmo,wfn1_param%nmo,"g_momo2",
     &                   -1,-1,g_momo2)) then
        call errquit(pname//" allocate g_momo2 failed",
     &               wfn1_ibd_proc(wfn1_param%nmo**2,nproc),GA_ERR)
      endif
c
c     == Find the linearly independent orbitals ==
c     == and orthonormalizing transformation    ==
c
      call wfn1_lindep(wfn1_mat%g_s,g_bfbf,wfn1_param%nmo,ndep_tol)
      if (ndep_nmo.gt.0) wfn1_param%nmo = min(wfn1_param%nbf,ndep_nmo)
c
c     == Allocate memory for the wavefunction ==
c
      call wfn1_allocate_wfn(wfn1_param,wfn1_inst,wfn1_wave)
      call wfn1_allocate_wfn(wfn1_param,wfn1_inst,wfn1_trial_wv)
c
c     == Allocate memory for the derivatives wrt the wavefunctions ==
c
      call wfn1_allocate_drv(wfn1_param,wfn1_deriv)
      call wfn1_allocate_drv(wfn1_param,wfn1_trial_drv)
c
c     == Populate the natural orbitals with the eigenfunctions of ==
c     == the S^{-1/2} matrix                                      == 
c
      call ga_copy_patch('N',
     &     g_bfbf,1,wfn1_param%nbf,
     &            wfn1_param%nbf-wfn1_param%nmo+1,wfn1_param%nbf,
     &     wfn1_wave%g_na,1,wfn1_param%nbf,1,wfn1_param%nmo)
      call ga_copy(wfn1_wave%g_na,wfn1_wave%g_nb)
c
      call guess_dens(wfn1_inst%rtdb,wfn1_inst%geom,wfn1_inst%basis,
     &                g_bfbf)
      allocate(eig(1:wfn1_param%nbf))
      call rhf_dens_to_mo(wfn1_inst%rtdb,wfn1_inst%geom,wfn1_inst%basis,
     &                    wfn1_param%nea+wfn1_param%neb,wfn1_param%nbf,
     &                    wfn1_param%nmo,1.0d-5,g_bfbf,wfn1_wave%g_na,
     &                    .false.,eig)
      deallocate(eig)
      call ga_copy(wfn1_wave%g_na,wfn1_wave%g_nb)
c
c     == Set up correlation functions ==
c
      call wfn1_guess_corr(wfn1_param,wfn1_wave)
c
c     == Set up electron configuration ==
c
      call wfn1_configuration(wfn1_param,wfn1_wave)
c
c     == Initial wavefunction complete ==
c
      oconverged = .false.
      curit      = 0
      call wfn1_e_g(wfn1_param,wfn1_inst,wfn1_mat,wfn1_wave,wfn1_deriv,
     &              wfn1_energy,g_bfmo,g_momo1)
      do while ((.not.oconverged).and.(curit.lt.maxit))
        curit = curit + 1
        call wfn1_linesearch(wfn1_param,wfn1_inst,wfn1_mat,wfn1_wave,
     &                       wfn1_deriv,wfn1_energy,
     &                       wfn1_trial_wv,wfn1_trial_drv,wfn1_trial_en,
     &                       g_bfmo,g_momo1,g_momo2)

        error = wfn1_trace(wfn1_trial_drv,wfn1_trial_drv)
        oconverged = (error < conv_tol)
c
c       swap trial wavefunctions, derivatives and energies with the 
c       regular ones.
c
        wfn1_t_wv      = wfn1_trial_wv
        wfn1_trial_wv  = wfn1_wave
        wfn1_wave      = wfn1_t_wv
c
        wfn1_t_drv     = wfn1_trial_drv
        wfn1_trial_drv = wfn1_deriv
        wfn1_deriv     = wfn1_t_drv
c
        wfn1_energy    = wfn1_trial_en

        if (iproc.eq.0) then
          write(luout,'(i5,f20.8,f16.6)')curit,wfn1_energy%E_tot,error
        endif
      enddo
      wfn1_task_energy = oconverged
c
c     == Store results ==
c
      if (.not.rtdb_put(rtdb,'wfn1:energy',MT_DBL,1,wfn1_energy%E_tot))
     &  call errquit(pname//' failed to store energy in rtdb', 
     &               0, RTDB_ERR)
      if (.not.rtdb_put(rtdb,'wfn1:converged',MT_LOG,1,oconverged))
     &  call errquit(pname//' failed to store converged in rtdb', 
     &               0, RTDB_ERR)
c
c     == Clean up ==
c
      call wfn1_deallocate_drv(wfn1_trial_drv)
      call wfn1_deallocate_drv(wfn1_deriv)
      call wfn1_deallocate_wfn(wfn1_trial_wv)
      call wfn1_deallocate_wfn(wfn1_wave)
      call wfn1_deallocate_mat(wfn1_mat)
c
      if (.not.ga_destroy(g_momo2)) then
        call errquit(pname//" deallocate g_momo2 failed",0,GA_ERR)
      endif
      if (.not.ga_destroy(g_momo1)) then
        call errquit(pname//" deallocate g_momo1 failed",0,GA_ERR)
      endif
      if (.not.ga_destroy(g_bfmo)) then
        call errquit(pname//" deallocate g_bfmo failed",0,GA_ERR)
      endif
      if (.not.ga_destroy(g_bfbf)) then
        call errquit(pname//" deallocate g_bfbf failed",0,GA_ERR)
      endif
c
c     call int_1e_uncache_ga()
c     call schwarz_tidy
      call int_terminate
      call util_print_pop
c
      if (.not.bas_destroy(wfn1_inst%basis))
     &  call errquit(pname//' Failed to destroy basis',0,BASIS_ERR)
      if (.not.geom_destroy(wfn1_inst%geom))
     &  call errquit(pname//' Failed to destroy geometry',0,GEOM_ERR)
c
      end function wfn1_task_energy
