c
c     ------------------------------------------------------------------
c     Main routine for calculating the gradients of the TDDFT excitation
c     energies with respect to geometric perturbations.
c     ------------------------------------------------------------------
c
      logical function tddft_gradient(rtdb)
c
c     Notation:  Xia (i=occupied,a=virtual) 
c     Occupied index runs faster than the virtual index.
c
c     The notation closely follows the Furche, Ahlrichs paper
c     "Adiabatic time-dependent density functional methods for excited state properties"
c     F. Furche, R. Ahlrichs, J. Chem. Phys. 117, 7433 (2002)
c
c     Please cite the following paper for the NWChem implementation:
c
c     Simulating One-Photon Absorption and Resonance Raman Scattering Spectra
c     Using Analytical Excited State Energy Gradients within
c     Time-Dependent Density Functional Theory
c
c     D.W. Silverstein, N. Govind, H.J.J. van Dam, L. Jensen
c     J. Chem. Theory Comput., 2013, 9 (12), pp 5490–5503 DOI: 10.1021/ct4007772
c
      implicit none
c
#include "mafdecls.fh"
#include "global.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "xc.fh"
#include "rtdb.fh"
#include "cosmo.fh"
#include "tddft_grad_util.fh"
c
      character*255 title ! the job title
c
      integer rtdb        ! the runtime database handle
      integer ipol        ! =1 (restricted), =2 (unrestricted)
      logical tda         ! True if Tamm-Dancoff approximation
      logical oskel       ! Use of symmetry in Fock build
      logical osinglet    ! Compute singlet excited states?
      logical otriplet    ! Compute triplet excited states?
      logical lhascd
      integer nsroots     ! the number of singlet roots
      integer ihdl_sroots ! the singlet roots handle
      integer iptr_sroots ! the singlet roots list
      integer ntroots     ! the number of triplet roots
      integer ihdl_troots ! the triplet roots handle
      integer iptr_troots ! the triplet roots list
      integer nroots      ! the number of roots
      integer ihdl_roots  ! the roots handle
      integer iptr_roots  ! the roots list
      integer ntotroots   ! the total number of roots
      integer ihdl_bfao   ! the AO basis set handle
      integer ihdl_geom   ! the geometry handle
c
      integer nao         ! the number of AOs
      integer nmo(2)      ! the number of MOs
      integer nfc(2)      ! the number of frozen core orbitals
      integer nfv(2)      ! the number of frozen virtual orbitals
      integer nocc(2)     ! the number of occupied orbitals
      integer naoc(2)     ! the number of active occupied orbitals
      integer nav(2)      ! the number of active virtual orbitals
      integer nov(2)      ! nov(i) = naocc(i)*nav(i) 
      integer nat         ! the number of atoms
c
      integer ihdl_mo_e   ! mo energy memory handle
      integer iptr_mo_e   ! mo energy memory index
      integer ihdl_root_e ! excitation energy memory handle
      integer iptr_root_e ! excitation energy memory index
c
      double precision tol2e ! 2-electron integral tolerance
c
      double precision solve_thresh ! convergence threshold for solver
c
      integer solve_maxit  ! the maximum number of iterations in solver
      integer solve_maxsub ! the maximum number expansion vectors
                           ! per root
c
      character*255 fn_solve   ! filename for linear system solutions
      character*6   solve_when ! when to store solutions
c
      character*255 fn_civecs         ! filename for unrestricted vectors
      character*255 fn_civecs_singlet ! filename for singlet state vectors
      character*255 fn_civecs_triplet ! filename for triplet state vectors
      character*255 fn_movecs         ! filename for MO vectors
c
      character*40 xname(numfunc)       ! Names of X functionals
      character*40 cname(numfunc)       ! Names of C functionals
      character*40 xcname(numfunc)      ! Names of XC functionals
      double precision xfac(numfunc)    ! Prefactors for X functionals
      double precision cfac(numfunc)    ! Prefactors for C functionals
      logical lxfac(numfunc)            ! True if local X functional
      logical nlxfac(numfunc)           ! True if gc X functional
      logical lcfac(numfunc)            ! True if local C functional
      logical nlcfac(numfunc)           ! True if gc C functional
      logical xccomb(numfunc)           ! True if composite XC functional
c
      logical lhf                       ! True if HF (no DFT)
      logical ldft                      ! True if pure DFT (no HF exch)
      logical lhybrid                   ! True if hybrid DFT and HF
      logical lhasdft                   ! True if it has DFT
      logical lhashf                    ! True if it has HF
c
c     Global array for various quantities 
c     g_something(1) = alpha component (also in RHF case)
c     g_something(2) = beta  component
c
      integer g_mo(2)      ! global arrays to the MOs
      integer g_xpy(2)     ! global arrays to hold X+Y
      integer g_xmy(2)     ! global arrays to hold X-Y
      integer g_tp(2)      ! global arrays to hold T, and later P
      integer g_r(2)       ! global arrays to hold right-hand-sides
      integer g_z(2)       ! global arrays to hold solutions
      integer g_w(2)       ! global arrays to hold energy weighted density matrices
      integer g_g          ! global array to hold the gradients
c
      double precision kfac ! Prefactor for HF exchange
c
      double precision tddft_grad_wall ! wall clock time
      double precision tddft_grad_cpu  ! cpu time
c
      logical tddft_energy
      external tddft_energy
      logical oroot
      logical odebug
c
      character*32 pname
c
      parameter (oskel=.false.)
      logical task_gradient
      external task_gradient
      logical task_energy
      external task_energy
c
c MN solvation models -->
c
      integer do_cosmo_vem, istep_cosmo_vem
      integer g_vem(3)
c
      if (.not. rtdb_get
     $   (rtdb,'cosmo:do_cosmo_vem',mt_int,1,do_cosmo_vem))
     $    do_cosmo_vem = 0
c
c <-- MN solvation models
c
      pname = "tddft_gradient: "
c
      oroot = ga_nodeid().eq.0
c
c     Debug flag
      odebug = .false.
      if (rtdb_get(rtdb,'tddft_grad:debug', mt_log, 1, odebug))
     &     odebug = .true.
c
c     TDDFT gradient has not been computed successfully yet.
      tddft_gradient = .false.
c
c     Make sure we have the excited states first
      call ga_sync
      if (.not.tddft_energy(rtdb)) return
      call ga_sync
c
c     Get the start times
      tddft_grad_wall = -util_wallsec()
      tddft_grad_cpu  = -util_cpusec()
c
c     Init pstat
      call dft_pstat_init(rtdb)
c      
      call util_print_push
      call util_print_rtdb_load(rtdb,'tddft')
c
c     Print header
      if (oroot) call tddft_grad_print_header(luout)
c
c     Parameter initialization
      call tddft_grad_init(rtdb,ipol,tda,nat,
     +  nao,nocc,nmo,nfc,nfv,naoc,nav,nov,fn_movecs,numfunc,
     +  xname,cname,xcname,xfac,cfac,lxfac,nlxfac,lcfac,
     +  nlcfac,xccomb,lhf,ldft,lhybrid,lhasdft,lhashf,
     +  nroots, iptr_roots, ihdl_roots, fn_civecs,
     +  nsroots,iptr_sroots,ihdl_sroots,fn_civecs_singlet,
     +  ntroots,iptr_troots,ihdl_troots,fn_civecs_triplet,
     +  osinglet,otriplet, ntotroots,ihdl_bfao,ihdl_geom,tol2e,
     +  solve_thresh,solve_maxit,solve_maxsub,
     +  fn_solve,solve_when,kfac,title,lhascd)
c
c     Cosmo initialization, if necessary
      cosmo_on = .false.
      if (rtdb_get(rtdb,'slv:cosmo',mt_log,1,cosmo_on)) then
       if(cosmo_on) then
c
c MN solvation models -->
c
c Tell hnd_coschg in src/property/cosmo.F that we have an excited-state
c density (istep_cosmo_vem = 1) when it needs to know it (i.e., only when
c a VEM calculation is called). After calling hnd_coschg via cosmo_initialize 
c istep_cosmo_vem changes either to 2 (further iteration) or 
c to 3 (exit the calculation upon convergence).
c
         istep_cosmo_vem = 0
         if (do_cosmo_vem.ne.0) istep_cosmo_vem = 1
         if (.not. rtdb_put
     $   (rtdb,'cosmo:istep_cosmo_vem',mt_int,1,istep_cosmo_vem))
     $   call errquit(
     $ 'tddft_gradient: cannot put istep_cosmo_vem in rtdb',
     $   0,rtdb_err)
c
c <-- MN solvation models
c
        call cosmo_initialize(rtdb,ihdl_geom,ihdl_bfao,.false.)
        cosmo_phase = 2
        if (.not. rtdb_put(rtdb,'cosmo_phase',mt_int,1,cosmo_phase))
     $   call errquit('cosmo_input: put cosmo_phase failed',
     $           911, RTDB_ERR)
       endif
      endif
c
c     Integral initialization
      call tddft_grad_init_ints(rtdb,ihdl_bfao,ihdl_geom,tol2e)
c
c     Print calculation information
      call tddft_grad_print_info(luout,nroots,int_mb(iptr_roots),
     +     nsroots,int_mb(iptr_sroots),ntroots,int_mb(iptr_troots),
     +     title)
c
c     Print start times
      if (oroot) then
        if (util_print('tddft grad time',print_default)) then
          write(luout,*)
          write(luout,"(' Start at time  cpu:',f11.1,'s  wall:',f11.1,
     +                  's')")-tddft_grad_cpu,-tddft_grad_wall
          write(luout,*)
        endif
      end if
c
c     Allocate orbital energies
      if (oroot.and.odebug) write(luout,*) "call tddft_grad_create_e"
      call tddft_grad_create_e(ipol,nao,ihdl_mo_e,iptr_mo_e) 
      if (oroot.and.odebug) write(luout,*) "done tddft_grad_create_e"
c
c     Allocate MO coefficients.
      if (oroot.and.odebug) write(luout,*) "call tddft_grad_create_mo"
      call tddft_grad_create_mo(ipol,ihdl_geom,ihdl_bfao,g_mo)
      if (oroot.and.odebug) write(luout,*) "done tddft_grad_create_mo"
c
c     Allocate (X+Y) & (X-Y)
      if (oroot.and.odebug) write(luout,*) "call tddft_grad_create_xy"
      call tddft_grad_create_xy(tda,ipol,naoc,nav,ntotroots,g_xpy,g_xmy)
      if (oroot.and.odebug) write(luout,*) "done tddft_grad_create_xy"
c
c     Allocate T, T+(X+Y), or T+V(=P)
      if (oroot.and.odebug) write(luout,*) "call tddft_grad_create_tp"
      call tddft_grad_create_tp(ipol,naoc,nav,ntotroots,g_tp)
      if (oroot.and.odebug) write(luout,*) "done tddft_grad_create_tp"
c
c     Allocate excitation energies
      if (oroot.and.odebug) 
     &   write(luout,*)"call tddft_grad_create_root_e"
      call tddft_grad_create_root_e(ntotroots,ihdl_root_e,iptr_root_e)
      if (oroot.and.odebug) 
     &   write(luout,*)"done tddft_grad_create_root_e"
c
c     Load the MO vectors from the file specified by fn_movecs
      if (oroot.and.odebug) write(luout,*) "call tddft_grad_load_mo"
      call tddft_grad_load_mo(ipol,nao,fn_movecs,dbl_mb(iptr_mo_e),g_mo)
      if (oroot.and.odebug) write(luout,*) "done tddft_grad_load_mo"
c
c     Load CI vectors from TDDFT calculation: (X+Y) and (X-Y)
      if (oroot.and.odebug) write(luout,*) "call tddft_grad_load_xy"
      call tddft_grad_load_xy(tda,ipol,naoc,nav,nov,nroots,
     +     int_mb(iptr_roots),fn_civecs,nsroots,int_mb(iptr_sroots),
     +     fn_civecs_singlet,ntroots,int_mb(iptr_troots),
     +     fn_civecs_triplet,g_xpy,g_xmy,dbl_mb(iptr_root_e))
      if (oroot.and.odebug) write(luout,*) "done tddft_grad_load_xy"
c
c     Compute the TDDFT energy for each root using the solution vectors
      if (oroot.and.odebug)
     &    write(luout,*)"call tddft_grad_compute_energy"
      call tddft_grad_compute_energy(rtdb,ihdl_geom,ihdl_bfao,tol2e,tda,
     +   ipol,nao,nfc,naoc,nocc,nav,nfv,ntotroots,
     +   dbl_mb(iptr_mo_e),g_mo,g_xpy,g_xmy,kfac,lhashf,otriplet,lhascd)
      if (oroot.and.odebug)
     &    write(luout,*)"done tddft_grad_compute_energy"
c
c     Compute the unrelaxed difference density T
      if (oroot.and.odebug) write(luout,*) "call tddft_grad_compute_t"
      call tddft_grad_compute_t(tda,ipol,naoc,nav,ntotroots,g_xpy,g_xmy,
     +     g_tp)
      if (oroot.and.odebug) write(luout,*) "done tddft_grad_compute_t"
c
c     Allocate right hand side (R: g_r) and solution arrays (Z: g_z)
      if (oroot.and.odebug) write(luout,*) "call tddft_grad_create_r"
      call tddft_grad_create_r(ipol,naoc,nav,ntotroots,g_r)
      call tddft_grad_create_r(ipol,naoc,nav,ntotroots,g_z)
      if (oroot.and.odebug) write(luout,*) "done tddft_grad_create_r"
c
c     Assemble right hand side: R
      if (oroot.and.odebug) write(luout,*) "call tddft_grad_compute_r"
      call tddft_grad_compute_r(rtdb,ihdl_geom,ihdl_bfao,tol2e,tda,ipol,
     +     nao,nfc,naoc,nocc,nav,nfv,ntotroots,g_mo,g_xpy,g_xmy,g_tp,
     +     g_r,kfac,lhashf,otriplet,lhascd)
      if (oroot.and.odebug) write(luout,*) "done tddft_grad_compute_r"
c
c     Solve the Z equation: (A+B)Z = -R
c     Note: g_z has the solutions
      if (oroot.and.odebug) write(luout,*) "call tddft_grad_solve"
      call tddft_grad_solve(rtdb,ihdl_geom,ihdl_bfao,ipol,ntotroots,
     +     nao,nocc,naoc,nav,nfc,nfv,g_mo,iptr_mo_e,tda,oskel,kfac,
     +     lhashf,otriplet,solve_thresh,solve_maxsub,solve_maxit,
     +     fn_solve,solve_when,tol2e,g_r,g_z,lhascd)
      if (oroot.and.odebug) write(luout,*) "done tddft_grad_solve"
c
c     Deallocate right hand side: R
      if (oroot.and.odebug) write(luout,*) "call tddft_grad_destroy_r"
      call tddft_grad_destroy_r(ipol,g_r)
      if (oroot.and.odebug) write(luout,*) "done tddft_grad_destroy_r"
c
c     Compute the relaxed density matrix: P=T+Z
      if (oroot.and.odebug) write(luout,*) "call tddft_grad_compute_p"
      call tddft_grad_compute_p(ipol,naoc,nav,ntotroots,g_tp,g_z)
      if (oroot.and.odebug) write(luout,*) "done tddft_grad_compute_p"
c
c     Create Lagrange multipliers: W
      if (oroot.and.odebug) write(luout,*) "call tddft_grad_create_tp"
      call tddft_grad_create_tp(ipol,naoc,nav,ntotroots,g_w)
      if (oroot.and.odebug) write(luout,*) "done tddft_grad_create_tp"
c
c     Compute Lagrange multipliers: W
c     W is the energy-weighted difference density matrix
      if (oroot.and.odebug) write(luout,*) "call tddft_grad_compute_w"
      call tddft_grad_compute_w(rtdb,ihdl_geom,ihdl_bfao,tol2e,tda,ipol,
     +     ntotroots,nfc,naoc,nocc,nav,nfv,nao,g_mo,g_tp,g_z,
     +     g_xpy,g_xmy,dbl_mb(iptr_mo_e),dbl_mb(iptr_root_e),g_w,
     +     kfac,lhashf,otriplet,lhascd)
      if (oroot.and.odebug) write(luout,*) "done tddft_grad_compute_w"
c
c     Finalize integrals
      if (oroot.and.odebug)write(luout,*)"call tddft_grad_finalize_ints"
      call tddft_grad_finalize_ints(rtdb)
      if (oroot.and.odebug)write(luout,*)"done tddft_grad_finalize_ints"
c
c     Create gradients array
      if (oroot.and.odebug) write(luout,*) "call tddft_grad_create_g"
      call tddft_grad_create_g(nat,ntotroots,g_g)
      if (oroot.and.odebug) write(luout,*) "done tddft_grad_create_g"
c
c     Compute the gradients
      if (oroot.and.odebug) write(luout,*) "call tddft_grad_compute_g"
      call tddft_grad_compute_g(rtdb,ihdl_geom,ihdl_bfao,tda,nat,nao,
     +  nfc,naoc,nocc,nav,nfv,ipol,ntotroots,oskel,tol2e,g_mo,
     +  g_tp,g_w,g_xpy,g_xmy,g_g,kfac,lhashf,otriplet)
      if (oroot.and.odebug) write(luout,*) "done tddft_grad_compute_g"
c
c     Deallocate arrays
      if (oroot.and.odebug) write(luout,*) "call deallocate arrays"
      call tddft_grad_destroy_tp(ipol,g_w)
      call tddft_grad_destroy_r(ipol,g_z)
      call tddft_grad_destroy_tp(ipol,g_tp)
      call tddft_grad_destroy_xy(tda,ipol,g_xpy,g_xmy)
      call tddft_grad_destroy_mo(ipol,g_mo)
      call tddft_grad_destroy_root_e(ihdl_root_e)
      call tddft_grad_destroy_e(ihdl_mo_e)
      if (oroot.and.odebug) write(luout,*) "done deallocate arrays"
c MN solvation models -->
c
c We don't need to print out the excited-state gradients in the case of VEM
c
      if (do_cosmo_vem.ne.0) goto 100
c
c <-- MN solvation models
c
c     Print gradients of the excited states of interest
      if (oroot.and.odebug) 
     &     write(luout,*) "call tddft_grad_print_result"
      if (ipol.eq.2) then
        call tddft_grad_print_result(ntotroots,int_mb(iptr_roots),
     +                             nat,ihdl_geom,g_g,1.0d0)
      else
        if (osinglet) then
          call tddft_grad_print_result(ntotroots,int_mb(iptr_sroots),
     +                               nat,ihdl_geom,g_g,1.0d0)
        else if (otriplet) then
          call tddft_grad_print_result(ntotroots,int_mb(iptr_troots),
     +                               nat,ihdl_geom,g_g,1.0d0)
        endif
      endif
      if (oroot.and.odebug) 
     &      write(luout,*) "done tddft_grad_print_result"
c
c MN solvation models -->
c
 100  continue
c
c <-- MN solvation models
c
c     Store gradients
      if (oroot.and.odebug)
     &      write(luout,*)"call tddft_grad_store_gradient"
      call tddft_grad_store_gradient(rtdb,nat,ntotroots,
     +                         int_mb(iptr_sroots),g_g)
      if (oroot.and.odebug)
     &      write(luout,*)"done tddft_grad_store_gradient"
c
c     Deallocate gradients
      if (oroot.and.odebug)
     &      write(luout,*) "call tddft_grad_destroy_g"
      call tddft_grad_destroy_g(g_g)
      if (oroot.and.odebug)
     &      write(luout,*) "done tddft_grad_destroy_g"
c
c     Finalize gradients
      if (oroot.and.odebug) 
     &      write(luout,*) "call tddft_grad_finalize"
      call tddft_grad_finalize(nsroots,ihdl_sroots,ntroots,ihdl_troots,
     +       nroots, ihdl_roots, ihdl_geom, ihdl_bfao)
      if (oroot.and.odebug) 
     &      write(luout,*) "done tddft_grad_finalize"
c
c     Clean up cosmo stuff, if necessary
      if (oroot.and.odebug.and.cosmo_on)
     &      write(luout,*) "call cosmo_tidy"
      if (cosmo_on) then
         call cosmo_tidy(rtdb)
         cosmo_on = .false.
      endif
      if (oroot.and.odebug.and.cosmo_on)
     &      write(luout,*) "done cosmo_tidy"
c
c     Clean up
      call ga_sync()
      call grid_cleanup(.true.)
      tddft_gradient = .true.
c
c     Print end times and wrap up
      tddft_grad_cpu  = tddft_grad_cpu  + util_cpusec()
      tddft_grad_wall = tddft_grad_wall + util_wallsec()
      if (oroot) then
        if (util_print('tddft grad time',print_default)) then
          write(LuOut,*)
          write(LuOut,"(' TDDFT Gradient time  cpu:',f11.1,'s  wall:',
     +                  f11.1,'s')")tddft_grad_cpu,tddft_grad_wall
        endif
      endif
c
c MN solvation models -->
c
c cosmo-vem: Do one more TDDFT iteration. The whole SCF--TDDFT energy--TDDFT density
c calculation will be repeated as many times as needed to achieve 
c convergence on the VEM energy as determined by hnd_coschg in src/property/cosmo.F.
c When such a convergence is achieved istep_cosmo_vem is set to 3 and the program exits the loop.  
c
      cosmo_on = .false.
      if (rtdb_get(rtdb,'slv:cosmo',mt_log,1,cosmo_on)) then
       if(cosmo_on) then
        istep_cosmo_vem = 0
        if (.not. rtdb_get
     $(rtdb,'cosmo:istep_cosmo_vem',mt_int,1,istep_cosmo_vem))
     $  call errquit('tddft_grad: cannot get istep_cosmo_vem',
     $ 0,rtdb_err)
        if (istep_cosmo_vem.eq.2.and.do_cosmo_vem.ne.0) then
         if (.not.task_gradient(rtdb)) call errquit
     $('error: task_gradient(rtdb) failed for cosmo-vem',
     & 0,geom_err)
        endif
c
c cosmo-vem: Do an additional DFT energy calculation if emission is requested overwise do nothing. 
c
        if (istep_cosmo_vem.eq.3.and.do_cosmo_vem.eq.2) then
         if (.not. rtdb_cput(rtdb, 'task:theory', 1, 'dft'))
     $ call errquit('error: failed changing theory for for cosmo-vem',
     1 0,rtdb_err)
         if (.not. rtdb_put(rtdb,'dft:converged', mt_log, 1, .false.))
     $ call errquit("task: could not write to rtdb",
     $ 0,rtdb_err)
         if (.not.task_energy(rtdb)) call errquit
     $('error: task_energy(rtdb) failed for cosmo-vem',
     & 0,geom_err)
        endif
c
c cosmo-vem: Destroy global arrays
c
        if (istep_cosmo_vem.eq.3) then
         if(rtdb_get(rtdb,'cosmo:g_vem',mt_int,3,g_vem)) then
          if (.not. ga_destroy(g_vem(1))) call errquit
     $        ('hnd_coschg: could not destroy g_vem(1)', 0, ga_err)
          if (.not. ga_destroy(g_vem(2))) call errquit
     $        ('hnd_coschg: could not destroy g_vem(2)', 0, ga_err)
          if (.not. ga_destroy(g_vem(3))) call errquit
     $        ('hnd_coschg: could not destroy g_vem(3)', 0, ga_err)
          call ga_sync()
         endif
        endif
       endif
      endif
c
c <-- MN solvation models 
c
      call dft_pstat_print
      call util_print_pop
      return
      end
c $Id$
