C> @file rt_tddft_cs_prop.F
C> @ingroup rt_tddft

C--------------------------------------------------------------------
C>  Closed-shell propagation driver.
C>
C>  Propagates closed shell density matrix using von Neumann
C>  dynamics.
C--------------------------------------------------------------------
      subroutine rt_tddft_cs_prop (params)
      implicit none

#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "util.fh"
#include "cdft.fh"
#include "matutils.fh"
#include "zora.fh"
#include "rt_tddft.fh"

      
C     == Inputs ==
      type(rt_params_t), intent(inout) :: params  !< struct containing parameters


C     == Parameters ==
      character(*), parameter     :: pname       = "rt_tddft_cs_prop: "
      

C     == External ==
      logical, external  :: rt_tddft_atpoint
      logical, external  :: rt_tddft_at_snap      
      logical, external  :: rt_tddft_cs_tdfock
      logical, external  :: rt_tddft_checkpoint


      
C     == Variables ==
      integer g_zdens_ao, g_zfock_ao !in AO basis
      integer g_zdens_mo, g_zfock_mo !in MO basis
      integer g_densre_ao, g_densim_ao
      integer g_densre_mo, g_densim_mo
      integer g_zfock_mo_n12dt  !stored previous F'(t-0.5*dt) values
      integer g_zfock_ao_init
      integer g_movecs_mo_gs
      integer g_zscr_ao, g_zscr_mo, g_zscr2_mo
      
      double precision charge, charge_elec, charge_nuc

      type(rt_energies_t) energies

      integer n, m              !aliases for ns_ao, ns_mo
      integer me
      logical ok
      integer it, rest_it, ig, iorb
      integer lmoocc, imoocc

      character*60 outstr
      logical lrestart

      double precision zjunk
      integer i, j
      

C     (time-dependent system properties)
      double precision tt       ! current time
      double precision cpusec
      double precision Exc(2)   !at most we have two Vxc


      type(rt_vector_t) dipmom 
      type(rt_quad_t) quadmom

      
C     (profiling)
      double precision elapsed, elapsed_checks, elapsed_props


C     (estimating time remaining)
      double precision est_cputime_start
      integer est_it_start


      call rt_tddft_cs_confirm (params,'rt_tddft_cs_prop.F')
      
      me = ga_nodeid()


C     == Initialization ==
      

      if (me.eq.0) call util_print_centered (luout,
     $     "Closed shell propagation", 40, .true.)
      

      
C
C     Memory allocation of time-varying properties.  Recall that all
C     quantites in AO basis are ns_ao x ns_ao, while those in MO basis are
C     ns_mo x ns_mo.
C
      n = params%ns_ao            !local alias for clarity
      m = params%ns_mo            !  ""       ""

      ok = .true.

C     (Real and imaginary parts of dens mat, used for analysis)
CXXX  [KAL]: remove or unify these???
      ok=ok.and.ga_create(mt_dbl,n,n,"Re[P]" , 0, 0, g_densre_ao)
      ok=ok.and.ga_create(mt_dbl,n,n,"Im[P]" , 0, 0, g_densim_ao)
      ok=ok.and.ga_create(mt_dbl,m,m,"Re[P']", 0, 0, g_densre_mo)
      ok=ok.and.ga_create(mt_dbl,m,m,"Im[P']", 0, 0, g_densim_mo)

C     (Fock and dens mats in AO and MO basis)
      ok=ok.and.ga_create(mt_dcpl,n,n,"P" , 0,0, g_zdens_ao)
      ok=ok.and.ga_create(mt_dcpl,m,m,"P'", 0,0, g_zdens_mo)
      ok=ok.and.ga_create(mt_dcpl,n,n,"F" , 0,0, g_zfock_ao)
      ok=ok.and.ga_create(mt_dcpl,m,m,"F'", 0,0, g_zfock_mo)

      ok=ok.and.ga_create(mt_dcpl,n,n,"AO scr", 0,0, g_zscr_ao)
      ok=ok.and.ga_create(mt_dcpl,m,m,"MO scr", 0,0, g_zscr_mo)
      ok=ok.and.ga_create(mt_dcpl,m,m,"MO scr2", 0,0, g_zscr2_mo)


C     (ground state movecs, MO basis, pure real)
      ok=ok.and.ga_create(mt_dbl,m,m,"movecs_mo_gs", 0,0,g_movecs_mo_gs)


C     (Fock mats at previous times (for Magnus) in MO basis)
      ok=ok.and.ga_create(mt_dcpl,m,m,"F'(t-1/2*dt)",
     $     0,0, g_zfock_mo_n12dt)

      
C     (static Fock matrix) !XXXXX
      if (params%static) then
         ok=ok.and.ga_create(mt_dcpl,n,n,"F(t=0)",
     $        0,0, g_zfock_ao_init)
      endif



      if (.not.ok) call errquit(pname//"GA allocation failed", 0,GA_ERR)
      
      if (.not.ma_push_get(MT_Dbl,m,'moocc',lmoocc,imoocc))
     &     call errquit(pname//'cannot allocate moocc',0, MA_ERR)

      call ga_zero (g_zdens_ao)
      call ga_zero (g_zdens_mo)
      call ga_zero (g_zfock_ao)
      call ga_zero (g_zfock_mo)


      if (do_zora) call rt_tddft_init_zora_scalar (params)


C
C     Check that canonical orthogonalization is set up correctly.
C
      call canorg_check (params, params%g_s)
      

C
C     Get starting ("ground state") movecs and density matrix in AO
C     basis (note 1 for closed shell). Then convert to dens mat to MO
C     basis.
C
C      call rt_tddft_movecs_zdens (params, 1, g_movecs_ao_gs, g_zdens_ao)
      call rt_tddft_movecs_zdens (params, 1, g_zdens_ao)
      call canorg_trans(params,"P","AO->MO", g_zdens_ao, g_zdens_mo)



C
C     Build ground state Fock matrix from ground state AO SCF density
C     and convert to MO basis.  Check that the imaginary part of F MO is
C     pure zero (as it comes from the SCF) and store ground state
C     energy.
C
C     We also store the Fock matrix if we are doing a "static"
C     calculation.
C
      params%skip_fock = .false.

      call zfock_cs_build (params, g_zdens_ao, energies, g_zfock_ao)
      call canorg_trans (params,"F","AO->MO", g_zfock_ao, g_zfock_mo)

      if ( (params%mocap_active).or.(params%cap_active) ) then
         call rt_tddft_print_notice (
     $        "MO and/or spatial CAP active: "//
     $        "Disabled initial MO Fock check")
      else
         if (.not. mat_is_purereal (g_zfock_mo, params%tol_zero))
     $        call errquit (pname//
     $        "initial Fock matrix in MO basis is not pure real", 0, 0)
      endif


      call rt_tddft_calc_Etot (params, energies)
      energies%gs = energies%tot


      if (params%static) then
         call ga_copy (g_zfock_ao, g_zfock_ao_init)
         params%g_zfock_ao_init(1) = g_zfock_ao_init
         params%skip_fock = .true.
      endif

      
C
C     Check that this Fock matrix is consistent with SCF Fock matrix.
C
      if (params%matrix_checks) then
         call zmat_compare_and_print (params%ns_ao, g_zfock_ao,
     $        "fock_cs_gs", params%tol_zero,
     $        "GS and RT-TDDFT Fock matrices")

         call zmat_compare_and_print (params%ns_ao, g_zdens_ao,
     $        "dens_cs_gs", params%tol_zero,
     $        "GS and RT-TDDFT density matrices")
      endif

C
C     If we need the eigenvectors of the ground state Fock matrix (for
C     now only when print MO occupations or dipcontribs), diagonalize in
C     MO basis. Again, pure real.  If we are using a CAP, F will be
C     complex but we give the user a warning and continue anyways.
C
      if (params%print_moocc .or. params%print_dipcontribs) then

         if ( (params%mocap_active).or.(params%cap_active) ) then
            call rt_tddft_print_notice (
     $           "MO and/or spatial CAP active: "//
     $           "Disabled initial MO Fock check "//
     $           "for moocc and/or dipcontribs")
         else

            if (.not. mat_is_purereal (g_zfock_mo, params%tol_zero))
     $           call errquit (pname//
     $           "initial Fock matrix in MO basis is not pure real",0,0)
            
         endif
         
         call rt_tddft_movecs_gs (params, g_zfock_mo, g_movecs_mo_gs)
      endif


C     ======================
C     == Time propagation ==
C     ======================

C
C     If user specified "load restart" then look for restart
C     file by name runname.rt_restart.
C     
      if (.not.rtdb_get (params%rtdb, "rt_tddft:restart",
     $     mt_log, 1, lrestart))
     $     lrestart = .false.

      if (lrestart) then
         call rt_tddft_restart_load (params, it, tt, 1, g_zdens_ao) ! 1 for rks
         call canorg_trans (params, "P", "AO->MO",
     $        g_zdens_ao, g_zdens_mo)
         call rt_tddft_print_msg (params, tt,
     $        "Propagation successfully restarted")
      else
         it = 1
         tt = params%tmin
         call rt_tddft_print_msg(params, tt, "Propagation started")
      endif

      
C
C     Build an initial TD Fock (includes applied field, unlike F
C     computed above), and clone it into F'(t-0.5*dt) for Magnus.
C
      if (.not. rt_tddft_cs_tdfock (params, tt, g_zdens_ao,
     $     energies, g_zfock_ao))
     $     call errquit (pname//"tdfock failed", 0, 0)
      
      call canorg_trans (params,"F","AO->MO", g_zfock_ao, g_zfock_mo)
      call ga_copy (g_zfock_mo, g_zfock_mo_n12dt)



C     (store start times and index for estimaing time remaining)
      est_cputime_start = util_cpusec ()
      est_it_start = it

      
      do while (it .lt. params%nt) !start time loop

         if (params%prof) call prof_start (elapsed)

C     
C     Skip time propagation if asked for.
C     
         if (params%noprop) then
            call rt_tddft_print_notice ("Propagation skipped")
            exit
         endif


         

C
C     Convert density matrix to AO basis and extract real and im parts.
C         
CXXX  [KAL]: this is perhaps wasteful;
C     - if we are not printing or doing checks there should be no need to do this here.
C
         call canorg_trans(params,"P","MO->AO",g_zdens_mo,g_zdens_ao)
         call convert_z2d (g_zdens_ao, g_densre_ao, g_densim_ao)
         call convert_z2d (g_zdens_mo, g_densre_mo, g_densim_mo)

         
C
C     == Occasionally do some checks. == 
C
C     For the check we build new Fock matrix in AO basis and then
C     convert to MO basis (wasteful but we dont do checks very often).
C     

         if (rt_tddft_atpoint(it, params%nt, params%nchecks)) then
            
            if (params%prof) call prof_start (elapsed_checks)
            
            if (.not. rt_tddft_cs_tdfock (params, tt, g_zdens_ao,
     $           energies, g_zfock_ao))
     $           call errquit (pname//"tdfock failed", 0, 0)

            call ga_zero (g_zfock_mo)
            call canorg_trans(params,"F","AO->MO",g_zfock_ao,g_zfock_mo)  !note should be same as g_zfock_mo already stored (from Magnus)


            if (.not. rt_tddft_checkpoint (params, 1, g_zfock_ao,
     $           g_zdens_ao, g_zfock_mo, g_zdens_mo)) then
               call rt_tddft_print_msg (params, tt,
     $              "WARNING: CHECKS FAILED")
            else
               call rt_tddft_print_msg (params, tt, "Checks passed")
            endif
            
            if (params%prof) call prof_end (elapsed_checks, "Checks")
         endif

         
C     (prints estimated time remaining)
         if (rt_tddft_atpoint(it, params%nt, params%nestims)) then
            if (it .ne. est_it_start)
     $           call rt_tddft_estimate_rem (params, est_cputime_start,
     $           est_it_start, it, tt)
         endif
            


C
C     == Occasionally dump real part of AO density matrix. ==
C
C         if (rt_tddft_atpoint(it, params%nt, params%nsnapshots)) then
         if (rt_tddft_at_snap (params, tt)) then
            call rt_tddft_snapshot_save (params, it, tt, g_densre_ao)
            call rt_tddft_print_msg (params, tt,
     $           "Saved density matrix snapshot")
         endif


C
C     == Occasionally calc and print time-dependent system properties. ==
C     
         if (rt_tddft_atpoint(it, params%nt, params%nprints)) then
            if (params%prof) call prof_start (elapsed_props)

C     (external applied field, also calcd in tdfock)
            if (params%print_field)
     $           call rt_tddft_print_fields (params, tt)
            

C     (compute MO occupations if we are printing them) 
C     PREVIOUS USED WHEN DOING ZORA
C            if (params%print_moocc .or. do_zora) then
            if (params%print_moocc) then
               call rt_tddft_moocc_calc (params, g_densre_mo,
     $              g_movecs_mo_gs, dbl_mb(imoocc))
            endif


C            call ga_print (g_densre_mo)


C     (cpu time calculation and output)
            if (params%print_cputime) then
               cpusec = util_cpusec ()
               call rt_tddft_print_val (params, tt, cpusec, "CPU time")
            endif


C     (system energies, parts already computed in fock build)
            if (params%print_energy) then
               call rt_tddft_calc_Etot (params, energies)
               call rt_tddft_print_energies(params, tt, energies)
            endif


C     (total charge on each geom)
            if (params%print_charge) call rt_tddft_print_charge (params,
     $           tt, g_densre_ao, "")


C     (molecular orbital occupations, calc in MO basis)
            if (params%print_moocc) then
               call rt_tddft_moocc_print (params,
     $              tt, dbl_mb(imoocc), "# MO Occupations")

C               call rt_tddft_complex_evals (params, g_zfock_mo,
C     $              "Fock MO")
               
            endif

C     (Transition Dipole Contributions [ASB])
            if (params%print_dipcontribs) then
               call rt_tddft_calc_dip_contribs (params,
     $              g_movecs_mo_gs, g_zdens_ao, tt, "CLOSED")
            endif


CXXX  [KAL]: FIX
C     (now, scale real part of AO dens mat if doing ZORA--this will
C     affect all subsequent properties for this time step)
            
C            if (do_zora)
C     $           call rt_tddft_zora_scaled_densre (params,
C     $           g_movecs_ao_gs, dbl_mb(imoocc), g_densre_ao)


C     (dipole moment on each geom, calc in AO basis)
            if (params%print_dipole) then
               do ig = 1, params%ngeoms
                  if (params%geom(ig)%measure) then
                     
                     outstr = "Dipole moment [" //
     $                    trim (params%geom(ig)%name) // "]"

                     call rt_tddft_calc_dip_geom (params,
     $                    g_densre_ao, ig, dipmom)
                     call rt_tddft_print_vector (params, tt,
     $                    dipmom, trim(outstr))
                     
                  endif
               enddo
               
            endif

C     (quadrupole moment, calc in AO basis)
            if (params%print_quad) then
               call rt_tddft_calc_quad (params, g_densre_ao, quadmom)
               call rt_tddft_print_quad (params, tt,
     $              quadmom, "Quadrupole moment")
            endif

            if (params%prof) call prof_end (elapsed_props, "Properties")
         endif  !end print point

         
C
C     Dump restart information to disk.
C
         if (rt_tddft_atpoint(it, params%nt, params%nrestarts))
     $        call rt_tddft_restart_save (params, it, tt, 1, g_zdens_ao)


C     
C     Step MO density matrix forward dt in time.
C     
         if (params%prop_method .eq. 1) then
            call prop_euler (params, tt, rt_tddft_cs_tdfock,
     $           1, g_zdens_mo, energies)

         elseif (params%prop_method .eq. 2) then
            call prop_rk4 (params, tt, rt_tddft_cs_tdfock,
     $           1, g_zdens_mo, energies)
            
         elseif (params%prop_method .eq. 3) then
            call prop_magnus_ord2_interpol (params, tt,
     $           rt_tddft_cs_tdfock, 1,
     $           g_zfock_mo_n12dt, g_zfock_mo, energies, g_zdens_mo)

         else
            call errquit (pname//"bad prop_method", 0, 0)
         endif


C
C     We are now at t + dt.
C
         it = it + 1
         tt = tt + params%dt
C         tt = params%tmin + (it-1)*params%dt  


         
         if (params%prof) call prof_end (elapsed, "Complete time step")

      enddo                     !end time loop


C
C     I suppose technically we should print the time-dependent system
C     properties here, as we have yet to that for the current (advanced)
C     time.  We're only losing one value (the last time) so it's not a
C     big deal, though....
C

      
C     (NOTE: uses it, tt values *outside* loop but not important)
      call rt_tddft_print_msg (params, tt, "Propagation finished")


C     == Cleanup GAs for time-dependent quantities and params ==
      ok = .true.
      ok = ok .and. ga_destroy(g_zdens_ao)
      ok = ok .and. ga_destroy(g_zfock_ao)
      ok = ok .and. ga_destroy(g_zdens_mo)
      ok = ok .and. ga_destroy(g_zfock_mo)
      ok = ok .and. ga_destroy(g_zfock_mo_n12dt)
      ok = ok .and. ga_destroy(g_densre_ao)
      ok = ok .and. ga_destroy(g_densim_ao)
      ok = ok .and. ga_destroy(g_densre_mo)
      ok = ok .and. ga_destroy(g_densim_mo)
      ok = ok .and. ga_destroy(g_movecs_mo_gs)
      ok = ok .and. ga_destroy(g_zscr_ao)
      ok = ok .and. ga_destroy(g_zscr_mo)
      ok = ok .and. ga_destroy(g_zscr2_mo)


C     (static Fock matrix)
      if (params%static) then
         ok=ok.and.ga_destroy (g_zfock_ao_init)
      endif


      if (.not.ok) call errquit(pname//"failed to destroy GAs", 0,0)

      if (.not.ma_chop_stack(lmoocc))
     &   call errquit(pname//'cannot chop stack',0, MA_ERR)


C
C     Clean up charge density fitting if needed.
C
C      if (cdfit) call zfock_cs_coulcdfit_clean (params)


      end subroutine
c $Id$
