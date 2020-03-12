C
C     rt_tddft_os_prop.F
C
C     Propagation of open shell density matrix using von Neumann
C     dynamics.
C
      subroutine rt_tddft_os_prop (params)
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
      type(rt_params_t), intent(inout) :: params  ! struct containing parameters


C     == Parameters ==
      character(*), parameter     :: pname       = "rt_tddft_os_prop: "
      integer, parameter          :: nspin = 2 !open shell

      
C     == External ==
      logical, external  :: rt_tddft_atpoint
      logical, external  :: rt_tddft_at_snap
      logical, external  :: rt_tddft_os_tdfock
      logical, external  :: rt_tddft_checkpoint

      
C     == Variables ==
      integer g_zdens_ao(2), g_zfock_ao(2) !in AO basis
      integer g_zdens_mo(2), g_zfock_mo(2) !in MO basis
      integer g_densre_ao(2), g_densim_ao(2)
      integer g_densre_mo(2), g_densim_mo(2)
      integer g_zfock_mo_n12dt(2) !stored previous F'(t-0.5*dt) values
      integer g_zfock_ao_init(2)
      integer g_densre_ao_tot
      integer g_movecs_mo_gs(2)
      integer g_zscr_mo(2)
      type(rt_energies_t) energies


      integer is ! spin counter
      integer ig ! geom counter
      integer n, m              !aliases for ns_ao, ns_mo
      integer me
      logical ok
      integer lmoocc(2), imoocc(2)

      double complex zjunk
      
      character*60 outstr
      logical lrestart
      
      integer it, rest_it


C     (time-dependent system properties)
      double precision tt       ! current time
      double precision dtrace
      double precision cpusec
      double precision s2val

      type(rt_vector_t) field !external (applied) uniform elec field
      type(rt_vector_t) dipmom_alpha, dipmom_beta, dipmom_tot
      type(rt_quad_t) qmom_alpha, qmom_beta, qmom_tot


      
C     (profiling)
      double precision elapsed, elapsed_checks, elapsed_props


C     (estimating time remaining)
      double precision est_cputime_start
      integer est_it_start


      call rt_tddft_os_confirm (params,'rt_tddft_os_prop.F')

      
      me = ga_nodeid()


C     == Initialization ==

      if (me.eq.0) call util_print_centered (luout,
     $     "Open shell propagation", 40, .true.)

      
C
C     Memory allocation of time-varying properties.  Recall that all
C     quantites in AO basis are ns_ao x ns_ao, while those in MO basis are
C     ns_mo x ns_mo.
C
      n = params%ns_ao            !local alias for clarity
      m = params%ns_mo            !  ""       ""

      ok = .true.

C     (Real and imaginary parts of dens mat, used for analysis)

      ok=ok.and.ga_create(mt_dbl,n,n,"Re[Pa]", 0, 0, g_densre_ao(1))
      ok=ok.and.ga_create(mt_dbl,n,n,"Re[Pb]", 0, 0, g_densre_ao(2))

      ok=ok.and.ga_create(mt_dbl,n,n,"Im[Pa]" , 0, 0, g_densim_ao(1))
      ok=ok.and.ga_create(mt_dbl,n,n,"Im[Pb]" , 0, 0, g_densim_ao(2))

      ok=ok.and.ga_create(mt_dbl,m,m,"Re[P'a]", 0, 0, g_densre_mo(1))
      ok=ok.and.ga_create(mt_dbl,m,m,"Re[P'b]", 0, 0, g_densre_mo(2))

      ok=ok.and.ga_create(mt_dbl,m,m,"Im[P'a]", 0, 0, g_densim_mo(1))
      ok=ok.and.ga_create(mt_dbl,m,m,"Im[P'b]", 0, 0, g_densim_mo(2))

      ok=ok.and.ga_create(mt_dbl,n,n,"Re[Ptot]", 0, 0, g_densre_ao_tot)

C     (Real and imaginary dens and Fock mats)
      ok=ok.and.ga_create(mt_dcpl,n,n,"Pa" , 0,0, g_zdens_ao(1))
      ok=ok.and.ga_create(mt_dcpl,n,n,"Pb" , 0,0, g_zdens_ao(2))

      ok=ok.and.ga_create(mt_dcpl,m,m,"P'a", 0,0, g_zdens_mo(1))
      ok=ok.and.ga_create(mt_dcpl,m,m,"P'b", 0,0, g_zdens_mo(2))
      
      ok=ok.and.ga_create(mt_dcpl,n,n,"Fa" , 0,0, g_zfock_ao(1))
      ok=ok.and.ga_create(mt_dcpl,n,n,"Fb" , 0,0, g_zfock_ao(2))

      ok=ok.and.ga_create(mt_dcpl,m,m,"F'a", 0,0, g_zfock_mo(1))
      ok=ok.and.ga_create(mt_dcpl,m,m,"F'b", 0,0, g_zfock_mo(2))

      ok=ok.and.ga_create(mt_dcpl,m,m,"MO scr 1", 0,0, g_zscr_mo(1))
      ok=ok.and.ga_create(mt_dcpl,m,m,"MO scr 2", 0,0, g_zscr_mo(2))


      
C     (ground state movecs, pure real)
      ok=ok.and.ga_create(mt_dbl,m,m,"movecs_gs a", 0,0,
     $     g_movecs_mo_gs(1))
      ok=ok.and.ga_create(mt_dbl,m,m,"movecs_gs b", 0,0,
     $     g_movecs_mo_gs(2))


C     (Fock mats at previous times (for Magnus) in MO basis)
      ok=ok.and.ga_create(mt_dcpl,m,m,"F'a(t-1/2*dt)",
     $     0,0, g_zfock_mo_n12dt(1))
      ok=ok.and.ga_create(mt_dcpl,m,m,"F'b(t-1/2*dt)",
     $     0,0, g_zfock_mo_n12dt(2))
      
C     (static Fock matrix)
      if (params%static) then
         ok=ok.and.ga_create(mt_dcpl,n,n,"F(t=0) alpha",
     $        0,0, g_zfock_ao_init(1))
         ok=ok.and.ga_create(mt_dcpl,n,n,"F(t=0) beta",
     $        0,0, g_zfock_ao_init(2))
      endif


      if (.not.ok) call errquit(pname//"GA allocation failed", 0,GA_ERR)
      
      if (.not.ma_push_get(MT_Dbl,m,'moocc alpha',lmoocc(1),imoocc(1)))
     &     call errquit(pname//'cannot allocate moocc(1)',0, MA_ERR)

      if (.not.ma_push_get(MT_Dbl,m,'moocc alpha',lmoocc(2),imoocc(2)))
     &     call errquit(pname//'cannot allocate moocc(2)',0, MA_ERR)


      do is = 1, nspin
         call ga_zero (g_zdens_ao(is))
         call ga_zero (g_zdens_mo(is))
         call ga_zero (g_zfock_ao(is))
         call ga_zero (g_zfock_mo(is))
      enddo

      
      if (do_zora) call rt_tddft_init_zora_scalar (params)
      

C
C     Check that canonical orthogonalization is set up correctly.
C
      call canorg_check (params, params%g_s)


C
C     Get starting ("ground state") density matrix in AO basis.  Then
C     convert to dens mat to MO basis.
C
      
      call rt_tddft_os_movecs_zdens (params, g_zdens_ao) !old way XXX should unify with closedshell way
C      call rt_tddft_movecs_zdens (params, 2, g_zdens_ao) !doesnt work?  I forget--check


      call canorg_trans(params,"P","AO->MO",g_zdens_ao(1),g_zdens_mo(1))
      call canorg_trans(params,"P","AO->MO",g_zdens_ao(2),g_zdens_mo(2))


C
C     Build ground state Fock matrix from ground state SCF density and
C     convert to MO basis.  Check that the imaginary part of F MO is
C     pure zero (as it comes from the SCF).  Also store ground state
C     energy.
C
C     We also store the Fock matrix if we are doing a "static"
C     calculation.
C
      params%skip_fock = .false.
      
      call zfock_os_build (params, g_zdens_ao, energies, g_zfock_ao) 
      call canorg_trans(params,"F","AO->MO",g_zfock_ao(1),g_zfock_mo(1))
      call canorg_trans(params,"F","AO->MO",g_zfock_ao(2),g_zfock_mo(2))

      if ( (params%mocap_active).or.(params%cap_active) ) then
         call rt_tddft_print_notice (
     $        "MO and/or spatial CAP active: "//
     $        "Disabled initial MO Fock check")
      else

         if (.not. mat_is_purereal (g_zfock_mo(1), params%tol_zero))
     $        call errquit (pname//
     $        "initial Fock alpha matrix in MO basis not pure real",0,0)

         if (.not. mat_is_purereal (g_zfock_mo(2), params%tol_zero))
     $        call errquit (pname//
     $        "initial Fock beta matrix in MO basis not pure real",0,0)
      endif

      
      call rt_tddft_calc_Etot (params, energies)
      energies%gs = energies%tot

      if (params%static) then
         call ga_copy (g_zfock_ao(1), g_zfock_ao_init(1))
         call ga_copy (g_zfock_ao(2), g_zfock_ao_init(2))
         params%g_zfock_ao_init = g_zfock_ao_init
         params%skip_fock = .true.
      endif


C
C     Check that this Fock matrix is consistent with SCF Fock matrix.
C
      if (params%matrix_checks) then
         call zmat_compare_and_print (params%ns_ao, g_zfock_ao(1),
     $        "fock_alpha_gs", params%tol_zero,
     $        "GS and RT-TDDFT alpha Fock matrices")
         
         call zmat_compare_and_print (params%ns_ao, g_zfock_ao(2),
     $        "fock_beta_gs", params%tol_zero,
     $        "GS and RT-TDDFT beta Fock matrices")
         
         call zmat_compare_and_print (params%ns_ao, g_zdens_ao(1),
     $        "dens_alpha_gs", params%tol_zero,
     $        "GS and RT-TDDFT alpha density matrices")
         
         call zmat_compare_and_print (params%ns_ao, g_zdens_ao(2),
     $        "dens_beta_gs", params%tol_zero,
     $        "GS and RT-TDDFT beta density matrices")
      endif


C
C     If we need the eigenvectors of the ground state Fock matrix (for
C     now only when print MO occupations or dipcontribs), diagonalize in
C     MO basis. Again, pure real.  If we are using a CAP, F will be
C     complex but we give the user a warning and continue anyways.
C
c$$$      if (params%print_dipcontribs) then
c$$$         if (.not. mat_is_purereal (g_zfock_mo(1), params%tol_zero))
c$$$     $        call errquit (pname//
c$$$     $        "initial Fock matrix in MO basis is not pure real", 0, 0)
c$$$
c$$$         if (.not. mat_is_purereal (g_zfock_mo(2), params%tol_zero))
c$$$     $        call errquit (pname//
c$$$     $        "initial Fock matrix in MO basis is not pure real", 0, 0)
c$$$
c$$$         call rt_tddft_movecs_gs (params, g_zfock_mo(1),
c$$$     $        g_movecs_mo_gs(1))
c$$$         call rt_tddft_movecs_gs (params, g_zfock_mo(2),
c$$$     $        g_movecs_mo_gs(2))
c$$$      endif

      if (params%print_moocc .or. params%print_dipcontribs) then

         if ( (params%mocap_active).or.(params%cap_active) ) then
            call rt_tddft_print_notice (
     $           "MO and/or spatial CAP active: "//
     $           "Disabled initial MO Fock check "//
     $           "for moocc and/or dipcontribs")
         else
            if (.not. mat_is_purereal (g_zfock_mo(1), params%tol_zero))
     $           call errquit (pname//
     $           "initial Fock matrix in MO basis is not pure real",0,0)
            
            if (.not. mat_is_purereal (g_zfock_mo(2), params%tol_zero))
     $           call errquit (pname//
     $           "initial Fock matrix in MO basis is not pure real",0,0)
         endif

         call rt_tddft_movecs_gs (params, g_zfock_mo(1),
     $        g_movecs_mo_gs(1))
         call rt_tddft_movecs_gs (params, g_zfock_mo(2),
     $        g_movecs_mo_gs(2))
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
         call rt_tddft_restart_load (params, it, tt, 2, g_zdens_ao) ! 2 for uks
         call canorg_trans (params, "P", "AO->MO",
     $        g_zdens_ao(1), g_zdens_mo(1))
         call canorg_trans (params, "P", "AO->MO",
     $        g_zdens_ao(2), g_zdens_mo(2))
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
      if (.not. rt_tddft_os_tdfock (params, tt, g_zdens_ao,
     $     energies, g_zfock_ao))
     $     call errquit (pname//"tdfock failed", 0, 0)
      
      call canorg_trans (params, "F","AO->MO",
     $     g_zfock_ao(1), g_zfock_mo(1))
      call canorg_trans (params, "F","AO->MO",
     $     g_zfock_ao(2), g_zfock_mo(2))
      
      call ga_copy (g_zfock_mo(1), g_zfock_mo_n12dt(1))
      call ga_copy (g_zfock_mo(2), g_zfock_mo_n12dt(2))
         


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
C     Convert alpha,beta density matrices to AO basis and unpack real
C     and im parts.
C         
CXXX  [KAL]: this is perhaps wasteful; if we are not printing or doing checks there should be no need to do this here.
         call canorg_trans (params, "P","MO->AO",
     $        g_zdens_mo(1), g_zdens_ao(1))
         call canorg_trans (params, "P","MO->AO",
     $        g_zdens_mo(2), g_zdens_ao(2))
         
         call convert_z2d(g_zdens_ao(1), g_densre_ao(1), g_densim_ao(1))
         call convert_z2d(g_zdens_ao(2), g_densre_ao(2), g_densim_ao(2))

         call convert_z2d(g_zdens_mo(1), g_densre_mo(1), g_densim_mo(1))
         call convert_z2d(g_zdens_mo(2), g_densre_mo(2), g_densim_mo(2))


C
C     == Occasionally do some checks. == 
C
C     For the check we build new Fock matrix in AO basis and then
C     convert to MO basis (wasteful but we dont do checks very often).
C     

         if (rt_tddft_atpoint(it, params%nt, params%nchecks)) then
            
            if (params%prof) call prof_start (elapsed_checks)
            
            if (.not. rt_tddft_os_tdfock (params, tt, g_zdens_ao,
     $           energies, g_zfock_ao))
     $           call errquit (pname//"tdfock failed", 0, 0)
            
            call ga_zero (g_zfock_mo)

            call canorg_trans (params, "F","AO->MO",
     $           g_zfock_ao(1), g_zfock_mo(1))
            call canorg_trans (params, "F","AO->MO",
     $           g_zfock_ao(2), g_zfock_mo(2))


            if (.not. rt_tddft_checkpoint (params, 2, g_zfock_ao,
     $           g_zdens_ao, g_zfock_mo, g_zdens_mo)) then
               call rt_tddft_print_msg (params, tt,
     $              "WARNING: CHECKS FAILED")
            else
               call rt_tddft_print_msg (params, tt, "Checks passed")
            endif
            
c$$$            call rt_tddft_os_checkpoint (params, g_zfock_ao, g_zdens_ao,
c$$$     $           g_zfock_mo, g_zdens_mo)

            
            if (params%prof) call prof_end (elapsed_checks, "Checks")

C     (prints estimated time remaining)
            if (it .ne. est_it_start)
     $           call rt_tddft_estimate_rem (params, est_cputime_start,
     $           est_it_start, it, tt)
         endif


C
C     == Occasionally dump real part of *total* a+b AO density matrix. ==
C
C
C         if (rt_tddft_atpoint(it, params%nt, params%nsnapshots)) then
         if (rt_tddft_at_snap (params, tt)) then
            
            call ga_zero (g_densre_ao_tot)
            call ga_add (1d0, g_densre_ao(1), 1d0, g_densre_ao(2),
     $           g_densre_ao_tot)
            
            call rt_tddft_snapshot_save(params, it, tt, g_densre_ao_tot)
            call rt_tddft_print_msg (params, tt,
     $           "Saved density matrix snapshot")
         endif



C
C     == Occasionally calc and print time-dependent system properties. ==
C     
         if (rt_tddft_atpoint(it, params%nt, params%nprints)) then
            if (params%prof) call prof_start (elapsed_props)

            call ga_zero (g_densre_ao_tot)
            call ga_add (1d0, g_densre_ao(1), 1d0, g_densre_ao(2),
     $           g_densre_ao_tot)  !done before?


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


C     (charge on each geom)
            if (params%print_charge) then 
               call rt_tddft_print_charge (params, tt,
     $              g_densre_ao(1), "(alpha spin)")
               call rt_tddft_print_charge (params, tt,
     $              g_densre_ao(2), "(beta spin)")
               call rt_tddft_print_charge (params, tt,
     $              g_densre_ao_tot, "(total spin)")
            endif


C     (molecular orbital occupations, calc in MO basis)
CXXX  [KAL]: I guess really no need to allocate two moooc arrays...
            if (params%print_moocc) then

               call rt_tddft_moocc_calc (params, g_densre_mo(1),
     $              g_movecs_mo_gs(1), dbl_mb(imoocc(1)))
               call rt_tddft_moocc_print (params,
     $              tt, dbl_mb(imoocc(1)),
     $              "# MO Occupations (alpha spin)")

               call rt_tddft_moocc_calc (params, g_densre_mo(2),
     $              g_movecs_mo_gs(2), dbl_mb(imoocc(2)))        !was incorrectly gs(1) before
               call rt_tddft_moocc_print (params,
     $              tt, dbl_mb(imoocc(2)),
     $              "# MO Occupations (beta spin)")

c$$$               call rt_tddft_complex_evals (params, tt, g_zfock_mo(1),
c$$$     $              "Fock MO (alpha)")
c$$$               call rt_tddft_complex_evals (params, tt, g_zfock_mo(2),
c$$$     $              "Fock MO (beta)")
               
            endif

CXXX [ASB]: Transition Dipole Contributions 
            if (params%print_dipcontribs) then
               call rt_tddft_calc_dip_contribs (params,
     $              g_movecs_mo_gs(1), g_zdens_ao(1), tt, "ALPHA")
               call rt_tddft_calc_dip_contribs (params,
     $              g_movecs_mo_gs(2), g_zdens_ao(2), tt, "BETA")
            endif


C     (<S^2> of entire active geom, calc in MO basis)
            if (params%print_s2) then
               call rt_tddft_calc_s2 (params, g_densre_mo, s2val)
               call rt_tddft_print_val (params, tt, s2val, "<S^2>")
            endif

            
C     (dipole moment, calc in AO basis)
C
C     note: I dont think dipmom_alpha + dipmom_beta = dip_total (note I
C     use the total density matrix to compute dipmom_tot)
C     
            if (params%print_dipole) then
               do ig = 1, params%ngeoms
                  
                  if (params%geom(ig)%measure) then

                     outstr = "Dipole moment [" //
     $                    trim (params%geom(ig)%name) //
     $                    "] (alpha spin)"
                     call rt_tddft_calc_dip_geom (params,
     $                    g_densre_ao(1), ig, dipmom_alpha)
                     call rt_tddft_print_vector (params, tt,
     $                    dipmom_alpha, outstr)

                     outstr = "Dipole moment [" //
     $                    trim (params%geom(ig)%name) //
     $                    "] (beta spin)"
                     call rt_tddft_calc_dip_geom (params,
     $                    g_densre_ao(2), ig, dipmom_beta)
                     call rt_tddft_print_vector (params, tt,
     $                    dipmom_beta, outstr)

                     outstr = "Dipole moment [" //
     $                    trim (params%geom(ig)%name) //
     $                    "] (total spin)"
                     call rt_tddft_calc_dip_geom (params,
     $                    g_densre_ao_tot, ig, dipmom_tot)
                     call rt_tddft_print_vector (params, tt,
     $                    dipmom_tot, outstr)
                     
                  endif         !end "measure" this geom
               enddo            !loop over geoms
            endif
            
C     (external applied field, also calcd in tdfock, note assume applied
C     equally to both spins)
            if (params%print_field)
     $           call rt_tddft_print_fields (params, tt)

            if (params%prof) call prof_end (elapsed_props, "Properties")
         endif  !end print point

         
C
C     Dump restart information to disk.
C
         if (rt_tddft_atpoint(it, params%nt, params%nrestarts))
     $        call rt_tddft_restart_save (params, it, tt, 2, g_zdens_ao)


C     
C     Step MO density matrix forward dt in time.
C     
         if (params%prop_method .eq. 1) then
            call prop_euler (params, tt, rt_tddft_os_tdfock,
     $           2, g_zdens_mo, energies)

         elseif (params%prop_method .eq. 2) then
            call prop_rk4 (params, tt, rt_tddft_os_tdfock,
     $           2, g_zdens_mo, energies)
            
         elseif (params%prop_method .eq. 3) then
            call prop_magnus_ord2_interpol (params, tt,
     $           rt_tddft_os_tdfock, 2,
     $           g_zfock_mo_n12dt, g_zfock_mo, energies, g_zdens_mo)

         else
            call errquit (pname//"bad prop_method", 0, 0)
         endif

         
C
C     We are now at t + dt.
C
         it = it + 1
C         tt = params%tmin + (it-1)*params%dt
         tt = tt + params%dt

         
         if (params%prof) call prof_end (elapsed, "Complete time step")

      enddo

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

      do is = 1, nspin
         ok = ok .and. ga_destroy(g_zdens_ao(is))
         ok = ok .and. ga_destroy(g_zfock_ao(is))
         ok = ok .and. ga_destroy(g_zdens_mo(is))
         ok = ok .and. ga_destroy(g_zfock_mo(is))
         ok = ok .and. ga_destroy(g_zfock_mo_n12dt(is))
         ok = ok .and. ga_destroy(g_densre_ao(is))
         ok = ok .and. ga_destroy(g_densim_ao(is))
         ok = ok .and. ga_destroy(g_densre_mo(is))
         ok = ok .and. ga_destroy(g_densim_mo(is))
         ok = ok .and. ga_destroy(g_movecs_mo_gs(is))
         ok = ok .and. ga_destroy(g_zscr_mo(is))
         if (params%static) ok=ok.and.ga_destroy (g_zfock_ao_init(is))
      enddo
      
      ok = ok .and. ga_destroy(g_densre_ao_tot)
      
      if (.not.ok) call errquit(pname//"failed to destroy GAs", 0,0)

      if (.not.ma_pop_stack(lmoocc(2)))
     &   call errquit(pname//'cannot pop moocc(2)',0, MA_ERR)

      if (.not.ma_pop_stack(lmoocc(1)))
     &   call errquit(pname//'cannot pop moocc(1)',0, MA_ERR)
      
      
      end subroutine
c $Id$
