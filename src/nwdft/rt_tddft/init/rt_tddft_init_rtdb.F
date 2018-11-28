C     
C     Read rt-tddft parameters from rtdb and load into params struct.
C     
      subroutine rt_tddft_init_rtdb (params)
      implicit none

#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "rtdb.fh"
#include "rt_tddft.fh"


C     == In/out ==
      type(rt_params_t) params

C     == Parameters ==
      character(*), parameter :: pname = "rt_tddft_init_rtdb: "

C     == Variables ==
      integer ifield, iexcite

      
      if (.not. rtdb_get (params%rtdb, "rt_tddft:dt",
     $     mt_dbl, 1, params%dt)) call errquit (pname//
     $     "failed to read dt from rtdb", 0, RTDB_ERR)
      
      if (.not. rtdb_get (params%rtdb, "rt_tddft:tmin",
     $     mt_dbl, 1, params%tmin)) call errquit (pname//
     $     "failed to read tmin from rtdb", 0, RTDB_ERR)

      if (.not. rtdb_get (params%rtdb, "rt_tddft:tmax",
     $     mt_dbl, 1, params%tmax)) call errquit (pname//
     $     "failed to read tmax from rtdb", 0, RTDB_ERR)

      if (.not.rtdb_get(params%rtdb, "rt_tddft:nchecks",
     $     mt_int, 1, params%nchecks)) call errquit (pname//
     $     "failed to read nchecks from rtdb", 0, RTDB_ERR)

      if (.not.rtdb_get(params%rtdb, "rt_tddft:checklvl",
     $     mt_int, 1, params%checklvl)) call errquit (pname//
     $     "failed to read nchecks from rtdb", 0, RTDB_ERR)
      
      if (.not.rtdb_get(params%rtdb, "rt_tddft:nprints",
     $     mt_int, 1, params%nprints)) call errquit (pname//
     $     "failed to read nprints from rtdb", 0, RTDB_ERR)

c$$$      if (.not.rtdb_get(params%rtdb, "rt_tddft:nsnapshots",
c$$$     $     mt_int, 1, params%nsnapshots)) call errquit (pname//
c$$$     $     "failed to read nsnapshots from rtdb", 0, RTDB_ERR)

      if (.not.rtdb_get(params%rtdb, "rt_tddft:nrestarts",
     $     mt_int, 1, params%nrestarts)) call errquit (pname//
     $     "failed to read nrestarts from rtdb", 0, RTDB_ERR)

      if (.not. rtdb_get (params%rtdb, "rt_tddft:prop_method",
     $     mt_int, 1, params%prop_method)) call errquit (pname//
     $     "failed to read prop_method from rtdb", 0, RTDB_ERR)

      if (.not.rtdb_get(params%rtdb, "rt_tddft:prof",
     $     mt_log, 1, params%prof)) call errquit (pname//
     $     "failed to read prof from rtdb", 0, RTDB_ERR)

      if (.not.rtdb_get(params%rtdb, "rt_tddft:matrix_checks",
     $     mt_log, 1, params%matrix_checks)) call errquit (pname//
     $     "failed to read matrix_checks from rtdb", 0, RTDB_ERR)

      if (.not.rtdb_get(params%rtdb, "rt_tddft:nodisk",
     $     mt_log, 1, params%nodisk)) call errquit (pname//
     $     "failed to read nodisk from rtdb", 0, RTDB_ERR)

      if (.not.rtdb_get(params%rtdb, "rt_tddft:noprop",
     $     mt_log, 1, params%noprop)) call errquit (pname//
     $     "failed to read noprop from rtdb", 0, RTDB_ERR)
      
      if (.not.rtdb_get(params%rtdb, "rt_tddft:static",
     $     mt_log, 1, params%static)) call errquit (pname//
     $     "failed to read static from rtdb", 0, RTDB_ERR)

      if (.not.rtdb_get(params%rtdb, "rt_tddft:use_dmat",
     $     mt_log, 1, params%use_dmat)) call errquit (pname//
     $     "failed to read use_dmat from rtdb", 0, RTDB_ERR)

      if (.not.rtdb_get(params%rtdb, "rt_tddft:viz_active",
     $     mt_log, 1, params%viz_active)) call errquit (pname//
     $     "failed to read viz_active from rtdb", 0, RTDB_ERR)

      if (params%viz_active) then
         if (.not.rtdb_get(params%rtdb, "rt_tddft:viz_dplot",
     $        mt_log, 1, params%viz_dplot)) call errquit (pname//
     $        "failed to read viz_dplot from rtdb", 0, RTDB_ERR)

         if (.not.rtdb_get(params%rtdb, "rt_tddft:viz_subgs",
     $        mt_log, 1, params%viz_subgs)) call errquit (pname//
     $        "failed to read viz_subgs from rtdb", 0, RTDB_ERR)
         
         if (.not. rtdb_get (params%rtdb, "rt_tddft:viz_tstart",
     $        mt_dbl, 1, params%viz_tstart)) call errquit (pname//
     $        "failed to read viz_tstart from rtdb", 0, RTDB_ERR)
         
         if (.not. rtdb_get (params%rtdb, "rt_tddft:viz_tend",
     $        mt_dbl, 1, params%viz_tend)) call errquit (pname//
     $        "failed to read viz_tend from rtdb", 0, RTDB_ERR)

         if (.not. rtdb_get (params%rtdb, "rt_tddft:viz_treference",
     $        mt_dbl, 1, params%viz_treference)) call errquit (pname//
     $        "failed to read viz_treference from rtdb", 0, RTDB_ERR)

      endif
         
      if (.not.rtdb_get(params%rtdb, "rt_tddft:print_dipole",
     $     mt_log, 1, params%print_dipole)) call errquit (pname//
     $     "failed to read print_dipole from rtdb", 0, RTDB_ERR)

      if (.not.rtdb_get(params%rtdb, "rt_tddft:print_quad",
     $     mt_log, 1, params%print_quad)) call errquit (pname//
     $     "failed to read print_quad from rtdb", 0, RTDB_ERR)

      if (.not.rtdb_get(params%rtdb, "rt_tddft:print_field",
     $     mt_log, 1, params%print_field)) call errquit (pname//
     $     "failed to read print_field from rtdb", 0, RTDB_ERR)

      if (.not.rtdb_get(params%rtdb, "rt_tddft:print_moocc",
     $     mt_log, 1, params%print_moocc)) call errquit (pname//
     $     "failed to read print_moocc from rtdb", 0, RTDB_ERR)

      if (.not.rtdb_get(params%rtdb, "rt_tddft:print_dipcontribs",
     $     mt_log, 1, params%print_dipcontribs)) call errquit (pname//
     $     "failed to read print_dipcontribs from rtdb", 0, RTDB_ERR)

      if (.not.rtdb_get(params%rtdb, "rt_tddft:print_energy",
     $     mt_log, 1, params%print_energy)) call errquit (pname//
     $     "failed to read print_energy from rtdb", 0, RTDB_ERR)

      if (.not.rtdb_get(params%rtdb, "rt_tddft:print_cputime",
     $     mt_log, 1, params%print_cputime)) call errquit (pname//
     $     "failed to read print_cputime from rtdb", 0, RTDB_ERR)

      if (.not.rtdb_get(params%rtdb, "rt_tddft:print_charge",
     $     mt_log, 1, params%print_charge)) call errquit (pname//
     $     "failed to read print_charge from rtdb", 0, RTDB_ERR)

      if (.not.rtdb_get(params%rtdb, "rt_tddft:print_convergence",
     $     mt_log, 1, params%print_convergence)) call errquit (pname//
     $     "failed to read print_convergence from rtdb", 0, RTDB_ERR)

      if (.not.rtdb_get(params%rtdb, "rt_tddft:print_s2",
     $     mt_log, 1, params%print_s2)) call errquit (pname//
     $     "failed to read print_s2 from rtdb", 0, RTDB_ERR)

      if (.not. rtdb_get (params%rtdb, "rt_tddft:exp_method",
     $     mt_int, 1, params%exp_method))
     $     call errquit (pname//"failed to read exp_method from rtdb",
     $     0, RTDB_ERR)

      if (.not. rtdb_get (params%rtdb, "rt_tddft:tol_zero",
     $     mt_dbl, 1, params%tol_zero)) call errquit (pname//
     $     "failed to read tol_zero from rtdb", 0, RTDB_ERR)

      if (.not. rtdb_get (params%rtdb, "rt_tddft:tol_series",
     $     mt_dbl, 1, params%tol_series)) call errquit (pname//
     $     "failed to read tol_series from rtdb", 0, RTDB_ERR)

      if (.not. rtdb_get (params%rtdb, "rt_tddft:tol_interpol",
     $     mt_dbl, 1, params%tol_interpol)) call errquit (pname//
     $     "failed to read tol_interpol from rtdb", 0, RTDB_ERR)

      if (.not.rtdb_get(params%rtdb, "rt_tddft:terms_series",
     $     mt_int, 1, params%terms_series)) call errquit (pname//
     $     "failed to read terms_series from rtdb", 0, RTDB_ERR)

      if (.not.rtdb_get(params%rtdb, "rt_tddft:terms_interpol",
     $     mt_int, 1, params%terms_interpol)) call errquit (pname//
     $     "failed to read interpol from rtdb", 0, RTDB_ERR)

      if (.not.rtdb_cget(params%rtdb,"rt_tddft:tag",1,params%tag))
     $     call errquit(pname//'Read failed for tag from rtdb',
     $     0,RTDB_ERR)


C
C     Read in fields
C
      if (.not.rtdb_get(params%rtdb, "rt_tddft:nfields",
     $     mt_int, 1, params%nfields)) call errquit (pname//
     $     "failed to read nfields from rtdb", 0, RTDB_ERR)

      if (params%nfields .gt. rt_max_fields)
     $     call errquit (pname//"nfields > rt_max_fields", 0, 0)

      do ifield = 1, params%nfields
         call rt_tddft_field_rtdb_get (params%rtdb, ifield,
     $        params%field(ifield))
      enddo


C
C     Read in excitation rules
C
      if (.not.rtdb_get(params%rtdb, "rt_tddft:nexcites",
     $     mt_int, 1, params%nexcites)) call errquit (pname//
     $     "failed to read nexcites from rtdb", 0, RTDB_ERR)

      if (params%nexcites .gt. rt_max_excites)
     $     call errquit (pname//"nexcites > rt_max_excites", 0, 0)

      do iexcite = 1, params%nexcites
         call rt_tddft_excite_rtdb_get (params%rtdb, iexcite,
     $        params%excite(iexcite))
      enddo


C     
C     MO CAP
C
      if (.not.rtdb_get(params%rtdb, "rt_tddft:mocap_active",
     $     mt_log, 1, params%mocap_active)) call errquit (pname//
     $     "failed to read mocap_active from rtdb", 0, RTDB_ERR)

      if (params%mocap_active) then
         if (.not. rtdb_get (params%rtdb, "rt_tddft:mocap_emin",
     $        mt_dbl, 1, params%mocap_emin)) call errquit (pname//
     $        "failed to read mocap_emin from rtdb", 0, RTDB_ERR)

         if (.not. rtdb_get (params%rtdb, "rt_tddft:mocap_maxval",
     $        mt_dbl, 1, params%mocap_maxval)) call errquit (pname//
     $        "failed to read mocap_maxval from rtdb", 0, RTDB_ERR)

         if (.not. rtdb_get (params%rtdb, "rt_tddft:mocap_prefac",
     $        mt_dbl, 1, params%mocap_prefac)) call errquit (pname//
     $        "failed to read mocap_prefac from rtdb", 0, RTDB_ERR)

         if (.not. rtdb_get (params%rtdb, "rt_tddft:mocap_expconst",
     $        mt_dbl, 1, params%mocap_expconst)) call errquit (pname//
     $        "failed to read mocap_expconst from rtdb", 0, RTDB_ERR)

         if (.not.rtdb_get(params%rtdb, "rt_tddft:mocap_checks",
     $        mt_log, 1, params%mocap_checks)) call errquit (pname//
     $        "failed to read mocap_checks from rtdb", 0, RTDB_ERR)

         if (.not.rtdb_get(params%rtdb, "rt_tddft:mocap_print",
     $        mt_log, 1, params%mocap_print)) call errquit (pname//
     $        "failed to read mocap_print from rtdb", 0, RTDB_ERR)
      endif


C     
C     Spatial CAP
C
      if (.not.rtdb_get(params%rtdb, "rt_tddft:cap_active",
     $     mt_log, 1, params%cap_active)) call errquit (pname//
     $     "failed to read cap_active from rtdb", 0, RTDB_ERR)

C     if (params%cap_active) then
C     XXX LOAD PARAMS

      
      end subroutine
c $Id$
