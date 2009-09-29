c
c     == Ehrenfest dynamics - SODFT driver ==
c
      logical function ehrenfest_sodft
     &                 (rtdb, Etold, Enuc, iVcoul_opt, iVxc_opt, 
     &                  iter, g_dens, g_dens_at, g_movecs, g_vxc, 
     &                  g_fock, g_svecs, svals, g_xcinv, g_s,
     &                  g_dipx, g_dipy, g_dipz, g_ddotE, g_vel)
c     
c     $Id: ehrenfest_sodft.F,v 1.54 2008-12-31 01:37:01 niri Exp $
c     
      implicit none
c
#include "errquit.fh"
c     
      integer rtdb              ! [input]
      double precision Etold, Enuc, trace
      integer iVcoul_opt
      integer iVxc_opt
      integer iter, swap(20),no_swap, no_det, idet
      integer g_dens(2), g_movecs(2), g_vxc(4), 
     &     g_fock, g_svecs, g_xcinv,  g_scr
c
c     For time propagation
      integer g_dipx    ! x-dipole integrals
      integer g_dipy    ! y-dipole integrals
      integer g_dipz    ! z-dipole integrals
      integer g_ddotE   ! D.E matrix
      integer g_vel     ! velocity integrals
      integer g_sp1     ! S^(1) matrix
      integer g_sm12    ! S^(-1/2) matrix
      integer g_sp12    ! S^(1/2) matrix
      integer g_sm1     ! S^(-1) matrix
      double precision  dipole(3),jcurr(3)
c
      integer nwork
      parameter ( nwork = 2 )
      integer g_work(nwork) ! work arrays
      integer g_dprdt, g_dpidt  ! real and imaginary parts of dP/dt
      integer g_temp
c
      logical rk4
      integer irk
      double precision xh
      double precision xh2
      double precision h2
      double precision h6
      integer g_k1(nwork) ! work arrays
      integer g_k2(nwork) ! work arrays
      integer g_k3(nwork) ! work arrays
      integer g_k4(nwork) ! work arrays
      integer g_kwork(nwork) ! work arrays
c
      double precision current(3)
c  
      integer g_dens_at(2)
      double precision  svals(*)
c
c     so
      integer g_densso(2), g_tmp2,   
     &     g_moso(2), g_old(2),
     &     g_fockso(2), g_scr2, g_damp_so(2), g_gmovecs(2)
c
      integer la, ia            ! complex*16 a(nbf_ao, nbf_ao)
      integer lw, iw            ! double precision w(nbf_ao)
      integer llwork 
      integer lwork, iwork      ! complex*16 work(3) 
      integer lrwork, irwork    ! double precision rwork
      integer info 
      integer lbuff, ibuff
      integer nbf_mo 
      integer g_s
c
c     declarations for fractional occupation 
c
      integer nel_fon, nmo_fon  ! number of fractionally occupied electrons and orbitals
      double precision avg_fon  ! fractional occupancy (averaged)
      integer nTotOcc           ! nTotOcc: no. of (occupied maybe f.o.) mo's 
      logical fon 
      double precision scale 
      integer kfon_occ, lfon_occ
c
c     so
      double precision rho_n, toll_s
c
      character*32 pname  ! name of this subroutine
c
#include "bas.fh"
#include "geom.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "cdft.fh"
#include "global.fh"
#include "msgids.fh"
#include "util.fh"
#include "../so_dft/consts.fh"
c
      Logical movecs_write_so, movecs_converged
      External movecs_write_so, movecs_converged
c     
      Logical movecs_read_header_so, movecs_read_so
      External movecs_read_header_so, movecs_read_so 
c
      Logical spinor_guess
      External spinor_guess
      integer  ga_create_atom_blocked
      external ga_create_atom_blocked
c     
      logical oprint, status
      double precision Exc(2), rms(2), derr(2)
      integer noc(2), nmo(2), icall(2)
c      integer nva
      integer n3c_dbl, n3c_int, n_batch
      integer iwhat_max
      integer l_3cwhat, k_3cwhat, l_3cERI, k_3cERI
      integer dft_n3cint, n_semi_bufs, fd
      external dft_n3cint
      double precision dft_n3cdbl
      external dft_n3cdbl
      Integer l_eval
      integer k_eval(2)
      integer natoms, nTotEl

      integer l_occ, k_occ
      integer i, j,  i1, jstart 
      integer me, nproc
      integer g_tmp, g_fockt, g_wght, g_xyz,g_nq
c     so
      integer g_so(3)
c     so
      integer nheap, nstack
      integer ispin, idone
      integer nexc
      integer iswitc
      integer itol_max, iaoacc_max
      integer itol_min, iAOacc_min
      double precision tol_rho_min, tol_rho_max
      integer npol
      integer leneval, lcd_coef, icd_coef
      integer lcntoce, icntoce, lcntobfr, icntobfr,
     &     lcetobfr, icetobfr, lrdens_atom, irdens_atom,
     &     nscr, lscr, iscr
      double precision start_wall, current_wall, elapsed_wall,
     &     save_wall, current_cpu, start_cpu,
     &     wall_time_reqd
      integer int_wall_time_reqd
      double precision Ecoul, Ecore, Edip
      double precision Ecoresf, Ecoreso
      double precision pp, delta
      double precision anucl_charg, anel
      double precision anoca, anocb, onempp
      double precision Etnew, tol2e, tol2e_sleazy,Etotal
c     convergence declarations
      double precision rlshift_input, rlshift_def
      integer ndamp_input, ndamp_def
c     
c     Note, damping, levelshifting, and diising logicals
c     are used to turn on/off these procedures per
c     iteration.  The alternative logicals nodamping, 
c     nolevelshifting, and nodiis are specified and held
c     for the entire convergence sequence.
c     
      logical diising, damping, levelshifting
      logical keep_damp_on,keep_levl_on, keep_diis_on
      Logical  IOLGC, mulliken
      logical converged, wght_GA
c      logical oconverged 
      logical oprint_parm, oprint_conv, oprint_vecs, 
     &     oprint_eval, oprint_syma, oprint_time, 
     &     oprint_info, oprint_tol, oprint_final_vecs, 
     &     oprint_energy_step, oprint_intermediate_fock,
     &     oprint_3c2e, oprint_interm_overlap, oprint_interm_S2,
     &     oprint_conv_details
c     
      integer ilo, ihi          ! For printing movecs analysis
      double precision eval_pr_tol_lo, eval_pr_tol_hi
      parameter (eval_pr_tol_lo = -1.5d0, eval_pr_tol_hi=0.5)
c     
c     early convergence tolerances
c     
      parameter(itol_min = 7, iAOacc_min = 12, tol_rho_min = 1.d-7)
c     
      double precision dft_dencvg, dft_time
      external dft_dencvg
      double precision homo, lumo, homo_lumo_gap
      integer l_ir, k_ir
      logical last_time_energy
      logical check_shift, lmaxov_sv
      character*7 name
      character*4 scftype
      character*255 basis_name, basis_trans
      integer nopen, nclosed
c     !!! BGJ
      logical cphf_poliz, do_poliz
      external cphf_poliz
c     !!! BGJ
      character*255 title1       ! Returns title of job that created vectors
      character*255 basis_name1  ! Returns name of basis set
      character*255 scftype1     ! Returns the SCF type of the vectors
      integer nbf1               ! Returns no. of functions in basis
      integer g_oep
c     integer ijk
c     
      integer ncanorg
      logical ldmix
      character*7 vecs_or_dens
      logical spinor
      logical dft_mem3c
      external dft_mem3c
c
c     Time propagation 
      integer itime,ntime
      double precision tmin,tmax,tt,dt
c
      ehrenfest_sodft = .false.
      nbf_mo = 2*nbf_ao
      lmaxov_sv = lmaxov
      oprint = util_print('information', print_low)
      oprint_info = util_print('common', print_debug)
      oprint_parm = util_print('parameters', print_default)
      oprint_3c2e = util_print('3c 2e integrals', print_default)
      oprint_conv = util_print('convergence', print_default)
      oprint_conv_details = util_print('convergence details', 
     &     print_high)
      oprint_vecs = util_print('intermediate vectors', print_high)
      oprint_eval = util_print('intermediate evals', print_high)
      oprint_syma = util_print('interm vector symm', print_high)
      oprint_time = util_print('dft timings', print_high)
      oprint_tol = util_print('screening parameters', print_high)
      oprint_energy_step = util_print('intermediate energy info',
     &     print_high)
      oprint_intermediate_fock = util_print('intermediate fock matrix',
     &     print_high)
      oprint_interm_S2 = util_print('intermediate S2',print_high)
      oprint_interm_overlap = util_print('intermediate overlap',
     &     print_high)
      oprint_final_vecs = util_print('final vectors', print_high)
c
      me = ga_nodeid()
      nproc = ga_nnodes()
      ispin=1
      pname = 'ehrenfest_sodft:'

      call int_1e_uncache_ga()      
c     !!! BGJ
c     Store SCF hamiltonian type as DFT for use in BGJ routines
      if (.not. rtdb_put(rtdb, 'bgj:scf_type', MT_INT, 1, 2))
     $  call errquit(pname//'put of bgj:scf_type failed',0, RTDB_ERR)
c     !!! BGJ
c     
c     see if levelshifting monitoring is desired
c     
      if (.not. rtdb_get(rtdb, 'dft:check_shift', mt_log, 1,
     &     check_shift))then
         check_shift = .false.      
      endif
c     
      if (.not. geom_ncent(geom, natoms))
     &  call errquit(pname//'geom_ncent failed',73,GEOM_ERR)
      if (.not. geom_nuc_charge(geom, anucl_charg))
     &  call errquit(pname//'geom_nuc_charge failed',0,GEOM_ERR)
c     
      anel = int(anucl_charg) - rcharge
c     
c     Pre-compute mapping vectors
c     
      if (.not.ma_push_get
     &  (mt_int,nshells_ao,'cntoce map',lcntoce,icntoce))
     &  call errquit(pname//'push_get failed', 13, MA_ERR)
      if (.not.ma_push_get
     &  (mt_int,nshells_ao*2,'cntoce map',lcntobfr,icntobfr))
     &  call errquit(pname//'push_get failed', 13, MA_ERR)
      if (.not.ma_push_get
     &  (mt_int,natoms*2,'cntoce map',lcetobfr,icetobfr))
     &  call errquit(pname//'push_get failed', 13, MA_ERR)
c     
      call build_maps(ao_bas_han, int_mb(icntoce), int_mb(icntobfr), 
     &     int_mb(icetobfr), natoms, nshells_ao)
c     
c     Set aside some memory for reduced density matrix
c     
      if (.not.MA_Push_Get(MT_Dbl,2*natoms*natoms,'rdens_atom',
     & lrdens_atom,irdens_atom))
     & call errquit(pname//'cannot allocate rdens_atom',0,MA_ERR)
c     
c     determine pattern of orbitals' occupancy
c     
      if (.not. rtdb_get(rtdb, 'dft:noc', mt_int, 2, noc))
     & call errquit(pname//'rtdb_get of noc failed',0,RTDB_ERR)
      if (.not. MA_Push_Get(MT_Dbl,nbf_ao*2,'mo occ',l_occ,k_occ))
     & call errquit(pname//'failed to alloc',999, MA_ERR)
c     
c     get orbital overlap tolerance
c     
      if (.not. rtdb_get(rtdb, 'dft:toll_s', MT_DBL, 1, toll_s))
     & call errquit(pname//'lost toll_s',0, RTDB_ERR)
c     
      nTotEl = noc(1) + noc(2)
      nmo(1) = nbf_ao
      nmo(2) = nbf_ao
c     
      anoca = noc(1)
      anocb = noc(2)
c     
c     UHF occupations
c     
      call dfill(nbf_mo, 0.0d0, dbl_mb(k_occ), 1)
      do i = 1, noc(1)
         dbl_mb(i-1+k_occ) = 1.0d0
      enddo
      do i = nbf_ao+1, nbf_ao+noc(2)
         dbl_mb(i-1+k_occ) = 1.0d0
      enddo
c     
      wght_GA = .false.
c     
c     Determine whether to fit the electronic charge density.
c     
      CDFIT = .FALSE.
      if (iVcoul_opt.eq.1)CDFIT = .TRUE.
      XCFIT = .FALSE.
      if (iVxc_opt.eq.1)XCFIT = .TRUE.
c
c     Define various constants.
c     
      npol = (ipol*(ipol+1))/2
c     
      itol_max = itol2e
      iaoacc_max = iaoacc
      tol_rho_max = tol_rho
      mulliken = .false.
      if (imull.eq.1)mulliken = .true.
      IOLGC = .TRUE.
      if (noio.eq.1)IOLGC = .FALSE.
c     
c     Energy decomposition switch
c     
      nExc    = idecomp + 1
      Etnew = 0.d0
c     
c     SCF energy convergence criterion. 
c     
      g_tmp = ga_create_atom_blocked(geom, AO_bas_han, 'ga_temp')
      g_fockt = ga_create_atom_blocked(geom, AO_bas_han, 'fock tr')
c     
c     Set up local convergence parameters
c     
      diising = diis
      damping = damp
      levelshifting = levelshift
      keep_damp_on = .false.
      keep_levl_on = .false.
      keep_diis_on = .false.
      ndamp_input = ndamp
      rlshift_input = rlshift
      ndamp_def = 0
      rlshift_def = 0.0
      rlshift = rlshift_def
c     
      if (nodamping)damping = .false.
      if (nolevelshifting) then 
         levelshifting = .false.
         rlshift = rlshift_def
      endif
      if (nodiis) then
         diising = .false.
      endif
      if (ncydp.ne.0)then
         damping = .true. 
         ndamp = ndamp_input
      endif
      if (ncysh.ne.0)then
         levelshifting = .true.
         rlshift = rlshift_input
      endif
      if (ncyds.ne.0)then
         diising = .true.
      endif
c     
c     Initialize DIIS call counter.
c     
      icall(1) = 0
      icall(2) = 0
c     
c     Begin the SCF cycle.
c     
c     allocate eigenvalue array
c     
      leneval = 4*nbf_ao 
      if (.not.MA_Push_Get(MT_Dbl,leneval,'eval',l_eval,k_eval(1)))
     & call errquit(pname//'Cannot allocate eval',0, MA_ERR)
      k_eval(2) = k_eval(1) + nbf_mo
c     
c     Dump DFT parameters (if debugging) to see if they make sense
c     
      if (me.eq.0.and.oprint_info)call dft_dump_info(me)
c     
c     Get initial density.
c     
      if (oprint_time)
     &     call dft_tstamp(' Before call to DFT_INIT.')
      scftype = 'UHF'
c     
c     allocate array for irreps
c     
      if (.not.MA_Push_Get(mt_int,2*nbf_ao,'dft:irreps',l_ir,k_ir))
     & call errquit(pname//'Cannot allocate irreps',0, MA_ERR)
      nopen = mult - 1
      nclosed = (nTotEl - nopen) / 2
c     
      if (.not. bas_name(ao_bas_han, basis_name, basis_trans))
     $ call errquit(pname//'bas_name?', 0, BASIS_ERR)
c     
c     get info for int2e_ and set sleazy tolerance
c     
      tol2e_sleazy = 1.d-3
      call scf_get_fock_param(rtdb, tol2e_sleazy)
c     
c     Force sleazy SCF into "direct" mode.
c     
      call fock_force_direct(rtdb)
c
c     Allocate Fock matrix, movecs, density matrices, scratch areas, etc
c
c     Real molecular orbital vectors
      if(.not.ga_create(mt_dbl,nbf_mo,nbf_mo,'Movecs Re',0,0,g_moso(1)))     
     & call errquit(pname//'Error creating Movecs Re',0, GA_ERR)
      call ga_zero(g_moso(1))
c
      if(.not.ga_create(mt_dbl,nbf_ao,nbf_ao,'Movecs Re',0,0, 
     &     g_gmovecs(1)))     
     & call errquit(pname//'Error creating Movecs Re',0, GA_ERR)
      call ga_zero(g_gmovecs(1))
c
c     Imaginary molecular orbital vectors
      if(.not.ga_create(mt_dbl,nbf_mo, nbf_mo,'Movecs Im',0,0, 
     &     g_moso(2)))
     & call errquit(pname//'Error creating Movecs Im',0, GA_ERR)
      call ga_zero(g_moso(2))
c
      if(.not.ga_create(mt_dbl,nbf_ao, nbf_ao,'Movecs Im',0,0, 
     &     g_gmovecs(2)))
     & call errquit(pname//'Error creating Movecs Im',0, GA_ERR)
      call ga_zero(g_gmovecs(2))
c
c     Real part of the fock matrix
      if(.not.ga_create(mt_dbl,nbf_mo,nbf_mo,'Fock Re',0,0, 
     & g_fockso(1)))
     & call errquit(pname//'Error creating Fock Re',0, GA_ERR)
      call ga_zero(g_fockso(1))
c
c     Imaginary part of the fock matrix
      if(.not.ga_create(mt_dbl,nbf_mo,nbf_mo,'Fock Im',0,0, 
     & g_fockso(2)))
     & call errquit(pname//'Error creating Fock Im',0, GA_ERR)
      call ga_zero(g_fockso(2))
c
c     Real part of the spin-orbit density matrix
      if(.not.ga_create(mt_dbl,nbf_mo,nbf_mo,'DenMx Re',0,0, 
     & g_densso(1)))
     & call errquit(pname//'Error creating DenMx Re',0, GA_ERR)
      call ga_zero(g_densso(1))
c
c     Imaginary part of the spin-orbit density matrix
      if(.not.ga_create(mt_dbl,nbf_mo,nbf_mo,'DenMx Im',0,0, 
     & g_densso(2)))
     & call errquit(pname//'Error creating DenMx Im',0, GA_ERR)
      call ga_zero(g_densso(2))
c
c     S^(1)
      if(.not.ga_create(mt_dbl, nbf_mo, nbf_mo,'S^(1)', 0, 0, g_sp1))
     & call errquit(pname//'Error creating g_sm1',0, GA_ERR)
      call ga_zero(g_sp1)
c
c     S^(-1)
      if(.not.ga_create(mt_dbl, nbf_mo, nbf_mo,'S^(-1)', 0, 0, g_sm1))
     & call errquit(pname//'Error creating g_sm1',0, GA_ERR)
      call ga_zero(g_sm1)
c
c     S^(1/2)
      if(.not.ga_create(mt_dbl, nbf_mo, nbf_mo,'S^(1/2)', 0, 0, g_sp12))
     & call errquit(pname//'Error creating g_sm1',0, GA_ERR)
      call ga_zero(g_sp12)
c
c     S^(-1/2)
      if(.not.ga_create(mt_dbl,nbf_mo,nbf_mo,'S^(-1/2)', 0, 0, g_sm12))
     & call errquit(pname//'Error creating g_sm1',0, GA_ERR)
      call ga_zero(g_sm12)
c
c     Spin-orbit matrices: 1->z, 2->y, 3->x
      if(.not.ga_create(mt_dbl,nbf_ao,nbf_ao,'so z',0,0, 
     &  g_so(1)))
     &  call errquit(pname//'Error creating so z',0, GA_ERR)
      call ga_zero(g_so(1))
      if(.not.ga_create(mt_dbl,nbf_ao,nbf_ao,'so y',0,0, 
     &  g_so(2)))
     &  call errquit(pname//'Error creating so y',0, GA_ERR)
      call ga_zero(g_so(2))
      if(.not.ga_create(mt_dbl,nbf_ao,nbf_ao,'so x',0,0, 
     &  g_so(3)))
     &  call errquit(pname//'Error creating so x',0, GA_ERR)
      call ga_zero(g_so(3))
c
c     Extra arrays
      if(.not.ga_create(mt_dbl,nbf_mo,nbf_mo,'old den', 0, 0, 
     & g_damp_so(1)))
     & call errquit(pname//'Error creating damp ga', 0, GA_ERR)
      call ga_zero(g_damp_so(1))
      if(.not.ga_create(mt_dbl,nbf_mo,nbf_mo,'old den', 0, 0, 
     & g_damp_so(2)))
     & call errquit(pname//'Error creating damp ga', 0, GA_ERR)
      call ga_zero(g_damp_so(2))
      if(.not.ga_create(mt_dbl,nbf_mo,nbf_mo,'scr2', 0, 0, g_scr2))
     & call errquit(pname//'Error creating scr2',0, GA_ERR)
      call ga_zero(g_scr2)
      if(.not.ga_create(mt_dbl,nbf_mo,nbf_mo,'old re',0,0, 
     & g_old(1)))
     & call errquit(pname//'Error creating Old Re',0, GA_ERR)
      call ga_zero(g_old(1))
      if(.not.ga_create(mt_dbl,nbf_mo,nbf_mo,'old Im',0,0, 
     & g_old(2)))
     & call errquit(pname//'Error creating Old Im',0, GA_ERR)
      call ga_zero(g_old(2))
      if(.not.ga_create(mt_dbl,nbf_mo,nbf_mo,'Tmp ReIm',0,0, 
     & g_tmp2))
     & call errquit(pname//'Error creating Tmp ReIm',0, GA_ERR)
      call ga_zero(g_tmp2)
      g_scr = ga_create_atom_blocked(geom, AO_bas_han, 'ga scr')
      call ga_zero(g_scr)
c
c     Workspaces for time propagation
c
        if(.not.ga_create(mt_dbl,nbf_ao,nbf_ao,'g_temp',0,0, g_temp))
     &  call errquit(pname//'Error creating g_temp',0, GA_ERR)
        call ga_zero(g_temp)

      do iwork = 1,nwork
        if(.not.ga_create(mt_dbl,nbf_mo,nbf_mo,'g_work',0,0,
     & g_work(iwork)))
     &  call errquit(pname//'Error creating g_work',0, GA_ERR)
        call ga_zero(g_work(iwork))
c
        if(.not.ga_create(mt_dbl,nbf_mo,nbf_mo,'g_k1',0,0,
     & g_k1(iwork)))
     &  call errquit(pname//'Error creating g_k1',0, GA_ERR)
        call ga_zero(g_k1(iwork))
c
        if(.not.ga_create(mt_dbl,nbf_mo,nbf_mo,'g_k2',0,0,
     & g_k2(iwork)))
     &  call errquit(pname//'Error creating g_k2',0, GA_ERR)
        call ga_zero(g_k2(iwork))
c
        if(.not.ga_create(mt_dbl,nbf_mo,nbf_mo,'g_k3',0,0,
     & g_k3(iwork)))
     &  call errquit(pname//'Error creating g_k3',0, GA_ERR)
        call ga_zero(g_k3(iwork))
c
        if(.not.ga_create(mt_dbl,nbf_mo,nbf_mo,'g_k4',0,0,
     & g_k4(iwork)))
     &  call errquit(pname//'Error creating g_k4',0, GA_ERR)
        call ga_zero(g_k4(iwork))
c
        if(.not.ga_create(mt_dbl,nbf_mo,nbf_mo,'g_kwork',0,0,
     & g_kwork(iwork)))
     &  call errquit(pname//'Error creating g_kwork',0, GA_ERR)
        call ga_zero(g_kwork(iwork))
c
      end do
c
c     dPr/dt
      if(.not.ga_create(mt_dbl,nbf_mo, nbf_mo,'g_dprdt',0,0,g_dprdt))
     &  call errquit(pname//'Error creating g_dprdt',0, GA_ERR)
      call ga_zero(g_dprdt)
c
c     dPi/dt
      if(.not.ga_create(mt_dbl,nbf_mo,nbf_mo,'g_dpidt',0,0,g_dpidt))
     &  call errquit(pname//'Error creating g_dpidt',0, GA_ERR)
      call ga_zero(g_dpidt)
c
      call ga_sync()
c
c     Calculate S powers: S,S^(-1/2),S^(1/2),S^(-1)
c     S^1:g_sp1, S^-1/2:g_sm12, S^+1/2:g_sp12 ,S^-1:g_sm1
c
      call calc_s_powers(g_scr,g_tmp,nbf_ao,toll_s,svals,g_svecs,
     &   g_sp1,g_sm12,g_sp12,g_sm1)
c
c     Vector Guess
c
      spinor = .false. 
      spinor=spinor_guess(movecs_in)
c
      if(.not.spinor)then 
         if (rtdb_get(rtdb,'sodft:fon',mt_log,1,fon)) then 
          if (.not.rtdb_get(rtdb,'sodft:nmo_fon',mt_int,1,nmo_fon)) then
             if (me.eq.0) then
                write(LuOut,*)"Warning: fractional occupation 
     &          calculation specified without nmo_fon"
             end if
          end if
          if (.not.rtdb_get(rtdb,'sodft:nel_fon',mt_int,1,nel_fon)) then
             if (me.eq.0) then
                write(LuOut,*)"Warning: fractional occupation 
     &          calculation specified without nel_fon"  
             end if
          end if
            anoca = noc(1)
            anocb = noc(2)
            noc(1)=(noc(1)+noc(2)-nel_fon)/2
            noc(2)=(noc(1)+noc(2)-nel_fon)/2
         endif
c
         call dft_guessin(movecs_in,ldmix,ncanorg,fon,
     &     vecs_or_dens, ipol,nbf_ao,g_movecs,g_gmovecs,
     &     toll_s,svals)
c
         call scf_vectors_guess(rtdb, tol2e_sleazy, geom, ao_bas_han, 
     &        basis_trans, movecs_in, movecs_out, 
     &        movecs_guess, scftype, nclosed, nopen, 
     &        nbf, nmo, noc(1), noc(2),  k_eval, k_occ, 
     &        k_ir, g_gmovecs, g_dens, 'density', 
     &        'dft', title, oskel, oadapt, 
     &        .true.) 
c
         call dft_guessout(nmo,nbf_ao,g_gmovecs,g_movecs,ipol)
c
         if (rtdb_get(rtdb,'sodft:fon',mt_log,1,fon))then 
            noc(1)=anoca
            noc(2)=anocb
         endif
c     
c     Spinor occupancies
c     
         call dfill(nbf_mo, 0.0d0, dbl_mb(k_occ), 1)
         do i = 1, nTotEl
            dbl_mb(i-1+k_occ) = 1.0d0
         enddo
c     
c     Map initial guess movecs from spin-free calculations g_moso(1) 
c     noc(1).ge.noc(2) is assumed
c         
         do i=1,min(noc(1),noc(2))
            call ga_dadd_patch(1.d0,g_movecs(1),1,nbf_ao,i,i, 
     $           0.d0,g_moso(1),1,nbf_ao,2*(i-1)+1,2*(i-1)+1,
     $           g_moso(1),1,nbf_ao,2*(i-1)+1,2*(i-1)+1) 
            call ga_dadd_patch(1.d0,g_movecs(2),1,nbf_ao,i,i, 
     $           0.d0,g_moso(1),1+nbf_ao,nbf_mo,2*(i-1)+2,2*(i-1)+2,
     $           g_moso(1),1+nbf_ao,nbf_mo,2*(i-1)+2,2*(i-1)+2)
         enddo
         do i=noc(2)+1,noc(1)
            call ga_dadd_patch(1.d0,g_movecs(1),1,nbf_ao,i,i, 
     $           0.d0,g_moso(1),1,nbf_ao,noc(2)+i,noc(2)+i,
     $           g_moso(1),1,nbf_ao,noc(2)+i,noc(2)+i) 
         enddo
         do i=noc(2)+1,noc(1)
            call ga_dadd_patch(1.d0,g_movecs(2),1,nbf_ao,i,i, 
     $           0.d0,g_moso(1),1+nbf_ao,nbf_mo,noc(1)+i,noc(1)+i,
     $           g_moso(1),1+nbf_ao,nbf_mo,noc(1)+i,noc(1)+i) 
         enddo
         do i=noc(1)+1,nbf_ao
            call ga_dadd_patch(1.d0,g_movecs(1),1,nbf_ao,i,i, 
     $           0.d0,g_moso(1),1,nbf_ao,2*(i-1)+1,2*(i-1)+1,
     $           g_moso(1),1,nbf_ao,2*(i-1)+1,2*(i-1)+1) 
            call ga_dadd_patch(1.d0,g_movecs(2),1,nbf_ao,i,i, 
     $           0.d0,g_moso(1),1+nbf_ao,nbf_mo,2*(i-1)+2,2*(i-1)+2,
     $           g_moso(1),1+nbf_ao,nbf_mo,2*(i-1)+2,2*(i-1)+2)
         enddo
      endif  !if not spinor

      if(spinor)then 
c     
c     read spinors from files 
c     
c     get MO vectors from file
c     
      status = movecs_read_header_so(movecs_in, title1, basis_name1,
     $        scftype1, nbf1)
c     
c     Should check much more info than just nbf for consistency
c     
c     get mo eigevectors
c     
         if (2*nbf_ao .ne. nbf1)then
          write(6,*)'ehrenfest_so movecs output = ',movecs_in
          call errquit(pname//'Could not read mo vectors',911,DISK_ERR)
         else 
          status = .true.
          status = status .and.
     $           movecs_read_so(movecs_in, dbl_mb(k_occ),
     $           dbl_mb(k_eval(1)), g_moso)
         endif
c     
         if (.not.status)then
          write(6,*)'ehrenfest_so movecs output = ',movecs_in
          call errquit(pname//'Could not read mo vectors',917,DISK_ERR)
         endif
c     
         call movecs_swap_so(rtdb,'dft',scftype,g_moso,
     &        dbl_mb(k_occ),dbl_mb(k_eval(1)))
      endif  !spinor
c     
c     Form Re and Im of density matrix
c     
      if (.not.rtdb_get(rtdb,'sodft:fon',mt_log,1,fon))
     &    fon = .false.
      if (.not.rtdb_get(rtdb,'sodft:nmo_fon',mt_int,1,nmo_fon))
     &    nmo_fon = 0
      if (.not.rtdb_get(rtdb,'sodft:nel_fon',mt_int,1,nel_fon))
     &    nel_fon = 0
      nTotOcc = (nTotEl-nel_fon) + nmo_fon 
      if(fon)then 
         avg_fon = dble(nel_fon)/dble(nmo_fon) 
         do i = (nTotEl-nel_fon)+1, nTotOcc  
            dbl_mb(i-1+k_occ) = avg_fon 
         enddo
      endif
c     
c     Fractionally occupied mo's are scaled by the sqrt of the fon
c     
      if(nmo_fon .ne. 0 .and. nel_fon .ne. 0 .and. fon)then 
         scale = sqrt(avg_fon)
         call ga_scale_patch(g_moso(1), 
     &        1, nbf_mo, (nTotEl-nel_fon)+1, nTotOcc, scale)
         call ga_scale_patch(g_moso(2), 
     &        1, nbf_mo, (nTotEl-nel_fon)+1, nTotOcc, scale)
      endif
c
c     Calculate the spin-orbit density matrix using the scaled mo's
c 
      call dft_densm_so(g_densso, g_moso, nbf_ao, nTotOcc)
c
c     Restore the scaled mo's
c
      if(nmo_fon .ne. 0 .and. nel_fon .ne. 0 .and. fon)then 
         scale = 1.0d0/scale 
         call ga_scale_patch(g_moso(1), 
     &        1, nbf_mo, (nTotEl-nel_fon)+1, nTotOcc, scale)
         call ga_scale_patch(g_moso(2), 
     &        1, nbf_mo, (nTotEl-nel_fon)+1, nTotOcc, scale)
      endif
c      
c     Calculate the spin-free density matrix from the spin-orbit density matrix
c
      call ga_zero(g_dens(1))
      call ga_zero(g_dens(2))
      !pack real part of dens mat g_densso(1) -> g_dens(1) and (2)
      call ga_dens_sf(g_dens, g_densso(1), nbf_ao)
c
c     Tidy up SCF
c     
      call fock_2e_tidy(rtdb)
c     
c     Set initial coulomb acc
c     
      if (movecs_guess.eq.'restart')ltight=.true.
c     
c     May not want levelshifting initially until sure that the
c     transformed Fock matrix will be diagonally dominant, or
c     alternatively shift the piss out of it.
c     
      if (movecs_guess.eq.'restart')then
         levelshifting = .true.
      else
         levelshifting = .false.
      endif
      iswitc = 0
      if (ltight)then
         itol2e = itol_max
         iAOacc = iAOacc_max
         tol_rho = tol_rho_max
         iswitc = 1
      else
         itol2e = min(itol_min,itol_max)
         iAOacc = min(iAOacc_min,iAOacc_max)
         tol_rho = max(tol_rho_min,tol_rho_max)
      endif
c     
      tol2e = 10.d0**(-itol_max)
c     
c     Restore SCF parameters
c     
      call scf_get_fock_param(rtdb, tol2e)
c
      iter = 1 
c 
cso   Memory allocation for the complex diagonalizer 
      if (.not.MA_Push_Get(MT_DCpl,nbf_mo*nbf_mo,'cpl a',la,ia))
     & call errquit(pname//'Cannot allocate cpl a',0, MA_ERR)
      if (.not.MA_Push_Get(MT_Dbl,nbf_mo,'cpl eval',lw,iw))
     & call errquit(pname//'Cannot allocate cpl eval',0, MA_ERR)
      llwork = max(1, 2*nbf_mo-1)
      if (.not.MA_Push_Get(MT_DCpl,llwork,'cpl work',lwork,iwork))
     & call errquit(pname//'Cannot allocate cpl work',0, MA_ERR)
      if (.not.MA_Push_Get(MT_Dbl,max(1,3*nbf_mo-2),'ws',lrwork,irwork))
     & call errquit(pname//'Cannot allocate ws',0, MA_ERR)
      if (.not.ma_push_get(mt_dbl,nbf_mo,'buff',lbuff,ibuff))
     & call errquit(pname//'push_get failed', 13, MA_ERR)
cso     
c     
c     Set up some sundries before the SCF loop
c
c     Write prep time required
c
      call ga_sync()
      if (me.eq.0.and.oprint)then
         current_cpu = util_cpusec()
         write(LuOut,20)current_cpu
   20    format(2x,' Time prior to 1st pass: ',f8.1)
      endif
c     
c     DFT_SCF timer
c     
      start_wall = util_wallsec()
      start_cpu = util_cpusec()
      dft_time = -start_cpu
c
      if (oprint_time) call dft_tstamp('   Before SCF iter loop. ')
c
      last_time_energy = .false.
c
      idet = 0 
c
c     START SCF LOOP
c
 1000 continue 

      if (me.eq.0 .and. oprint_conv_details)
     &   write(LuOut,124)damping, levelshifting, diising
 124  format(10x,' DAMPING=',l1,' LEVELSHIFTING=',l1,
     &           ' DIISING=',l1)
c
      if (me.eq.0.and.oprint_tol)write(LuOut,3234)itol2e,iAOacc,iXCacc
 3234 format(10x,'itol2e=',i2,' iAOacc=',i2,' iXCacc=',i2)
c
c     Initialize the Density matrix
c
      call ga_zero(g_dens(1))
      call ga_zero(g_dens(2))
c     Pack real part of dens mat g_densso(1) -> g_dens(1) and (2)
      call ga_dens_sf(g_dens, g_densso(1), nbf_ao)
c
c     Initialize and Build the Fock matrix
c
      call ga_zero(g_fockso(1))
      call ga_zero(g_fockso(2))
      call ga_zero(g_fock)
      Ecoul = 0.d0
      Exc(1) = 0.d0
      Exc(2) = 0.d0
      call build_fock_total(rtdb, ao_bas_han, nbf_ao, nbf_mo,
     &    Exc, Ecoul, nExc, iVxc_opt, g_xcinv, g_dens, g_vxc, IOLGC,
     &    g_wght, g_xyz, g_nq, wght_GA, rho_n, dbl_mb(irdens_atom),
     &    int_mb(icetobfr), natoms, ipol, nbf_ao_mxnbf_ce, g_tmp,
     &    g_fockso, g_densso, g_so, g_fock, oskel, Etnew)
c
      if (oprint_time)
     &     call dft_tstamp(' End of parallel region. ')
c
c     Last calculation of the energy
c
      if (last_time_energy) goto 2000
c     
      delta = abs(etold-etnew)
c     
      call ga_sync
      rms(1) = 0.d0
      rms(2) = 0.d0
      homo_lumo_gap = 200.0d0
c     
c     Symmetrize the Fock matrix
c
      if (oskel)
     &   call sym_symmetrize(geom, AO_bas_han, .false., g_fock)
c
c     DIIS step taken here.
c     
      if (diising)then
        call diis_driver_so(toll_s, derr, icall, nfock, 
     &           nbf_mo, g_fockso, g_densso, 
     &           g_svecs, svals, diising, nodiis)
        derr(2)=derr(1)
      endif
c     
c     Transform Fock matrix into the orthonormal basis
c     
      call ga_zero(g_tmp2)
      call trans_fock_to_ortho(g_tmp2,nbf_mo,g_sm12,g_fockso)
c     
c     Save the old vectors 
c     
      call ga_copy(g_moso(1), g_old(1))
      call ga_copy(g_moso(2), g_old(2))
c
c     Apply levelshifting if needed
c     
      if (levelshifting)then
        call levelshift_fock(nbf_mo,ntotocc,g_tmp2,g_sp12,g_moso,
     &   g_scr2,g_fockso,rlshift)
      else
        rlshift = 0.d0
      end if
c     
c     Solve for the eigenvalues and eigenvectors of the Hamiltonian.
c
      call diag_fock(nbf_mo,ia,g_fockso,ibuff,g_moso,iwork,irwork,
     &     k_eval,trace,llwork,info)
c
c     Back-transform eigenvectors with S^-1/2 (g_sm12)
c     
      call trans_vec_to_ao(nbf_mo,g_sm12,g_moso,g_fockso)
c
c     Keep orbital ordering according to principle
c     of maximum overlap with previous iteration.
c
      call ga_symmetrize(g_fock)
      if (lmaxov)
     .     call dft_mxovl(ao_bas_han, dbl_mb(k_eval(ispin)),
     &                  g_tmp, g_movecs(ispin), g_s,g_fock,
     .        noc,oprint_conv_details,homo,lumo)
c        
c     Determine the HOMO-LUMO gap 
c
      call calc_homolumogap(k_eval,nTotEl,rlshift,
     &  homo,lumo,homo_lumo_gap)
      if (me.eq.0 .and. oprint_conv_details) then
        write(LuOut,4224)homo,lumo,rlshift, homo_lumo_gap
      end if
c     
      call ga_sync
c     
c     Save previous density for convergence check.
c     
      call ga_copy(g_dens(1), g_movecs(1))
      call ga_copy(g_dens(2), g_movecs(2))
c
c     Symmetry adapt vectors?
c
      if (oadapt)then
         call scf_sym_adapt_so(ao_bas_han, g_moso,
     &        oprint_syma, 2*nbf_ao, name,
     &        .true., 
     &        int_mb(k_ir))
      endif      
c
c     Save the old density matrix for damping 
c
      call ga_copy(g_densso(1), g_damp_so(1))
      call ga_copy(g_densso(2), g_damp_so(2))
c     
c     Form a new density matrix.
c     
      call ga_sync 
      call ga_zero(g_densso(1))
      call ga_zero(g_densso(2)) 
c     
c     Fractionally occupied mo's are scaled by the sqrt of the fon
c     
      if(nmo_fon .ne. 0 .and. nel_fon .ne. 0 .and. fon)then 
         scale = sqrt(avg_fon)
         call ga_scale_patch(g_moso(1), 
     &        1, nbf_mo, (nTotEl-nel_fon)+1, nTotOcc, scale)
         call ga_scale_patch(g_moso(2), 
     &        1, nbf_mo, (nTotEl-nel_fon)+1, nTotOcc, scale)
      endif
      call dft_densm_so(g_densso,g_moso,nbf_ao,nTotOcc) 
c     
c     Restore the scaled mo's
c     
      if(nmo_fon .ne. 0 .and. nel_fon .ne. 0 .and. fon)then 
         scale = 1.0d0/scale 
         call ga_scale_patch(g_moso(1), 
     &        1, nbf_mo, (nTotEl-nel_fon)+1, nTotOcc, scale)
         call ga_scale_patch(g_moso(2), 
     &        1, nbf_mo, (nTotEl-nel_fon)+1, nTotOcc, scale)
      endif
c
c     g_dens(1)=Re(Daa)+Re(Dbb) and g_dens(2)=Re(Dbb)
c     For coulomb and xc potentials only the alpha, Re(Daa), and  
c     beta, Re(Dbb) densities are needed 
c
      call ga_zero(g_dens(1))
      call ga_zero(g_dens(2))
      !Pack real part of dens mat g_densso(1) -> g_dens(1) and (2)
      call ga_dens_sf(g_dens, g_densso(1), nbf_ao)
      call ga_sync
c     
c     Check convergence on Density.
c     
      rms(1) = dft_dencvg(g_dens(1), g_movecs(1), nbf_ao)
      rms(2) = dft_dencvg(g_dens(2), g_movecs(2), nbf_ao)
      call ga_sync
c     
      if (oprint_conv.and.iter.eq.1.and.me.eq.0)then
         nheap = MA_Inquire_Heap(MT_Dbl)
         nstack = MA_Inquire_Stack(MT_Dbl)
         write(LuOut,21)
         write(LuOut,'(10x,a,f10.2,i20)')
     &        ' Heap Space remaining (MW):  ',dble(nheap)*1.D-06,nheap
         write(LuOut,'(10x,a,f10.2,i20)')
     &        'Stack Space remaining (MW):  ',dble(nstack)*1.D-06,nstack
         call util_flush(LuOut)
         write(LuOut,1)
      endif
 21   format(/,10x,' Memory utilization after 1st SCF pass: ')
    1 format(/,
     &     1x,'  convergence    iter        energy       DeltaE   ',
     &     'RMS-Dens  Diis-err    time'/
     &     1x,'---------------- ----- ----------------- --------- ',
     &     '--------- ---------  ------')
      if (oprint_conv.and.me.eq.0)then
         current_cpu = util_cpusec()
         if (diising)then
            write(LuOut,2)ndamp,rlshift,
     &           iter, Etnew+Enuc,
     &           -etold+etnew,sqrt(rms(1)),derr(1),current_cpu
            if (ipol.eq.2)write(LuOut,3)sqrt(rms(2)),derr(2)
         else
            write(LuOut,22)ndamp,rlshift,
     &           iter, Etnew+Enuc,
     &           -etold+etnew,sqrt(rms(1)), current_cpu
            if (ipol.eq.2)write(LuOut,23)sqrt(rms(2))
         endif
         call util_flush(LuOut)
      endif
    2 format(1x,'d=',i2,',ls=',f3.1,',diis',1x,i5,f18.10,
     &     1p,3d10.2,0p,f8.1)
    3 format(51x,1p,2d10.2)
 22   format(1x,'d=',i2,',ls=',f3.1,6x,i5,f18.10,
     &     1p,2d10.2,10x,0p,f8.1)
 23   format(51x,1p,1d10.2)
c
      call ga_sync
c     
c     Save eigenvectors to movecs file
c     
      if (.not.movecs_write_so
     $     (rtdb, ao_bas_han, movecs_out, 'sodft', title,
     &     nbf_mo, dbl_mb(k_occ), dbl_mb(k_eval(1)), g_moso))
     &     call errquit(pname//'movec_write failed', 0, DISK_ERR)
c     
      call ga_sync
c     
      if (me .eq. 0.and.oprint_eval)then
         if (util_print('intermediate evals', print_default))then
            call util_print_centered(LuOut,'eigenvalues',
     &           20,.true.)
            call output(dbl_mb(k_eval(1)), 1, min(nTotEl+10,nbf_mo),
     &           1, 1, nbf_mo, 1, 1)
            call util_flush(6)
         endif
      endif
      if (oprint_vecs)then
         if (me .eq. 0)then
            write(LuOut,*)
            call util_print_centered(LuOut,
     &           'Intermediate MO vectors',40,.true.)
            write(LuOut,*)
            call util_flush(LuOut)
         endif
      endif
c     
c     If open shell compute overlap of alpha orbitals with beta 
c     orbitals.
c     
      if ((ipol.gt.1).and.(oprint_interm_overlap)) then
         call dft_mxspin_ovlp(nbf_ao,nmo,ao_bas_han, noc, g_movecs(1), 
     &       g_movecs(2),  g_tmp)
      endif
c     
c     <S2> for open shell
c     
      if ((ipol.gt.1).and.(oprint_interm_S2)) then
         call dft_s2_value(geom, AO_bas_han, .false., noc(1), noc(2),
     &        nbf_ao, g_dens(1), g_dens(2))
      endif
c     
      call ga_sync
c     
c     Check for SCF convergence.
c     
      call ga_sync
      call dft_scfcvg(rms, derr, Etold, Etnew, 
     &     e_conv, d_conv, g_conv, ipol, 
     &     iter, iterations, idone, rtdb,
     &     converged, diising)
      if (delta.lt.1.d-3)then
c     
c     Set coulomb acc to max (e.g., input parameter).
c     (note, may also require re-initializing DIIS)
c     
         itol2e = itol_max
         iAOacc = iAOacc_max
         tol_rho = tol_rho_max
         iswitc = iswitc+1
      endif
c     
c     Damping implemented here.
c    
      if (damping)then
         pp = dble(ndamp)*1.d-2
         onempp = 1.0d0 - pp
         call ga_dadd(pp, g_damp_so(1), onempp, g_densso(1), 
     &        g_densso(1))
         call ga_dadd(pp, g_damp_so(2), onempp, g_densso(2), 
     &        g_densso(2))
         call ga_zero(g_dens(1))
         call ga_zero(g_dens(2))
         call ga_dens_sf(g_dens, g_densso(1), nbf_ao)
         call ga_sync
      else
         ndamp = 0
      endif
      call ga_sync
      iter = iter + 1    ! update the iteration number
c     
c     Check convergence parameters.
c     
      if ((delta.lt.dampon.and.delta.gt.dampoff).or.iter.le.ncydp)then
         damping = .true.
         ndamp = ndamp_input
      else
         damping = .false.
         ndamp = ndamp_def
      endif
c     
      if ((delta.lt.levlon.and.delta.gt.levloff).or.
     &     (iter.le.ncysh))then
         if (homo_lumo_gap.lt.hl_tol)then
            levelshifting = .true.
            rlshift = rlshift_input
            if (check_shift)then
               if (lumo .lt. homo)then
                  levelshifting = .false.
                  if (me.eq.0 .and. oprint_conv_details)
     &                 write(LuOut,2224)homo, lumo
               endif
            endif
         else
            levelshifting = .false.
            rlshift = rlshift_def
         endif 
      else
         levelshifting = .false.
         rlshift = rlshift_def
      endif
c
 2224 format(10x,'HOMO = ',f6.2,' LUMO (with shift) = ',f6.2,
     &     /,10x,'Unshifted LUMO is less than HOMO.',
     &     /,10x,'Turning levelshifting OFF this iteration.')
c     
      if ((delta.lt.diison.and.delta.gt.diisoff).or.
     &     iter.le.ncyds.or.keep_diis_on)then
         diising = .true.
c     
c     Once started, keep DIIS on until diisoff threshold.
c     
         keep_diis_on = .true.
      else
         diising = .false.
      endif
      if (delta.lt.diisoff.or.(ncyds.gt.0.and.iter.gt.ncyds))then
         diising = .false.
         keep_diis_on = .false.
      endif
c     
      if (nodamping)damping = .false.
      if (nolevelshifting) then 
         levelshifting = .false.
         rlshift=rlshift_def
      endif        
      if (nodiis)diising = .false.
c     
      Etold = Etnew
c     
      lmaxov = lmaxov_sv
      if ((lumo - homo).lt.-hl_tol.and.lmaxov)then
         lmaxov = .false.
         if (me.eq.0 .and. oprint_conv_details)
     &        write(LuOut,224)homo, lumo
c
 224     format(10x,' HOMO = ',f6.2,' LUMO = ',f6.2,
     &        /,10x,'Significant orbital reordering with',
     &        ' maximum overlap',
     &        /,10x,'turned ON.  Turning max_ovl OFF.')
c
      endif
c
      if (oprint_energy_step.and.me.eq.0)then         
       current_cpu = util_cpusec()
       call print_energies(etnew,enuc,ecore,ecoul,exc,nexc,
     &             rho_n,current_cpu)
      endif
c     
c     Check for remaining time to exit "gracefully"
c     
      current_wall = util_wallsec()
      if ((iter-1).gt.1)then
         elapsed_wall = current_wall - save_wall
         save_wall = current_wall
      else
         elapsed_wall = current_wall - start_wall
         save_wall = current_wall
      endif
c     
      if (converged)then
c     
c     If converged probably need a few seconds to clean things up 
c     and calculate a few properties.
c     
         wall_time_reqd = 5.0
c     
      else
c     
c     If not converged probably need at least the amount time
c     required for previous iteration (multiply by 1.2 to be on the safe side).
c     
         wall_time_reqd = elapsed_wall*1.2d0
      endif
      int_wall_time_reqd = wall_time_reqd
      if (.not.util_test_time_remaining(rtdb, int_wall_time_reqd))then
         if (me.eq.0)then
            write(LuOut,*)
            call util_print_centered(LuOut,
     &           'Exiting due to time limitations.', 20, .true.)
            write(LuOut,*)
            goto 2000
         endif
      endif
      if (idone.eq.0.or.(iswitc.lt.2.and.iter.lt.iterations))
     &     go to 1000           ! BEGIN NEW SCF ITERATION
      if (idone.eq.1.and.(.not.last_time_energy))then
         last_time_energy = .true.
         go to 1000             ! BUILD FINAL TOTAL ENERGIES
      endif
c     
 2000 continue
c    
c     Print energies and some information
c
      if (me.eq.0.and.oprint)then
         if (.not.converged)then
            write(LuOut,*)
            call util_print_centered(LuOut,
     &           'Calculation failed to converge', 20, .true.)
            write(LuOut,*)
         endif
         dft_time = dft_time+util_cpusec()
         call print_energies(etnew,enuc,ecore,ecoul,exc,enuc,nexc,
     &                      rho_n,dft_time)
         call util_flush(LuOut)
      endif
c
c     Propagate the density matrix in time
c
      ! TODO: need ntime, tmin as part of the input
      tmin = 0.d0
      tmax = 0.06 !200.0
      dt = 0.02d0
      ntime = max (nint (tmax/dt), 1)
      rk4 = .true.
      do itime = 1,ntime
c
c        Increment the time step
         tt = tmin + (itime-1)*dt 
c
c        Add call to FDTD driver here to get the electric field
c
c        Initialize the Density matrix
         call ga_zero(g_dens(1))
         call ga_zero(g_dens(2))
c     
c        Calculate the instantaneous dipole moment
c        Dipole moment depends on the symmetric part of the density matrix
         call ga_dens_sf(g_dens, g_densso(1), nbf_ao) ! g_densso(1) -> g_dens(1) and (2)
         call calc_dipole(rtdb, ao_bas_han, ipol, g_dens, dipole, tt)
         if (me.eq.0)  then
            write(luout,*) "Dipole:",tt,dipole(1),dipole(2),dipole(3)
         end if
c
c        Perform an instantaneous Mulliken population analysis
         !call calc_mullpop(rtdb,geom,ao_bas_han,g_s,g_dens) ! no need for now
c
c        Build the Fock matrix
         call ga_zero(g_fockso(1))
         call ga_zero(g_fockso(2))
         call ga_zero(g_fock)
         call build_fock_total(rtdb, ao_bas_han, nbf_ao, nbf_mo,
     &    Exc, Ecoul, nExc, iVxc_opt, g_xcinv, g_dens, g_vxc, IOLGC,
     &    g_wght, g_xyz, g_nq, wght_GA, rho_n, dbl_mb(irdens_atom), 
     &    int_mb(icetobfr), natoms, ipol, nbf_ao_mxnbf_ce, g_tmp, 
     &    g_fockso, g_densso, g_so, g_fock, oskel, Etnew)
c          
c        Add dipole contribution to the Fock matrix
         Edip = 0.d0
         call calc_ddotE(tt, g_dipx, g_dipy, g_dipz, g_ddotE)
         !call ga_print(g_ddotE)
         call energy_dipole(g_dens, g_ddotE, Edip)
         call ga_fock_sf(g_ddotE, g_fockso(1), nbf_ao)
c
c        Total energy
         Etotal = Etnew + Edip + Enuc
         if (me.eq.0) then
            write(luout,*) "Etotal:",tt, Etotal   ! print total energy at time tt
         end if
c
c        Do the propagation
         if (rk4) then ! 4th order Runge-Kutta 
           call prop_densmat_rk4(nbf_mo,g_sm1,g_fockso,g_densso,
     &            g_work,g_k1,g_k2,g_k3,g_k4,g_kwork,tt,dt)
         else ! Euler
           call prop_densmat_euler(nbf_mo,g_sm1,g_fockso,g_densso,
     &            g_work,g_dprdt,g_dpidt,tt,dt)
         end if
c
c        Calculate the instantaneous current
c        Current = Sum(mu,nu) A(mu,nu) < chi_mu | nabla | chi_nu >
c        where A(mu,nu) is the anti-symmetric part of the density matrix
         call ga_dens_sf(g_dens, g_densso(2), nbf_ao) ! g_densso(2) -> g_dens(1) and (2)
         call calc_current(ao_bas_han, nbf_ao, ipol, g_dens,
     &           g_vel, g_temp, jcurr, tt)
         if (me.eq.0) then
            write(luout,*) "Current:", tt,jcurr(1),jcurr(2),jcurr(3)
         end if
c
      end do   ! time propagation loop
c
c     Prepare for clean up and close down
c
      if(fon)then 
      if (.not.ma_pop_stack(kfon_occ))
     & call errquit('ehrenfest_so: cannot pop stack:lfon_occ',0, MA_ERR)
      endif
      if (.not. ga_destroy(g_damp_so(1))) 
     & call errquit(pname//'Could not destroy g_damp_so', 0, GA_ERR)
      if (.not. ga_destroy(g_damp_so(2))) 
     & call errquit(pname//'Could not destroy g_damp_so', 0, GA_ERR)
c
      if (util_print('final vectors analysis', print_default)) then
         do ilo = 1,max(1,nTotEl-10)
            if (dbl_mb(k_eval(1)+ilo-1) .ge. eval_pr_tol_lo) 
     &           goto 961
         enddo
 961     do ihi = min(nTotEl+10,nbf_mo), nbf_mo
            if (dbl_mb(k_eval(1)+ihi-1) .ge. eval_pr_tol_hi) 
     &           goto 9611
         enddo
         ihi = max(ihi-1,1)
 9611    continue
         if (util_print('final vectors analysis', print_high)) then
            ilo = 1
            ihi = nbf_mo
         endif
         call movecs_anal_so(ao_bas_han, ilo, ihi, 0.15d0, 
     &        g_moso, 
     &        'DFT Final Molecular Orbital Analysis', 
     &        .true., dbl_mb(k_eval(1)), oadapt, 
     &        int_mb(k_ir), .true., dbl_mb(k_occ))
      endif
c     
c     Store energy and convergence status ... must store before
c     write movecs since date of insertion is used.
c     
      if (.not. rtdb_put(rtdb, 'sodft:energy',MT_DBL,1,(Etnew+Enuc)))
     & call errquit(pname//'Failed to store energy in rtdb',0,RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'sodft:converged',MT_LOG,1,converged))
     & call errquit(pname//'Failed to store converged in rtdb',0,
     &       RTDB_ERR)
c
c     Output energies and eigenvectors to disk
c     
      if (.not.movecs_write_so
     $     (rtdb, ao_bas_han, movecs_out, 'sodft', title,
     &     nbf_mo, dbl_mb(k_occ), dbl_mb(k_eval(1)), g_moso))
     &     call errquit(pname//'movec_write failed', 0, DISK_ERR)
      call ga_sync()
c     
c     Shut down DIIS.
c     
      if (icall(1).gt.0)then
         icall(1) = -1
         call diis_driver_so(toll_s, derr, icall, nfock, 
     &        nbf_mo, g_fockso, g_densso, 
     &        g_svecs, svals, diising, nodiis)
      endif
c     
c     If open shell compute overlap of alpha orbitals with beta orbitals.
c     
      if (ipol.gt.1)then
         call dft_mxspin_ovlp(nbf_ao,nmo,ao_bas_han,noc,g_movecs(1), 
     &        g_movecs(2), g_tmp)
      endif
c
      if (wght_GA)then
         if (.not. ga_destroy(g_wght)) 
     &  call errquit(pname//'Could not destroy g_wght', 0, GA_ERR)
         if (.not. ga_destroy(g_xyz)) 
     &  call errquit(pname//'Could not destroy g_xyz', 0, GA_ERR)
         if (.not. ga_destroy(g_nq)) 
     &  call errquit(pname//'Could not destroy g_nq', 0, GA_ERR)
      endif
c     
c     == Calculate <S2> for open shell ==
c     
      if (ipol.gt.1)
     &   call dft_s2_value(geom,AO_bas_han,.false.,noc(1),noc(2),
     &        nbf_ao,g_dens(1),g_dens(2))
c     
c     == Calculate moments ==
c
      if (natoms .gt. 1)
     &   call dft_mpole(rtdb, ao_bas_han, ipol, g_dens(1), g_dens(2))
c     
c     == Print ==
c     
      if (util_print('schwarz',print_high).and.(.not.CDFIT))
     & call schwarz_print(natoms, nshells_ao)
c     
      if (me .eq. 0)then
cng         if (util_print('final evals', print_high))then
            call util_print_centered(LuOut,'Final alpha eigenvalues',
     &           20,.false.)
            call util_print_centered(LuOut,
     &           '(all occupied plus 10 virtual)',20,.true.)
            call output(dbl_mb(k_eval(1)),
     &           1, min(noc(1)+10,nbf_ao),
     &           1, 1, nbf_ao, 1, 1)
            if (ipol.eq.2)then
               write(LuOut,*)
               call util_print_centered(LuOut,'Final beta eigenvalues',
     &              20,.true.)
               call util_print_centered(LuOut,
     &              '(all occupied plus 10 virtual)',20,.true.)
               call output(dbl_mb(k_eval(2)),
     &              1, min(noc(1)+10,nbf_ao),
     &              1, 1, nbf_ao, 1, 1)
            endif
            call util_flush(6)
cng         endif
         if (oprint_final_vecs)then
            write(LuOut,*)
            call util_print_centered(
     &           LuOut,'Final MO vectors',40,.true.)
            write(LuOut,*)
            call util_flush(LuOut)
         endif
      endif
c
c     Start deallocate
c
      if (.not.ma_pop_stack(lbuff)) then
      call ma_summarize_allocated_blocks()
      call util_flush(6)
      call util_flush(0)
      call errquit(pname//'Cannot pop stack:lbuff',lbuff, MA_ERR)
      endif
      if (.not.ma_pop_stack(lrwork))
     &   call errquit(pname//'Cannot pop stacklrwork',0, MA_ERR)
      if (.not.ma_pop_stack(lwork))
     &   call errquit(pname//'Cannot pop stack:lwork',0, MA_ERR)
      if (.not.ma_pop_stack(lw))
     &   call errquit(pname//'Cannot pop stack:lw',0, MA_ERR)
      if (.not.ma_pop_stack(la))
     &   call errquit(pname//'Cannot pop stack:la',0, MA_ERR)
      if (.not.ma_pop_stack(l_ir))
     &   call errquit(pname//'Cannot pop stack:l_ir',0, MA_ERR)
c     
      if (ipol.gt.1)then
       if (.not. ga_destroy(g_fockt)) 
     &   call errquit(pname//'Could not destroy g_fockt', 0, GA_ERR)
      endif
      if (.not. ga_destroy(g_tmp)) 
     &  call errquit(pname//'Could not destroy g_tmp', 0, GA_ERR)
c
      call fock_2e_tidy(rtdb)
c     
      if (.not.ma_pop_stack(l_eval))
     & call errquit(pname//'Cannot pop stack:l_eval',0, MA_ERR)
      if (.not.ma_pop_stack(l_occ))
     & call errquit(pname//'Cannot pop stack:l_occ',0, MA_ERR)
      if (.not.ma_pop_stack(lrdens_atom))
     & call errquit(pname//'Cannot pop stack:lrdens_atom',0, MA_ERR)
      if (.not.ma_pop_stack(lcetobfr))
     & call errquit(pname//'Cannot pop stack:lcetobfr',0, MA_ERR)
      if (.not.ma_pop_stack(lcntobfr))
     & call errquit(pname//'Cannot pop stack:lcntobfr',0, MA_ERR)
      if (.not.ma_pop_stack(lcntoce))
     & call errquit(pname//'Cannot pop stack:lcntoce',0, MA_ERR)
c
      if(.not.ga_destroy(g_moso(1)))     
     & call errquit(pname//'Error destroy Movecs Re',0, GA_ERR)
      if(.not.ga_destroy(g_moso(2)))     
     & call errquit(pname//'Error destroy Movecs Im',0, GA_ERR)
      if(.not.ga_destroy(g_fockso(1)))     
     & call errquit(pname//'Error destroy Fock Re',0, GA_ERR)
      if(.not.ga_destroy(g_fockso(2)))     
     & call errquit(pname//'Error destroy Fock Im',0, GA_ERR)
      if(.not.ga_destroy(g_densso(1)))     
     & call errquit(pname//'Error destroy DenMx Re',0, GA_ERR)
      if(.not.ga_destroy(g_densso(2)))     
     & call errquit(pname//'Error destroy DenMx Im',0, GA_ERR)
      if(.not.ga_destroy(g_tmp2))     
     & call errquit(pname//'Error destroy old re',0, GA_ERR)
      if(.not.ga_destroy(g_old(1)))     
     & call errquit(pname//'Error destroy old im',0, GA_ERR)
      if(.not.ga_destroy(g_old(2)))     
     & call errquit(pname//'Error destroy Tmp ReIm',0, GA_ERR)
      if(.not.ga_destroy(g_so(1)))     
     & call errquit(pname//'Error destroy so z',0, GA_ERR)
      if(.not.ga_destroy(g_so(2)))     
     & call errquit(pname//'Error destroy so y',0, GA_ERR)
      if(.not.ga_destroy(g_so(3)))     
     & call errquit(pname//'Error destroy so x',0, GA_ERR)
      if (.not. ga_destroy(g_scr)) 
     & call errquit(pname//'Could not destroy g_scr', 0, GA_ERR)
      if (.not. ga_destroy(g_scr2)) 
     & call errquit(pname//'Could not destroy g_scr2', 0, GA_ERR)
c
c     == Time propagation arrays ==
c
      if(.not.ga_destroy(g_sp1))     
     & call errquit(pname//'Error destroy g_sp1',0, GA_ERR)
      if(.not.ga_destroy(g_sm1))     
     & call errquit(pname//'Error destroy g_sm1',0, GA_ERR)
      if(.not.ga_destroy(g_sp12))     
     & call errquit(pname//'Error destroy g_sp12',0, GA_ERR)
      if(.not.ga_destroy(g_sm12))     
     & call errquit(pname//'Error destroy g_sm12',0, GA_ERR)
c
      if(.not.ga_destroy(g_temp))     
     & call errquit(pname//'Error destroy g_temp',0, GA_ERR)
c
      do iwork = 1,nwork
        if(.not.ga_destroy(g_work(iwork)))     
     &   call errquit(pname//'Error destroy g_work',0, GA_ERR)
        if(.not.ga_destroy(g_k1(iwork)))     
     &   call errquit(pname//'Error destroy g_work',0, GA_ERR)
        if(.not.ga_destroy(g_k2(iwork)))     
     &   call errquit(pname//'Error destroy g_work',0, GA_ERR)
        if(.not.ga_destroy(g_k3(iwork)))     
     &   call errquit(pname//'Error destroy g_work',0, GA_ERR)
        if(.not.ga_destroy(g_k4(iwork)))     
     &   call errquit(pname//'Error destroy g_work',0, GA_ERR)
        if(.not.ga_destroy(g_kwork(iwork)))     
     &   call errquit(pname//'Error destroy g_work',0, GA_ERR)
      end do
c
      if(.not.ga_destroy(g_dprdt))     
     & call errquit(pname//'Error destroy g_dprdt',0, GA_ERR)
      if(.not.ga_destroy(g_dpidt))     
     & call errquit(pname//'Error destroy g_dpidt',0, GA_ERR)
c
      ehrenfest_sodft = converged
c
 4224 format(10x,' HOMO = ',f6.2,' LUMO = ',f6.2,
     &           ' RLSHIFT = ',f6.2,' GAP = ',f6.2)
c
      return
      end
