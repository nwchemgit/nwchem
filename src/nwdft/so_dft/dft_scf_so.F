c
c    Main spin-orbit DFT driver
c
      logical function dft_scf_so
     &                 (rtdb, Etold, Enuc, iVcoul_opt, iVxc_opt, 
     &                  iter, g_dens, g_dens_at, g_movecs, g_vxc, 
     &                  g_fock, g_svecs, svals, g_xcinv, g_s)
c     
c     $Id$
c     
      implicit none
c
#include "errquit.fh"
#include "bas.fh"
#include "geom.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "cdft.fh"
#include "global.fh"
#include "msgids.fh"
#include "util.fh"
#include "cosmo.fh"
#include "zora.fh" 
#include "case.fh"
#include "modelpotential.fh"
#ifdef SCALAPACK
#include "dft_scaleig.fh"
#endif
c     
      integer rtdb              ! [input]
      double precision Etold, Enuc, trace
      integer iVcoul_opt
      integer iVxc_opt
      integer iter, swap(20),nswap, ndet, idet
      integer g_dens(2), g_movecs(2), g_vxc(4), 
     &     g_fock, g_svecs,
     &     g_xcinv,  g_scr
      integer g_dens_at(2)
      double precision  svals(*)
c     so
      integer g_densso(2), g_tmp_ri,   
     &     g_moso(2), g_old(2),
     &     g_fockso(2), g_scr2, g_damp_so(2), g_gmovecs(2)
c
      integer la, ia            ! complex*16 a(nbf_ao, nbf_ao)
      integer lw, iw            ! double precision w(nbf_ao)
      integer llwork 
      integer lwork, iwork      ! complex*16 work(3) 
      integer lrwork, irwork    ! double precision rwork
      integer info 
c      integer lbuff, ibuff
      integer nbf_mo 
      integer g_s
c
c     declarations for fractional occupation 
c
      integer nmo_fon,nmo_fon2  ! number of fractionally occupied orbitals
      integer ncore_fon ! number of fully occupied orbitals
      integer nfilled, ilow, ihigh
      double precision nel_fon,nel_fon2  ! fractional electron number
      double precision avg_fon,avg_fon2  ! fractional occupancy (averaged)
      integer nTotOcc           ! nTotOcc: no. of (occupied maybe f.o.) mo's 
      logical fon 
      integer ntmp_fon(4)
      double precision rtmp_fon(4), pstrace
      double precision scale, ncheck
      integer kfon_occ, lfon_occ
      logical debug_fon, det_eng 
      integer i_swap,j_swap
      logical noscf
c
c     so
      double precision rho_n, toll_s
c     
c     zora related
      double precision ener_scal, ener_kin
      double precision numelecs
      integer ncanorg
      logical ldmix
      character*7 vecs_or_dens
      integer icalczora
      logical do_purescalar
      logical KNucl
      logical ofinite
      integer l_zetanuc, k_zetanuc
c
c     virtual spectrum
      logical lvspec
      integer owstart(2), owend(2)
      character*255 vspecfilename
c
c     vdw contrib
      double precision dum
      logical xc_chkdispauto
      external xc_chkdispauto
      logical disp
c
c     misc
      double precision tanalyze
c
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
      double precision xc_hfexch
      external xc_hfexch
c
      logical oprint, status
c
      double precision Edisp
      double precision Exc(2), rms(2), derr(2)
      double precision Exc2(2)
      double precision exch0,exch0r,exch0i
      double precision exchZ,exchZr,exchZi
      double precision exchX,exchXr,exchXi
      double precision exchY,exchYr,exchYi
      double precision edftxc(2),edftxcs,edftxcd
      double precision exchtot
c
      integer nmo(2), icall(2)
c     integer nva
      integer n3c_dbl, n3c_int, n_batch
      integer iwhat_max
      integer l_3cwhat, k_3cwhat, l_3cERI, k_3cERI
      integer dft_n3cint, n_semi_bufs, fd
      external dft_n3cint
      double precision dft_n3cdbl
      external dft_n3cdbl
      Integer l_eval
      integer k_eval(2)
      integer iatom,natoms, nTotEl
c
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
      integer leneval
      integer lcntoce, icntoce, lcntobfr, icntobfr,
     &     lcetobfr, icetobfr, lrdens_atom, irdens_atom,
     &     nscr, lscr, iscr
      double precision start_wall, current_wall, elapsed_wall,
     &     save_wall, current_cpu, start_cpu,
     &     wall_time_reqd
      integer int_wall_time_reqd
      double precision ecoul, ecore, noso
      double precision pp, delta
      double precision anucl_charg, anel
      double precision anoca, anocb, onempp
      double precision etnew, tol2e, tol2e_sleazy
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
      double precision zero, onem, one, mone
      parameter(zero = 0.d0, one = 1.d0, mone=-1.0d0, onem = -one)
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
c     zora related
      logical dft_zora_read_so, dft_zora_write_so
      external dft_zora_read_so, dft_zora_write_so
c
      logical dft_zora_inquire_file_so
      external dft_zora_inquire_file_so
c
      character*255 zorafilename
      integer g_zora_sf(2)
      integer g_zora_scale_sf(2)
      integer g_zora_so(3)
      integer g_zora_scale_so(3)
      double precision Ezora_sf
      integer switch_focc
      logical do_avg_fon
      logical do_core_fon
c
      logical do_zora_recalc, zora_recalc
      external zora_recalc
c
      integer nmoout
      logical spinor
      logical dft_mem3c
      external dft_mem3c
c
c     ----- cosmo -----
c
      logical odbug
      logical ocosmo_got_gasphase
      double precision egas
      double precision esol
      double precision ecosmo ! intra energy between cosmo charges
      logical cosmo_esp
      character*255 cosmo_file
c
      double precision gcds, gspol, espol, gspoldyn, gspolneq
      double precision gstote, estote
      integer do_cosmo_vem, istep_cosmo_vem
      logical do_cosmo_smd
c
c     Preliminaries
      me = ga_nodeid()
      nproc = ga_nnodes()
c
      npol = (ipol*(ipol+1))/2
      itol_max = itol2e
      iaoacc_max = iaoacc
      tol_rho_max = tol_rho
c
      call ecce_print_module_entry('dft')
      dft_scf_so = .false.
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
      odbug = util_print('cosmo', print_high)
c
c     Solvation related
      gcds = 0.0
      gspol = 0.0
      gstote = 0.0
      do_cosmo_smd = .false.
      if (cosmo_on) then
       if (.not. rtdb_get
     $ (rtdb,'cosmo:do_cosmo_smd',mt_log,1,do_cosmo_smd))
     $ call errquit('
     $ dft_scf: cannot get do_cosmo_smd from rtdb',
     $ 0,rtdb_err)
       if (do_cosmo_smd) then
         if (.not. rtdb_get
     $   (rtdb,'cosmo:gcds',mt_dbl,1,gcds))
     $    call errquit(
     $    'dft_scf: cannot get gcds from rtdb',
     $    0,rtdb_err)
       endif
      endif
c
c     Scalar relativistic calculations within the two-component framework
      do_purescalar = .false.
      if (.not.rtdb_get(rtdb,'sodft:scalar',mt_log,1,do_purescalar))
     &     do_purescalar = .false.
c
      if (me.eq.0.and.do_purescalar) then
         call util_print_centered(LuOut,
     $        'Neglecting spin-orbit terms', 23, .true.)
         write(LuOut,*)
       endif ! me
c
c     do average fractional occupation by default
      do_avg_fon = .true.
      if (.not.rtdb_get(rtdb,'dft:avg_fon',mt_log,1,do_avg_fon))
     &     do_avg_fon = .true.
c
c     read MO analysis threshold
      tanalyze = 0.025d0
      if (.not.rtdb_get(rtdb,'movecs:tanalyze',MT_DBL,1,tanalyze))
     &   tanalyze = 0.025d0
c
c     do average occupation starting with core orbitals
      do_core_fon = .false.
      if (.not.rtdb_get(rtdb,'dft:core_fon',mt_log,1,do_core_fon))
     &     do_core_fon = .false.
c
      ispin=1
      call int_1e_uncache_ga()      
c     !!! BGJ
c     Store SCF hamiltonian type as DFT for use in BGJ routines
      if (.not. rtdb_put(rtdb, 'bgj:scf_type', MT_INT, 1, 2))
     $     call errquit('dft_scf_so: put of bgj:scf_type failed',0,
     &       RTDB_ERR)
c     !!! BGJ
c     
c     see if levelshifting monitoring is desired
      if (.not. rtdb_get(rtdb, 'dft:check_shift', mt_log, 1,
     &     check_shift))then
         check_shift = .false.      
      endif
c     
      if (.not. geom_ncent(geom, natoms))
     &     call errquit('dft_scf_so: geom_ncent failed',73, GEOM_ERR)
      if (.not. geom_nuc_charge(geom, anucl_charg))
     &     call errquit('dft_scf_so: geom_nuc_charge failed', 
     & 0, GEOM_ERR)
c     
      anel = int(anucl_charg) - rcharge
c
c     Determine whether to fit the electronic charge density.
      k_3cERI  = 0
      k_3cwhat = 0
      CDFIT = .FALSE.
      if (iVcoul_opt.eq.1)CDFIT = .TRUE.
      XCFIT = .FALSE.
      if (iVxc_opt.eq.1)XCFIT = .TRUE.
c
      if (CDFIT)then
         if (oprint_time)
     &     call dft_tstamp(' Before 3c-2e initialize.')
c
         if(dft_mem3c(rtdb,natoms,npol,oprint_parm,oprint_3c2e,
     O                n3c_int,n3c_dbl,
     O                l_3ceri,k_3ceri, l_3cwhat,k_3cwhat)) then
            incore = .false.
            call dft_3cincor(n_batch, n3c_int, int_mb(k_3cwhat),
     &                       dbl_mb(k_3cERI), n3c_dbl)
            incore=.true.
         else
            if (me.eq.0 .and. oprint_3c2e)write(LuOut,3230)
            incore=.false.
         endif
      endif
 3230 format(/,10x,'Incore memory use for 3-center 2e- integrals is ',
     &     'turned off. ')
c     
c     Pre-compute mapping vectors
      if (.not.ma_push_get
     &     (mt_int,nshells_ao,'cntoce map',lcntoce,icntoce))
     &     call errquit('dft_scf_so:push_get failed', 13, MA_ERR)
      if (.not.ma_push_get
     &     (mt_int,nshells_ao*2,'cntoce map',lcntobfr,icntobfr))
     &     call errquit('dft_scf_so:push_get failed', 13, MA_ERR)
      if (.not.ma_push_get
     &     (mt_int,natoms*2,'cntoce map',lcetobfr,icetobfr))
     &     call errquit('dft_scf_so:push_get failed', 13, MA_ERR)
c     
      call build_maps(ao_bas_han, int_mb(icntoce), int_mb(icntobfr), 
     &     int_mb(icetobfr), natoms, nshells_ao)
c
c     Memory for reduced density matrix
      if (.not.ma_push_get(MT_Dbl,2*natoms*natoms,'rdens_atom',
     &     lrdens_atom,irdens_atom))
     &     call errquit('dft_scf_so: cannot allocate rdens_atom',
     & 0, MA_ERR)
c     
c     Memory for occupancy
      if (.not. ma_push_get(MT_Dbl,nbf_ao*2,'mo occ',l_occ,k_occ))
     &     call errquit('dft_scf_so: failed to alloc',999, MA_ERR)
c     
c     Allocate eigenvalue array
      leneval = 4*nbf_ao 
      if (.not.ma_push_get(MT_Dbl,leneval,'eval',l_eval,k_eval(1)))
     &     call errquit('dft_scf_so: cannot allocate eval',0, MA_ERR)
      k_eval(2) = k_eval(1) + nbf_mo
c
c     Allocate array for irreps
      if (.not.ma_push_get(mt_int,2*nbf_ao,'dft:irreps',l_ir,k_ir))
     &     call errquit('dft_scf_so: cannot allocate irreps',0, MA_ERR)

c     Memory for complex diagonalizer 
      if (.not.ma_push_get(MT_DCpl,nbf_mo*nbf_mo,'cpl a',la,ia))
     &     call errquit('dft_scf: cannot allocate cpl a',0, MA_ERR)
      if (.not.ma_push_get(MT_Dbl,nbf_mo,'cpl eval',lw,iw))
     &     call errquit('dft_scf: cannot allocate cpl eval',0, MA_ERR)
      llwork = max(1, 2*nbf_mo-1)
      if (.not.ma_push_get(MT_DCpl,llwork,'cpl work',lwork,iwork))
     &     call errquit('dft_scf: cannot allocate cpl work',0, MA_ERR)
      if (.not.ma_push_get(MT_Dbl,max(1,3*nbf_mo-2),'w.s',lrwork,
     &     irwork))
     &     call errquit('dft_scf: cannot allocate w.s',0, MA_ERR)
c
c     Get orbital overlap tolerance
      if (.not. rtdb_get(rtdb, 'dft:toll_s', MT_DBL, 1, toll_s))
     .     call errquit('dft_scf_so: lost toll_s ',0, RTDB_ERR)
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
      nopen = mult - 1
      nclosed = (nTotEl - nopen) / 2
c     
      if (.not. bas_name(ao_bas_han, basis_name, basis_trans))
     $     call errquit('dft_scf_so: bas_name?', 0, BASIS_ERR)
c     
c     get info for int2e_ and set sleazy tolerance
c     
      tol2e_sleazy = 1.d-3
      call scf_get_fock_param(rtdb, tol2e_sleazy)
c     
c     Force sleazy SCF into "direct" mode.
c     
      call fock_force_direct(rtdb)
cso
cso   allocate Fock matrix and movecs 
cso
c     real molecular orbital vectors
      if(.not.ga_create(mt_dbl,nbf_mo,nbf_mo,'Movecs Re',0,0, 
     &     g_moso(1)))     
     &     call errquit('dft_scf_so: error creating Movecs Re',0,
     &       GA_ERR)
c
c     imaginary molecular orbital vectors
      if(.not.ga_create(mt_dbl,nbf_mo, nbf_mo,'Movecs Im',0,0, 
     &     g_moso(2)))
     &     call errquit('dft_scf_so: error creating Movecs Im',0,
     &       GA_ERR)
c
      call ga_zero(g_moso(1))
      call ga_zero(g_moso(2))
c
      call ga_sync() 
c
c     real part of the fock matrix
      if(.not.ga_create(mt_dbl,nbf_mo,nbf_mo,'Fock Re',0,0, 
     &     g_fockso(1)))
     &     call errquit('dft_scf_so: error creating Fock Re',0, GA_ERR)
      call ga_zero(g_fockso(1))
c
c     imaginary part of the fock matrix
      if(.not.ga_create(mt_dbl,nbf_mo,nbf_mo,'Fock Im',0,0, 
     &     g_fockso(2)))
     &     call errquit('dft_scf_so: error creating Fock Im',0, GA_ERR)
      call ga_zero(g_fockso(2))
c
c     extra arrays
      if(.not.ga_create(mt_dbl,nbf_mo,nbf_mo,'old re',0,0, 
     &     g_old(1)))
     &     call errquit('dft_scf_so: error creating Old Re',0, GA_ERR)
      call ga_zero(g_old(1))
      if(.not.ga_create(mt_dbl,nbf_mo,nbf_mo,'old Im',0,0, 
     &     g_old(2)))
     &     call errquit('dft_scf_so: error creating Old Im',0, GA_ERR)
      call ga_zero(g_old(2))
c
c     real part of the spin-orbit density matrix
      if(.not.ga_create(mt_dbl,nbf_mo,nbf_mo,'DenMx Re',0,0, 
     &     g_densso(1)))
     &     call errquit('dft_scf_so: error creating DenMx Re',0, GA_ERR)
      call ga_zero(g_densso(1))
c
c     imaginary part of the spin-orbit density matrix
      if(.not.ga_create(mt_dbl,nbf_mo,nbf_mo,'DenMx Im',0,0, 
     &     g_densso(2)))
     &     call errquit('dft_scf_so: error creating DenMx Im',0, GA_ERR)
      call ga_zero(g_densso(2))
c
c     extra arrays
      if(.not.ga_create(mt_dbl,nbf_mo,nbf_mo,'Tmp ReIm',0,0, 
     &     g_tmp_ri))
     &     call errquit('dft_scf_so: error creating Tmp ReIm',0, GA_ERR)
      call ga_zero(g_tmp_ri)
c
c     spin-orbit matrices: 1->z, 2->y, 3->x
      if(.not.ga_create(mt_dbl,nbf_ao,nbf_ao,'so z',0,0, 
     &     g_so(1)))
     &     call errquit('dft_scf_so: error creating so z',0, GA_ERR)
      call ga_zero(g_so(1))
      if(.not.ga_create(mt_dbl,nbf_ao,nbf_ao,'so y',0,0, 
     &     g_so(2)))
     &     call errquit('dft_scf_so: error creating so y',0, GA_ERR)
      call ga_zero(g_so(2))
      if(.not.ga_create(mt_dbl,nbf_ao,nbf_ao,'so x',0,0, 
     &     g_so(3)))
     &     call errquit('dft_scf_so: error creating so x',0, GA_ERR)
      call ga_zero(g_so(3))
c
c     extra arrays
      if(.not.ga_create(mt_dbl, 2*nbf, 2*nbf,'old den', 0, 0, 
     &     g_damp_so(1)))
     &     call errquit('dft_scf_so: error creating damp ga', 0, GA_ERR)
      call ga_zero(g_damp_so(1))
      if(.not.ga_create(mt_dbl, 2*nbf, 2*nbf,'old den', 0, 0, 
     &     g_damp_so(2)))
     &     call errquit('dft_scf_so: error creating damp ga', 0, GA_ERR)
      call ga_zero(g_damp_so(2))
c
c     == get atomic zora contributions ==
      if (.not.rtdb_get(rtdb,'zora',mt_log,1,do_zora)) do_zora = .false.
c
c     == zora arrays ==
      if (do_zora) then
       if (me.eq.0) then
         call util_print_centered(LuOut,
     $        'Performing ZORA calculations', 23, .true.)
         write(LuOut,*)
       endif ! me
c
c      == get filename for zora data ==
       call util_file_name('zora_so',.false.,.false.,zorafilename)
c
c      == zora: scalar arrays ==
       if(.not.ga_duplicate(g_vxc(1),g_zora_sf(1),'g_zora_sf'))
     &    call errquit('dft_scf_so: error creating g_zora_sf',0, 
     &       GA_ERR)
       call ga_zero(g_zora_sf(1))
       if(.not.ga_duplicate(g_vxc(1),g_zora_sf(2),'g_zora_sf'))
     &    call errquit('dft_scf_so: error creating g_zora_sf',0, 
     &       GA_ERR)
       call ga_zero(g_zora_sf(2))
c
c      == zora: scalar energy scaling arrays ==
       if(.not.ga_duplicate(g_vxc(1),g_zora_scale_sf(1),
     G      'g_zora_scale_sf'))
     & call errquit('dft_scf_so: error creating g_zora_scale_sf',0, 
     &    GA_ERR)
       call ga_zero(g_zora_scale_sf(1))
       if(.not.ga_duplicate(g_vxc(1),g_zora_scale_sf(2),
     G      'g_zora_scale_sf'))
     & call errquit('dft_scf_so: error creating g_zora_scale_sf',0, 
     &    GA_ERR)
       call ga_zero(g_zora_scale_sf(2))
c
       if(.not.ga_duplicate(g_vxc(1),g_zora_so(1),'g_zora_so1'))
     & call errquit('dft_scf_so: error creating g_zora_so 1',0, GA_ERR)
       call ga_zero(g_zora_so(1))
       if(.not.ga_duplicate(g_vxc(1),g_zora_so(2),'g_zora_so2'))
     & call errquit('dft_scf_so: error creating g_zora_so 2',0, GA_ERR)
       call ga_zero(g_zora_so(2))
       if(.not.ga_duplicate(g_vxc(1),g_zora_so(3),'g_zora_so3'))
     & call errquit('dft_scf_so: error creating g_zora_so 3',0, GA_ERR)
       call ga_zero(g_zora_so(3))
c
c      == zora: spin-orbit energy scaling arrays ==
       if(.not.ga_duplicate(g_vxc(1),g_zora_scale_so(1),'g_zorasc_so1'))
     & call errquit('dft_scf_so: error creating g_zora_scale_so 1',0, 
     &    GA_ERR)
       call ga_zero(g_zora_scale_so(1))
       if(.not.ga_duplicate(g_vxc(1),g_zora_scale_so(2),'g_zorasc_so2'))
     & call errquit('dft_scf_so: error creating g_zora_scale_so 2',0, 
     &    GA_ERR)
       call ga_zero(g_zora_scale_so(2))
       if(.not.ga_duplicate(g_vxc(1),g_zora_scale_so(3),'g_zorasc_so3'))
     &   call errquit('dft_scf_so: error creating g_zora_scale_so',0, 
     &    GA_ERR)
       call ga_zero(g_zora_scale_so(3))
c
c      == generate an superposition of atomic densities ==
       call ga_zero(g_dens_at(1))
       if (ipol.gt.1) call ga_zero(g_dens_at(2))
       call guess_dens(rtdb, geom, ao_bas_han, g_dens_at)
       if (oskel) call ga_symmetrize(g_dens_at(1))
       if(ipol.gt.1) then
            call ga_copy(g_dens_at(1),g_dens_at(2))
            call ga_dscal(g_dens_at(1),dble(ntotel-nclosed)/(ntotel))
            call ga_dscal(g_dens_at(2),dble(nclosed)/(ntotel))
            if(oskel) call ga_symmetrize(g_dens_at(2))
       end if
c
c       == in case fon is used together with zora ==
c       == pstrace is queried in the grid code ==
c
       fon = .false.
       if (rtdb_get(rtdb,'dft:fon',mt_log,1,fon)) then 
         pstrace=ga_ddot(g_dens_at(1),g_s)
         pstrace=pstrace + ga_ddot(g_dens_at(2),g_s)
         if(ga_nodeid().eq.0) write (luout,'(5x,a,1x,e15.7)')
     &      'tr(P*S): ',pstrace
         if (.not. rtdb_put(rtdb, 'dft:pstrace', mt_dbl, 1, pstrace))
     &     call errquit('dft_scf: rtdb_put pstrace failed', 1, RTDB_ERR)
       end if ! fon check
c
c      == try reading the zora contributions from file ==
       icalczora = 0
       if (.not.dft_zora_read_so(zorafilename, nbf_ao, ipol, nmo, mult,
     &     g_zora_sf, g_zora_scale_sf, g_zora_so, g_zora_scale_so)) 
     &     icalczora = 1
c 
c       check if zora corrections have to recalculated
        do_zora_recalc = zora_recalc(rtdb)
        if (do_zora_recalc) then
           icalczora = 1
           if (me.eq.0) then
              write(luout,*) "Forcing recalculation of ZORA corrections"
           endif
        endif
c
       if (icalczora.eq.1) then
        if (me.eq.0) then
         call util_print_centered(LuOut,
     $        'Generating atomic ZORA corrections', 23, .true.)
         write(LuOut,*)
        end if ! me = 0
c
c       == read modelpotential flag ==
        use_modelpotential = .false.
        if (.not. rtdb_get(rtdb, 'modelpotential', mt_log, 1,
     &          use_modelpotential)) use_modelpotential = .false.
c
c       == calculate the zora atomic corrections ==
        Knucl=.false.
        ofinite=.false.
        status=rtdb_get(rtdb,'prop:ofinite',mt_log,1,ofinite) ! Check if gaussian nucl model requested
        if (.not.ma_alloc_get(mt_dbl,natoms,
     &                  'zetanuc',l_zetanuc,k_zetanuc))
     &    call errquit('dft_scf: ma failed',0,MA_ERR)
        call get_zetanuc_arr(geom,natoms,dbl_mb(k_zetanuc)) !  zetanuc_arr(i) i=1,natoms
        do iatom = 1,natoms ! == loop over the atoms ==
           dbl_mb(k_zetanuc+iatom-1)=dsqrt(dbl_mb(k_zetanuc+iatom-1)) !  Calc sqrt(zetanuc)
        enddo ! end-loop-iatom
        call zora_getv_so(rtdb, g_dens_at, g_zora_sf, g_zora_scale_sf,
     &                    g_zora_so, g_zora_scale_so, 
     &                    ofinite,dbl_mb(k_zetanuc),
     &                    Knucl,nexc)
         if (.not. MA_free_heap(l_zetanuc))  call errquit
     &           ('dft_scf_so:cannot popstack',111, MA_ERR)
c
c       == write out the atomic zora corrections to file ==
        if (.not.dft_zora_write_so(rtdb, ao_bas_han, zorafilename,
     &    nbf_ao, ipol, nmo, mult, g_zora_sf, g_zora_scale_sf, 
     &    g_zora_so, g_zora_scale_so))
     &    call errquit('dft_scf_so: dft_zora_write_so failed', 0, 
     &            DISK_ERR)
c
       end if ! icalczora
c
c       == for scalar calculations via the 2 component formalism
        if (do_purescalar) then 
            call ga_zero(g_zora_so(1))
            call ga_zero(g_zora_so(2))
            call ga_zero(g_zora_so(3))
            call ga_zero(g_zora_scale_so(1))
            call ga_zero(g_zora_scale_so(2))
            call ga_zero(g_zora_scale_so(3))
        end if
c
      end if  ! do_zora
c
c     check for fon input
c
      fon = .false.
      nmo_fon = 0
      nmo_fon2 = 0
      nel_fon = 0
      nel_fon2 = 0
      ncore_fon = 0
      if (rtdb_get(rtdb,'dft:fon',mt_log,1,fon)) then 
c
c       variable 'fon' should be true, otherwise there's 
c       something fishy going on:
        if (.not.fon) call errquit(
     &     'dft_scf_so: fon stored in RTDB but not .true.', 1,
     &       RTDB_ERR)

c       note: *_fon variables are read here not as arrays with
c       two elements, assuming that the input didn't specify
c       spects for alpha and beta separately

        if (.not.rtdb_get(rtdb,'dft:nmo_fon',mt_int,4,ntmp_fon)) then
          if (me.eq.0) then
            write(LuOut,*)"Error: fractional occupation 
     &         calculation specified without setting number of orbitals"
          end if
          call errquit('dft_scf_so: nmo_fon rtdb_get failed', 1,
     &       RTDB_ERR)
        else
          nmo_fon = ntmp_fon(1)
          nmo_fon2 = ntmp_fon(3)
        end if
        if (.not.rtdb_get(rtdb,'dft:nel_fon',mt_dbl,4,rtmp_fon)) then
          if (me.eq.0) then
            write(LuOut,*)"Error: fractional occupation 
     &         specified without setting number of electrons"  
          end if
          call errquit('dft_scf_so: nel_fon rtdb_get failed', 1,
     &       RTDB_ERR)
        else
          nel_fon = rtmp_fon(1)
          nel_fon2 = rtmp_fon(3)
        end if
        if (.not.rtdb_get(rtdb,'dft:ncore_fon',mt_int,2,ntmp_fon)) then
          if (me.eq.0) then
            write(LuOut,*)"Error: fractional occupation 
     &         calculation specified without setting filled levels"  
          end if
          call errquit('dft_scf_so:  nel_core rtdb_get failed', 1,
     &       RTDB_ERR)
        else
          ncore_fon = ntmp_fon(1)
        end if
        if (rtdb_get(rtdb, 'dft:debugfon', mt_log, 1,
     &     debug_fon)) continue

      else ! keyword not in RTDB
        fon = .false.
        debug_fon = .false.
      end if ! fon
c
      call ga_sync()
      call ga_zero(g_densso(1))
      call ga_zero(g_densso(2))
      spinor = .false. 
      spinor=spinor_guess(movecs_in)
c
      if(.not.spinor)then 
        
c       fractional occupations:
        if (fon) then
c         presumably good enough for initial guess. We won't tinker
c         with that. later, the fractional occupations
c         are caluclated explicitly
          anoca = noc(1)
          anocb = noc(2)
          noc(1)=(noc(1)+noc(2)-nel_fon-nel_fon2)/2
          noc(2)=(noc(1)+noc(2)-nel_fon-nel_fon2)/2
        endif                   ! fon
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
c        fon: undo temp setting of noc(:)
         if (fon)then 
           noc(1)=anoca
           noc(2)=anocb
         endif
c     
c     spinor occupancies
c     
         call dfill(nbf_mo, 0.0d0, dbl_mb(k_occ), 1)
         do i = 1, nTotEl
            dbl_mb(i-1+k_occ) = 1.0d0
         enddo
c     
c     map initial guess movecs from spin-free calculations g_moso(1) 
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
c
      endif  !if not spinor
c
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
            write(6,*)'dft_scf_so movecs output = ',movecs_in
            call errquit('dft_scf_so: could not read mo vectors',911,
     &       DISK_ERR)
         else 
            status = .true.
            status = status .and.
     $           movecs_read_so(movecs_in, dbl_mb(k_occ),
     $           dbl_mb(k_eval(1)), g_moso)
         endif
c     
         if (.not.status)then
            write(6,*)'dft_scf_so movecs output = ',movecs_in
            call errquit('dft_scf_so: could not read mo vectors',917,
     &       DISK_ERR)
         endif
c     
      endif  !spinor
c
c     Apply any swaps of orbitals
c
      call movecs_swap_so(rtdb,'dft',scftype,g_moso,
     &        dbl_mb(k_occ),dbl_mb(k_eval(1)))
c
c     Orthogonalize all input/guess vectors
c     Do not mix the occupied with the virtuals
c
      call ga_orthog_mos_so(ao_bas_han, g_moso)
c
c     Form Re and Im of density matrix
c     
c     applying fractional occupations

      nTotOcc = nTotEl  ! if there is no fractional occupation 
      if(fon)then 
        if (nmo_fon.lt.1) call errquit(
     &     'dft_scf_so:fon nmo_fon <1',
     &     1, INPUT_ERR)
        if (nel_fon.lt.0d0) call errquit(
     &     'dft_scf_so:fon nel_fon <0',
     &     1, INPUT_ERR)

        avg_fon = nel_fon/dble(nmo_fon)
        if (nmo_fon2.lt.1.or.nel_fon2.le.0d0) then
           avg_fon2 = 0.0d0
        else
           avg_fon2 = nel_fon2/dble(nmo_fon2)
        endif
        nTotOcc = ncore_fon + nmo_fon + nmo_fon2
c
c       initialization
        do i = 1,2*nbf_ao
          dbl_mb(i-1+k_occ) = 0.d0
        end do
c
        ncheck = 0d0
        nfilled = ncore_fon ! number of filled orbitals
c
c       partially fill molecular orbitals either from the core or valence
        if (do_core_fon) then       
c
c         partially filled core molecular orbitals
          do i = 1, nmo_fon
            if (i> 2*nbf_ao) call errquit(
     &       'dft_densm:fon focc index exceeds 2nbf error 1',
     &       i, INPUT_ERR)
            if (do_avg_fon) then
             dbl_mb(i-1+k_occ) = avg_fon
            else
             dbl_mb(i-1+k_occ) = nel_fon
            end if  ! do_avg_fon
            ncheck = ncheck + dbl_mb(i-1+k_occ)
          end do
c
c         fully filled outer orbitals
          do i = nmo_fon + 1, nmo_fon + nfilled
            if (i> 2*nbf_ao) call errquit(
     &       'dft_densm:fon focc index exceeds 2nbf error 2',
     &       i, INPUT_ERR)
            dbl_mb(i-1+k_occ) = 1d0
            ncheck = ncheck + dbl_mb(i-1+k_occ)
          end do
c
        else   ! don't do core_fon
c
c         fully filled core molecular orbitals
          do i = 1, nfilled
            if (i> 2*nbf_ao) call errquit(
     &       'dft_densm:fon focc index exceeds 2nbf error 1',
     &       i, INPUT_ERR)
            dbl_mb(i-1+k_occ) = 1d0
            ncheck = ncheck + dbl_mb(i-1+k_occ)
          end do
c
c         partially filled outer molecular orbitals
          do i = nfilled + 1, nfilled + nmo_fon
            if (i> 2*nbf_ao) call errquit(
     &       'dft_densm:fon focc index exceeds 2nbf error 2',
     &       i, INPUT_ERR)
c
            if (do_avg_fon) then
             dbl_mb(i-1+k_occ) = avg_fon
            else
             dbl_mb(i-1+k_occ) = nel_fon
            end if  ! do_avg_fon
            ncheck = ncheck + dbl_mb(i-1+k_occ)
          end do
c
c         second partially filled outer molecular orbitals
          do i = nfilled + nmo_fon + 1, nfilled + nmo_fon + nmo_fon2
            if (i> 2*nbf_ao) call errquit(
     &       'dft_densm:fon focc index exceeds 2nbf error 2',
     &       i, INPUT_ERR)
c
            if (do_avg_fon) then
             dbl_mb(i-1+k_occ) = avg_fon2
            else
             dbl_mb(i-1+k_occ) = nel_fon2
            end if  ! do_avg_fon
            ncheck = ncheck + dbl_mb(i-1+k_occ)
          end do
c
        endif ! do_core_fon
c
          if(abs(ncheck-dble(nTotEl)).gt.1d-3 .and. me.eq.0) then
            write(luout,*) ' frac. electrons ',ncheck,' vs ',nTotEl
          end if
      endif ! fon
c
c     the fractionally occupied mo's are scaled by the sqrt of the fon
c     
      if (fon) then 
         if (do_core_fon) then
          ilow = 1
          ihigh = nmo_fon
         else
          ilow = nfilled +1
          ihigh = nfilled + nmo_fon + nmo_fon2
         end if
c
         do i = ilow,ihigh
          if (i> nbf_mo) call errquit(
     &        'dft_densm:fon g_moso index exceeds nbf_mo',
     &        i, INPUT_ERR)
           scale = sqrt(dbl_mb(i-1+k_occ)) ! sqrt(occupation number)
           call ga_scale_patch(g_moso(1),1,nbf_mo,i,i,scale)
           call ga_scale_patch(g_moso(2),1,nbf_mo,i,i,scale)
         end do ! i loop
      end if ! fon
c
c     fractional occupations by orbital
      switch_focc = 0
      status=rtdb_get(rtdb,'focc:occup_switch', mt_int,1,switch_focc)
      if (switch_focc.eq.1 .and. status) then ! using specified occupations 
        call dft_frac_mo_so(rtdb,g_moso,dbl_mb(k_occ),g_densso,
     &                     nbf_mo,nTotOcc)
      else
        call dft_densm_so(g_densso, g_moso, nbf_ao, nTotOcc) ! standard way
      endif
c
c     restore the scaled mo's
c
      if (fon) then 
         if (do_core_fon) then
          ilow = 1
          ihigh = nmo_fon
         else
          ilow = nfilled + 1
          ihigh = nfilled + nmo_fon + nmo_fon2
         end if
         do i = ilow,ihigh
           if (i> nbf_mo) call errquit(
     &        'dft_densm:fon g_moso index exceeds nbf_mo',
     &        i, INPUT_ERR)
           if (dbl_mb(i-1+k_occ) < 1d-4) call errquit(
     &        'dft_densm:fon frac occup < 1E-4. Aborting',
     &        i, INPUT_ERR)           
           scale = 1d0/sqrt(dbl_mb(i-1+k_occ)) ! 1/sqrt(occupation number)
           call ga_scale_patch(g_moso(1),1,nbf_mo, i, i, scale)
           call ga_scale_patch(g_moso(2),1,nbf_mo, i, i, scale)
         end do
       end if ! fon
c      
c     calculate the spin-free density matrix from the spin-orbit density matrix
c
      call ga_zero(g_dens(1))
      call ga_zero(g_dens(2))
      call ga_dens_sf(g_dens, g_densso(1), nbf_ao)

c     check fon, calculate tr[P S] and print
c      if (debug_fon) call dft_pstrace(g_dens,ao_bas_han,nbf_ao,oskel)
      if (fon) then
        pstrace=ga_ddot(g_dens(1),g_s)
        pstrace=pstrace + ga_ddot(g_dens(2),g_s)
        if(ga_nodeid().eq.0) write (luout,'(5x,a,1x,e15.7)')
     &     'tr(P*S): ',pstrace
        if (.not. rtdb_put(rtdb, 'dft:pstrace', mt_dbl, 1, pstrace))
     &     call errquit('dft_scf: rtdb_put pstrace failed', 1, RTDB_ERR)
      end if
c
c     Tidy up SCF
c     
      call fock_2e_tidy(rtdb)
c     
c     set initial coulomb acc
c     
c     write(6,*)' movecs_guess = ',movecs_guess
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

c     rlshift = 2.0
      endif
      iswitc = 0
#ifdef LOOSE_START
      if (ltight)then
#endif
         itol2e = itol_max
         iAOacc = iAOacc_max
         tol_rho = tol_rho_max
         iswitc = 1
#ifdef LOOSE_START
      else
         itol2e = min(itol_min,itol_max)
         iAOacc = min(iAOacc_min,iAOacc_max)
         tol_rho = max(tol_rho_min,tol_rho_max)
      endif
#endif
c     
      tol2e = 10.d0**(-itol_max)
c     
c     Restore SCF parameters
c     
      call scf_get_fock_param(rtdb, tol2e)
c     
c     Put the total density matrix in g_dens(1)
c     
      call ga_dadd(one,g_dens(1),one,g_dens(2),g_dens(1))
c     
c     Call to Mulliken Pop Analysis for initial density
c     
      if (mulliken)then
         if (me.eq.0)call dft_header
     &        (' Total Density - Mulliken Population Analysis')
         call mull_pop(geom,ao_bas_han,g_dens(1),g_s,'total')
c     
c        analysis of spin density
c     
         if (me.eq.0) call dft_header
     &        (' Spin Density - Mulliken Population Analysis')
         call ga_dadd(one,g_dens(1),-2.d0,g_dens(2),g_dens(2))
         call mull_pop(geom,ao_bas_han,g_dens(2),g_s,'spin') 
c     
c        restore beta density in g_dens(2)
c     
         call ga_dadd(one,g_dens(1),-1.d0,g_dens(2),g_dens(2))
         call ga_dscal(g_dens(2),0.5d0)
      endif  ! mulliken
 
      iter = 1   ! iteration counter
c
c     dispersion bit
      Edisp = 0.d0
      if (.not.rtdb_get(rtdb, 'dft:disp', mt_log, 1, disp))
     &   disp=.false.
      if(disp.or.xc_chkdispauto())
     &     call xc_vdw(rtdb,geom,Edisp,dum,'energy')
c     
c     Top of infinite SCF iteration loop
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
c     start DFT_SCF timer
c     
      start_wall = util_wallsec()
      start_cpu = util_cpusec()
      dft_time = -start_cpu
c
      if (oprint_time)
     &     call dft_tstamp('   Before SCF iter loop. ')
c 
 3000 continue
c
      last_time_energy = .false.
c
      det_eng = .false. 
      idet = 0 
c
c     SCF loop
c
 1000 continue   ! top of the scf loop

      if (me.eq.0 .and. oprint_conv_details)
     &   write(LuOut,124)damping, levelshifting, diising
 124  format(10x,' DAMPING=',l1,' LEVELSHIFTING=',l1,
     &           ' DIISING=',l1)
c
      if (me.eq.0.and.oprint_tol)write(LuOut,3234)itol2e,iAOacc,iXCacc
 3234 format(10x,'itol2e=',i2,' iAOacc=',i2,' iXCacc=',i2)

      Ecoul  = ZERO
      Exc(1) = ZERO
      Exc(2) = ZERO
c
      call ga_zero(g_fockso(1))
      call ga_zero(g_fockso(2))
c     
c     Accumulate core hamiltonian into Fock matrix; 
c     compute core energy
c     
      call ga_zero(g_fock)
      call int_1e_ga(ao_bas_han, ao_bas_han, g_fock, 'kinetic', oskel)
      call int_1e_ga(ao_bas_han, ao_bas_han, g_fock, 'potential', oskel)
c
c     Re(Dsf)=Re(Daa)+Re(Dbb)=g_dens(1), Re(Dbb) = g_dens(2)
c     <Hsf> = Re(Dsf) dot T+Vsf 
      Ecore = ga_ddot(g_dens(1), g_fock)
      if (me.eq.0.and.oprint_info)
     &     write(luout,'(1x,a,e15.7)') "Ecore =", ecore
c
c     Restore alpha and beta densities 
c     g_dens(1) = Re(Daa), g_dens(2) = Re(Dbb)
      call ga_dadd(one, g_dens(1), onem, g_dens(2), g_dens(1))
c
c     Calculate cosmo charges
      if (cosmo_on.and.cosmo_phase.eq.2) then
       cosmo_file = "cosmo.xyz"
       call cosmo_charges_from_dmat(rtdb, ao_bas_han, geom,
     &          ecosmo, odbug, ipol, g_dens, cosmo_file)
       call cosmo_fock(ao_bas_han, g_fock)
       call ga_dadd(one, g_dens(1), one, g_dens(2), g_dens(1))
       Ecore = ga_ddot(g_dens(1), g_fock)
       call ga_dadd(one, g_dens(1), onem, g_dens(2), g_dens(1))
      end if ! cosmo check
c
c     Incorporate g_fock into supermatrix g_fockso
      call ga_fock_sf(g_fock, g_fockso(1), nbf_ao) 
c
c     Accumulate s.o. contribution to fock matrix 
c   
      call ga_zero(g_so(1))
      call ga_zero(g_so(2))
      call ga_zero(g_so(3))
c
c     Calculate the spin-orbit contributions from the ecp
      if( .not. do_zora .and. .not. do_purescalar) then
         call int_1e_ga(ao_bas_han, ao_bas_han, g_so, 'so', oskel)
         call ga_scale(g_so(1),dble(0.5d0))  ! z
         call ga_scale(g_so(2),dble(0.5d0))  ! y
         call ga_scale(g_so(3),dble(0.5d0))  ! x
c
c        Add in the s.o. contribution to the fock matrix
         call ga_fock_so(g_so, g_fockso, nbf_ao)
      end if
c
c     Calculate one-electron energy terms
c  
      call ga_zero(g_tmp) 
      call ga_dens_so(g_tmp, g_densso, nbf_ao, 'z') 
      Ecore = Ecore + ga_ddot(g_tmp, g_so(1)) 
c
      call ga_zero(g_tmp) 
      call ga_dens_so(g_tmp, g_densso, nbf_ao, 'y') 
      Ecore = Ecore + ga_ddot(g_tmp, g_so(2)) 
c
      call ga_zero(g_tmp) 
      call ga_dens_so(g_tmp, g_densso, nbf_ao, 'x') 
      Ecore = Ecore + ga_ddot(g_tmp, g_so(3))
c
      if (do_zora) then
       call ga_zero(g_tmp) 
       call ga_dens_so(g_tmp, g_densso, nbf_ao, 'z') 
       Ecore = Ecore + ga_ddot(g_tmp, g_zora_so(1))
c
       call ga_zero(g_tmp) 
       call ga_dens_so(g_tmp, g_densso, nbf_ao, 'y') 
       Ecore = Ecore + ga_ddot(g_tmp, g_zora_so(2))
c
       call ga_zero(g_tmp) 
       call ga_dens_so(g_tmp, g_densso, nbf_ao, 'x') 
       Ecore = Ecore + ga_ddot(g_tmp, g_zora_so(3))
c
       Ezora_sf = ga_ddot(g_dens(1),g_zora_sf(1))
     &     + ga_ddot(g_dens(2),g_zora_sf(2))
       Ecore = Ecore + Ezora_sf
       if (me.eq.0.and.oprint_info)
     &     write(luout,'(1x,a,e15.7)') "Ecore+SC+SO(z+y+x) =", ecore
c
       call ga_copy(g_zora_sf(1),g_vxc(1))
       call ga_copy(g_zora_sf(2),g_vxc(2))
c
c      == patch scalar zora parts into Fock matrix ==
       call ga_sync()
       call ga_dadd_patch( 1.d0, g_fockso(1), 1, nbf_ao,
     &     1, nbf_ao,
     &     1.0d0, g_vxc(1),  1, nbf_ao,
     &     1, nbf_ao,
     &     g_fockso(1), 1, nbf_ao,
     &     1, nbf_ao)
       call ga_dadd_patch( 1.d0, g_fockso(1), 1+nbf_ao, nbf_mo,
     &     1+nbf_ao, nbf_mo,
     &     1.0d0, g_vxc(2),  1, nbf_ao,
     &     1, nbf_ao,
     &     g_fockso(1), 1+nbf_ao, nbf_mo,
     &     1+nbf_ao, nbf_mo)
c
c      == patch in spin-orbit zora part ==
       if (.not.do_purescalar)
     &     call ga_fock_so(g_zora_so, g_fockso, nbf_ao)
c
      end if ! do_zora
c
c     Pre-compute reduced total density matrices over atoms
c 
      call dfill(ipol*natoms*natoms, 0.0d0, dbl_mb(irdens_atom), 1)
      nscr = nbf_ao_mxnbf_ce*nbf_ao_mxnbf_ce
      if (.not.ma_push_get(MT_Dbl,nscr,'scr',lscr,iscr)) ! allocate lscr
     &   call errquit('dft_scf: cannot allocate scr',0, MA_ERR) 
      call util_ga_mat_reduce(nbf_ao, natoms, int_mb(icetobfr), g_dens, 
     &                        ipol, dbl_mb(irdens_atom), 'absmax', 
     &                        dbl_mb(iscr), nbf_ao_mxnbf_ce,.true.)
c      write(*,'("irdens",5f10.7)')
c     &     (dbl_mb(irdens_atom+i),i=0,ipol*natoms*natoms-1)
      if (.not.ma_pop_stack(lscr))  ! free lscr
     &   call errquit('dft_scf: cannot pop stack:lscr',0, MA_ERR)
c
      g_vxc(1) = g_tmp
      call ga_zero(g_vxc(1))
      rho_n = 0.0d0
      call ga_zero(g_vxc(2))
c
      ecoul = 0.d0
      if (cam_exch) call case_setflags(.false.) ! unset attenuation for coulomb
      if (CDFIT)then
        call getv_coul_cd(rtdb, ecoul, g_densso, g_fockso, g_vxc,
     &     dbl_mb(k_3cERI), nTotEl, n_batch, n3c_int,
     &     int_mb(k_3cwhat), n3c_dbl, IOLGC, natoms, iVcoul_opt)
      endif  ! CDFIT
c
c     Note that g_dens(1) now contains the alpha density
c     matrix and g_dens(2) contains the beta
c     
c     Pre-compute reduced alpha and beta density matrices over atoms
c     
      call dfill(ipol*natoms*natoms, 0.0d0, dbl_mb(irdens_atom), 1)
      nscr = nbf_ao_mxnbf_ce*nbf_ao_mxnbf_ce
      if (.not.ma_push_get(MT_Dbl,nscr,'scr',lscr,iscr)) ! allocate lscr
     &     call errquit('dft_scf: cannot allocate scr',0, MA_ERR)
      call util_ga_mat_reduce(nbf_ao, natoms, int_mb(icetobfr), 
     &     g_dens, ipol, dbl_mb(irdens_atom), 
     &     'absmax', dbl_mb(iscr), nbf_ao_mxnbf_ce,.true.)
      if (.not.ma_pop_stack(lscr))  ! free lscr
     &  call errquit('dft_scf: cannot pop stacklscr:',0, MA_ERR)
c     
c     == Coulomb part via traditional 4c approach ==
c
      if (.not.CDFIT) 
     &    call getv_coul(rtdb, ecoul, g_densso, g_fockso, g_vxc)
      if (me.eq.0.and.oprint_info) write(luout,'(1x,a,e15.7)')
     &   "ECoul =", ecoul
c
c     Compute the XC potential and energy.
c
c     == HF Exchange part ==

c     == attenuation ==
      if (cam_exch) call case_setflags(.true.) ! set attenuation, if needed
c
c     == D0, V0(Daa + Dbb) ==
      exch0r = 0.d0
      exch0i = 0.d0
      exch0 = 0.d0
      if (me.eq.0.and.oprint_info) write(luout,*) "calling getv_exch0"
      call getv_exch0(rtdb, exch0r, g_densso, g_fockso, g_vxc, 'real')
      call getv_exch0(rtdb, exch0i, g_densso, g_fockso, g_vxc, 'imag')
      exch0 = exch0r + exch0i
      if (me.eq.0.and.oprint_info) write(luout,'(1x,a,3e15.7)')
     &     "Exch0Re, Exch0I, sum =", exch0r, exch0i, exch0
c
c     == DZ, i/2 Vz*(Daa - Dbb) ==
      exchZ = 0.d0
      exchZr = 0.d0
      exchZi = 0.d0
      if (me.eq.0.and.oprint_info) write(luout,*) "calling getv_exchZ"
      call getv_exchZ(rtdb, exchZr, g_densso, g_fockso, g_vxc, 'real')
      call getv_exchZ(rtdb, exchZi, g_densso, g_fockso, g_vxc, 'imag')
      exchZ = exchZr + exchZi
      if (me.eq.0.and.oprint_info) write(luout,'(1x,a,3e15.7)')
     &     "ExchZRe, ExchZI, sum =", exchZr, exchZi, exchZ
c
c     == DX, i/2 Vx*(Dab + Dba) ==
      exchX = 0.d0
      exchXr = 0.d0
      exchXi = 0.d0
      if (.not.do_purescalar) then 
        if (me.eq.0.and.oprint_info) write(luout,*) "calling getv_exchX"
        call getv_exchX(rtdb, exchXr, g_densso, g_fockso, g_vxc, 'real')
        call getv_exchX(rtdb, exchXi, g_densso, g_fockso, g_vxc, 'imag')
        exchX = exchXr + exchXi
        if (me.eq.0.and.oprint_info) write(luout,'(1x,a,3e15.7)')
     &     "ExchXRe, ExchXI, sum =", exchXr, exchXi, exchX
      end if ! .not.do_purescalar
c
c     == DY, 1/2 Vy*(Dab - Dba)
      exchY = 0.d0
      exchYr = 0.d0
      exchYi = 0.d0
      if (.not.do_purescalar) then 
        if (me.eq.0.and.oprint_info) write(luout,*) "calling getv_exchY"
        call getv_exchY(rtdb, exchYr, g_densso, g_fockso, g_vxc, 'real')
        call getv_exchY(rtdb, exchYi, g_densso, g_fockso, g_vxc, 'imag')
        exchY = exchYr + exchYi
        if (me.eq.0.and.oprint_info) write(luout,'(1x,a,3e15.7)')
     &     "ExchYRe, ExchYI, sum =", exchYr, exchYi, exchY
      end if ! .not.do_purescalar
c
c     == attenuation ==
      if (cam_exch) call case_setflags(.false.) ! turn off attenuation
c
c     == tally up exchange energy ==
      exchtot = exch0 + exchZ + exchX + exchY ! HF exchange
c
c     == DFT Exchange-Correlation ==
c
      edftxc(1) = 0.d0
      edftxc(2) = 0.d0
      if (me.eq.0.and.oprint_info)  write(luout,*) "calling getv_dftxc"
      call getv_dftxc(rtdb, nexc, edftxc, rho_n, g_densso,
     &         g_fockso, g_vxc, 'none')
c
      if (nexc.eq.1) then       ! don't decompose if nexc = 1
        Exc(1) = edftxc(1) + edftxc(2) ! dft exchange + correlation
        Exc(2) = zero
      else
        Exc(1) = edftxc(1)      ! exchange
        Exc(2) = edftxc(2)      ! correlation
      end if
c
      if (me.eq.0.and.oprint_info) write(luout,'(1x,a,3e15.7)')
     &   'Etot, Ecoul, Exch =',
     &   Ecore + Ecoul + exchtot + Exc(1) + Exc(2) + Edisp,
     &   Ecoul,
     &   exchtot  + Exc(1) + Exc(2)
c
      if (oprint_time)
     &     call dft_tstamp(' End of parallel region. ')
c     
c     Calculate the total electronic energy.
c     
      if (nexc.eq.1)then
         Etnew = Ecore + Ecoul + exchtot + Exc(1) + Edisp
         if(det_eng)goto 2001
      else
         Etnew = Ecore + Ecoul + exchtot + Exc(1) + Exc(2) + Edisp
         if(det_eng)goto 2001
      endif
c
c     Add in cosmo contributions
c
      if (cosmo_on.and.cosmo_phase.eq.2) then
        if (do_cosmo_smd) then
             etnew = etnew + ecosmo + gcds
        else
             etnew = etnew + ecosmo
        end if
      end if
c
c     == is this a oneshot "noscf" type calculation, no diagonalization ==
      noscf = .false.
      if (.not.rtdb_get(rtdb, 'dft:noscf', mt_log, 1, noscf)) 
     &   noscf=.false.
      if (noscf) then
        last_time_energy = .true.
        converged = .true.
        levelshifting = .false. ! turn off levelshifting for noscf
      end if
c
      if (last_time_energy)then
c     
c        If open shell put the total density matrix back in  g_dens(1) and quit.
c     
         call ga_dadd(one, g_dens(1), one, g_dens(2), g_dens(1))
         goto 2000
      endif
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
      call ga_symmetrize(g_fock)
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
      g_scr = ga_create_atom_blocked(geom, AO_bas_han, 'ga scr')
c     
c     Put s-1/2 in g_scr.
c     
      iw = 2
      call diis_bld12_so(toll_s, svals, g_svecs, g_scr, 
     &     g_tmp, nbf_ao, iw)
c     
c     map s-1/2 to the nbf_mo by nbf_mo g_scr2 
c    
      if(.not.ga_create(mt_dbl, 2*nbf, 2*nbf,'scr2', 0, 0, g_scr2))
     &     call errquit('dft_scf_so: error creating scr2',0, GA_ERR)
      call ga_zero(g_scr2)
      call ga_fock_sf(g_scr, g_scr2, nbf_ao)
c     
c     Transform Fock matrix.
c     
      call ga_zero(g_tmp_ri)   
      call ga_dgemm('T', 'N', nbf_mo, nbf_mo, nbf_mo, one, 
     &                 g_scr2, g_fockso(1), zero, g_tmp_ri)
      call ga_dgemm('N', 'N', nbf_mo, nbf_mo, nbf_mo, one, 
     &                 g_tmp_ri, g_scr2, zero, g_fockso(1))

      call ga_zero(g_tmp_ri)   
      call ga_dgemm('T', 'N', nbf_mo, nbf_mo, nbf_mo, one, 
     &                 g_scr2, g_fockso(2), zero, g_tmp_ri)
      call ga_dgemm('N', 'N', nbf_mo, nbf_mo, nbf_mo, one, 
     &                 g_tmp_ri, g_scr2, zero, g_fockso(2)) 
c     
c     Level shifting is implemented here (similarity 
c     transformation before standard eigensolver).  Note,
c     levelshifting is appropriate once a transformation
c     is available which makes the resulting Fock matrix 
c     diagonally dominant, e.g., in an approximate MO basis.  
c     Also note, there are many matrix multiplies with S^+-1/2 
c     which are redundant if one is sure that the former basis
c     is orthonormal.
c     
c     levelshifting = .false. 
      if (levelshifting)then
c     
c     save the old vectors 
c     
         call ga_copy(g_moso(1), g_old(1))
         call ga_copy(g_moso(2), g_old(2))
c     
c        Build a matrix which is diagonal in the "MO" rep,
c        back-transform, and shift the current Fock matrix
c     
c        Use S^+1/2 * old movecs (as a transform).
c     
         iw = 3
         call diis_bld12_so(toll_s, svals, g_svecs, g_scr, 
     &                   g_tmp, nbf_ao, iw)
         call ga_zero(g_scr2)
         call ga_fock_sf(g_scr, g_scr2, nbf_ao)
         call ga_zero(g_tmp_ri)
         call ga_dgemm('N', 'N', nbf_mo, nbf_mo, nbf_mo, one, 
     &                 g_scr2, g_moso(1), zero, g_tmp_ri)
         call ga_copy(g_tmp_ri,  g_moso(1)) 

         call ga_zero(g_tmp_ri) 
         call ga_dgemm('N', 'N', nbf_mo, nbf_mo, nbf_mo, one, 
     &                 g_scr2, g_moso(2), zero, g_tmp_ri)
         call ga_copy(g_tmp_ri,  g_moso(2)) 
c     
c        Build diagonal matrix.
c     
         call ga_zero(g_tmp_ri)
         do j = nTotOcc+1+me, nbf_mo, nproc
            call ga_put(g_tmp_ri, j, j, j, j, rlshift, 1)
         enddo
c     
c        Transform this into "AO" basis and add to current 
c        Fock matrix
c     
         call ga_dgemm('N', 'N', nbf_mo, nbf_mo, nbf_mo, one, 
     &                 g_moso(1), g_tmp_ri, zero, g_scr2)
         call ga_dgemm('N', 'T', nbf_mo, nbf_mo, nbf_mo, one, 
     &                 g_scr2, g_moso(1), one, g_fockso(1))
         call ga_dgemm('N', 'N', nbf_mo, nbf_mo, nbf_mo, one, 
     &                 g_moso(2), g_tmp_ri, zero, g_scr2)
         call ga_dgemm('N', 'T', nbf_mo, nbf_mo, nbf_mo, one, 
     &                 g_scr2, g_moso(2), one, g_fockso(1))

         call ga_dgemm('N', 'N', nbf_mo, nbf_mo, nbf_mo, one, 
     &                 g_moso(1), g_tmp_ri, zero, g_scr2)
         call ga_dgemm('N', 'T', nbf_mo, nbf_mo, nbf_mo, mone, 
     &                 g_scr2, g_moso(2), one, g_fockso(2))
         call ga_dgemm('N', 'N', nbf_mo, nbf_mo, nbf_mo, one, 
     &                 g_moso(2), g_tmp_ri, zero, g_scr2)
         call ga_dgemm('N', 'T', nbf_mo, nbf_mo, nbf_mo, one, 
     &                 g_scr2, g_moso(1), one, g_fockso(2))
      else
        rlshift = 0.0
      endif
c     
c     Solve for the eigenvalues and eigenvectors of the Hamiltonian.
c     
      call ga_symmetrize(g_fock)
      if (oprint_intermediate_fock)then     
      endif
cso#if defined(PARALLEL_DIAG)
cso      call ga_diag_std(g_fock, g_tmp, Dbl_MB(k_eval(ispin)))
cso#else 
cso      call ga_diag_std_seq(g_fock, g_tmp, Dbl_MB(k_eval(ispin)))
cso#endif
cso      call ga_diag_compl(g_fockso(1), g_fockso(2), g_moso(1), 
cso     &                   g_moso(2),  Dbl_MB(k_eval(1))) 
c      write(*,*)"compare"
      call ga_zero(g_moso(1))
      call ga_zero(g_moso(2))
#ifdef SCALAPACK
      nmoout=0
      call dft_scaleig_so(g_fockso,g_moso,dbl_mb(k_eval(1)),
     A      nmoout)
#else
      call dft_eig_so(nbf_mo, g_fockso, g_moso, dbl_mb(k_eval(1)))
#endif
c      write(*,'("before transform")')
c      write(*,*)(Dbl_mb(k_eval(1)+i),i=0,nbf_mo-1)
c      call ga_print(g_moso(1))
c      call ga_print(g_moso(2))
c     
c     Check HOMO/LUMO gap.
c     
      homo = Dbl_MB(k_eval(1)+nTotEl-1)
      lumo = Dbl_MB(k_eval(1)+nTotEl)
c     
c     If levelshifting then tidy up.
c  
      if (levelshifting)then
c     
c        Put S^-1/2 back in g_scr2 (use g_fock as temp scr).
c     
         iw = 2
         call diis_bld12_so(toll_s, svals, g_svecs, g_scr, 
     &                   g_fock, nbf_ao, iw)
         call ga_zero(g_scr2)
         call ga_fock_sf(g_scr, g_scr2, nbf_ao)
      endif
c     
c     Back-transform eigenvectors with S^-1/2.
c     
      call ga_dgemm('N', 'N', nbf_mo, nbf_mo, nbf_mo, one, 
     &     g_scr2, g_moso(1), zero, g_fockso(1))
      call ga_dgemm('N', 'N', nbf_mo, nbf_mo, nbf_mo, one, 
     &     g_scr2, g_moso(2), zero, g_fockso(2))
      if (.not. ga_destroy(g_scr)) call errquit
     &     ('dft_scf: could not destroy g_scr', 0, GA_ERR)
      if (.not. ga_destroy(g_scr2)) call errquit
     &     ('dft_scf_so: could not destroy g_scr2', 0, GA_ERR)
c     
c     Keep orbital ordering according to principle
c     of maximum overlap with previous iteration.
c
      if (lmaxov)
     .     call dft_mxovl(ao_bas_han, dbl_mb(k_eval(ispin)),
     &                  g_tmp, g_movecs(ispin), g_s,g_fock,
     .        noc,oprint_conv_details,homo,lumo)

      call ga_copy(g_fockso(1), g_moso(1))
      call ga_copy(g_fockso(2), g_moso(2))

c        
c     determine homo-lumo gap 
c
      homo_lumo_gap = min(homo_lumo_gap, (lumo-homo-rlshift))
      if (me.eq.0 .and. oprint_conv_details)
     &   write(LuOut,4224)homo,lumo,rlshift, homo_lumo_gap
 4224 format(10x,' HOMO = ',f6.2,' LUMO = ',f6.2,
     &           ' RLSHIFT = ',f6.2,' HL_GAP = ',f6.2)
c     
      call ga_sync
c     
c     Save previous density for convergence check.
c     
      call ga_copy(g_dens(1), g_movecs(1))
      call ga_copy(g_dens(2), g_movecs(2))
c
c     symmetry adapt vectors?
c
      if (oadapt)then
         call scf_sym_adapt_so(ao_bas_han, g_moso,
     &        oprint_syma, 2*nbf_ao, name,
     &        .true., 
     &        int_mb(k_ir))
      endif      
c
c     save the old density matrix for damping 
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
c     
c     the fractionally occupied mo's are scaled by the sqrt of the fon
c     
      if(fon)then 
         if (do_core_fon) then
          ilow = 1
          ihigh = nmo_fon
         else
          ilow = nfilled + 1
          ihigh = nfilled + nmo_fon + nmo_fon2
         end if
         do i = ilow, ihigh
           if (i> nbf_mo) call errquit(
     &        'dft_densm:fon g_moso index exceeds nbf_mo',
     &        i, INPUT_ERR)
           scale = sqrt(dbl_mb(i-1+k_occ)) ! sqrt(occupation number)
           call ga_scale_patch(g_moso(1),1,nbf_mo,i,i,scale)
           call ga_scale_patch(g_moso(2),1,nbf_mo,i,i,scale)
         end do
      end if ! fon
c
c     fractional occupations by orbital
      switch_focc = 0
      status=rtdb_get(rtdb,'focc:occup_switch', mt_int,1,switch_focc)
      if (switch_focc.eq.1 .and. status) then ! using specified occupations
        call dft_frac_mo_so(rtdb,g_moso,dbl_mb(k_occ),g_densso,
     &                     nbf_mo,nTotOcc)
      else
        call dft_densm_so(g_densso, g_moso, nbf_ao, nTotOcc) ! standard way
      endif
c     
c     restore the scaled mo's
c     
      if (fon) then 
         if (do_core_fon) then
          ilow = 1
          ihigh = nmo_fon
         else
          ilow = nfilled + 1
          ihigh = nfilled + nmo_fon + nmo_fon2
         end if
         do i = ilow,ihigh
           if (i> nbf_mo) call errquit(
     &        'dft_densm:fon g_moso index exceeds nbf_mo',
     &        i, INPUT_ERR)
           if (dbl_mb(i-1+k_occ) < 1d-4) call errquit(
     &        'dft_densm:fon frac occup < 1E-4. Aborting',
     &        i, INPUT_ERR)           
           scale = 1d0/sqrt(dbl_mb(i-1+k_occ)) ! 1/sqrt(occupation number)
           call ga_scale_patch(g_moso(1),1,nbf_mo,i,i,scale)
           call ga_scale_patch(g_moso(2),1,nbf_mo,i,i,scale)
         end do
       end if ! fon
cso
cso   g_dens(1)=Re(Daa)+Re(Dbb) and g_dens(2)=Re(Dbb)
cso   For coulomb and xc potentials only the alpha, Re(Daa), and  
cso   beta, Re(Dbb) densities are needed 
cso
      call ga_zero(g_dens(1))
      call ga_zero(g_dens(2))
      call ga_dens_sf(g_dens, g_densso(1), nbf_ao)

c     check g_dens, calculate tr[P S] and print
c     if (debug_fon) call dft_pstrace(g_dens,ao_bas_han,nbf_ao,oskel)
      if (fon) then
        pstrace=ga_ddot(g_dens(1),g_s)
        pstrace=pstrace + ga_ddot(g_dens(2),g_s)
        if(ga_nodeid().eq.0) write (luout,'(5x,a,1x,e15.7)')
     &       'tr(P*S): ',pstrace 
        if (.not. rtdb_put(rtdb, 'dft:pstrace', mt_dbl, 1, pstrace))
     &     call errquit('dft_scf: rtdb_put pstrace failed', 1, RTDB_ERR)
      end if
c
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
c
c     ----- cosmo message -----
c
         if(cosmo_on.and.cosmo_phase.eq.1) then
            write(LuOut,909)
         elseif(cosmo_on) then
            write(LuOut,910)
         endif
 909     format(2x,'   COSMO gas phase')
 910     format(2x,'   COSMO solvation phase')
c
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
c     ecce ouput
c
      call ecce_print1 ('iteration counter', mt_int, iter, 1)
      call ecce_print1 ('iterative total energy difference', 
     &                  mt_dbl, -etold+etnew, 1)
      call ecce_print1 ('iterative total density difference', 
     &                  mt_dbl, sqrt(rms(1)), 1)
c
      call ga_sync
c     
c     save eigenvectors to movecs file
c     
      if (.not.movecs_write_so
     $     (rtdb, ao_bas_han, movecs_out, 'sodft', title,
     &     nbf_mo, dbl_mb(k_occ), dbl_mb(k_eval(1)), g_moso))
     &     call errquit('dft_scf_so: movec_write failed', 0, DISK_ERR)
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
         call dft_mxspin_ovlp(nbf_ao,nmo,ao_bas_han,g_movecs(1), 
     &       g_movecs(2), g_tmp)
      endif
c     
c     computation of <S2> for open shell
c     
      if ((ipol.gt.1).and.(oprint_interm_S2)) then
         
         call dft_s2_value(geom, AO_bas_han, .false., noc(1), noc(2),
     &        nbf_ao, g_dens(1), g_dens(2))
      endif
c     
c     Form the total density matrix.
c     
      call ga_dadd(one, g_dens(1), one, g_dens(2), g_dens(1))
      call ga_sync
c     
c     Check for SCF convergence.
c     
      call ga_sync
c
      call dft_scfcvg(rms, derr, Etold, Etnew, 
     &     e_conv, d_conv, g_conv, ipol, 
     &     iter, iterations, idone, rtdb,
     &     converged, diising, 0)
c
#ifdef LOOSE_START
      if (delta.lt.1.d-3)then
#endif
c     
c     Set coulomb acc to max (e.g., input parameter).
c     (note, may also require re-initializing DIIS)
c     
         itol2e = itol_max
         iAOacc = iAOacc_max
         tol_rho = tol_rho_max
         iswitc = iswitc+1
#ifdef LOOSE_START
      endif
#endif
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
         call ga_dadd(one, g_dens(1), one, g_dens(2), g_dens(1))
         call ga_sync
      else
         ndamp = 0
      endif
      call ga_sync
      iter = iter + 1
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
      Etold = Etnew  ! store new energy
c     
      lmaxov = lmaxov_sv
      if ((lumo - homo).lt.-hl_tol.and.lmaxov)then
         lmaxov = .false.
         if (me.eq.0 .and. oprint_conv_details)
     &        write(LuOut,224)homo, lumo
 224     format(10x,' HOMO = ',f6.2,' LUMO = ',f6.2,
     &        /,10x,'Significant orbital reordering with',
     &        ' maximum overlap',
     &        /,10x,'turned ON.  Turning max_ovl OFF.')
      endif
c
      if (oprint_energy_step.and.me.eq.0)then         
         current_cpu = util_cpusec()
         if (nexc.le.1)then
            write(LuOut,222)etnew+enuc, ecore, Ecoul, Exc(1), enuc, 
     &          rho_n, current_cpu
         else
            write(LuOut,223)etnew+enuc, ecore, Ecoul, Exc(1), Exc(2),
     &           enuc, rho_n, current_cpu
         endif
         if (abs(Edisp).gt.0.0d0) 
     &            write(LuOut,225)Edisp
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
c     ----- if cosmo ... after gas_phase, do sol_phase ... -----
c
        if(cosmo_on.and.cosmo_phase.eq.1) then
         converged           = .false.
         egas                = etnew+enuc
         cosmo_phase         = 2
         ocosmo_got_gasphase = .true.
         if(odbug) write(LuOut,*)
     &     'gas_phase done, do sol-phase now ...'
c
c     ----- reset convergence aids -----
c
         iter = 0
         if(diising) then
            call diis_driver_so(toll_s, derr, icall, nfock, 
     &           nbf_mo, g_fockso, g_densso, 
     &           g_svecs, svals, diising, nodiis)
         endif
         damping=.true.
         if(damping) then
         endif
c
         call ga_sync()
         go to 3000
        else if(cosmo_on) then
         esol      =etnew+enuc
        endif
c
c     
c     If converged probably need a few seconds to clean things up 
c     and calculate a few properties.
c     
         wall_time_reqd = 5.0
c     
c        == scale zora eigenvalues and energy ==
         ener_scal = 0.d0
         if (do_zora) then
            call dft_zora_scale_so(
     &                   rtdb,g_dens_at,nexc, ! Added by FA
     &                   geom,
     &                   ao_bas_han,
     &                   nbf,
     &                   nbf_ao,
     &                   nbf_mo,
     &                   g_dens,
     &                   g_s,
     &                   g_moso,
     &                   g_zora_scale_sf,
     &                   g_zora_scale_so,
     &                   dbl_mb(k_eval(1)),
     &                   dbl_mb(k_occ),
     &                   nTotOcc,
     &                   noc, ! FA
     &                   ipol,
     &                   ener_scal)
         end if
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
c
      if (idone.eq.0.or.(iswitc.lt.2.and.iter.lt.iterations))
     &     go to 1000           ! begin new iteration
      if (idone.eq.1.and.(.not.last_time_energy))then
         last_time_energy = .true.
         go to 1000             ! build final total energies
      endif
c     
 2000 continue
c    
      if(cosmo_on.and.cosmo_phase.eq.2) then
          gstote = esol - gcds
          if (do_cosmo_smd) then
           if(.not.rtdb_get(rtdb,'cosmo:gspol',mt_dbl,1,gspol))
     $     call errquit(
     $     'dft_scf: cannot get gspol from rtdb',
     $     0,rtdb_err)
          endif
      endif
c
      if (me.eq.0.and.oprint)then
         if (.not.converged)then
            write(LuOut,*)
            call util_print_centered(LuOut,
     &           'Calculation failed to converge', 20, .true.)
            write(LuOut,*)
         endif
         dft_time = dft_time+util_cpusec()
c
c        Tally up last energy
         if (nexc.eq.1)then
            Etnew = Ecore + Ecoul + exchtot + Exc(1) + Edisp
         else
            Etnew = Ecore + Ecoul + exchtot + Exc(1) + Exc(2) + Edisp
         endif
c
c        Add in cosmo contributions
c
         if (cosmo_on.and.cosmo_phase.eq.2) then
           if (do_cosmo_smd) then
             etnew = etnew + ecosmo + gcds
           else
             etnew = etnew + ecosmo
           end if
         end if
c
         if (nexc.le.1)then
          write(LuOut,222)etnew+enuc,
     &                      ecore,
     &                      ecoul,
     &                      exchtot + exc(1),
     &                      enuc
         else
          write(LuOut,223)etnew+enuc,
     &                      ecore,
     &                      ecoul,
     &                      exchtot + exc(1),
     &                      exc(2),
     &                      enuc
         end if
         if (abs(Edisp).gt.0.0d0) 
     &            write(LuOut,225)Edisp
c
         if (cosmo_on.and.cosmo_phase.eq.2) then
          if (do_cosmo_smd) then
           write(LuOut,226) ecosmo+gcds
          else
           write(LuOut,226) ecosmo
          end if
         endif
c
         if (do_zora) write(luout,2221) ener_scal
         write(luout,2222) rho_n
         write(luout,2223) dft_time
c
         if (me.eq.0.and.oprint) then
          if(cosmo_on.and.cosmo_phase.eq.2) then
           call cosmo_smd_output(egas, esol, do_cosmo_smd, gcds,
     $          gstote, gspol)
          endif
         endif
c
 222  format(//
     &     '      Total SO-DFT energy =', f22.12/
     &     '      One electron energy =', f22.12/
     &     '           Coulomb energy =', f22.12/
     &     '    Exchange-Corr. energy =', f22.12/
     &     ' Nuclear repulsion energy =', f22.12/)
c
 223  format(//
     &     '      Total SO-DFT energy =', f22.12/
     &     '      One electron energy =', f22.12/
     &     '           Coulomb energy =', f22.12/
     &     '          Exchange energy =', f22.12/
     &     '       Correlation energy =', f22.12/
     &     ' Nuclear repulsion energy =', f22.12/)
c
 225  format('    Dispersion correction =', f22.12/)
c
 226  format('             COSMO energy =', f22.12/)
c
 2221 format('       Scaling correction =', f22.12/)
 2222 format(' Numeric. integr. density =', f22.12/)
 2223 format('     Total iterative time =', f9.1,'s'//)
c
         call util_flush(LuOut)
      endif
c
c     print out the determinantal energies 
c
 2001 continue 
      if (me.eq.0.and.oprint.and.det_eng)then
         write(LuOut,*)
c         call util_print_centered(LuOut,
c     &        'Calculation failed to converge', 20, .true.)
         write(LuOut,*)
         dft_time = dft_time+util_cpusec()
         if (nexc.le.1)then
            write(LuOut,232)etnew+enuc, ecore, Ecoul, Exc(1), enuc, 
     &           rho_n, dft_time
            if (do_zora) write(luout,2221) ener_scal
         else
            write(LuOut,233)etnew+enuc, ecore, Ecoul, Exc(1), Exc(2),
     &           enuc, rho_n, dft_time
            if (do_zora) write(luout,2221) ener_scal
         endif
 232     format(//
     &        '       Determinant Energy'/
     &        '       Determinant Energy'/
     &        '      Total SO-DFT energy =', f22.12/
     &        '      One electron energy =', f22.12/
     &        '           Coulomb energy =', f22.12/
     &        '    Exchange-Corr. energy =', f22.12/
     &        ' Nuclear repulsion energy =', f22.12//
     &        ' Numeric. integr. density =', f22.12//
     &        '     Total iterative time =', f9.1,'s'//)
 233     format(//
     &        '       Determinant Energy'/
     &        '      Total SO-DFT energy =', f22.12/
     &        '      One electron energy =', f22.12/
     &        '           Coulomb energy =', f22.12/
     &        '          Exchange energy =', f22.12/
     &        '       Correlation energy =', f22.12/
     &        ' Nuclear repulsion energy =', f22.12//
     &        ' Numeric. integr. density =', f22.12//
     &        '     Total iterative time =', f9.1,'s'//)
         call util_flush(LuOut)
      endif 
c
c     calculate the determinantal energy
c     
      if(.not.fon) goto 2002 
      if (.not.rtdb_get(rtdb,'sodft:ndet',mt_int,1,ndet)) ndet = 0
c
      if(idet .eq. 0 .and. ndet. gt. 0)then 
           if(.not.ma_push_get(mt_int,(nmo_fon+nmo_fon2)*ndet,'det',
     &        kfon_occ,lfon_occ))
     &        call errquit('cannot alloctate lfon_occ',0, MA_ERR)
           if (.not.rtdb_get(rtdb,'sodft:occupancy',mt_int,
     &        (nmo_fon+nmo_fon2)*ndet,int_mb(lfon_occ)))
     &        call errquit('no occupancy specified',0, RTDB_ERR)
       endif  ! (idet.eq.0) and (ndet.gt.0)
c
      if(idet.eq.ndet) goto 2002    ! all determinant energies have been calculated
c
c     Switch the mo order to match the determinantal occupancy 
c
      i_swap = ncore_fon + 1  ! first partially filled orbital
      j_swap = i_swap + idet   ! forms the pair (i,j) to be swapped
      if (j_swap.ge.i_swap) then
         swap(1) = i_swap  ! assign the pair
         swap(2) = j_swap  ! assign the pair
         if(.not.rtdb_put(rtdb,'sodft:swap',mt_int,2,swap))  ! swap is queried in movecs_swap_so
     &        call errquit('swap: failed to put nelem in rtdb', 0,
     &       RTDB_ERR)
         call movecs_swap_so(rtdb,'dft',scftype,g_moso,   ! do the swap
     &        dbl_mb(k_occ),dbl_mb(k_eval(1)))
      end if
c     
c     form a new density matrix and calculate the corresponding energy 
c     
      call dft_densm_so(g_densso,g_moso,nbf_ao,nTotEl) 
      call ga_zero(g_dens(1))
      call ga_zero(g_dens(2))
      call ga_dens_sf(g_dens, g_densso(1), nbf_ao)
      call ga_dadd(one, g_dens(1), one, g_dens(2), g_dens(1))
      call ga_sync
      det_eng = .true. 
c     
c     Restore the mo order so that we have a fixed reference
c     
      if (me.eq.0) write(luout,4443)
 4443 format(/,1x,'Restoring orbital order')
c
      if (j_swap.ge.i_swap) then
         swap(1) = i_swap  ! assign the pair
         swap(2) = j_swap  ! assign the pair
         if(.not.rtdb_put(rtdb,'sodft:swap',mt_int,2,swap))  ! swap is queried in movecs_swap_so
     &        call errquit('swap: failed to put nelem in rtdb', 0,
     &       RTDB_ERR)
         call movecs_swap_so(rtdb,'dft',scftype,g_moso,  ! do the swap
     &        dbl_mb(k_occ),dbl_mb(k_eval(1)))
      end if
c
c     Write determinant configuration     
c
      if(me .eq. 0)then 
         write(Luout,4444)
     &    (int_mb(lfon_occ+idet*(nmo_fon+nmo_fon2)+i-1), 
     &    i=1,nmo_fon+nmo_fon2)
      endif
      idet = idet + 1   ! determinant counter
      go to 1000 
 4444 format(/,1x,'Determinant Occupancy: ',8(1x,i3,1x))
c
      if(fon.and.(ndet.gt.0))then 
       if (.not.ma_pop_stack(kfon_occ))
     &   call errquit('dft_scf: cannot pop stack:lfon_occ',0,MA_ERR)
      endif
c     
 2002 continue

      if (.not. ga_destroy(g_damp_so(1))) call errquit
     &     ('dft_scf_so: could not destroy g_damp_so', 0, GA_ERR)
      if (.not. ga_destroy(g_damp_so(2))) call errquit
     &     ('dft_scf_so: could not destroy g_damp_so', 0, GA_ERR)
c     
c     symmetry adapt vectors last time print symmetries, etc.
c     
c     if (oadapt)then
c     call scf_movecs_sym_adapt(ao_bas_han, g_movecs, oprint, 
c     &                             nbf_ao, '- alpha', .true., 
c     &                             int_mb(k_ir))
c         if (ipol.eq.2)
c     &      call scf_movecs_sym_adapt(ao_bas_han, g_movecs(2), oprint, 
c     &                                nbf_ao, '- beta', .true., 
c     &                                int_mb(k_ir+nbf_ao))
c      endif      
c
c     Vector analysis stolen from rohf.F
c
      if (util_print('final vectors analysis', print_default)) then
         ilo = 1
         ihi = nbf_mo
         call movecs_anal_so(ao_bas_han, ilo, ihi, tanalyze, 
     &        g_moso, 
     &        'DFT Final Molecular Orbital Analysis', 
     &        .true., dbl_mb(k_eval(1)), oadapt, 
     &        int_mb(k_ir), .true., dbl_mb(k_occ))
      endif
c     
c     call to Mulliken Pop Ananlysis
c     
      if (mulliken)then
         if (me.eq.0)
     &      call dft_header
     &      (' Total Density - Mulliken Population Analysis')
         call mull_pop(geom,ao_bas_han,g_dens(1),g_s, 'total')
         if (ipol.eq.2)then
c     
c           analysis of spin density
c     
            if (me.eq.0)call dft_header
     &         (' Spin Density - Mulliken Population Analysis')
            call ga_dadd(one,g_dens(1),-2.d0,g_dens(2),g_dens(2))
            call mull_pop(geom,ao_bas_han,g_dens(2),g_s,'spin')
            call ga_dadd(one,g_dens(1),-1.d0,g_dens(2),g_dens(2))
            call ga_dscal(g_dens(2),0.5d0)
         endif
      endif
c     
c     end infinite loop for SCF iterations
c     
c     Store energy and convergence status ... must store before
c     write movecs since date of insertion is used.
c     
      if (.not. rtdb_put(rtdb, 'sodft:energy',MT_DBL,1,(Etnew+Enuc)))
     &   call errquit('dft_scf: failed to store energy in rtdb', 0,
     &       RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'sodft:converged',MT_LOG,1,converged))
     &   call errquit('dft_scf: failed to store converged in rtdb',0,
     &       RTDB_ERR)
c      if (rtdb_get(rtdb, 'sodft:converged', mt_log, 1, oconverged))
c     &     write(*,*)"converged=", converged
c
c     output energies and eigenvectors to disk
c     
      if (.not.movecs_write_so
     $     (rtdb, ao_bas_han, movecs_out, 'sodft', title,
     &     nbf_mo, dbl_mb(k_occ), dbl_mb(k_eval(1)), g_moso))
     &     call errquit('dft_scf_so: movec_write failed', 0, DISK_ERR)
      call ga_sync()
c
c     Write out the alpha/beta split eigenvalues
c
      call print_evals_ab(nbf_ao, dbl_mb(k_occ), dbl_mb(k_eval(1)))
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
         call dft_mxspin_ovlp(nbf_ao,nmo,ao_bas_han,g_movecs(1), 
     &        g_movecs(2), g_tmp)
      endif
c
      if (wght_GA)then
         if (.not. ga_destroy(g_wght)) call errquit
     &      ('dft_scf: could not destroy g_wght', 0, GA_ERR)
         if (.not. ga_destroy(g_xyz)) call errquit
     &      ('dft_scf: could not destroy g_xyz', 0, GA_ERR)
         if (.not. ga_destroy(g_nq)) call errquit
     &      ('dft_scf: could not destroy g_nq', 0, GA_ERR)
      endif
c     
c     Restore alpha and beta densities.
c
      if (ipol .gt. 1)
     &   call ga_dadd(one,g_dens(1),onem,g_dens(2),g_dens(1))
c     
c     computation of <S2> for open shell
c     
      if (ipol.gt.1)then
         call dft_s2_value(geom,AO_bas_han,.false.,noc(1),noc(2),
     &        nbf_ao,g_dens(1),g_dens(2))
      endif
c     
c     computation of multipole moments
c
      if (natoms .gt. 1)
     &   call dft_mpole(rtdb, ao_bas_han, ipol, g_dens(1), g_dens(2))
c
c     calculate virtual spectrum
c
      lvspec = .false.
      if (.not.rtdb_get(rtdb,'dft:lvspec',mt_log,1,lvspec))
     &         lvspec=.false.
      if (.not. rtdb_get(rtdb,'dft:owstart',mt_int,2,owstart))
     &         lvspec=.false.
      if (.not. rtdb_get(rtdb,'dft:owend',mt_int,2,owend))
     &         lvspec=.false.
      if (lvspec) then
       call util_file_name('vspec',.false.,.false.,vspecfilename)
       call dft_vspec_so(rtdb,ao_bas_han,vspecfilename,nbf_ao,nbf_mo, 
     &    noc, nTotEl, g_moso, k_eval, owstart, owend)
      endif
c     
c     print stolen for uhf.F
c     
      if (util_print('schwarz',print_high).and.(.not.CDFIT))then
         call schwarz_print(natoms, nshells_ao)
      endif
c     
      if (me .eq. 0)then
         if (util_print('final evals', print_high))then
            call util_print_centered(LuOut,'Final eigenvalues',
     &           20,.false.)
            call util_print_centered(LuOut,
     &           '(all occupied plus 10 virtual)',20,.true.)
            call output(dbl_mb(k_eval(1)),
     &           1, min(noc(1)+10,nbf_ao),
     &           1, 1, nbf_ao, 1, 1)
            call util_flush(6)
         endif  ! util_print
c
         if (oprint_final_vecs)then
            write(LuOut,*)
            call util_print_centered(
     &           LuOut,'Final MO vectors',40,.true.)
            write(LuOut,*)
            call util_flush(LuOut)
         endif
      endif
c
c     ECCE printout
c     
      call movecs_ecce(nbf_ao, nmo, 1, nmo(1), dbl_mb(k_eval(1)),
     &                 dbl_mb(k_occ), int_mb(k_ir), 
     &                 g_movecs(1), 'dft', 'alpha')
      if (ipol.eq.2)then ! spin-unrestricted
         call movecs_ecce(nbf_ao, nmo, 1, nmo(2), dbl_mb(k_eval(2)),
     &                    dbl_mb(k_occ+nbf_ao), int_mb(k_ir+nbf_ao), 
     &                    g_movecs(2), 'dft', 'beta')
      endif
      call ecce_print1 ('total energy', mt_dbl, (Etold+Enuc), 1)
      call ecce_print1 ('nuclear repulsion energy', mt_dbl, Enuc, 1)
      call ecce_print1 ('coulomb energy', mt_dbl, Ecoul, 1)
      call ecce_print1 ('exchange energy', mt_dbl, Exc(1), 1)
      if (nexc.gt.1) then
         call ecce_print1 ('correlation energy', mt_dbl, Exc(2), 1)
      endif
c
c     free memory 
      if (ipol.gt.1)then
         if (.not. ga_destroy(g_fockt)) call errquit
     &      ('dft_scf: could not destroy g_fockt', 0, GA_ERR)
      endif
      if (.not. ga_destroy(g_tmp)) call errquit
     &   ('dft_scf: could not destroy g_tmp', 0, GA_ERR)
c
      call fock_2e_tidy(rtdb)
c     
      if (converged)then
         call ecce_print_module_exit('dft', 'ok')
      else
         call ecce_print_module_exit('dft', 'failed')
      endif
c     
c     call ma_summarize_allocated_blocks()
c
      if (.not.ma_pop_stack(lrwork))
     &   call errquit('dft_scf: cannot pop stack:lrwork',0, MA_ERR)
      if (.not.ma_pop_stack(lwork))
     &   call errquit('dft_scf: cannot pop stack:lwork',0, MA_ERR)
      if (.not.ma_pop_stack(lw))
     &   call errquit('dft_scf: cannot pop stack:lw',0, MA_ERR)
      if (.not.ma_pop_stack(la))
     &   call errquit('dft_scf: cannot pop stack:la',0, MA_ERR)
      if (.not.ma_pop_stack(l_ir))
     &   call errquit('dft_scf: cannot pop stack:l_ir',0, MA_ERR)
      if (.not.ma_pop_stack(l_eval))
     &   call errquit('dft_scf: cannot pop stack:l_eval',0, MA_ERR)
      if (.not.ma_pop_stack(l_occ))
     &   call errquit('dft_scf: cannot pop stack:l_occ',0, MA_ERR)
      if (.not.ma_pop_stack(lrdens_atom))
     &   call errquit('dft_scf: cannot pop stack:lrdens_atom',0, MA_ERR)
      if (.not. ma_verify_allocator_stuff())
     $     call errquit('dft_scf: ma corrupt',0, MA_ERR)
      if (.not.ma_pop_stack(lcetobfr))
     &   call errquit('dft_scf: cannot pop stack:lcetobfr',0, MA_ERR)
      if (.not.ma_pop_stack(lcntobfr))
     &   call errquit('dft_scf: cannot pop stack:lcntobfr',0, MA_ERR)
      if (.not.ma_pop_stack(lcntoce))
     &   call errquit('dft_scf: cannot pop stack:lcntoce',0, MA_ERR)
c
      if (CDFIT) then
       if (.not.ma_pop_stack(l_3cwhat))
     &      call errquit('dft_scf: cannot pop stack:l_3cwhat',0, MA_ERR)
       if (.not.ma_pop_stack(l_3cERI))
     &      call errquit('dft_scf: cannot pop stack:l_3cERI',0, MA_ERR)
      endif
c
#if 0
C---------------------------------------------------------------------------
C     KAL
C
C     Dump matrices to file for later comparison.  Used in the real-time
C     TDDFT code (with "matrix_checks" option) for checking purposes.
C
C     NOTE: place just before destroying g_focks in dft_scf.F
C     
      call zmat_compare_dump (g_fockso(1), g_fockso(2), "fock_so_gs")
      call zmat_compare_dump (g_densso(1), g_densso(2), "dens_so_gs")

      if (do_zora) then
         call zmat_compare_dump(g_zora_sf(1), g_zora_sf(2),"zora_sf_gs")
         call zmat_compare_dump(g_zora_so(1), g_zora_so(2),"zora_so_gs")
      endif
C
C     KAL
C---------------------------------------------------------------------------
#endif      

      if(.not.ga_destroy(g_moso(1)))     
     &     call errquit('dft_scf_so: error destroy Movecs Re',0, GA_ERR)
      if(.not.ga_destroy(g_moso(2)))     
     &     call errquit('dft_scf_so: error destroy Movecs Im',0, GA_ERR)
      if(.not.ga_destroy(g_fockso(1)))     
     &     call errquit('dft_scf_so: error destroy Fock Re',0, GA_ERR)
      if(.not.ga_destroy(g_fockso(2)))     
     &     call errquit('dft_scf_so: error destroy Fock Im',0, GA_ERR)
      if(.not.ga_destroy(g_densso(1)))     
     &     call errquit('dft_scf_so: error destroy DenMx Re',0, GA_ERR)
      if(.not.ga_destroy(g_densso(2)))     
     &     call errquit('dft_scf_so: error destroy DenMx Im',0, GA_ERR)
      if(.not.ga_destroy(g_tmp_ri))     
     &     call errquit('dft_scf_so: error destroy old re',0, GA_ERR)
      if(.not.ga_destroy(g_old(1)))     
     &     call errquit('dft_scf_so: error destroy old im',0, GA_ERR)
      if(.not.ga_destroy(g_old(2)))     
     &     call errquit('dft_scf_so: error destroy Tmp ReIm',0, GA_ERR)
      if(.not.ga_destroy(g_so(1)))     
     &     call errquit('dft_scf_so: error destroy so z',0, GA_ERR)
      if(.not.ga_destroy(g_so(2)))     
     &     call errquit('dft_scf_so: error destroy so y',0, GA_ERR)
      if(.not.ga_destroy(g_so(3)))     
     &     call errquit('dft_scf_so: error destroy so x',0, GA_ERR)
c
c     == deallocate zora arrays ==
      if (do_zora) then
c
c      == spin-free parts ==
       if (.not. ga_destroy(g_zora_sf(1))) call errquit(
     &   'dft_scf_so: ga_destroy failed zora sf 1',0, GA_ERR)
       if (.not. ga_destroy(g_zora_sf(2))) call errquit(
     &   'dft_scf_so: ga_destroy failed zora sf 2',0, GA_ERR)
c
       if (.not. ga_destroy(g_zora_scale_sf(1))) call errquit(
     &   'dft_scf_so: ga_destroy failed zora scale sf 1',0, GA_ERR)
       if (.not. ga_destroy(g_zora_scale_sf(2))) call errquit(
     &   'dft_scf_so: ga_destroy failed zora scale sf 2',0, GA_ERR)
c
c      == spin-orbit parts ==
       if(.not.ga_destroy(g_zora_so(1))) call errquit(
     &   'dft_scf_so: ga_destroy failed zora so z',0, GA_ERR)
       if(.not.ga_destroy(g_zora_so(2))) call errquit(
     &   'dft_scf_so: ga_destroy failed zora so y',0, GA_ERR)
       if(.not.ga_destroy(g_zora_so(3))) call errquit(
     &   'dft_scf_so: ga_destroy failed zora so x',0, GA_ERR)
c
       if(.not.ga_destroy(g_zora_scale_so(1))) call errquit(
     &   'dft_scf_so: ga_destroy failed zora so scale z',0, GA_ERR)
       if(.not.ga_destroy(g_zora_scale_so(2))) call errquit(
     &   'dft_scf_so: ga_destroy failed zora so scale y',0, GA_ERR)
       if(.not.ga_destroy(g_zora_scale_so(3))) call errquit(
     &   'dft_scf_so: ga_destroy failed zora so scale x',0, GA_ERR)
c
      end if  !do_zora
c
      dft_scf_so = converged
c
c !!! BGJ
      if (.not. rtdb_get(rtdb, 'bgj:poliz', mt_log,
     &     1, do_poliz)) then
         do_poliz = .false.
      endif
      if (do_poliz) then
c         write(*,*)'*** dft_scf: calling cphf_poliz'
         if (.not. cphf_poliz(rtdb)) ! Never executed.
     $        call errquit(' cphf_poliz: failed from dft_scf !',0,
     &       CALC_ERR)
      endif
c !!! BGJ
      return
c     
 1111 format(15x,'Core Energy:              ',f22.10)
c     
      end
