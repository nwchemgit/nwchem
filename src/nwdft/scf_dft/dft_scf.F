
      logical function dft_scf(rtdb, Etold, Enuc, iVcoul_opt, iVxc_opt, 
     &                         iter, g_dens, g_movecs, g_vxc, g_fock, 
     &                         g_svecs, svals, g_xcinv, g_s)
c     
c     $Id: dft_scf.F,v 1.193 2001-10-18 20:10:19 edo Exp $
c     
      implicit none
c
      integer rtdb              ! [input]
      double precision Etold, Enuc
      integer iVcoul_opt
      integer iVxc_opt
      integer iter
      integer g_dens(2), g_movecs(2), g_vxc(4), 
     &        g_fock, g_svecs, 
     &        g_xcinv
      double precision rho_n, toll_s
      double precision svals(*) ! [in] S evals
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
#ifdef BETTIS
#include "screened_op.fh"
#endif
c     
#if defined(DFT_TIMER)
#include "timer.fh"
#endif
      Logical movecs_write, movecs_converged
      External movecs_write, movecs_converged
      external dft_lindep
      integer dft_lindep
c
      integer  ga_create_atom_blocked
      external ga_create_atom_blocked
c
      integer g_s ! [in]
      logical oprint
      double precision Exc(2), rms(2), derr(2)
      integer noc(2), nmo(2), icall(2)
      integer n3c_dbl, n3c_int, n_batch
      integer iwhat_max
      integer l_3cwhat, k_3cwhat, l_3cERI, k_3cERI
      integer n_semi_bufs, fd
cc    Jorge Garza/Begin
c
c     Variables used to obtain the exact exchange-correlation potential
c
      integer g_rho_exact(2), g_tmp_exact(2), g_addit_exact(2),
     &        g_movecs_aux(2)
      integer incre, exact_pot, split, n_levels(2), act_levels
      double precision lamda, p_lamda, dif_lamda, lamda_old,
     &                 Ecoul_aux(2)
c
      integer test_sic, condfukui, l_degen, i_degen(2)
      double precision totsic
cc    Jorge Garza/End
      Integer l_eval
      integer k_eval(2)
      integer natoms, nTotEl
      integer l_occ, k_occ
      integer i
      integer me, nproc
      integer g_tmp, g_fockt,  g_wght, g_xyz,g_nq
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
     &        lcetobfr, icetobfr, lrdens_atom, irdens_atom
      double precision start_wall, current_wall, elapsed_wall,
     &                 save_wall, current_cpu, start_cpu,
     &                 wall_time_reqd
      integer int_wall_time_reqd
      double precision ecoul, ecore
      double precision pp, delta
      double precision anucl_charg, anel
      double precision  onempp
      double precision etnew, tol2e, tol2e_sleazy,tdots
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
      logical diising, damping, levelshifting,ldmix
      logical keep_damp_on,keep_levl_on, keep_diis_on
      Logical  IOLGC, mulliken
      logical converged, wght_GA
      logical oprint_parm, oprint_conv, oprint_vecs, 
     &        oprint_eval, oprint_syma, oprint_time, 
     &        oprint_info, oprint_tol, oprint_final_vecs, 
     &        oprint_energy_step, oprint_intermediate_fock,
     &        oprint_3c2e, oprint_interm_overlap, oprint_interm_S2,
     &        oprint_conv_details, oprint_sic,oprintinertia
      double precision zero, onem, one
      parameter(zero = 0.d0, one = 1.d0, onem = -one)
c
      integer ilo, ihi          ! For printing movecs analysis
      integer iter_by_5, iter_ck
      double precision eval_pr_tol_lo, eval_pr_tol_hi
      parameter (eval_pr_tol_lo = -1.5d0, eval_pr_tol_hi=0.5)
C
c
c     early convergence tolerances
c
      parameter(itol_min = 7, iAOacc_min = 12, tol_rho_min = 1.d-7)
c
      double precision dft_dencvg, dft_time
      external dft_dencvg
      double precision homo, lumo, homo_lumo_gap
      double precision ssmear
      integer l_ir, k_ir
      logical last_time_energy
      logical check_shift, lmaxov_sv,fon, reorth
      character*7 name
      character*4 scftype
      character*7 vecs_or_dens
      character*255 basis_name, basis_trans,blob

      integer nopen, nclosed, g_denso(2),k_tmp,l_tmp,nstdis,
     .     l_gdiis,k_gdiis,ipolmod
      logical dft_mem3c
      external dft_mem3c
c !!! BGJ
      logical cphf_poliz, do_poliz
      external cphf_poliz
c !!! BGJ
c
c     ----- cosmo -----
c
      logical odbug
      logical ocosmo
      logical ogas_phase
      logical osol_phase
      double precision egas
      double precision esol
      double precision ecosmo,dft_cosmoe,dft_cosmoadde
      external dft_cosmoe,dft_cosmoadde


      integer g_focks(2)
      integer nel_fon(2),nmo_fon(2),g_movad,g_xij,
     ,     ncanorg
c
      call ecce_print_module_entry('dft')
      dft_scf = .false.
      lmaxov_sv = lmaxov
      oprint = util_print('information', print_low)
      oprint_info = util_print('common', print_debug)
      oprint_parm = util_print('parameters', print_default)
      oprint_3c2e = util_print('3c 2e integrals', print_default)
      oprint_conv = util_print('convergence', print_default)
      oprint_conv_details = util_print('convergence details', 
     &                                 print_high)
      oprint_vecs = util_print('intermediate vectors', print_high)
      oprint_eval = util_print('intermediate evals', print_high)
      oprint_syma = util_print('interm vector symm', print_high)
      oprint_time = util_print('dft timings', print_high)
      oprint_tol = util_print('screening parameters', print_high)
      oprint_energy_step = util_print('intermediate energy info',
     &                                print_high)
      oprint_intermediate_fock = util_print('intermediate fock matrix',
     &                                      print_high)
      oprint_interm_S2 = util_print('intermediate S2',print_high)
      oprint_interm_overlap = util_print('intermediate overlap',
     &                                      print_high)
      oprint_final_vecs = util_print('final vectors', print_high)
      oprint_sic = util_print('SIC information', print_high)
      oprintinertia = util_print('inertia', print_default).and.
     . ga_nodeid().eq.0
      odbug = util_print('cosmo', print_high)

      me = ga_nodeid()
      nproc = ga_nnodes()
c
c !!! BGJ
c     Store SCF hamiltonian type as DFT for use in BGJ routines
      if (.not. rtdb_put(rtdb, 'bgj:scf_type', MT_INT, 1, 2))
     $     call errquit('dft_scf: put of bgj:scf_type failed',0)
c !!! BGJ
c
c     check if smearing is on
c
      ssmear=0.d0
      nstdis=1
      if (rtdb_get(rtdb, 'dft:smear_sigma', mt_dbl, 1,
     &   ssmear))then
      endif
      if (rtdb_get(rtdb, 'dft:fon', mt_log, 1,
     &   fon))then
      if (.not.rtdb_get(rtdb, 'dft:nmo_fon', mt_int, 2,
     & nmo_fon)) call errquit('dft_scf: rtdb_gett failed', 2300)
      if (.not.rtdb_get(rtdb, 'dft:nel_fon', mt_int, 2,
     & nel_fon)) call errquit('dft_scfp: rtdb_get failed', 2301)
      else
         fon=.false.
      endif
      if (.not.rtdb_get(rtdb, 'dft:reorth', mt_log, 1,
     &   reorth))reorth=.false.
c
c     see if levelshifting monitoring is desired
c
      if (.not. rtdb_get(rtdb, 'dft:check_shift', mt_log, 1,
     &   check_shift))then
         check_shift = .false.      
      endif
c     
      if (.not. geom_ncent(geom, natoms))
     &   call errquit('dft_scf: geom_ncent failed',73)     
      if (.not. geom_nuc_charge(geom, anucl_charg))
     &   call errquit('dft_scf: geom_nuc_charge failed', 0)
c     
      anel = int(anucl_charg) - rcharge
c
c     Pre-compute mapping vectors
c 
      if (.not.ma_push_get
     &   (mt_int,nshells_ao,'cntoce map',lcntoce,icntoce))
     &   call errquit('dft_scf:push_get failed', 13)
      if (.not.ma_push_get
     &   (mt_int,nshells_ao*2,'cntoce map',lcntobfr,icntobfr))
     &   call errquit('dft_scf:push_get failed', 13)
      if (.not.ma_push_get
     &   (mt_int,natoms*2,'cntoce map',lcetobfr,icetobfr))
     &   call errquit('dft_scf:push_get failed', 13)
c     
      call build_maps(ao_bas_han, int_mb(icntoce), int_mb(icntobfr), 
     &                int_mb(icetobfr), natoms, nshells_ao)

      if (.not. rtdb_get(rtdb, 'dft:noc', mt_int, 2, noc))
     &   call errquit('dft_scf: rtdb_get of noc failed', 0)
      if (.not. MA_Push_Get(MT_Dbl,nbf_ao*ipol,'mo occ',l_occ,k_occ))
     &   call errquit('dft_scf: failed to alloc',999)
c
c     get/set orbital overlap tolerance
c
      if (.not. rtdb_get(rtdb, 'dft:toll_s', MT_DBL, 1, toll_s))
     .     call errquit('dftscf: lost toll_s ',0)
      nmo(1)=dft_lindep(nbf_ao,toll_s,
     .     svals,me.eq.0,luout)
      if (ipol.eq.1)then
c     
c        noc(1) = No. of electrons alpha
c     
         nTotEl = 2*noc(1)
      else
c     
         nTotEl = noc(1) + noc(2)
         nmo(2) = nmo(1)
      endif
c
c    Jorge Garza/Begin
      call xc_sicinit(rtdb, test_sic, condfukui, exact_pot,
     &                l_degen, i_degen, noc, act_levels)
c     Jorge Garza/End
c
c     Set aside some memory for reduced density matrix
c
      ipolmod=ipol
      if(test_sic.ne.0) ipolmod=2
      if (.not.MA_Push_Get(MT_Dbl,ipolmod*natoms*natoms,'rdens_atom',
     &   lrdens_atom,irdens_atom))
     &   call errquit('dft_scf: cannot allocate rdens_atom',0)
c     
c     determine pattern of orbitals' occupancy
c     
c      
c
      if (ipol .eq. 1)then
c
c        RHF occupations
c
         call dfill(nbf_ao, 0.0d0, dbl_mb(k_occ), 1)
         do i = 1, noc(1)
            dbl_mb(i-1+k_occ) = 2.0d0
         enddo
      else  
c
c        UHF occupations
c
         call dfill(2*nbf_ao, 0.0d0, dbl_mb(k_occ), 1)
         do i = 1, noc(1)
            dbl_mb(i-1+k_occ) = 1.0d0
         enddo
         do i = nbf_ao+1, nbf_ao+noc(2)
            dbl_mb(i-1+k_occ) = 1.0d0
         enddo
      endif 
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
      if (oprint_time)
     &     call dft_tstamp(' Before 3c-2e initialize.')
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
      g_focks(1)=g_fock
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
      if (nodiis)then
         diising = .false.
      else
         if (.not.MA_alloc_Get(MT_int,nfock*ipol*2,
     ,        'ga handle for diis',l_gdiis,k_gdiis))
     ,        call errquit('dft_scf: cannot allocate g_diis',0)
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
c     allocate eigenvalue array, including second pointer to beta block
c
      leneval = nbf_ao * ipol
      if (.not.MA_Push_Get(MT_Dbl,leneval,'eval',l_eval,k_eval(1)))
     &     call errquit('dft_scf: cannot allocate eval',0)
      if (ipol .eq. 2)then  
         k_eval(2) = k_eval(1) + nbf_ao
      endif 
c     
c     Dump DFT parameters (if debugging) to see if they make sense
c     
      if (me.eq.0.and.oprint_info)call dft_dump_info(me)

      if (ipol.eq.1)then
         scftype = 'RHF'
      elseif (ipol.eq.2)then
         scftype = 'UHF'
      endif
c
c     allocate array for irreps
c
      if (.not.MA_Push_Get(mt_int,ipol*nbf_ao,'dft:irreps',l_ir,k_ir))
     &   call errquit('dft_scf: cannot allocate irreps',0)
      nopen = mult - 1
      nclosed = (nTotEl - nopen) / 2
c
      if (.not. bas_name(ao_bas_han, basis_name, basis_trans))
     $     call errquit('dft_scf: bas_name?', 0)
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
      vecs_or_dens='density'
      ldmix=.false.
      ncanorg=0
      if((movecs_in.ne.'atomic'.and.
     *     movecs_in.ne.'fragment'.and.movecs_in.ne.'hcore')
     .     .or.fon) then
         vecs_or_dens='vectors'
         ldmix=.true.
         ncanorg=1
      else
         if (.not. ga_create(mt_dbl, nbf_ao, nbf_ao, 'Xij canorg', 
     &        nbf, 0, g_xij)) 
     .        call errquit(' ddftdiagn: ga_create failed ',0)
      endif
      call scf_vectors_guess(rtdb, tol2e_sleazy, geom, ao_bas_han, 
     &                       basis_trans, movecs_in, movecs_out, 
     &                       movecs_guess, scftype, nclosed, nopen, 
     &                       nbf, nbf, noc(1), noc(2),  k_eval, k_occ, 
     &                       k_ir, g_movecs, g_dens, vecs_or_dens, 
     &                       'dft', title, oskel, oadapt, 
     &                       .true.) 
      if(oskel) then
        call ga_symmetrize(g_dens(1))
        if(ipol.eq.2) call ga_symmetrize(g_dens(2))
      endif


      if(ldmix) then
         tdots=0.d0
         call dft_densm(g_dens, g_movecs, 
     &        nbf_ao, ipol, 
     &        geom, AO_bas_han,noc,ntotel,
     &        Dbl_MB(k_eval(1)),dbl_mb(k_occ),
     &        ssmear,tdots,iter,.false.,
     .        fon, nel_fon,nmo_fon)
c
c  Jorge Garza/Begin/Exact_pot
c
         if (exact_pot.eq.1) then
           call open_xc_exact_pot(g_rho_exact, g_tmp_exact,
     &                            g_addit_exact, g_movecs_aux,
     &                            g_dens, dif_lamda, lamda,
     &                            p_lamda, lamda_old, Ecoul_aux)
         end if
c
c Jorge Garza/End/Exact_pot

      call dft_prevals('inter',me,oprint_eval,oprint_vecs,noc,
     ,     g_movecs,dbl_mb(k_eval(1)),dbl_mb(k_occ))

         endif
c
c     Tidy up SCF
c
      call fock_2e_tidy(rtdb)
#if defined(DFT_TIMER)
      call dft_second(tp3)
#endif
c
      if (CDFIT)then
         if(dft_mem3c(
     I     natoms,npol,oprint_parm,oprint_3c2e,
     O     n3c_int,n3c_dbl,n_semi_bufs,
     O     l_3ceri,k_3ceri, l_3cwhat,k_3cwhat)) then
            call dft_3cincor(n_batch, n3c_int, int_mb(k_3cwhat), 
     &                       dbl_mb(k_3cERI), n3c_dbl, iwhat_max, 
     &                       n_semi_bufs, fd)
            incore=.true.
         else
            if (me.eq.0 .and. oprint_3c2e)write(LuOut,3230)
            incore=.false.
         endif
      endif
 3230 format(/,10x,'Incore memory use for 3-center 2e- integrals is ',
     &     'turned off. ')
#if defined(DFT_TIMER)
      call dft_second(tp4)
      pgmsect(4) = ' 3c-2e Ints'
      call tcllct(t_pgm(1,4),tp4,tp3)
      call dft_second(tp3)
#endif
c
c     set initial coulomb acc
c
c      write(LuOut,*)' movecs_guess = ',movecs_guess
      if (movecs_guess.eq.'restart')ltight=.true.
c     
c     May not want levelshifting initially until sure that the
c     transformed Fock matrix will be diagonally dominant, or
c     alternatively shift the piss out of it.
c     
      if (movecs_guess.eq.'restart'.or.ncanorg.eq.1)then
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
#ifdef BETTIS
c     Screened exchange works only with direct SCF
      if (screened_exchange()) then
         write(*,*)'Forcing direct SCF for screened exchange'
         call fock_force_direct(rtdb)
      endif
#endif
      tol2e = 10.d0**(-itol_max)
c 
c     Restore SCF parameters
c
      call scf_get_fock_param(rtdb, tol2e)
c     
c     If open shell put the total density matrix in g_dens(1)
c
      if(ipol.eq.2)then             
         call ga_dadd(one,g_dens(1),one,g_dens(2),g_dens(1))
      endif
c     
c     Call to Mulliken Pop Ananlysis for initial density
c     
      if (mulliken) call dft_mulwrap(me,g_dens,g_s)
c
      call dft_cosmoinit(rtdb,ocosmo,odbug)
c
c    The loop associated with the label 135 is to find the exact exchange-
c    correlation potential.
c
 135  continue
c     
c     Top of infinite SCF iteration loop
c
c     Write prep time required
c
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
      call dft_cosmochk(rtdb,ocosmo,odbug,ogas_phase,osol_phase)
c
      iter = 1
c
      if (ipol.gt.1)then
         g_fockt = ga_create_atom_blocked(geom, AO_bas_han, 'fock tr')
      endif
 3000 continue
      if (ipol.gt.1)then
         if(.not.ga_duplicate(g_fock, g_focks(2), ' beta fock'))
     .       call errquit(' dftscf:failed duplicate ',g_fock)
      endif
c
c     ----- dft scf proper -----
c
      last_time_energy = .false.
c
      g_denso(1) = 
     =     ga_create_atom_blocked(geom,AO_bas_han,'oldDMa')
      if(ipol.eq.2)
     +     g_denso(2) = 
     =     ga_create_atom_blocked(geom,AO_bas_han,'oldDMb')
#if defined(DFT_TIMER)
      call dft_second(tp4)
      pgmsect(5) = ' initial guess'
      call tcllct(t_pgm(1,5),tp4,tp3)
#endif
c
      if(iterations.eq.0) then
         if (.not. ga_destroy(g_denso(1))) call errquit
     &        ('dft_scf: could not destroy g_denso', 0)
         if(ipol.eq.2) then
            if (.not. ga_destroy(g_focks(2))) call errquit
     &           ('dft_scf: could not destroy g_focks2', 0)
            if (.not. ga_destroy(g_denso(2))) call errquit
     &           ('dft_scf: could not destroy g_dens2', 0)
         if (.not. MA_free_heap(l_gdiis))
     .        call errquit(' dft_scf:cannot popstack',111)
         endif
         goto 1970
      endif
 1000 continue

c
      if (me.eq.0 .and. oprint_conv_details)
     &   write(LuOut,124)damping, levelshifting, diising
 124  format(10x,' DAMPING=',l1,' LEVELSHIFTING=',l1,
     &           ' DIISING=',l1)
c
      if (me.eq.0.and.oprint_tol)write(LuOut,3234)itol2e,iAOacc,iXCacc
 3234 format(10x,'itol2e=',i2,' iAOacc=',i2,' iXCacc=',i2)

      Ecoul  = ZERO
      Exc(1)= ZERO
      Exc(2) = ZERO
      rms(1) = 0.d0
      rms(2) = 0.d0
      ecosmo=dft_cosmoe(rtdb,ocosmo,odbug,
     .     geom,iter,g_movecs,ao_bas_han)

c
c     build DFT Fock matrix
c
      call dft_fockbld(rtdb,natoms,ntotel,
     ,     g_fock,g_dens,g_vxc,g_movecs,g_xcinv,
     ,     g_tmp,g_nq,g_wght,g_xyz,
     ,     iwhat_max,ivcoul_Opt,nexc,ivxc_opt,wght_GA,
     ,     n_batch,n3c_int,n3c_dbl,n_semi_bufs,fd,IOLGC,
     ,     k_eval, dbl_mb(irdens_atom),
     ,     dbl_mb(k_3ceri),int_mb(k_3cwhat),
     ,     int_mb(icetobfr),
     ,     iter, ecore, ecoul,exc,noc,rho_n,ldmix,
     ,     test_sic,exact_pot,ecoul_aux,dif_lamda,totsic,
     ,     i_degen,g_rho_exact,g_tmp_exact, g_addit_exact,
     ,     n_levels, act_levels, oprint_time)
c     
c     Calculate the total electronic energy.
c     
      if (nExc.eq.1)then
         Etnew = Ecore + Ecoul + Exc(1)
      else
         Etnew = Ecore + Ecoul + Exc(1) + Exc(2)
      endif
      if (oprint_conv_details.and.me.eq.0)then
         write(LuOut,*)'Etnew, Ecore, Ecoul, Exc(1), Exc(2), Totsic: ',
     &                  Etnew, Ecore, Ecoul, Exc(1), Exc(2), Totsic
      endif
      etnew=dft_cosmoadde(rtdb,ocosmo,odbug,
     .     ecosmo,etnew)
c     
      if (last_time_energy)then
         etnew=etnew+tdots
c     
c        If open shell put the total density matrix back in 
c        g_dens(1) and quit.
c     
         if (ipol.eq.2)then
            call ga_dadd(one, g_dens(1), one, g_dens(2), g_dens(1))
         endif
         goto 2000
      endif
      homo_lumo_gap = 200.0d0
      if (ipol.gt.1) then
        call ga_copy(g_fock, g_fockt)
      endif
      do ispin = 1, ipol
c     
c        Add V xc to Fock Matrix
c     
         if (ispin.eq.2) then
            call ga_copy(g_fockt, g_focks(ispin))
            call ga_dadd(one, g_vxc(ispin), one, g_focks(ispin), 
     .           g_focks(ispin))
         endif
c     
#if defined(DFT_TIMER)
         call dft_second(tp3)
#endif
c     
c
c        Symmetrize the Fock matrix
c
         if (oskel)
     &      call sym_symmetrize(geom, AO_bas_han, .false., 
     .        g_focks(ispin))
c
         call ga_symmetrize(g_focks(ispin))
      enddo
         if (diising)then
c     
c           DIIS step taken here.
c     
            if (oprint_time)
     &           call dft_tstamp(' calling diis driver ')
            nstdis=nstdis+1
            call diis_driver(toll_s, derr,  ipol, 
     &           icall, nfock, nbf_ao, geom, 
     &           ao_bas_han, g_focks, g_dens, 
     &           g_tmp, g_svecs, svals, diising, 
     &           nodiis,int_mb(k_gdiis),lkeeps,
     .           g_s12m,g_s12p)
            if (oprint_time)
     &           call dft_tstamp(' called diis driver ')
         endif
c     
#if defined(DFT_TIMER)
         call dft_second(tp4)
         pgmsect(9) = '   DIIS   '
         call tcllct(t_pgm(1,9),tp4,tp3)
#endif
c     
      do ispin = 1, ipol
#if defined(DFT_TIMER)
         call dft_second(tp3)
#endif
         call dft_diagn(levelshifting,
     ,        nmo(ispin),ncanorg,g_xij,
     .        g_focks(ispin),g_s,g_movecs(ispin),g_tmp,g_svecs,
     ,        dbl_mb(k_eval(ispin)),svals,noc(ispin),
     ,        homo,lumo,toll_s,
     ,        oprint_intermediate_fock,oprint_time,
     ,        oprint_conv_details)

c        
c        determine homo-lumo gap 
c
         homo_lumo_gap = min(homo_lumo_gap, (lumo-homo-rlshift))
         if (me.eq.0 .and. oprint_conv_details)
     &      write(LuOut,4224)homo,lumo,rlshift, homo_lumo_gap
c     
         call ga_sync
#if defined(DFT_TIMER)
         call dft_second(tp4)
         pgmsect(11) = '   Diag   '
         call tcllct(t_pgm(1,11),tp4,tp3)
#endif
c     
c        Save previous density for convergence check.
c     
         call ga_copy(g_dens(ispin), g_denso(ispin))
c
c        symmetry adapt vectors?
c
         if (oadapt)then
            if(nmo(ispin).ne.0) then
               if (.not. ga_create(mt_dbl, nbf_ao, nmo(ispin),
     ,              'movad', nbf_ao, 0, g_movad)) call errquit(
     '              'dft_scf: gacreate failed ',0)
               call ga_copy_patch('n', g_movecs(ispin), 1, nbf_ao, 
     .              1, nmo(ispin),
     ,              g_movad, 1,nbf_ao,1,nmo(ispin))
            else
               g_movad=g_movecs(ispin)
            endif
            if(ispin.eq.1)name = '- alpha'
            if(ispin.eq.2)name = '- beta'
            call scf_movecs_sym_adapt(ao_bas_han, g_movad,
     &           oprint_syma, nmo, name,
     &           .true., 
     &                                int_mb(k_ir+nbf_ao*(ispin-1)))
            if(nmo(ispin).ne.0) then
               call ga_copy_patch('n', 
     ,              g_movad, 1,nbf_ao,1,nmo(ispin),
     .              g_movecs(ispin), 1, nbf_ao, 1, nmo(ispin))
               if (.not. ga_destroy(g_movad)) call errquit(
     &              'dft_scf: could not destroy g_movad', 0)
            endif
         endif      
c
c        occasionally (every 5 iterations) clean up (by orthogonalizing
c        away any accumulated noise) the mo vectors.
c
         iter_by_5 = iter/5
         iter_ck = iter - 5*iter_by_5
         if (iter_ck.eq.0.and.reorth)
     &      call ga_orthog_mos(ao_bas_han, g_movecs(ispin))

      enddo                     ! end big loop over ispin

      call ga_sync
c     
c     Form a new density matrix.
c     
      tdots=0.d0
      if (oprint_time)
     &     call dft_tstamp(' calling densm ')
      call dft_densm(g_dens, g_movecs, 
     &     nbf_ao, ipol, 
     &     geom, AO_bas_han,noc,ntotel,
     &      Dbl_MB(k_eval(1)),dbl_mb(k_occ),
     &     ssmear,tdots,iter,.true.,
     .        fon, nel_fon,nmo_fon)

      if (oprint_time)
     &     call dft_tstamp(' called densm ')
      
      if(ssmear.gt.1d-9.and.ipol.eq.2.and.(noc(2).gt.noc(1))) then
c
c     swap alpha and beta
c
      if (.not.MA_Push_Get(MT_Dbl,nbf_ao,'ssccrr',l_tmp,k_tmp))
     &     call errquit('dft_scf: cannot allocate ktmp',0)
         call dft_swapab(g_dens, g_movecs,nbf_ao,noc,
     ,        dbl_mb(k_eval(1)),dbl_mb(k_occ),g_tmp,dbl_mb(k_tmp))
      if (.not.ma_pop_stack(l_tmp))
     &   call errquit('dft_scf: cannot pop stack',33)
            if(me.eq.0) write(LuOut,*) ' Swapped AB - resetting Diis'
            if(diising) then
            icall(1) = -1
            nstdis=1
            call diis_driver(toll_s, derr,  ipol, icall, nfock, 
     &                       nbf_ao, geom, ao_bas_han, g_focks, 
     &                       g_dens, g_tmp, g_svecs, svals, 
     &                       diising, nodiis,int_mb(k_gdiis),
     ,           lkeeps,g_s12m,g_s12p)
            icall(1) = 0
            icall(2) = 0
            endif
      endif
         etnew=etnew+tdots
c
      delta = abs(etold-etnew)
c
c     check for symmetry breaking
c
      if(oskel.and.ipol.eq.2) then
          call dft_symbr(geom,ao_bas_han,
     .     g_dens(1),g_dens(2), g_tmp)
       endif
c

      do ispin = 1, ipol
c     
c     
c     Check convergence on Density.
c     
         rms(ispin) = dft_dencvg(g_dens(ispin), g_denso(ispin), nbf_ao)
c     
c        Damping implemented here.
c     
         if (damping)then
            pp = dble(ndamp)*1.d-2
            onempp = 1.0d0 - pp
            call ga_dadd(pp, g_denso(ispin),onempp, g_dens(ispin), 
     &                   g_dens(ispin))
         else
            ndamp = 0
         endif
      enddo                     ! end loop over ispin
      call ga_sync
c     
      if (oprint_conv.and.iter.eq.1.and.me.eq.0)then
c
         if(ogas_phase.or.osol_phase) then
c
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
         endif
c
c     ----- cosmo message -----
c
         if(ogas_phase.and.osol_phase) then
            write(LuOut,909)
         elseif(.not.ogas_phase.and..not.osol_phase) then
            write(LuOut,910)
         endif
 909     format(2x,'   gas_phase')
 910     format(2x,'   sol_phase')
c
      endif
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
      if (.not.movecs_write(rtdb, ao_bas_han, movecs_out, 'dft', title,
     &                      nbf_ao, ipol, nmo, dbl_mb(k_occ), nbf_ao, 
     &                      dbl_mb(k_eval(1)), nbf_ao, g_movecs))
     &   call errquit('dft_scf: movec_write failed', 0)
      call dft_prevals('inter',me,oprint_eval,oprint_vecs,noc,
     ,     g_movecs,dbl_mb(k_eval(1)),dbl_mb(k_occ))
c     
c     If open shell compute overlap of alpha orbitals with beta 
c     orbitals.
c     
      if ((ipol.gt.1).and.(oprint_interm_overlap)) then
         call dft_mxspin_ovlp(ao_bas_han, g_movecs(1), g_movecs(2), 
     &                        g_tmp)
      endif
c     
c     computation of <S2> for open shell
c     
      if ((ipol.gt.1).and.(oprint_interm_S2)) then

         call dft_s2_value(geom, AO_bas_han, .false., noc(1), noc(2),
     &                     nbf_ao, g_dens(1), g_dens(2))
      endif
c     
c     
c     Form the total density matrix.
c     
      if (ipol.eq.2)then
         call ga_dadd(one, g_dens(1), one, g_dens(2), g_dens(1))
      endif
c     
c     Check for SCF convergence.
c     
      call ga_sync
      if (oprint_time)
     &     call dft_tstamp(' calling scfcvg ')
      call dft_scfcvg(rms, derr, Etold, Etnew,
     &                e_conv, d_conv, g_conv, ipol, 
     &                iter, iterations, idone, rtdb,
     &                converged, diising)
      if (oprint_time)
     &     call dft_tstamp(' called scfcvg ')
      if (abs(delta).lt.1.d-3)then
c     
c        Set coulomb acc to max (e.g., input parameter).
c        (note, may also require re-initializing DIIS)
c     
         itol2e = itol_max
         iAOacc = iAOacc_max
         tol_rho = tol_rho_max
         iswitc = iswitc+1
      endif
c     
      iter = iter + 1
cedo      if(ssmear.gt.1d-9.and.mod(nstdis,30).eq.0) then
cedo            if(me.eq.0) write(LuOut,*) ' resetting Diis',nstdis
cedo            icall(1) = -1
cedo            nstdis=1
cedo            call diis_driver(toll_s, derr(1), 1, ipol, icall, nfock, 
cedo     &                       nbf_ao, geom, ao_bas_han, g_fock, 
cedo     &                       g_dens(1), g_tmp, g_svecs, isvals, 
cedo     &                       diising, nodiis,int_mb(k_gdiis),
cedo     ,           lkeeps,g_s12m,g_s12p)
cedo            icall(1) = 0
cedo            icall(2) = 0
cedo         endif
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
     &    (iter.le.ncysh))then
         if (homo_lumo_gap.lt.hl_tol)then
            levelshifting = .true.
            rlshift = rlshift_input
            if (check_shift)then
               if (lumo .lt. homo)then
                  levelshifting = .false.
                  if (me.eq.0 .and. oprint_conv_details)
     &               write(LuOut,2224)homo, lumo
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
      if ((delta.lt.diison.and.delta.gt.diisoff).or.
     &     iter.le.ncyds.or.keep_diis_on)then
         diising = .true.
c     
c        Once started, keep DIIS on until diisoff threshold.
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
     &      write(LuOut,224)homo, lumo
      endif
c     
      if (oprint_energy_step.and.me.eq.0)then         
         current_cpu = util_cpusec()
         if (nexc.le.1)then
            write(LuOut,222)etnew+enuc, ecore, Ecoul, Exc(1), enuc, 
     &                      rho_n, current_cpu
         else
            write(LuOut,223)etnew+enuc, ecore, Ecoul, Exc(1), Exc(2),
     &                      enuc, rho_n, current_cpu
         endif
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
c        If converged probably need a few seconds to clean things up 
c        and calculate a few properties.
c
         wall_time_reqd = 5.0
      else
c
c        If not converged probably need at least the amount time
c        required for previous iteration (multiply by 1.2 to be on the safe side).
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
      if (idone.eq.0.or.(iswitc.lt.2.and.iter.lt.iterations).or.
     .     (iterations.eq.0))
     &   go to 1000 ! begin new iteration
      if (idone.eq.1.and.(.not.last_time_energy))then
         last_time_energy = .true.
         go to 1000             ! build final total energies
      endif
c     
 2000 continue
c
c     get rid of levelshifting in evals
c
         call dft_remshft(levelshifting,ipol,nbf_ao,noc,
     ,        rlshift,dbl_mb(k_eval(1)))
      
      Etnew = Ecore + Ecoul + Exc(1) +tdots
      etnew=dft_cosmoadde(rtdb,ocosmo,odbug,
     .     ecosmo,etnew)
c    
      if (.not. ga_destroy(g_denso(1))) call errquit
     &   ('dft_scf: could not destroy g_denso', 0)
      if(ipol.eq.2) then
          if (.not. ga_destroy(g_denso(2))) call errquit
     &     ('dft_scf: could not destroy g_denso', 0)
          if (.not. ga_destroy(g_focks(2))) call errquit
     &     ('dft_scf: could not destroy g_focks2', 0)
      end if
c
      if (exact_pot.eq.1) then
        call parlam_xc_exact_pot(lamda, lamda_old, dif_lamda,
     &                               p_lamda, incre, g_tmp_exact,
     &                               g_addit_exact, g_movecs_aux,
     &                               g_movecs, iter, k_eval,
     &                               Ecoul_aux, noc, split)
        if (lamda.le.900.0) go to 135
        call close_xc_exact_pot(g_rho_exact, g_tmp_exact,
     &                          g_addit_exact, g_movecs_aux)
      end if
c
c     Jorge Garza / Begin
      if (test_sic.eq.1) then
         call xc_sicdeg(i_degen, k_eval, noc, n_levels,act_levels)
         call localize_sic(g_movecs, noc, k_eval)
         call xc_sic_pert(rtdb, nExc, iVxc_opt,
     &                    g_wght, g_xyz, g_nq,
     &                    wght_GA, dbl_mb(irdens_atom),
     &                    int_mb(icetobfr), natoms,
     &                    g_movecs, noc, totsic, i_degen, n_levels)
         if (me.eq.0.and.oprint_sic)then
            write(LuOut,*) ' SIC after localization',totsic
         endif
         Etnew = Etnew + totsic
         Exc(1) = Exc(1) + totsic
      end if
c
c     Jorge Garza / End
c     

c
c     ----- if cosmo ... after gas_phase, do sol_phase ... -----
c
      if(ogas_phase.and.osol_phase) then
         converged =.false.
         egas      =etnew+enuc
         ogas_phase=.false.
         osol_phase=.false.
         if (.not. rtdb_put(rtdb,'slv:cosmo',mt_log,1,.true.))
     $       call errquit('dft_scf: rtdb_put failed',0)
         if(odbug) then
            write(LuOut,*) 'gas_phase done, do sol-phase now ...'
         endif
c
c     ----- reset convergence aids -----
c
         iter  =1
         if(diising) then
            icall(1) = -1
            nstdis=1
            call diis_driver(toll_s, derr,  ipol, icall, nfock,
     &                       nbf_ao, geom, ao_bas_han, g_focks,
     &                       g_dens, g_tmp, g_svecs, svals,
     &                       diising, nodiis,int_mb(k_gdiis),
     ,           lkeeps,g_s12m,g_s12p)
            icall(1) = 0
            icall(2) = 0
         endif
         damping=.true.
         if(damping) then
         endif
c
         call ga_sync()
         go to 3000
      elseif(.not.ogas_phase.and..not.osol_phase) then
         esol      =etnew+enuc
      endif
      call ga_sync()
c
      if (me.eq.0.and.oprint)then
         if (.not.converged)then
            write(LuOut,*)
            call util_print_centered(LuOut,
     &           'Calculation failed to converge', 20, .true.)
            write(LuOut,*)
         endif
         dft_time = dft_time+util_cpusec()
         if (nexc.le.1)then
            write(LuOut,222)etnew+enuc, ecore, Ecoul, Exc(1), enuc, 
     &           rho_n, dft_time
         else
            write(LuOut,223)etnew+enuc, ecore, Ecoul, Exc(1), Exc(2),
     &           enuc, rho_n, dft_time
         endif
         if(ocosmo) call dft_cosmosum(egas,esol)
      endif
 1970 continue
c
c     Vector analysis stolen from rohf.F
c
      do ispin = 1, ipol
         if (util_print('final vectors analysis', print_default)) then
            do ilo = 1,max(1,nclosed-10)
               if (dbl_mb(k_eval(ispin)+ilo-1) .ge. eval_pr_tol_lo) 
     &            goto 961
            enddo
 961        do ihi = min(nclosed+nopen+10,nbf_ao), nbf_ao
               if (dbl_mb(k_eval(ispin)+ihi-1) .ge. eval_pr_tol_hi) 
     &            goto 9611
            enddo
            ihi = max(ihi-1,1)
 9611       continue
            if (util_print('final vectors analysis', print_high)) then
               ilo = 1
               ihi = nbf_ao
            endif
            if (ipol.eq.1) then
               blob='DFT Final Molecular Orbital Analysis'
            else if(ipol.eq.2.and.ispin.eq.1) then
               blob='DFT Final Alpha Molecular Orbital Analysis' 
            else
               blob='DFT Final Beta Molecular Orbital Analysis' 
            endif
            call movecs_print_anal(ao_bas_han, ilo, ihi, 0.15d0, 
     &           g_movecs(ispin), 
     &           blob, 
     &           .true., dbl_mb(k_eval(ispin)), oadapt, 
     &           int_mb(k_ir+(ispin-1)*nbf_ao),
     &           .true., dbl_mb(k_occ+(ispin-1)*nbf_ao))
         endif
      enddo
c
c Jorge Garza/Begin
c
      if (condfukui.eq.1) then
         call fukui(g_movecs, k_eval, tol2e, rtdb, nExc, iVxc_opt, 
     &              g_xcinv, IOLGC, g_wght, g_xyz, g_nq, wght_GA,
     &              rho_n, irdens_atom,
     &              icetobfr, natoms)
      endif
c
c Jorge Garza/End
c
      if (mulliken) call dft_mulwrap(me,g_dens,g_s)
c     
c     end infinite loop for SCF iterations
c     
c     Store energy and convergence status ... must store before
c     write movecs since date of insertion is used.
c     
      if (.not. rtdb_put(rtdb, 'dft:energy', MT_DBL, 1, (Etnew+Enuc)))
     &   call errquit('dft_scf: failed to store energy in rtdb', 0)
      if (.not. rtdb_put(rtdb, 'dft:converged', MT_LOG, 1, converged))
     &   call errquit('dft_scf: failed to store converged in rtdb', 0)
c
c     output energies and eigenvectors to disk
c     
      if (.not.movecs_write(rtdb, ao_bas_han, movecs_out, 'dft', title,
     &                      nbf_ao, ipol, nmo, dbl_mb(k_occ), nbf_ao, 
     &                      dbl_mb(k_eval(1)), nbf_ao, g_movecs))
     &                      call errquit('dft_scf: movec_write failed',
     &                      0)
c     
c     Shut down DIIS.
c     
      if (icall(1).gt.0)then
         icall(1) = -1
         call diis_driver(toll_s, derr(1),  ipol, icall, nfock, 
     &                 nbf_ao, geom, ao_bas_han, g_focks, g_dens(1), 
     &                    g_tmp, g_svecs, svals, diising, nodiis,
     .        int_mb(k_gdiis),lkeeps,g_s12m,g_s12p)
         if (.not. MA_free_heap(l_gdiis))
     .        call errquit(' dft_scf:cannot popstack',111)
      endif
c     
c     If open shell compute overlap of alpha orbitals with beta orbitals.
c     
      if (ipol.gt.1 .and. util_print('final overlap',print_default))then
         call dft_mxspin_ovlp(ao_bas_han, g_movecs(1), g_movecs(2), 
     &                        g_tmp)
      endif
c
      if (wght_GA)then
         if (.not. ga_destroy(g_wght)) call errquit
     &      ('dft_scf: could not destroy g_wght', 0)
         if (.not. ga_destroy(g_xyz)) call errquit
     &      ('dft_scf: could not destroy g_xyz', 0)
         if (.not. ga_destroy(g_nq)) call errquit
     &      ('dft_scf: could not destroy g_nq', 0)
      endif
c     
c     Restore alpha and beta densities.
c
      if (ipol .gt. 1)
     &   call ga_dadd(one,g_dens(1),onem,g_dens(2),g_dens(1))
c     
c     computation of <S2> for open shell
c     
      if (ipol.gt.1 .and. util_print('final s2', print_low))then

         call dft_s2_value(geom,AO_bas_han,.false.,noc(1),noc(2),
     &        nbf_ao,g_dens(1),g_dens(2))

      endif
c
c     write to rtdb that fractional occupancy exists
c
      if(fon.or.ssmear.gt.1d-9) then
      if (.not.rtdb_put(rtdb, 'dft:fractional_occup', mt_log, 1,
     &   .true.))call errquit('dftscf: cannot write rtdb',1)
      endif

c     
c     computation of moments
c
      if (natoms .gt. 1) then
         if(oprintinertia)call geom_momint(geom)
         call dft_mpole(rtdb, ao_bas_han, ipol, g_dens(1), g_dens(2))
      endif
c     
c     print stolen for uhf.F
c     
      if (util_print('schwarz',print_high).and.(.not.CDFIT))then
         call schwarz_print(natoms, nshells_ao)
      endif
c     
      call dft_prevals('final',me,oprint_eval,oprint_vecs,noc,
     ,     g_movecs,dbl_mb(k_eval(1)),dbl_mb(k_occ))
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
      if (nexc.gt. 1)then
         call ecce_print1 ('correlation energy', mt_dbl, Exc(2), 1)
      endif
      if (.not.ma_chop_stack(l_ir))
     &   call errquit('dft_scf: cannot chop stack',98)
c     
      if (ipol.gt.1)then
         if (.not. ga_destroy(g_fockt)) call errquit
     &      ('dft_scf: could not destroy g_fockt', 0)
      endif
      if (.not. ga_destroy(g_tmp)) call errquit
     &   ('dft_scf: could not destroy g_tmp', 0)
      if(lkeeps) then
      if (.not. ga_destroy(g_s12p)) call errquit
     &   ('dft_scf: could not destroy g_s12p', 0)
      if (.not. ga_destroy(g_s12m)) call errquit
     &   ('dft_scf: could not destroy g_s12m', 0)
      endif
     
c
      call fock_2e_tidy(rtdb)
c     
      if (converged)then
         call ecce_print_module_exit('dft', 'ok')
      else
         call ecce_print_module_exit('dft', 'failed')
      endif
c     
c     eval deallocation moved here from inside iteration loop
c     
      if (.not.ma_chop_stack(lcntoce))
     &   call errquit('dft_scf: cannot chop stack',99)
c
      dft_scf = converged
c
c !!! BGJ
      if (.not. rtdb_get(rtdb, 'bgj:poliz', mt_log,
     &     1, do_poliz)) then
         do_poliz = .false.
      endif
      if (do_poliz) then
         write(LuOut,*)'*** dft_scf: calling cphf_poliz'
         if (.not. cphf_poliz(rtdb)) ! Never executed.
     $        call errquit(' cphf_poliz: failed from dft_scf !',0)
      endif
c !!! BGJ

      return
c     
 21   format(/,10x,' Memory utilization after 1st SCF pass: ')
    1 format(/,
     &     1x,'  convergence    iter        energy       DeltaE   ',
     &     'RMS-Dens  Diis-err    time'/
     &     1x,'---------------- ----- ----------------- --------- ',
     &     '--------- ---------  ------')
    2 format(1x,'d=',i2,',ls=',f3.1,',diis',1x,i5,f18.10,
     &     1p,3d10.2,0p,f8.1)
    3 format(51x,1p,2d10.2)
 22   format(1x,'d=',i2,',ls=',f3.1,6x,i5,f18.10,
     &     1p,2d10.2,10x,0p,f8.1)
 23   format(51x,1p,1d10.2)
 1111 format(15x,'Core Energy:              ',f20.10)
 222  format(//
     &     '         Total DFT energy =', f20.12/
     &     '      One electron energy =', f20.12/
     &     '           Coulomb energy =', f20.12/
     &     '    Exchange-Corr. energy =', f20.12/
     &     ' Nuclear repulsion energy =', f20.12//
     &     ' Numeric. integr. density =', f20.12//
     &     '     Total iterative time =', f9.1,'s'//)
 223  format(//
     &     '         Total DFT energy =', f20.12/
     &     '      One electron energy =', f20.12/
     &     '           Coulomb energy =', f20.12/
     &     '          Exchange energy =', f20.12/
     &     '       Correlation energy =', f20.12/
     &     ' Nuclear repulsion energy =', f20.12//
     &     ' Numeric. integr. density =', f20.12//
     &     '     Total iterative time =', f9.1,'s'//)
 224  format(10x,' HOMO = ',f6.2,' LUMO = ',f6.2,
     &     /,10x,'Significant orbital reordering with',
     &     ' maximum overlap',
     &     /,10x,'turned ON.  Turning max_ovl OFF.')
 2224 format(10x,'HOMO = ',f6.2,' LUMO (with shift) = ',f6.2,
     &     /,10x,'Unshifted LUMO is less than HOMO.',
     &     /,10x,'Turning levelshifting OFF this iteration.')
 4224 format(10x,' HOMO = ',f6.2,' LUMO = ',f6.2,
     &              ' RLSHIFT = ',f6.2,' HL_GAP = ',f6.2)
c     
      end

      subroutine dft_dump_info(me)
      implicit none
      integer i, me
c      
#include "cdft.fh"
c
      write(*,*)' me = ',me
c
c     Logicals
c
      write(*,*)' LOGICALS '
      write(*,*)' nodamping = ', nodamping
      write(*,*)' nolevelshifting = ', nolevelshifting
      write(*,*)' nodiis = ', nodiis
      write(*,*)' ltight = ', ltight
      write(*,*)' lmaxov = ', lmaxov
      write(*,*)' incore = ', incore
      write(*,*)' ldelley = ', ldelley
      write(*,*)' store_wght = ', store_wght
      write(*,*)' XCFIT = ', XCFIT
      write(*,*)' CDFIT = ', CDFIT
      write(*,*)' levelshift = ', levelshift
      write(*,*)' damp = ', damp
      write(*,*)' diis = ', diis
      write(*,*)' direct = ', direct
      write(*,*)' oskel = ', oskel
      write(*,*)' oadapt = ', oadapt
      write(*,*)' lcfac(i) = ', (lcfac(i),i=1,15)
      write(*,*)' nlcfac(i) = ', (nlcfac(i),i=1,15)
      write(*,*)' lxfac(i) = ', (lxfac(i),i=1,15)
      write(*,*)' nlxfac(i) = ', (nlxfac(i),i=1,15)
      write(*,*)' xccomb(i) = ', (xccomb(i),i=1,15)
c
c     Integers
c
      write(*,*)' INTEGERS '
      write(*,*)' iAOacc = ', iAOacc
      write(*,*)' iCDacc = ', iXCacc
      write(*,*)' iXCacc = ', iXCacc
      write(*,*)' irqmax = ', irqmax
      write(*,*)' itol2e = ', itol2e
      write(*,*)' imull = ', imull
      write(*,*)' iangquad = ', iangquad
      write(*,*)' iterations = ', iterations
      write(*,*)' iguess = ', iguess
      write(*,*)' mult = ', mult
      write(*,*)' ipol = ', ipol
      write(*,*)' nfock = ', nfock
      write(*,*)' ncydp = ', ncydp
      write(*,*)' ncysh = ', ncysh
      write(*,*)' ncyds = ', ncyds
      write(*,*)' idirect = ', idirect
      write(*,*)' noio = ', noio
      write(*,*)' ndamp = ', ndamp
      write(*,*)' idecomp = ', idecomp
      write(*,*)' nagrid = ', nagrid
      write(*,*)' nquad_task = ', nquad_task
      write(*,*)' nrinc = ', nrinc
      write(*,*)' nrmax = ', nrmax
      write(*,*)' geom = ', geom
      write(*,*)' ncenters = ', ncenters
      write(*,*)' nbf = ', nbf
      write(*,*)' nradpts = ', nradpts
      write(*,*)' nang_leb_pts = ', nang_leb_pts
      write(*,*)' nbf_ao = ', nbf_ao
      write(*,*)' nbf_cd = ', nbf_cd
      write(*,*)' nbf_xc = ', nbf_xc
      write(*,*)' MXqshells = ', MXqshells
      write(*,*)' MXrad = ', MXrad
      write(*,*)' MXagrid = ', MXagrid
      write(*,*)' MXtheta = ', MXtheta
      write(*,*)' MXphi = ', MXphi
      write(*,*)' IAC_MAX = ', IAC_MAX
      write(*,*)' MAX_CHUNK = ', MAX_CHUNK
      write(*,*)' nqshells = ', nqshells
      write(*,*)' iqrot = ', iqrot
      write(*,*)' ntypes = ', ntypes
      write(*,*)' nquad_task = ', nquad_task
      write(*,*)' dft_ntags_bsmx = ', dft_ntags_bsmx
      write(*,*)' dft_max_cent = ', dft_max_cent
      write(*,*)' AO_bas_han = ', AO_bas_han
      write(*,*)' CD_bas_han = ', CD_bas_han
      write(*,*)' XC_bas_han = ', XC_bas_han
      write(*,*)' nrad: ',(nrad(i), i = 1, dft_ntags_bsmx)
      write(*,*)' nang: ',(nang(i), i = 1, dft_ntags_bsmx)
c
c     Doubles
c
      write(*,*)' DOUBLES '
      write(*,*)' rlshift = ', rlshift
      write(*,*)' dampon = ', dampon
      write(*,*)' dampoff = ', dampoff
      write(*,*)' levlon = ', levlon
      write(*,*)' levloff = ', levloff
      write(*,*)' diison = ', diison
      write(*,*)' diisoff = ', diisoff
      write(*,*)' rcharge = ', rcharge
      write(*,*)' rm = ', rm
      write(*,*)' tol_rho = ', tol_rho
      write(*,*)' radius = ', radius
      write(*,*)' rscale = ', rscale
      write(*,*)' cfac = ', (cfac(i),i=1,15)
      write(*,*)' xfac =', (xfac(i),i=1,15)
      write(*,*)' hl_tol = ', hl_tol
      write(*,*)' e_conv = ', e_conv
      write(*,*)' d_conv = ', d_conv
      write(*,*)' g_conv = ', g_conv
c
c     Characters
c
      write(*,*)' CHARACTERS '
      write(*,*)' title = ', title
      write(*,*)' gridtype = ', gridtype
      write(*,*)' xname = ', xname
      write(*,*)' cname = ', cname
c
      return
      end
