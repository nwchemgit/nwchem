      logical function dft_scf(rtdb, Etold, Enuc, iVcoul_opt, iVxc_opt, 
     &                         iter, g_dens, g_movecs, g_vxc, g_fock, 
     &                         g_svecs, isvals, g_cdinv, g_xcinv, 
     &                         g_2ceri)
c     
c     $Id: dft_scf.F,v 1.103 1998-08-20 23:52:22 d3e129 Exp $
c     
      implicit none
c
      integer rtdb              ! [input]
      double precision Etold, Enuc
      integer iVcoul_opt
      integer iVxc_opt
      integer iter
      integer g_dens(2), g_movecs(2), g_vxc(4), 
     &        g_fock, g_svecs, isvals, 
     &        g_cdinv, g_xcinv, g_2ceri, g_scr
      double precision rho_n
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
c     
#if defined(DFT_TIMER)
#include "timer.fh"
#endif
      Logical movecs_write, movecs_converged
      External movecs_write, movecs_converged
c
      integer  ga_create_atom_blocked
      external ga_create_atom_blocked
c
      logical oprint, oprint_sym
      double precision Exc(2), rms(2), derr(2), n3c_dbl_r, availm_r
      integer noc(2), nmo(2), icall(2)
      integer n3c_dbl, n3c_int, n_batch
      integer avail, availm
      integer iwhat_max
      integer l_3cwhat, k_3cwhat, l_3cERI, k_3cERI
      integer dft_n3cint, n_semi_bufs, fd
      external dft_n3cint
      double precision dft_n3cdbl
      external dft_n3cdbl
      Integer l_eval
      integer k_eval(2)
      integer lmax              ! CMM
      integer finest            ! CMM
      integer natoms, nTotEl
      integer l_occ, k_occ
      integer i, j, iw
      integer me, nproc
      integer g_tmp, g_fockt, g_s, g_wght, g_xyz
      integer nheap, nstack
      integer ispin, idone
      integer nexc
      integer iswitc
      integer itol_max, iaoacc_max
      integer itol_min, iAOacc_min
      double precision tol_rho_min, tol_rho_max
      integer nao2_max, batch_max, icount
      integer npol, nprint
      integer leneval, lcd_coef, icd_coef
      integer lcntoce, icntoce, lcntobfr, icntobfr,
     &        lcetobfr, icetobfr, lrdens_atom, irdens_atom,
     &        nscr, lscr, iscr, max_component, max_elem_ang_scr
      double precision start_wall, current_wall, elapsed_wall,
     &                 save_wall, current_cpu, start_cpu,
     &                 wall_time_reqd
      integer int_wall_time_reqd
      double precision ecoul, ecore
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
      Logical LResult, IOLGC, mulliken, lSPM
      logical converged, wght_GA
      logical oprint_parm, oprint_conv, oprint_vecs, 
     &        oprint_eval, oprint_syma, oprint_time, 
     &        oprint_info, oprint_tol, oprint_final_vecs, 
     &        oprint_energy_step, oprint_intermediate_fock,
     &        oprint_3c2e, oprint_interm_overlap, oprint_interm_S2
      double precision zero, onem, one
      parameter(zero = 0.d0, one = 1.d0, onem = -one)
c
      integer ilo, ihi          ! For printing movecs analysis
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
      integer stat1, stat2, l_ir, k_ir
      logical last_time_energy
      logical check_shift
      character*7 name
      character*4 scftype
      character*255 basis_name, basis_trans
      integer nopen, nclosed
c
      call ecce_print_module_entry('dft')
c
      dft_scf = .false.
      oprint = util_print('information', print_low)
      oprint_info = util_print('common', print_debug)
      oprint_parm = util_print('parameters', print_default)
      oprint_3c2e = util_print('3c 2e integrals', print_default)
      oprint_conv = util_print('convergence', print_default)
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
      me = ga_nodeid()
      nproc = ga_nnodes()
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
c
c     Set aside some memory for reduced density matrix
c
      if (.not.MA_Push_Get(MT_Dbl,ipol*natoms*natoms,'rdens_atom',
     &   lrdens_atom,irdens_atom))
     &   call errquit('dft_scf: cannot allocate rdens_atom',0)
c     
c     determine pattern of orbitals' occupancy
c     
      if (.not. rtdb_get(rtdb, 'dft:noc', mt_int, 2, noc))
     &   call errquit('dft_scf: rtdb_get of noc failed', 0)
      if (.not. MA_Push_Get(MT_Dbl,nbf_ao*ipol,'mo occ',l_occ,k_occ))
     &   call errquit('dft_scf: failed to alloc',999)
c      
      if (ipol.eq.1)then
c     
c        noc(1) = No. of electrons alpha
c     
         nTotEl = 2*noc(1)
         nmo(1) = nbf_ao
      else
c     
         nTotEl = noc(1) + noc(2)
         nmo(1) = nbf_ao
         nmo(2) = nbf_ao
      endif
      anoca = noc(1)
      anocb = noc(2)
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
#if defined(DFT_TIMER)
      call dft_second(tp3)
#endif
c
      if (CDFIT)then
c     
c        Determine how big a buffer can be allocated to 3-center
c        2e- integrals.
c     
c        - amount needed for all incore:
c     
         n3c_dbl_r = dft_n3cdbl(AO_bas_han,CD_bas_han,itol_max)
         n3c_int = dft_n3cint(AO_bas_han,CD_bas_han,itol_max)
c     
c        find - (minimum)amount local available memory on all nodes
c     
         call ga_sync
         avail = MA_inquire_avail(mt_dbl)
         call ga_igop(msg_min_stack_avail, avail, 1, 'min')
c     
c        estimate and subtract off amount needed for DIIS
c
         availm = avail - ((nfock+4)*nbf_ao*nbf_ao)/nproc
c
c        estimate and subtract off amount needed for XC numerical integration
c        in xc_quadv0
c
         availm = availm - (natoms*(natoms+1)/2 + 13*natoms + 
     &                      nqmax*(7*ipol +  npol + natoms + 
     &                             nbf_ao_mxnbf_ce + 4) + 
     &                      nbf_ao_mxcont + nbf_ao_mxprim + 
     &                      2*nbf_ao_mxnbf_ce*nbf_ao_mxnbf_ce)
c
c        estimate and subtract off amount needed for XC numerical integration
c        in xc_quadv0_a
c     
         max_component = (nbf_ao_mxang+1)*(nbf_ao_mxang+2)/2 
         if (nbf_ao_mxang .eq. 0)then
            max_elem_ang_scr = max_component * 3
         elseif (nbf_ao_mxang .le. 3)then
            max_elem_ang_scr = max_component * 9
         else ! general case
            max_elem_ang_scr = max_component * 28
         endif 
c
         nscr = 3*nqmax*nbf_ao_mxcont + 
     &          max(3*nqmax*nbf_ao_mxprim,nqmax*max_elem_ang_scr) + 1
c
c        The big chunk is the memory needed for new_eval_gbsets 
c        which is roughly 4*nqmax*nbf_ao.  This is reduced by
c        screening (and chunking up the angular grid) and is 
c        computed at the end of xc_setquad to be 4*max_pr_mbfnq.
c
         availm = availm - (nqmax*(natoms + 3*nbf_ao_mxnbf_ce + 1) +
     &                      4*max_pr_mbfnq + 
     &                      nbf_ao + nscr)
c
c        Subtract off a few extra bits
c
         availm = availm - 100000
c
         if(availm.lt.0)then
           availm = 0
         endif
         availm_r = dble(availm)
         availm_r = min(availm_r,n3c_dbl_r)
         if (availm_r.lt.n3c_dbl_r)then
c
c           cannot all be done incore so let's restrict memory use to a
c           max of 3 MW for semi-direct (I/O)
c
            availm_r = min(availm_r,3d+6)
         endif
         n3c_dbl = nint(availm_r)
         if (me.eq.0 .and. oprint_parm)
     &      write(LuOut,3228)avail, availm, n3c_dbl
c     
c        Loops are parallelized over the products of AO shells
c        (check for zero ... must be at least 1).
c     
         icount = (nshells_ao + mod(nshells_ao,2))/(2*nproc)
         icount = max(icount,1)
         if (nproc.eq.1)icount = 1
c     
         nao2_max = nbf_ao_mxnbf_cn**2
c
c        compute size of largest block of integrals computed
c
         batch_max = nao2_max*nbf_cd*icount
c     
         if ((n3c_dbl.lt.batch_max).or.direct)then
            incore = .false.
            n3c_dbl = 1
            n3c_int = 1
         endif
         if (availm_r.lt.n3c_dbl_r)then
c
c           cannot all be done incore, compute number of buffers 
c           required if semi-direct (I/O)
c
            n_semi_bufs = nint(n3c_dbl_r/availm_r) + 1
         else
            n_semi_bufs = 0
         endif
         if (.not.MA_Push_Get
     &      (MT_Dbl,n3c_dbl,'3c ERI space',l_3cERI,k_3cERI))
     &      call errquit('dft_scf: push_get failed', 12)
c     
         if (.not.MA_Push_Get
     &      (MT_int,n3c_int,'3c what space',l_3cwhat,k_3cwhat))
     &      call errquit('dft_scf:push_get failed', 13)
c     
         if (incore)then
            if (me.eq.0 .and. oprint_3c2e)
     &         write(LuOut,3229)n3c_dbl*1.d-6
            call dft_3cincor(rtdb, n_batch, n3c_int, int_mb(k_3cwhat), 
     &                       dbl_mb(k_3cERI), n3c_dbl, iwhat_max, 
     &                       n_semi_bufs, fd, int_mb(icntoce), 
     &                       int_mb(icntobfr), int_mb(icetobfr), natoms)
         else
            if (me.eq.0 .and. oprint_3c2e)write(LuOut,3230)
         endif
      endif
 3228 format(10x,'Minimum dble words available (all nodes) is: ',i15,
     &     /,10x,'         This is reduced (for later use) to: ',i15,
     &     /,10x,'                   Suggested buffer size is: ',i15)
 3229 format(/,10x,f6.3,' MW buffer allocated for incore 3-center '
     &     /,10x,'2e- integral storage on stack. ')
 3230 format(/,10x,'Incore memory use for 3-center 2e- integrals is ',
     &     'turned off. ')
#if defined(DFT_TIMER)
      call dft_second(tp4)
      pgmsect(5) = ' 3c-2e Ints'
      call tcllct(t_pgm(1,5),tp4,tp3)
#endif
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
      if (ipol.gt.1)then
         g_fockt = ga_create_atom_blocked(geom, AO_bas_han, 'fock tr')
      endif
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
      if (nodiis)diising = .false.
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
c     Do not allow levelshifting initially until sure that the
c     transformed Fock matrix will be diagonally dominant.
c     
      levelshifting = .false.
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
c     
c     Get initial density.
c     
      if (oprint_time)
     &     call dft_tstamp(' Before call to DFT_INIT.')
c      call dft_init_guess(rtdb, g_dens, g_movecs, g_fock)
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
      call scf_vectors_guess(rtdb, tol2e_sleazy, geom, ao_bas_han, 
     &                       basis_trans, movecs_in, movecs_out, 
     &                       movecs_guess, scftype, nclosed, nopen, 
     &                       nbf, nbf, noc(1), noc(2),  k_eval, k_occ, 
     &                       k_ir, g_movecs, g_dens, 'density', 
     &                       'dft', title, oskel, oadapt, 
     &                       .true.) 
c
c     Tidy up SCF
c
      call fock_2e_tidy(rtdb)
c
c     set initial coulomb acc
c
c      write(6,*)' movecs_guess = ',movecs_guess
      if (movecs_guess.eq.'restart')ltight=.true.
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
c     If open shell put the total density matrix in g_dens(1)
c
      if(ipol.eq.2)then             
         call ga_dadd(one,g_dens(1),one,g_dens(2),g_dens(1))
      endif
c     
c     Call to Mulliken Pop Ananlysis for initial density
c     
      if (mulliken)then
         g_s = ga_create_atom_blocked(geom, AO_bas_han, 'AO ovl')
         call ga_zero(g_s)
         call int_1e_ga(ao_bas_han,ao_bas_han,g_s,'overlap',.false.)
         if (me.eq.0)call dft_header
     &      (' Total Density - Mulliken Population Analysis')
         call mull_pop(geom,ao_bas_han,g_dens(1),g_s,'total')
         if (ipol.eq.2)then
c         
c           Analysis of spin density
c         
            if (me.eq.0) call dft_header
     &         (' Spin Density - Mulliken Population Analysis')
            call ga_dadd(one,g_dens(1),-2.d0,g_dens(2),g_dens(2))
            call mull_pop(geom,ao_bas_han,g_dens(2),g_s,'spin') 
c
c           Restore beta density in g_dens(2)
c
            call ga_dadd(one,g_dens(1),-1.d0,g_dens(2),g_dens(2))
            call ga_dscal(g_dens(2),0.5d0)
         endif
c         call gadest(g_s)
         if (.not. ga_destroy(g_s)) call errquit
     &      ('dft_scf: could not destroy g_s', 0)
      endif
      iter = 1
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
      last_time_energy = .false.
c
 1000 continue
c
      if (me.eq.0.and.oprint_tol)write(LuOut,3234)itol2e,iAOacc,iXCacc
 3234 format(10x,'itol2e=',i2,' iAOacc=',i2,' iXCacc=',i2)

      Ecoul  = ZERO
      Exc(1)= ZERO
      Exc(2) = ZERO
c     
c     Accumulate core hamiltonian into Fock matrix; compute core energy
c     
      call ga_zero(g_fock)
      call int_1e_ga(ao_bas_han, ao_bas_han, g_fock, 'kinetic', oskel)
      call int_1e_ga(ao_bas_han, ao_bas_han, g_fock, 'potential', oskel)
c
      Ecore = ga_ddot(g_dens(1), g_fock)
c     
c     Perform SPM analysis
c     
      LSPM=.true.
      if (.not. rtdb_get(rtdb, 'dft:cmm', mt_int, 1, lmax))
     &                   LSPM = .false.
      if (LSPM)then
         lresult = geom_charge_center(geom)
         if (.not. rtdb_get(rtdb, 'dft:finest', mt_int, 1, finest))
     &        finest = 4
         if (oprint_time)call dft_tstamp(' Before call to CMM.   ')
         call dft_getspm(geom,lmax,finest,g_dens,ao_bas_han)
         call errquit('dft_scf: end of CMM ',1234)
      endif
c
c     Pre-compute reduced total density matrices over atoms
c 
      call dfill(ipol*natoms*natoms, 0.0d0, dbl_mb(irdens_atom), 1)
      nscr = nbf_ao_mxnbf_ce*nbf_ao_mxnbf_ce
      if (.not.MA_Push_Get(MT_Dbl,nscr,'scr',lscr,iscr))
     &   call errquit('dft_scf: cannot allocate scr',0)
      call util_ga_mat_reduce(nbf_ao, natoms, int_mb(icetobfr), g_dens, 
     &                        ipol, dbl_mb(irdens_atom), 'rms', 
     &                        dbl_mb(iscr), nbf_ao_mxnbf_ce)
      if (.not.ma_pop_stack(lscr))
     &   call errquit('dft_scf: cannot pop stack',0)
c
#if defined(DFT_TIMER)
      call dft_second(tp3)
#endif
c
      if (CDFIT)then
c     
c        Fit the electron charge density.
c     
         if (.not.MA_Push_Get(MT_Dbl,nbf_cd,'cd_coef',lcd_coef,
     &        icd_coef))
     &        call errquit('dft_scf: cannot allocate cd_coef',0)
         if (oprint_time)
     &        call dft_tstamp(' Before call to FITCD.   ')
         call dft_fitcd(Dbl_MB(icd_coef), dbl_mb(k_3cERI), Ecoul, rtdb,
     &                  g_cdinv, g_dens, nTotEl, n_batch, n3c_int,
     &                  int_mb(k_3cwhat), n3c_dbl, iwhat_max, 
     &                  n_semi_bufs, fd, IOLGC, 
     &                  dbl_mb(irdens_atom), int_mb(icntoce), 
     &                  int_mb(icntobfr), int_mb(icetobfr), natoms)
      endif
c     
#if defined(DFT_TIMER)
      call dft_second(tp4)
      pgmsect(6) = '   FitCD  '
      call tcllct(t_pgm(1,6),tp4,tp3)
      call dft_second(tp3)
#endif
      if (oprint_time)
     &     call dft_tstamp(' Before call to GETVCOUL.')
      call dft_getvc(Dbl_MB(icd_coef), dbl_mb(k_3cERI), Ecoul, rtdb,
     &               g_tmp, iVcoul_opt, g_2ceri, g_dens, n_batch, 
     &               n3c_int, int_mb(k_3cwhat), n3c_dbl, iwhat_max,
     &               n_semi_bufs, fd, IOLGC,
     &               dbl_mb(irdens_atom), int_mb(icntoce), 
     &               int_mb(icntobfr), int_mb(icetobfr), natoms)
c     
c     Add V coul to Fock Matrix
c     
      call ga_dadd(one, g_tmp, one, g_fock, g_fock)
      if (CDFIT)then
         if (.not.ma_pop_stack(lcd_coef))
     &        call errquit('dft_scf: cannot pop stack',0)
      endif
c     
#if defined(DFT_TIMER)
      call dft_second(tp4)
      pgmsect(7) = '   GtVcoul'
      call tcllct(t_pgm(1,7),tp4,tp3)
      call dft_second(tp3)
#endif
      if (ipol.gt.1)then
c     
c        Restore alpha and beta densities.
c     
         call ga_dadd(one, g_dens(1), onem, g_dens(2), g_dens(1))
c     
c        Note that g_dens(1) now contains the alpha density
c        matrix and g_dens(2) contains the beta
c
c        Pre-compute reduced alpha and beta density matrices over atoms
c 
         call dfill(ipol*natoms*natoms, 0.0d0, dbl_mb(irdens_atom), 1)
         nscr = nbf_ao_mxnbf_ce*nbf_ao_mxnbf_ce
         if (.not.MA_Push_Get(MT_Dbl,nscr,'scr',lscr,iscr))
     &      call errquit('dft_scf: cannot allocate scr',0)
         call util_ga_mat_reduce(nbf_ao, natoms, int_mb(icetobfr), 
     &                           g_dens, ipol, dbl_mb(irdens_atom), 
     &                           'rms', dbl_mb(iscr), nbf_ao_mxnbf_ce)
         if (.not.ma_pop_stack(lscr))
     &      call errquit('dft_scf: cannot pop stack',0)
      endif
c     
c     Compute the XC potential and energy.
c     
      g_vxc(1) = g_tmp
      call ga_zero(g_vxc(1))
      rho_n = 0.0d0
      if (ipol.eq.2)call ga_zero(g_vxc(2))
      if (oprint_time)call dft_tstamp(' Before call to GETVXC.  ')
      call xc_getv(rtdb, Exc, nExc, iVxc_opt, g_xcinv, g_dens, 
     &             g_vxc, IOLGC, g_wght, g_xyz, wght_GA, rho_n,
     &             dbl_mb(irdens_atom), int_mb(icntoce), 
     &             int_mb(icntobfr), int_mb(icetobfr), natoms)
      call ga_sync
      call ga_dadd(1.d0, g_vxc(1), 1.d0, g_fock, g_fock)
      if (ipol.gt.1)then
c     
c     Put (-alpha+beta) in g_vxc(2)
c     
         call ga_dadd(-1.d0, g_vxc(1), 1.d0, g_vxc(2), g_vxc(2))
      endif
c     
      if (oprint_time)
     &     call dft_tstamp(' End of parallel region. ')
c     
#if defined(DFT_TIMER)
      call dft_second(tp4)
      pgmsect(8) = '   GetVXC '
      call tcllct(t_pgm(1,8),tp4,tp3)
#endif
c     
c     Calculate the total electronic energy.
c     
      if (nExc.eq.1)then
         Etnew = Ecore + Ecoul + Exc(1)
      else
         Etnew = Ecore + Ecoul + Exc(1) + Exc(2)
      endif
c      write(6,*)'Etnew, Ecore, Ecoul, Exc(1), Exc(2): ',
c     &           Etnew, Ecore, Ecoul, Exc(1), Exc(2)
c     
      if (last_time_energy)then
c     
c        If open shell put the total density matrix back in 
c        g_dens(1) and quit.
c     
         if (ipol.eq.2)then
            call ga_dadd(one, g_dens(1), one, g_dens(2), g_dens(1))
         endif
         goto 2000
      endif
c
      delta = abs(etold-etnew)
c
      call ga_sync
      rms(1) = 0.d0
      rms(2) = 0.d0
      homo_lumo_gap = 200.0d0
      do ispin = 1, ipol
c     
c        Add V xc to Fock Matrix
c     
         if (ispin.eq.2) 
     &        call ga_copy(g_fockt, g_fock)
         if (ispin.eq.2)then
            call ga_dadd(one, g_vxc(ispin), one, g_fock, g_fock)
         endif
c     
#if defined(DFT_TIMER)
         call dft_second(tp3)
#endif
c     
         if (ispin.eq.1.and.ipol.gt.1)call ga_copy(g_fock, g_fockt)
c
c        Symmetrize the Fock matrix
c
         if (oskel)
     &      call sym_symmetrize(geom, AO_bas_han, .false., g_fock)
c
         if (diising)then
c     
c           DIIS step taken here.
c     
            if (diising)then
               call diis_driver(derr(ispin), ispin, ipol, icall, nfock, 
     &                          nbf_ao, geom, ao_bas_han, g_fock, 
     &                          g_dens(ispin), g_tmp, g_svecs, isvals, 
     &                          diising, nodiis)
            endif
         endif
c     
#if defined(DFT_TIMER)
         call dft_second(tp4)
         pgmsect(9) = '   DIIS   '
         call tcllct(t_pgm(1,9),tp4,tp3)
#endif
c     
#if defined(DFT_TIMER)
         call dft_second(tp3)
#endif
c
         g_scr = ga_create_atom_blocked(geom, AO_bas_han, 'ga scr')
c     
c        Put s-1/2 in g_scr.
c     
         iw = 2
         call diis_bld12(dbl_mb(isvals), g_svecs, g_scr, 
     &                   g_tmp, nbf_ao, iw)
c     
c        Transform Fock matrix.
c     
         call ga_dgemm('T', 'N', nbf_ao, nbf_ao, nbf_ao, one, 
     &                 g_scr, g_fock, zero, g_tmp)
         call ga_dgemm('N', 'N', nbf_ao, nbf_ao, nbf_ao, one, 
     &                 g_tmp, g_scr, zero, g_fock)
c     
c        Level shifting is implemented here (similarity 
c        transformation before standard eigensolver).  Note,
c        levelshifting is appropriate once a transformation
c        is available which makes the resulting Fock matrix 
c        diagonally dominant, e.g., in an approximate MO basis.  
c        Also note, there are many matrix multiplies with S^+-1/2 
c        which are redundant if one is sure that the former basis
c        is orthonormal.
c     
         if (levelshifting)then
c     
c           Build a matrix which is diagonal in the "MO" rep,
c           back-transform, and shift the current Fock matrix
c     
c           Use S^+1/2 * old movecs (as a transform).
c     
            iw = 3
            call diis_bld12(dbl_mb(isvals), g_svecs, g_scr, 
     &                      g_tmp, nbf_ao, iw)
            call ga_dgemm('N', 'N', nbf_ao, nbf_ao, nbf_ao, one, 
     &                    g_scr, g_movecs(ispin), zero, g_tmp)
            call ga_copy(g_tmp,  g_movecs(ispin))
c     
c           Build diagonal matrix.
c     
            call ga_zero(g_tmp)
            do j = noc(ispin)+1+me, nbf_ao, nproc
               call ga_put(g_tmp, j, j, j, j, rlshift, 1)
            enddo
c     
c           Transform this into "AO" basis and add to current 
c           Fock matrix
c     
            call ga_dgemm('N', 'N', nbf_ao, nbf_ao, nbf_ao, one, 
     &                    g_movecs(ispin), g_tmp, zero, g_scr)
            call ga_dgemm('N', 'T', nbf_ao, nbf_ao, nbf_ao, one, 
     &                    g_scr, g_movecs(ispin), one, g_fock)
         else
            rlshift = 0.0
         endif
c     
c        Solve for the eigenvalues and eigenvectors of the Hamiltonian.
c     
         if (oprint_intermediate_fock)then     
            call ga_print(g_fock)
         endif
#if defined(PARALLEL_DIAG)
         call ga_diag_std(g_fock, g_tmp, Dbl_MB(k_eval(ispin)))
#else 
         call ga_diag_std_seq(g_fock, g_tmp, Dbl_MB(k_eval(ispin)))
#endif
c     
c        Check HOMO/LUMO gap.
c     
         if(noc(ispin).gt.0) then
           homo = Dbl_MB(k_eval(ispin)+noc(ispin)-1)
           lumo = Dbl_MB(k_eval(ispin)+noc(ispin))
         else
           homo=-100.0d0
           lumo=100.0d0
         endif
c     
c        If levelshifting then tidy up.
c     
         if (levelshifting)then
c     
c           Put S^-1/2 back in g_scr (use g_fock as temp scr).
c     
            iw = 2
            call diis_bld12(dbl_mb(isvals), g_svecs, g_scr, 
     &                      g_fock, nbf_ao, iw)
c     
c           Put old-movecs back.
c     
            call ga_dgemm('N', 'N', nbf_ao, nbf_ao, nbf_ao, one, 
     &                    g_scr, g_movecs(ispin), zero, g_fock)
            call ga_copy(g_fock, g_movecs(ispin))
         endif
c     
c        Back-transform eigenvectors with S^-1/2.
c     
         call ga_dgemm('N', 'N', nbf_ao, nbf_ao, nbf_ao, one, 
     &                 g_scr, g_tmp, zero, g_fock)
         if (.not. ga_destroy(g_scr)) call errquit
     &      ('dft_scf: could not destroy g_scr', 0)
c     
c        Keep orbital ordering according to principle
c        of maximum overlap with previous iteration.
c     
         if (lmaxov)then
            call dft_mxovl(ao_bas_han, dbl_mb(k_eval(ispin)),
     &                     g_fock, g_movecs(ispin), g_tmp)
            homo = -100.0d0
            lumo =  100.0d0
            if(noc(ispin).gt.0) then
              do i = 1, noc(ispin)
                 homo = max(homo,(Dbl_MB(k_eval(ispin)+i-1)))
              enddo
              do i = noc(ispin)+1, nbf
                 lumo = min(lumo,(Dbl_MB(k_eval(ispin)+i-1)))
              enddo
            endif
         endif
         call ga_copy(g_fock, g_movecs(ispin))
c        
c        determine homo-lumo gap 
c
         homo_lumo_gap = min(homo_lumo_gap, (lumo-homo-rlshift))
c     
         call ga_sync
#if defined(DFT_TIMER)
         call dft_second(tp4)
         pgmsect(10) = '   Diag   '
         call tcllct(t_pgm(1,10),tp4,tp3)
#endif
c     
c        Save previous density for convergence check.
c     
         call ga_copy(g_dens(ispin), g_tmp)
c
c        symmetry adapt vectors?
c
         if (oadapt)then
            if(ispin.eq.1)name = '- alpha'
            if(ispin.eq.2)name = '- beta'
            call scf_movecs_sym_adapt(ao_bas_han, g_movecs(ispin),
     &                                oprint_syma, nbf_ao, name,
     &                                .true., 
     &                                int_mb(k_ir+nbf_ao*(ispin-1)))
c23456789012345678901234567890123456789012345678901234567890123456789012

         endif      
c     
c        Form a new density matrix.
c     
         call dft_densm(g_dens(ispin), g_movecs(ispin), 
     &                  nbf_ao, ipol, noc(ispin))
c     
         call ga_sync
c     
c        Check convergence on Density.
c     
         rms(ispin) = dft_dencvg(g_dens(ispin), g_tmp, nbf_ao)
         call ga_sync
c     
c        Damping implemented here.
c     
         if (damping)then
            pp = dble(ndamp)*1.d-2
            onempp = 1.0d0 - pp
            call ga_dadd(pp, g_tmp, onempp, g_dens(ispin), 
     &                   g_dens(ispin))
         else
            ndamp = 0
         endif
         call ga_sync
      enddo                     ! end loop over ispin
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
c     
      call ga_sync
c     
      if (me .eq. 0.and.oprint_eval)then
         if (util_print('intermediate evals', print_default))then
            call util_print_centered(LuOut,'alpha eigenvalues',
     &                               20,.true.)
            call output(dbl_mb(k_eval(1)), 1, min(noc(1)+10,nbf_ao),
     &                  1, 1, nbf_ao, 1, 1)
            if (ipol.eq.2)then
               write(LuOut,*)
               call util_print_centered(LuOut,'beta eigenvalues',
     &                                  20,.true.)
               call output(dbl_mb(k_eval(2)), 1, min(noc(1)+10,nbf_ao),
     &                     1, 1, nbf_ao, 1, 1)
            endif
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
         call ga_print(g_movecs)
         if (ipol.eq.2)call ga_print(g_movecs(2))
      endif
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
      call ga_sync
c     
c     Check for SCF convergence.
c     
      call ga_sync
      call dft_scfcvg(rms, derr, Etold, Etnew, Enuc, 
     &                e_conv, d_conv, g_conv, ipol, 
     &                iter, iterations, idone, rtdb,
     &                converged, diising)
      if (delta.lt.1.d-3)then
c     
c        Set coulomb acc to max (e.g., input parameter).
c        (note, may also require re-initializing DIIS)
c     
         itol2e = itol_max
         iAOacc = iAOacc_max
         tol_rho = tol_rho_max
         iswitc = iswitc+1
c     
c        Re-initialize DIIS
c     
c         if (icall(1).gt.0)then
c            icall(1) = -1
c            call diis_driver(derr(1), 1, ipol, icall, nfock, nbf_ao, 
c     &                       geom, ao_bas_han, g_fock, g_dens(1), 
c     &                       g_tmp, g_svecs, isvals, diising, nodiis)
c            icall(1) = 0
c            icall(2) = 0
c         endif
      endif
c     
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
     &    (iter.le.ncysh))then
         if (homo_lumo_gap.lt.hl_tol)then
            levelshifting = .true.
            rlshift = rlshift_input
c            if (check_shift)then
c               if (lumo .lt. homo)then
c                  levelshifting = .false.
c                  if (me.eq.0 .and. oprint_conv)
c     &               write(LuOut,2224)homo, lumo
c               endif
c            endif
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
c      if ((lumo - homo).lt.-hl_tol.and.lmaxov)then
c         lmaxov = .false.
c         if (me.eq.0 .and. oprint_conv)write(LuOut,224)homo, lumo
c 224     format(10x,' HOMO = ',f6.2,' LUMO = ',f6.2,
c     &        /,10x,'Significant orbital reordering with',
c     &        ' maximum overlap',
c     &        /,10x,'turned ON.  Turning max_ovl OFF.')
c      endif
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
      if (idone.eq.0.or.(iswitc.lt.2.and.iter.lt.iterations))
     &   go to 1000 ! begin new iteration
      if (idone.eq.1.and.(.not.last_time_energy))then
         last_time_energy = .true.
         go to 1000             ! build final total energies
      endif
c     
 2000 continue
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
 222     format(//
     &        '         Total DFT energy =', f20.12/
     &        '      One electron energy =', f20.12/
     &        '           Coulomb energy =', f20.12/
     &        '    Exchange-Corr. energy =', f20.12/
     &        ' Nuclear repulsion energy =', f20.12//
     &        ' Numeric. integr. density =', f20.12//
     &        '     Total iterative time =', f9.1,'s'//)
 223     format(//
     &        '         Total DFT energy =', f20.12/
     &        '      One electron energy =', f20.12/
     &        '           Coulomb energy =', f20.12/
     &        '          Exchange energy =', f20.12/
     &        '       Correlation energy =', f20.12/
     &        ' Nuclear repulsion energy =', f20.12//
     &        ' Numeric. integr. density =', f20.12//
     &        '     Total iterative time =', f9.1,'s'//)
         call util_flush(LuOut)
      endif
c
c     symmetry adapt vectors last time print symmetries, etc.
c
c      if (oadapt)then
c         call scf_movecs_sym_adapt(ao_bas_han, g_movecs, oprint, 
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
            call movecs_print_anal(ao_bas_han, ilo, ihi, 0.15d0, 
     &           g_movecs(ispin), 
     &           'DFT Final Molecular Orbital Analysis', 
     &           .true., dbl_mb(k_eval(ispin)), oadapt, 
     &           int_mb(k_ir+(ispin-1)*nbf_ao),
     &           .true., dbl_mb(k_occ+(ispin-1)*nbf_ao))
         endif
      enddo
c     
c     call to Mulliken Pop Ananlysis
c     
      if (mulliken)then
         g_s = ga_create_atom_blocked(geom, AO_bas_han, 'AO ovl')
         call ga_zero(g_s)
         call int_1e_ga(ao_bas_han,ao_bas_han,g_s,'overlap',.false.)
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
         if (.not. ga_destroy(g_s)) call errquit
     &      ('dft_scf: could not destroy g_s', 0)
      endif
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
      call ga_sync()
c     
c     Shut down DIIS.
c     
      if (icall(1).gt.0)then
         icall(1) = -1
         call diis_driver(derr(1), 1, ipol, icall, nfock, nbf_ao, geom, 
     &                    ao_bas_han, g_fock, g_dens(1), g_tmp, g_svecs, 
     &                    isvals, diising, nodiis)
      endif
c     
c     If open shell compute overlap of alpha orbitals with beta orbitals.
c     
      if (ipol.gt.1)then
         call dft_mxspin_ovlp(ao_bas_han, g_movecs(1), g_movecs(2), 
     &                        g_tmp)
      endif
c
      if (wght_GA)then
         if (.not. ga_destroy(g_wght)) call errquit
     &      ('dft_scf: could not destroy g_wght', 0)
         if (.not. ga_destroy(g_xyz)) call errquit
     &      ('dft_scf: could not destroy g_xyz', 0)
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
c     computation of moments
c
      if (natoms .gt. 1)
     &   call dft_mpole(rtdb, ao_bas_han, ipol, g_dens(1), g_dens(2))
c     
c     print stolen for uhf.F
c     
      if (util_print('schwarz',print_high).and.(.not.CDFIT))then
         call schwarz_print(natoms, nshells_ao)
      endif
c     
      if (me .eq. 0)then
         if (util_print('final evals', print_high))then
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
         endif
         if (oprint_final_vecs)then
            write(LuOut,*)
            call util_print_centered(
     &           LuOut,'Final MO vectors',40,.true.)
            write(LuOut,*)
            call util_flush(LuOut)
         endif
      endif
      if (oadapt)then
         nprint = min(noc(1)+10,nbf_ao)
         if (util_print('all vector symmetries', print_high))then
            nprint = nbf_ao
            oprint_sym = util_print('final vector symmetries',
     &           print_default)
            call scf_movecs_sym_adapt(ao_bas_han, g_movecs,
     &           oprint_sym, nprint, '- alpha', .false., 1)
            if (ipol.eq.2)
     &         call scf_movecs_sym_adapt(ao_bas_han, g_movecs(2), 
     &         oprint_sym, nprint, '- beta', .false., 1)
         endif
      endif
      if (oprint_final_vecs)then
         call ga_print(g_movecs)
         if (ipol.eq.2)call ga_print(g_movecs(2))
      endif
c     
c     ECCE printout
c     
c     requested?
      stat1 = 0
c
      call ecce_print_control(stat1, stat2) ! stat2 set to print status
      stat1 = stat2
      call ecce_print_control(stat1, stat2) ! restore previous setting
c
c     we need the value from process 0
c
      call ga_brdcst(12321, stat1, ma_sizeof(mt_int,1,mt_byte), 0)
      if (stat1 .eq. 1)then
c
         if (oadapt)then
            call scf_movecs_sym_adapt(ao_bas_han, g_movecs(1),
     &                                .false., nbf_ao, '- alpha',
     &                                .true., int_mb(k_ir))
         else
            call ifill (nbf_ao, 1, int_mb(k_ir), 1)
         endif
         call movecs_ecce(nbf_ao, nmo, 1, min(noc(1)+20,nmo(1)), 
     $                    dbl_mb(k_eval(1)),
     &                    dbl_mb(k_occ), int_mb(k_ir), 
     &                    g_movecs(1), 'dft', 'alpha')
         if (ipol.eq.2)then ! spin-unrestricted
            if (oadapt)then
               call scf_movecs_sym_adapt(ao_bas_han, g_movecs(2),
     &                                   .false., nbf_ao, '- alpha',
     &                                   .true., int_mb(k_ir))
            endif
            call movecs_ecce(nbf_ao, nmo, 1, min(noc(1)+20,nmo(1)),
     $                       dbl_mb(k_eval(2)),
     &                       dbl_mb(k_occ+nbf_ao), int_mb(k_ir), 
     &                       g_movecs(2), 'dft', 'alpha')
         endif
         call ecce_print1 ('total energy', mt_dbl, (Etold+Enuc), 1)
         call ecce_print1 ('nuclear repulsion energy', mt_dbl, Enuc, 1)
         call ecce_print1 ('coulomb energy', mt_dbl, Ecoul, 1)
         call ecce_print1 ('exchange energy', mt_dbl, Exc(1), 1)
         if (nexc.gt. 1)then
            call ecce_print1 ('correlation energy', mt_dbl, Exc(2), 1)
         endif
      endif
      if (.not.ma_pop_stack(l_ir))
     &   call errquit('dft_scf: cannot pop stack',0)
c     
      if (ipol.gt.1)then
         if (.not. ga_destroy(g_fockt)) call errquit
     &      ('dft_scf: could not destroy g_fockt', 0)
      endif
      if (.not. ga_destroy(g_tmp)) call errquit
     &   ('dft_scf: could not destroy g_tmp', 0)
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
      if (.not.ma_pop_stack(l_eval))
     &   call errquit('dft_scf: cannot pop stack',0)
      if (CDFIT)then
         if (.not.ma_pop_stack(l_3cwhat))
     &      call errquit('dft_scf: cannot pop stack',0)
         if (.not.ma_pop_stack(l_3cERI))
     &      call errquit('dft_scf: cannot pop stack',0)
      endif
      if (.not.ma_pop_stack(l_occ))
     &   call errquit('dft_scf: cannot pop stack',0)
      if (.not.ma_pop_stack(lrdens_atom))
     &   call errquit('dft_scf: cannot pop stack',0)
      if (.not.ma_pop_stack(lcetobfr))
     &   call errquit('dft_scf: cannot pop stack',0)
      if (.not.ma_pop_stack(lcntobfr))
     &   call errquit('dft_scf: cannot pop stack',0)
      if (.not.ma_pop_stack(lcntoce))
     &   call errquit('dft_scf: cannot pop stack',0)
c
      dft_scf = converged
c
      return
c     
 1111 format(15x,'Core Energy:              ',f20.10)
c     
      end

      subroutine build_maps(basis, cntoce, cntobfr, cetobfr, natoms, 
     &                      nshells)
      implicit none
      integer basis, natoms, nshells
      integer cntoce(nshells), cntobfr(2,nshells), cetobfr(2,natoms)
c      
#include "bas.fh"
c
      integer ish, iat
c     
c     Build maps (for speed).
c     
      do ish = 1, nshells
         if (.not. bas_cn2ce(basis, ish, cntoce(ish)))
     &        call errquit('build_maps: bad basis', 0)
      end do
c
      do ish = 1, nshells
         if (.not. bas_cn2bfr(basis,ish,cntobfr(1,ish),cntobfr(2,ish)))
     &        call errquit('build_maps: bad basis', 0)
      end do
c
      do iat = 1, natoms
         if (.not. bas_ce2bfr(basis,iat,cetobfr(1,iat),cetobfr(2,iat)))
     &        call errquit('build_maps: bad basis', 0)
      end do
c
      return
      end
      subroutine util_ga_mat_reduce(n, nr, map, g_a, n_a, r, op,
     &                              scr, mxmap)
      implicit none
c     
      integer n                         ! Original size [input]
      integer nr                        ! Reduced size  [input]
      integer map(2,nr)                 ! map(1,*)=lo, map(2,*)=hi [input]
      integer n_a                       ! number of GAs to reduce
      integer g_a(n_a)                  ! Original GA handle(s) [input]
      double precision r(nr,nr,n_a)     ! Reduced matrix [output]
      integer mxmap                     ! max size of map vector
      double precision scr(mxmap,mxmap) ! scratch
      character*(*) op                  ! Reduction operation
c     
c     R(i,j) <= R(i,j) op A(map(1,i):map(2,i),map(1,j):map(2,j))
c     
c     where op is one of 'abssum', 'absmax', 'rms' (extend as necessary)
c     
      integer ir, jr, i, j, k, ielems, jelems
      double precision sum
c     
      do k = 1, n_a
c
c         write(6,*) ' util_mat_reduce: input matrix '
c         call ga_print(g_a(k))
c     
         do jr = 1, nr
            do ir = 1, nr
               sum = 0.0d0
               ielems = map(2,ir) - map(1,ir) + 1
               jelems = map(2,jr) - map(1,jr) + 1
               if (ielems.gt.0 .and. jelems.gt.0) then
                  call ga_get(g_a(k), map(1,ir), map(2,ir), map(1,jr), 
     &                 map(2,jr), scr, mxmap)
                  if (op .eq. 'abssum') then
                     do j = 1, jelems
                        do i = 1, ielems
                           sum = sum + abs(scr(i,j))
                        end do
                     end do
                  else if (op .eq. 'absmax') then
                     do j = 1, jelems
                        do i = 1, ielems
                           sum = max(sum, abs(scr(i,j)))
                        end do
                     end do
                  else if (op .eq. 'rms') then
                     do j = 1, jelems
                        do i = 1, ielems
                           sum = sum + scr(i,j)*scr(i,j)
                        end do
                     enddo
                     sum = sqrt(sum)
                  else
                     call errquit('util_ga_mat_reduce: unknown op', 0)
                  end if
                  r(ir,jr,k) = max(r(ir,jr,k),sum)
               endif
            end do
         end do
c         write(6,*) ' util_mat_reduce: reduced matrix '
c         call output(r(1,1,k), 1, nr, 1, nr, nr, nr, 1)
      enddo
c     
      end
      subroutine util_irreg_mat_reduce(n_row, n_col, nr_row, nr_col, 
     &                                 row_map, col_map, a, r, op)
      implicit none
c     
      integer n_row                      ! Original row size [input]
      integer n_col                      ! Original col size [input]
      integer nr_row                     ! Reduced row size  [input]
      integer nr_col                     ! Reduced col size  [input]
      integer row_map(2,nr_row)          ! map(1,*)=lo, map(2,*)=hi [input]
      integer col_map(2,nr_col)          ! map(1,*)=lo, map(2,*)=hi [input]
      double precision a(n_row, n_col)   ! Original matrix [input]
      double precision r(nr_row, nr_col) ! Reduced matrix [output]
      character*(*) op                   ! Reduction operation
c     
c     R(i,j) <= R(i,j) op A(map(1,i):map(2,i),map(1,j):map(2,j))
c     
c     where op is one of 'abssum', 'absmax', 'rms' (extend as necessary)
c     
      integer ir, jr, i, j
      double precision sum
c     
      do jr = 1, nr_col
         do ir = 1, nr_row
            sum = 0.0d0
            if (op .eq. 'abssum') then
               do j = col_map(1,jr), col_map(2,jr)
                  do i = row_map(1,ir), row_map(2,ir)
                     if (i.ne.0.and.j.ne.0)sum = sum + abs(a(i,j))
                  end do
               end do
            else if (op .eq. 'absmax') then
               do j = col_map(1,jr), col_map(2,jr)
                  do i = row_map(1,ir), row_map(2,ir)
                     if (i.ne.0.and.j.ne.0)sum = max(sum, abs(a(i,j)))
                  end do
               end do
            else if (op .eq. 'rms') then
               do j = col_map(1,jr), col_map(2,jr)
                  do i = row_map(1,ir), row_map(2,ir)
                     if (i.ne.0.and.j.ne.0)sum = sum + a(i,j)*a(i,j)
                  end do
               enddo
               sum = sqrt(sum)
            else
               call errquit('util_irreg_mat_reduce: unknown op', 0)
            end if
            r(ir,jr) = max(r(ir,jr),sum)
         end do
      end do
c     
c      write(6,*) ' util_irreg_mat_reduce: input matrix '
c      call output(a, 1, n_row, 1, n_col, n_row, n_col, 1)
c
c      write(6,*) ' Row map begin: ',(row_map(1,ir),ir = 1,nr_row)
c      write(6,*) ' Row map end: ',(row_map(2,ir),ir = 1,nr_row)
c      write(6,*) ' Col map begin: ',(col_map(1,ir),ir = 1,nr_col)
c      write(6,*) ' Col map end: ',(col_map(2,ir),ir = 1,nr_col)
c     
c      write(6,*) ' util_irreg_mat_reduce: reduced matrix '
c      call output(r, 1, nr_row, 1, nr_col, nr_row, nr_col, 1)
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
