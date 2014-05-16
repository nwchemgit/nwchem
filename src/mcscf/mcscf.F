      logical function mcscf( rtdb )
*
* $Id$
*
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "util.fh"
#include "sym.fh"
#include "pstat.fh"
#include "mcscfprof.fh"
#include "inp.fh"
c     
      integer rtdb
c     
      integer geom, basis
      integer nbf, nclosed, nvir, nact, noper, nsym
      integer nactel, nela, nelb, multip, orlen, symstate
      character*8 state
      integer l_occ, k_occ, l_sym, k_sym
      integer l_evals, k_evals
      integer l_dm1, k_dm1, l_dm2, k_dm2
      integer l_tmp, k_tmp
      integer g_movecs, g_coul, g_exch
      integer g_conj, g_gradold, g_hcore
      integer mo_lo, mo_hi
      integer i, j, iter, maxiter, nmixed, clo, chi
      integer icgreset                ! Number of iterations between CG resets
      integer blen
c     
      double precision pfac
      double precision thresh, lshift0
      double precision eone, etwo, enrep, energy
      double precision ecore, eci, e0
      double precision citol, citol_now
      double precision tol2e, gnorm, gnrmv(100), cjgtol
      double precision cputime, walltime
      logical oe2approx, ostepci, oconj, ocanonical, onatural
      logical ocjrst, ocjrst1, ocanon, ocirestart
      logical ocimicro, ocivecprint
      logical oconverged
      logical oskel
      logical oprinfo, oprhigh, oprconv, oprciconv
      logical ohalf, oaodisk, oblk, oaokeep, oaoreuse      ! toggles for 4-index
      logical ocisolve
      integer ciiterlo, ciiterhi      ! Iteration range for CI solution
      integer k_tsym, l_tsym, k_evl, l_evl, k_ord, l_ord, iord_in, iorb
c
      double precision time_iter_start
      integer time_prev_iter
c     
      character*255 title, movecs_in
c     
      integer maxop, maxireps, nci
      parameter (maxop = 120, maxireps=20)
      integer nop, nir,  class_dim(maxireps)
      character*8 zir(maxireps), zclass(maxireps), zname
      double precision chars(maxireps*maxireps)
c
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
      integer ga_create_JKblocked
      external ga_create_JKblocked
      logical movecs_converged
      external movecs_converged
      logical int_normalize
      external int_normalize
      logical moints_aodisk_init, status, detci_memsiz
      external moints_aodisk_init
c
      character*1 mchar(9)
c     
      data mchar/'1','2','3','4','5','6','7','8','9'/
      data ohalf/.true./
      data blen/25/
c
c     This is a stupid kludge because of Benny breaking encapsulation rules.
c     Store SCF hamiltonian type as HF for use in BGJ routines
      if (.not. rtdb_put(rtdb, 'bgj:scf_type', MT_INT, 1, 1))
     $     call errquit('mcscf: put of bgj:scf_type failed',0,
     &       RTDB_ERR)
c     
      cputime = util_cpusec()
      walltime= util_wallsec()
c     
      mcscf = .false.
      oconverged = .false.
c     
      call util_print_push
      call util_print_rtdb_load(rtdb, 'mcscf')
      call ecce_print_module_entry('mcscf')
c     
      oprinfo = (util_print('information',print_low).and.
     $     (ga_nodeid().eq.0))
      oprconv = (util_print('conv',print_default).and.
     $     (ga_nodeid().eq.0))
      oprciconv = (util_print('ciconv',print_high).and.
     $     (ga_nodeid().eq.0))
c     
c     Get w.f. parameters
c     
      call mcscf_getinfo(rtdb, geom, basis, nbf, nclosed, nact )
c
      if (.not. sym_abelian_group(geom)) call errquit
     $     ('mcscf: the molecular point group is not abelian',0,
     &       GEOM_ERR)
c
      nvir = nbf - nclosed - nact
      orlen = (nclosed*nvir) + (nact*nvir) + (nclosed*nact)
      nsym = sym_number_ops(geom)+1
c
*      if (.not.rtdb_get(rtdb,'scf:skeleton',MT_LOG,1,oskel))
*     $     oskel = sym_number_ops(geom).gt.0
      oskel = sym_number_ops(geom).gt.0
c
      if (.not.geom_nuc_rep_energy( geom, enrep ))
     $     call errquit('mcscf: cannot retrieve nuclear repulsion',0,
     &       GEOM_ERR)
c
      if (.not.rtdb_get(rtdb,'mcscf:tol2e',MT_DBL,1,tol2e))        
     $     tol2e = 1.d-9       ! Redundant recovered later
c
      call sym_group_name(geom,zname)
      if (.not. sym_char_table(zname, nop, nir, class_dim,
     $     zir, zclass, chars)) call errquit
     $     ('mcscf: failed to get character table info',0,
     &       GEOM_ERR)
c     
c     Get electron and spin multiplicity (NB: for info only) and symm
c     Active elec and multiplicity must be set --- no defaults
c     
      if (rtdb_cget(rtdb,'mcscf:state',1,state)) then
        do multip = 1,9
	   if (state(1:1).eq.mchar(multip)) goto 656
        enddo
	call errquit('mcscf: failed to parse multiplicity in state',0,
     &       INPUT_ERR)
656     do symstate = 1,nir
	  if (inp_compare(.false.,state(2:),zir(symstate))) goto 657
        enddo
	call errquit('mcscf: failed to recognize irrep in state',0,
     &       INPUT_ERR)
 657    if (.not.rtdb_put(rtdb,'mcscf:statesymmetry',MT_INT,1,symstate))
     $     call errquit('mcscf: rtdb put of state sym?',0, RTDB_ERR)
        if (.not.rtdb_put(rtdb,'mcscf:multiplicity',MT_INT,1,multip))
     $     call errquit('mcsf: rtdb put of spin multiplicity ?',0,
     &       RTDB_ERR)
      endif
      if (.not.rtdb_get(rtdb,'mcscf:multiplicity',MT_INT,1,multip))
     $   call errquit('mcsf: rtdb get of spin multiplicity ?',0,
     &       RTDB_ERR)
      if (.not.rtdb_get(rtdb,'mcscf:statesymmetry',MT_INT,1,symstate))
     $   symstate = 1
      if (.not.rtdb_get(rtdb,'mcscf:nactelec',MT_INT,1,nactel))
     $     call errquit('number of active electrons not set',0,
     &       RTDB_ERR)
      nela = (nactel + multip - 1)/2
      nelb = nactel - nela
      if ((mod((nactel + multip - 1),2).ne.0).or.
     $     (nela.lt.0).or.(nelb.lt.0))
     $     call errquit('mcscf: incompatible nelec and spin',0,
     &       INPUT_ERR)
c     
c     Get convergence info
c     
      if (.not.rtdb_get(rtdb,'mcscf:maxiter',MT_INT,1,maxiter))
     $     maxiter = 20
      if (.not.rtdb_get(rtdb,'mcscf:cgreset',MT_INT,1,icgreset))  
     $     icgreset = 20
      if (.not.rtdb_get(rtdb,'mcscf:thresh',MT_DBL,1,thresh))
     $     thresh = 1.d-4
      if (.not. rtdb_get(rtdb,'mcscf:lshift',MT_DBL,1,lshift0))
     $     lshift0 = 0.1d0
      if (.not.rtdb_get(rtdb,'mcscf:e2approx',MT_LOG,1,oe2approx))
     $     oe2approx = .false.
      if (.not.rtdb_get(rtdb,'mcscf:cisolve',MT_LOG,1,ocisolve))
     $     ocisolve = .true.
      if (.not.rtdb_get(rtdb,'mcscf:conjugacy',MT_LOG,1,oconj))
     $     oconj = .true.
      if (.not.rtdb_get(rtdb,'mcscf:microci',MT_LOG,1,ostepci))
     $     ostepci = .true.
      if (.not.rtdb_get(rtdb,'mcscf:citol',MT_DBL,1,citol))
     $     citol = min(1.d-6,thresh)
      if (.not. rtdb_get(rtdb,'mcscf:canonical',MT_LOG,1,ocanonical))
     $     ocanonical = .true.
      if (.not. rtdb_get(rtdb,'mcscf:natural',MT_LOG,1,onatural))
     $     onatural = .true.
      if (.not. rtdb_get(rtdb,'mcscf:ciiterlo',MT_INT,1,ciiterlo))
     $     ciiterlo = 0
      if (.not. rtdb_get(rtdb,'mcscf:ciiterhi',MT_INT,1,ciiterhi))
     $     ciiterhi = maxiter + 1
      if (.not. rtdb_cget(rtdb, 'title', 1, title))
     $     title = ' '
c     
c     Profiling and statistics
c     
      call mcscf_profile_init(rtdb)
c     
c     Print info
c     
      oprhigh = util_print('extra info',print_high)
      if (oprinfo) then
         write(6,*)
         call util_print_centered(6, 'NWChem Direct MCSCF Module',
     $        40, .true.)
         write(6,*)
         write(6,*)
         if (title .ne. ' ') then
            call util_print_centered(6, title, 40, .false.)
            write(6,*)
            write(6,*)
         endif
         write(6,901) nbf, nclosed, nact, nactel,zir(symstate),
     $        multip, orlen
         if (oprhigh) write(6,902) (nclosed*nvir),(nact*nvir),
     $        (nclosed*nact)
 901     format(5x,35('-'),/,
     $        7x,'Basis functions:',10x,i5,/,
     $        7x,'Inactive shells:',10x,i5,/,
     $        7x,'Active shells:',12x,i5,/,
     $        7x,'Active electrons:',9x,i5,/,
     $        7x,'Symmetry:    ',16x,a,/,
     $        7x,'Multiplicity:',13x,i5,/,
     $        7x,'Orbital rotations:',8x,i5)
 902     format(12x,'Inact - Virt',9x,i5,/,
     $        12x,'Act - Virt',11x,i5,/,
     $        12x,'Inact - Act',10x,i5,/)
         write(6,903)
 903     format(5x,35('-'))
         write(6,904) oe2approx, oconj, ostepci, ocanonical,
     $        maxiter, thresh, icgreset, lshift0
 904     format(/,5x,35('-'),/,
     $        7x,'E2 extrapolation:',9x,l5,/,
     $        7x,'Conjugacy:',16x,l5,/,
     $        7x,'Micro CI relax:',11x,l5,/,
     $        7x,'Canonical orbitals:',7x,l5,/,
     $        7x,'Maximum iterations:',7x,i5,/,
     $        7x,'Converge threshold:',1p,e12.2,0p/,
     $        7x,'Conjugacy resets:',9x,i5,/,
     $        7x,'Initial levelshift:',f12.2,/,
     $        5x,35('-'))
      endif
c     
      if (nvir.lt.0) then
        call errquit("mcscf: negative number of virtual orbs",
     &               nvir,INPUT_ERR)
      endif
      if (.not.int_normalize(rtdb,basis))
     &    call errquit('mcscf:int_normalize failed',0, INT_ERR)
c     
c     See if we are already converged

      if (.not.rtdb_cget(rtdb,'mcscf:input vectors',1,movecs_in)) then
         if (.not.rtdb_cget(rtdb,'scf:input vectors',1,movecs_in))
     $        movecs_in = 'atomic'
      endif
      if (.not. (movecs_in.eq.'atomic'.or.movecs_in.eq.'hcore')) then
         call util_file_name_resolve(movecs_in, .false.)
         if (movecs_converged(rtdb, basis, 'mcscf', movecs_in)) then
            if (.not. rtdb_get(rtdb, 'mcscf:energy', mt_dbl, 1, energy))
     $           call errquit('mcscf:failed getting converged energy',0,
     &       RTDB_ERR)
            if (ga_nodeid().eq.0 .and. oprinfo) then
               write(6,1101) energy
 1101          format(/'  The MCSCF is already converged '//,
     $              '       Total MCSCF energy =', f20.12/)
               call util_flush(6)
            endif
            oconverged = .true.
            mcscf = .true.
            goto 401
         endif
      endif
c     
c     Create occupation and orbital energy vectors
c     
      if (.not.ma_push_get(MT_DBL, nbf, 'eps', l_evals, k_evals))
     $     call errquit('mcscf: cannot allocate',0, MA_ERR)
      if (.not.ma_push_get(MT_DBL, nbf, 'MO occ', l_occ, k_occ))
     $     call errquit('mcscf: cannot allocate',0, MA_ERR)
c     
c     Allocate 1- & 2-PDM
c     
      if (.not.ma_push_get(MT_DBL, (nact*nact*nact*nact),
     $     '2P density', l_dm2, k_dm2))
     $     call errquit('mcscf: cannot allocate MO density',0, MA_ERR)
      if (.not.ma_push_get(MT_DBL, (nact*nact),
     $     '1P density', l_dm1, k_dm1))
     $     call errquit('mcscf: cannot allocate MO density',0, MA_ERR)
c     
c     Initialize stuff
c     
      call int_init(rtdb, 1, basis)
      call schwarz_init(geom, basis)
      call scf_get_fock_param( rtdb, tol2e )
c     
c     Create and get MO vectors 
c     
*ga:1:0
      if (.not.ga_create(MT_DBL,nbf,nbf,'MO vectors',nbf,0,g_movecs))
     $     call errquit('mcscf: cannot allocate MO vectors',0, GA_ERR)

      call mcscf_getmovecs( rtdb, geom, basis, nbf, dbl_mb(k_occ),
     $     dbl_mb(k_evals), g_movecs )
c 
c     Orbital irreps
c     
      if (.not.ma_push_get(MT_INT, nbf, 'MO sym', l_sym, k_sym))
     $     call errquit('mcscf: cannot allocate symmetry',0, MA_ERR)
      call sym_movecs_adapt( basis, 1.d-8, g_movecs,
     $     int_mb(k_sym), nmixed )
      if (nmixed .ne. 0) call errquit(
     $     'mcscf: symmetry contamination in starting MOs', nmixed,
     &       GEOM_ERR)
      call ga_orthog_mos(basis, g_movecs)
c
      if (.not.ga_duplicate(g_movecs,g_hcore,'1e hamiltonian'))
     $     call errquit('mcscf: cannot dup MO vectors',0, GA_ERR)
c
      status = detci_memsiz(nact, nela, nelb, nsym,
     $     int_mb(k_sym+nclosed), nci )
      if (oprinfo) write(6,7761) nci
 7761 format(7x,'No. of determinants:',6x,i8/)
c
c     Order orbitals based in CAS input
c
      if (.not.ma_push_get(MT_INT, nact, 'order', l_ord, k_ord))
     $     call errquit('mcscf: cannot allocate order',0, MA_ERR)
      iord_in=-1
      do i=1,nir
         do j=nclosed+1,nclosed+nact
            if (int_mb(k_sym+j-1).eq.i) then
               iord_in = iord_in + 1
               int_mb(k_ord+iord_in) = j
            endif
         enddo
      enddo
      call ga_copy(g_movecs,g_hcore)
      if (.not.ma_push_get(MT_DBL, nact,'tmp evals', l_evl, k_evl))
     $   call errquit('mcscf_reorder: cannot allocate tmp evals',0,
     &                 MA_ERR)
      call dcopy(nact,dbl_mb(k_evals+nclosed),1,dbl_mb(k_evl),1)
      if (.not.ma_push_get(MT_INT, nact,'tmp sym', l_tsym, k_tsym))
     $    call errquit('mcscf_reorder: cannot allocate tmp sym',0,
     &                 MA_ERR)
      call icopy(nact,int_mb(k_sym+nclosed),1,int_mb(k_tsym),1)
      do i=1,nact 
         iorb=int_mb(k_ord+i-1)
         call ga_copy_patch('N',g_hcore,1,nbf,iorb,iorb,
     &                          g_movecs,1,nbf,nclosed+i,nclosed+i)
         dbl_mb(k_evals+nclosed+i-1) = dbl_mb(k_evl+iorb-nclosed-1)
         int_mb(k_sym+nclosed+i-1)   = int_mb(k_tsym+iorb-nclosed-1)
      enddo
      if (.not.ma_pop_stack(l_tsym))
     $     call errquit('mcscf: cannot pop temp sym',0, MA_ERR)
      if (.not.ma_pop_stack(l_evl))
     $     call errquit('mcscf: cannot pop tmp evls',0, MA_ERR)
      if (.not.ma_pop_stack(l_ord))
     $     call errquit('mcscf: cannot pop order',0, MA_ERR)
c     
c     Print orbital info
c     
      if (util_print('orbinfo',print_default).and.ga_nodeid().eq.0) then
         write(6,550)
 550     format(/,2x,'Starting Orbital Energies')
         write(6,551) (dbl_mb(k_evals+i-1),i=1,nbf)
 551     format(6f13.6)
         write(6,887) 
 887     format(/,2x,'Orbital Symmetry Irreps')
         write(6,888) (zir(int_mb(k_sym+i)),i=0,nbf-1)
 888     format(15(a4,1x))
      endif
c     
c     Create conjugacy vectors
c     
c     - If we use the MCSCF as a FullCI program orlen will be 0
c       causing ga_create to fail. So we need to be at bit creative...
*ga:1:0
      if (.not.ga_create(MT_DBL,max(1,orlen),1,'Conj',0,0,g_conj))
     $     call errquit('mcscf: cannot allocate',0, GA_ERR)
*ga:1:0
      if (.not.ga_create(MT_DBL,max(1,orlen),1,'Old Gradient',0,0,
     $                   g_gradold))
     $     call errquit('mcscf: cannot allocate',0, GA_ERR)
      call ga_zero(g_conj)
      call ga_zero(g_gradold)
      pfac = 1.d0
      ocjrst = .false.
c     
c     Allocate operator matrices
c     Memory test required here!
c     
      mo_lo = nclosed + 1
      mo_hi = nclosed + nact
      noper = (nact*(nact+1))/2
      g_coul = ga_create_JKblocked(noper,nbf,nbf,'Coulomb Oper')
      g_exch = ga_create_JKblocked(noper,nbf,nbf,'X Oper')
      if (.not.rtdb_get(rtdb,'mcscf:aoblock',MT_LOG,1,oblk)) then
         if (.not.rtdb_get(rtdb,'fourindex:aoblock',MT_LOG,1,oblk))
     $        oblk = .true.
      endif
c
c  4-Index disk caching info
c
      if (.not.rtdb_get(rtdb, 'mcscf:aodisk', MT_LOG, 1, oaodisk)) then
        if (.not.rtdb_get(rtdb,'fourindex:aodisk', MT_LOG,
     $                    1, oaodisk)) oaodisk = .true.
      endif
      if (.not.rtdb_get(rtdb, 'mcscf:aoreuse',MT_LOG,1,oaoreuse)) then
        if (.not.rtdb_get(rtdb,'fourindex:oaoreuse', MT_LOG,
     $                    1, oaoreuse)) oaoreuse = .false.
      endif
      if (.not.rtdb_get(rtdb, 'mcscf:aokeep', MT_LOG, 1, oaokeep)) then
        if (.not.rtdb_get(rtdb,'fourindex:oaokeep', MT_LOG,
     $                    1, oaokeep)) oaokeep = .false.
      endif
      if (.not.moints_aodisk_init( oaodisk, oaoreuse ))
     $  call errquit('mcscf: cannot initialize 4-index AO cache',0,
     &       INT_ERR)
c     
c     Initial 4-Index Tranformation
c     
      call moints_build_2x( basis, ohalf, oskel,
     $     mo_lo, mo_lo, mo_hi, 1, nbf,
     $     g_movecs, g_coul, .true.,
     $     g_exch, .true., blen, oblk )
c     
c     Core energy
c     
      call dfill((nact*nact),0.d0,dbl_mb(k_dm1),1)
      call dfill((nact*nact*nact*nact),0.d0,dbl_mb(k_dm2),1)
      call mcscf_etrace( geom, basis, nbf, nclosed, nact,
     $     .false., oskel, tol2e, dbl_mb(k_dm1),
     $     dbl_mb(k_dm2), g_movecs, g_coul,
     $     eone, etwo, ecore )
c     
c     
c     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c     
c     MCSCF Macroiterations
c     
c     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c     
c     
      cjgtol = 1.d-3
      gnorm = 4.d0
      e0 = ecore + enrep
      ocanon = .false.
      ocivecprint = util_print('civec',print_debug)
      do iter=1,maxiter
         time_iter_start = util_wallsec()
         if (oprconv) write(6,996) iter
 996     format(///,10x,'%%%%% MacroIteration ',i3,' %%%%%',/)
         if (mod(iter,5) .eq. 0) call ga_orthog_mos(basis, g_movecs)
c     
c     CI logic
c     
         ocirestart = (iter.gt.1).and.(.not.(ocanon))
         ocisolve = ((iter.ge.ciiterlo).and.(iter.le.ciiterhi)).or.
     $        ((iter.eq.1).and.(ciiterlo.ne.0)).or.
     $        (ocanon)
c     
c     CI solver
c     
         if (ocisolve) then
c$$$  call mcscf_ifock( geom, basis, nbf, nclosed, nact,              
c$$$  $                       oskel, tol2e, g_movecs, eone, etwo,
c$$$  $                       ecore, g_hcore )                          ! this will give correct CI total energy
            call mcscf_fcore( basis, nbf, nclosed, nact, g_movecs,
     $           g_coul, g_exch, g_hcore ) ! wrong frozen core energy but saves 1 fock build
***   e0 = ecore + enrep  ! RJH ... report only valence CI number
            if (omcscfprof) call pstat_on(ps_ci)
            citol_now = min(citol, (gnorm*0.05d0))
            call mcscf_cisolve( rtdb, geom, basis, nbf, nclosed, nact,
     $           nsym, int_mb(k_sym), 0.0d0, ! Was e0 here 
     $           dbl_mb(k_evals), g_hcore, g_coul, 
     $           citol_now, ocirestart, oprciconv,
     $           ocivecprint,
     $           dbl_mb(k_dm1), dbl_mb(k_dm2), eci )
            if (util_print('ci energy',print_default).and.
     $           (ga_nodeid().eq.0)) write(6,882) eci
 882        format(/,11x,'>>>| Valence CI:',6x,f22.12)
            if (omcscfprof) call pstat_off(ps_ci)
         else
            if (oprinfo) write(6,834) 
 834        format(5x,'--- Skipped CI component ---')
         endif
c     
c     Conjugacy reset
c     
         ocjrst1 = (ocjrst).or.(mod(iter,icgreset).eq.0).or.
     $        ((gnorm.lt.cjgtol).and.(.not.(ocanonical)))
         if ((oconj).and.(ocjrst1)) then
            call ga_zero(g_conj)
            call ga_zero(g_gradold)
            cjgtol = gnorm*gnorm
         endif
         ocjrst = ocjrst1
c     
c     Microiterative CI relaxation logic
c     
         ocimicro = ostepci.and.(gnorm.lt.5.d0).and.
     $        (iter.ge.ciiterlo).and.(iter.le.ciiterhi)
         if ((ostepci).and.(.not.(ocimicro)).and.
     $        (oprinfo)) write(6,812)
 812     format(5x,'--- Micro CI disabled ---')
c     
c     Orbital solver
c     
         if (omcscfprof) call pstat_on(ps_orb)
         ocanon = ocanonical.and.(mod((iter-1),10).eq.0)
         call mcscf_orbsolve( rtdb, geom, basis, nbf, nclosed, nact,
     $        oskel, thresh, int_mb(k_sym),
     $        dbl_mb(k_dm1), dbl_mb(k_dm2), g_coul,
     $        g_exch, oconj, ocjrst, ocanon, ocimicro,
     $        g_conj, g_gradold, g_hcore, eone, etwo,
     $        ecore, gnorm, dbl_mb(k_evals),
     $        g_movecs, oconverged, oprconv,
     $        dbl_mb(k_occ))
c     
         energy = eone + etwo + enrep
         call ecce_print1('total energy', mt_dbl, energy, 1)
         call ecce_print1('orbital gradient norm', mt_dbl, gnorm, 1)
c     
         gnrmv(iter) =  gnorm
         ocjrst = ocjrst.or.ocanon
         if (omcscfprof) call pstat_off(ps_orb)
         if (oconverged) goto 27
c
         time_prev_iter = util_wallsec() - time_iter_start
         time_iter_start = util_wallsec()
         if (.not. util_test_time_remaining(rtdb,time_prev_iter*2))
     $        goto 27
c         
      enddo
c     
c     Print convergence status
c     
 27   continue
      if (oprinfo) then
         if (gnorm.lt.thresh) then
            write(6,928)
 928        format(////,5x,30('='),/,12x,'MCSCF Converged',/,5x,
     $           30('='))
         else
            write(6,929)
 929        format(//,5x,'--- maximum iterations reached ---')
         endif
      endif
c     
      oconverged = gnorm .le. thresh
      mcscf = oconverged
c     
c     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c     
c     Analysis and Printing Section
c     
c     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c     
c     
c     Canonical Orbitals
c     
      if ((oprinfo).and.(util_print('movecs',print_debug))) then
         if (ocanonical) then
            write(6,810)
 810        format(///,10x,'Canonical Orbitals and Eigenvalues')
         else
            write(6,811)
 811        format(///,10x,'Orbitals and Eigenvalues')
         endif
         if (.not.ma_push_get(MT_DBL, (nbf*(nclosed+nact)), 'tmp',
     $        l_tmp, k_tmp))
     $        call errquit('mcscf: cannot allocate local MO',0, MA_ERR)
         call ga_get(g_movecs, 1, nbf, 1, (nclosed+nact),
     $        dbl_mb(k_tmp), nbf)
         chi = 0
 32      clo = chi + 1
         chi = min((clo + 7),(nact+nclosed))
         write(6,*)
         write(6,881) (dbl_mb(k_evals+i-1),i=clo,chi)
         write(6,*)
         do i=1,nbf
            write(6,881) (dbl_mb(k_tmp+(j-1)*nbf+i-1),j=clo,chi)
         enddo
         if (chi.ne.(nact+nclosed)) goto 32
         if (.not.ma_pop_stack(l_tmp))
     $        call errquit('mcscf: cannot pop local MO',0, MA_ERR)
      else if ((oprinfo) .and. util_print('evals',print_default)) then
         if (ocanonical) then
            write(6,8101)
 8101       format(///,10x,'Canonical Eigenvalues')
         else
            write(6,8111)
 8111       format(///,10x,'Non-canonical Eigenvalues')
         endif
         call output(dbl_mb(k_evals), 1, nbf, 1, 1, nbf, 1, 1)
         write(6,*)
      endif
c     
c     Solve CI for Canonical orbitals
c     
      call moints_build_2x( basis, ohalf, oskel,
     $     mo_lo, mo_lo, mo_hi, 1, nbf,
     $     g_movecs, g_coul, .true.,
     $     g_exch, .true., blen, oblk )
      call mcscf_fcore( basis, nbf, nclosed, nact, g_movecs,
     $     g_coul, g_exch, g_hcore )
      ocivecprint = util_print('converged civec',print_default)
      call mcscf_cisolve( rtdb, geom, basis, nbf, nclosed, nact,
     $     nsym, int_mb(k_sym), e0, dbl_mb(k_evals),
     $     g_hcore, g_coul, 
     $     citol, .true., .false., ocivecprint,
     $     dbl_mb(k_dm1), dbl_mb(k_dm2), eci )
c
c     Save convergence status ... must do this before saving
c     final orbitals since date of entry is examined at restart
c     
         if (.not. rtdb_put(rtdb, 'mcscf:converged', MT_LOG, 1,
     $        oconverged))
     $        call errquit('mcscf: failed to store converged in rtdb',0,
     &       RTDB_ERR)

c     
c     Natural Orbitals
c     
      if (onatural) then
         call mcscf_naturalorb( nbf, nclosed, nact, dbl_mb(k_dm1),
     $        dbl_mb(k_occ), dbl_mb(k_evals),
     $        g_movecs )
c     
c     Again, for degenerate irreps that we must split, symmetry
c     may be contaminated.  Adapt and orthog one more time
c     
         call sym_movecs_adapt(basis, 1.d-8, g_movecs,
     $        int_mb(k_sym), nmixed)
         call ga_orthog_mos(basis, g_movecs)
         call movecs_fix_phase(g_movecs)
c
         call mcscf_save_vec(rtdb,geom,basis, nbf, 
     $        dbl_mb(k_evals), dbl_mb(k_occ), g_movecs)
c     
         if ((oprinfo).and.(util_print('movecs',print_debug))) then
            write(6,880)
 880        format(///,10x,'Natural Orbitals and Occupation')
            if (.not.ma_push_get(MT_DBL, (nbf*(nclosed+nact)), 'tmp',
     $           l_tmp, k_tmp))
     $           call errquit('mcscf: cannot allocate local MO',0,
     &       MA_ERR)
            call ga_get(g_movecs, 1, nbf, 1, (nclosed+nact),
     $           dbl_mb(k_tmp), nbf)
            chi = 0
 33         clo = chi + 1
            chi = min((clo + 7),(nact+nclosed))
            write(6,*)
            write(6,881) (dbl_mb(k_occ+i-1),i=clo,chi)
 881        format(8f12.6)
            write(6,*)
            do i=1,nbf
               write(6,881) (dbl_mb(k_tmp+(j-1)*nbf+i-1),j=clo,chi)
            enddo
            if (chi.ne.(nact+nclosed)) goto 33
            if (.not.ma_pop_stack(l_tmp))
     $           call errquit('mcscf: cannot pop local MO',0,
     &       MA_ERR)
         else if ((oprinfo) .and. util_print('nocc',print_default))then
            write(6,8801)
 8801       format(///,10x,'Natural orbital occupation numbers')
            call output(dbl_mb(k_occ), 1, nact+nclosed, 1, 1, nbf, 1, 1)
            write(6,*)
         endif
c     
c     Resolve CI for Natural Orbitals
c     
         call moints_build_2x( basis, ohalf, oskel,
     $        mo_lo, mo_lo, mo_hi, 1, nbf,
     $        g_movecs, g_coul, .true.,
     $        g_exch, .true., blen, oblk )
         call mcscf_fcore( basis, nbf, nclosed, nact, g_movecs,
     $        g_coul, g_exch, g_hcore )
         ocivecprint = util_print('civec',print_debug)
         call mcscf_cisolve( rtdb, geom, basis, nbf, nclosed, nact,
     $        nsym, int_mb(k_sym), e0, dbl_mb(k_evals),
     $        g_hcore, g_coul, 
     $        citol, .false., .false., ocivecprint,
     $        dbl_mb(k_dm1), dbl_mb(k_dm2), eci )
      endif
c     
c     Energy trace
c     
      call mcscf_etrace( geom, basis, nbf, nclosed, nact,
     $     .true., oskel, tol2e, dbl_mb(k_dm1),
     $     dbl_mb(k_dm2), g_movecs, g_coul,
     $     eone, etwo, ecore )
      energy = eone + etwo + enrep
c     
      cputime  = util_cpusec() - cputime
      walltime = util_wallsec() - walltime
      if (util_print('trace energy',print_low).and.
     $     (ga_nodeid().eq.0)) write(6,674) energy, gnorm, 
     $     cputime, walltime
 674  format(//
     $     '       Final MCSCF results '/
     $     '       ------------------- '//
     $     '       Total MCSCF energy =', f20.12/
     $     '            Gradient norm =', 1p,d13.1,0p//
     $     '    CPU time for solution =', f9.1,'s'/
     $     '   Wall time for solution =', f9.1,'s'//)
c     
      call ecce_print1('total energy', mt_dbl, energy, 1)
      call ecce_print1('orbital gradient norm', mt_dbl, gnorm, 1)
      call ecce_print1('two-electron energy', mt_dbl, etwo, 1)
      call ecce_print1('one-electron energy', mt_dbl, eone, 1)
      call ecce_print1('nuclear repulsion energy', mt_dbl, enrep, 1)
      call movecs_ecce(nbf, nbf, 1, min(nclosed+nact+20,nbf),
     $     dbl_mb(k_evals),
     $     dbl_mb(k_occ), int_mb(k_sym), 
     $     g_movecs, 'mcscf', ' ')
c
      if (util_print('final vectors analysis', print_default))
     $     call movecs_print_anal(basis, 
     $     max(1,nclosed-10), min(nbf,nclosed+nact+10),
     $     0.15d0, g_movecs, 'Analysis of MCSCF natural orbitals',
     $     .true., dbl_mb(k_evals), .true., int_mb(k_sym), 
     $     .true., dbl_mb(k_occ))
c     
      call mcscf_dump_for_grad(rtdb, geom, basis, nbf, nclosed, nact,
     $     oskel, tol2e, dbl_mb(k_dm1), dbl_mb(k_dm2), 
     $     g_movecs,  g_coul)

#ifdef MCSCF_DEBUGGER
c     
c     Debugging routine
c     
      call mcscf_debugger( rtdb, basis, geom, nbf, nclosed, nact,
     $     g_movecs, dbl_mb(k_evals) )
#endif
c     
c     Store energy and convergence status
c     
      if (.not. rtdb_put(rtdb, 'mcscf:energy', MT_DBL, 1, energy)) 
     $     call errquit('mcscf: failed to store energy in rtdb', 0,
     &       RTDB_ERR)

      call moints_aodisk_close(oaokeep)
c     
c     Cleanup
c     
      call scf_tidy(rtdb)
      call schwarz_tidy()
      call int_terminate()
      if (.not.ga_destroy(g_hcore))
     $     call errquit('mcscf: cannot destroy 1e hamiltonian',0,
     &       GA_ERR)
      if (.not.ga_destroy(g_conj))
     $     call errquit('mcscf: cannot destroy conjugacy',0,
     &       GA_ERR)
      if (.not.ga_destroy(g_gradold))
     $     call errquit('mcscf: cannot destroy old gradient',0, GA_ERR)
      if (.not.ma_pop_stack(l_sym))
     $     call errquit('mcscf: cannot pop stack?',0, MA_ERR)
      if (.not.ma_pop_stack(l_dm1))
     $     call errquit('mcscf: cannot pop stack?',0, MA_ERR)
      if (.not.ma_pop_stack(l_dm2))
     $     call errquit('mcscf: cannot pop stack?',0, MA_ERR)
      if (.not.ma_pop_stack(l_occ))
     $     call errquit('mcscf: cannot pop stack?',0, MA_ERR)
      if (.not.ma_pop_stack(l_evals))
     $     call errquit('mcscf: cannot pop stack?',0, MA_ERR)
      if (.not.ga_destroy(g_exch))
     $     call errquit('mcscf: cannot destroy exchange',0, GA_ERR)
      if (.not.ga_destroy(g_coul))
     $     call errquit('mcscf: cannot destroy Coulomb',0, GA_ERR)
      if (.not.ga_destroy(g_movecs))
     $     call errquit('mcscf: cannot destroy MO vectors',0, GA_ERR)
c
 401  call mcscf_profile_tidy()
      if (.not.bas_destroy(basis))
     $     call errquit('mcscf: cannot free basis handle',0, BASIS_ERR)
      if (.not.geom_destroy(geom))
     $     call errquit('mcscf: cannot destory geom handle', 0,
     &       GEOM_ERR)
c     
      if (oconverged) then
         call ecce_print_module_exit('mcscf', 'ok')
      else
         call ecce_print_module_exit('mcscf', 'failed')
      endif
c     
      call util_print_pop
      return
      end







       subroutine mcscf_getinfo( rtdb, geom, basis, nbf, nclosed, nact )
       implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "geom.fh"
c
       integer rtdb
       integer geom
       integer basis
       integer nbf
       integer nclosed
       integer nact
c
       integer nactel
       integer nelec
       double precision charge
       double precision nuclear_charge
c
c  Create geom and basis handles
c
       if (.not. geom_create(geom, 'geometry'))
     $     call errquit('get_db_handles: geom_create?', 0, GEOM_ERR)
       if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     $     call errquit('get_db_handles: no geometry ', 0, RTDB_ERR)
       if (.not. bas_create(basis, 'ao basis'))
     $     call errquit('get_db_handles: bas_create?', 0, BASIS_ERR)
       if (.not. bas_rtdb_load(rtdb, geom, basis, 'ao basis'))
     $     call errquit('get_db_handles: no ao basis ', 0,
     &       RTDB_ERR)
c
c  Basis info
c
      if (.not. bas_numbf(basis, nbf)) call errquit
     $     ('scf_init: basis info',0, BASIS_ERR)
c
c  Find number of electrons
c
      if (.not. rtdb_get(rtdb, 'charge', MT_DBL, 1, charge))
     $     charge = 0.0d0
      if (.not. geom_nuc_charge(geom, nuclear_charge))
     $     call errquit('scf: geom_nuc_charge failed', 0, GEOM_ERR)
      nelec = nint(nuclear_charge - charge)
      if (nelec .le. 0) call errquit
     $     ('scf: negative no. of electrons ?', nelec, INPUT_ERR)
      if (abs(nuclear_charge - charge - dble(nelec)) .gt. 1.0d-8)
     $     call errquit('mcscf: charge neglected?', 0, INPUT_ERR)
c
c  Determine no. of closed and active shells
c  **** Have to correctly determine number of active from MCSCF input ****
c
       if (.not. rtdb_get(rtdb, 'mcscf:nact', MT_INT, 1, nact))
     $    call errquit('mcscf:no active orbitals',0, RTDB_ERR)
       if (.not. rtdb_get(rtdb, 'mcscf:nactelec', MT_INT, 1, nactel))
     $    call errquit('mcscf:no active electrons',0, RTDB_ERR)
       if (nactel.gt.(2*nact))
     $    call errquit('mcscf: actel > 2*actorb',nactel, UNKNOWN_ERR)
       nclosed = (nelec-nactel) / 2
       if (.not. rtdb_put(rtdb, 'mcscf:nclosed', mt_int, 1, nclosed))
     $      call errquit('mcscf: failed storing nclosed', 0, RTDB_ERR)
c

       end










      subroutine mcscf_profile_init(rtdb)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "pstat.fh"
#include "mcscfprof.fh"
c
      integer rtdb
      integer junk
c
c  Check if General profiling enabled
c
      if (.not.rtdb_get(rtdb,'mcscf:prof',MT_LOG,1,omcscfprof))
     $     omcscfprof = .false.
      if (.not.(omcscfprof)) return
c
c  If Fock profiling enabled call SCF initialization
c
      if (.not.rtdb_get(rtdb,'mcscf:fockprof',MT_LOG,1,ofockprof)) then
        if (.not.rtdb_get(rtdb,'scf:pstat',MT_LOG,1,ofockprof))
     $     ofockprof = .false.
      endif
      if (ofockprof) then
        if (.not.rtdb_put(rtdb,'scf:pstat',MT_LOG,1,.true.))
     $     call errquit('mcscf_profile: rtdb_put of scf:stat failed',0,
     &       RTDB_ERR)
        call scf_pstat_init(rtdb)
      endif
c
c  Overall MCSCF profiling
c
      if (.not. pstat_allocate('mcscf:orbsolve', pstat_qstat, 0,junk,
     $     ps_orb)) call errquit('mcscf: scf_pstat_init', 0, MEM_ERR)
      if (.not. pstat_allocate('mcscf:ci', pstat_qstat, 0,junk,
     $     ps_ci)) call errquit('mcscf: scf_pstat_init', 0, MEM_ERR)
      if (.not. pstat_allocate('mcscf:4index', pstat_qstat, 0,junk,
     $     ps_4index)) call errquit('mcscf: scf_pstat_init', 0, MEM_ERR)
      if (.not. pstat_allocate('mcscf:prcnd', pstat_qstat, 0,junk,
     $     ps_prcnd)) call errquit('mcscf: scf_pstat_init', 0, MEM_ERR)
      if (.not. pstat_allocate('mcscf:hess',
     $     pstat_qstat, 0, junk, ps_hv))
     $     call errquit('mcscf: scf_pstat_init', 0, MEM_ERR)
      if (.not. pstat_allocate('mcscf:hess:1e',
     $     pstat_qstat, 0, junk, ps_hv1e))
     $     call errquit('mcscf: scf_pstat_init', 0, MEM_ERR)
      if (.not. pstat_allocate('mcscf:hess:2e',
     $     pstat_qstat, 0, junk, ps_hv2e))
     $     call errquit('mcscf: scf_pstat_init', 0, MEM_ERR)
      if (.not. pstat_allocate('mcscf:hess:2e-ao',
     $     pstat_qstat, 0, junk, ps_hv2eao))
     $     call errquit('mcscf: scf_pstat_init', 0, MEM_ERR)
      if (.not. pstat_allocate('mcscf:hess:2e-mo',
     $     pstat_qstat, 0, junk, ps_hv2emo))
     $     call errquit('mcscf: scf_pstat_init', 0, MEM_ERR)
      if (.not. pstat_allocate('mcscf:fock',
     $     pstat_qstat, 0, junk, ps_gfock))
     $     call errquit('mcscf: scf_pstat_init', 0, MEM_ERR)
c
c Detailed CI profiling
c
      if (.not.rtdb_get(rtdb,'mcscf:ciprof',MT_LOG,1,ociprof))
     $   ociprof = .false.
      if (ociprof) then
        if (.not. pstat_allocate('mcscf:cidiag',
     $     pstat_qstat, 0, junk, ps_cidiag))
     $     call errquit('mcscf: scf_pstat_init', 0, MEM_ERR)
        if (.not. pstat_allocate('mcscf:cidens',
     $     pstat_qstat, 0, junk, ps_cidens))
     $     call errquit('mcscf: scf_pstat_init', 0, MEM_ERR)
      endif

      
      return
      end




      subroutine mcscf_profile_tidy()
      implicit none
#include "errquit.fh"
#include "pstat.fh"
#include "mcscfprof.fh"
c
      if (.not.(omcscfprof)) return
      if (ofockprof) call scf_pstat_print()
      if (.not.pstat_free(ps_orb)) call errquit('mcscf_pstat?',0,
     &       MEM_ERR)
      if (.not.pstat_free(ps_ci)) call errquit('mcscf_pstat?',0,
     &       MEM_ERR)
      if (.not.pstat_free(ps_4index)) call errquit('mcscf_pstat?',0,
     &       MEM_ERR)
      if (.not.pstat_free(ps_prcnd)) call errquit('mcscf_pstat?',0,
     &       MEM_ERR)
      if (.not.pstat_free(ps_hv)) call errquit('mcscf_pstat?',0,
     &       MEM_ERR)
      if (.not.pstat_free(ps_hv1e)) call errquit('mcscf_pstat?',0,
     &       MEM_ERR)
      if (.not.pstat_free(ps_hv2e)) call errquit('mcscf_pstat?',0,
     &       MEM_ERR)
      if (.not.pstat_free(ps_hv2eao)) call errquit('mcscf_pstat?',0,
     &       MEM_ERR)
      if (.not.pstat_free(ps_hv2emo)) call errquit('mcscf_pstat?',0,
     &       MEM_ERR)
      if (ociprof) then
        if (.not.pstat_free(ps_cidiag)) call errquit('mcscf_pstat?',0,
     &       MEM_ERR)
        if (.not.pstat_free(ps_cidens)) call errquit('mcscf_pstat?',0,
     &       MEM_ERR)
      endif
      return
      end



