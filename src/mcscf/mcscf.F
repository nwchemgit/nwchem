       subroutine mcscf( rtdb )
       implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "context.fh"
#include "util.fh"
#include "sym.fh"
#include "pstat.fh"
#include "mcscfprof.fh"
c
       integer rtdb
c
       integer geom, basis
       integer nbf, nclosed, nvir, nact, noper, nsym
       integer nactel, nela, nelb, multip, orlen
       integer l_occ, k_occ, l_sym, k_sym
       integer l_evals, k_evals
       integer l_dm1, k_dm1, l_dm2, k_dm2
       integer l_tmp, k_tmp
       integer g_movecs, g_coul, g_exch
       integer g_conj, g_gradold, g_hcore
       integer mo_lo, mo_hi
       integer i, j, iter, maxiter, nmixed, clo, chi
       integer icgreset                                   ! Number of iterations between CG resets
       integer blen
       double precision pfac
       double precision thresh, lshift0
       double precision eone, etwo, enrep, etrace
       double precision ecore, eci, e0
       double precision citol, citol_now
       double precision tol2e, gnorm, gnrmv(100), cjgtol
       logical oe2approx, ostepci, oconj, ocanonical, onatural
       logical ocjrst, ocjrst1, ocanon, ocirestart
       logical ocimicro, ocivecprint
       logical oconverged
       logical oskel
       logical oprinfo, oprhigh
       logical ohalf
       logical ocisolve
       logical oblk
       integer ciiterlo, ciiterhi                         ! Iteration range for CI solution
c
       integer ga_create_atom_blocked
       external ga_create_atom_blocked
       integer ga_create_JKblocked
       external ga_create_JKblocked
c
       data ohalf/.true./
       data blen/25/
c
c
c
       if (.not. context_push('mcscf'))
     $   call errquit('mcscf: context_push failed', 0)
       call util_print_push
       call util_print_rtdb_load(rtdb, 'mcscf')
       oprinfo = (util_print('information',print_default).and.
     $           (ga_nodeid().eq.0))
c
c  Get w.f. parameters
c
       call mcscf_getinfo(rtdb, geom, basis, nbf, nclosed, nact )
       nvir = nbf - nclosed - nact
       orlen = (nclosed*nvir) + (nact*nvir) + (nclosed*nact)
       nsym = sym_number_ops(geom)+1
       if (.not.rtdb_get(rtdb,'scf:skeleton',MT_LOG,1,oskel))
     $   oskel = sym_number_ops(geom).gt.0
       if (.not.geom_nuc_rep_energy( geom, enrep ))
     $   call errquit('mcscf: cannot retrieve nuclear repulsion',0)
       if (.not.rtdb_get(rtdb,'mcscf:tol2e',MT_DBL,1,tol2e))        
     $      tol2e = 1.d-12                                                   ! Redundant recovered later
c
c  Get electron and spin multiplicity (NB: for info only)
c  Active elec and multiplicity must be set --- no defaults
c
       if (.not.rtdb_get(rtdb,'mcscf:nactelec',MT_INT,1,nactel))
     $   call errquit('number of active electrons not set',0)
       if (.not.rtdb_get(rtdb,'mcscf:multiplicity',MT_INT,1,multip))
     $   call errquit('spin multiplicity not set',0)
       nela = (nactel + multip - 1)/2
       nelb = nactel - nela
       if ((mod((nactel + multip - 1),2).ne.0).or.
     $    (nela.lt.0).or.(nelb.lt.0))
     $    call errquit('mcscf: incompatible elec and spin',0)
c
c  Get convergence info
c
       if (.not.rtdb_get(rtdb,'mcscf:maxiter',MT_INT,1,maxiter)) then
         if (.not.rtdb_get(rtdb,'scf:maxiter',MT_INT,1,maxiter))
     $     maxiter = 5
       endif
       if (.not.rtdb_get(rtdb,'mcscf:cgreset',MT_INT,1,icgreset))  
     $   icgreset = 20
       if (.not.rtdb_get(rtdb,'mcscf:thresh',MT_DBL,1,thresh)) then
         if (.not.rtdb_get(rtdb,'scf:thresh',MT_DBL,1,thresh))
     $     thresh = 1.d-6
       endif
       if (.not. rtdb_get(rtdb,'mcscf:lshift',MT_DBL,1,lshift0)) then
         if (.not. rtdb_get(rtdb,'scf:lshift',MT_DBL,1,lshift0))
     $     lshift0 = 0.1d0
       endif
       if (.not.rtdb_get(rtdb,'mcscf:e2approx',MT_LOG,1,oe2approx))
     $      oe2approx = .true.
       if (.not.rtdb_get(rtdb,'mcscf:cisolve',MT_LOG,1,ocisolve))
     $      ocisolve = .true.
       if (.not.rtdb_get(rtdb,'mcscf:conjugacy',MT_LOG,1,oconj))
     $      oconj = .true.
       if (.not.rtdb_get(rtdb,'mcscf:microci',MT_LOG,1,ostepci))
     $      ostepci = .true.
       if (.not.rtdb_get(rtdb,'mcscf:citol',MT_DBL,1,citol))
     $      citol = 1.d-5
       if (.not. rtdb_get(rtdb,'mcscf:canonical',MT_LOG,1,ocanonical))
     $      ocanonical = .true.
       if (.not. rtdb_get(rtdb,'mcscf:natural',MT_LOG,1,onatural))
     $      onatural = .true.
       if (.not. rtdb_get(rtdb,'mcscf:ciiterlo',MT_INT,1,ciiterlo))
     $      ciiterlo = 0
       if (.not. rtdb_get(rtdb,'mcscf:ciiterhi',MT_INT,1,ciiterhi))
     $      ciiterhi = maxiter + 1
c
c  Profiling and statistics
c
       call mcscf_profile_init(rtdb)
c
c  Print info
c
       oprhigh = util_print('extra info',print_high)
       if (oprinfo) write(6,900)
 900   format(/,27x,14('='),/
     $          28x,'Direct MCSCF',/
     $          27x,14('='),/)
       if (oprinfo) write(6,901) nbf, nclosed, nact, nactel,
     $                           multip, orlen
       if (oprhigh) write(6,902) (nclosed*nvir),(nact*nvir),
     $                           (nclosed*nact)
 901   format(18x,35('-'),/,
     $        20x,'Basis functions:',10x,i5,/,
     $        20x,'Inactive shells:',10x,i5,/,
     $        20x,'Active shells:',12x,i5,/,
     $        20x,'Active electrons:',9x,i5,/,
     $        20x,'Multiplicity:',13x,i5,/,
     $        20x,'Orbital rotations:',8x,i5)
 902   format(25x,'Inact - Virt',9x,i5,/,
     $        25x,'Act - Virt',11x,i5,/,
     $        25x,'Inact - Act',10x,i5,/)
       if (oprinfo) write(6,903)
 903   format(18x,35('-'))
       if (oprinfo) write(6,904) oe2approx, oconj, ostepci, ocanonical,
     $                           maxiter, thresh, icgreset, lshift0
 904   format(/,18x,35('-'),/,
     $          20x,'E2 extrapolation:',9x,l5,/,
     $          20x,'Conjugacy:',16x,l5,/,
     $          20x,'Micro CI relax:',11x,l5,/,
     $          20x,'Canonical orbitals:',7x,l5,/,
     $          20x,'Maximum iterations:',7x,i5,/,
     $          20x,'Converge threshold:',e12.2,/,
     $          20x,'Conjugacy resets:',9x,i5,/,
     $          20x,'Initial levelshift:',f12.2,/,
     $          18x,35('-'))
c
c  Create occupation and orbital energy vectors
c
       if (.not.ma_push_get(MT_DBL, nbf, 'eps', l_evals, k_evals))
     $     call errquit('mcscf: cannot allocate',0)
       if (.not.ma_push_get(MT_DBL, nbf, 'MO occ', l_occ, k_occ))
     $     call errquit('mcscf: cannot allocate',0)
c
c  Allocate 1- & 2-PDM
c
       if (.not.ma_push_get(MT_DBL, (nact*nact*nact*nact),
     $                      '2P density', l_dm2, k_dm2))
     $     call errquit('mcscf: cannot allocate MO density',0)
       if (.not.ma_push_get(MT_DBL, (nact*nact),
     $                      '1P density', l_dm1, k_dm1))
     $     call errquit('mcscf: cannot allocate MO density',0)
c
c  Initialize stuff
c
       call int_normalize(basis)
       call int_init(rtdb, 1, basis)
       call schwarz_init(geom, basis)
       call scf_get_fock_param( rtdb, tol2e )
c
c  Create and get MO vectors 
c
       if (.not.ga_create(MT_DBL,nbf,nbf,'MO vectors',nbf,1,g_movecs))
     $   call errquit('mcscf: cannot allocate MO vectors',0)
       call mcscf_getmovecs( rtdb, geom, basis, nbf, dbl_mb(k_occ),
     $                       dbl_mb(k_evals), g_movecs )
c
c  Orbital irreps
c 
       if (.not.ma_push_get(MT_INT, nbf, 'MO sym', l_sym, k_sym))
     $     call errquit('mcscf: cannot allocate symmetry',0)
       call sym_movecs_adapt( basis, 1.d-8, g_movecs,
     $                        int_mb(k_sym), nmixed )
       if (nmixed .ne. 0) call errquit(
     $   'mcscf: symmetry contamination in starting MOs', nmixed)
c
c  Print orbital info
c
       if (oprinfo) then
         write(6,550)
 550     format(/,2x,'Starting Orbital Energies')
         write(6,551) (dbl_mb(k_evals+i-1),i=1,nbf)
 551     format(7f12.6)
         write(6,887) 
 887     format(/,2x,'Orbital Symmetry Irreps')
         write(6,888) (int_mb(k_sym+i),i=0,nbf-1)
 888     format(16i3)
       endif
       if (.not.ga_duplicate(g_movecs,g_hcore,'1e hamiltonian'))
     $   call errquit('mcscf: cannot dup MO vectors',0)
c
c   Create conjugacy vectors
c
       if (.not.ga_create(MT_DBL,orlen,1,'Conj',1,1,g_conj))
     $      call errquit('mcscf: cannot allocate',0)
       if (.not.ga_create(MT_DBL,orlen,1,'Old Gradient',1,1,g_gradold))
     $      call errquit('mcscf: cannot allocate',0)
       call ga_zero(g_conj)
       call ga_zero(g_gradold)
       pfac = 1.d0
       ocjrst = .false.
c
c   Allocate operator matrices
c   Memory test required here!
c
       mo_lo = nclosed + 1
       mo_hi = nclosed + nact
       noper = (nact*(nact+1))/2
       g_coul = ga_create_JKblocked(noper,nbf,nbf,'Coulomb Oper')
       g_exch = ga_create_JKblocked(noper,nbf,nbf,'X Oper')
       if (.not.rtdb_get(rtdb,'mcscf:aoblock',MT_LOG,1,oblk)) then
         if (.not.rtdb_get(rtdb,'fourindex:aoblock',MT_LOG,1,oblk))
     $      oblk = .false.
       endif
c
c   Initial 4-Index Tranformation
c
       call moints_build_2x( basis, ohalf, oskel,
     $                       mo_lo, mo_lo, mo_hi, 1, nbf,
     $                       g_movecs, g_coul, .true.,
     $                       g_exch, .true., blen, oblk )
c
c   Core energy
c
       call dfill((nact*nact),0.d0,dbl_mb(k_dm1),1)
       call dfill((nact*nact*nact*nact),0.d0,dbl_mb(k_dm2),1)
       call mcscf_etrace( geom, basis, nbf, nclosed, nact,
     $                    .false., oskel, tol2e, dbl_mb(k_dm1),
     $                    dbl_mb(k_dm2), g_movecs, g_coul,
     $                    eone, etwo, ecore )
c
c
c  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c                    MCSCF Macroiteractions
c
c  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c
       cjgtol = 1.d-3
       gnorm = 4.d0
       e0 = ecore + enrep
       ocanon = .false.
       ocivecprint = util_print('civec',print_debug)
       do iter=1,maxiter
         if (oprinfo) write(6,996) iter
 996     format(///,10x,'%%%%% MacroIteration ',i3,' %%%%%',/)
c
c  CI logic
c
         ocirestart = (iter.gt.1).and.(.not.(ocanon))
         ocisolve = ((iter.ge.ciiterlo).and.(iter.le.ciiterhi)).or.
     $              ((iter.eq.1).and.(ciiterlo.ne.0)).or.
     $              (ocanon)
c
c  CI solver
c
         if (ocisolve) then
c$$$           call mcscf_ifock( geom, basis, nbf, nclosed, nact,              
c$$$     $                       oskel, tol2e, g_movecs, eone, etwo,
c$$$     $                       ecore, g_hcore )                          ! this will give correct CI total energy
           call mcscf_fcore( basis, nbf, nclosed, nact, g_movecs,
     $                       g_coul, g_exch, g_hcore )                 ! wrong CI total energy but saves 1 fock build
           e0 = ecore + enrep
           if (omcscfprof) call pstat_on(ps_ci)
           citol_now = min(citol, (gnorm*0.05d0))
           call mcscf_cisolve( rtdb, geom, basis, nbf, nclosed, nact,
     $                         nsym, int_mb(k_sym), e0,
     $                         dbl_mb(k_evals), g_hcore, g_coul, 
     $                         citol_now, ocirestart, .true.,
     $                         ocivecprint,
     $                         dbl_mb(k_dm1), dbl_mb(k_dm2), eci )
           if (util_print('ci energy',print_default).and.
     $         (ga_nodeid().eq.0)) write(6,882) eci
 882       format(/,11x,'>>>| CI energy:',7x,f22.12)
           if (omcscfprof) call pstat_off(ps_ci)
         else
           if (oprinfo) write(6,834) 
 834       format(5x,'--- Skipped CI component ---')
         endif
c
c  Conjugacy reset
c
         ocjrst1 = (ocjrst).or.(mod(iter,icgreset).eq.0).or.
     $             ((gnorm.lt.cjgtol).and.(.not.(ocanonical)))
         if ((oconj).and.(ocjrst1)) then
           call ga_zero(g_conj)
           call ga_zero(g_gradold)
           cjgtol = gnorm*gnorm
         endif
         ocjrst = ocjrst1
c
c  Microiterative CI relaxation logic
c
         ocimicro = ostepci.and.(gnorm.lt.5.d0).and.
     $              (iter.ge.ciiterlo).and.(iter.le.ciiterhi)
         if ((ostepci).and.(.not.(ocimicro)).and.
     $       (oprinfo)) write(6,812)
 812     format(5x,'--- Micro CI disabled ---')
c
c  Orbital solver
c
         if (omcscfprof) call pstat_on(ps_orb)
         ocanon = ocanonical.and.(mod((iter-1),10).eq.0)
         call mcscf_orbsolve( rtdb, geom, basis, nbf, nclosed, nact,
     $                        oskel, thresh, int_mb(k_sym),
     $                        dbl_mb(k_dm1), dbl_mb(k_dm2), g_coul,
     $                        g_exch, oconj, ocjrst, ocanon, ocimicro,
     $                        g_conj, g_gradold, g_hcore, eone, etwo,
     $                        ecore, gnorm, dbl_mb(k_evals),
     $                        g_movecs, oconverged )
         gnrmv(iter) =  gnorm
         ocjrst = ocjrst.or.ocanon
         if (omcscfprof) call pstat_off(ps_orb)
         if (oconverged) goto 27
       enddo
c
c  Print convergence status
c
 27    continue
       if (oprinfo) then
         if (gnorm.lt.thresh) then
           write(6,928)
 928       format(////,5x,30('='),/,12x,'MCSCF Converged',/,5x,30('='))
         else
           write(6,929)
 929       format(//,5x,'--- maximum iterations reached ---')
         endif
       endif
c
c
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c             Analysis and Printing Section
c
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c
c  Canonical Orbitals
c
       if ((oprinfo).and.(util_print('movecs',print_default))) then
         if (ocanonical) then
           write(6,810)
 810       format(///,10x,'Canonical Orbitals and Eigenvalues')
         else
           write(6,811)
 811       format(///,10x,'Orbitals and Eigenvalues')
         endif
         if (.not.ma_push_get(MT_DBL, (nbf*(nclosed+nact)), 'tmp',
     $                        l_tmp, k_tmp))
     $        call errquit('mcscf: cannot allocate local MO',0)
         call ga_get(g_movecs, 1, nbf, 1, (nclosed+nact),
     $               dbl_mb(k_tmp), nbf)
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
     $     call errquit('mcscf: cannot pop local MO',0)
       endif
c
c  Solve CI for Canonical orbitals
c
       call moints_build_2x( basis, ohalf, oskel,
     $                       mo_lo, mo_lo, mo_hi, 1, nbf,
     $                       g_movecs, g_coul, .true.,
     $                       g_exch, .true., blen, oblk )
       call mcscf_fcore( basis, nbf, nclosed, nact, g_movecs,
     $                   g_coul, g_exch, g_hcore )
       ocivecprint = util_print('converged civec',print_default)
       citol = min(citol,gnorm*0.1d0)
       call mcscf_cisolve( rtdb, geom, basis, nbf, nclosed, nact,
     $                     nsym, int_mb(k_sym), e0, dbl_mb(k_evals),
     $                     g_hcore, g_coul, 
     $                     citol, .true., .false., ocivecprint,
     $                     dbl_mb(k_dm1), dbl_mb(k_dm2), eci )
c
c  Natural Orbitals
c
       if (onatural) then
         call mcscf_naturalorb( nbf, nclosed, nact, dbl_mb(k_dm1),
     $                          dbl_mb(k_occ), dbl_mb(k_evals),
     $                          g_movecs )
         if ((oprinfo).and.(util_print('movecs',print_default))) then
           write(6,880)
 880       format(///,10x,'Natural Orbitals and Occupation')
           if (.not.ma_push_get(MT_DBL, (nbf*(nclosed+nact)), 'tmp',
     $                          l_tmp, k_tmp))
     $          call errquit('mcscf: cannot allocate local MO',0)
           call ga_get(g_movecs, 1, nbf, 1, (nclosed+nact),
     $                dbl_mb(k_tmp), nbf)
           chi = 0
 33        clo = chi + 1
           chi = min((clo + 7),(nact+nclosed))
           write(6,*)
           write(6,881) (dbl_mb(k_occ+i-1),i=clo,chi)
 881       format(8f12.6)
           write(6,*)
           do i=1,nbf
             write(6,881) (dbl_mb(k_tmp+(j-1)*nbf+i-1),j=clo,chi)
           enddo
           if (chi.ne.(nact+nclosed)) goto 33
           if (.not.ma_pop_stack(l_tmp))
     $          call errquit('mcscf: cannot pop local MO',0)
         endif
c
c  Resolve CI for Natural Orbitals
c
         call moints_build_2x( basis, ohalf, oskel,
     $                         mo_lo, mo_lo, mo_hi, 1, nbf,
     $                         g_movecs, g_coul, .true.,
     $                         g_exch, .true., blen, oblk )
         call mcscf_fcore( basis, nbf, nclosed, nact, g_movecs,
     $                     g_coul, g_exch, g_hcore )
         ocivecprint = util_print('civec',print_debug)
         call mcscf_cisolve( rtdb, geom, basis, nbf, nclosed, nact,
     $                       nsym, int_mb(k_sym), e0, dbl_mb(k_evals),
     $                       g_hcore, g_coul, 
     $                       citol, .false., .false., ocivecprint,
     $                       dbl_mb(k_dm1), dbl_mb(k_dm2), eci )
       endif
c
c  Energy trace
c     
       call mcscf_etrace( geom, basis, nbf, nclosed, nact,
     $                    .true., oskel, tol2e, dbl_mb(k_dm1),
     $                    dbl_mb(k_dm2), g_movecs, g_coul,
     $                    eone, etwo, ecore )
       etrace = eone + etwo + enrep
       if (util_print('trace energy',print_default).and.
     $     (ga_nodeid().eq.0)) write(6,674) etrace, gnorm
 674   format(//,20x,'Energy:',11x,f20.14,
     $         /,20x,'Gradient Norm:',4x,e20.4)

#ifdef MCSCF_DEBUGGER
c
c  Debugging routine
c
       call mcscf_debugger( rtdb, basis, geom, nbf, nclosed, nact,
     $                      g_movecs, dbl_mb(k_evals) )
#endif
c
c  Cleanup
c
 401   continue
       call mcscf_profile_tidy()
       call scf_tidy()
       call schwarz_tidy()
       call int_terminate()
       if (.not.ga_destroy(g_hcore))
     $      call errquit('mcscf: cannot destroy 1e hamiltonian',0)
       if (.not.ga_destroy(g_conj))
     $      call errquit('mcscf: cannot destroy conjugacy',0)
       if (.not.ga_destroy(g_gradold))
     $      call errquit('mcscf: cannot destroy old gradient',0)
       if (.not.bas_destroy(basis))
     $   call errquit('mcscf: cannot free basis handle',0)
       if (.not.geom_destroy(geom))
     $   call errquit('mcscf: cannot destory geom handle', 0)
       if (.not.ma_pop_stack(l_sym))
     $      call errquit('mcscf: cannot pop stack?',0)
       if (.not.ma_pop_stack(l_dm1))
     $      call errquit('mcscf: cannot pop stack?',0)
       if (.not.ma_pop_stack(l_dm2))
     $      call errquit('mcscf: cannot pop stack?',0)
       if (.not.ma_pop_stack(l_occ))
     $      call errquit('mcscf: cannot pop stack?',0)
       if (.not.ma_pop_stack(l_evals))
     $      call errquit('mcscf: cannot pop stack?',0)
       if (.not.ga_destroy(g_exch))
     $      call errquit('mcscf: cannot destroy exchange',0)
       if (.not.ga_destroy(g_coul))
     $      call errquit('mcscf: cannot destroy Coulomb',0)
       if (.not.ga_destroy(g_movecs))
     $      call errquit('mcscf: cannot destroy MO vectors',0)
c
c
c
       call util_print_pop
       if (.not. context_pop('mcscf'))
     $   call errquit('context_pop failed',0)
       return
       end







       subroutine mcscf_getinfo( rtdb, geom, basis, nbf, nclosed, nact )
       implicit none
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
       character*255 stsym
       double precision charge
       double precision nuclear_charge
c
c  Create geom and basis handles
c
       if (.not. geom_create(geom, 'geometry'))
     $     call errquit('get_db_handles: geom_create?', 0)
       if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     $     call errquit('get_db_handles: no geometry ', 0)
       if (.not. bas_create(basis, 'ao basis'))
     $     call errquit('get_db_handles: bas_create?', 0)
       if (.not. bas_rtdb_load(rtdb, geom, basis, 'ao basis'))
     $     call errquit('get_db_handles: no ao basis ', 0)
c
c  Basis info
c
      if (.not. bas_numbf(basis, nbf)) call errquit
     $     ('scf_init: basis info',0)
c
c  Find number of electrons
c
      if (.not. rtdb_get(rtdb, 'charge', MT_DBL, 1, charge))
     $     charge = 0.0d0
      if (.not. geom_nuc_charge(geom, nuclear_charge))
     $     call errquit('scf: geom_nuc_charge failed', 0)
      nelec = nint(nuclear_charge - charge)
      if (nelec .le. 0) call errquit
     $     ('scf: negative no. of electrons ?', nelec)
      if (abs(nuclear_charge - charge - dble(nelec)) .gt. 1.0d-8)
     $     call errquit('mcscf: charge neglected?', 0)
c
c  Determine no. of closed and active shells
c  **** Have to correctly determine number of active from MCSCF input ****
c
       if (.not. rtdb_get(rtdb, 'mcscf:nact', MT_INT, 1, nact))
     $    call errquit('mcscf:no active orbitals',0)
       if (.not. rtdb_get(rtdb, 'mcscf:nactelec', MT_INT, 1, nactel))
     $    call errquit('mcscf:no active electrons',0)
       if (nactel.gt.(2*nact))
     $    call errquit('mcscf: actel > 2*actorb',nactel)
       nclosed = (nelec-nactel) / 2
c
c  Reconcile symmetry info
c
       if (.not.rtdb_cget(rtdb, 'mcscf:statesymmlabel', 1, stsym))
     $   stsym = 'a'

       return
       end










      subroutine mcscf_profile_init(rtdb)
      implicit none
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
     $     call errquit('mcscf_profile: rtdb_put of scf:stat failed',0)
        call scf_pstat_init(rtdb)
      endif
c
c  Overall MCSCF profiling
c
      if (.not. pstat_allocate('mcscf:orbsolve', pstat_qstat, 0,junk,
     $     ps_orb)) call errquit('mcscf: scf_pstat_init', 0)
      if (.not. pstat_allocate('mcscf:ci', pstat_qstat, 0,junk,
     $     ps_ci)) call errquit('mcscf: scf_pstat_init', 0)
      if (.not. pstat_allocate('mcscf:4index', pstat_qstat, 0,junk,
     $     ps_4index)) call errquit('mcscf: scf_pstat_init', 0)
      if (.not. pstat_allocate('mcscf:prcnd', pstat_qstat, 0,junk,
     $     ps_prcnd)) call errquit('mcscf: scf_pstat_init', 0)
      if (.not. pstat_allocate('mcscf:hess',
     $     pstat_qstat, 0, junk, ps_hv))
     $     call errquit('mcscf: scf_pstat_init', 0)
      if (.not. pstat_allocate('mcscf:hess:1e',
     $     pstat_qstat, 0, junk, ps_hv1e))
     $     call errquit('mcscf: scf_pstat_init', 0)
      if (.not. pstat_allocate('mcscf:hess:2e',
     $     pstat_qstat, 0, junk, ps_hv2e))
     $     call errquit('mcscf: scf_pstat_init', 0)
      if (.not. pstat_allocate('mcscf:hess:2e-ao',
     $     pstat_qstat, 0, junk, ps_hv2eao))
     $     call errquit('mcscf: scf_pstat_init', 0)
      if (.not. pstat_allocate('mcscf:hess:2e-mo',
     $     pstat_qstat, 0, junk, ps_hv2emo))
     $     call errquit('mcscf: scf_pstat_init', 0)
      if (.not. pstat_allocate('mcscf:fock',
     $     pstat_qstat, 0, junk, ps_gfock))
     $     call errquit('mcscf: scf_pstat_init', 0)
c
c Detailed CI profiling
c
      if (.not.rtdb_get(rtdb,'mcscf:ciprof',MT_LOG,1,ociprof))
     $   ociprof = .false.
      if (ociprof) then
        if (.not. pstat_allocate('mcscf:cidiag',
     $     pstat_qstat, 0, junk, ps_cidiag))
     $     call errquit('mcscf: scf_pstat_init', 0)
        if (.not. pstat_allocate('mcscf:cidens',
     $     pstat_qstat, 0, junk, ps_cidens))
     $     call errquit('mcscf: scf_pstat_init', 0)
      endif

      
      return
      end




      subroutine mcscf_profile_tidy()
      implicit none
#include "pstat.fh"
#include "mcscfprof.fh"
c
      if (.not.(omcscfprof)) return
      if (ofockprof) call scf_pstat_print()
      if (.not.pstat_free(ps_orb)) call errquit('mcscf_pstat?',0)
      if (.not.pstat_free(ps_ci)) call errquit('mcscf_pstat?',0)
      if (.not.pstat_free(ps_4index)) call errquit('mcscf_pstat?',0)
      if (.not.pstat_free(ps_prcnd)) call errquit('mcscf_pstat?',0)
      if (.not.pstat_free(ps_hv)) call errquit('mcscf_pstat?',0)
      if (.not.pstat_free(ps_hv1e)) call errquit('mcscf_pstat?',0)
      if (.not.pstat_free(ps_hv2e)) call errquit('mcscf_pstat?',0)
      if (.not.pstat_free(ps_hv2eao)) call errquit('mcscf_pstat?',0)
      if (.not.pstat_free(ps_hv2emo)) call errquit('mcscf_pstat?',0)
      if (ociprof) then
        if (.not.pstat_free(ps_cidiag)) call errquit('mcscf_pstat?',0)
        if (.not.pstat_free(ps_cidens)) call errquit('mcscf_pstat?',0)
      endif
      return
      end



