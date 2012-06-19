c  $Id$
c
c
c  Orbital MCSCF macroiteration
c
c      i)  Energy and gradient
c     ii)  Precondition using exact or approx. Hessian
c    iii)  Apply conjugacy to search vector
c     iv)  Line search on E2(T) surface / optional CI relaxation
c      v)  Update orbitals / optional canonicalization
c
       subroutine mcscf_orbsolve( rtdb, geom, basis, nbf, nclosed,
     $                           nact, oskel, gnorm_tol, orbsym,
     $                           dm1, dm2, g_coul, g_exch, oconj,
     $                           ocjrst, ocanon, ocirelax, g_conj,
     $                           g_gradold, g_hcore, eone, etwo, ecore,
     $                           gnorm, orbe, g_movecs, oconverged,
     $                           oprconv, occ )
       implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "sym.fh"
#include "util.fh"
#include "msgids.fh"
#include "mcscfprof.fh"
c
       integer rtdb, geom, basis                       ! [input] Handles
       integer nbf, nclosed, nact                      ! [input] Orbital parameters
       logical oskel                                   ! [input] Symmetry toggle
       double precision gnorm_tol                      ! [input] Required covergence tolerance
       integer orbsym(*)                               ! [input] Orbital irreps
       double precision dm1(*)                         ! [input] 1-PDM
       double precision dm2(*)                         ! [input] 2-PDM
       integer g_coul                                  ! [input] Coulomb operator
       integer g_exch                                  ! [input] Exchange operator
       logical oconj                                   ! [input] Conjugacy toggle
       logical ocjrst                                  ! [input/output] Conjugacy reset
       logical ocanon                                  ! [input] Canonicalization toggle
       logical ocirelax                                ! [input] Microiterative CI relax toggle
       integer g_conj                                  ! [input] Conjugacy vector
       integer g_gradold                               ! [input] Old gradient
       integer g_hcore                                 ! [input] 1e hamiltonian
       double precision eone                           ! [output] 1e-energy
       double precision etwo                           ! [output] 2e-energy
       double precision ecore                          ! [output] clsoed-shell core energy
       double precision gnorm                          ! [output] gradient norm
       double precision orbe(nbf)                      ! [output] MCSCF orbital energies
       double precision occ(nbf)                       ! [input] MCSCF orbital occupation (junk until natorb)
       integer g_movecs                                ! [input/output] MO vectors
       logical oconverged                              ! [output] convergence flag
       logical oprconv                                 ! [input] if true then print
c
       integer maxpiter
       integer maxlinsrch
       parameter(maxpiter=15)
       parameter(maxlinsrch=16)
c
       double precision tol_floor,tol_ceil
       parameter(tol_floor=1.d-12)
       parameter(tol_ceil=1.d-1)
c
       integer iter, vlen, pflg
       integer molo, mohi, i, iter1, nstep
       double precision tol2e, enrep, energy, lshift0
       double precision snorm, alpha, gsscale, gnorm0
       double precision eprec, p_acc, ndigit
       double precision eev(0:maxlinsrch), gsv(0:maxlinsrch)
       double precision av(0:maxlinsrch)
       double precision hmax, hmin, maxx
       double precision lsrch_tol
       integer g_afock, g_ifock, g_gfock
       integer g_newvecs, g_grad, g_kvec, g_srch, g_hdiag
       integer blen, nmixed, nflip, nsym, nnact
       integer g_u, g_b, g_tmp, g_actints
       logical ocoul, oexch, ohalf
       logical oe2approx
       logical odelta_core
       logical olsrchprint, ociprint
       logical ocivecprint, oprcndprint
       logical oblk, obshift
       double precision micro_citol, citol
       double precision lshift, bshift
       double precision e0, eci
       save lshift
       double precision mone, knorm
c
       logical line_search
       external line_search
c
       data lshift/1.0d0/
       data mone/-1.d0/
       data blen/25/
       data nflip/0/
       data ohalf/.true./
       data ocoul/.true./
c
c  Retrieve info from the database
c
       vlen = (nclosed+nact)*(nbf-nclosed-nact) + nclosed*nact
       if (.not. rtdb_get(rtdb,'mcscf:tol2e',MT_DBL,1,tol2e)) then
         if (.not. rtdb_get(rtdb,'scf:tol2e',MT_DBL,1,tol2e))
     $     tol2e = 1.d-10
       endif
       if (.not. rtdb_get(rtdb,'mcscf:lshift',MT_DBL,1,lshift0)) then
         if (.not. rtdb_get(rtdb,'scf:lshift',MT_DBL,1,lshift0))
     $     lshift0 = 0.1d0
       endif
       if (.not. rtdb_get(rtdb,'mcscf:hessiantype',MT_INT,1,pflg))
     $   pflg = 2
       if (.not. rtdb_get(rtdb,'mcscf:e2approx',MT_LOG,1,oe2approx))
     $   oe2approx = .false.
       if (.not. rtdb_get(rtdb,'mcscf:citol',MT_DBL,1,citol))
     $      citol = 1.d-8
       if (.not. rtdb_get(rtdb,'mcscf:line_search_tol',MT_DBL,
     $                    1,lsrch_tol)) lsrch_tol = 1.d-1
       if (.not.geom_nuc_rep_energy( geom, enrep ))
     $   call errquit('mcscf_orb: cannot retrieve nuclear repulsion',0,
     &       GEOM_ERR)
       molo = nclosed + 1
       mohi = nclosed + nact
       nsym = sym_number_ops(geom)+1
       if (.not.rtdb_get(rtdb,'mcscf:aoblock',MT_LOG,1,oblk)) then
         if (.not.rtdb_get(rtdb,'fourindex:aoblock',MT_LOG,1,oblk))
     $      oblk = .true.
       endif
       if (.not. rtdb_get(rtdb,'mcscf:obshift',mt_log, 1, obshift))
     $      obshift = .true.
c
c  Allocate MO matrices
c
*ga:1:0
       if (.not.ga_create(MT_DBL,nbf,nbf,'Act Fock',nbf,0,g_afock))
     $      call errquit('mcscf_orb: cannot allocate active Fock',0,
     &       GA_ERR)
*ga:1:0
       if (.not.ga_create(MT_DBL,nbf,nbf,'In Fock',nbf,0,g_ifock))
     $      call errquit('mcscf_orb: cannot allocate inactive Fock',0,
     &       GA_ERR)
*ga:1:0
       if (.not.ga_create(MT_DBL,nbf,nbf,'Gen Fock',nbf,0,g_gfock))
     $      call errquit('mcscf_orb: cannot allocate general Fock',0,
     &       GA_ERR)
       if (oe2approx) then
         nnact = (nact*(nact+1))/2
*ga:1:0
         if (.not.ga_create(MT_DBL,nbf,nbf,'Unitrary',nbf,0,g_u))
     $        call errquit('mcscf_orb: cannot allocate unitrary',0,
     &       GA_ERR)
*ga:1:0
         if (.not.ga_create(MT_DBL,nbf,nbf,'B matrix',nbf,0,g_b))
     $        call errquit('mcscf_orb: cannot allocate B',0,
     &       GA_ERR)
*ga:1:0
         if (.not.ga_create(MT_DBL,nbf,nbf,'Temp',nbf,0,g_tmp))
     $        call errquit('mcscf_orb: cannot allocate B',0,
     &       GA_ERR)
*ga:1:0
         if (.not.ga_create(MT_DBL,(nact*nact),nnact,'Act Ints',
     $                      (nact*nact), 0, g_actints ))
     $        call errquit('mcscf_orb: cannot allocate active ints',0,
     &       GA_ERR)
       endif
       if (.not.ga_duplicate(g_movecs,g_newvecs,'New MO vectors'))
     $      call errquit('mcscf_orb: cannot duplicate movecs',0,
     &       GA_ERR)
c
c  Allocate gradient vectors
c
*ga:1:0
       if (.not.ga_create(MT_DBL,max(1,vlen),1,'Gradient',0,0,g_grad))
     $      call errquit('mcscf_orbsolve: cannot allocate',0,
     &       GA_ERR)
*ga:1:0
       if (.not.ga_create(MT_DBL,max(1,vlen),1,'K vector',0,0,g_kvec))
     $      call errquit('mcscf_orbsolve: cannot allocate',0,
     &       GA_ERR)
*ga:1:0
       if (.not.ga_create(MT_DBL,max(1,vlen),1,'Search',0,0,g_srch))
     $      call errquit('mcscf_orbsolve: cannot allocate',0,
     &       GA_ERR)
*ga:1:0
       if (.not.ga_create(MT_DBL,max(1,vlen),1,'Hess Diag',0,0,g_hdiag))
     $      call errquit('mcscf_orbsolve: cannot allocate',0, GA_ERR)
       call ga_zero(g_kvec)
       call ga_zero(g_srch)
c
c  Compute Fock matrices and gradients
c
       call mcscf_fock( geom, basis, nbf, nclosed, nact,
     $                  oskel, tol2e, dm1, dm2,
     $                  g_movecs, g_coul, eone, etwo, ecore,
     $                  g_ifock, g_afock, g_gfock )
       call mcscf_gfock2grad( nbf, nclosed, nact, g_gfock, g_grad )
       call mcscf_vsymm_chk( nbf, nclosed, nact, orbsym,
     $                       g_grad, .true., .false., maxx )
       call ga_screen(g_grad,1.d-12)
       gnorm = sqrt(ga_ddot(g_grad,g_grad))
       gnorm0 = gnorm
       energy = enrep + eone + etwo
       if (oprconv) write(6,922) energy,gnorm
 922   format(11x,'>>>| MCSCF energy:',4x,f22.12,/,
     $        11x,'>>>| Gradient norm:',3x,d22.9)
       oconverged = (gnorm.lt.gnorm_tol)
       if (oconverged) goto 555                              ! jump to clean-up and return
       e0 = ecore + enrep
c
c  Precompute 1e-Hessian diagonal
c
       call mcscf_hdiag( nbf, nclosed, nact, pflg, dm1, dm2,
     $                   g_coul, g_exch, g_ifock,
     $                   g_afock, g_gfock, g_hdiag )
c
c  RHS = -Gradient
c
       call ga_dadd( -1.d0, g_grad, 0.d0, g_kvec, g_kvec )
c
c  Adjust level-shift
c
       if (gnorm .gt. 1.d-1) then
         lshift = lshift0
       else if (gnorm .gt. 1.d-2) then
         lshift = lshift0/5.d0
       else
         lshift = 0.d0
       endif
*       write(6,*) ' bshift ', bshift
       bshift = 0.d0
c
c  Check on vanishing/negative Hessian diagonals
c
       call mcscf_chk_hdiag( pflg, vlen, g_hdiag, lshift, hmin, hmax,
     $      oprconv)
c       CALL GA_ROWPRINT(' Scaled Hessian Diag ', G_HDIAG )
c
c  Precondition accuracy (same as rohf code)
c     rjh ... it's only first order convergent
c
       ndigit = 2.d0
       if (pflg .eq. 2) then
         ndigit = 1.5d0
*         if (lshift .eq. 0.0d0) ndigit =
*     $     min((1+abs(log10(gnorm/gnorm_tol))),abs(log10(gnorm)))
       endif
c
c  Precondition
c
       p_acc = max(min(10.d0**(-ndigit), tol_ceil), tol_floor)
       oprcndprint = util_print('precondition',print_high)
 1010  call mcscf_precond( geom, basis, nclosed, nact,
     $                     oskel, orbsym, tol2e, pflg,
     $                     lshift, maxpiter, oprcndprint,
     $                     dm1, dm2,
     $                     g_coul, g_exch, g_movecs,
     $                     g_ifock, g_afock, g_gfock,
     $                     g_hdiag, g_kvec, g_srch, p_acc )
c
       call mcscf_vsymm_chk( nbf, nclosed, nact, orbsym,
     $                       g_srch, .true., .false., maxx )
       if (maxx .gt. 1d-6) call errquit('mcscf_orb: large symmetry'//
     $      ' contamination in search direction',0, GEOM_ERR)
C       call ga_screen(g_srch, 1.d-10 )
       snorm = sqrt(ga_ddot(g_srch,g_srch))
c
c  Apply conjugacy if previous iteration did not force reset
c  and if gradient is not too large
c
       if ((gnorm0.lt.0.2d0).and.(oconj).and.(.not.(ocjrst))) then
         call mcscf_conj( g_grad, g_gradold, g_conj, g_srch )
       else if (ocjrst) then
         if (oprconv) write(6,723)
 723     format(5x,'--- Conjugacy reset ---')
       endif
       call mcscf_vsymm_chk( nbf, nclosed, nact, orbsym,
     $                       g_srch, .true., .false., maxx )
       if (maxx .gt. 1d-6) call errquit('mcscf_orb: large symmetry'//
     $      ' contamination in search direction after conjugacy',0,
     &       GEOM_ERR)
c
c  Info print
c
       olsrchprint = oprconv
       ociprint = util_print('microci',print_high)
       ocivecprint = util_print('civec',print_debug).and.ociprint
       if ((ga_nodeid().eq.0).and.
     $     (util_print('converge info',print_high)))
     $   write(6,966) hmin,hmax,(snorm/gnorm), lshift, bshift
 966   format(/,'Min/Max Hessian diags:',2x,2e10.1,
     $        /,'Scaled search norm:',5x,e12.2,
     $        /,'Level shift:',15x,f12.3,
     $        /,'E2 shift:',18x,f12.3)
c
c  Begin line search...
c
c$$$       CALL GA_ROWPRINT('step vector',G_KVEC)
*       if (oprconv) write(6,973)
* 973   format(/,10x,'step',8x,'energy',9x,'p.grad',6x,'ratio',
*     $          10x,'time',/,7x,63('-'))
       alpha = 1.0d0
       eev(0) = energy
       gsv(0) = ga_ddot(g_grad, g_srch)
       if (gsv(0) .gt. 0) then
c
c     Oops ... the line search is not down hill which means that
c     the preconditioning has failed, or there is an inconsistency
c     in the energy/gradient, or the conjugacy is messed up.
c     Ideally want to redo the preconditioning with a larger level shift
c     and also reset conjugacy.  Eventually we'll end up on the right track.
c     While we're at it we could also disable the seocnd order algorithm since
c     things are obviously messed up but this is not done now.
c
          lshift = 2.0d0*lshift + 0.5d0
          ocjrst = .true.
          if (ga_nodeid() .eq. 0) then
             write(6,'(a,f8.2)') 
     $            ' !! warning - +ve gradient along search '//
     $            'direction - Increasing shift ',  lshift
             call util_flush(6)
          endif
          goto 1010
       endif
c
       av(0) = 0.d0
*       micro_citol = gnorm0*0.1d0*lsrch_tol
       micro_citol = min(citol, (gnorm*0.05d0))
       eprec = max(gnorm0*0.01d0, tol2e*100.0d0,micro_citol*10.0d0)
       do iter=1,maxlinsrch
         iter1 = iter - 1
         nstep = iter
*         if (oprconv) write(6,974) av(iter1), eev(iter1), gsv(iter1),
*     $                                abs(gsv(iter1)/gsv(0)),
*     $                                util_cpusec()
* 974     format(5x,f10.4,2x,f15.10,2x,d10.2,2x,f10.4,2x,f10.1)
         if (line_search( nstep, maxlinsrch, eev, gsv, av, eprec,
     $                    lsrch_tol, olsrchprint )) goto 100
         alpha = av(nstep-1)
         call ga_dadd( alpha, g_srch, 0.d0, g_kvec, g_kvec )
         if (oe2approx) then
           call rohf_k2u( basis, nbf, nbf, nclosed, nact, g_kvec, g_u )
           call ga_dgemm('n', 'n', nbf, nbf, nbf, 1.d0, g_movecs,
     $                    g_u, 0.d0, g_newvecs )
         else
           call rohf_k2cf( basis, nbf, nbf, nclosed, nact, g_kvec,
     $                     g_movecs, g_newvecs )
         endif
c$$$         PRINT*,'New MO vectors'
c$$$         CALL GA_PRINT_X(G_NEWVECS)
c
c   Regenerate J and K 
c   Optionally relax CI
c
         if (oe2approx) then
           if (ocirelax) then
             odelta_core = (gnorm0.gt.0.1d0)                         ! toggle response of inactive orbitals in core
             odelta_core = .true.
             call mcscf_trfJKu( nbf, nclosed, nact, g_u, g_coul,
     $                          g_exch, g_actints )
             call mcscf_b_trf_ifock( geom, basis, nbf, nclosed, nact,
     $                               oskel, tol2e, g_movecs, g_ifock,
     $                               g_u, odelta_core, g_tmp )
             call mcscf_cisolve( rtdb, geom, basis, nbf, nclosed, nact,
     $                           nsym, orbsym, e0, orbe, g_tmp,
     $                           g_actints, micro_citol, .true.,
     $                           ociprint, ocivecprint, dm1, dm2, eci )
           endif
         else
           oexch = ocirelax
           if (omcscfprof) call pstat_on(ps_4index)
           call moints_build_2x( basis, ohalf, oskel,
     $                           molo, molo, mohi, 1, nbf,
     $                           g_newvecs, g_coul, ocoul,
     $                           g_exch, oexch, blen, oblk )
           if (omcscfprof) call pstat_off(ps_4index)
           if (ocirelax) then
             call mcscf_fcore( basis, nbf, nclosed, nact, g_newvecs,
     $                         g_coul, g_exch, g_hcore )
             call mcscf_cisolve( rtdb, geom, basis, nbf, nclosed, nact,
     $                           nsym, orbsym, e0, orbe, g_hcore,
     $                           g_coul, micro_citol, .true.,
     $                           ociprint, ocivecprint, dm1, dm2, eci )
           endif
         endif
c
c   New energy and gradient
c
         if (oe2approx) then
            bshift = 0.0d0
            if (obshift) bshift = 0.2d0 * lshift
           call mcscf_b( geom, basis, nbf, nclosed, nact, orbsym,
     $                   dm1, dm2, oskel, tol2e, bshift, g_movecs,
     $                   g_ifock, g_afock, g_coul, g_exch, g_u, g_b )
           call ga_dgemm( 't', 'n', nbf, nbf, nbf, 1.d0, g_b, g_u,
     $                    0.d0, g_tmp )
c
           call mcscf_gfock2grad(nbf, nclosed, nact, g_tmp, g_grad )
c
           do i=1+ga_nodeid(),nbf,ga_nnodes()
             call ga_acc(g_u, i, i, i, i, mone, 1, 1.d0 )
           enddo
c
c     Apply level shift of 0.5*bshift*T^2 to the energy, bshift*T to gradient
c     Don't have T in a convenient form but we have K.
c
**           write(6,*) ' BSHIFT ', bshift
           call ga_dadd( bshift, g_kvec, 1.d0, g_grad, g_grad )
           knorm = ga_ddot(g_kvec, g_kvec)
           call ga_transpose( g_b, g_tmp )
           call ga_dadd( 1.d0, g_gfock, 1.d0, g_tmp, g_b )
           call ga_transpose( g_b, g_tmp )
           energy = eev(0) + ga_ddot( g_tmp, g_u ) + 0.5d0*bshift*knorm
         else
           call mcscf_fock( geom, basis, nbf, nclosed, nact,
     $                      oskel, tol2e, dm1, dm2,
     $                      g_newvecs, g_coul, eone, etwo, ecore,
     $                      g_ifock, g_afock, g_gfock )
           e0 = ecore + enrep                               ! update core energy valid for previous macro
c$$$           PRINT*,' Generalized Fock matrix'
c$$$           CALL GA_PRINT_X(G_GFOCK)
           call mcscf_gfock2grad( nbf, nclosed, nact, g_gfock, g_grad )
           energy = enrep + eone + etwo
         endif
c$$$         CALL GA_ROWPRINT('Gradient',G_GRAD)
         av(iter) = alpha
         eev(iter) = energy
         gsv(iter) = ga_ddot(g_grad, g_srch)
         gsscale = abs(gsv(iter)/gsv(0))
       enddo
c
c  Line search exhausted...take last successful step
c  (alpha0) for new MO-vectors
c
       if (ga_nodeid().eq.0) write(6,927) alpha
 927   format(10x,'Line search exhausted --- step length is:',f10.2)
       call ga_dadd( alpha, g_srch, 0.d0, g_kvec, g_kvec )
       call rohf_k2cf( basis, nbf, nbf, nclosed, nact, g_kvec,
     $                 g_movecs, g_newvecs )
 100   continue
c
c  MO Info
c
       gnorm = sqrt(ga_ddot(g_grad, g_grad))
       if (util_print('new movecs',print_debug)) then
         if (ga_nodeid().eq.0) write(6,330) 
 330     format('New MO vectors')
         call ga_print(g_newvecs)
       endif
c  
c  Safe to copy new MOs into result
c  Canonicalization and symmetry adaption
c
       ocjrst = .false.
       if (ocanon) then
         call mcscf_canonical( nbf, nclosed, nact, g_ifock, g_afock,
     $                         g_gfock, orbe, g_newvecs, nflip )
c
c     Canonicalization can mess up the symmetry for degenerate irreps
c     that we have to split. So symm adapt here and ignore
c     any contamination.  Contamination elsewhere is still bad.
c
         call sym_movecs_adapt(basis, 1.d-8, g_newvecs, orbsym, nmixed)
         call ga_orthog_mos(basis, g_newvecs)
c
         call mcscf_gfock2grad( nbf, nclosed, nact, g_gfock, g_grad )
         if (oprconv) write(6,331)
 331     format(5x,'--- Orbitals canonicalized ---')
         if (util_print('canonical',print_high)) then
           if (ga_nodeid().eq.0) then
             write(6,332)
 332         format('Canonicalization - orbital energies')
             write(6,333) (orbe(i),i=1,nbf)
 333         format(7f12.6)
             write(6,334) nflip
 334         format(i5,' orbitals flipped')
             ocjrst = (nflip.gt.0)                   ! Reset conjugacy if flipped orbitals
           endif
         endif
       endif
c
c  Safe to copy MO vectors
c
       call ga_copy(g_newvecs, g_movecs )
       call sym_movecs_adapt( basis, 1.d-8, g_movecs, orbsym, nmixed )
*       write(6,*) ' MOVECS after sym adaption'
*       call ga_print(g_movecs)
       if (nmixed .ne. 0)
     $   call errquit('mcscf_orb: symmetry contamination', nmixed,
     &       GEOM_ERR)
c
c  Save updated MOs
c  Orbitals maybe reordered on save (movecs_lock)
c  so pass temporary MOs
c
       call ga_copy(g_movecs, g_newvecs) ! rjh ... save sym adapted vectors
       call mcscf_save_vec(rtdb, geom, basis, nbf, orbe, occ, g_newvecs)
c
c  Update Coulomb and Exchange operator
c     rjh ... must do it if have canonicalized even if exact energy LS
c
       if ((ocanon) .or. (oe2approx).or.(.not.(ocirelax))) then
         if (omcscfprof) call pstat_on(ps_4index)
         call moints_build_2x( basis, ohalf, oskel,
     $                         molo, molo, mohi, 1, nbf,
     $                         g_movecs, g_coul, .true.,
     $                         g_exch, .true., blen, oblk )
         if (omcscfprof) call pstat_off(ps_4index)
       endif
c
c  Clean-up
c
 555   continue
       if (oe2approx) then
         if (.not.ga_destroy(g_u))
     $        call errquit('mcscf_orbsolve: cannot destroy unitrary',0,
     &       GA_ERR)
         if (.not.ga_destroy(g_b))
     $        call errquit('mcscf_orbsolve: cannot destroy B',0,
     &       GA_ERR)
         if (.not.ga_destroy(g_tmp))
     $        call errquit('mcscf_orbsolve: cannot destroy B',0, GA_ERR)
         if (.not.ga_destroy(g_actints))
     $        call errquit('mcscf_orbsolve: cannot destroy B',0, GA_ERR)
       endif
       if (.not.ga_destroy(g_newvecs))
     $      call errquit('mcscf_orbsolve: cannot destroy MO matrix',0,
     &       GA_ERR)
       if (.not.ga_destroy(g_afock))
     $      call errquit('mcscf_orbsolve: cannot destroy act fock',0,
     &       GA_ERR)
       if (.not.ga_destroy(g_ifock))
     $      call errquit('mcscf_orbsolve: cannot destroy inact fock',0,
     &       GA_ERR)
       if (.not.ga_destroy(g_gfock))
     $      call errquit('mcscf_orbsolve: cannot destroy gen fock',0,
     &       GA_ERR)
       if (.not.ga_destroy(g_grad))
     $      call errquit('mcscf_orbsolve: cannot destroy gradient',0,
     &       GA_ERR)
       if (.not.ga_destroy(g_kvec))
     $      call errquit('mcscf_orbsolve: cannot destroy k-vector',0,
     &       GA_ERR)
       if (.not.ga_destroy(g_srch))
     $      call errquit('mcscf_orbsolve: cannot destroy search',0,
     &       GA_ERR)
       if (.not.ga_destroy(g_hdiag))
     $      call errquit('mcscf_orbsolve: cannot destroy search',0,
     &       GA_ERR)
       return
       end






c
c  Apply conjugacy to preconditioned vector
c
c
c          k+1 -1T  k+1   k          k+1  k+1   k
c         g  .A  .(g   - g )       -s  .(g   - g )            k+1
c  beta = ---------------------  =  --------------       NB. s   is preconditioned
c            k   k+1   k              k   k+1   k            search direction before
c           s .(g   - g )            s .(g   - g )           conjugacy is applied
c
c
c   k+1     k+1         k
c  s     = s    + beta.s
c
c

       subroutine mcscf_conj( g_grad, g_gradold, g_conj, g_srch )
       implicit none
#include "global.fh"
       integer g_grad, g_gradold, g_conj, g_srch
       double precision numer, denom, beta

       numer = ga_ddot(g_srch, g_grad) - ga_ddot(g_srch,g_gradold)
       denom = ga_ddot(g_conj, g_grad) - ga_ddot(g_conj,g_gradold)
       if (abs(denom).gt.1.d-12) then
         beta = -numer/denom
       else
         beta = 0.d0
       endif
       call ga_dadd( beta, g_conj, 1.d0, g_srch, g_srch )
       call ga_copy( g_srch, g_conj )
       call ga_copy( g_grad, g_gradold )

       return
       end






c
c  Approximate transformation of the active-space 
c  MO integrals, accurate for E2(T) surface.
c
c
       subroutine mcscf_trfJKu( nbf, nclosed, nact, g_u, g_coul,
     $                          g_exch, g_actints )
       implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
       integer nbf, nclosed, nact
       integer g_u, g_coul, g_exch
       integer g_actints
c
       integer l_s, k_s, l_z, k_z, k_x, l_u, k_u, ld
       integer clo, chi, aoff, aend, zoff
       integer t, u, tt, tu, tus, uts, v, w, vw, vwoff
       integer tv, uw, tvs, vts, uws, wus, wtop
       double precision qq, pp
       logical ga_check_JKblocked
       integer ijkl
       external ijkl
*       INTEGER IQ, VV
*       DOUBLE PRECISION XYZ,XXX(10000)
c
       zoff = nbf*nclosed + nclosed
       aoff = nclosed + 1
       aend = nclosed + nact
       if (.not.ma_push_get(MT_DBL,(nbf*nbf),'t',l_s,k_s))
     $   call errquit('mcscf_trfJKu: cannot allocate s',0, MA_ERR)
       if (.not.ma_push_get(MT_DBL,(nbf*nbf),'t',l_z,k_z))
     $   call errquit('mcscf_trfJKu: cannot allocate t',0, MA_ERR)
       if (.not.ma_push_get(MT_DBL,(nbf*nbf),'u',l_u,k_u))
     $   call errquit('mcscf_trfJKu: cannot allocate t',0, MA_ERR)
       if (.not.ga_check_JKblocked( g_coul, nact, nbf, clo, chi ))
     $   call errquit('mcscf_trfJKu: wrong distrib operator',0,
     &       GA_ERR)
       call ga_get(g_u,1,nbf,aoff,aend,dbl_mb(k_u),nbf)
       call ga_zero(g_actints)
c
c  
c         (2)     tu       tu           vw
c   (tu|vw)   = -J   + (U J  U)   + (U J  U)
c                 vw           vw           tu
c
       do t=1,nact
         do u=1,t
           tu = (t*(t-1))/2 + u
           tus = (t-1)*nact + u
           uts = (u-1)*nact + t
           qq = 1.d0
           if (t.eq.u) qq = 0.5d0
           if ((tu.ge.clo).and.(tu.le.chi)) then
             call ga_access(g_coul, 1, (nbf*nbf), tu, tu, k_x, ld )
c$$$             DO V=1,NACT
c$$$               CALL DCOPY(NACT,DBL_MB(K_X+(NCLOSED+V-1)*NBF+NCLOSED),1,
c$$$     $                         XXX((V-1)*NACT+1),1)
c$$$             ENDDO
c$$$             WRITE(6,'(10x,2I4)') T,U
c$$$             CALL MOINTS_MATPRINT(NACT,NACT,XXX)
             call dgemm( 't', 'n', nact, nbf, nbf, 1.d0, dbl_mb(k_u),
     $                   nbf, dbl_mb(k_x), nbf, 0.d0,
     $                   dbl_mb(k_z), nact )
             call dgemm( 'n', 'n', nact, nact, nbf, 1.d0, dbl_mb(k_z),
     $                   nact, dbl_mb(k_u), nbf, 0.d0,
     $                   dbl_mb(k_s), nact )
             do v=1,nact
               do w=1,v
                 vw = (v*(v-1))/2 + w
                 vwoff = (v-1)*nact + w - 1
                 call ga_acc(g_actints, tus, tus, vw, vw,
     $                       dbl_mb(k_s+vwoff), 1, qq )
                 call ga_acc(g_actints, uts, uts, vw, vw,
     $                       dbl_mb(k_s+vwoff), 1, qq )
               enddo
             enddo
             do v=1,nact
               call daxpy( nact, -1.d0, dbl_mb(k_x+zoff+(v-1)*nbf), 1,
     $                     dbl_mb(k_s+(v-1)*nact), 1 )
             enddo
             call ga_acc( g_actints, 1, (nact*nact), tu, tu,
     $                    dbl_mb(k_s), (nact*nact), 1.d0 )
             call ga_release(g_coul, 1, (nbf*nbf), tu, tu )
           endif
         enddo
       enddo
c$$$       CALL GA_ZERO(G_ACTINTS)
c
c         (2)        tu
c   (tv|uw)   += (T K  T)
c                        vw
c
c$$$       CALL DFILL(1000,0.d0,XXX,1)

       do t=1,nact
         tt = (t-1)*nbf + nclosed + t - 1
         dbl_mb(k_u+tt) = dbl_mb(k_u+tt) - 1.d0
       enddo
       do t=1,nact
         do u=1,t
           tu = (t*(t-1))/2 + u
           tus = (t-1)*nact + u
           uts = (u-1)*nact + t
           if ((tu.ge.clo).and.(tu.le.chi)) then
             call ga_access(g_exch, 1, (nbf*nbf), tu, tu, k_x, ld )
             call dgemm( 't', 'n', nact, nbf, nbf, 1.d0, dbl_mb(k_u),
     $                   nbf, dbl_mb(k_x), nbf, 0.d0,
     $                   dbl_mb(k_z), nact )
             call dgemm( 'n', 'n', nact, nact, nbf, 1.d0, dbl_mb(k_z),
     $                   nact, dbl_mb(k_u), nbf, 0.d0,
     $                   dbl_mb(k_s), nact )
             call ga_release(g_exch, 1, (nbf*nbf), tu, tu )
c$$$             WRITE(6,'(10x,2I4)') T,U
c$$$             CALL MOINTS_MATPRINT(NACT,NACT,DBL_MB(K_S))
             do v=1,nact
               tv = (max(t,v)*(max(t,v)-1))/2 + min(t,v)
               tvs = (t-1)*nact + v
               vts = (v-1)*nact + t
               pp = 1.d0
               if (t.eq.v) pp = 2.d0
               wtop = nact
               if (t.eq.u) wtop = v
               do w=1,wtop
                 uw = (max(u,w)*(max(u,w)-1))/2 + min(u,w)
                 uws = (u-1)*nact + w
                 wus = (w-1)*nact + u
                 vwoff = (v-1)*nact + w - 1
c
c  Logic is messy but irreducible (?)
c
                 qq = pp
                 if (u.eq.w) qq = qq*2.d0
                 if ((t.eq.w).and.(u.eq.v).and.(t.ne.u)) qq = qq*2.d0
                 call ga_acc(g_actints, tvs, tvs, uw, uw,
     $                       dbl_mb(k_s+vwoff), 1, qq )
                 if (t.ne.v)
     $             call ga_acc(g_actints, vts, vts, uw, uw,
     $                         dbl_mb(k_s+vwoff), 1, qq )

                 if (uw.ne.tv) then
                   call ga_acc(g_actints, uws, uws, tv, tv,
     $                         dbl_mb(k_s+vwoff), 1, qq )
                   if (u.ne.w)
     $                call ga_acc(g_actints, wus, wus, tv, tv,
     $                            dbl_mb(k_s+vwoff), 1, qq )
                 endif
c
c  Debug
c
c$$$                 IQ = IJKL(T,V,U,W)
c$$$                 XYZ = DBL_MB(K_S+VWOFF)
c$$$                 XXX(IQ) = XXX(IQ) + DBL_MB(K_S+VWOFF)*QQ
               enddo
             enddo
           endif
         enddo
       enddo
c
c
c        
c$$$       PRINT*, ' FINAL ACTIVE INTS '
c$$$       DO T=1,NACT
c$$$         DO U=1,T
c$$$           WRITE(6,'(10x,2I4)') T,U
c$$$           TU = (T*(T-1))/2 + U
c$$$           CALL GA_GET(G_ACTINTS, 1, (NACT*NACT), TU, TU,
c$$$     $                 DBL_MB(K_S), 1 )
c$$$           CALL MOINTS_MATPRINT(NACT,NACT,DBL_MB(K_S))
c$$$         ENDDO
c$$$       ENDDO
c
c
c$$$
c$$$       DO T=1,NACT
c$$$         DO U=1,T
c$$$           DO V=1,T
c$$$             TU = V
c$$$             IF (V.EQ.T) TU=U
c$$$             DO W=1,TU
c$$$               IQ = IJKL(T,U,V,W)
c$$$               WRITE(6,'(4I3,3X,I5,3X,F20.8)') T,U,V,W,IQ,XXX(IQ)
c$$$             ENDDO
c$$$           ENDDO
c$$$         ENDDO
c$$$       ENDDO

       if (.not.ma_pop_stack(l_u))
     $   call errquit('mcscf_trfJKu: failed pop stack',0, MA_ERR)
       if (.not.ma_pop_stack(l_z))
     $   call errquit('mcscf_trfJKu: failed pop stack',0, MA_ERR)
       if (.not.ma_pop_stack(l_s))
     $   call errquit('mcscf_trfJKu: failed pop stack',0, MA_ERR)
       return
       end








      integer function ijkl( i, j, k, l )
      implicit none
      integer i, j, k, l
      integer ij, kl
      
      ij = (max(i,j)*(max(i,j)-1))/2 + min(i,j)
      kl = (max(k,l)*(max(k,l)-1))/2 + min(k,l)
      ijkl = (max(ij,kl)*(max(ij,kl)-1))/2 + min(ij,kl)
      return
      end





      subroutine mcscf_pr_j( nbf, nclosed, nact, g_coul )
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      integer nbf, nclosed, nact
      integer g_coul
      double precision xxx(1000)
      integer clo, chi
      integer t,u,tu,k_x,ld,v
      logical ga_check_JKblocked

      if (.not.ga_check_JKblocked( g_coul, nact, nbf, clo, chi ))
     $   call errquit('mcscf_trfJKu: wrong distrib operator',0, GA_ERR)
      do t=1,nact
        do u=1,t
          tu = (t*(t-1))/2 + u
          if ((tu.ge.clo).and.(tu.le.chi)) then
            write(6,'(10X,2I4)') t,u
            call ga_access(g_coul, 1, (nbf*nbf), tu, tu, k_x, ld )
            do v=1,nact
              call dcopy(nact,dbl_mb(k_x+(nclosed+v-1)*nbf+nclosed),1,
     $                        xxx((v-1)*nact+1),1)
            enddo
            call ga_release(g_coul, 1, (nbf*nbf), tu, tu )
            call moints_matprint( nact, nact, xxx )
          endif
        enddo
      enddo
      return
      end

             





c
c  Ensure approximate Hessian diagonals are
c  positive definite
c

       subroutine mcscf_chk_hdiag( pflg, vlen, g_hdiag, lshift,
     $                             hmin, hmax, oprconv )
       implicit none
#include "global.fh"
#include "msgids.fh"
       integer pflg
       integer vlen
       integer g_hdiag
       double precision lshift
       double precision hmin, hmax
       
       integer i, ismall
       double precision yy

       logical oprconv          ! [input] if true then print

       if (pflg.gt.1) then
         hmax = -10.d0
         hmin = 10.d0
         ismall = 0
         do i=ga_nodeid()+1,vlen,ga_nnodes()
           call ga_get(g_hdiag,i,i,1,1,yy,1)
           yy = yy + lshift
           hmax = max(hmax,yy)
           hmin = min(hmin,yy)
           if (abs(yy).lt.0.05d0) ismall = ismall + 1
         enddo
         call ga_sync()
         call ga_dgop(Msg_HessMin,hmin,1,'min')
         call ga_dgop(Msg_HessMax,hmax,1,'max')
         if (hmin.lt.0.1d0) then
           do i=ga_nodeid()+1,vlen,ga_nnodes()
             call ga_get(g_hdiag,i,i,1,1,yy,1)
             yy = yy + lshift - hmin + 0.1d0
             call ga_put(g_hdiag,i,i,1,1,yy,1)
           enddo
         endif
         call ga_sync()
         call ga_igop(msg_mcscf_chkhess,ismall,1,'+')
         if ((ga_nodeid().eq.0).and.(ismall.gt.0) .and. oprconv) then
           write(6,112) ismall
 112       format(5x,'--- Warning:',i5,' Hessian diagonals < 0.05 ---')
         endif
       else
         do i=ga_nodeid()+1,vlen,ga_nnodes()
           call ga_get(g_hdiag,i,i,1,1,yy,1)
           yy = yy + lshift
           call ga_put(g_hdiag,i,i,1,1,yy,1)
         enddo
       endif
       return
       end
