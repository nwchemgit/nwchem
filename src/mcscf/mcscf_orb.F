c
c  Solve orbital MCSCF equations using PCG/NR
c
c  $Id: mcscf_orb.F,v 1.5 1996-02-13 23:51:18 rg240 Exp $
c
       subroutine mcscf_orbsolve( rtdb, geom, basis, nbf, nclosed,
     $                           nact, osym, orbsym, dm1, dm2,
     $                           g_coul, g_exch, oconj, pfac,
     $                           g_conj, g_gradold, g_hcore, orbe,
     $                           eone, etwo, gnorm, g_movecs )
       implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "sym.fh"
#include "msgids.fh"
c
       integer rtdb, geom, basis                       ! [input] Handles
       integer nbf, nclosed, nact                      ! [input] Orbital parameters
       logical osym                                    ! [input] Symmetry toggle
       integer orbsym                                  ! [input] Orbital irreps
       double precision dm1(*)                         ! [input] 1-PDM
       double precision dm2(*)                         ! [input] 2-PDM
       integer g_coul                                  ! [input] Coulomb operator
       integer g_exch                                  ! [input] Exchange operator
       logical oconj                                   ! [input] Conjugacy toggle
       double precision pfac                           ! [input/output] Conj projection factor
       integer g_conj                                  ! [input] Conjugacy vector
       integer g_gradold                               ! [input] Old gradient
       integer g_hcore                                 ! [input] 1e hamiltonian
       double precision orbe(*)                        ! [input] ROHF orbital energies
       double precision eone                           ! [output] 1e-energy + nuclear repulsion
       double precision etwo                           ! [output] 2e-energy
       double precision gnorm                          ! [output] gradient norm
       integer g_movecs                                ! [input/output] MO vectors
c
       integer maxpiter
       integer maxlinsrch
       parameter(maxpiter=15)
       parameter(maxlinsrch=16)
       integer iter, vlen, maxiter, pflg
       integer molo, mohi, i, ismall, nstep
       double precision tol2e, enrep, energy
       double precision snorm, lshift, lshift0, p_acc
       double precision alpha, gsscale, gnorm0
       double precision eev(0:maxlinsrch), gsv(0:maxlinsrch)
       double precision av(0:maxlinsrch)
       double precision xx, yy, eold, hmax, hmin, maxx
       logical oskel
       integer g_afock, g_ifock, g_gfock
       integer g_newvecs, g_grad, g_kvec, g_srch, g_hdiag
c
       double precision linesrch
       external linesrch
       logical line_search
       external line_search
       double precision linesrch1
       external linesrch1
       data pflg/2/
c
c
c
       vlen = (nclosed+nact)*(nbf-nclosed-nact) + nclosed*nact
       if (.not. rtdb_get(rtdb,'mcscf:tol2e',MT_DBL,1,tol2e)) then
         if (.not. rtdb_get(rtdb,'scf:tol2e',MT_DBL,1,tol2e))
     $     tol2e = 1.0d-12
       endif
       if (.not. rtdb_get(rtdb,'mcscf:lshift',MT_DBL,1,lshift0)) then
         if (.not. rtdb_get(rtdb,'scf:lshift',MT_DBL,1,lshift0))
     $     lshift0 = 0.0d0
       endif
       if (.not.rtdb_get(rtdb,'scf:skeleton',MT_LOG,1,oskel)) then
         oskel = sym_number_ops(geom).gt.0
       endif
       if (.not.geom_nuc_rep_energy( geom, enrep ))
     $   call errquit('mcscf: cannot retrieve nuclear repulsion',0)
       molo = nclosed + 1
       mohi = nclosed + nact
c
c  Allocate MO matrices
c
       if (.not.ga_create(MT_DBL,nbf,nbf,'Act Fock',nbf,1,g_afock))
     $      call errquit('mcscf: cannot allocate active Fock',0)
       if (.not.ga_create(MT_DBL,nbf,nbf,'In Fock',nbf,1,g_ifock))
     $      call errquit('mcscf: cannot allocate inactive Fock',0)
       if (.not.ga_create(MT_DBL,nbf,nbf,'Gen Fock',nbf,1,g_gfock))
     $      call errquit('mcscf: cannot allocate general Fock',0)
       if (.not.ga_duplicate(g_movecs,g_newvecs,'New MO vectors'))
     $      call errquit('mcscf: cannot duplicate movecs',0)
c
c  Allocate gradient vectors
c
       if (.not.ga_create(MT_DBL,vlen,1,'Gradient',1,1,g_grad))
     $      call errquit('mcscf_orbsolve: cannot allocate',0)
       if (.not.ga_create(MT_DBL,vlen,1,'K vector',1,1,g_kvec))
     $      call errquit('mcscf_orbsolve: cannot allocate',0)
       if (.not.ga_create(MT_DBL,vlen,1,'Search',1,1,g_srch))
     $      call errquit('mcscf_orbsolve: cannot allocate',0)
       if (.not.ga_create(MT_DBL,vlen,1,'Hess Diag',1,1,g_hdiag))
     $      call errquit('mcscf_orbsolve: cannot allocate',0)
       call ga_zero(g_kvec)
       call ga_zero(g_srch)
c
c  Compute Fock matrices and gradients
c
       call mcscf_fock( geom, basis, nbf, nclosed, nact,
     $                  oskel, tol2e, dm1, dm2,
     $                  g_movecs, g_coul, eone, etwo,
     $                  g_ifock, g_afock, g_gfock )
       call mcscf_gfock2grad( nbf, nclosed, nact, g_gfock, g_grad )
       call mcscf_vsymm_chk( nbf, nclosed, nact, orbsym,
     $                       g_grad, .true., .false., maxx )
       call ga_screen(g_grad,1.d-10)
       gnorm = sqrt(ga_ddot(g_grad,g_grad))
       gnorm0 = gnorm
       eone = eone + enrep
       energy = eone + etwo
c
c  Precompute 1e-Hessian diagonal
c
       call mcscf_hdiag( nbf, nclosed, nact, dm1, dm2,
     $                   g_coul, g_exch, g_ifock,
     $                   g_afock, g_gfock, g_hdiag )
c
c  RHS = -Gradient
c
       call ga_dadd( -1.d0, g_grad, 0.d0, g_kvec, g_kvec )
c
c  Adjust level-shift
c
       if (gnorm.gt.1.d-2) then
         lshift = lshift0
       else if (lshift.eq.lshift0) then
         lshift = lshift/5.d0
       else
         lshift = 0.d0
       endif
c
c  Check on vanishing/negative Hessian diagonals
c
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
       if (hmin.lt.0.1d0) then
         do i=ga_nodeid()+1,vlen,ga_nnodes()
           call ga_get(g_hdiag,i,i,1,1,yy,1)
           yy = yy + lshift - hmin + 0.1d0
           call ga_put(g_hdiag,i,i,1,1,yy,1)
         enddo
       endif
       call ga_sync()
       call ga_igop(Msg_ChkMO,ismall,1,'+')
       if ((ga_nodeid().eq.0).and.(ismall.gt.0)) then
         write(6,112) ismall
 112     format('Warning:',i5,' Hessian diagonals < 0.05d0')
       endif
       if (ga_nodeid().eq.0) write(6,113) hmax,hmin
 113   format('Hessian spectrum:',2e12.2)
c$$$       CALL GA_ROWPRINT(' HESSIAN DIAG ', g_hdiag )
c
c  Precondition
c
       p_acc = 1.d-1
       call mcscf_precond( geom, basis, nclosed, nact,
     $                     oskel, orbsym, tol2e, pflg,
     $                     lshift, maxpiter, dm1, dm2,
     $                     g_coul, g_exch, g_movecs,
     $                     g_ifock, g_afock, g_gfock,
     $                     g_hdiag, g_kvec, g_srch, p_acc )
c$$$       CALL GA_ROWPRINT(' RHS ', g_kvec )
c$$$       CALL GA_ROWPRINT(' SEARCH DIRECTION ', g_srch )
       call mcscf_vsymm_chk( nbf, nclosed, nact, orbsym,
     $                       g_srch, .true., .false., maxx )
       call ga_screen(g_srch, 1.d-10 )
       snorm = sqrt(ga_ddot(g_srch,g_srch))
c
c  Conjugacy
c
       if (oconj)
     $   call mcscf_conj( g_grad, g_gradold, pfac, g_conj, g_srch )
c
c  Info
c
       if (ga_nodeid().eq.0) write(6,966) (snorm/gnorm), lshift
 966   format(/,'Scaled search norm:',5x,e12.2,
     $        /,'Level shift:',15x,f12.3)
c$$$       IF (GA_NODEID().EQ.0) WRITE(6,967)
c$$$ 967   FORMAT(/,15X,'ALPHA',15X,'ENERGY',12X,'GPROJ',4X,
c$$$     $              'GPROJSCALED',/,14X,60('-'))
c
c  Line search...
c
       alpha = 1.0d0
       eold = energy
       eev(0) = energy
       gsv(0) = ga_ddot(g_grad, g_srch)
       av(0) = 0.d0
       do iter=1,maxlinsrch
         nstep = iter
         if (line_search( nstep, maxlinsrch, eev, gsv, av, 1.d-8,
     $       2.d-1, .true. )) goto 100
         alpha = av(nstep-1)
         call ga_dadd( alpha, g_srch, 0.d0, g_kvec, g_kvec )
         call rohf_k2cf( basis, nbf, nbf, nclosed, nact, g_kvec,
     $                   g_movecs, g_newvecs )
c
c   Resolve secular eqn for each step?
c   Have to regenerate J and K for STEPCI
c   only J required otherwise
c
#ifdef STEPCI
         if (gnorm0.gt.1.d-2) then
           write(6,771) 
 771       format(20x,'CI resolved')
           call moints_build_6x( basis, osym,
     $                           molo, molo, mohi, 1, nbf,
     $                           g_newvecs, g_coul, .true.,
     $                           g_exch, .true., 16 )
           call mcscf_fcore( basis, nbf, nclosed, nact, g_newvecs,
     $                       g_coul, g_exch, g_hcore )
           call mcscf_cisolve( rtdb, geom, basis, nbf, nclosed, nact,
     $                         orbe, g_hcore, g_coul, dm1, dm2,
     $                         .false., .false. )
        else
         call moints_build_6x( basis, osym,
     $                         molo, molo, mohi, 1, nbf,
     $                         g_newvecs, g_coul, .true.,
     $                         g_exch, .false., 16 )
        endif
#else
         call moints_build_6x( basis, osym,
     $                         molo, molo, mohi, 1, nbf,
     $                         g_newvecs, g_coul, .true.,
     $                         g_exch, .false., 16 )
#endif
c
c
         call mcscf_fock( geom, basis, nbf, nclosed, nact,
     $                    oskel, tol2e, dm1, dm2,
     $                    g_newvecs, g_coul, eone, etwo,
     $                    g_ifock, g_afock, g_gfock )
         call mcscf_gfock2grad( nbf, nclosed, nact, g_gfock, g_grad )
         av(iter) = alpha
         eev(iter) = eone + enrep + etwo
         gsv(iter) = ga_ddot(g_grad, g_srch)
         gsscale = abs(gsv(iter)/gsv(0))
c$$$         if (ga_nodeid().eq.0)
c$$$     $     write(6,957) '---',av(iter),eev(iter),gsv(iter),gsscale
c$$$ 957     format(5x,a3,2x,f12.5,5x,f20.12,2e12.2)
c$$$         if (gsscale.lt.2.d-1) goto 100
c$$$         alpha = linesrch1( (iter+1), av, eev, gsv )
       enddo
c
c  line search exhausted...take last successful step
c  (alpha0) for new MO-vectors
c
       if (ga_nodeid().eq.0) write(6,927) alpha
 927   format(10x,'Line search exhausted --- step length is:',f10.2)
       call ga_dadd( alpha, g_srch, 0.d0, g_kvec, g_kvec )
       call rohf_k2cf( basis, nbf, nclosed, nact, g_kvec,
     $                 g_movecs, g_newvecs )
c  
c  Safe to copy new MOs into result
c
 100   continue
       call ga_copy(g_newvecs, g_movecs )
       gnorm = sqrt(ga_ddot(g_grad, g_grad))
c
c  Clean-up
c
       if (.not.ga_destroy(g_newvecs))
     $      call errquit('mcscf_orbsolve: cannot destroy MO matrix',0)
       if (.not.ga_destroy(g_afock))
     $      call errquit('mcscf_orbsolve: cannot destroy MO matrix',0)
       if (.not.ga_destroy(g_ifock))
     $      call errquit('mcscf_orbsolve: cannot destroy MO matrix',0)
       if (.not.ga_destroy(g_gfock))
     $      call errquit('mcscf_orbsolve: cannot destroy MO matrix',0)
       if (.not.ga_destroy(g_grad))
     $      call errquit('mcscf_orbsolve: cannot destroy gradient',0)
       if (.not.ga_destroy(g_kvec))
     $      call errquit('mcscf_orbsolve: cannot destroy k-vector',0)
       if (.not.ga_destroy(g_srch))
     $      call errquit('mcscf_orbsolve: cannot destroy search',0)
       if (.not.ga_destroy(g_hdiag))
     $      call errquit('mcscf_orbsolve: cannot destroy search',0)
       return
       end






       subroutine lsrch_proj( g_grad, g_srch, gs0, alpha1, alpha0 )
C$Id: mcscf_orb.F,v 1.5 1996-02-13 23:51:18 rg240 Exp $
       implicit none
#include "global.fh"
       integer g_grad, g_srch
       double precision gs0, alpha1, alpha0
       double precision gs1, atmp

       gs1 = ga_ddot(g_grad,g_srch)
       atmp = (gs1*alpha0 - gs0*alpha1)/(gs1 - gs0)
       alpha0 = alpha1
       alpha1 = atmp
       gs0 = gs1
       return
       end








       subroutine mcscf_conj( g_grad, g_gradold, pfac, g_conj, g_srch )
       implicit none
#include "global.fh"
       integer g_grad, g_gradold, g_conj, g_srch
       double precision pfac
       double precision pfac1, beta

       pfac1 = ga_ddot(g_grad,g_srch)
       beta = (pfac1 - ga_ddot(g_gradold,g_srch))/pfac
       call ga_dadd(beta,g_conj,1.d0,g_srch,g_srch)
       call ga_copy(g_srch,g_conj)
       call ga_copy(g_grad,g_gradold)
       pfac = pfac1
       return
       end
