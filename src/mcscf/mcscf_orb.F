c
c  Solve orbital MCSCF equations using PCG/NR
c
c  $Id: mcscf_orb.F,v 1.3 1995-12-15 11:59:15 d3g681 Exp $
c
       subroutine mcscf_orbsolve( rtdb, geom, basis, nbf, nclosed,
     $                           nact, osym, dm1, dm2,
     $                           g_coul, g_exch, oconj, pfac,
     $                           g_conj, g_gradold,
     $                           eone, etwo, g_movecs )
       implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "sym.fh"
c
       integer rtdb, geom, basis                       ! [input] Handles
       integer nbf, nclosed, nact                      ! [input] Orbital parameters
       logical osym                                    ! [input] Symmetry toggle
       double precision dm1(*)                         ! [input] 1-PDM
       double precision dm2(*)                         ! [input] 2-PDM
       integer g_coul                                  ! [input] Coulomb operator
       integer g_exch                                  ! [input] Exchange operator
       logical oconj                                   ! [input] Conjugacy toggle
       double precision pfac                           ! [input/output] Conj projection factor
       integer g_conj                                  ! [input] Conjugacy vector
       integer g_gradold                               ! [input] Old gradient
       double precision eone                           ! [output] 1e-energy + nuclear repulsion
       double precision etwo                           ! [output] 2e-energy
       integer g_movecs                                ! [input/output] MO vectors
c
       integer maxpiter
       integer maxlinsrch
       parameter(maxpiter=15)
       parameter(maxlinsrch=16)
       integer iter, ngrad, maxiter, pflg
       integer molo, mohi
       double precision tol2e, enrep, energy
       double precision gnorm, lshift, lshift0, p_acc
       double precision alpha, gsscale
       double precision eev(0:maxlinsrch), gsv(0:maxlinsrch)
       double precision av(0:maxlinsrch)
       double precision xx, yy, eold
       logical oskel
       integer g_afock, g_ifock, g_gfock
       integer g_newvecs, g_grad, g_kvec, g_srch, g_hdiag
c
       double precision linesrch
       external linesrch
       double precision linesrch1
       external linesrch1
       data pflg/2/
c
c
c
       ngrad = (nclosed+nact)*(nbf-nclosed-nact) + nclosed*nact
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
       if (.not.ga_create(MT_DBL,nGrad,1,'Gradient',1,1,g_grad))
     $      call errquit('mcscf_orbsolve: cannot allocate',0)
       if (.not.ga_create(MT_DBL,nGrad,1,'K vector',1,1,g_kvec))
     $      call errquit('mcscf_orbsolve: cannot allocate',0)
       if (.not.ga_create(MT_DBL,nGrad,1,'Search',1,1,g_srch))
     $      call errquit('mcscf_orbsolve: cannot allocate',0)
       if (.not.ga_create(MT_DBL,nGrad,1,'Hess Diag',1,1,g_hdiag))
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
       gnorm = sqrt(ga_ddot(g_grad,g_grad))
       eone = eone + enrep
       energy = eone + etwo
       if (ga_nodeid().eq.0) write(6,933) energy,gnorm
 933   format(10x,'ENERGY:',f20.10,5x, 'Gradient Norm:',e12.3)
c
c  Precompute 1e-Hessian diagonal
c
       call mcscf_hdiag( nbf, nclosed, nact, dm1, dm2,
     $                   g_coul, g_exch, g_ifock,
     $                   g_afock, g_gfock, g_hdiag )
C       CALL GA_ROWPRINT(' HESSIAN DIAG ', g_hdiag )
c
c  RHS = -G 
c
       call ga_dadd( -1.d0, g_grad, 0.d0, g_kvec, g_kvec )
       lshift = lshift0
       if (gnorm.le.1.d-2) lshift = 0.d0
c
c  Precondition
c
       p_acc = 1.d-3
       call mcscf_precond( geom, basis, nclosed, nact,
     $                     oskel, tol2e, pflg, lshift,
     $                     maxpiter, dm1, dm2,
     $                     g_coul, g_exch, g_movecs,
     $                     g_ifock, g_afock, g_gfock,
     $                     g_hdiag, g_kvec, g_srch, p_acc )
C       CALL GA_ROWPRINT(' SEARCH DIRECTION ', g_srch )
c
c  Conjugacy
c
       if (oconj) then
         call mcscf_conj( g_grad, g_gradold, pfac, g_conj, g_srch )
         IF (GA_NODEID().EQ.0)
     $     WRITE(6,'(A,5X,E10.2)') 'Applied conjugacy',pfac
       endif
c
c
c  Line search...
c    **Approximate** because J integrals 
c    are not regenerated
c
       alpha = 0.8d0
       eold = energy
       eev(0) = energy
       gsv(0) = ga_ddot(g_grad, g_srch)
       av(0) = 0.d0
       if (ga_nodeid().eq.0) write(6,966) lshift
 966   format(/,'Level shift:',f12.2,
     $        /,15x,'Alpha',15x,'Energy',12x,'GProj',4x,
     $         'GProjScaled',/,14x,60('-'))
       do iter=1,maxlinsrch
         call ga_dadd( alpha, g_srch, 0.d0, g_kvec, g_kvec )
         call rohf_k2cf( basis, nbf, nclosed, nact, g_kvec,
     $                   g_movecs, g_newvecs )
         if (gnorm.gt.1.d-3) then
           call moints_build_6x( basis, osym,
     $                           molo, molo, mohi, 1, nbf,
     $                           g_newvecs, g_coul, .true.,
     $                           g_exch, .false., 16 )
         endif
         call mcscf_fock( geom, basis, nbf, nclosed, nact,
     $                    oskel, tol2e, dm1, dm2,
     $                    g_newvecs, g_coul, eone, etwo,
     $                    g_ifock, g_afock, g_gfock )
         call mcscf_gfock2grad( nbf, nclosed, nact, g_gfock, g_grad )
         av(iter) = alpha
         eev(iter) = eone + enrep + etwo
         gsv(iter) = ga_ddot(g_grad, g_srch)
         gsscale = abs(gsv(iter)/gsv(0))
         if (ga_nodeid().eq.0)
     $     write(6,957) '---',av(iter),eev(iter),gsv(iter),gsscale
 957     format(5x,a3,2x,f12.5,5x,f20.12,2e12.2)
         if (gsscale.lt.1.d-3) goto 100
         alpha = linesrch1( (iter+1), av, eev, gsv )
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
C$Id: mcscf_orb.F,v 1.3 1995-12-15 11:59:15 d3g681 Exp $
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
