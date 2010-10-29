       subroutine mcscf_precond( geom, basis, nclosed, nact,
     $                           oskel, orbsym, tol2e, pflg,
     $                           lshift, maxiter, oprcndprint,
     $                           dm1, dm2,
     $                           g_coul, g_exch, g_movecs,
     $                           g_ifock, g_afock, g_gfock,
     $                           g_hdiag, g_rhs, g_soln, acc )
*
* $Id$
*
       implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "sym.fh"
#include "mcscfprof.fh"
c
c
       integer geom, basis                              ! [input] Handles
       integer nclosed, nact                            ! [input] Parameters
       logical oskel                                    ! [input] Symmetry toggle
       integer orbsym(*)                                ! [input] Orbital irreps
       double precision tol2e                           ! [input] Hessian accuracy
       integer pflg                                     ! [input] Hessian-type selection
       double precision lshift                          ! [input] Level-shift
       integer maxiter                                  ! [input] Maximum iterations
       logical oprcndprint                              ! [input] Print toggle
       double precision dm1(*), dm2(*)                  ! [input] 1- & 2-PDM                  
       integer g_coul, g_exch                           ! [input] Coulomb & exchange operators
       integer g_movecs                                 ! [input] MO vectors
       integer g_ifock, g_afock, g_gfock                ! [input] Inactive, active & gen. Fock matrices
       integer g_hdiag                                  ! [input] Hessian diagonal
       integer g_rhs                                    ! [input] Negative gradient (RHS)
       integer g_soln                                   ! [output] Solution
       double precision acc                             ! [input/output] Solution accuracy
c
c
c
       integer maxbas
       parameter(maxbas=7)
       integer nbf, nvir, vlen
       integer g_basis, g_sigma
       integer g_xv, g_axv
       integer nbas, ibas, jbas, iter
       double precision pconvg
       double precision aa(maxbas,maxbas), bb(maxbas)
       double precision basis_norm(maxbas)
       integer itemp(maxbas)
       double precision rnorm, bnorm, proj, scale, deltasoln
       double precision rnorm_old, ratio
       double precision maxx
c$$$       DOUBLE PRECISION XX
       integer info
c
c
c
       if (omcscfprof) call pstat_on(ps_prcnd)
       if (.not. bas_numbf( basis, nbf))
     $    call errquit('mcscf_precond: bas_numbf failed',0,
     &       BASIS_ERR)
       nvir = nbf - nclosed - nact
       vlen = nvir*(nclosed+nact) + nclosed*nact
c
c
c
*ga:1:0
       if (.not.ga_create(MT_DBL,vlen,maxbas,'basis',0,maxbas,g_basis))
     $    call errquit('mcscf_precond: failed to allocate',0, GA_ERR)
*ga:1:0
       if (.not.ga_create(MT_DBL,vlen,maxbas,'sigma',0,maxbas,g_sigma))
     $    call errquit('mcscf_precond: failed to allocate',0, GA_ERR)
*ga:1:0
       if (.not.ga_create(MT_DBL,vlen,1,'x vector',0,0,g_xv))
     $    call errquit('mcscf_precond: failed to allocate',0, GA_ERR)
*ga:1:0
       if (.not.ga_create(MT_DBL,vlen,1,'prod vector',0,0,g_axv))
     $    call errquit('mcscf_precond: failed to allocate',0, GA_ERR)
1010       nbas = 0
       if ((oprcndprint).and.(ga_nodeid().eq.0)) write(6,937) 
 937   format(/,20x,'Iter',7x,'RNorm',8x,'Accuracy',5x,'DeltaSoln',2x,
     $        'Convg Rate',/,20x,59('-'))
c
c
c
       rnorm_old = 1.d0
       call ga_zero(g_soln)
       bnorm = sqrt(ga_ddot(g_rhs,g_rhs))
       call ga_copy(g_rhs,g_xv)
       call mcscf_hdiag_scale( nbf, nclosed, nact, lshift, dm1,
     $                         g_ifock, g_afock, g_gfock,
     $                         g_hdiag, g_xv )
c$$$       PRINT*,'LSHIFT: ',lshift
c$$$       CALL GA_ROWPRINT( 'RHS ', G_RHS )
c$$$       CALL GA_ROWPRINT( 'GUESS SOLUTION', G_XV )
C       GOTO 120
c
c  Start iterative solver
c
       do iter=1,maxiter
         nbas = nbas + 1
         basis_norm(nbas) = sqrt(ga_ddot(g_xv,g_xv))
         call mcscf_hessv( geom, basis, nbf, nclosed, nact,
     $                     tol2e, oskel, pflg, lshift, dm1, dm2,
     $                     g_movecs, g_ifock, g_afock, g_gfock,
     $                     g_coul, g_exch, g_xv, g_axv )
         call mcscf_vsymm_chk( nbf, nclosed, nact, orbsym,
     $                         g_axv, .true., .false., maxx )
c
         if (maxx .gt. 1d-6) call errquit ! rjh ... always expect the worst
     $        ('mcscf_prcnd: large symmetry contamination', 0,
     &       GEOM_ERR)
c
C         CALL GA_ROWPRINT( 'HESSIAN VECTOR PRODUCT', G_AXV )

         call ga_copy_patch('n',g_axv,1,vlen,1,1,
     $                       g_sigma,1,vlen,nbas,nbas)
         call ga_copy_patch('n',g_xv,1,vlen,1,1,
     $                       g_basis,1,vlen,nbas,nbas)

         do ibas=1,nbas
           do jbas=1,ibas
             aa(jbas,ibas) = ga_ddot_patch(
     $                              g_sigma, 'n', 1, vlen, ibas, ibas,
     $                              g_sigma, 'n', 1, vlen, jbas, jbas)
             aa(ibas,jbas) = aa(jbas,ibas)
           enddo
           bb(ibas)= ga_ddot_patch( g_sigma, 'n', 1, vlen, ibas, ibas,
     $                              g_rhs, 'n', 1, vlen, 1, 1 )
         enddo
c$$$         PRINT*,'SUBSPACE MATRIX'
c$$$         CALL MAT_PRINT2(AA,1,NBAS,1,NBAS,MAXBAS,MAXBAS)
c$$$         PRINT*,'SUBSPACE RHS'
c$$$         CALL MAT_PRINT2(BB,1,NBAS,1,1,MAXBAS,1)
c
c
         call dgesv(nbas, 1, aa, maxbas, itemp, bb, maxbas, info)
c
c
c$$$         PRINT*,'INFO:',info
c$$$         PRINT*,'SUBSPACE SOLUTION'
c$$$         CALL MAT_PRINT2(BB,1,NBAS,1,1,MAXBAS,1,1)
c
c Residue & Current solution vector
c
         call ga_copy(g_rhs,g_xv)
         call ga_copy(g_soln,g_axv)
         call ga_zero(g_soln)
         do ibas=1,nbas
           call ga_dadd_patch( -bb(ibas), g_sigma, 1, vlen, ibas, ibas,
     $                         1.d0,      g_xv, 1, vlen, 1, 1,
     $                                    g_xv, 1, vlen, 1, 1 )
           call ga_dadd_patch( bb(ibas), g_basis, 1, vlen, ibas, ibas,
     $                         1.d0,     g_soln, 1, vlen, 1, 1,
     $                                   g_soln, 1, vlen, 1, 1 )
         enddo
         call ga_dadd(1.d0,g_axv,-1.d0,g_soln,g_axv)
         deltasoln = ga_ddot(g_axv,g_axv)
         rnorm = sqrt(ga_ddot(g_xv,g_xv))
         pconvg = rnorm/bnorm
         ratio = rnorm/rnorm_old
         rnorm_old = rnorm
         if ((oprcndprint).and.(ga_nodeid().eq.0))
     $     write(6,955) iter,rnorm,pconvg,deltasoln,ratio
 955     format(20x,i3,2x,3e14.3,f10.4)
         if (pconvg.lt.acc) goto 119
c
c  Form new basis vector from residue
c
         call mcscf_hdiag_scale( nbf, nclosed, nact, lshift, dm1,
     $                           g_ifock, g_afock, g_gfock,
     $                           g_hdiag, g_xv )
c
c  Gram-Schmidt
c
         call ga_copy(g_xv,g_axv)
         do ibas=1,nbas
           proj = ga_ddot_patch( g_axv, 'n', 1, vlen, 1, 1,
     $                           g_basis, 'n', 1, vlen, ibas, ibas )
           scale = -proj/basis_norm(ibas)
           call ga_dadd_patch( 1.d0,  g_xv, 1, vlen, 1, 1,
     $                         scale, g_basis, 1, vlen, ibas, ibas,
     $                                g_xv, 1, vlen, 1, 1)
         enddo
c
c  Contract basis if necessary
c
         if (nbas.eq.maxbas) then
           call ga_copy_patch('n' ,g_soln, 1, vlen, 1, 1,
     $                             g_basis, 1, vlen, 1, 1 )
           call ga_zero(g_axv)
           do ibas=1,nbas
             call ga_dadd_patch( bb(ibas), g_sigma, 1, vlen,ibas,ibas,
     $                           1.d0,  g_axv, 1, vlen, 1, 1,
     $                                  g_axv, 1, vlen, 1, 1 )
           enddo
           call ga_copy_patch('n' ,g_axv, 1, vlen, 1, 1,
     $                             g_sigma, 1, vlen, 1, 1 )
           nbas = 1
         endif
c
c$$$         PROJ = GA_DDOT_PATCH( G_AXV, 'N', 1, VLEN, 1, 1,
c$$$     $                         G_AXV, 'N', 1, VLEN, 1, 1 )
c$$$         XX = GA_DDOT_PATCH( G_XV, 'N', 1, VLEN, 1, 1,
c$$$     $                       G_XV, 'N', 1, VLEN, 1, 1 )
c$$$         WRITE(6,702) SQRT(PROJ),SQRT(XX)
c$$$ 702     FORMAT('ORIG NORM:',5X,E15.4,'NEW NORM:',5X,E15.4)
c$$$         PROJ = GA_DDOT_PATCH( G_XV, 'N', 1, VLEN, 1, 1,
c$$$     $                         G_AXV, 'N', 1, VLEN, 1, 1 )
c$$$         WRITE(6,701) PROJ
c$$$ 701     FORMAT('OVERLAP WITH ORIG:',5X,E15.4)
c$$$         DO IBAS=1,NBAS
c$$$           PROJ = GA_DDOT_PATCH( G_XV, 'N', 1, VLEN, 1, 1,
c$$$     $                           G_BASIS, 'N', 1, VLEN, IBAS, IBAS )
c$$$           WRITE(6,700) IBAS,PROJ,BASIS_NORM(IBAS)
c$$$ 700       FORMAT('ORTHOG:',I5,5X,2E15.4)
c$$$         ENDDO
c$$$         XX = 1.D0/SQRT(GA_DDOT(G_XV,G_XV))
c$$$         CALL GA_DSCAL(G_XV,XX)
C         CALL GA_ROWPRINT( 'NEW BASIS', G_XV )
       enddo
c
c Maximum iterations exhausted
c
       if (ga_nodeid().eq.0) write(6,961) pconvg, acc
 961   format(5x,'--- Warning: Precondition failed to converge:',
     $        'Residue: current=',d8.1,' required=',d8.1,'---')

       if (abs(pconvg/acc) .lt. 2.0d0) goto 120  ! Close enuf

c
c      Increase the level shift and try again.
c      Failure is not acceptable since otherwise the line search
c      will get lost.  C.f. Sean Connery ... "Only losers whine
c      about their best; winners go home and **** the prom queen."

      lshift = lshift*2.0d0 + 0.2

      if (ga_nodeid() .eq. 0) then
        write(6,'(a,f8.2)') ' Increase level shift to ', lshift
        call util_flush(6)
      endif
      goto 1010
c
c
c
 119   continue
C$$$ CHECK SOLUTION ACCURACY
c$$$       WRITE(6,*)
c$$$       CALL MCSCF_HESSV( GEOM, BASIS, NBF, NCLOSED, NACT,
c$$$     $                   OSKEL, DM1, DM2, G_MOVECS,
c$$$     $                   G_IFOCK, G_AFOCK, G_GFOCK,
c$$$     $                   G_COUL, G_EXCH, G_SOLN, G_AXV )
c$$$       CALL GA_DADD(1.D0,G_RHS,-1.D0,G_AXV,G_AXV)
c$$$       PCONVG = SQRT(GA_DDOT(G_AXV,G_AXV))
c$$$       WRITE(6,981) PCONVG
c$$$ 981   FORMAT(/,30X,'SOLUTION ACCURACY:',E12.3)
       if ((oprcndprint).and.(ga_nodeid().eq.0)) write(6,771)
 771   format(5x,'--- converged ---')
c
c
c
 120   continue
       acc = pconvg
C       CALL GA_ROWPRINT( 'SOLUTION', G_AXV )
       if (.not.ga_destroy(g_basis))
     $    call errquit('mcscf_precond: failed to destroy',0, GA_ERR)
       if (.not.ga_destroy(g_sigma))
     $    call errquit('mcscf_precond: failed to destroy',0, GA_ERR)
       if (.not.ga_destroy(g_xv))
     $    call errquit('mcscf_precond: failed to destroy',0, GA_ERR)
       if (.not.ga_destroy(g_axv))
     $    call errquit('mcscf_precond: failed to destroy',0, GA_ERR)
       if (omcscfprof) call pstat_off(ps_prcnd)

       return
       end
