c
c  Davidson diagonalization with 
c  3-term recursion
c
      subroutine detci_dvddiag( n, mxv, precon, maxiter, tol,
     $                          g_x, g_x0, g_q, g_s, g_s0, g_sq, e )
      implicit none
#include "global.fh"
#include "util.fh"
      integer n                               ! [input] dimension
      external mxv                            ! [input] matrix-vector product routine
      external precon                         ! [input] preconditioner routine
      integer maxiter                         ! [input] maximum allowed iterations
      double precision tol                    ! [input] convergence tolerance
      integer g_x                             ! [input/output] starting guess/soln vector
      integer g_x0                            ! [scratch] previous iteration solution
      integer g_q                             ! [scratch] new basis vector
      integer g_s                             ! [scratch] Vector product of soln
      integer g_s0                            ! [scratch] Vector product of prev. soln
      integer g_sq                            ! [scratch] Vector product of new basis
      double precision e                      ! [output] Solution eigenvalue

c$$$      double precision x(n)                   ! [input/output] starting guess/soln vector
c$$$      double precision x0(n)                  ! [scratch] previous iteration solution
c$$$      double precision q(n)                   ! [scratch] new basis vector
c$$$      double precision s(n)                   ! [scratch] Vector product of soln
c$$$      double precision s0(n)                  ! [scratch] Vector product of prev. soln
c$$$      double precision sq(n)                  ! [scratch] Vector product of new basis
c
c
c
      integer nscr, n3
      parameter(nscr=100)
      parameter(n3=3)
      double precision cc(n3,n3)               ! Subspace solution coeff.
      double precision ov(n3,n3)               ! Subspace overlap matrix 
      double precision ee(n3)                  ! Subspace eigenvalues
      double precision scr(nscr)
      double precision rnorm, rnormscal, xx, eold
      integer info, iter, nz
      REAL TX, TY
c
c
      double precision dnrm2
      external dnrm2
c
c
      TY = util_cpusec()
      WRITE(6,906)
 906  FORMAT(/,'Iter',8x,'Energy',10x,'DeltaE',
     $       6x,'Rnorm',7x,'Coeff',7x,'Time')
      CALL FLUSH(6)
c
c
c
      eold = 0.d0
      iter = 0
      call ga_copy( g_x, g_q )
      call ga_zero( g_x )
      call ga_zero( g_x0 )
      call ga_zero( g_s ) 
      call ga_zero( g_s0 )
      call ga_zero( g_sq )
      call dfill((n3*n3),0.d0,cc,1)
      call dfill((n3*n3),0.d0,ov,1)
c$$$      call dcopy(n,x,1,q,1)
c$$$      call dfill(n,0.d0,x,1)
c$$$      call dfill(n,0.d0,x0,1)
c$$$      call dfill(n,0.d0,s,1)
c$$$      call dfill(n,0.d0,s0,1)
c      
c
c
 100  continue
      iter = iter + 1
      nz = min(iter,3)
c
c  Matrix-vector product
c
      call ga_zero( g_sq )
c$$$      call dfill(n,0.d0,sq,1)
      call mxv(g_q,g_sq)
c$$$      CALL DETCI_PRINT( sq, 1.D-2 )
c
c  Subspace terms
c      
      call dfill((n3*n3),0.d0,cc,1)
      call dfill((n3*n3),0.d0,ov,1)
      cc(1,1) = ga_ddot(g_sq,g_q)
      cc(2,1) = ga_ddot(g_sq,g_x)
      cc(3,1) = ga_ddot(g_sq,g_x0)
      cc(2,2) = ga_ddot(g_s,g_x)
      cc(3,2) = ga_ddot(g_s,g_x0)
      cc(3,3) = ga_ddot(g_s0,g_x0)
      ov(1,1) = ga_ddot(g_q,g_q)
      ov(2,1) = ga_ddot(g_q,g_x)
      ov(3,1) = ga_ddot(g_q,g_x0)
      ov(2,2) = ga_ddot(g_x,g_x)
      ov(3,2) = ga_ddot(g_x,g_x0)
      ov(3,3) = ga_ddot(g_x0,g_x0)

c$$$      cc(1,1) = ddot(n,sq,1,q,1)
c$$$      cc(2,1) = ddot(n,sq,1,x,1)
c$$$      cc(3,1) = ddot(n,sq,1,x0,1)
c$$$      cc(2,2) = ddot(n,s,1,x,1)
c$$$      cc(3,2) = ddot(n,s,1,x0,1)
c$$$      cc(3,3) = ddot(n,s0,1,x0,1)
c$$$      ov(1,1) = ddot(n,q,1,q,1)
c$$$      ov(2,1) = ddot(n,q,1,x,1)
c$$$      ov(3,1) = ddot(n,q,1,x0,1)
c$$$      ov(2,2) = ddot(n,x,1,x,1)
c$$$      ov(3,2) = ddot(n,x,1,x0,1)
c$$$      ov(3,3) = ddot(n,x0,1,x0,1)
C
c$$$      WRITE(6,'(/,A,/)') 'HAMILTONIAN & OVERLAP'
c$$$      CALL MAT_PRINT(N3,N3,CC)
c$$$      WRITE(6,*)
c$$$      CALL MAT_PRINT(N3,N3,OV)
c$$$      WRITE(6,*)
c$$$      IF (ITER.GE.10) RETURN
c
c  Generalized eigenvalue
c
      call dsygv(1,'V','L',nz,cc,n3,ov,n3,ee,scr,nscr,info)
c
c$$$      WRITE(6,'(/,A,5x,I5)') 'INFO',info
c$$$      WRITE(6,'(/,A)') 'EIGENVECTORS'
c$$$      CALL MAT_PRINT(N3,N3,CC)
c$$$      WRITE(6,'(/,A)') 'EIGENVALUES'
c$$$      CALL MAT_PRINT(N3,1,EE)
c
c
c                  k+1               k         k    k-1
c  New solution:  x   = c .q  +  c .x  +  c .(x  - x   )
c                        1        2        3
c
      call ga_dadd( cc(1,1), g_q,  cc(2,1), g_x, g_q )
      call ga_dadd( cc(3,1), g_x0, 1.d0,    g_q, g_q )
      call ga_copy( g_q, g_x0 )
      call ga_dadd( -1.d0, g_x, 1.d0, g_x0, g_x0 )
      call ga_copy( g_q, g_x )

c$$$      call dscal(n,cc(1,1),q,1)
c$$$      call daxpy(n,cc(2,1),x,1,q,1)
c$$$      call daxpy(n,cc(3,1),x0,1,q,1)
c$$$      call dcopy(n,q,1,x0,1)
c$$$      call daxpy(n,-1.d0,x,1,x0,1)
c$$$      call dcopy(n,q,1,x,1)
c
c
c  New sigma vector (analogous to above)
c
      call ga_dadd( cc(1,1), g_sq, cc(2,1), g_s,  g_sq )
      call ga_dadd( cc(3,1), g_s0, 1.d0,    g_sq, g_sq )
      call ga_copy( g_sq, g_s0 )
      call ga_dadd( -1.d0, g_s, 1.d0, g_s0, g_s0 )
      call ga_copy( g_sq, g_s )

c$$$      call dscal(n,cc(1,1),sq,1)
c$$$      call daxpy(n,cc(2,1),s,1,sq,1)
c$$$      call daxpy(n,cc(3,1),s0,1,sq,1)
c$$$      call dcopy(n,sq,1,s0,1)
c$$$      call daxpy(n,-1.d0,s,1,s0,1)
c$$$      call dcopy(n,sq,1,s,1)
c
c  Residue  q = s - e.x
c

      call ga_dadd( -ee(1), g_x, 1.d0, g_s, g_q )
      rnorm = sqrt(ga_ddot( g_q, g_q ))

c$$$      call dcopy(n,s,1,q,1)
c$$$      call daxpy(n,-ee(1),x,1,q,1)
c$$$      rnorm = dnrm2(n,q,1)

      rnormscal = rnorm/sqrt(abs(ee(1)))
      TX = util_cpusec() - TY
      TY = util_cpusec()
      if (ga_nodeid().eq.0) then
        write(6,881) iter,ee(1),(ee(1)-eold),rnormscal,abs(cc(1,1)),tx
 881    format(i3,f20.10,3e12.2,f10.2)
        call flush(6)
      endif
      eold = ee(1)
      if ((rnormscal.lt.tol).or.(iter.ge.maxiter)) goto 101
c
c  Preconditioning with normalization
c      

      call precon( g_q, g_sq )
      xx = sqrt(ga_ddot( g_q, g_q ))
      call ga_dscal( g_q, (1.d0/xx) )

c$$$      call precon(q,sq)
c$$$      xx = dnrm2(n,q,1)
c$$$      call dscal(n,(1.d0/xx),q,1)
c
c  Normalize old vector?
c
c$$$      xx = dnrm2(n,x0,1)
c$$$      call dscal(n,(1.d0/xx),x0,1)
c$$$      call dscal(n,(1.d0/xx),s0,1)
c
c
c
      goto 100
c
c  End main loop
c
 101  continue
      xx = sqrt( ga_ddot(g_x,g_x) )
      call ga_dscal( g_x, (1.d0/xx) )
      e = ee(1)
      return
      end

      









c
c
c Simple single root Davidson diag
c
      subroutine detci_dvddiag0( mxv, precon, nw, tol, eval, g_w, g_x)
      implicit none
#include "global.fh"
#include "util.fh"
      external mxv                            ! [input] matrix-vector product routine
      external precon                         ! [input] preconditioner routine
      integer nw                              ! [input] vectors in workspace
      double precision tol                    ! [input] convergence tolerance
      double precision eval                   ! [output] solution eigenvalue
      integer g_w(nw)                         ! [scratch] workspace
      integer g_x                             ! [input/output] starting guess/soln. vector
c$$$      double precision w(n,nw)                ! [scratch] workspace
c$$$      double precision x(n)                   ! [input/output] starting guess/soln. vector
c
      integer mxsubsp
      integer nz
      integer maxiter
      parameter(mxsubsp=20)
      parameter(nz=mxsubsp*32)
      parameter(maxiter=12)
      double precision e(mxsubsp)
      double precision b(mxsubsp*mxsubsp)
      double precision z(nz)
      integer mxbas, ibas, info, iter
      double precision rnorm, rnormscal, xx, eold
      integer i,j,ij
c
      double precision dnrm2
      external dnrm2
c
c  Initialize
c
      iter = 0
      eold = 0.d0
      mxbas = min((nw/2),mxsubsp)
      ibas = 1
      call ga_copy( g_x, g_w(1) )
c$$$      call dcopy(n,x,1,w(1,1),1)
c
c  Main iterative loop
c
 101  continue
      iter = iter + 1
      call ga_zero( g_w(mxbas+ibas) )
      call mxv( g_w(ibas), g_w(mxbas+ibas) )
c$$$      call mxv( w(1,ibas), w(1,mxbas+ibas) )
      call dfill((ibas*ibas),0.d0,b,1)
      do i=1,ibas
        do j=1,i
          ij = (j-1)*ibas + i
          b(ij) = ga_ddot( g_w(i), g_w(mxbas+j) )
c$$$          b(ij) = ddot(n,w(1,i),1,w(1,mxbas+j),1)
        enddo
      enddo
      call dsyev( 'V', 'L', ibas, b, ibas, e, z, nz, info )
c$$$      PRINT*,'INFO:',INFO
c
c  Residue
c
      call ga_zero(g_x)
c$$$      call dfill(n,0.d0,x,1)
      do i=1,ibas
        call ga_dadd( b(i), g_w(mxbas+i), 1.d0, g_x, g_x )
        call ga_dadd( (-b(i)*e(1)), g_w(i), 1.d0, g_x, g_x )
c$$$        call daxpy(n,b(i),w(1,mxbas+i),1,x,1)
c$$$        call daxpy(n,(-b(i)*e(1)),w(1,i),1,x,1)
      enddo
c$$$      IF (GA_NODEID().EQ.0) PRINT*,'<<< Residue >>>'
c$$$      CALL DETCI_PRINT(G_X,1.D-3)
c
c
c
      rnorm = sqrt(ga_ddot(g_x,g_x))
c$$$      rnorm = dnrm2(n,x,1)
      rnormscal = rnorm/sqrt(abs(e(1)))
      if (ga_nodeid().eq.0) then
        write(6,881) iter,e(1),(e(1)-eold),rnorm,rnormscal,abs(b(ibas))
 881    format(i3,f20.10,4e12.2)
      endif
      if (rnormscal.lt.tol) goto 199
c
c  Contract basis (if necessary)
c
      if (ibas.eq.mxbas) goto 199
c$$$      if (ibas.eq.mxbas) then
c$$$        call dscal(n,b(1),w(1,1),1)
c$$$        call dscal(n,b(1),w(1,mxbas+1),1)
c$$$        do i=2,ibas
c$$$          call daxpy(n,b(i),w(1,i),1,w(1,1),1)
c$$$          call daxpy(n,b(i),w(1,mxbas+i),1,w(1,mxbas+1),1)
c$$$        enddo
c$$$        ibas = 1
c$$$      endif
c
c  Preconditioner
c
      ibas = ibas + 1
      call precon( g_x, g_w(ibas) )
c$$$      IF (GA_NODEID().EQ.0) PRINT*,'<<< Preconditioned Residue >>>'
c$$$      CALL DETCI_PRINT(G_X,1.D-3)
c$$$      call precon(x)
c
c  Add new basis vector and orthonormalize
c
      call ga_copy( g_x, g_w(ibas) )
      do i=1,ibas-1
        xx = ga_ddot( g_x, g_w(i) )
        call ga_dadd( -xx, g_w(i), 1.d0, g_w(ibas), g_w(ibas) )
      enddo
      xx = sqrt( ga_ddot( g_w(ibas), g_w(ibas) ))
      call ga_dscal( g_w(ibas), (1.d0/xx) )
      call ga_zero( g_w(mxbas+ibas) )
c$$$      IF (GA_NODEID().EQ.0) THEN
c$$$        WRITE(6,'(A,5X,E10.2)') 'Normalization:',XX
c$$$        PRINT*,'<<< New basis >>>'
c$$$      ENDIF
c$$$      CALL DETCI_PRINT( G_W(IBAS), 1.D-6 )
c$$$      call dcopy(n,x,1,w(1,ibas),1)
c$$$      do i=1,ibas-1
c$$$        xx = ddot(n,x,1,w(1,i),1)
c$$$        call daxpy(n,-xx,w(1,i),1,w(1,ibas),1)
c$$$      enddo
c$$$      xx = dnrm2(n,w(1,ibas),1)
c$$$      call dscal(n,(1.d0/xx),w(1,ibas),1)
c$$$      call dfill(n,0.d0,w(1,mxbas+ibas),1)
      eold = e(1)
      if (iter.lt.maxiter) goto 101
      ibas = ibas - 1
c
c  Contract basis into solution vector
c
 199  continue
      call ga_zero(g_x)
      do i=1,ibas
        call ga_dadd( b(i), g_w(i), 1.d0, g_x, g_x )
      enddo
c$$$      call dfill(n,0.d0,x,1)
c$$$      do i=1,ibas
c$$$        call daxpy(n,b(i),w(1,i),1,x,1)
c$$$      enddo
      eval = e(1)
c
c  Check:  E = x.H.x
c
c$$$      CALL DFILL(N,0.D0,W(1,1),1)
c$$$      CALL MXV( X, W(1,1) )
c$$$      XX = DDOT(N,X,1,W(1,1),1)
c$$$      WRITE(6,'(10x,A,5X,F20.10)')  '2e energy (expection value)',XX
c
c
c

      return
      end











