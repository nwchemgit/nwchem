      subroutine shake_chain2(n,nb,
     >                       tol,maxit,
     >                       dsq,mass,
     >                       r2,r1)
      implicit none
      integer n,nb
      real*8  tol
      integer maxit
      real*8  dsq(*),mass(*)
      real*8  r2(3,*), r1(3,*)
      
c        REAL*8       RXI(n), RYI(n), RZI(n)
c        REAL*8       PXI(n), PYI(n), PZI(n)
        LOGICAL     MOVING(n)
        LOGICAL     MOVED(n)

        LOGICAL     DONE
        INTEGER     IT, A, B, I,J
        REAL*8      PXAB, PYAB, PZAB, PABSQ
        REAL*8      RXAB, RYAB, RZAB, RABSQ, DIFFSQ, RPAB
        REAL*8      GAB, DX, DY, DZ, TOL2
        REAL*8      RPTOL, RMA, RMB
        PARAMETER ( RPTOL = 1.0E-6 )
        
        TOL2   = 2.0 * tol

        do A = 1, n
c              RXI(A) = r1(1,indx(A))
c              RYI(A) = r1(2,indx(A))
c              RZI(A) = r1(3,indx(A))
c
c              PXI(A) = r2(1,indx(A))
c              PYI(A) = r2(2,indx(A))
c              PZI(A) = r2(3,indx(A))

              MOVING(A) = .FALSE.
              MOVED(A)  = .TRUE.
        end do

           IT = 0
           DONE = .FALSE.

C       ** BEGIN ITERATIVE LOOP **

1000       IF ( ( .NOT. DONE ) .AND. ( IT .LE. MAXIT ) ) THEN

              DONE = .TRUE.

              DO 300 A = 1, nb

                 B = A + 1
                 IF ( B .GT. n ) B = 1

                 IF ( MOVED(A) .OR. MOVED(B) ) THEN

                    PXAB = r2(1,A) - r2(1,B)  !PXI(A) - PXI(B)
                    PYAB = r2(2,A) - r2(2,B)  !PYI(A) - PYI(B)
                    PZAB = r2(3,A) - r2(3,B)  !PZI(A) - PZI(B)
                    call lattice_min_difference(PXAB,PYAB,PZAB)

                    PABSQ  = PXAB ** 2 + PYAB ** 2 + PZAB ** 2
                    RABSQ  = DSQ(A)
                    DIFFSQ = RABSQ - PABSQ


                    IF ( ABS(DIFFSQ) .GT. ( RABSQ * TOL2 ) ) THEN

                       RXAB = r1(1,A) - r1(1,B) !RXI(A) - RXI(B)
                       RYAB = r1(2,A) - r1(2,B) !RYI(A) - RYI(B)
                       RZAB = r1(3,A) - r1(3,B) !RZI(A) - RZI(B)
                       call lattice_min_difference(RXAB,RYAB,RZAB)

                       RPAB = RXAB * PXAB + RYAB * PYAB + RZAB * PZAB

                       IF ( RPAB .LT. ( RABSQ * RPTOL ) ) THEN

                          STOP 'CONSTRAINT FAILURE '

                       ENDIF

                       RMA = 1.0 / mass(A)
                       RMB = 1.0 / mass(B)
                       GAB = DIFFSQ / ( 2.0 * ( RMA + RMB ) * RPAB )
                       DX  = RXAB * GAB
                       DY  = RYAB * GAB
                       DZ  = RZAB * GAB

                       r2(1,A) = r2(1,A) + RMA*DX !PXI(A) = PXI(A) + RMA * DX
                       r2(2,A) = r2(2,A) + RMA*DY !PYI(A) = PYI(A) + RMA * DY
                       r2(3,A) = r2(3,A) + RMA*DZ !PZI(A) = PZI(A) + RMA * DZ

                       r2(1,B) = r2(1,B) - RMB*DX !PXI(B) = PXI(B) - RMB * DX
                       r2(2,B) = r2(2,B) - RMB*DY !PYI(B) = PYI(B) - RMB * DY
                       r2(3,B) = r2(3,B) - RMB*DZ !PZI(B) = PZI(B) - RMB * DZ

                       MOVING(A) = .TRUE.
                       MOVING(B) = .TRUE.
                       DONE = .FALSE.

                    ENDIF

                 ENDIF

300           CONTINUE

              DO 400 A = 1, n
                 MOVED(A) = MOVING(A)
                 MOVING(A) = .FALSE.
400           CONTINUE

              IT = IT + 1
              GOTO 1000

           ENDIF

C       ** END ITERATIVE LOOP **

           IF ( .NOT. DONE ) THEN

              WRITE(*,'('' TOO MANY CONSTRAINT ITERATIONS '')')
              STOP

           ENDIF

c           do a = 1, n
c              r2(1,indx(a)) = PXI(a)
c              r2(2,indx(a)) = PYI(a)
c              r2(3,indx(a)) = PZI(a)
c              r1(1,indx(a)) = RXI(a)
c              r1(2,indx(a)) = RYI(a)
c              r1(3,indx(a)) = RZI(a)
c           end do

        RETURN
        END
c $Id$
