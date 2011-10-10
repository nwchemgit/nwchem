      subroutine shake_chain(n,indx,nb,
     >                       tol,maxit,
     >                       dsq,mass,
     >                       r2,r1)
      implicit none
      integer n,indx(*),nb
      real*8  tol
      integer maxit
      real*8  dsq(*),mass(*)
      real*8  r2(3,*), r1(3,*)
      



        REAL*8       RXI(n), RYI(n), RZI(n)
        REAL*8       PXI(n), PYI(n), PZI(n)
        LOGICAL     MOVING(n)
        LOGICAL     MOVED(n)

        LOGICAL     DONE
        INTEGER     IT, A, B, I
        REAL*8      PXAB, PYAB, PZAB, PABSQ
        REAL*8      RXAB, RYAB, RZAB, RABSQ, DIFFSQ, RPAB
        REAL*8      GAB, DX, DY, DZ, TOL2
        REAL*8      RPTOL, RMA, RMB
        PARAMETER ( RPTOL = 1.0E-6 )

        real*8 da(3),ua(3,3)
        real*8 c1,c2,c3

*     **** external functions ****
      real*8   lattice_unita
      external lattice_unita


*     ***** Determine the unit lattice vectors and distances ******
      do i=1,3
         da(i)   = lattice_unita(1,i)**2
     >           + lattice_unita(2,i)**2
     >           + lattice_unita(3,i)**2
         da(i)   = da(i)**(0.5)
         ua(1,i) = lattice_unita(1,i)/da(i)
         ua(2,i) = lattice_unita(2,i)/da(i)
         ua(3,i) = lattice_unita(3,i)/da(i)
      end do



        
        TOL2   = 2.0 * TOL



        do a = 1, n
              RXI(a) = r1(1,indx(a))
              RYI(a) = r1(2,indx(a))
              RZI(a) = r1(3,indx(a))
              PXI(a) = r2(1,indx(a))
              PYI(a) = r2(2,indx(a))
              PZI(a) = r2(3,indx(a))

              MOVING(A) = .FALSE.
              MOVED(A)  = .TRUE.
        end do

           IT = 0
           DONE = .FALSE.

C       ** BEGIN ITERATIVE LOOP **

1000       IF ( ( .NOT. DONE ) .AND. ( IT .LE. MAXIT ) ) THEN

              DONE = .TRUE.

              DO 300 A = 1, NB

                 B = A + 1
                 IF ( B .GT. N ) B = 1

                 IF ( MOVED(A) .OR. MOVED(B) ) THEN

                    PXAB = PXI(A) - PXI(B)
                    PYAB = PYI(A) - PYI(B)
                    PZAB = PZI(A) - PZI(B)
                    c1 = PXAB*ua(1,1) + PYAB*ua(2,1) + PZAB*ua(3,1)
                    c2 = PXAB*ua(1,2) + PYAB*ua(2,2) + PZAB*ua(3,2)
                    c3 = PXAB*ua(1,3) + PYAB*ua(2,3) + PZAB*ua(3,3)
                    c1 = c1 - ANINT(c1/da(1))*da(1)
                    c2 = c2 - ANINT(c2/da(2))*da(2)
                    c3 = c3 - ANINT(c3/da(3))*da(3)
                    PXAB = ua(1,1)*c1 + ua(1,2)*c2 + ua(1,3)*c3
                    PYAB = ua(2,1)*c1 + ua(2,2)*c2 + ua(2,3)*c3
                    PZAB = ua(3,1)*c1 + ua(3,2)*c2 + ua(3,3)*c3
                 

                    PABSQ  = PXAB ** 2 + PYAB ** 2 + PZAB ** 2
                    RABSQ  = DSQ(A)
                    DIFFSQ = RABSQ - PABSQ

                    IF ( ABS(DIFFSQ) .GT. ( RABSQ * TOL2 ) ) THEN

                       RXAB = RXI(A) - RXI(B)
                       RYAB = RYI(A) - RYI(B)
                       RZAB = RZI(A) - RZI(B)
                       c1 = RXAB*ua(1,1) + RYAB*ua(2,1) + RZAB*ua(3,1)
                       c2 = RXAB*ua(1,2) + RYAB*ua(2,2) + RZAB*ua(3,2)
                       c3 = RXAB*ua(1,3) + RYAB*ua(2,3) + RZAB*ua(3,3)
                       c1 = c1 - ANINT(c1/da(1))*da(1)
                       c2 = c2 - ANINT(c2/da(2))*da(2)
                       c3 = c3 - ANINT(c3/da(3))*da(3)
                       RXAB = ua(1,1)*c1 + ua(1,2)*c2 + ua(1,3)*c3
                       RYAB = ua(2,1)*c1 + ua(2,2)*c2 + ua(2,3)*c3
                       RZAB = ua(3,1)*c1 + ua(3,2)*c2 + ua(3,3)*c3

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

                       PXI(A) = PXI(A) + RMA * DX
                       PYI(A) = PYI(A) + RMA * DY
                       PZI(A) = PZI(A) + RMA * DZ
                       PXI(B) = PXI(B) - RMB * DX
                       PYI(B) = PYI(B) - RMB * DY
                       PZI(B) = PZI(B) - RMB * DZ

                       MOVING(A) = .TRUE.
                       MOVING(B) = .TRUE.
                       DONE = .FALSE.

                    ENDIF

                 ENDIF

300           CONTINUE

              DO 400 A = 1, N
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

           do a = 1, n
              r2(1,indx(a)) = PXI(a)
              r2(2,indx(a)) = PYI(a)
              r2(3,indx(a)) = PZI(a)
              r1(1,indx(a)) = RXI(a)
              r1(2,indx(a)) = RYI(a)
              r1(3,indx(a)) = RZI(a)
           end do





        RETURN
        END
c $Id$
