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
        INTEGER     IT, A, B, I,J
        REAL*8      PXAB, PYAB, PZAB, PABSQ
        REAL*8      RXAB, RYAB, RZAB, RABSQ, DIFFSQ, RPAB
        REAL*8      GAB, DX, DY, DZ, TOL2
        REAL*8      RPTOL, RMA, RMB
        PARAMETER ( RPTOL = 1.0E-6 )

        real*8 ua(3,3),ub(3,3),volume
        real*8 c1,c2,c3

*     **** external functions ****
      real*8   lattice_unita
      external lattice_unita


*     ***** Determine the unit lattice vectors and distances ******
      do j=1,3
      do i=1,3
         ua(i,j) = lattice_unita(i,j)
      end do
      end do
      ub(1,1) = ua(2,2)*ua(3,3) - ua(3,2)*ua(2,3)
      ub(2,1) = ua(3,2)*ua(1,3) - ua(1,2)*ua(3,3)
      ub(3,1) = ua(1,2)*ua(2,3) - ua(2,2)*ua(1,3)
      ub(1,2) = ua(2,3)*ua(3,1) - ua(3,3)*ua(2,1)
      ub(2,2) = ua(3,3)*ua(1,1) - ua(1,3)*ua(3,1)
      ub(3,2) = ua(1,3)*ua(2,1) - ua(2,3)*ua(1,1)
      ub(1,3) = ua(2,1)*ua(3,2) - ua(3,1)*ua(2,2)
      ub(2,3) = ua(3,1)*ua(1,2) - ua(1,1)*ua(3,2)
      ub(3,3) = ua(1,1)*ua(2,2) - ua(2,1)*ua(1,2)
      volume = ua(1,1)*ub(1,1)
     >       + ua(2,1)*ub(2,1)
     >       + ua(3,1)*ub(3,1)
      volume = 1.0d0/volume
      call dscal(9,volume,ub,1)


        
        TOL2   = 2.0 * tol



        do A = 1, n
              RXI(A) = r1(1,indx(A))
              RYI(A) = r1(2,indx(A))
              RZI(A) = r1(3,indx(A))
              PXI(A) = r2(1,indx(A))
              PYI(A) = r2(2,indx(A))
              PZI(A) = r2(3,indx(A))

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
                    c1 = PXAB*ub(1,1) + PYAB*ub(2,1) + PZAB*ub(3,1)
                    c2 = PXAB*ub(1,2) + PYAB*ub(2,2) + PZAB*ub(3,2)
                    c3 = PXAB*ub(1,3) + PYAB*ub(2,3) + PZAB*ub(3,3)
                    c1 = c1 - ANINT(c1)
                    c2 = c2 - ANINT(c2)
                    c3 = c3 - ANINT(c3)
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
                       c1 = RXAB*ub(1,1) + RYAB*ub(2,1) + RZAB*ub(3,1)
                       c2 = RXAB*ub(1,2) + RYAB*ub(2,2) + RZAB*ub(3,2)
                       c3 = RXAB*ub(1,3) + RYAB*ub(2,3) + RZAB*ub(3,3)
                       c1 = c1 - ANINT(c1)
                       c2 = c2 - ANINT(c2)
                       c3 = c3 - ANINT(c3)
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
