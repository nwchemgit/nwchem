      Subroutine xc_dirac(tol_rho, fac, lfac, nlfac, rho, Amat, nq, 
     &                    ipol, Ex, qwght, ldew, func)
c
C$Id: xc_dirac.F,v 1.8 1997-12-08 23:48:50 d3h449 Exp $
c
      Implicit none
c
      integer nq, ipol
      double precision fac, Ex
      logical ldew, lfac, nlfac
      double precision func(*)  ! value of the functional [output]
c
c     Charge Density
c
      double precision rho(nq,(ipol*(ipol+1))/2)
c
c     Quadrature Weights
c
      double precision qwght(nq)
c
c     Partial Derivatives of the Exchange Energy Functional
c
      double precision Amat(nq,ipol)
c
c     Compute the partial derivatives of the exchange functional of Dirac.
c
      double precision P1, P2, P3, P4, tol_rho
c
c     P1 =       -(3/PI)**(1/3)
c     P2 = -(3/4)*(3/PI)**(1/3)
c     P3 =       -(6/PI)**(1/3)
c     P4 = -(3/4)*(6/PI)**(1/3)
c
      Parameter (P1 = -0.9847450218426959D+00)
      Parameter (P2 = -0.7385587663820219D+00)
      Parameter (P3 = -0.1240700981798799D+01)
      Parameter (P4 = -0.9305257363490993D+00)
      double precision rho13, rho32, rho33, one_third, Exp
      Parameter (one_third = 1.d0/3.d0)
      integer n
c
      Exp = 0.d0
c
      if (ipol.eq.1)then
c
c        ======> SPIN-RESTRICTED <======
c
         do 10 n = 1, nq
            if (abs(rho(n,1)).gt.tol_rho)then
               rho13 = rho(n,1)**one_third
               Amat(n,1) = Amat(n,1) + rho13
               Exp = Exp + rho(n,1)*rho13*qwght(n)
               if(ldew)func(n) = func(n) + rho(n,1)*rho13*fac*P2
            endif
   10    continue
         call dscal(nq, P1*fac, Amat(1,1), 1)
c
         Ex = Ex + P2*Exp*fac
c
      else
c
c        ======> SPIN-UNRESTRICTED <======
c
         do 20 n = 1,nq
            if (abs(rho(n,1)).gt.tol_rho)then
               if (rho(n,2).lt.0.d0.or.rho(n,3).lt.0.d0)then
                  write(6,*)' n, nq, rho(n,1), rho(n,2), rho(n,3): ',
     &                        n, nq, rho(n,1), rho(n,2), rho(n,3)
                  call errquit(' negative density ',911)
               endif
               rho32=rho(n,2)**one_third
               rho33=rho(n,3)**one_third
               Amat(n,1) = Amat(n,1) + P3*rho32*fac
               Amat(n,2) = Amat(n,2) + P3*rho33*fac
               if (ldew)func(n) = func(n) + ( rho32*rho(n,2) + 
     &                                        rho33*rho(n,3)   )*P4*fac
               Exp = Exp + (rho32*rho(n,2)+rho33*rho(n,3))*qwght(n)
c      write(6,*)'from xc_dirac ; rho32, rho33, qwght(n): ',
c     &                           rho32, rho33, qwght(n)
            endif
   20    continue
c 
         Ex = Ex + P4*Exp*fac
c
      endif
c      write(6,*)'from xc_dirac ; p4, exp, fac, Ex: ', 
c     &                           p4, exp, fac, Ex
      return
      end
