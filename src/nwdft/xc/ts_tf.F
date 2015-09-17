cc AJL/Begin/FDE
      Subroutine ts_tf(tol_rho, fac, lfac, nlfac, rho, Amat, nq, 
     &                    ipol, Ts, qwght, ldew, func)
c
      Implicit none
#include "errquit.fh"
c
#include "stdio.fh"
c
      integer nq, ipol
      double precision fac, Ts
      logical ldew, lfac, nlfac
      double precision func(*)  ! value of the functional [output]
      double precision two_rho32, two_rho33
c
c     Charge Density
c
      double precision rho(nq,(ipol*(ipol+1))/2)
c
c     Quadrature Weights
c
      double precision qwght(nq)
c
c     Partial First Derivatives of the Kinetic Energy Functional
c
      double precision Amat(nq,ipol)
c
c     Compute the partial derivatives of the Kinetic functional of Dirac.
c
      double precision P1, P2, P3, P4, tol_rho
c
c Daniel (5-23-12): P1 and P2 are for restricted calculations, P3 
c and P4 are for unrestricted calculations.
c
c     P1 =            (1/2) *(3*PI^2)**(1/3)
c     P2 =           (3/10) *(3*PI^2)**(1/3)
c     P3 = 2**(2/3)*  (1/2) *(3*PI^2)**(1/3)
c     P4 = 2**(2/3)* (3/10) *(3*PI^2)**(1/3) 
c
c Calculated 14/10/2014 with Wolfram Alpha
      Parameter (P1 = 0.15468338631400680D+01)
      Parameter (P2 = 0.92810031788404078D+00)
      Parameter (P3 = 0.24554457015685778D+01)
      Parameter (P4 = 0.14732674209411467D+01)
      double precision rho13, rhom23, two_third, term
      double precision rho32, rho33
      Parameter (two_third = 2.d0/3.d0)
      integer n
c
      term = 0.d0
c
      if (ipol.eq.1)then
c
c        ======> SPIN-RESTRICTED <======
c
         do n = 1, nq
c            write(luout,*) rho(n,1)
            if (rho(n,1).gt.tol_rho)then
               rho13=rho(n,1)**two_third
               Amat(n,1) = Amat(n,1) + rho13
               term = term + rho(n,1)*rho13*qwght(n)
               if(ldew) func(n) = func(n) + rho(n,1)*rho13*fac*P2
            endif
         enddo
         call dscal(nq, P1*fac, Amat(1,1), 1)
c 
         Ts = Ts + P2*term*fac
c
      else
c
c        ======> SPIN-UNRESTRICTED <======
c
         do n = 1,nq
c            write(luout,*) rho(n,1), rho(n,2), rho(n,3)
            if (rho(n,1).gt.tol_rho)then

               rho32 = 0.d0
               rho33 = 0.d0

               if (rho(n,2).gt.tol_rho) rho32 = rho(n,2)**two_third
               if (rho(n,3).gt.tol_rho) rho33 = rho(n,3)**two_third

               Amat(n,1) = Amat(n,1) + P3*rho32*fac
               Amat(n,2) = Amat(n,2) + P3*rho33*fac

               if (ldew)func(n) = func(n) + ((rho(n,2)*rho32)
     &              + (rho(n,3)*rho33))*P4*fac

               term = term+((rho(n,2)*rho32)
     &              + (rho(n,3)*rho33))*qwght(n)

c      write(LuOut,*)'from ts_tf ; rho32, rho33, qwght(n): ',
c     &                           rho32, rho33, qwght(n)
            endif
         enddo
c 
         Ts = Ts + P4*term*fac
      endif
c      write(LuOut,*)'from ts_tf ; p4, term, fac, Ex: ', 
c     &                           P4, term, fac, Ts 
      return
      end
c
cc AJL/End
