      Subroutine xc_gill96(tol_rho, fac, lfac, nlfac, rho, delrho, 
     &                      Amat, Cmat, nq, ipol, Ex, qwght,ldew,func)
c
C$Id: xc_gill96.F,v 1.1 1999-06-08 17:35:52 edo Exp $
c
      implicit none
c      
c
      double precision tol_rho, fac, Ex
      integer nq, ipol
      logical lfac, nlfac,ldew
      double precision func(*)  ! value of the functional [output]
c
c     Charge Density
c
      double precision rho(nq,ipol*(ipol+1)/2)
c
c     Charge Density Gradient
c
      double precision delrho(nq,3,ipol)
c
c     Quadrature Weights
c
      double precision qwght(nq)
c
c     Sampling Matrices for the XC Potential
c
      double precision Amat(nq,ipol), Cmat(nq,*)
c
c
      double precision BETA
      Parameter (BETA = 1d0/137d0)
c
c References:
c
c    Gill , Mol. Phys. 89, 433 (1996)
c
c***************************************************************************
c
      integer n
      double precision C, rho13, rho43, gamma, x, g,  dg,
     &      t,y ,f,ff
       f(y)=-BETA*y**(1.5d0)
       ff(y)= -1.5d0*BETA*sqrt(y)
c
c
c     Uniform electron gas constant
c
      C = -(1.5d0)*(0.75d0/acos(-1d0))**(1d0/3d0)
c
      if (ipol.eq.1) then
c
c        ======> SPIN-RESTRICTED <======
c
         do 10 n = 1, nq
            if (rho(n,1).lt.tol_rho) goto 10
c
c           Spin alpha:
c
            rho13 = (0.5d0*rho(n,1))**(1.d0/3.d0)
            rho43 = rho13**4 
            gamma = delrho(n,1,1)*delrho(n,1,1) +
     &              delrho(n,2,1)*delrho(n,2,1) +
     &              delrho(n,3,1)*delrho(n,3,1)
            if (dsqrt(gamma).gt.tol_rho)then
               gamma = 0.25d0 * gamma
               x = sqrt(gamma) / rho43
            else
               x = 0d0
            endif
c
            g = f(x)

            dg = ff(x)
c
            if (lfac) then
               Ex = Ex + 2d0*rho43*C*qwght(n)*fac
               if(ldew)func(n) = func(n) + 2.d0*rho43*C*fac
               Amat(n,1) = Amat(n,1) + (4d0/3d0)*rho13*C*fac
            endif
c
            if (nlfac) then
               Ex = Ex + 2d0*rho43*g*qwght(n)*fac
               if(ldew)func(n) = func(n) + 2.d0*rho43*g*fac
               Amat(n,1) = Amat(n,1) + (4d0/3d0)*rho13*(g-x*dg)*fac
            endif
c
            if (x.gt.tol_rho) then
               t = 0.5d0 * dg / sqrt(gamma) * fac
               Cmat(n,1) = Cmat(n,1) + t
            endif
c
c
 10      continue
c
      else
c
c        ======> SPIN-UNRESTRICTED <======
c
         do 20 n = 1, nq
            if (dabs(rho(n,1)).lt.tol_rho) goto 20
            if (dabs(rho(n,2)).lt.tol_rho) goto 25
c
c           Spin alpha:
c
            rho13 = rho(n,2)**(1.d0/3.d0)
            rho43 = rho13**4 
            gamma = delrho(n,1,1)*delrho(n,1,1) +
     &              delrho(n,2,1)*delrho(n,2,1) +
     &              delrho(n,3,1)*delrho(n,3,1)
            if (dsqrt(gamma).gt.tol_rho)then
               x = sqrt(gamma) / rho43
            else
               x = 0d0
            endif
c
            g = f(x)
            dg = ff(x)
c
            if (lfac) then
               Ex = Ex + rho43*C*qwght(n)*fac
               if (ldew)func(n) = func(n) + rho43*C*fac
               Amat(n,1) = Amat(n,1) + (4d0/3d0)*rho13*C*fac
            endif
c
            if (nlfac) then
               Ex = Ex + rho43*g*qwght(n)*fac
               if (ldew)func(n) = func(n) + rho43*g*fac
               Amat(n,1) = Amat(n,1) + (4d0/3d0)*rho13*(g-x*dg)*fac
            endif
c
            if (x.gt.tol_rho) then
               t = dg / sqrt(gamma) * fac
               Cmat(n,1) = Cmat(n,1) + t * 0.5d0

            endif
c
c
 25         continue
c
c           Spin beta:
c
            if (dabs(rho(n,3)).lt.tol_rho) goto 20
c
            rho13 = rho(n,3)**(1.d0/3.d0)
            rho43 = rho13**4 
            gamma = delrho(n,1,2)*delrho(n,1,2) +
     &              delrho(n,2,2)*delrho(n,2,2) +
     &              delrho(n,3,2)*delrho(n,3,2)
            if (dsqrt(gamma).gt.tol_rho)then
               x = sqrt(gamma) / rho43
            else
               x = 0d0
            endif
c
            g = f(x)
            dg = ff(x)
c
            if (lfac) then
               Ex = Ex + rho43*C*qwght(n)*fac
               if (ldew)func(n) = func(n) + rho43*C*fac
               Amat(n,2) = Amat(n,2) + (4d0/3d0)*rho13*C*fac
            endif
c
            if (nlfac) then
               Ex = Ex + rho43*g*qwght(n)*fac
               if (ldew)func(n) = func(n) +rho43*g*fac
               Amat(n,2) = Amat(n,2) + (4d0/3d0)*rho13*(g-x*dg)*fac
            endif
c
            if (x.gt.tol_rho) then
               t = dg / sqrt(gamma) * fac
               Cmat(n,3) = Cmat(n,3) + t * 0.5d0

            endif
c
c
 20      continue
c
      endif
c
      return
      end
