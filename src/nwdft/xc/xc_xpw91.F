      Subroutine xc_xpw91(tol_rho, fac, lfac, nlfac, rho, delrho, 
     &                      Amat, Cmat, nq, ipol, Ex, qwght,ldew,func)     
c
C$Id: xc_xpw91.F,v 1.4 1999-02-15 16:47:18 bjohnson Exp $
c
      implicit none
c
#include "dft2drv.fh"
c      
      double precision fac, Ex
      integer nq, ipol
      logical lfac, nlfac,ldew
      double precision func(*)  ! value of the functional [output]
c
c     Charge Density & Its Cube Root
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
c     Sampling Matrices for the XC Potential & Energy
c
      double precision Amat(nq,ipol), Cmat(nq,*)
c
c
c     Compute the partial derivatives of the exchange functional of Perdew91.
c
c     Becke & Perdew  Parameters
c
      integer DPOW
      double precision BETA, TOLL, tol_rho, AX, pi,
     &     ONE, TWO, SIX, ONEHALF, FOURTHIRDS,CPW91,BETAPW91
      parameter (TOLL=1.E-20)

      Parameter (BETA = 0.0042D0, CPW91=1.6455D0,DPOW=4)

      Parameter (ONE=1.D0, TWO=2.D0, SIX=6.D0)
      Parameter (ONEHALF=0.5D0, FOURTHIRDS=4.D0/3.D0)
c
c References:
c
c
c***************************************************************************

      integer n
      double precision x, sinhm1
      double precision rho13, rho43, Gna, Xa,  Ha, denom, num,  
     &                  fprimex, 
     &      g,  gamma,fx,x2a,ten6xd,expo,bbx2
c
      sinhm1(x)=log(x+dsqrt(1+x*x))
      pi=acos(-1.d0)
      BETAPW91=(pi*36.d0)**(-5.d0/3.d0)*5.d0
      AX=-(0.75d0/pi)**(1.d0/3.d0)*1.5d0
      if (ipol.eq.1 )then
c
c        ======> SPIN-RESTRICTED <======
c
         do 10 n = 1, nq
            if (rho(n,1).lt.tol_rho) goto 10

            rho13 = (0.5d0*rho(n,1))**(1.d0/3.d0)
            rho43 = rho13**4 
            gamma = delrho(n,1,1)*delrho(n,1,1) +
     &           delrho(n,2,1)*delrho(n,2,1) +
     &           delrho(n,3,1)*delrho(n,3,1)
            if (dsqrt(gamma).gt.tol_rho)then
               GNa = 0.5d0*dsqrt(gamma )
               xa = GNa/rho43
               x2a=xa*xa
               ten6xd=Xa**DPOW*1.d-6
               expo=exp(-CPW91*x2a)
               bbx2=(BETA-BETAPW91)*x2a
               Ha = sinhm1(Xa)
               denom = 1.d0/(1.d0 + 6d0*(beta*xa)*ha-ten6xd/ax)
               num = BETA*x2a-bbx2*expo-ten6xd
               fx=num*denom
               fprimex =((2.d0*beta*xa-2.d0*bbx2*expo*(1d0/xa+CPW91)-
     -              ten6xd/xa*dpow) -
     -              (6.d0*beta*(ha + xa/sqrt(1.d0+x2a)) -
     -              ten6xd/ax/xa*dpow)*fx)*denom
            else
               GNa = 0.d0
               Xa = 0.d0
               fx=0.d0
               fprimex=0.d0
            endif
c
            if (lfac)then
               Ex = Ex + 2d0*rho43*AX*qwght(n)*fac
               if(ldew)func(n) = func(n) + TWO*rho43*AX*fac
               Amat(n,1) = Amat(n,1) + FOURTHIRDS*rho13*AX*fac
            endif
c
            if (nlfac)then
               Ex = Ex - 2d0*rho43*fx*qwght(n)*fac
               if(ldew)func(n) = func(n) - TWO*rho43*fx*fac
               Amat(n,1) = Amat(n,1) - 
     -              FOURTHIRDS*rho13*(fx - xa*fprimex)*fac
            endif
c
            if (Gna.gt.tol_rho)  then
               g = -.5d0*fprimex/Gna*fac
               Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + g
#if 0
               Cmat(n,1,1) = Cmat(n,1,1) + g*delrho(n,1,1)
               Cmat(n,2,1) = Cmat(n,2,1) + g*delrho(n,2,1)
               Cmat(n,3,1) = Cmat(n,3,1) + g*delrho(n,3,1)
#endif
            endif
c
   10    continue
c
      else
c
c        ======> SPIN-UNRESTRICTED <======
c
         do 20 n = 1, nq
           if (rho(n,1).lt.tol_rho) goto 20
c
c           Spin alpha:
c
           if (rho(n,2).gt.tol_rho) then

            rho13 = rho(n,2)**(1.d0/3.d0)
            rho43 = rho13**4 
            gamma =    delrho(n,1,1)*delrho(n,1,1) +
     &                    delrho(n,2,1)*delrho(n,2,1) +
     &                    delrho(n,3,1)*delrho(n,3,1)  
            if (gamma.gt.tol_rho**2)then
               GNa = dsqrt( gamma )
               xa = GNa/rho43
               x2a=xa*xa
               ten6xd=Xa**DPOW*1.d-6
               expo=exp(-CPW91*x2a)
               bbx2=(BETA-BETAPW91)*x2a
               Ha = sinhm1(Xa)
               denom = 1.d0/(1.d0 + 6d0*(beta*xa)*ha-ten6xd/ax)
               num = BETA*x2a-bbx2*expo-ten6xd
               fx=num*denom
               fprimex =((2.d0*beta*xa-2.d0*bbx2*expo*(1d0/xa+CPW91)-
     -              ten6xd/xa*dpow) -
     -              (6.d0*beta*(ha + xa/sqrt(1.d0+x2a)) -
     -              ten6xd/ax/xa*dpow)*fx)*denom
            else
               GNa = 0.d0
               Xa = 0.d0
               fx=0.d0
               fprimex=0.d0
            endif
c
            if (lfac)then
               Ex = Ex + rho43*AX*qwght(n)*fac
               if(ldew)func(n) = func(n) + rho43*AX*fac
               Amat(n,1) = Amat(n,1) + FOURTHIRDS*rho13*AX*fac
            endif
c
            if (nlfac)then
               Ex = Ex - rho43*fx*qwght(n)*fac
               if(ldew)func(n) = func(n) - rho43*fx*fac
               Amat(n,1) = Amat(n,1) - 
     -              FOURTHIRDS*rho13*(fx - xa*fprimex)*fac
            endif
c
            if (Gna.gt.tol_rho)  then
               g = -fprimex/Gna*fac
               Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + g*0.5d0
#if 0
               Cmat(n,1,1) = Cmat(n,1,1) + g*delrho(n,1,1)
               Cmat(n,2,1) = Cmat(n,2,1) + g*delrho(n,2,1)
               Cmat(n,3,1) = Cmat(n,3,1) + g*delrho(n,3,1)
#endif
            endif
            endif
c
c           Spin beta:
c
            if (rho(n,3).gt.tol_rho) then

            rho13=rho(n,3)**(1.d0/3.d0)
            rho43 = rho13**4 
c
            gamma =   delrho(n,1,2)*delrho(n,1,2) +
     &           delrho(n,2,2)*delrho(n,2,2) +
     &           delrho(n,3,2)*delrho(n,3,2)   
            if (gamma.gt.tol_rho**2)then
               GNa = sqrt(gamma)
               xa = GNa/rho43
               x2a=xa*xa
               ten6xd=Xa**DPOW*1.d-6
               expo=exp(-CPW91*x2a)
               bbx2=(BETA-BETAPW91)*x2a
               Ha = sinhm1(Xa)
               denom = 1.d0/(1.d0 + 6d0*(beta*xa)*ha-ten6xd/ax)
               num = BETA*x2a-bbx2*expo-ten6xd
               fx=num*denom
               fprimex =((2.d0*beta*xa-2.d0*bbx2*expo*(1d0/xa+CPW91)-
     -              ten6xd/xa*dpow) -
     -              (6.d0*beta*(ha + xa/sqrt(1.d0+x2a)) -
     -              ten6xd/ax/xa*dpow)*fx)*denom
            else
               GNa = 0.d0
               Xa = 0.d0
               fx=0.d0
               fprimex=0.d0
            endif
c
            if (lfac)then
               Ex = Ex + rho43*AX*qwght(n)*fac
               if(ldew)func(n) = func(n) + rho43*AX*fac
               Amat(n,2) = Amat(n,2) + FOURTHIRDS*rho13*AX*fac
            endif
c
            if (nlfac)then
               Ex = Ex - rho43*fx*qwght(n)*fac
               if(ldew)func(n) = func(n) - rho43*fx*fac
               Amat(n,2) = Amat(n,2) - 
     -              FOURTHIRDS*rho13*(fx - xa*fprimex)*fac
            endif
c
            if (Gna.gt.tol_rho)  then
               g = -fprimex/Gna*fac
               Cmat(n,D1_GBB) = Cmat(n,D1_GBB) + g*0.5d0
#if 0
               Cmat(n,1,2) = Cmat(n,1,2) + g*delrho(n,1,1)
               Cmat(n,2,2) = Cmat(n,2,2) + g*delrho(n,2,1)
               Cmat(n,3,2) = Cmat(n,3,2) + g*delrho(n,3,1)
#endif
            endif
          endif
c        
   20    continue
c
      endif
c
      return
      end
