#ifndef SECOND_DERIV
      Subroutine xc_xmpw91(tol_rho, fac, lfac, nlfac, rho, delrho, 
     &                      Amat, Cmat, nq, ipol, Ex, qwght,ldew,func)     
#else
      Subroutine xc_xmpw91_d2(tol_rho, fac, lfac, nlfac, rho, delrho, 
     &                         Amat, Amat2, Cmat, Cmat2, nq, ipol, Ex,
     &                         qwght,ldew,func)
#endif
c
C$Id: xc_xmpw91.F,v 1.3 2003-02-19 02:18:05 edo Exp $
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
      double precision DPOW
      double precision BETA, TOLL, tol_rho, AX, pi,
     &     ONE, TWO, SIX, ONEHALF, FOURTHIRDS,CPW91,BETAPW91
      parameter (TOLL=1.E-20)

!      Parameter (BETA = 0.0042D0, CPW91=1.6455D0,DPOW=4) ! pw91 paper
!      Parameter (BETA = 0.0046D0, CPW91=1.6455D0,DPOW=3.73D0)!mpw91 paper
      Parameter (BETA = 0.00426D0, CPW91=1.6455D0,DPOW=3.73D0)!errata 
C
C     C. Adamo confirmed that there is a typo in the JCP paper
c     BETA is 0.00426 instead of 0.0046
C     adamo@ext.jussieu.fr
C

      Parameter (ONE=1.D0, TWO=2.D0, SIX=6.D0)
      Parameter (ONEHALF=0.5D0, FOURTHIRDS=4.D0/3.D0)
#ifdef SECOND_DERIV
c
c     Second Derivatives of the Exchange Energy Functional
c
      double precision Amat2(nq,NCOL_AMAT2), Cmat2(nq,NCOL_CMAT2)
      double precision rhom23, f2x, d2den,d2num
#endif
c
c References:
c
c
c***************************************************************************

      integer n
      double precision x, sinhm1,dsinhm1
      double precision rho13, rho43,  Xa,  Ha, denom, num,  
     &                  fprimex, d1num,d1den,
     &        gamma,gna,fx,x2a,ten6xd,expo,bbx2
c
      sinhm1(x)=log(x+dsqrt(1d0+x*x))
      dsinhm1(x)=1d0/dsqrt(1d0+x*x)
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
            if (gamma.gt.tol_rho**2)then
               gamma = 0.25d0*gamma
               xa = sqrt(gamma)/rho43
               x2a=xa*xa
               ten6xd=Xa**DPOW*1.d-6
               expo=exp(-CPW91*x2a)
               bbx2=(BETA-BETAPW91)*x2a
               Ha = sinhm1(Xa)
               denom = 1.d0/(1.d0 + 6d0*(beta*xa)*ha-ten6xd/ax)
               num = -BETA*x2a+bbx2*expo+ten6xd
               fx=num*denom
               d1num=-2.d0*xa*(beta+bbx2*expo*(1d0/x2a-CPW91))+
     +              ten6xd/xa*dpow
               d1den=6.d0*beta*(ha + xa*dsinhm1(xa)) -
     -              ten6xd/ax/xa*dpow
               fprimex=(d1num - d1den*fx)*denom
            else
               gamma = 0.d0
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
               Ex = Ex + 2d0*rho43*fx*qwght(n)*fac
               if(ldew)func(n) = func(n) + TWO*rho43*fx*fac
               Amat(n,1) = Amat(n,1) +
     +              FOURTHIRDS*rho13*(fx-xa*fprimex)*fac
               if (xa.gt.tol_rho)  then
                  Cmat(n,D1_GAA)=Cmat(n,D1_GAA)+
     +                 .5d0*fprimex/sqrt(gamma)*fac
               endif
            endif
c
#ifdef SECOND_DERIV
            rhom23 = rho13 / (0.5d0*rho(n,1))
            d2num=-2d0*(beta+(BETA-BETAPW91)*expo*
     *           (1d0-2d0*CPW91*xa*(1d0+xa-CPW91*x2a)))+
     +           ten6xd/x2a*dpow*(dpow-1d0)
            d2den=6.d0*beta*dsinhm1(xa)*(2d0-x2a/(1d0+x2a)) -
     -              ten6xd/ax/x2a*dpow*(dpow-1d0)
            f2x=(d2num-2d0*fprimex*d1num-fx*d2den)*denom
            f2x = d2num*denom + 4d0*BETA*xa*d1den*denom**2 -
     &            fx*d2den*denom +
     &            2d0*fx*(d1den**2)*denom**2
            f2x = d2num*denom - 2d0*d1num*d1den*denom**2 -
     &            fx*d2den*denom +
     &            2d0*fx*(d1den**2)*denom**2
            f2x = denom*(d2num -fx*d2den - 2d0*denom*d1den*
     *           (d1num - fx*d1den))
            f2x = denom*(d2num -fx*d2den - 2d0*d1den*fprimex)
            if (lfac)then
               Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) +
     &              (4d0/9d0)*rhom23*Ax*fac
            endif
            if (nlfac)then
               Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA)
     &              + (4d0/9d0)*rhom23*(fx-xa*fprimex+4d0*x2a*f2x)*fac
               Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA)
     &              - (4d0/3d0)*(rhom23**2/rho(n,1))*f2x*fac
               if (xa.gt.tol_rho) then
                  Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA)
     &                 - 0.25d0*gamma**(-1.5d0)*(fprimex-xa*f2x)*fac
               endif
            endif
#endif
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
               fx=-num*denom
               fprimex =-((2.d0*beta*xa-2.d0*bbx2*expo*(1d0/xa+CPW91)-
     -              ten6xd/xa*dpow) +
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
               Ex = Ex + rho43*fx*qwght(n)*fac
               if(ldew)func(n) = func(n) + rho43*fx*fac
               Amat(n,1) = Amat(n,1) +
     -              FOURTHIRDS*rho13*(fx - xa*fprimex)*fac
            endif
c
            if (Gna.gt.tol_rho)  then
               Cmat(n,D1_GAA) = Cmat(n,D1_GAA)+0.5d0*fprimex/Gna*fac
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
               fx=-num*denom
               fprimex =-((2.d0*beta*xa-2.d0*bbx2*expo*(1d0/xa+CPW91)-
     -              ten6xd/xa*dpow) +
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
               Ex = Ex + rho43*fx*qwght(n)*fac
               if(ldew)func(n) = func(n) + rho43*fx*fac
               Amat(n,2) = Amat(n,2) + 
     +              FOURTHIRDS*rho13*(fx - xa*fprimex)*fac
            endif
c
            if (Gna.gt.tol_rho)  then
               Cmat(n,D1_GBB) = Cmat(n,D1_GBB) +0.5d0*fprimex/Gna*fac
            endif
          endif
c        
   20    continue
c
      endif
c
      return
      end
#ifndef SECOND_DERIV
#define SECOND_DERIV
c
c     Compile source again for the 2nd derivative case
c
#include "xc_xmpw91.F"
#endif
