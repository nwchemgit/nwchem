      Subroutine xc_xpw91(tol_rho, fac, lfac, nlfac, rho, delrho, 
     &                      Amat, Cmat, nq, ipol, Ex, qwght,ldew,func)     
c
C$Id: xc_xpw91.F,v 1.2 1999-01-05 22:06:34 edo Exp $
c
      implicit none
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
      double precision Amat(nq,ipol), Cmat(nq,3,ipol)
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
      double precision rho13a, rho43a, Gna, Xa, Sa, Ha, denom, num,  
     &                  fprimex, GNb, Xb, Sb, Hb, Tb, Numb, Gb, Pb, 
     &     rho13b, rho43b, gaa, g, gbb, gamma,fx,x2a,ten6xd,expo,bbx2
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
c
c           Spin alpha:
c
            rho13a = (0.5d0*rho(n,1))**(1.d0/3.d0)
            rho43a = rho13a**4 
            gamma = delrho(n,1,1)*delrho(n,1,1) +
     &           delrho(n,2,1)*delrho(n,2,1) +
     &           delrho(n,3,1)*delrho(n,3,1)
            if (dsqrt(gamma).gt.tol_rho)then
               GNa = dsqrt( 0.25D0*gamma )
               xa = GNa/rho43a
               x2a=xa*xa
               ten6xd=Xa**DPOW*1.d-6
               expo=exp(-CPW91*x2a)
               bbx2=(BETA-BETAPW91)*x2a
            else
               GNa = 0.d0
               Xa = 0.d0
               x2a=0.d0
               ten6xd=0.d0
               expo=1.d0
               bbx2=0.d0
            endif
            Ha = sinhm1(Xa)
            denom = 1.d0/(1.d0 + 6d0*(beta*xa)*ha-ten6xd/ax)
            num = BETA*x2a-bbx2*expo-ten6xd
            fx=num*denom
            fprimex = ((2.d0*beta*xa-2.d0*bbx2/xa*expo+
     +           2.d0*bbx2*xa*expo*CPW91-ten6xd/xa*dpow) -
     -           (6.d0*beta*ha +6.d0*beta*xa/sqrt(1.d0+x2a) -
     -           ten6xd/ax/xa*dpow)*fx)*denom
c
            if (lfac)then
               Ex = Ex + TWO*rho43a*AX*qwght(n)*fac
               if(ldew)func(n) = func(n) + TWO*rho43a*AX*fac
               Amat(n,1) = Amat(n,1) + FOURTHIRDS*rho13a*AX*fac
            endif
c
            if (nlfac)then
               Ex = Ex - TWO*rho43a*fx*qwght(n)*fac
               if(ldew)func(n) = func(n) - TWO*rho43a*fx*fac
               Amat(n,1) = Amat(n,1) - 
     -              FOURTHIRDS*rho13a*(fx - xa*fprimex)*fac
            endif
c
            if (Gna.gt.tol_rho)  then
               g = -.5d0*fprimex/Gna*fac
               Cmat(n,1,1) = Cmat(n,1,1) + g*delrho(n,1,1)
               Cmat(n,2,1) = Cmat(n,2,1) + g*delrho(n,2,1)
               Cmat(n,3,1) = Cmat(n,3,1) + g*delrho(n,3,1)
            endif
c
   10    continue
c
      else
c
c        ======> SPIN-UNRESTRICTED <======
c
         do 20 n = 1, nq
            if (dabs(rho(n,1)).lt.tol_rho) goto 20
c
c           Spin alpha:
c
            rho13a = rho(n,2)**(1.d0/3.d0)
            rho43a = rho13a**4 
            GNa = sqrt(   delrho(n,1,1)*delrho(n,1,1) +
     &                    delrho(n,2,1)*delrho(n,2,1) +
     &                    delrho(n,3,1)*delrho(n,3,1)   )
            if (GNa.gt.tol_rho)then
               Xa = GNa/rho43a
            else
               Xa = 0.d0
            endif
            Sa = sqrt(Xa*Xa + ONE)
            Ha = log(Xa + Sa)
            denom = 1/(ONE + SIX*(BETA*Xa)*Ha)
            num = -(BETA*Xa)*Xa*denom
            fprimex = ((SIX*(BETA*Xa)**2)*(Xa/Sa - Ha) - TWO*(BETA*Xa))*
     *           denom**2
c
c           Spin beta:
c
            rho13b = rho(n,3)**(1.d0/3.d0)
            rho43b = rho13b**4 
c
            GNb = sqrt(   delrho(n,1,2)*delrho(n,1,2) +
     &                    delrho(n,2,2)*delrho(n,2,2) +
     &                    delrho(n,3,2)*delrho(n,3,2)   )
            if (GNb.gt.tol_rho)then
               Xb = GNb/rho43b
            else
               Xb = 0.d0
            endif
            Sb = sqrt(Xb*Xb + ONE)
            Hb = log(Xb + Sb)
            Tb = 1/(ONE + SIX*(BETA*Xb)*Hb)
            Numb = -(BETA*Xb)*Xb*Tb
            Gb = AX + Numb
            Pb = ((SIX*(BETA*Xb)**2)*(Xb/Sb - Hb) - TWO*(BETA*Xb))*Tb**2
c
            if (lfac)then
               Amat(n,1) = Amat(n,1) + FOURTHIRDS*rho13a*AX*fac
               Amat(n,2) = Amat(n,2) + FOURTHIRDS*rho13b*AX*fac
            endif
c
            if (nlfac)then
               Amat(n,1) = Amat(n,1) + 
     +             FOURTHIRDS*rho13a*(num - Xa*fprimex)*fac
               Amat(n,2) = Amat(n,2) + 
     +              FOURTHIRDS*rho13b*(Numb - Xb*Pb)*fac
            endif
c
            if (lfac)then
               Ex = Ex + AX*(rho43a + rho43b)*qwght(n)*fac
               if (ldew)func(n) = func(n) +  
     *              AX*(rho43a + rho43b)*fac
            endif
c
            if (nlfac)then
               Ex = Ex + (rho43a*num + rho43b*Numb)*qwght(n)*fac
               if (ldew)func(n) = func(n) +  
     *              (rho43a*num + rho43b*Numb)*fac
            endif
            if (GNa.gt.tol_rho)then
               gaa = fprimex/GNa*fac
               Cmat(n,1,1) = Cmat(n,1,1) + gaa*delrho(n,1,1) 
               Cmat(n,2,1) = Cmat(n,2,1) + gaa*delrho(n,2,1) 
               Cmat(n,3,1) = Cmat(n,3,1) + gaa*delrho(n,3,1) 
            endif
            if (GNb.gt.tol_rho)then
               gbb = Pb/GNb*fac
               Cmat(n,1,2) = Cmat(n,1,2) + gbb*delrho(n,1,2) 
               Cmat(n,2,2) = Cmat(n,2,2) + gbb*delrho(n,2,2) 
               Cmat(n,3,2) = Cmat(n,3,2) + gbb*delrho(n,3,2) 
            endif
c        
   20    continue
c
      endif
c
      return
      end
