      Subroutine xc_xpw91(tol_rho, fac, lfac, nlfac, rho, delrho, 
     &                      Amat, Cmat, nq, ipol, Ex, qwght,ldew,func)     
c
C$Id: xc_xpw91.F,v 1.1 1999-01-05 00:44:02 edo Exp $
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
      double precision BETA, TOLL, tol_rho, AX, 
     &     ONE, TWO, SIX, ONEHALF, FOURTHIRDS,TEN6,CPW91,BETAPW91
      parameter (TOLL=1.E-20)

      Parameter (BETA = 0.0042D0, CPW91=1.6455D0,DPOW=4,
     ,     TEN6=1.d-6)
c
c     AX = -(3/4)*(6/PI)**(1/3)
c
      Parameter (AX = -0.9305257363490993D+00)

      Parameter (ONE=1.D0, TWO=2.D0, SIX=6.D0)
      Parameter (ONEHALF=0.5D0, FOURTHIRDS=4.D0/3.D0)
c
c References:
c
c
c***************************************************************************

      integer n
      double precision x, sinhm1
      double precision CD13a, CD43a, Gna, Xa, Sa, Ha, denom, Ca,  
     &                 Ga, Pa, GNb, Xb, Sb, Hb, Tb, Cb, Gb, Pb, CD13b, 
     &                 CD43b, gaa, g, gbb, gamma
c
      sinhm1(x)=log(x+dsqrt(1+x*x))
      BETAPW91=(acos(-1.d0)*36.d0)**(-5.d0/3.d0)*5.d0
      if (ipol.eq.1 )then
c
c        ======> SPIN-RESTRICTED <======
c
         do 10 n = 1, nq
            if (rho(n,1).lt.tol_rho) goto 10
c
c           Spin alpha:
c
            CD13a = (0.5d0*rho(n,1))**(1.d0/3.d0)
            CD43a = CD13a**4 
            gamma = delrho(n,1,1)*delrho(n,1,1) +
     &           delrho(n,2,1)*delrho(n,2,1) +
     &           delrho(n,3,1)*delrho(n,3,1)
            if (dsqrt(gamma).gt.tol_rho)then
               GNa = dsqrt( 0.25D0*gamma )
               Xa = GNa/CD43a
            else
               GNa = 0.d0
               Xa = 0.d0
            endif
            Sa = sqrt(Xa*Xa + ONE)
            Ha = sinhm1(Xa)
            denom = 1.d0/(ONE + SIX*(BETA*Xa)*Ha-Xa**DPOW*TEN6/AX)
            Ca = -(BETA*Xa)*Xa*denom
            Ga = AX + Ca
            Pa = ((SIX*(BETA*Xa)**2)*(Xa/Sa-Ha)-TWO*(BETA*Xa))*
     *           denom**2
c
            if (lfac)then
               Amat(n,1) = Amat(n,1) + FOURTHIRDS*CD13a*AX*fac
            endif
c
            if (nlfac)then
               Amat(n,1) = Amat(n,1) + FOURTHIRDS*CD13a*(Ca - Xa*Pa)*fac
            endif
c
            if (lfac)then
               Ex = Ex + TWO*CD43a*AX*qwght(n)*fac
               if(ldew)func(n) = func(n) + TWO*CD43a*AX*fac
            endif
c
            if (nlfac)then
               Ex = Ex + TWO*CD43a*Ca*qwght(n)*fac
               if(ldew)func(n) = func(n) + TWO*CD43a*Ca*fac
            endif
c
            if (Gna.gt.tol_rho)  then
               g = .5d0*Pa/Gna*fac
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
            CD13a = rho(n,2)**(1.d0/3.d0)
            CD43a = CD13a**4 
            GNa = sqrt(   delrho(n,1,1)*delrho(n,1,1) +
     &                    delrho(n,2,1)*delrho(n,2,1) +
     &                    delrho(n,3,1)*delrho(n,3,1)   )
            if (GNa.gt.tol_rho)then
               Xa = GNa/CD43a
            else
               Xa = 0.d0
            endif
            Sa = sqrt(Xa*Xa + ONE)
            Ha = log(Xa + Sa)
            denom = 1/(ONE + SIX*(BETA*Xa)*Ha)
            Ca = -(BETA*Xa)*Xa*denom
            Ga = AX + Ca
            Pa = ((SIX*(BETA*Xa)**2)*(Xa/Sa - Ha) - TWO*(BETA*Xa))*
     *           denom**2
c
c           Spin beta:
c
            CD13b = rho(n,3)**(1.d0/3.d0)
            CD43b = CD13b**4 
c
            GNb = sqrt(   delrho(n,1,2)*delrho(n,1,2) +
     &                    delrho(n,2,2)*delrho(n,2,2) +
     &                    delrho(n,3,2)*delrho(n,3,2)   )
            if (GNb.gt.tol_rho)then
               Xb = GNb/CD43b
            else
               Xb = 0.d0
            endif
            Sb = sqrt(Xb*Xb + ONE)
            Hb = log(Xb + Sb)
            Tb = 1/(ONE + SIX*(BETA*Xb)*Hb)
            Cb = -(BETA*Xb)*Xb*Tb
            Gb = AX + Cb
            Pb = ((SIX*(BETA*Xb)**2)*(Xb/Sb - Hb) - TWO*(BETA*Xb))*Tb**2
c
            if (lfac)then
               Amat(n,1) = Amat(n,1) + FOURTHIRDS*CD13a*AX*fac
               Amat(n,2) = Amat(n,2) + FOURTHIRDS*CD13b*AX*fac
            endif
c
            if (nlfac)then
               Amat(n,1) = Amat(n,1) + FOURTHIRDS*CD13a*(Ca - Xa*Pa)*fac
               Amat(n,2) = Amat(n,2) + FOURTHIRDS*CD13b*(Cb - Xb*Pb)*fac
            endif
c
            if (lfac)then
               Ex = Ex + AX*(CD43a + CD43b)*qwght(n)*fac
               if (ldew)func(n) = func(n) +  
     *              AX*(CD43a + CD43b)*fac
            endif
c
            if (nlfac)then
               Ex = Ex + (CD43a*Ca + CD43b*Cb)*qwght(n)*fac
               if (ldew)func(n) = func(n) +  
     *              (CD43a*Ca + CD43b*Cb)*fac
            endif
            if (GNa.gt.tol_rho)then
               gaa = Pa/GNa*fac
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
