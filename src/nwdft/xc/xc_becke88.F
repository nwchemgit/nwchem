      Subroutine xc_becke88(tol_rho, fac, lfac, nlfac, rho, delrho, 
     &                      Amat, Cmat, nq, ipol, Ex, qwght)     
c
C$Id: xc_becke88.F,v 1.9 1997-12-08 23:48:50 d3h449 Exp $
c
      implicit none
c      
      double precision fac, Ex
      integer nq, ipol
      logical lfac, nlfac
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
c     Compute the partial derivatives of the exchange functional of Becke.
c
c     Becke Empirical Parameter
c
      double precision BETA, TOLL, tol_rho, P2, 
     &     ONE, TWO, SIX, ONEHALF, FOURTHIRDS

      Parameter (BETA = 0.0042D0, TOLL=1.E-20)
c
c     P2 = -(3/4)*(6/PI)**(1/3)
c
      Parameter (P2 = -0.9305257363490993D+00)

      Parameter (ONE=1.D0, TWO=2.D0, SIX=6.D0)
      Parameter (ONEHALF=0.5D0, FOURTHIRDS=4.D0/3.D0)
c
c References:
c
c    Becke, Phys. Rev. A 38, 3098 (1988)
c    Johnson, Gill & Pople, J. Chem. Phys. 98, 5612 (1993)
c
c***************************************************************************

      integer n
      double precision x, sinhm1
      double precision CD13a, CD43a, Gna, Xa, Sa, Ha, Ta, Ca, P2LDA, 
     &                 Ga, Pa, GNb, Xb, Sb, Hb, Tb, Cb, Gb, Pb, CD13b, 
     &                 CD43b, gaa, g, gbb, gamma
c
      sinhm1(x)=log(x+dsqrt(1+x*x))
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
     &              delrho(n,2,1)*delrho(n,2,1) +
     &              delrho(n,3,1)*delrho(n,3,1)
            if (dsqrt(gamma).gt.tol_rho)then
               GNa = dsqrt( 0.25D0*gamma )
               Xa = GNa/CD43a
            else
               GNa = 0.d0
               Xa = 0.d0
            endif
            Sa = sqrt(Xa*Xa + ONE)
            Ha = sinhm1(Xa)
            Ta = 1/(ONE + SIX*(BETA*Xa)*Ha)
            Ca = -(BETA*Xa)*Xa*Ta
            P2LDA = P2
            Ga = P2LDA + Ca
            Pa = ((SIX*(BETA*Xa)**2)*(Xa/Sa - Ha) - TWO*(BETA*Xa))*Ta**2
c
            if (lfac)then
               Amat(n,1) = Amat(n,1) + FOURTHIRDS*CD13a*P2LDA*fac
            endif
c
            if (nlfac)then
               Amat(n,1) = Amat(n,1) + FOURTHIRDS*CD13a*(Ca - Xa*Pa)*fac
            endif
c
            if (lfac)then
               Ex = Ex + TWO*CD43a*P2LDA*qwght(n)*fac
            endif
c
            if (nlfac)then
               Ex = Ex + TWO*CD43a*Ca*qwght(n)*fac
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
c
            Xa = GNa/CD43a
            Sa = sqrt(Xa*Xa + ONE)
            Ha = log(Xa + Sa)
            Ta = 1/(ONE + SIX*(BETA*Xa)*Ha)
            Ca = -(BETA*Xa)*Xa*Ta
            P2LDA = P2
            Ga = P2LDA + Ca
            Pa = ((SIX*(BETA*Xa)**2)*(Xa/Sa - Ha) - TWO*(BETA*Xa))*Ta**2
c
c           Spin beta:
c
            CD13b = rho(n,3)**(1.d0/3.d0)
            CD43b = CD13b**4 
c
            GNb = sqrt(   delrho(n,1,2)*delrho(n,1,2) +
     &                    delrho(n,2,2)*delrho(n,2,2) +
     &                    delrho(n,3,2)*delrho(n,3,2)   )

            Xb = GNb/CD43b
            Sb = sqrt(Xb*Xb + ONE)
            Hb = log(Xb + Sb)
            Tb = 1/(ONE + SIX*(BETA*Xb)*Hb)
            Cb = -(BETA*Xb)*Xb*Tb
            Gb = P2LDA + Cb
            Pb = ((SIX*(BETA*Xb)**2)*(Xb/Sb - Hb) - TWO*(BETA*Xb))*Tb**2
c
            if (lfac)then
               Amat(n,1) = Amat(n,1) + FOURTHIRDS*CD13a*P2LDA*fac
               Amat(n,2) = Amat(n,2) + FOURTHIRDS*CD13b*P2LDA*fac
            endif
c
            if (nlfac)then
               Amat(n,1) = Amat(n,1) + FOURTHIRDS*CD13a*(Ca - Xa*Pa)*fac
               Amat(n,2) = Amat(n,2) + FOURTHIRDS*CD13b*(Cb - Xb*Pb)*fac
            endif
c
            if (lfac)then
               Ex = Ex + P2LDA*(CD43a + CD43b)*qwght(n)*fac
            endif
c
            if (nlfac)then
               Ex = Ex + (CD43a*Ca + CD43b*Cb)*qwght(n)*fac
            endif
cwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww
c            write(6,*)' ex, cd13a, cd43a, cd13b, cd43b, ca, cb: ',
c     &                  ex, cd13a, cd43a, cd13b, cd43b, ca, cb
c
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
