      Subroutine xc_becke88(fac,
     *     rho,delrho,Amat,Cmat,nq,ipol,Ex,qwght)

C$Id: xc_becke88.F,v 1.4 1996-05-04 21:08:34 og845 Exp $
      implicit none
      
      double precision fac,Ex
      integer nq,ipol

      logical grad_corr
c
c     logical grad_corr true if only GCA part needed
c
c--> Charge Density & Its Cube Root

      double precision rho(ipol*(ipol+1)/2,nq)

c--> Charge Density Gradient

      double precision delrho(ipol,3,nq)

c--> Quadrature Weights

      double precision qwght(nq)

c--> Sampling Matrices for the XC Potential & Energy

      double precision Amat(ipol,nq),Cmat(ipol,3,nq)

c
c Compute the partial derivatives of the exchange functional of Becke.
c
c--> Becke Empirical Parameter

      double precision BETA,TOLL,P2,
     ,     ONE,TWO,SIX,ONEHALF,FOURTHIRDS

      Parameter (BETA = 0.0042D0,TOLL=1.E-20)

c--> P2 = -(3/4)*(6/PI)**(1/3)

      Parameter (P2 = -0.9305257363490993D+00)

      Parameter (ONE=1.D0,TWO=2.D0,SIX=6.D0)
      Parameter (ONEHALF=0.5D0,FOURTHIRDS=4.D0/3.D0)

c References:
c
c    Becke, Phys. Rev. A 38, 3098 (1988)
c    Johnson, Gill & Pople, J. Chem. Phys. 98, 5612 (1993)
c
c***************************************************************************

      integer n
      double precision x,sinhm1
      double precision CD13a,CD43a,Gna,Xa,Sa,Ha,Ta,Ca,P2LDA,Ga,Pa,
     ,     GNb,Xb,Sb,Hb,Tb,Cb,Gb,Pb,CD13b,CD43b,
     ,     gaa,g,gbb,gamma
c

      sinhm1(x)=log(x+dsqrt(1+x*x))
      grad_corr=.false.
      if(abs(fac-1.d0).gt.1.e-8) grad_corr=.true.
      if( ipol.eq.1 )then

c                      ======> SPIN-RESTRICTED <======

       do 10 n = 1,nq
       if(rho(1,n).lt.TOLL) goto 10
c Spin alpha:
        CD13a = (0.5d0*rho(1,n))**(1.d0/3.d0)
        CD43a = CD13a**4 
        gamma=delrho(1,1,n)*delrho(1,1,n)+
     +       delrho(1,2,n)*delrho(1,2,n)+
     +       delrho(1,3,n)*delrho(1,3,n)
        if(dsqrt(gamma).gt.TOLL) then
          GNa = dsqrt( 0.25D0*gamma )
          Xa = GNa/CD43a
        else
          GNa=0.d0
          Xa=0.d0
        endif
        Sa = sqrt(Xa*Xa + ONE)
        Ha = sinhm1(Xa)
        Ta = 1/(ONE + SIX*(BETA*Xa)*Ha)
        Ca = -(BETA*Xa)*Xa*Ta
        P2LDA=P2
        if(grad_corr) P2LDA=0.d0
        Ga = P2LDA + Ca
        Pa = ( (SIX*(BETA*Xa)**2)*(Xa/Sa - Ha) - TWO*(BETA*Xa) )*Ta**2

c ...

        Amat(1,n) = Amat(1,n) + FOURTHIRDS*CD13a*(Ga - Xa*Pa)*fac
        if(Gna.gt.TOLL)  then
          g         = .5d0*Pa/Gna*fac
          Cmat(1,1,n) = Cmat(1,1,n) + g*delrho(1,1,n)
          Cmat(1,2,n) = Cmat(1,2,n) + g*delrho(1,2,n)
          Cmat(1,3,n) = Cmat(1,3,n) + g*delrho(1,3,n)
        endif
        Ex      = Ex + TWO*CD43a*Ga*qwght(n)*fac


   10  continue

      else

c                      ======> SPIN-UNRESTRICTED <======

       do 20 n = 1,nq
         if(dabs(rho(1,n)).lt.TOLL) goto 20

c Spin alpha:

        CD13a = rho(2,n)**(1.d0/3.d0)
        CD43a = CD13a**4 
        GNa = sqrt(
     (  delrho(1,1,n)*delrho(1,1,n)+
     +       delrho(1,2,n)*delrho(1,2,n)+
     +       delrho(1,3,n)*delrho(1,3,n))

        Xa = GNa/CD43a
        Sa = sqrt(Xa*Xa + ONE)
        Ha = log(Xa + Sa)
        Ta = 1/(ONE + SIX*(BETA*Xa)*Ha)
        Ca = -(BETA*Xa)*Xa*Ta
        P2LDA=P2
        if(grad_corr) P2LDA=0.d0
        Ga = P2LDA + Ca
        Pa = ( (SIX*(BETA*Xa)**2)*(Xa/Sa - Ha) - TWO*(BETA*Xa) )*Ta**2

c Spin beta:

        CD13b = rho(3,n)**(1.d0/3.d0)
        CD43b = CD13b**4 

        GNb = sqrt(
     (  delrho(2,1,n)*delrho(2,1,n)+
     +       delrho(2,2,n)*delrho(2,2,n)+
     +       delrho(2,3,n)*delrho(2,3,n))

        Xb = GNb/CD43b
        Sb = sqrt(Xb*Xb + ONE)
        Hb = log(Xb + Sb)
        Tb = 1/(ONE + SIX*(BETA*Xb)*Hb)
        Cb = -(BETA*Xb)*Xb*Tb
        Gb = P2LDA + Cb
        Pb = ( (SIX*(BETA*Xb)**2)*(Xb/Sb - Hb) - TWO*(BETA*Xb) )*Tb**2

c ...

        Amat(1,n) = Amat(1,n) +FOURTHIRDS*CD13a*(Ga - Xa*Pa)*fac
        Amat(2,n) = Amat(2,n) + FOURTHIRDS*CD13b*(Gb - Xb*Pb)*fac
        if(GNa.gt.TOLL) then
          gaa = Pa/GNa*fac
          Cmat(1,1,n) = Cmat(1,1,n) + gaa*delrho(1,1,n) 
          Cmat(1,2,n) = Cmat(1,2,n) + gaa*delrho(1,2,n) 
          Cmat(1,3,n) = Cmat(1,3,n) + gaa*delrho(1,3,n) 
        endif
        if(GNb.gt.TOLL) then
          gbb = Pb/GNb*fac
          Cmat(2,1,n) = Cmat(2,1,n) + gbb*delrho(2,1,n) 
          Cmat(2,2,n) = Cmat(2,2,n) + gbb*delrho(2,2,n) 
          Cmat(2,3,n) = Cmat(2,3,n) + gbb*delrho(2,3,n) 
        endif
        
        Ex      = Ex + (CD43a*Ga + CD43b*Gb)*qwght(n)*fac

   20  continue

      end if

      end
