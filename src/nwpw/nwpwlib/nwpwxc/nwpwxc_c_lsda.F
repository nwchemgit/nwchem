C> \ingroup nwpwxc
C> @{
C> 
C> \file nwpwxc_c_lsda.F
C> Another implementation of the PW91 LDA correlation functional
C>
C> @}
C>
C> \ingroup nwpwxc_priv
C> @{
C>
C> \brief An alternative implementation of the functional 
C> `nwpwxc_c_pw91lda`
C>
C> This stems from the original authors of the functional and hence
C> it follows different calling conventions from the other routines
C> in the library. For this reason it is also only meant to be used
C> as an internal component of PW91 related functionals.
C>
c
c Perdew 91 local correlation functional at one grid point 
c

      Subroutine nwpwxc_c_lsda(tol_rho,
     R     rs,zeta,PotLC,dLdS,dLdZ,d2LdSS,d2LdSZ,d2LdZZ)

C     Evaluate the Perdew 91 local correlation functional and its
C     derivatives at one point. 
      implicit none
      double precision tol_rho
      double precision rs
      double precision PotLC,dLdS, dLdZ, d2dSS, d2dSZ, d2dZZ
      double precision eps0c(6), eps1c(6), epsc(6)
      double precision F1, F2, F3, F4, F6, F8, F9, F12, F24, F36
      double precision GammaI,Zeta,FZeta,dfZdz,d2fZdz 
      double precision EU,dEUdRS,d2UdRS
      double precision EP,dEPdRS,d2PdRS
      double precision AlphaM,dAMdRS,d2AdRS
      double precision GZ, HZ, dGZ, dHZ, d2GZ, d2HZ
      double precision d2LdSS,d2LdSZ,d2LdZZ,fzzi,
     Z     z2,z3,z4,third,pi
 
      data eps0c/0.03109070D0,0.21370D0, 7.5957D0,3.5876D0,1.6382D0,
     &         0.49294D0/
      data eps1c/0.01554535D0,0.20548D0,14.1189D0,6.1977D0,3.3662D0,
     &         0.62517D0/
      data epsc/0.01688690D0,0.11125D0,10.3570D0,3.6231D0,0.88026D0,
     &         0.49671D0/
      data F1/1.0d0/, F2/2.0d0/, F3/3.0d0/, F4/4.0d0/, F6/6.0d0/,
     &  F8/8.0d0/, F9/9.0d0/, F12/12.0d0/, F24/24.0d0/, F36/36.0d0/

      Pi = F4*ATan(F1)
      Third = F1 / F3


      epsc(1) = F1 / (F6*Pi*Pi)
      FZZI = F9*(F2**Third-F1) / F4
      GammaI = F1 / (F2*F2**Third-F2)
      
      Call nwpwxc_EvFZet(tol_rho,GammaI,Zeta,FZeta,dfZdz,d2fZdz)
      Call nwpwxc_EvPWLC(eps0c(1),eps0c(2),eps0c(3),eps0c(4),eps0c(5),
     &    eps0c(6), RS,EU,dEUdRS,d2UdRS)
      Call nwpwxc_EvPWLC(eps1c(1),eps1c(2),eps1c(3),eps1c(4),eps1c(5),
     &    eps1c(6),RS,EP,dEPdRS,d2PdRS)
      Call nwpwxc_EvPWLC(epsc(1),epsc(2),epsc(3),epsc(4),epsc(5),
     &    epsc(6),RS,AlphaM,dAMdRS,d2AdRS)
      Z2 = Zeta*Zeta
      Z3 = Zeta*Z2
      Z4 = Zeta*Z3
      GZ = FZeta*Z4
      HZ = FZZI*(FZeta-GZ)
      PotLC = EU*(F1-GZ) + EP*GZ - AlphaM*HZ
      dLdS = dEUdRS*(F1-GZ) + dEPdRS*GZ - dAMdRS*HZ
      dGZ = dfZdz*Z4 + F4*FZeta*Z3
      dHZ = FZZI*(dFZdz-dGZ)
      dLdz = (EP-EU)*dGZ - AlphaM*dHZ
      d2GZ = d2fZdz*Z4 + F8*Z3*dfZdz + F12*FZeta*Z2
      d2HZ = FZZI*(d2FZdz-d2GZ)
      d2LdSS = d2UdRS*(F1-GZ) + d2PdRS*GZ - d2AdRS*HZ
      d2LdSZ = (dEPdRS-dEUdRS)*dGZ - dAMdRS*dHZ
      d2LdZZ = (EP-EU)*d2GZ - AlphaM*d2HZ
      Return
      End

c
c   f(zeta)
c
      Subroutine nwpwxc_EvFZet(tol_rho,S,Zeta,FZeta,dfZdz,d2fZdz)
      Implicit none
c
c     evaluate f(Zeta) and its derivatives for lsdac.
c
      double precision tol_rho
      double precision Small
      double precision S, Zeta, FZeta,dfZdz,d2fZdz
      double precision Zero, One, Two, Three, Four, Nine, F8, F27
      double precision OMZ, OPZ, OMZ2, OPZ2, OMZ3, OPZ3
      double precision F13, F43, F49, F827
      data Zero/0.0d0/, One/1.0d0/, Two/2.0d0/, Three/3.0d0/,
     $  Four/4.0d0/, Nine/9.0d0/, F8/8.0D0/, F27/27.0D0/
C
c      Small = 1.0d-14
      small=tol_rho
      FZeta = -Two
      dfZdz = Zero
      d2fZdz = Zero
      OMZ = One - Zeta
      OPZ = One + Zeta
      OMZ2 = OMZ**2
      OPZ2 = OPZ**2
      F13 = One / Three
      F43 = Four / Three
      F49 = Four / Nine
      F827 = F8 / F27
      If(OMZ.gt.Small) then
        OMZ3 = OMZ ** F13
        fZeta = fZeta + OMZ*OMZ3
        dfZdz = dfZdz - OMZ3
        d2fZdz = d2fZdz + OMZ3/OMZ
        endIf
      If(OPZ.gt.Small) then
        OPZ3 = OPZ ** F13
        fZeta = fZeta + OPZ*OPZ3
        dfZdz = dfZdz + OPZ3
        d2fZdz = d2fZdz + OPZ3/OPZ
        endIf
      fZeta = fZeta * S
      dfZdz = dfZdz * F43 * S
      d2fZdz = d2fZdz * F49 * S
      Return
      End

c
c  pw91 local correlation
c 
      Subroutine nwpwxc_EvPWLC(A,A1,B1,B2,B3,B4,RS,V,dVdRS,d2VdRS)
      Implicit none
C
C     Evaluate the interpolation function for PW91 local correlation.
C
      double precision A,A1,B1,B2,B3,B4,RS,V,dVdRS,d2VdRS
      double precision F1,F2, F3, F4
      double precision Q0, RS12, RS32,Q1,Q2
      double precision dQ0dRS,dQ1dRS,dQ2dRS
      double precision d2Q1dS, d2Q2dS
      data F1/1.0d0/, F2/2.0d0/, F3/3.0d0/, F4/4.0d0/
C
      Q0 = -F2*A*(F1+A1*RS)
      RS12 = Sqrt(RS)
      RS32 = RS*RS12
      Q1 = F2*A*(B1*RS12+B2*RS+B3*RS32+B4*RS*RS)
      Q2 = Log(F1+F1/Q1)
      V = Q0*Q2
C
      dQ0dRS = -F2*A*A1
      dQ1dRS = A*(B1/RS12+F2*B2+F3*B3*RS12+F4*B4*RS)
      dQ2dRS = -dQ1dRS/(Q1+Q1**2)
      dVdRS = dQ0dRS*Q2 + Q0*dQ2dRS
C
      d2Q1dS = A*(F3*B3/(RS12*F2)-B1/(RS32*F2)+F4*B4)
      d2Q2dS = (F2*Q1+F1)*(dQ1dRS/(Q1+Q1**2))**2 - d2Q1dS/(Q1+Q1**2)
      d2VdRS = F2*dQ0dRS*dQ2dRS + Q0*d2Q2dS
      Return
      End
C> @}
c $Id$
