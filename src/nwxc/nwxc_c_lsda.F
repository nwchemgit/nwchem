#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
C> \ingroup nwxc
C> @{
C> 
C> \file nwxc_c_lsda.F
C> Another implementation of the PW91 LDA correlation functional
C>
C> @}
#endif
C>
C> \ingroup nwxc_priv
C> @{
C>
C> \brief An alternative implementation of the functional 
C> `nwxc_c_pw91lda`
C>
C> This stems from the original authors of the functional and hence
C> it follows different calling conventions from the other routines
C> in the library. For this reason it is also only meant to be used
C> as an internal component of PW91 related functionals.
C>
c
c Perdew 91 local correlation functional at one grid point 
c
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
      Subroutine nwxc_c_lsda_p(tol_rho,
     R     rs,zeta,PotLC,dLdS,dLdZ,d2LdSS,d2LdSZ,d2LdZZ)
#else
      Subroutine nwxc_c_lsda(tol_rho,
     R     rs,zeta,PotLC,dLdS,dLdZ,d2LdSS,d2LdSZ,d2LdZZ)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      Subroutine nwxc_c_lsda_d2(tol_rho,
     R     rs,zeta,PotLC,dLdS,dLdZ,d2LdSS,d2LdSZ,d2LdZZ)
#else
      Subroutine nwxc_c_lsda_d3(tol_rho,
     R     rs,zeta,PotLC,dLdS,dLdZ,d2LdSS,d2LdSZ,d2LdZZ)
#endif

C     Evaluate the Perdew 91 local correlation functional and its
C     derivatives at one point. 
c
#include "nwad.fh"
c
      implicit none
c
#include "intf_nwxc_evfzet.fh"
#include "intf_nwxc_evpwlc.fh"
c
      double precision tol_rho
      type(nwad_dble)::rs
      type(nwad_dble)::PotLC
      double precision dLdS, dLdZ, d2dSS, d2dSZ, d2dZZ
      double precision eps0c(6), eps1c(6), epsc(6)
      double precision F1, F2, F3, F4, F6, F8, F9, F12, F24, F36
      type(nwad_dble)::Zeta,FZeta,z4,GZ,HZ
      double precision GammaI,dfZdz,d2fZdz 
      type(nwad_dble)::EU,EP
      double precision dEUdRS,d2UdRS
      double precision dEPdRS,d2PdRS
      type(nwad_dble)::AlphaM
      double precision dAMdRS,d2AdRS
      double precision dGZ, dHZ, d2GZ, d2HZ
      double precision d2LdSS,d2LdSZ,d2LdZZ,fzzi,
     Z     third,pi
 
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
      
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
      Call nwxc_EvFZet_p(tol_rho,GammaI,Zeta,FZeta,dfZdz,d2fZdz)
      Call nwxc_EvPWLC_p(eps0c(1),eps0c(2),eps0c(3),eps0c(4),eps0c(5),
     &    eps0c(6),RS,EU,dEUdRS,d2UdRS)
      Call nwxc_EvPWLC_p(eps1c(1),eps1c(2),eps1c(3),eps1c(4),eps1c(5),
     &    eps1c(6),RS,EP,dEPdRS,d2PdRS)
      Call nwxc_EvPWLC_p(epsc(1),epsc(2),epsc(3),epsc(4),epsc(5),
     &    epsc(6),RS,AlphaM,dAMdRS,d2AdRS)
#else
      Call nwxc_EvFZet(tol_rho,GammaI,Zeta,FZeta,dfZdz,d2fZdz)
      Call nwxc_EvPWLC(eps0c(1),eps0c(2),eps0c(3),eps0c(4),eps0c(5),
     &    eps0c(6),RS,EU,dEUdRS,d2UdRS)
      Call nwxc_EvPWLC(eps1c(1),eps1c(2),eps1c(3),eps1c(4),eps1c(5),
     &    eps1c(6),RS,EP,dEPdRS,d2PdRS)
      Call nwxc_EvPWLC(epsc(1),epsc(2),epsc(3),epsc(4),epsc(5),
     &    epsc(6),RS,AlphaM,dAMdRS,d2AdRS)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      Call nwxc_EvFZet_d2(tol_rho,GammaI,Zeta,FZeta,dfZdz,d2fZdz)
      Call nwxc_EvPWLC_d2(eps0c(1),eps0c(2),eps0c(3),eps0c(4),eps0c(5),
     &    eps0c(6),RS,EU,dEUdRS,d2UdRS)
      Call nwxc_EvPWLC_d2(eps1c(1),eps1c(2),eps1c(3),eps1c(4),eps1c(5),
     &    eps1c(6),RS,EP,dEPdRS,d2PdRS)
      Call nwxc_EvPWLC_d2(epsc(1),epsc(2),epsc(3),epsc(4),epsc(5),
     &    epsc(6),RS,AlphaM,dAMdRS,d2AdRS)
#else
      Call nwxc_EvFZet_d3(tol_rho,GammaI,Zeta,FZeta,dfZdz,d2fZdz)
      Call nwxc_EvPWLC_d3(eps0c(1),eps0c(2),eps0c(3),eps0c(4),eps0c(5),
     &    eps0c(6),RS,EU,dEUdRS,d2UdRS)
      Call nwxc_EvPWLC_d3(eps1c(1),eps1c(2),eps1c(3),eps1c(4),eps1c(5),
     &    eps1c(6),RS,EP,dEPdRS,d2PdRS)
      Call nwxc_EvPWLC_d3(epsc(1),epsc(2),epsc(3),epsc(4),epsc(5),
     &    epsc(6),RS,AlphaM,dAMdRS,d2AdRS)
#endif
c     Z2 = Zeta*Zeta
c     Z3 = Zeta*Z2
      Z4 = Zeta**4
      GZ = FZeta*Z4
      HZ = FZZI*(FZeta-GZ)
      PotLC = EU*(F1-GZ) + EP*GZ - AlphaM*HZ
c     dLdS = dEUdRS*(F1-GZ) + dEPdRS*GZ - dAMdRS*HZ
c     dGZ = dfZdz*Z4 + F4*FZeta*Z3
c     dHZ = FZZI*(dFZdz-dGZ)
c     dLdz = (EP-EU)*dGZ - AlphaM*dHZ
c     d2GZ = d2fZdz*Z4 + F8*Z3*dfZdz + F12*FZeta*Z2
c     d2HZ = FZZI*(d2FZdz-d2GZ)
c     d2LdSS = d2UdRS*(F1-GZ) + d2PdRS*GZ - d2AdRS*HZ
c     d2LdSZ = (dEPdRS-dEUdRS)*dGZ - dAMdRS*dHZ
c     d2LdZZ = (EP-EU)*d2GZ - AlphaM*d2HZ
      Return
      End

c
c   f(zeta)
c
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
      Subroutine nwxc_EvFZet_p(tol_rho,S,Zeta,FZeta,dfZdz,d2fZdz)
#else
      Subroutine nwxc_EvFZet(tol_rho,S,Zeta,FZeta,dfZdz,d2fZdz)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      Subroutine nwxc_EvFZet_d2(tol_rho,S,Zeta,FZeta,dfZdz,d2fZdz)
#else
      Subroutine nwxc_EvFZet_d3(tol_rho,S,Zeta,FZeta,dfZdz,d2fZdz)
#endif
#include "nwad.fh"
      Implicit none
c
c     evaluate f(Zeta) and its derivatives for lsdac.
c
      double precision tol_rho
      double precision Small
      type(nwad_dble)::Zeta, FZeta
      double precision S, dfZdz,d2fZdz
      double precision Zero, One, Two, Three, Four, Nine, F8, F27
      type(nwad_dble)::OMZ, OPZ
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
c     OMZ2 = OMZ**2.0d0
c     OPZ2 = OPZ**2.0d0
      F13 = One / Three
      F43 = Four / Three
      F49 = Four / Nine
      F827 = F8 / F27
      If(OMZ.gt.Small) then
c       OMZ3 = OMZ ** F13
        fZeta = fZeta + OMZ**F43
c       dfZdz = dfZdz - OMZ3
c       d2fZdz = d2fZdz + OMZ3/OMZ
      endIf
      If(OPZ.gt.Small) then
c       OPZ3 = OPZ ** F13
        fZeta = fZeta + OPZ**F43
c       dfZdz = dfZdz + OPZ3
c       d2fZdz = d2fZdz + OPZ3/OPZ
      endIf
      fZeta = fZeta * S
c     dfZdz = dfZdz * F43 * S
c     d2fZdz = d2fZdz * F49 * S
      Return
      End

c
c  pw91 local correlation
c 
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
      Subroutine nwxc_EvPWLC_p(A,A1,B1,B2,B3,B4,RS,V,dVdRS,d2VdRS)
#else
      Subroutine nwxc_EvPWLC(A,A1,B1,B2,B3,B4,RS,V,dVdRS,d2VdRS)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      Subroutine nwxc_EvPWLC_d2(A,A1,B1,B2,B3,B4,RS,V,dVdRS,d2VdRS)
#else
      Subroutine nwxc_EvPWLC_d3(A,A1,B1,B2,B3,B4,RS,V,dVdRS,d2VdRS)
#endif
#include "nwad.fh"
      Implicit none
C
C     Evaluate the interpolation function for PW91 local correlation.
C
      type(nwad_dble)::RS,V,RS12,RS32
      double precision A,A1,B1,B2,B3,B4,dVdRS,d2VdRS
      double precision F1,F2, F3, F4
      type(nwad_dble)::Q0,Q1,Q2
      double precision dQ0dRS,dQ1dRS,dQ2dRS
      double precision d2Q1dS,d2Q2dS
      data F1/1.0d0/, F2/2.0d0/, F3/3.0d0/, F4/4.0d0/
C
      Q0 = -F2*A*(F1+A1*RS)
      RS12 = Sqrt(RS)
      RS32 = RS*RS12
      Q1 = F2*A*(B1*RS12+B2*RS+B3*RS32+B4*RS**2)
      Q2 = Log(F1+F1/Q1)
      V = Q0*Q2
C
c     dQ0dRS = -F2*A*A1
c     dQ1dRS = A*(B1/RS12+F2*B2+F3*B3*RS12+F4*B4*RS)
c     dQ2dRS = -dQ1dRS/(Q1+Q1**2)
c     dVdRS = dQ0dRS*Q2 + Q0*dQ2dRS
C
c     d2Q1dS = A*(F3*B3/(RS12*F2)-B1/(RS32*F2)+F4*B4)
c     d2Q2dS = (F2*Q1+F1)*(dQ1dRS/(Q1+Q1**2))**2 - d2Q1dS/(Q1+Q1**2)
c     d2VdRS = F2*dQ0dRS*dQ2dRS + Q0*d2Q2dS
      Return
      End
#ifndef NWAD_PRINT
#define NWAD_PRINT
c
c     Compile source again for the 2nd derivative case
c
#include "nwxc_c_lsda.F"
#endif
#ifndef SECOND_DERIV
#define SECOND_DERIV
c
c     Compile source again for the 2nd derivative case
c
#include "nwxc_c_lsda.F"
#endif
#ifndef THIRD_DERIV
#define THIRD_DERIV
c
c     Compile source again for the 3rd derivative case
c
#include "nwxc_c_lsda.F"
#endif
#undef NWAD_PRINT
C> @}
c $Id$
