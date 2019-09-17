#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
C> \ingroup nwxc
C> @{
C>
C> \file nwxc_c_b95.F
C> The B95 correlation functional
C>
C> @}
#endif
C>
C> \ingroup nwxc_priv
C> @{
C>
C> \brief Evaluate the B95 correlation functional
C>
C> Evaluate the B95 meta-GGA [1]. This routine is also used to 
C> evaluate the PW6B95 and PWB6K functionals [2].
C>
C> ### References ###
C>
C> [1] A.D. Becke,
C> "Density‐functional thermochemistry. IV. A new dynamical correlation
C> functional and implications for exact‐exchange mixing"
C> J. Chem. Phys. <b>104</b>, 1040-1046 (1996), DOI:
C> <a href="https://doi.org/10.1063/1.470829">
C> 10.1063/1.470829</a>.
C>
C> [2] Y. Zhao, D.G. Truhlar,
C> "Design of density functionals that are broadly accurate for
C> thermochemistry, thermochemical kinetics, and nonbonded
C> interactions", J. Phys. Chem. A <b>109</b> 5656-5667 (2005), DOI:
C> <a href="https://doi.org/10.1021/jp050536c">
C> 10.1021/jp050536c</a>.
C>
c    Bc95 correlation functional          
c           META GGA
C         utilizes ingredients:
c                              rho   -  density
c                              delrho - gradient of density
c                              tau (tauN)- K.S kinetic energy density
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
      Subroutine nwxc_c_b95_p(param, tol_rho, ipol, nq, wght, 
     &                        rho, rgamma, tau, func)
#else
      Subroutine nwxc_c_b95(param, tol_rho, ipol, nq, wght, 
     &                      rho, rgamma, tau, func)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      Subroutine nwxc_c_b95_d2(param, tol_rho, ipol, nq, wght, 
     &                         rho, rgamma, tau, func)
#else
      Subroutine nwxc_c_b95_d3(param, tol_rho, ipol, nq, wght, 
     &                         rho, rgamma, tau, func)
#endif
c
c$Id$
c
c  Reference
c    Becke, A. D. J. Chem. Phys. 1996, 104, 1040.
c
#include "nwad.fh"
c
      implicit none
#include "intf_nwxc_c_lsda.fh"
#include "intf_nwxc_c_b95css.fh"
c
#include "nwxc_param.fh"
c
#if defined(NWAD_PRINT)
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      type(nwad_dble)::param(2)
#else
      double precision param(2)
#endif
#else
      double precision param(2)!< [Input] Parameters of the functional
                               !< - param(1): \f$ C_{opp} \f$
                               !< - param(2): \f$ C_{\sigma\sigma} \f$
#endif
      double precision tol_rho !< [Input] The lower limit on the density
      integer nq               !< [Input] The number of points
      integer ipol             !< [Input] The number of spin channels
      double precision wght    !< [Input] The weight of the functional
c
c     Charge Density
c
      type(nwad_dble)::rho(nq,*) !< [Input] The density
c
c     Charge Density Gradient
c
      type(nwad_dble)::rgamma(nq,*) !< [Input] The norm of the density gradients
c
c     Kinetic Energy Density
c
      type(nwad_dble)::tau(nq,*) !< [Input] The kinetic energy density
c      
c     The functional
c
      type(nwad_dble)::func(*)  !< [Output] The value of the functional
c
c     Sampling Matrices for the XC Potential & Energy
c
c     double precision Amat(nq,*) !< [Output] The derivative wrt rho
c     double precision Cmat(nq,*) !< [Output] The derivative wrt rgamma
c     double precision Mmat(nq,*) !< [Output] The derivative wrt tau
c
c     Threshold parameters
c
      double precision DTol,F1, F2, F3, F4
#if defined(NWAD_PRINT)
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      type(nwad_dble) :: COpp
#else
      double precision COpp 
#endif
#else
      double precision COpp 
#endif
      Data F1/1.0d0/,F2/2.0d0/,
     &     F3/3.0d0/,F4/4.0d0/ 

      integer n

c    call to the bc95css subroutine
      type(nwad_dble)::PA,GAA,TA,FA,EUA,ChiA
      type(nwad_dble)::PB,GBB,TB,FB,EUB,ChiB
      double precision EUPA,ChiAP,ChiAG,EUPB,ChiBP,ChiBG
      double precision FPA,FGA,FTA
      double precision FPB,FGB,FTB
c
      type(nwad_dble)::P,RS,RSP,Zeta,EUEG,Denom,PotLC
      double precision  sop
      double precision Pi, F6, F43, Pi34, F13, 
     &dZdA,dZdB,dLdS,dLdZ,d2LdSS,d2LdSZ,d2LdZZ
      double precision DenPA, DenPB, DenGA, DenGB
      double precision EUEGPA,EUEGPB

      
c
c     ======> BOTH SPIN-RESTRICETED AND UNRESTRICTED <======
c
      DTol=tol_rho
      sop=1.0d0
      Copp = param(1)
c     if (ijmswitch.eq.1) then
c     Parameters for PW6B95 Correlation
c       COpp=0.00262d0 
c     elseif (ijmswitch.eq.2) then
c     Parameters for PWB6K Correlation
c       COpp=0.00353d0
c     endif
      Pi = F4*ATan(F1)
      F6=6.0d0
      F43 = F4 / F3
      Pi34 = F3 / (F4*Pi)
      F13 = F1 / F3

      do 20 n = 1, nq
       if (ipol.eq.1) then
         if (rho(n,R_T).lt.DTol) goto 20
       else
         if (rho(n,R_A)+rho(n,R_B).lt.DTol) goto 20
       endif
       if (ipol.eq.1) then
c
c    get the density, gradient, and tau for the alpha spin from the total 
c
         PA = rho(n,R_T)/F2
c        GAA = (    delrho(n,1,1)*delrho(n,1,1) +
c    &                 delrho(n,2,1)*delrho(n,2,1) +
c    &                 delrho(n,3,1)*delrho(n,3,1))/4.0d0
         GAA = rgamma(n,G_TT)/F4
c  In the bc95css subroutine, we use 2*TA as the tau, so we do not divide 
c  the tau by 2 here

         TA = tau(n,T_T)
!         TA=0.0005d0
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)                  
#if defined(NWAD_PRINT)
         Call nwxc_c_b95ss_p(param,dtol,PA,GAA,TA,FA,FPA,FGA,FTA,EUA,
     &                       ChiA,EUPA,ChiAP,ChiAG)
#else
         Call nwxc_c_b95ss(param,dtol,PA,GAA,TA,FA,FPA,FGA,FTA,EUA,
     &                     ChiA,EUPA,ChiAP,ChiAG)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)                  
         Call nwxc_c_b95ss_d2(param,dtol,PA,GAA,TA,FA,FPA,FGA,FTA,EUA,
     &                        ChiA,EUPA,ChiAP,ChiAG)
#else
         Call nwxc_c_b95ss_d3(param,dtol,PA,GAA,TA,FA,FPA,FGA,FTA,EUA,
     &                        ChiA,EUPA,ChiAP,ChiAG)
#endif
         PB = PA
         GBB = GAA
         TB = TA
         FB = FA
         FPB = FPA
         FGB = FGA
         FTB = FTA
         EUB = EUA
         ChiB = ChiA
         EUPB = EUPA
         ChiBP = ChiAP
         ChiBG = ChiAG

         func(n)=func(n)+ 2.0d0*FA*wght
c        Amat(n,D1_RA)=Amat(n,D1_RA)+ FPA*wght
c        Cmat(n,D1_GAA)=  Cmat(n,D1_GAA) + FGA*wght
c        Mmat(n,D1_TA)=  Mmat(n,D1_TA) + FTA*wght
#if 0
      write (0,'(A,3F20.6)') " Amat Cmat Mmat",FPA,FGA,FTA
#endif
 
 
c UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUnrestricted
      else  ! ipol=2
c
c        ======> SPIN-UNRESTRICTED <======
c
c
c       alpha
c
         
         PA = rho(n,R_A)
         if (PA.le.0.5d0*DTol) go to 25
c        GAA =   delrho(n,1,1)*delrho(n,1,1) +
c    &           delrho(n,2,1)*delrho(n,2,1) +
c    &          delrho(n,3,1)*delrho(n,3,1)
         GAA = rgamma(n,G_AA)
c
c  In the bc95css subroutine, we use 2*TA as the tau 
c
         TA = tau(n,T_A)*2.0d0

#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)                  
#if defined(NWAD_PRINT)
         Call nwxc_c_b95ss_p(param,dtol,PA,GAA,TA,FA,FPA,FGA,FTA,EUA,
     &                       ChiA,EUPA,ChiAP,ChiAG)
#else
         Call nwxc_c_b95ss(param,dtol,PA,GAA,TA,FA,FPA,FGA,FTA,EUA,
     &                     ChiA,EUPA,ChiAP,ChiAG)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)                  
         Call nwxc_c_b95ss_d2(param,dtol,PA,GAA,TA,FA,FPA,FGA,FTA,EUA,
     &                        ChiA,EUPA,ChiAP,ChiAG)
#else
         Call nwxc_c_b95ss_d3(param,dtol,PA,GAA,TA,FA,FPA,FGA,FTA,EUA,
     &                        ChiA,EUPA,ChiAP,ChiAG)
#endif
         func(n)=func(n)+ FA*wght
c        Amat(n,D1_RA)=Amat(n,D1_RA)+ FPA*wght
c        Cmat(n,D1_GAA)=  Cmat(n,D1_GAA) + FGA*wght
c      2*0.5=1.0 for Mmat
c        Mmat(n,D1_TA)=  Mmat(n,D1_TA) + FTA*wght
#if 0
      write (0,'(A,3F20.6)') "AAmat Cmat Mmat",FPA,FGA,FTA
#endif
c
c  In the bc95css subroutine, we use 2*TA as the tau, 
c
c
c       Beta 
c
 25       continue
         PB = rho(n,R_B)
         if(PB.le.0.5d0*DTol) go to 30
c        GBB =   delrho(n,1,2)*delrho(n,1,2) +
c    &           delrho(n,2,2)*delrho(n,2,2) +
c    &          delrho(n,3,2)*delrho(n,3,2)
         GBB = rgamma(n,G_BB)

         TB = tau(n,T_B)*2.0d0

#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)                  
#if defined(NWAD_PRINT)
         Call nwxc_c_b95ss_p(param,dtol,PB,GBB,TB,FB,FPB,FGB,FTB,EUB,
     &                       ChiB,EUPB,ChiBP,ChiBG)
#else
         Call nwxc_c_b95ss(param,dtol,PB,GBB,TB,FB,FPB,FGB,FTB,EUB,
     &                     ChiB,EUPB,ChiBP,ChiBG)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)                  
         Call nwxc_c_b95ss_d2(param,dtol,PB,GBB,TB,FB,FPB,FGB,FTB,EUB,
     &                        ChiB,EUPB,ChiBP,ChiBG)
#else
         Call nwxc_c_b95ss_d3(param,dtol,PB,GBB,TB,FB,FPB,FGB,FTB,EUB,
     &                        ChiB,EUPB,ChiBP,ChiBG)
#endif
         func(n)=func(n)+ FB*wght
c        Amat(n,D1_RB)= Amat(n,D1_RB)+ FPB*wght
c        Cmat(n,D1_GBB)=  Cmat(n,D1_GBB) + FGB*wght
c        Mmat(n,D1_TB)=  Mmat(n,D1_TB) + FTB*wght
#if 0
      write (0,'(A,3F20.6)') "BAmat Cmat Mmat",FPB,FGB,FTB
#endif
      endif
 30   continue
      P = PA+PB
      If(PA.gt.0.5d0*DTol.and.PB.gt.0.5d0*DTol) then
          RS = (Pi34/P) ** F13 
          RSP = -RS/(F3*P)
          Zeta = (PA-PB)/P
c         dZdA = (F1-Zeta)/P
c         dZdB = (-F1-Zeta)/P
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)                  
#if defined(NWAD_PRINT)
          Call nwxc_c_lsda_p(dtol,
     D         RS,Zeta,PotLC,dLdS,dLdZ,d2LdSS,d2LdSZ,
     $         d2LdZZ)
#else
          Call nwxc_c_lsda(dtol,
     D         RS,Zeta,PotLC,dLdS,dLdZ,d2LdSS,d2LdSZ,
     $         d2LdZZ)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)                  
          Call nwxc_c_lsda_d2(dtol,
     D         RS,Zeta,PotLC,dLdS,dLdZ,d2LdSS,d2LdSZ,
     $         d2LdZZ)
#else
          Call nwxc_c_lsda_d3(dtol,
     D         RS,Zeta,PotLC,dLdS,dLdZ,d2LdSS,d2LdSZ,
     $         d2LdZZ)
#endif
          EUEG = P*PotLC - EUA - EUB
          Denom = F1 + COpp*(ChiA+ChiB)
          func(n)=func(n)+ sop*EUEG/Denom*wght
c         DenPA = COpp*ChiAP
c         DenPB = COpp*ChiBP
c         DenGA = COpp*ChiAG
c         DenGB = COpp*ChiBG
c         EUEGPA = PotLC + P*dLdS*RSP + P*dLdZ*dZdA - EUPA
c         EUEGPB = PotLC + P*dLdS*RSP + P*dLdZ*dZdB - EUPB
c         if (ipol.eq.1) then 
c           Amat(n,D1_RA)  = Amat(n,D1_RA) + 
c    &               sop*(EUEGPA/Denom - EUEG*DenPA/Denom**2)*wght
c           Cmat(n,D1_GAA) =  Cmat(n,D1_GAA) 
c    &                     - sop*(EUEG*DenGA/Denom**2)*wght
c         else
c           Amat(n,D1_RA) = Amat(n,D1_RA) + 
c    &                 sop*(EUEGPA/Denom - EUEG*DenPA/Denom**2) *wght
c           Amat(n,D1_RB) = Amat(n,D1_RB) +
c    &                 sop*(EUEGPB/Denom - EUEG*DenPB/Denom**2) *wght
c           Cmat(n,D1_GAA) = Cmat(n,D1_GAA)
c    &                     - sop*EUEG*DenGA/Denom**2*wght
c           Cmat(n,D1_GBB) = Cmat(n,D1_GBB)
c    &                     - sop*EUEG*DenGB/Denom**2*wght
c         endif
      endIf
c      write (*,*) "Amat(n,1),Cmat(n,1),Mmat(n,1)",Amat(n,1),Cmat(n,1)
c     & ,Mmat(n,1)
c      stop
20    continue
      end

#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
      Subroutine nwxc_c_b95ss_p(param,dtol,PX,GX,TX,F,FP,FG,FT,EUEG,
     &                          Chi,EUEGP,ChiP,ChiG)
#else
      Subroutine nwxc_c_b95ss(param,dtol,PX,GX,TX,F,FP,FG,FT,EUEG,
     &                        Chi,EUEGP,ChiP,ChiG)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      Subroutine nwxc_c_b95ss_d2(param,dtol,PX,GX,TX,F,FP,FG,FT,EUEG,
     &                        Chi,EUEGP,ChiP,ChiG)
#else
      Subroutine nwxc_c_b95ss_d3(param,dtol,PX,GX,TX,F,FP,FG,FT,EUEG,
     &                        Chi,EUEGP,ChiP,ChiG)
#endif
c
#include "nwad.fh"
c
      Implicit none
#include "intf_nwxc_c_lsda.fh"
C
C     Compute the same-spin part of the bc95 correlation functional for one grid
C     point and one spin-case.
C
C
c     integer ijmswitch
#if defined(NWAD_PRINT)
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      type(nwad_dble)::param(2)
#else
      double precision param(2)
#endif
#else
      double precision param(2)!< [Input] Parameters of the functional
                               !< - param(1): \f$ C_{opp} \f$
                               !< - param(2): \f$ C_{\sigma\sigma} \f$
#endif
      double precision dtol
      type(nwad_dble)::PX, GX, TX, F, E
      type(nwad_dble)::D, RS, RSP, EUEG, DUEG, Chi, Denom
      type(nwad_dble)::Zeta,PotLC
      double precision FP, FG, FT
      double precision EUEGP, ChiP, ChiG
      double precision Zero, Pt25, F1, F2, F3, F4, F5, F6, F8, F11
      double precision Pi, Pi34, F13, F23, F43, F53, F83, F113
      double precision FDUEG
      double precision DenomG, DenomP, DUEGP, DP, DG, DT  
      double precision d2LdSS,d2LdSZ,d2LdZZ,dLdS,dLdZ
#if defined(NWAD_PRINT)
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      type(nwad_dble)::Css
#else
      double precision Css
#endif
#else
      double precision Css
#endif


      Data Zero/0.0d0/, Pt25/0.25d0/, F1/1.0d0/, F2/2.0d0/, F3/3.0d0/,
     $  F4/4.0d0/, F5/5.0d0/, F6/6.0d0/, F8/8.0d0/, F11/11.0d0/
C
c     if (ijmswitch.eq.1) then
C     Parameters for PW6B95 Correlation
c      Css=0.03668d0
c     elseif (ijmswitch.eq.2) then
C     Parameters for PWB6K Correlation
c      Css=0.04120d0
c     endif
      Css = param(2)
c     DTol =1.0d-6 
      If(PX.le.DTol) then
        EUEG = Zero
        Chi = Zero
        EUEGP = Zero
        ChiP = Zero
        ChiG = Zero
        PX = Zero
        GX = Zero 
        TX = Zero
        F  = Zero
        FP = Zero
        FG = Zero
        FT = Zero
      else
        Pi = F4*ATan(F1)
        Pi34 = F3 / (F4*Pi)
        F13 = F1 / F3
        F23 = F2 / F3
        F43 = F2 * F23
        F53 = F5 / F3
        F83 = F8 / F3
        F113 = F11 / F3
        FDUEG = (F3/F5)*(F6*Pi*Pi)**F23
        RS = (Pi34/PX) ** F13
        Zeta = F1
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
        Call nwxc_c_lsda_p(dtol,
     D       RS,Zeta,PotLC,dLdS,dLdZ,d2LdSS,d2LdSZ,d2LdZZ)
#else
        Call nwxc_c_lsda(dtol,
     D       RS,Zeta,PotLC,dLdS,dLdZ,d2LdSS,d2LdSZ,d2LdZZ)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
        Call nwxc_c_lsda_d2(dtol,
     D       RS,Zeta,PotLC,dLdS,dLdZ,d2LdSS,d2LdSZ,d2LdZZ)
#else
        Call nwxc_c_lsda_d3(dtol,
     D       RS,Zeta,PotLC,dLdS,dLdZ,d2LdSS,d2LdSZ,d2LdZZ)
#endif
        EUEG = PX*PotLC
        D = TX - Pt25*GX/PX
        DUEG = FDUEG*PX**F53
        Chi = GX/PX**F83
        Denom = F1 + Css*Chi
        E = D*EUEG/(DUEG*Denom*Denom)
c        write (*,*) "ijmswitch, Css, E= ",ijmswitch, Css, E
c        stop
        F = E 
c
c       RSP = -RS/(F3*Px)
c       ChiG = F1/PX**F83
c       ChiP = -F83*Chi/PX
c       DenomG = Css*ChiG
c       DenomP = Css*ChiP
c       DUEGP = F53*DUEG/PX
c       DP = Pt25*GX/PX**2
c       DG = -Pt25/PX
c       DT = F1
c       EUEGP = PotLC + PX*dLdS*RSP
c       FP = DP*EUEG/(DUEG*Denom*Denom) +
c    $      D*EUEGP/(DUEG*Denom*Denom)
c    $      - D*EUEG*DUEGP/(DUEG*Denom)**2 -
c    $      F2*D*EUEG*DenomP/(DUEG*Denom*Denom*Denom)
c       FG =DG*EUEG/(DUEG*Denom*Denom) -
c    $      F2*D*EUEG*DenomG/(DUEG*Denom*Denom*Denom)
c       FT =DT*EUEG/(DUEG*Denom*Denom)
       Endif
       Return
       End

#ifndef NWAD_PRINT
#define NWAD_PRINT
c
c     Compile source again for the 2nd derivative case
c
#include "nwxc_c_b95.F"
#endif
#ifndef SECOND_DERIV
#define SECOND_DERIV
c
c     Compile source again for the 2nd derivative case
c
#include "nwxc_c_b95.F"
#endif
#ifndef THIRD_DERIV
#define THIRD_DERIV
c
c     Compile source again for the 3rd derivative case
c
#include "nwxc_c_b95.F"
#endif
#undef NWAD_PRINT
C>
C> @}
