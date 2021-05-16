C> \ingroup nwpwxc
C> @{
C>
C> \file nwpwxc_c_b95.F
C> The B95 correlation functional
C>
C> @}
C>
C> \ingroup nwpwxc_priv
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

      Subroutine nwpwxc_c_b95(param, tol_rho, ipol, nq, wght, 
     &                      rho, rgamma, tau, func, Amat, Cmat, Mmat)
c
c$Id$
c
c  Reference
c    Becke, A. D. J. Chem. Phys. 1996, 104, 1040.
c
      implicit none
c
#include "nwpwxc_param.fh"
c
      double precision param(2)!< [Input] Parameters of the functional
                               !< - param(1): \f$ C_{opp} \f$
                               !< - param(2): \f$ C_{\sigma\sigma} \f$
      double precision tol_rho !< [Input] The lower limit on the density
      integer nq               !< [Input] The number of points
      integer ipol             !< [Input] The number of spin channels
      double precision wght    !< [Input] The weight of the functional
c
c     Charge Density
c
      double precision rho(nq,*) !< [Input] The density
c
c     Charge Density Gradient
c
      double precision rgamma(nq,*) !< [Input] The norm of the density gradients
c
c     Kinetic Energy Density
c
      double precision tau(nq,*) !< [Input] The kinetic energy density
c      
c     The functional
c
      double precision func(*)  !< [Output] The value of the functional
c
c     Sampling Matrices for the XC Potential & Energy
c
      double precision Amat(nq,*) !< [Output] The derivative wrt rho
      double precision Cmat(nq,*) !< [Output] The derivative wrt rgamma
      double precision Mmat(nq,*) !< [Output] The derivative wrt tau
c
c     Threshold parameters
c
      double precision DTol,F1, F2, F3, F4,COpp 
      Data COpp/0.0031d0/,F1/1.0d0/,F2/2.0d0/,
     & F3/3.0d0/,F4/4.0d0/ 

      integer n

c    call to the bc95css subroutine
      double precision PA,GAA,TA,FA,FPA,FGA,FTA,EUA,ChiA,EUPA
     &,ChiAP,ChiAG
      double precision PB,GBB,TB,FB,FPB,FGB,FTB,EUB,ChiB,EUPB
     &,ChiBP,ChiBG
c
      double precision  sop
      double precision Pi, F6, F43, Pi34, F13, 
     &RS,RSP,Zeta,dZdA,dZdB,PotLC,dLdS,dLdZ,d2LdSS,d2LdSZ,d2LdZZ
      double precision P, EUEG,Denom, DenPA, DenPB, DenGA, DenGB
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
         GAA = rgamma(n,G_TT)/4.0d0
c  In the bc95css subroutine, we use 2*TA as the tau, so we do not divide 
c  the tau by 2 here

         TA = tau(n,T_T) 
!         TA=0.0005d0
                  
         Call nwpwxc_c_b95ss(param,dtol,PA,GAA,TA,FA,FPA,FGA,FTA,EUA,
     &                     ChiA,EUPA,ChiAP,ChiAG)
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
         Amat(n,D1_RA)=Amat(n,D1_RA)+ FPA*wght
         Cmat(n,D1_GAA)=  Cmat(n,D1_GAA) + FGA*wght
         Mmat(n,D1_TA)=  Mmat(n,D1_TA) + FTA*wght
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

         Call nwpwxc_c_b95ss(param,dtol,PA,GAA,TA,FA,FPA,FGA,FTA,EUA,
     &                     ChiA,EUPA,ChiAP,ChiAG)
         func(n)=func(n)+ FA*wght
         Amat(n,D1_RA)=Amat(n,D1_RA)+ FPA*wght
         Cmat(n,D1_GAA)=  Cmat(n,D1_GAA) + FGA*wght
c      2*0.5=1.0 for Mmat
         Mmat(n,D1_TA)=  Mmat(n,D1_TA) + FTA*wght
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

         Call nwpwxc_c_b95ss(param,dtol,PB,GBB,TB,FB,FPB,FGB,FTB,EUB,
     &                     ChiB,EUPB,ChiBP,ChiBG)
         func(n)=func(n)+ FB*wght
         Amat(n,D1_RB)= Amat(n,D1_RB)+ FPB*wght
         Cmat(n,D1_GBB)=  Cmat(n,D1_GBB) + FGB*wght
         Mmat(n,D1_TB)=  Mmat(n,D1_TB) + FTB*wght
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
          dZdA = (F1-Zeta)/P
          dZdB = (-F1-Zeta)/P
          Call nwpwxc_c_lsda(dtol,
     D         RS,Zeta,PotLC,dLdS,dLdZ,d2LdSS,d2LdSZ,
     $      d2LdZZ)
          EUEG = P*PotLC - EUA - EUB
          Denom = F1 + COpp*(ChiA+ChiB)
          func(n)=func(n)+ sop*EUEG/Denom*wght
          DenPA = COpp*ChiAP
          DenPB = COpp*ChiBP
          DenGA = COpp*ChiAG
          DenGB = COpp*ChiBG
          EUEGPA = PotLC + P*dLdS*RSP + P*dLdZ*dZdA - EUPA
          EUEGPB = PotLC + P*dLdS*RSP + P*dLdZ*dZdB - EUPB
          if (ipol.eq.1) then 
            Amat(n,D1_RA)  = Amat(n,D1_RA) + 
     &               sop*(EUEGPA/Denom - EUEG*DenPA/Denom**2)*wght
            Cmat(n,D1_GAA) =  Cmat(n,D1_GAA) 
     &                     - sop*(EUEG*DenGA/Denom**2)*wght
          else
            Amat(n,D1_RA) = Amat(n,D1_RA) + 
     &                 sop*(EUEGPA/Denom - EUEG*DenPA/Denom**2) *wght
            Amat(n,D1_RB) = Amat(n,D1_RB) +
     &                 sop*(EUEGPB/Denom - EUEG*DenPB/Denom**2) *wght
            Cmat(n,D1_GAA) = Cmat(n,D1_GAA)
     &                     - sop*EUEG*DenGA/Denom**2*wght
            Cmat(n,D1_GBB) = Cmat(n,D1_GBB)
     &                     - sop*EUEG*DenGB/Denom**2*wght
          endif
      endIf
c      write (*,*) "Amat(n,1),Cmat(n,1),Mmat(n,1)",Amat(n,1),Cmat(n,1)
c     & ,Mmat(n,1)
c      stop
20    continue
      end

      Subroutine nwpwxc_c_b95_d2()
      call errquit(' bc95: d2 not coded ',0,0)
      return
      end




      Subroutine nwpwxc_c_b95ss(param,dtol,PX,GX,TX,F,FP,FG,FT,EUEG,
     &                        Chi,EUEGP,ChiP,ChiG)
      Implicit none
C
C     Compute the same-spin part of the bc95 correlation functional for one grid
C     point and one spin-case.
C
C
c     integer ijmswitch
      double precision param(2)!< [Input] Parameters of the functional
                               !< - param(1): \f$ C_{opp} \f$
                               !< - param(2): \f$ C_{\sigma\sigma} \f$
      double precision PX, GX, TX, F, FP, FG, FT, DTol
      double precision EUEG, Chi, EUEGP, ChiP, ChiG,Css
      double precision Zero, Pt25, F1, F2, F3, F4, F5, F6, F8, F11
      double precision Pi, Pi34, F13, F23, F43, F53, F83, F113
      double precision RS, FDUEG, D,  RSP,DUEG, Denom, PotLC
      double precision E, DenomG, DenomP, DUEGP, DP, DG, DT  
      double precision d2LdSS,d2LdSZ,d2LdZZ,dLdS,dLdZ
     


      Data Zero/0.0d0/, Pt25/0.25d0/, F1/1.0d0/, F2/2.0d0/, F3/3.0d0/,
     $  F4/4.0d0/, F5/5.0d0/, F6/6.0d0/, F8/8.0d0/, F11/11.0d0/,
     $  Css/0.038d0/
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
        Call nwpwxc_c_lsda(dtol,
     D       RS,F1,PotLC,dLdS,dLdZ,d2LdSS,d2LdSZ,d2LdZZ)
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
        RSP = -RS/(F3*Px)
        ChiG = F1/PX**F83
        ChiP = -F83*Chi/PX
        DenomG = Css*ChiG
        DenomP = Css*ChiP
        DUEGP = F53*DUEG/PX
        DP = Pt25*GX/PX**2
        DG = -Pt25/PX
        DT = F1
        EUEGP = PotLC + PX*dLdS*RSP
        FP = DP*EUEG/(DUEG*Denom*Denom) +
     $      D*EUEGP/(DUEG*Denom*Denom)
     $      - D*EUEG*DUEGP/(DUEG*Denom)**2 -
     $      F2*D*EUEG*DenomP/(DUEG*Denom*Denom*Denom)
        FG =DG*EUEG/(DUEG*Denom*Denom) -
     $      F2*D*EUEG*DenomG/(DUEG*Denom*Denom*Denom)
        FT =DT*EUEG/(DUEG*Denom*Denom)
       Endif
       Return
       End
C>
C> @}
