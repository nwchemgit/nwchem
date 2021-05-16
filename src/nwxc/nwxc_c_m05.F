#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
C> \ingroup nwxc
C> @{
C>
C> \file nwxc_c_m05.F
C> Implementation of the M05 correlation functional
C>
C> @}
#endif
C>
C> \ingroup nwxc_priv
C> @{
C>
C> \brief The M05 correlation functional
C>
C> The M05 functional [1,2] is a meta-GGA of which this evaluates
C> the correlation component.
C>
C> Due to the form of the meta-GGAs we need to screen on the kinetic
C> energy density to ensure that LDA will be obtained when the kinetic
C> energy density goes to zero [3].
C>
C> ### References ###
C>
C> [1] Y Zhao, NE Schultz, DG Truhlar,
C>     "Exchange-correlation functional with broad accuracy for 
C>     metallic and nonmetallic compounds, kinetics, and 
C>     noncovalent interactions",
C>     J.Chem.Phys. <b>123</b>, 161103-161106 (2005), DOI:
C>     <a href="https://doi.org/10.1063/1.2126975">
C>     10.1063/1.2126975</a>.
C>
C> [2] Y Zhao, NE Schultz, DG Truhlar,
C>     "Design of density functionals by combining the method of
C>     constraint satisfaction parametrization for thermochemistry,
C>     thermochemical kinetics, and noncovalent interactions",
C>     J.Chem.Theory.Comput. <b>2</b>, 364-382 (2006), DOI:
C>     <a href="https://doi.org/10.1021/ct0502763">
C>     10.1021/ct0502763</a>.
C>
C> [3] J. Gr&auml;fenstein, D. Izotov, D. Cremer,
C>     "Avoiding singularity problems associated with meta-GGA exchange
C>     and correlation functionals containing the kinetic energy
C>     density", J. Chem. Phys. <b>127</b>, 214103 (2007), DOI:
C>     <a href="https://doi.org/10.1063/1.2800011">
C>     10.1063/1.2800011</a>.
C>
c    M05 and M05-2X correlation functional 
c           META GGA
C         utilizes ingredients:
c                              rho   -  density
c                              delrho - gradient of density
c                              tau (tauN)- K.S kinetic energy density
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
      Subroutine nwxc_c_m05_p(param, tol_rho, ipol, nq, wght, rho,
     &                        rgamma, tau, func)
#else
      Subroutine nwxc_c_m05(param, tol_rho, ipol, nq, wght, rho, rgamma,
     &                      tau, func)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      Subroutine nwxc_c_m05_d2(param, tol_rho, ipol, nq, wght, rho,
     &                         rgamma, tau, func)
#else
      Subroutine nwxc_c_m05_d3(param, tol_rho, ipol, nq, wght, rho,
     &                         rgamma, tau, func)
#endif
c
c$Id$
c
c
c     [a]       Zhao, Y.; Schultz, N. E.; Truhlar, D. G. J. Chem. Phys. 2005, 123, 161103;
c      Note that in this communication we interchanged cCab,i and cCss,i in Table 1. 
c     [b]       Zhao, Y.; Schultz, N. E.; Truhlar, D. G. J. Chem. Theory Comput. 2006, in press.    

#include "nwad.fh"
      implicit none
c
#include "intf_nwxc_c_lsda.fh"
#include "intf_nwxc_m05css.fh"
c
#include "nwxc_param.fh"
c
c     Input and other parameters
c
#if defined(NWAD_PRINT)
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      type(nwad_dble)::param(*)
      type(nwad_dble):: sopp0, sopp1,sopp2, sopp3, sopp4
#else
      double precision param(*)
      double precision  sopp0, sopp1,sopp2, sopp3, sopp4
#endif
#else
      double precision param(*) !< [Input] Parameters of functional
                                !< - param(1): \f$ c_{C\alpha\beta,0} \f$
                                !< - param(2): \f$ c_{C\alpha\beta,1} \f$
                                !< - param(3): \f$ c_{C\alpha\beta,2} \f$
                                !< - param(4): \f$ c_{C\alpha\beta,3} \f$
                                !< - param(5): \f$ c_{C\alpha\beta,4} \f$
                                !< - param(6): \f$ c_{C\sigma\sigma,0} \f$
                                !< - param(7): \f$ c_{C\sigma\sigma,1} \f$
                                !< - param(8): \f$ c_{C\sigma\sigma,2} \f$
                                !< - param(9): \f$ c_{C\sigma\sigma,3} \f$
                                !< - param(10): \f$ c_{C\sigma\sigma,4} \f$
      double precision  sopp0, sopp1,sopp2, sopp3, sopp4
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
c     Charge Density Gradient Norm
c
      type(nwad_dble)::rgamma(nq,*) !< [Input] The density gradient norm
c
c     Kinetic Energy Density
c
      type(nwad_dble)::tau(nq,*) !< [Input] The kinetic energy density
c
c     Functional values
c
      type(nwad_dble)::func(*) !< [Output] The functional value
c
c     Sampling Matrices for the XC Potential
c
c     double precision Amat(nq,*) !< [Output] Derivative wrt density
c     double precision Cmat(nq,*) !< [Output] Derivative wrt rgamma
c     double precision Mmat(nq,*) !< [Output] Derivative wrt tau
c
c     Threshold parameters
c
      double precision DTol,F1, F2, F3, F4,COpp 
      Data COpp/0.0031d0/,F1/1.0d0/,F2/2.0d0/,
     & F3/3.0d0/,F4/4.0d0/ 
c
      integer n

c    call to the m05css subroutine
      type(nwad_dble)::PA,GAA,TA,FA,EUA,ChiA
      double precision FPA,FGA,FTA,EUEGA,EUPA,ChiAP,ChiAG
      type(nwad_dble)::PB,GBB,TB,FB,EUB,ChiB
      double precision FPB,FGB,FTB,EUEGB,EUPB,ChiBP,ChiBG
c
      type(nwad_dble)::RS,Zeta,PotLC,P,U,W,EUEG
      double precision  sop
      double precision Pi, F6, F43, Pi34, F13, 
     &RSP,dZdA,dZdB,dLdS,dLdZ,d2LdSS,d2LdSZ,d2LdZZ
      double precision dUdChiA,dUdChiB,dUdPA,dUdPB,dUdGA,dUdGB,
     &dWdU,dWdPA,dWdPB, dWdGA,dWdGB,EUEGPA,EUEGPB

      
c
c     ======> BOTH SPIN-RESTRICETED AND UNRESTRICTED <======
c
      DTol=tol_rho
      sop=1.0d0
      sopp0= param(1)
      sopp1= param(2)
      sopp2= param(3)
      sopp3= param(4)
      sopp4= param(5)
c     if (ijzy.eq.1) then
C     Parameters for M05 Correlation
c        sopp0= 1.00000d0
c        sopp1= 3.78569d0
c        sopp2= -14.15261d0
c        sopp3= -7.46589d0
c        sopp4= 17.94491d0
c     elseif (ijzy.eq.2) then
C     Parameters for M05-2X Correlation
c        sopp0= 1.00000
c        sopp1= 1.09297
c        sopp2= -3.79171
c        sopp3= 2.82810
c        sopp4= -10.58909
c     endif

      Pi = F4*ATan(F1)
      F6=6.0d0
      F43 = F4 / F3
      Pi34 = F3 / (F4*Pi)
      F13 = F1 / F3

      do 20 n = 1, nq

       FA = 0.0d0
       EUA = 0.0d0
       ChiA = 0.0d0
       FB = 0.0d0
       EUB = 0.0d0
       ChiB = 0.0d0

       if (ipol.eq.1) then
         if (rho(n,R_T).le.DTol) goto 20
c
c    get the density, gradient, and tau for the alpha spin from the total 
c
         PA = rho(n,R_T)/F2
c        GAA = (    delrho(n,1,1)*delrho(n,1,1) +
c    &                 delrho(n,2,1)*delrho(n,2,1) +
c    &                 delrho(n,3,1)*delrho(n,3,1))/4
         GAA = rgamma(n,G_TT)/4.0d0
c  In the m05css subroutine, we use 2*TA as the tau, so we do not divide 
c  the tau by 2 here

         TA = tau(n,T_T)
         if (TA.gt.DTol) then
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)                  
#if defined(NWAD_PRINT)
           Call nwxc_m05css_p(param,DTol,PA,GAA,TA,FA,FPA,FGA,FTA,EUA,
     &                  ChiA,EUPA,ChiAP,ChiAG)
#else
           Call nwxc_m05css(param,DTol,PA,GAA,TA,FA,FPA,FGA,FTA,EUA,
     &                  ChiA,EUPA,ChiAP,ChiAG)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)                  
           Call nwxc_m05css_d2(param,DTol,PA,GAA,TA,FA,FPA,FGA,FTA,EUA,
     &                  ChiA,EUPA,ChiAP,ChiAG)
#else
           Call nwxc_m05css_d3(param,DTol,PA,GAA,TA,FA,FPA,FGA,FTA,EUA,
     &                  ChiA,EUPA,ChiAP,ChiAG)
#endif
         endif
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

c        Ec = Ec + 2.d0*FA*qwght(n)            !factor of 2 account for both spin
         func(n)        = func(n)+ FA*2d0*wght
c        Amat(n,D1_RA)  = Amat(n,D1_RA)+ FPA*wght
c        Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + FGA*wght
c        Mmat(n,D1_TA)  = Mmat(n,D1_TA) + FTA*wght
c         write (*,*) "PA,GAA,TA",PA,GAA,TA
c         write (*,*) "FPA,FGA,FTA",FPA,FGA,FTA 
c UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUnrestricted
      else  ! ipol=2
c
c        ======> SPIN-UNRESTRICTED <======
c
         PA  = 0.0d0
         PB  = 0.0d0
         GAA = 0.0d0
         GBB = 0.0d0
         TA  = 0.0d0
         TB  = 0.0d0
c
c       alpha
c

         if (rho(n,R_A).le.0.5d0*DTol) go to 25
         PA = rho(n,R_A)
c        GAA =   delrho(n,1,1)*delrho(n,1,1) +
c    &           delrho(n,2,1)*delrho(n,2,1) +
c    &          delrho(n,3,1)*delrho(n,3,1)
         GAA = rgamma(n,G_AA)
c
c  In the m05css subroutine, we use 2*TA as the tau 
c
         TA = 2.0d0*tau(n,T_A) 
         if (TA.gt.DTol) then
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)                  
#if defined(NWAD_PRINT)
           Call nwxc_m05css_p(param,DTol,PA,GAA,TA,FA,FPA,FGA,FTA,EUA,
     &                        ChiA,EUPA,ChiAP,ChiAG)
#else
           Call nwxc_m05css(param,DTol,PA,GAA,TA,FA,FPA,FGA,FTA,EUA,
     &                      ChiA,EUPA,ChiAP,ChiAG)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)                  
           Call nwxc_m05css_d2(param,DTol,PA,GAA,TA,FA,FPA,FGA,FTA,EUA,
     &                         ChiA,EUPA,ChiAP,ChiAG)
#else
           Call nwxc_m05css_d3(param,DTol,PA,GAA,TA,FA,FPA,FGA,FTA,EUA,
     &                         ChiA,EUPA,ChiAP,ChiAG)
#endif
         endif
c        Ec = Ec + FA*qwght(n)     
         func(n)=func(n)+ FA*wght
c        Amat(n,D1_RA)=Amat(n,D1_RA)+ FPA*wght
c        Cmat(n,D1_GAA)=  Cmat(n,D1_GAA) + FGA*wght
c        Mmat(n,D1_TA)=  Mmat(n,D1_TA) + FTA*wght
c
c  In the m05css subroutine, we use 2*TB as the tau, 
c
c
c       Beta 
c
 25      continue
         if (rho(n,R_B).le.0.5d0*DTol) go to 30
         PB = rho(n,R_B)
c        GBB =   delrho(n,1,2)*delrho(n,1,2) +
c    &           delrho(n,2,2)*delrho(n,2,2) +
c    &          delrho(n,3,2)*delrho(n,3,2)
         GBB = rgamma(n,G_BB)

         TB = 2.0d0*tau(n,T_B)
         if (TB.gt.DTol) then
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)                  
#if defined(NWAD_PRINT)
           Call nwxc_m05css_p(param,DTol,PB,GBB,TB,FB,FPB,FGB,FTB,EUB,
     &                        ChiB,EUPB,ChiBP,ChiBG)
#else
           Call nwxc_m05css(param,DTol,PB,GBB,TB,FB,FPB,FGB,FTB,EUB,
     &                      ChiB,EUPB,ChiBP,ChiBG)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)                  
           Call nwxc_m05css_d2(param,DTol,PB,GBB,TB,FB,FPB,FGB,FTB,EUB,
     &                         ChiB,EUPB,ChiBP,ChiBG)
#else
           Call nwxc_m05css_d3(param,DTol,PB,GBB,TB,FB,FPB,FGB,FTB,EUB,
     &                         ChiB,EUPB,ChiBP,ChiBG)
#endif
         endif
c        Ec = Ec + FB*qwght(n)          
         func(n)=func(n)+ FB*wght
c        Amat(n,2)= Amat(n,2)+ FPB*wght
c        Cmat(n,3)=  Cmat(n,3) + FGB*wght
c        Mmat(n,2)=  Mmat(n,2) + FTB*wght
      endif

 30   continue
      P = PA + PB
      
      If((PA.gt.0.5d0*DTol).or.(PB.gt.0.5d0*DTol)) then
          RS = (Pi34/P) ** F13
c         RSP = -RS/(F3*P)
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
          U = COpp*(ChiA+ChiB)/(F1 + COpp*(ChiA+ChiB))
          W = sopp0+U*(sopp1+U*(sopp2+U*(sopp3+U*sopp4)))
c         Ec = Ec + sop*EUEG*W*qwght(n)
          func(n)=func(n)+ sop*EUEG*W*wght
c         dUdChiA =COpp/(F1 + COpp*(ChiA+ChiB))**2
c         dUdChiB =COpp/(F1 + COpp*(ChiA+ChiB))**2
c         dUdPA= dUdChiA*ChiAP
c         dUdPB= dUdChiB*ChiBP
c         dUdGA= dUdChiA*ChiAG
c         dUdGB= dUdChiB*ChiBG
c         dWdU =sopp1+U*(F2*sopp2+U*(F3*sopp3+U*F4*sopp4))
c         dWdPA= dWdU*dUdPA
c         dWdPB= dWdU*dUdPB
c         dWdGA= dWdU*dUdGA
c         dWdGB= dWdU*dUdGB
c         EUEGPA = PotLC + P*dLdS*RSP + P*dLdZ*dZdA - EUPA
c         EUEGPB = PotLC + P*dLdS*RSP + P*dLdZ*dZdB - EUPB
c         if (ipol.eq.1) then 
c           Amat(n,D1_RA) = Amat(n,D1_RA)
c    &                    + sop*(EUEGPA*W + EUEG*dWdPA)*wght
c           Cmat(n,D1_GAA)=  Cmat(n,D1_GAA)
c    &                    + sop*(EUEG*dWdGA)*wght
c         else
c           Amat(n,D1_RA) = Amat(n,D1_RA)
c    &                    + sop*(EUEGPA*W + EUEG*dWdPA)*wght
c           Amat(n,D1_RB) = Amat(n,D1_RB)
c    &                    + sop*(EUEGPB*W + EUEG*dWdPB)*wght
c           Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + sop*EUEG*dWdGA*wght
c           Cmat(n,D1_GBB) = Cmat(n,D1_GBB) + sop*(EUEG*dWdGB)*wght
c         endif
      endIf
c      write (*,*) "PA, PB, GAA, GBB,ipol",PA, PB, GAA, GBB,ipol
c      write (*,*) "FA, FB,FGA, FGB",FA, FB,FGA, FGB
c      Stop
20    continue
      end


#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)                  
#if defined(NWAD_PRINT)
      Subroutine nwxc_m05css_p(param,DTol,PX,GX,TX,F,FP,FG,FT,EUEG,Chi,
     &                       EUEGP,ChiP,ChiG)
#else
      Subroutine nwxc_m05css(param,DTol,PX,GX,TX,F,FP,FG,FT,EUEG,Chi,
     &                       EUEGP,ChiP,ChiG)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)                  
      Subroutine nwxc_m05css_d2(param,DTol,PX,GX,TX,F,FP,FG,FT,EUEG,Chi,
     &                          EUEGP,ChiP,ChiG)
#else
      Subroutine nwxc_m05css_d3(param,DTol,PX,GX,TX,F,FP,FG,FT,EUEG,Chi,
     &                          EUEGP,ChiP,ChiG)
#endif
#include "nwad.fh"
      Implicit none
c
#include "intf_nwxc_c_lsda.fh"
C
C     Compute the same-spin part of the m05 correlation functional for one grid
C     point and one spin-case.
C
C
#ifdef NWAD_PRINT
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      type(nwad_dble)::param(*)
      type(nwad_dble)::sss0,sss1, sss2, sss3, sss4
#else
      double precision param(*)
      double precision sss0,sss1, sss2, sss3, sss4
#endif
#else
      double precision param(*) !< [Input] Parameters of functional
                                !< - param(1): \f$ c_{C\alpha\beta,0} \f$
                                !< - param(2): \f$ c_{C\alpha\beta,1} \f$
                                !< - param(3): \f$ c_{C\alpha\beta,2} \f$
                                !< - param(4): \f$ c_{C\alpha\beta,3} \f$
                                !< - param(5): \f$ c_{C\alpha\beta,4} \f$
                                !< - param(6): \f$ c_{C\sigma\sigma,0} \f$
                                !< - param(7): \f$ c_{C\sigma\sigma,1} \f$
                                !< - param(8): \f$ c_{C\sigma\sigma,2} \f$
                                !< - param(9): \f$ c_{C\sigma\sigma,3} \f$
                                !< - param(10): \f$ c_{C\sigma\sigma,4} \f$
      double precision sss0,sss1, sss2, sss3, sss4
#endif
      integer ijzy
      type(nwad_dble)::PX, GX, TX, F, RS, PotLC
      type(nwad_dble)::EUEG,D,E,U,W,Chi,Fscc,Zeta
      double precision FP, FG, FT, DTol
      double precision EUEGP, ChiP, ChiG
      double precision Zero, Pt25, F1, F2, F3, F4, F5, F6, F8, F11
      double precision ss, Css
      double precision Pi, Pi34, F13, F23, F43, F53, F83, F113
      double precision FDUEG, RSP, dFsccP, dFsccG
      double precision dFsccT, dUdChi, dWdU, dWdP, dWdG 
      double precision d2LdSS,d2LdSZ,d2LdZZ,dLdS,dLdZ
     


      Data Zero/0.0d0/, Pt25/0.25d0/, F1/1.0d0/, F2/2.0d0/, F3/3.0d0/,
     $  F4/4.0d0/, F5/5.0d0/, F6/6.0d0/, F8/8.0d0/, F11/11.0d0/,
     $  Css/0.06d0/
C
c      DTol=1.0D-7
c      write (*,*) DTol
      ss=1.0
      sss0=  param(6)
      sss1=  param(7)
      sss2=  param(8)
      sss3=  param(9)
      sss4=  param(10)
c     if (ijzy.eq.1) then
C     Parameters for M05 Correlation
c        sss0=  1.00000d0
c        sss1=  3.77344d0
c        sss2=  -26.04463d0
c        sss3=  30.69913d0
c        sss4=  -9.22695d0
c     elseif (ijzy.eq.2) then
C     Parameters for M05-2X Correlation
c        sss0=  1.00000d0
c        sss1=  -3.05430d0
c        sss2=  7.61854d0
c        sss3=  1.47665d0
c        sss4=  -11.92365d0
c     endif
      
      If ((PX.le.DTol).or.(TX.le.dtol)) then
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
C        DUEG = FDUEG*PX**F53
        Chi = GX/PX**F83
        U = Css*Chi/(F1 + Css*Chi)
        W = sss0+U*(sss1+U*(sss2+U*(sss3+U*sss4)))
        Fscc=D/TX
        E = Fscc*W*EUEG
        F = E*ss
c       RSP = -RS/(F3*Px)
c       ChiG = F1/PX**F83
c       ChiP = -F83*Chi/PX
c       dFsccP=Pt25*GX/(TX*PX**2)
c       dFsccG=-Pt25/(TX*PX)
c       dFsccT=Pt25*GX/(PX*TX**2)
c       dUdChi=Css/((F1+Css*Chi)**2)
c       dWdU=sss1+U*(F2*sss2+U*(F3*sss3+U*F4*sss4))
c       dWdP=dWdU*dUdChi*ChiP
c       dWdG=dWdU*dUdChi*ChiG 
c       EUEGP = PotLC + PX*dLdS*RSP
c       FP = ss*(dFsccP*W*EUEG 
c    $                 + Fscc*dWdP*EUEG
c    $                 + Fscc*W*EUEGP)
c       FG = ss*(dFsccG*W*EUEG
c    $                 + Fscc*dWdG*EUEG)

c       FT = ss*(dFsccT*W*EUEG)
       Endif

       Return
       End


c------------------------------------------------------------------------
C> \brief Calculate the dlDF correlation functional
C>
C> Calculate the dlDF correlation functional [1].
C>
C> ### References ###
C>
C> [1] K Pernal, R Podeszwa, K Patkowski, K Szalewicz,
C> "Dispersionless density functional theory",
C> Phys.Rev.Lett. <b>103</b>, 263201-263204 (2009), DOI:
C> <a href="https://doi.org/10.1103/PhysRevLett.103.263201">
C> 10.1103/PhysRevLett.103.263201</a>.
C>
c     dlDF correlation functional 
c           META GGA
C         utilizes ingredients:
c                              rho   -  density
c                              delrho - gradient of density
c                              tau (tauN)- K.S kinetic energy density
c     References: 
c     [a]	Pernal,Podeszwa,Patkowski,Szalewicz, PRL 103 263201 (2009)
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
      Subroutine nwxc_c_dldf_p(tol_rho, ipol, nq, wght, rho, rgamma,
     &                         tau, func)
#else
      Subroutine nwxc_c_dldf(tol_rho, ipol, nq, wght, rho, rgamma, tau,
     &                       func)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      Subroutine nwxc_c_dldf_d2(tol_rho, ipol, nq, wght, rho, rgamma,
     &                          tau, func)
#else
      Subroutine nwxc_c_dldf_d3(tol_rho, ipol, nq, wght, rho, rgamma,
     &                          tau, func)
#endif

c
#include "nwad.fh"
      implicit none
c
#include "intf_nwxc_c_lsda.fh"
#include "intf_nwxc_dldfcss.fh"
c
#include "nwxc_param.fh"
c
c     Input and other parameters
c
      double precision tol_rho !< [Input] The lower limit on the density
      integer nq               !< [Input] The number of points
      integer ipol             !< [Input] The number of spin channels
      double precision wght    !< [Input] The weight of the functional
c
c     Charge Density 
c
      type(nwad_dble)::rho(nq,*) !< [Input] The density
c
c     Charge Density Gradient Norm
c
      type(nwad_dble)::rgamma(nq,*) !< [Input] The density gradient norm
c
c     Kinetic Energy Density
c
      type(nwad_dble)::tau(nq,*) !< [Input] The kinetic energy density
c
c     Functional values
c
      type(nwad_dble)::func(*) !< [Output] The functional value
c
c     Sampling Matrices for the XC Potential
c
c     double precision Amat(nq,*) !< [Output] Derivative wrt density
c     double precision Cmat(nq,*) !< [Output] Derivative wrt rgamma
c     double precision Mmat(nq,*) !< [Output] Derivative wrt tau
c
c     Threshold parameters
c
      double precision DTol,F1, F2, F3, F4,COpp 
      Data COpp/0.0031d0/,F1/1.0d0/,F2/2.0d0/,
     & F3/3.0d0/,F4/4.0d0/ 

      integer n

c    call to the dldfcss subroutine
      type(nwad_dble)::PA,GAA,TA,FA,EUA,ChiA
      double precision FPA,FGA,FTA,EUEGA,EUPA,ChiAP,ChiAG
      type(nwad_dble)::PB,GBB,TB,FB,EUB,ChiB
      double precision FPB,FGB,FTB,EUEGB,EUPB,ChiBP,ChiBG
c
      type(nwad_dble)::P,RS,Zeta,EUEG,U,W,PotLC
      double precision  sop, sopp0, sopp1,sopp2, sopp3, sopp4
      double precision Pi, F6, F43, Pi34, F13, 
     &RSP,dZdA,dZdB,dLdS,dLdZ,d2LdSS,d2LdSZ,d2LdZZ
      double precision dUdChiA,dUdChiB,dUdPA,dUdPB,dUdGA,dUdGB,
     &dWdU,dWdPA,dWdPB, dWdGA,dWdGB,EUEGPA,EUEGPB

      
c
c     ======> BOTH SPIN-RESTRICETED AND UNRESTRICTED <======
c
      DTol=tol_rho
      sop=1.0d0
C
C     Parameters for dlDF Correlation
         sopp0= 1.0d0
         sopp1= 5.9515308d0
         sopp2= -11.1602877d0
         sopp3= 0.0d0
         sopp4= 0.0d0


      Pi = F4*ATan(F1)
      F6=6.0d0
      F43 = F4 / F3
      Pi34 = F3 / (F4*Pi)
      F13 = F1 / F3

      do 20 n = 1, nq

       FA = 0.0d0
       EUA = 0.0d0
       ChiA = 0.0d0
       FB = 0.0d0
       EUB = 0.0d0
       ChiB = 0.0d0

       if (ipol.eq.1) then
         if (rho(n,R_T).lt.DTol) goto 20
c
c    get the density, gradient, and tau for the alpha spin from the total 
c
         PA = rho(n,R_T)/F2
c        GAA = (    delrho(n,1,1)*delrho(n,1,1) +
c    &                 delrho(n,2,1)*delrho(n,2,1) +
c    &                 delrho(n,3,1)*delrho(n,3,1))/4
         GAA = rgamma(n,G_TT)/4.0d0
c  In the dldfcss subroutine, we use 2*TA as the tau, so we do not divide 
c  the tau by 2 here

         TA = tau(n,T_T)
         if (TA.ge.DTol) then
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
           Call nwxc_dldfcss_p(DTol,PA,GAA,TA,FA,FPA,FGA,FTA,EUA,
     &                         ChiA,EUPA,ChiAP,ChiAG)
#else
           Call nwxc_dldfcss(DTol,PA,GAA,TA,FA,FPA,FGA,FTA,EUA,
     &                       ChiA,EUPA,ChiAP,ChiAG)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
           Call nwxc_dldfcss_d2(DTol,PA,GAA,TA,FA,FPA,FGA,FTA,EUA,
     &                          ChiA,EUPA,ChiAP,ChiAG)
#else
           Call nwxc_dldfcss_d3(DTol,PA,GAA,TA,FA,FPA,FGA,FTA,EUA,
     &                          ChiA,EUPA,ChiAP,ChiAG)
#endif
         endif
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

c        Ec = Ec + 2.d0*FA*qwght(n)            !factor of 2 account for both spin
         func(n)=func(n)+ FA*2d0*wght
c        Amat(n,D1_RA)  = Amat(n,D1_RA)+ FPA*wght
c        Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + FGA*wght
c        Mmat(n,D1_TA)  = Mmat(n,D1_TA) + FTA*wght
c         write (*,*) "PA,GAA,TA",PA,GAA,TA
c         write (*,*) "FPA,FGA,FTA",FPA,FGA,FTA 
c UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUnrestricted
      else  ! ipol=2
c
c        ======> SPIN-UNRESTRICTED <======
c
         PA  = 0.0d0
         PB  = 0.0d0
         GAA = 0.0d0
         GBB = 0.0d0
         TA  = 0.0d0
         TB  = 0.0d0
c
c       alpha
c

         if (rho(n,R_A).le.0.5d0*DTol) go to 25
         PA = rho(n,R_A)
c        GAA =   delrho(n,1,1)*delrho(n,1,1) +
c    &           delrho(n,2,1)*delrho(n,2,1) +
c    &          delrho(n,3,1)*delrho(n,3,1)
         GAA = rgamma(n,G_AA)
c
c  In the dldfcss subroutine, we use 2*TA as the tau 
c
         TA = 2.0d0*tau(n,T_A)
         if (TA.gt.DTol) then
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
           Call nwxc_dldfcss_p(DTol,PA,GAA,TA,FA,FPA,FGA,FTA,EUA,
     &                         ChiA,EUPA,ChiAP,ChiAG)
#else
           Call nwxc_dldfcss(DTol,PA,GAA,TA,FA,FPA,FGA,FTA,EUA,
     &                       ChiA,EUPA,ChiAP,ChiAG)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
           Call nwxc_dldfcss_d2(DTol,PA,GAA,TA,FA,FPA,FGA,FTA,EUA,
     &                          ChiA,EUPA,ChiAP,ChiAG)
#else
           Call nwxc_dldfcss_d3(DTol,PA,GAA,TA,FA,FPA,FGA,FTA,EUA,
     &                          ChiA,EUPA,ChiAP,ChiAG)
#endif
         endif
c        Ec = Ec + FA*qwght(n)     
         func(n)=func(n)+ FA*wght
c        Amat(n,D1_RA)=Amat(n,D1_RA)+ FPA*wght
c        Cmat(n,D1_GAA)=  Cmat(n,D1_GAA) + FGA*wght
c        Mmat(n,D1_TA)=  Mmat(n,D1_TA) + FTA*wght
c
c  In the dldfcss subroutine, we use 2*TB as the tau, 
c
c
c       Beta 
c
 25      continue
         if (rho(n,R_B).le.0.5d0*DTol) go to 30
         PB = rho(n,R_B)
c        GBB =   delrho(n,1,2)*delrho(n,1,2) +
c    &           delrho(n,2,2)*delrho(n,2,2) +
c    &          delrho(n,3,2)*delrho(n,3,2)
         GBB = rgamma(n,G_BB)

         TB = 2.0d0*tau(n,T_B) 
         if (TB.gt.DTol) then
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
           Call nwxc_dldfcss_p(DTol,PB,GBB,TB,FB,FPB,FGB,FTB,EUB,
     &                         ChiB,EUPB,ChiBP,ChiBG)
#else
           Call nwxc_dldfcss(DTol,PB,GBB,TB,FB,FPB,FGB,FTB,EUB,
     &                       ChiB,EUPB,ChiBP,ChiBG)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
           Call nwxc_dldfcss_d2(DTol,PB,GBB,TB,FB,FPB,FGB,FTB,EUB,
     &                          ChiB,EUPB,ChiBP,ChiBG)
#else
           Call nwxc_dldfcss_d3(DTol,PB,GBB,TB,FB,FPB,FGB,FTB,EUB,
     &                          ChiB,EUPB,ChiBP,ChiBG)
#endif
         endif
c        Ec = Ec + FB*qwght(n)          
         func(n)=func(n)+ FB*wght
c        Amat(n,D1_RB)= Amat(n,D1_RB)+ FPB*wght
c        Cmat(n,D1_GBB)=  Cmat(n,D1_GBB) + FGB*wght
c        Mmat(n,D1_TB)=  Mmat(n,D1_TB) + FTB*wght
      endif

 30   continue
      P = PA + PB
      
      If((PA.gt.0.5d0*DTol).or.(PB.gt.0.5d0*DTol)) then
          RS = (Pi34/P) ** F13
c         RSP = -RS/(F3*P)
          Zeta = (PA-PB)/P
c         dZdA = (F1-Zeta)/P
c         dZdB = (-F1-Zeta)/P
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
          Call nwxc_c_lsda_p(dtol,
     &         RS,Zeta,PotLC,dLdS,dLdZ,d2LdSS,d2LdSZ,d2LdZZ)
#else
          Call nwxc_c_lsda(dtol,
     &         RS,Zeta,PotLC,dLdS,dLdZ,d2LdSS,d2LdSZ,d2LdZZ)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
          Call nwxc_c_lsda_d2(dtol,
     &         RS,Zeta,PotLC,dLdS,dLdZ,d2LdSS,d2LdSZ,d2LdZZ)
#else
          Call nwxc_c_lsda_d3(dtol,
     &         RS,Zeta,PotLC,dLdS,dLdZ,d2LdSS,d2LdSZ,d2LdZZ)
#endif
          EUEG = P*PotLC - EUA - EUB
          U = COpp*(ChiA+ChiB)/(F1 + COpp*(ChiA+ChiB))
          W = sopp0+U*(sopp1+U*(sopp2+U*(sopp3+U*sopp4)))
c         Ec = Ec + sop*EUEG*W*qwght(n)
          func(n)=func(n)+ sop*EUEG*W*wght
c         dUdChiA =COpp/(F1 + COpp*(ChiA+ChiB))**2
c         dUdChiB =COpp/(F1 + COpp*(ChiA+ChiB))**2
c         dUdPA= dUdChiA*ChiAP
c         dUdPB= dUdChiB*ChiBP
c         dUdGA= dUdChiA*ChiAG
c         dUdGB= dUdChiB*ChiBG
c         dWdU =sopp1+U*(F2*sopp2+U*(F3*sopp3+U*F4*sopp4))
c         dWdPA= dWdU*dUdPA
c         dWdPB= dWdU*dUdPB
c         dWdGA= dWdU*dUdGA
c         dWdGB= dWdU*dUdGB
c         EUEGPA = PotLC + P*dLdS*RSP + P*dLdZ*dZdA - EUPA
c         EUEGPB = PotLC + P*dLdS*RSP + P*dLdZ*dZdB - EUPB
c         if (ipol.eq.1) then 
c           Amat(n,D1_RA) = Amat(n,D1_RA)
c    &                    + sop*(EUEGPA*W + EUEG*dWdPA)*wght
c           Cmat(n,D1_GAA)=  Cmat(n,D1_GAA)
c    &                    + sop*(EUEG*dWdGA)*wght
c         else
c           Amat(n,D1_RA) = Amat(n,D1_RA)
c    &                    + sop*(EUEGPA*W + EUEG*dWdPA)*wght
c           Amat(n,D1_RB) = Amat(n,D1_RB)
c    &                    + sop*(EUEGPB*W + EUEG*dWdPB)*wght
c           Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + sop*EUEG*dWdGA*wght
c           Cmat(n,D1_GBB) = Cmat(n,D1_GBB) + sop*(EUEG*dWdGB)*wght
c         endif
      endIf
c      write (*,*) "PA, PB, GAA, GBB,ipol",PA, PB, GAA, GBB,ipol
c      write (*,*) "FA, FB,FGA, FGB",FA, FB,FGA, FGB
c      Stop
20    continue
      end


#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
      Subroutine nwxc_dldfcss_p(DTol,PX,GX,TX,F,FP,FG,FT,EUEG,Chi,EUEGP,
     &                   ChiP,ChiG)
#else
      Subroutine nwxc_dldfcss(DTol,PX,GX,TX,F,FP,FG,FT,EUEG,Chi,EUEGP,
     &                   ChiP,ChiG)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      Subroutine nwxc_dldfcss_d2(DTol,PX,GX,TX,F,FP,FG,FT,EUEG,Chi,
     &                   EUEGP, ChiP,ChiG)
#else
      Subroutine nwxc_dldfcss_d3(DTol,PX,GX,TX,F,FP,FG,FT,EUEG,Chi,
     &                   EUEGP, ChiP,ChiG)
#endif
#include "nwad.fh"
      Implicit none
c
#include "intf_nwxc_c_lsda.fh"
c
C
C     Compute the same-spin part of the dldf correlation functional for one grid
C     point and one spin-case.
C
C
      type(nwad_dble)::PX, GX, TX, F, RS, D, E, W, U, Fscc
      type(nwad_dble)::EUEG, Chi, Zeta
      type(nwad_dble)::PotLC
      double precision FP, FG, FT, DTol
      double precision EUEGP, ChiP, ChiG
      double precision Zero, Pt25, F1, F2, F3, F4, F5, F6, F8, F11
      double precision ss, sss0,sss1, sss2, sss3, sss4, Css
      double precision Pi, Pi34, F13, F23, F43, F53, F83, F113
      double precision FDUEG, RSP, dFsccP, dFsccG
      double precision dFsccT, dUdChi, dWdU, dWdP, dWdG 
      double precision d2LdSS,d2LdSZ,d2LdZZ,dLdS,dLdZ
     


      Data Zero/0.0d0/, Pt25/0.25d0/, F1/1.0d0/, F2/2.0d0/, F3/3.0d0/,
     $  F4/4.0d0/, F5/5.0d0/, F6/6.0d0/, F8/8.0d0/, F11/11.0d0/,
     $  Css/0.06d0/
C
c      DTol=1.0D-7
c      write (*,*) DTol
      ss=1.0
C
C     Parameters for dlDF Correlation
         sss0=  1.0d0
         sss1=  -2.5960897d0
         sss2=  2.2233793d0
         sss3=  0.0d0
         sss4=  0.0d0

      
      If ((PX.le.DTol).or.(TX.le.dtol)) then
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
C        DUEG = FDUEG*PX**F53
        Chi = GX/PX**F83
        U = Css*Chi/(F1 + Css*Chi)
        W = sss0+U*(sss1+U*(sss2+U*(sss3+U*sss4)))
        Fscc=D/TX
        E = Fscc*W*EUEG
        F = E*ss
c       RSP = -RS/(F3*Px)
c       ChiG = F1/PX**F83
c       ChiP = -F83*Chi/PX
c       dFsccP=Pt25*GX/(TX*PX**2)
c       dFsccG=-Pt25/(TX*PX)
c       dFsccT=Pt25*GX/(PX*TX**2)
c       dUdChi=Css/((F1+Css*Chi)**2)
c       dWdU=sss1+U*(F2*sss2+U*(F3*sss3+U*F4*sss4))
c       dWdP=dWdU*dUdChi*ChiP
c       dWdG=dWdU*dUdChi*ChiG 
c       EUEGP = PotLC + PX*dLdS*RSP
c       FP = ss*(dFsccP*W*EUEG 
c    $                 + Fscc*dWdP*EUEG
c    $                 + Fscc*W*EUEGP)
c       FG = ss*(dFsccG*W*EUEG
c    $                 + Fscc*dWdG*EUEG)

c       FT = ss*(dFsccT*W*EUEG)
       Endif

       Return
       End


#ifndef NWAD_PRINT
#define NWAD_PRINT
c
c     Compile source again for the 2nd derivative case
c
#include "nwxc_c_m05.F"
#endif
#ifndef SECOND_DERIV
#define SECOND_DERIV
c
c     Compile source again for the 2nd derivative case
c
#include "nwxc_c_m05.F"
#endif
#ifndef THIRD_DERIV
#define THIRD_DERIV
c
c     Compile source again for the 3rd derivative case
c
#include "nwxc_c_m05.F"
#endif
#undef NWAD_PRINT
C> @}
