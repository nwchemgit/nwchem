C> \ingroup nwpwxc
C> @{
C>
C> \file nwpwxc_c_m05.F
C> Implementation of the M05 correlation functional
C>
C> @}
C>
C> \ingroup nwpwxc_priv
C> @{
C>
C> \brief The M05 correlation functional
C>
C> The M05 functional [1,2] is a meta-GGA of which this evaluates
C> the correlation component.
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
c    M05 and M05-2X correlation functional 
c           META GGA
C         utilizes ingredients:
c                              rho   -  density
c                              delrho - gradient of density
c                              tau (tauN)- K.S kinetic energy density

      Subroutine nwpwxc_c_m05(param,tol_rho,ipol,nq,wght,rho,rgamma,
     &                      tau, func, Amat, Cmat, Mmat)
c
c$Id$
c
c
c     [a]       Zhao, Y.; Schultz, N. E.; Truhlar, D. G. J. Chem. Phys. 2005, 123, 161103;
c      Note that in this communication we interchanged cCab,i and cCss,i in Table 1. 
c     [b]       Zhao, Y.; Schultz, N. E.; Truhlar, D. G. J. Chem. Theory Comput. 2006, in press.    

      implicit none
#include "nwpwxc_param.fh"
c
c     Input and other parameters
c
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
      double precision tol_rho !< [Input] The lower limit on the density
      integer nq               !< [Input] The number of points
      integer ipol             !< [Input] The number of spin channels
      double precision wght    !< [Input] The weight of the functional
c
c     Charge Density 
c
      double precision rho(nq,*) !< [Input] The density
c
c     Charge Density Gradient Norm
c
      double precision rgamma(nq,*) !< [Input] The density gradient norm
c
c     Kinetic Energy Density
c
      double precision tau(nq,*) !< [Input] The kinetic energy density
c
c     Functional values
c
      double precision func(*) !< [Output] The functional value
c
c     Sampling Matrices for the XC Potential
c
      double precision Amat(nq,*) !< [Output] Derivative wrt density
      double precision Cmat(nq,*) !< [Output] Derivative wrt rgamma
      double precision Mmat(nq,*) !< [Output] Derivative wrt tau
c
c     Threshold parameters
c
      double precision DTol,F1, F2, F3, F4,COpp 
      Data COpp/0.0031d0/,F1/1.0d0/,F2/2.0d0/,
     & F3/3.0d0/,F4/4.0d0/ 
c
      integer n

c    call to the m05css subroutine
      double precision PA,GAA,TA,FA,FPA,FGA,FTA,EUA,EUEGA,ChiA,EUPA
     &,ChiAP,ChiAG
      double precision PB,GBB,TB,FB,FPB,FGB,FTB,EUB,EUEGB,ChiB,EUPB
     &,ChiBP,ChiBG
c
      double precision  sop, sopp0, sopp1,sopp2, sopp3, sopp4
      double precision Pi, F6, F43, Pi34, F13, 
     &RS,RSP,Zeta,dZdA,dZdB,PotLC,dLdS,dLdZ,d2LdSS,d2LdSZ,d2LdZZ
      double precision P, EUEG, U, W
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
       if (rho(n,R_T).lt.DTol) goto 20
       if (ipol.eq.1) then
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
                  
         Call nwpwxc_m05css(param,DTol,PA,GAA,TA,FA,FPA,FGA,FTA,EUA,
     &                ChiA,EUPA,ChiAP,ChiAG)
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
         Amat(n,D1_RA)  = Amat(n,D1_RA)+ FPA*wght
         Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + FGA*wght
         Mmat(n,D1_TA)  = Mmat(n,D1_TA) + FTA*wght
c         write (*,*) "PA,GAA,TA",PA,GAA,TA
c         write (*,*) "FPA,FGA,FTA",FPA,FGA,FTA 
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
c  In the m05css subroutine, we use 2*TA as the tau 
c
         TA = 2*tau(n,T_A) 

         Call nwpwxc_m05css(param,DTol,PA,GAA,TA,FA,FPA,FGA,FTA,EUA,
     &                ChiA,EUPA,ChiAP,ChiAG)
c        Ec = Ec + FA*qwght(n)     
         func(n)=func(n)+ FA*wght
         Amat(n,D1_RA)=Amat(n,D1_RA)+ FPA*wght
         Cmat(n,D1_GAA)=  Cmat(n,D1_GAA) + FGA*wght
         Mmat(n,D1_TA)=  Mmat(n,D1_TA) + FTA*wght
c
c  In the m05css subroutine, we use 2*TB as the tau, 
c
c
c       Beta 
c
 25      continue
         PB = rho(n,R_B)
         if (PB.le.0.5d0*DTol) go to 30
c        GBB =   delrho(n,1,2)*delrho(n,1,2) +
c    &           delrho(n,2,2)*delrho(n,2,2) +
c    &          delrho(n,3,2)*delrho(n,3,2)
         GBB = rgamma(n,G_BB)

         TB = 2*tau(n,T_B)

         Call nwpwxc_m05css(param,DTol,PB,GBB,TB,FB,FPB,FGB,FTB,EUB,
     &                ChiB,EUPB,ChiBP,ChiBG)
c        Ec = Ec + FB*qwght(n)          
         func(n)=func(n)+ FB*wght
         Amat(n,2)= Amat(n,2)+ FPB*wght
         Cmat(n,3)=  Cmat(n,3) + FGB*wght
         Mmat(n,2)=  Mmat(n,2) + FTB*wght
      endif

 30   continue
      P = PA + PB
      
      If((PA.gt.0.5d0*DTol).and.(PB.gt.0.5d0*DTol)) then
          RS = (Pi34/P) ** F13
          RSP = -RS/(F3*P)
          Zeta = (PA-PB)/P
          dZdA = (F1-Zeta)/P
          dZdB = (-F1-Zeta)/P
          Call nwpwxc_c_lsda(dtol,
     D         RS,Zeta,PotLC,dLdS,dLdZ,d2LdSS,d2LdSZ,
     $         d2LdZZ)
          EUEG = P*PotLC - EUA - EUB
          U = COpp*(ChiA+ChiB)/(F1 + COpp*(ChiA+ChiB))
          W = sopp0+U*(sopp1+U*(sopp2+U*(sopp3+U*sopp4)))
c         Ec = Ec + sop*EUEG*W*qwght(n)
          func(n)=func(n)+ sop*EUEG*W*wght
          dUdChiA =COpp/(F1 + COpp*(ChiA+ChiB))**2
          dUdChiB =COpp/(F1 + COpp*(ChiA+ChiB))**2
          dUdPA= dUdChiA*ChiAP
          dUdPB= dUdChiB*ChiBP
          dUdGA= dUdChiA*ChiAG
          dUdGB= dUdChiB*ChiBG
          dWdU =sopp1+U*(F2*sopp2+U*(F3*sopp3+U*F4*sopp4))
          dWdPA= dWdU*dUdPA
          dWdPB= dWdU*dUdPB
          dWdGA= dWdU*dUdGA
          dWdGB= dWdU*dUdGB
          EUEGPA = PotLC + P*dLdS*RSP + P*dLdZ*dZdA - EUPA
          EUEGPB = PotLC + P*dLdS*RSP + P*dLdZ*dZdB - EUPB
          if (ipol.eq.1) then 
            Amat(n,D1_RA) = Amat(n,D1_RA)
     &                    + sop*(EUEGPA*W + EUEG*dWdPA)*wght
            Cmat(n,D1_GAA)=  Cmat(n,D1_GAA)
     &                    + sop*(EUEG*dWdGA)*wght
          else
            Amat(n,D1_RA) = Amat(n,D1_RA)
     &                    + sop*(EUEGPA*W + EUEG*dWdPA)*wght
            Amat(n,D1_RB) = Amat(n,D1_RB)
     &                    + sop*(EUEGPB*W + EUEG*dWdPB)*wght
            Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + sop*EUEG*dWdGA*wght
            Cmat(n,D1_GBB) = Cmat(n,D1_GBB) + sop*(EUEG*dWdGB)*wght
          endif
      endIf
c      write (*,*) "PA, PB, GAA, GBB,ipol",PA, PB, GAA, GBB,ipol
c      write (*,*) "FA, FB,FGA, FGB",FA, FB,FGA, FGB
c      Stop
20    continue
      end

      Subroutine nwpwxc_c_m05_d2()
      implicit none
      call errquit(' cm05: d2 not coded ',0,0)
      return
      end




      Subroutine nwpwxc_m05css(param,DTol,PX,GX,TX,F,FP,FG,FT,EUEG,Chi,
     &                       EUEGP,ChiP,ChiG)
      Implicit none
C
C     Compute the same-spin part of the m05 correlation functional for one grid
C     point and one spin-case.
C
C
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
      integer ijzy
      double precision PX, GX, TX, F, FP, FG, FT, DTol
      double precision EUEG, Chi, EUEGP, ChiP, ChiG
      double precision Zero, Pt25, F1, F2, F3, F4, F5, F6, F8, F11
      double precision ss, sss0,sss1, sss2, sss3, sss4, Css
      double precision Pi, Pi34, F13, F23, F43, F53, F83, F113
      double precision RS, FDUEG, D, Fscc, RSP, dFsccP, dFsccG
      double precision E, W, U, dFsccT, dUdChi, dWdU, dWdP, dWdG 
      double precision d2LdSS,d2LdSZ,d2LdZZ,PotLC,dLdS,dLdZ
     


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
        Call nwpwxc_c_lsda(dtol,
     D       RS,F1,PotLC,dLdS,dLdZ,d2LdSS,d2LdSZ,d2LdZZ)
        EUEG = PX*PotLC
        D = TX - Pt25*GX/PX
C        DUEG = FDUEG*PX**F53
        Chi = GX/PX**F83
        U = Css*Chi/(F1 + Css*Chi)
        W = sss0+U*(sss1+U*(sss2+U*(sss3+U*sss4)))
        Fscc=D/TX
        E = Fscc*W*EUEG
        F = E*ss
        RSP = -RS/(F3*Px)
        ChiG = F1/PX**F83
        ChiP = -F83*Chi/PX
        dFsccP=Pt25*GX/(TX*PX**2)
        dFsccG=-Pt25/(TX*PX)
        dFsccT=Pt25*GX/(PX*TX**2)
        dUdChi=Css/((F1+Css*Chi)**2)
        dWdU=sss1+U*(F2*sss2+U*(F3*sss3+U*F4*sss4))
        dWdP=dWdU*dUdChi*ChiP
        dWdG=dWdU*dUdChi*ChiG 
        EUEGP = PotLC + PX*dLdS*RSP
        FP = ss*(dFsccP*W*EUEG 
     $                 + Fscc*dWdP*EUEG
     $                 + Fscc*W*EUEGP)
        FG = ss*(dFsccG*W*EUEG
     $                 + Fscc*dWdG*EUEG)

        FT = ss*(dFsccT*W*EUEG)
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

      Subroutine nwpwxc_c_dldf(tol_rho,ipol,nq,wght,rho,rgamma,tau,
     &                       func, Amat, Cmat, Mmat)

c
      implicit none
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
      double precision rho(nq,*) !< [Input] The density
c
c     Charge Density Gradient Norm
c
      double precision rgamma(nq,*) !< [Input] The density gradient norm
c
c     Kinetic Energy Density
c
      double precision tau(nq,*) !< [Input] The kinetic energy density
c
c     Functional values
c
      double precision func(*) !< [Output] The functional value
c
c     Sampling Matrices for the XC Potential
c
      double precision Amat(nq,*) !< [Output] Derivative wrt density
      double precision Cmat(nq,*) !< [Output] Derivative wrt rgamma
      double precision Mmat(nq,*) !< [Output] Derivative wrt tau
c
c     Threshold parameters
c
      double precision DTol,F1, F2, F3, F4,COpp 
      Data COpp/0.0031d0/,F1/1.0d0/,F2/2.0d0/,
     & F3/3.0d0/,F4/4.0d0/ 

      integer n

c    call to the dldfcss subroutine
      double precision PA,GAA,TA,FA,FPA,FGA,FTA,EUA,EUEGA,ChiA,EUPA
     &,ChiAP,ChiAG
      double precision PB,GBB,TB,FB,FPB,FGB,FTB,EUB,EUEGB,ChiB,EUPB
     &,ChiBP,ChiBG
c
      double precision  sop, sopp0, sopp1,sopp2, sopp3, sopp4
      double precision Pi, F6, F43, Pi34, F13, 
     &RS,RSP,Zeta,dZdA,dZdB,PotLC,dLdS,dLdZ,d2LdSS,d2LdSZ,d2LdZZ
      double precision P, EUEG, U, W
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
       if (rho(n,R_T).lt.DTol) goto 20
       if (ipol.eq.1) then
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
                  
         Call nwpwxc_dldfcss(DTol,PA,GAA,TA,FA,FPA,FGA,FTA,EUA,
     &                ChiA,EUPA,ChiAP,ChiAG)
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
         Amat(n,D1_RA)  = Amat(n,D1_RA)+ FPA*wght
         Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + FGA*wght
         Mmat(n,D1_TA)  = Mmat(n,D1_TA) + FTA*wght
c         write (*,*) "PA,GAA,TA",PA,GAA,TA
c         write (*,*) "FPA,FGA,FTA",FPA,FGA,FTA 
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
c  In the dldfcss subroutine, we use 2*TA as the tau 
c
         TA = 2*tau(n,T_A)

         Call nwpwxc_dldfcss(DTol,PA,GAA,TA,FA,FPA,FGA,FTA,EUA,
     &                ChiA,EUPA,ChiAP,ChiAG)
c        Ec = Ec + FA*qwght(n)     
         func(n)=func(n)+ FA*wght
         Amat(n,D1_RA)=Amat(n,D1_RA)+ FPA*wght
         Cmat(n,D1_GAA)=  Cmat(n,D1_GAA) + FGA*wght
         Mmat(n,D1_TA)=  Mmat(n,D1_TA) + FTA*wght
c
c  In the dldfcss subroutine, we use 2*TB as the tau, 
c
c
c       Beta 
c
 25      continue
         PB = rho(n,R_B)
         if (PB.le.0.5d0*DTol) go to 30
c        GBB =   delrho(n,1,2)*delrho(n,1,2) +
c    &           delrho(n,2,2)*delrho(n,2,2) +
c    &          delrho(n,3,2)*delrho(n,3,2)
         GBB = rgamma(n,G_BB)

         TB = 2*tau(n,T_B) 

         Call nwpwxc_dldfcss(DTol,PB,GBB,TB,FB,FPB,FGB,FTB,EUB,
     &                ChiB,EUPB,ChiBP,ChiBG)
c        Ec = Ec + FB*qwght(n)          
         func(n)=func(n)+ FB*wght
         Amat(n,D1_RB)= Amat(n,D1_RB)+ FPB*wght
         Cmat(n,D1_GBB)=  Cmat(n,D1_GBB) + FGB*wght
         Mmat(n,D1_TB)=  Mmat(n,D1_TB) + FTB*wght
      endif

 30   continue
      P = PA + PB
      
      If((PA.gt.0.5d0*DTol).and.(PB.gt.0.5d0*DTol)) then
          RS = (Pi34/P) ** F13
          RSP = -RS/(F3*P)
          Zeta = (PA-PB)/P
          dZdA = (F1-Zeta)/P
          dZdB = (-F1-Zeta)/P
          Call nwpwxc_c_lsda(dtol,
     D         RS,Zeta,PotLC,dLdS,dLdZ,d2LdSS,d2LdSZ,
     $      d2LdZZ)
          EUEG = P*PotLC - EUA - EUB
          U = COpp*(ChiA+ChiB)/(F1 + COpp*(ChiA+ChiB))
          W = sopp0+U*(sopp1+U*(sopp2+U*(sopp3+U*sopp4)))
c         Ec = Ec + sop*EUEG*W*qwght(n)
          func(n)=func(n)+ sop*EUEG*W*wght
          dUdChiA =COpp/(F1 + COpp*(ChiA+ChiB))**2
          dUdChiB =COpp/(F1 + COpp*(ChiA+ChiB))**2
          dUdPA= dUdChiA*ChiAP
          dUdPB= dUdChiB*ChiBP
          dUdGA= dUdChiA*ChiAG
          dUdGB= dUdChiB*ChiBG
          dWdU =sopp1+U*(F2*sopp2+U*(F3*sopp3+U*F4*sopp4))
          dWdPA= dWdU*dUdPA
          dWdPB= dWdU*dUdPB
          dWdGA= dWdU*dUdGA
          dWdGB= dWdU*dUdGB
          EUEGPA = PotLC + P*dLdS*RSP + P*dLdZ*dZdA - EUPA
          EUEGPB = PotLC + P*dLdS*RSP + P*dLdZ*dZdB - EUPB
          if (ipol.eq.1) then 
            Amat(n,D1_RA) = Amat(n,D1_RA)
     &                    + sop*(EUEGPA*W + EUEG*dWdPA)*wght
            Cmat(n,D1_GAA)=  Cmat(n,D1_GAA)
     &                    + sop*(EUEG*dWdGA)*wght
          else
            Amat(n,D1_RA) = Amat(n,D1_RA)
     &                    + sop*(EUEGPA*W + EUEG*dWdPA)*wght
            Amat(n,D1_RB) = Amat(n,D1_RB)
     &                    + sop*(EUEGPB*W + EUEG*dWdPB)*wght
            Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + sop*EUEG*dWdGA*wght
            Cmat(n,D1_GBB) = Cmat(n,D1_GBB) + sop*(EUEG*dWdGB)*wght
          endif
      endIf
c      write (*,*) "PA, PB, GAA, GBB,ipol",PA, PB, GAA, GBB,ipol
c      write (*,*) "FA, FB,FGA, FGB",FA, FB,FGA, FGB
c      Stop
20    continue
      end


      Subroutine nwpwxc_dldfcss(DTol,PX,GX,TX,F,FP,FG,FT,EUEG,Chi,EUEGP,
     &                   ChiP,ChiG)
      Implicit none
C
C     Compute the same-spin part of the dldf correlation functional for one grid
C     point and one spin-case.
C
C
      double precision PX, GX, TX, F, FP, FG, FT, DTol
      double precision EUEG, Chi, EUEGP, ChiP, ChiG
      double precision Zero, Pt25, F1, F2, F3, F4, F5, F6, F8, F11
      double precision ss, sss0,sss1, sss2, sss3, sss4, Css
      double precision Pi, Pi34, F13, F23, F43, F53, F83, F113
      double precision RS, FDUEG, D, Fscc, RSP, dFsccP, dFsccG
      double precision E, W, U, dFsccT, dUdChi, dWdU, dWdP, dWdG 
      double precision d2LdSS,d2LdSZ,d2LdZZ,PotLC,dLdS,dLdZ
     


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
        Call nwpwxc_c_lsda(dtol,
     D       RS,F1,PotLC,dLdS,dLdZ,d2LdSS,d2LdSZ,d2LdZZ)
        EUEG = PX*PotLC
        D = TX - Pt25*GX/PX
C        DUEG = FDUEG*PX**F53
        Chi = GX/PX**F83
        U = Css*Chi/(F1 + Css*Chi)
        W = sss0+U*(sss1+U*(sss2+U*(sss3+U*sss4)))
        Fscc=D/TX
        E = Fscc*W*EUEG
        F = E*ss
        RSP = -RS/(F3*Px)
        ChiG = F1/PX**F83
        ChiP = -F83*Chi/PX
        dFsccP=Pt25*GX/(TX*PX**2)
        dFsccG=-Pt25/(TX*PX)
        dFsccT=Pt25*GX/(PX*TX**2)
        dUdChi=Css/((F1+Css*Chi)**2)
        dWdU=sss1+U*(F2*sss2+U*(F3*sss3+U*F4*sss4))
        dWdP=dWdU*dUdChi*ChiP
        dWdG=dWdU*dUdChi*ChiG 
        EUEGP = PotLC + PX*dLdS*RSP
        FP = ss*(dFsccP*W*EUEG 
     $                 + Fscc*dWdP*EUEG
     $                 + Fscc*W*EUEGP)
        FG = ss*(dFsccG*W*EUEG
     $                 + Fscc*dWdG*EUEG)

        FT = ss*(dFsccT*W*EUEG)
       Endif

       Return
       End


      Subroutine nwpwxc_c_dldf_d2()
      implicit none
      call errquit(' cdldf: d2 not coded ',0,0)
      return
      end

C> @}
