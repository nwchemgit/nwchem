C> \ingroup nwpwxc
C> @{
C>
C> \file nwpwxc_c_m06.F
C> Implementation of the M06 correlation functional
C>
C> @}
C>
C> \ingroup nwpwxc_priv
C> @{
C>
C> \brief The M06 correlation functional
C>
C> The M06 functional [1,2] is a meta-GGA of which this evaluates
C> the correlation component.
C>
C> ### References ###
C>
C> [1] Y Zhao, DG Truhlar, 
C> "A new local density functional for main-group thermochemistry,
C> transition metal bonding, thermochemical kinetics, and noncovalent
C> interactions",
C> J. Chem. Phys. <b>125</b>, 194101 (2006), DOI:
C> <a href="https://doi.org/10.1063/1.2370993">
C> 10.1063/1.2370993</a>.
C>
C> [2] Y Zhao, DG Truhlar,
C> "Density functional for spectroscopy: No long-range self-interaction
C> error, good performance for Rydberg and charge-transfer states,
C> and better performance on average than B3LYP for ground states",
C> J. Phys. Chem. A <b>110</b>, 13126-13130 (2006), DOI:
C> <a href="https://doi.org/10.1021/jp066479k">
C> 10.1021/jp066479k</a>.
C>
c    M06 suite correlation functional 
c           META GGA
C         utilizes ingredients:
c                              rho   -  density
c                              delrho - gradient of density
c                              tau (tauN)- K.S kinetic energy density
c                              ijzy - 1  M06-L  
c                              ijzy - 2  M06-HF
c                              ijzy - 3  M06
c                              ijzy - 4  M06-2X

      Subroutine nwpwxc_c_m06(param,tol_rho,ipol,nq,wght,rho,rgamma,
     &                      tau, func, Amat, Cmat, Mmat)
c
c$Id$
c
c     [a]   Zhao, Y. and  Truhlar, D. G. J. Chem. Phys. 2006, 125, 194101;
c     [b]   Zhao, Y. and  Truhlar, D. G. J. Phys. Chem. A (2006),110(49),13126-13130.
c
      implicit none
c
#include "nwpwxc_param.fh"
c
c     Input and other parameters
c
      double precision param(*) !< [Input] Parameters of functional
                                !< - param(1): \f$ d_{C\alpha\beta,0} \f$
                                !< - param(2): \f$ d_{C\alpha\beta,1} \f$
                                !< - param(3): \f$ d_{C\alpha\beta,2} \f$
                                !< - param(4): \f$ d_{C\alpha\beta,3} \f$
                                !< - param(5): \f$ d_{C\alpha\beta,4} \f$
                                !< - param(6): \f$ d_{C\alpha\beta,5} \f$
                                !< - param(7): \f$ d_{C\sigma\sigma,0} \f$
                                !< - param(8): \f$ d_{C\sigma\sigma,1} \f$
                                !< - param(9): \f$ d_{C\sigma\sigma,2} \f$
                                !< - param(10): \f$ d_{C\sigma\sigma,3} \f$
                                !< - param(11): \f$ d_{C\sigma\sigma,4} \f$
                                !< - param(12): \f$ d_{C\sigma\sigma,5} \f$
                                !< - param(13): \f$ c_{C\alpha\beta,0} \f$
                                !< - param(14): \f$ c_{C\alpha\beta,1} \f$
                                !< - param(15): \f$ c_{C\alpha\beta,2} \f$
                                !< - param(16): \f$ c_{C\alpha\beta,3} \f$
                                !< - param(17): \f$ c_{C\alpha\beta,4} \f$
                                !< - param(18): \f$ c_{C\sigma\sigma,0} \f$
                                !< - param(19): \f$ c_{C\sigma\sigma,1} \f$
                                !< - param(20): \f$ c_{C\sigma\sigma,2} \f$
                                !< - param(21): \f$ c_{C\sigma\sigma,3} \f$
                                !< - param(22): \f$ c_{C\sigma\sigma,4} \f$
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
      double precision F1, F2, F3, F4,COpp 
      Data COpp/0.0031d0/,F1/1.0d0/,F2/2.0d0/,
     & F3/3.0d0/,F4/4.0d0/ 

      integer n

c    call to the m06css subroutine
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
      sop=1.0d0
      sopp0= param(13)
      sopp1= param(14)
      sopp2= param(15)
      sopp3= param(16)
      sopp4= param(17)
c     if (ijzy.eq.1) then
C     Parameters for M06-L Correlation
c        sopp0= 6.042374D-01
c        sopp1= 1.776783D+02
c        sopp2= -2.513252D+02
c        sopp3= 7.635173D+01
c        sopp4= -1.255699D+01
c     elseif (ijzy.eq.2) then
c     Parameters for M06-HF Correlation
c        sopp0= 1.674634D+00
c        sopp1= 5.732017D+01
c        sopp2= 5.955416D+01
c        sopp3= -2.311007D+02
c        sopp4= 1.255199D+02
c     elseif (ijzy.eq.3) then
c     Parameters for M06 Correlation
c        sopp0= 3.741539D+00
c        sopp1= 2.187098D+02
c        sopp2= -4.531252D+02
c        sopp3= 2.936479D+02
c        sopp4= -6.287470D+01
c     elseif (ijzy.eq.4) then
c     Parameters for M06-2X Correlation
c        sopp0= 8.833596D-01
c        sopp1= 3.357972D+01
c        sopp2= -7.043548D+01
c        sopp3= 4.978271D+01
c        sopp4= -1.852891D+01
c     endif

      call nwpwxc_c_vs98(param, tol_rho, ipol, nq, wght, rho,
     &                 rgamma, tau, func, Amat, Cmat, Mmat)
      
      Pi = F4*ATan(F1)
      F6=6.0d0
      F43 = F4 / F3
      Pi34 = F3 / (F4*Pi)
      F13 = F1 / F3

      do 20 n = 1, nq
       if (rho(n,R_T).lt.Tol_Rho) goto 20
       if (ipol.eq.1) then
c
c    get the density, gradient, and tau for the alpha spin from the total 
c
         PA = rho(n,R_T)/F2
c        GAA = (    delrho(n,1,1)*delrho(n,1,1) +
c    &                 delrho(n,2,1)*delrho(n,2,1) +
c    &                 delrho(n,3,1)*delrho(n,3,1))/F4
         GAA = rgamma(n,G_TT)/F4
         if(sqrt(gaa).lt.tol_rho) goto 20
c  In the m06css subroutine, we use 2*TA as the tau, so we do not divide 
c  the tau by 2 here

         TA = tau(n,T_T)
         if(ta.lt.tol_rho) goto 20
                  
         Call nwpwxc_m06css(param,Tol_Rho,PA,GAA,TA,FA,FPA,FGA,FTA,EUA,
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
         if (PA.le.0.5d0*Tol_Rho) go to 25
c        GAA =   delrho(n,1,1)*delrho(n,1,1) +
c    &           delrho(n,2,1)*delrho(n,2,1) +
c    &          delrho(n,3,1)*delrho(n,3,1)
         GAA = rgamma(n,G_AA)
c
c  In the m06css subroutine, we use 2*TA as the tau 
c
         TA = 2*tau(n,T_A) 

         Call nwpwxc_m06css(param,Tol_Rho,PA,GAA,TA,FA,FPA,FGA,FTA,EUA,
     &                ChiA,EUPA,ChiAP,ChiAG)
c        Ec = Ec + FA*qwght(n)     
         func(n)=func(n)+ FA*wght
         Amat(n,D1_RA)  = Amat(n,D1_RA)+ FPA*wght
         Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + FGA*wght
         Mmat(n,D1_TA)  = Mmat(n,D1_TA) + FTA*wght
c
c  In the m06css subroutine, we use 2*TB as the tau, 
c
c
c       Beta 
c
 25      continue
         PB = rho(n,R_B)
         if (PB.le.0.5d0*Tol_Rho) go to 30
c        GBB =   delrho(n,1,2)*delrho(n,1,2) +
c    &           delrho(n,2,2)*delrho(n,2,2) +
c    &          delrho(n,3,2)*delrho(n,3,2)
         GBB = rgamma(n,G_BB)

         TB = 2*tau(n,T_B)

         Call nwpwxc_m06css(param,Tol_Rho,PB,GBB,TB,FB,FPB,FGB,FTB,EUB,
     &                ChiB,EUPB,ChiBP,ChiBG)
c        Ec = Ec + FB*qwght(n)          
         func(n)=func(n)+ FB
         Amat(n,D1_RB)  = Amat(n,D1_RB)+ FPB
         Cmat(n,D1_GBB) = Cmat(n,D1_GBB) + FGB
         Mmat(n,D1_TB)  = Mmat(n,D1_TB) + FTB
      endif

 30   continue
      P = PA + PB
      
      If((PA.gt.0.5d0*Tol_Rho).and.(PB.gt.0.5d0*Tol_Rho)) then
          RS = (Pi34/P) ** F13
          RSP = -RS/(F3*P)
          Zeta = (PA-PB)/P
          dZdA = (F1-Zeta)/P
          dZdB = (-F1-Zeta)/P
          Call nwpwxc_c_lsda(tol_rho,
     R         RS,Zeta,PotLC,dLdS,dLdZ,d2LdSS,d2LdSZ,
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
            Amat(n,D1_RA)  = Amat(n,D1_RA) 
     +                     + sop*(EUEGPA*W + EUEG*dWdPA)*wght
            Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + sop*(EUEG*dWdGA)*wght
          else
            Amat(n,D1_RA)  = Amat(n,D1_RA) 
     +                     + sop*(EUEGPA*W + EUEG*dWdPA)*wght
            Amat(n,D1_RB)  = Amat(n,D1_RB) 
     +                     + sop*(EUEGPB*W + EUEG*dWdPB)*wght
            Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + sop*EUEG*dWdGA*wght
            Cmat(n,D1_GBB) = Cmat(n,D1_GBB) + sop*(EUEG*dWdGB)*wght
          endif
      endIf
c      write (*,*) "PA, PB, GAA, GBB,ipol",PA, PB, GAA, GBB,ipol
c      write (*,*) "FA, FB,FGA, FGB",FA, FB,FGA, FGB
c      Stop
20    continue
      end

      Subroutine nwpwxc_c_m06_d2()
      implicit none
      call errquit(' xc06: d2 not coded ',0,0)
      return
      end




      Subroutine nwpwxc_m06css(param,Tol_Rho,PX,GX,TX,F,FP,FG,FT,EUEG,
     &                       Chi,EUEGP,ChiP,ChiG)
      Implicit none
C
C     Compute the same-spin part of the m06 correlation functional for one grid
C     point and one spin-case.
C
C
      double precision param(22)
      double precision PX, GX, TX, F, FP, FG, FT, Tol_Rho
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
c      Tol_Rho=1.0D-7
c      write (*,*) Tol_Rho
      ss=1.0
      sss0= param(18)
      sss1= param(19)
      sss2= param(20)
      sss3= param(21)
      sss4= param(22)
c     if (ijzy.eq.1) then
C     Parameters for M06-L Correlation
c        sss0=  5.349466D-01
c        sss1=  5.396620D-01
c        sss2=  -3.161217D+01
c        sss3=  5.149592D+01
c        sss4=  -2.919613D+01
c     elseif (ijzy.eq.2) then
c     Parameters for M06-HF Correlation
c        sss0=  1.023254D-01
c        sss1=  -2.453783D+00
c        sss2=  2.913180D+01
c        sss3=  -3.494358D+01
c        sss4=  2.315955D+01
c     elseif (ijzy.eq.3) then
c     Parameters for M06 Correlation
c        sss0=  5.094055D-01
c        sss1=  -1.491085D+00
c        sss2=  1.723922D+01
c        sss3=  -3.859018D+01
c        sss4=  2.845044D+01
c     elseif (ijzy.eq.4) then
c     Parameters for M06-2X Correlation
c        sss0=  3.097855D-01
c        sss1=  -5.528642D+00
c        sss2=  1.347420D+01
c        sss3=  -3.213623D+01
c        sss4=  2.846742D+01
c     endif
      
      If ((PX.le.Tol_Rho))  then
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
        Call nwpwxc_c_lsda(tol_rho,
     R       RS,F1,PotLC,dLdS,dLdZ,d2LdSS,d2LdSZ,d2LdZZ)
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

C> @}
