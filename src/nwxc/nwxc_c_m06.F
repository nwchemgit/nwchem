#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if !defined(NWAD_PRINT)
C> \ingroup nwxc
C> @{
C>
C> \file nwxc_c_m06.F
C> Implementation of the M06 correlation functional
C>
C> @}
#endif
#endif
C>
C> \ingroup nwxc_priv
C> @{
C>
C> \brief The M06 correlation functional
C>
C> The M06 functional [1,2] is a meta-GGA of which this evaluates
C> the correlation component.
C>
C> Due to the form of the meta-GGAs we need to screen on the kinetic
C> energy density to ensure that LDA will be obtained when the kinetic
C> energy density goes to zero [3].
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
C> [3] J. Gr&auml;fenstein, D. Izotov, D. Cremer,
C>     "Avoiding singularity problems associated with meta-GGA exchange
C>     and correlation functionals containing the kinetic energy
C>     density", J. Chem. Phys. <b>127</b>, 214103 (2007), DOI:
C>     <a href="https://doi.org/10.1063/1.2800011">
C>     10.1063/1.2800011</a>.
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

#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
      Subroutine nwxc_c_m06_p(param, tol_rho, ipol, nq, wght, rho,
     &                        rgamma, tau, func)
#else
      Subroutine nwxc_c_m06(param, tol_rho, ipol, nq, wght, rho, rgamma,
     &                      tau, func)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      Subroutine nwxc_c_m06_d2(param, tol_rho, ipol, nq, wght, rho,
     &                         rgamma, tau, func)
#else
      Subroutine nwxc_c_m06_d3(param, tol_rho, ipol, nq, wght, rho,
     &                         rgamma, tau, func)
#endif
c
c$Id$
c
c     [a]   Zhao, Y. and  Truhlar, D. G. J. Chem. Phys. 2006, 125, 194101;
c     [b]   Zhao, Y. and  Truhlar, D. G. J. Phys. Chem. A (2006),110(49),13126-13130.
c
#include "nwad.fh"
c
      implicit none
c
#include "intf_nwxc_c_lsda.fh"
#include "intf_nwxc_c_vs98.fh"
#include "intf_nwxc_m06css.fh"
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
      double precision F1, F2, F3, F4,COpp 
      Data COpp/0.0031d0/,F1/1.0d0/,F2/2.0d0/,
     & F3/3.0d0/,F4/4.0d0/ 

      integer n

c    call to the m06css subroutine
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

#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
      call nwxc_c_vs98_p(param, tol_rho, ipol, nq, wght, rho,
     &                 rgamma, tau, func)
#else
      call nwxc_c_vs98(param, tol_rho, ipol, nq, wght, rho,
     &                 rgamma, tau, func)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      call nwxc_c_vs98_d2(param, tol_rho, ipol, nq, wght, rho,
     &                    rgamma, tau, func)
#else
      call nwxc_c_vs98_d3(param, tol_rho, ipol, nq, wght, rho,
     &                    rgamma, tau, func)
#endif
      
      Pi = F4*ATan(F1)
      F6=6.0d0
      F43 = F4 / F3
      Pi34 = F3 / (F4*Pi)
      F13 = F1 / F3

      do 20 n = 1, nq
       EUA  = 0.0d0
       EUB  = 0.0d0
       ChiA = 0.0d0
       ChiB = 0.0d0
       if (ipol.eq.1) then
         if (rho(n,R_T).lt.Tol_Rho) goto 20
c
c    get the density, gradient, and tau for the alpha spin from the total 
c
         PA = rho(n,R_T)/F2
c        GAA = (    delrho(n,1,1)*delrho(n,1,1) +
c    &                 delrho(n,2,1)*delrho(n,2,1) +
c    &                 delrho(n,3,1)*delrho(n,3,1))/F4
         PB = PA
         GAA = rgamma(n,G_TT)/F4
c        if(sqrt(gaa).lt.tol_rho) goto 20
c  In the m06css subroutine, we use 2*TA as the tau, so we do not divide 
c  the tau by 2 here

         TA = tau(n,T_T)
         if(ta.lt.tol_rho) goto 30
                  
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
         Call nwxc_m06css_p(param,Tol_Rho,PA,GAA,TA,FA,FPA,FGA,FTA,EUA,
     &                ChiA,EUPA,ChiAP,ChiAG)
#else
         Call nwxc_m06css(param,Tol_Rho,PA,GAA,TA,FA,FPA,FGA,FTA,EUA,
     &                ChiA,EUPA,ChiAP,ChiAG)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
         Call nwxc_m06css_d2(param,Tol_Rho,PA,GAA,TA,FA,FPA,FGA,FTA,EUA,
     &                ChiA,EUPA,ChiAP,ChiAG)
#else
         Call nwxc_m06css_d3(param,Tol_Rho,PA,GAA,TA,FA,FPA,FGA,FTA,EUA,
     &                ChiA,EUPA,ChiAP,ChiAG)
#endif
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
         if (rho(n,R_A).le.0.5d0*Tol_Rho) go to 25
         PA = rho(n,R_A)
c        GAA =   delrho(n,1,1)*delrho(n,1,1) +
c    &           delrho(n,2,1)*delrho(n,2,1) +
c    &          delrho(n,3,1)*delrho(n,3,1)
c
c  In the m06css subroutine, we use 2*TA as the tau 
c
         if (tau(n,T_A).le.0.5d0*Tol_Rho) go to 25
         GAA = rgamma(n,G_AA)
         TA = 2.0d0*tau(n,T_A) 
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
         Call nwxc_m06css_p(param,Tol_Rho,PA,GAA,TA,FA,FPA,FGA,FTA,EUA,
     &                ChiA,EUPA,ChiAP,ChiAG)
#else
         Call nwxc_m06css(param,Tol_Rho,PA,GAA,TA,FA,FPA,FGA,FTA,EUA,
     &                ChiA,EUPA,ChiAP,ChiAG)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
         Call nwxc_m06css_d2(param,Tol_Rho,PA,GAA,TA,FA,FPA,FGA,FTA,EUA,
     &                ChiA,EUPA,ChiAP,ChiAG)
#else
         Call nwxc_m06css_d3(param,Tol_Rho,PA,GAA,TA,FA,FPA,FGA,FTA,EUA,
     &                ChiA,EUPA,ChiAP,ChiAG)
#endif
c        Ec = Ec + FA*qwght(n)     
         func(n)=func(n)+ FA*wght
c        Amat(n,D1_RA)  = Amat(n,D1_RA)+ FPA*wght
c        Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + FGA*wght
c        Mmat(n,D1_TA)  = Mmat(n,D1_TA) + FTA*wght
c
c  In the m06css subroutine, we use 2*TB as the tau, 
c
c
c       Beta 
c
 25      continue
         if (rho(n,R_B).le.0.5d0*Tol_Rho) go to 30
         PB = rho(n,R_B)
c        GBB =   delrho(n,1,2)*delrho(n,1,2) +
c    &           delrho(n,2,2)*delrho(n,2,2) +
c    &          delrho(n,3,2)*delrho(n,3,2)

         if (tau(n,T_B).le.0.5d0*Tol_Rho) go to 30
         GBB = rgamma(n,G_BB)
         TB = 2.0d0*tau(n,T_B)
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
         Call nwxc_m06css_p(param,Tol_Rho,PB,GBB,TB,FB,FPB,FGB,FTB,EUB,
     &                ChiB,EUPB,ChiBP,ChiBG)
#else
         Call nwxc_m06css(param,Tol_Rho,PB,GBB,TB,FB,FPB,FGB,FTB,EUB,
     &                ChiB,EUPB,ChiBP,ChiBG)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
         Call nwxc_m06css_d2(param,Tol_Rho,PB,GBB,TB,FB,FPB,FGB,FTB,EUB,
     &                ChiB,EUPB,ChiBP,ChiBG)
#else
         Call nwxc_m06css_d3(param,Tol_Rho,PB,GBB,TB,FB,FPB,FGB,FTB,EUB,
     &                ChiB,EUPB,ChiBP,ChiBG)
#endif
c        Ec = Ec + FB*qwght(n)          
         func(n)=func(n)+ FB*wght
c        Amat(n,D1_RB)  = Amat(n,D1_RB)+ FPB
c        Cmat(n,D1_GBB) = Cmat(n,D1_GBB) + FGB
c        Mmat(n,D1_TB)  = Mmat(n,D1_TB) + FTB
      endif

 30   continue
      P = PA + PB
      
      If((PA.gt.0.5d0*Tol_Rho).or.(PB.gt.0.5d0*Tol_Rho)) then
          RS = (Pi34/P) ** F13
c         RSP = -RS/(F3*P)
          Zeta = (PA-PB)/P
c         dZdA = (F1-Zeta)/P
c         dZdB = (-F1-Zeta)/P
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
          Call nwxc_c_lsda_p(tol_rho,
     R         RS,Zeta,PotLC,dLdS,dLdZ,d2LdSS,d2LdSZ,
     $         d2LdZZ)
#else
          Call nwxc_c_lsda(tol_rho,
     R         RS,Zeta,PotLC,dLdS,dLdZ,d2LdSS,d2LdSZ,
     $         d2LdZZ)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
          Call nwxc_c_lsda_d2(tol_rho,
     R         RS,Zeta,PotLC,dLdS,dLdZ,d2LdSS,d2LdSZ,
     $         d2LdZZ)
#else
          Call nwxc_c_lsda_d3(tol_rho,
     R         RS,Zeta,PotLC,dLdS,dLdZ,d2LdSS,d2LdSZ,
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
c           Amat(n,D1_RA)  = Amat(n,D1_RA) 
c    +                     + sop*(EUEGPA*W + EUEG*dWdPA)*wght
c           Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + sop*(EUEG*dWdGA)*wght
c         else
c           Amat(n,D1_RA)  = Amat(n,D1_RA) 
c    +                     + sop*(EUEGPA*W + EUEG*dWdPA)*wght
c           Amat(n,D1_RB)  = Amat(n,D1_RB) 
c    +                     + sop*(EUEGPB*W + EUEG*dWdPB)*wght
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
      Subroutine nwxc_m06css_p(param,Tol_Rho,PX,GX,TX,F,FP,FG,FT,EUEG,
     &                       Chi,EUEGP,ChiP,ChiG)
#else
      Subroutine nwxc_m06css(param,Tol_Rho,PX,GX,TX,F,FP,FG,FT,EUEG,
     &                       Chi,EUEGP,ChiP,ChiG)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      Subroutine nwxc_m06css_d2(param,Tol_Rho,PX,GX,TX,F,FP,FG,FT,EUEG,
     &                       Chi,EUEGP,ChiP,ChiG)
#else
      Subroutine nwxc_m06css_d3(param,Tol_Rho,PX,GX,TX,F,FP,FG,FT,EUEG,
     &                       Chi,EUEGP,ChiP,ChiG)
#endif
#include "nwad.fh"
      Implicit none
c
#include "intf_nwxc_c_lsda.fh"
c
C
C     Compute the same-spin part of the m06 correlation functional for one grid
C     point and one spin-case.
C
C
#if defined(NWAD_PRINT)
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      type(nwad_dble)::param(22)
      type(nwad_dble)::sss0, sss1, sss2, sss3, sss4
#else
      double precision param(22)
      double precision sss0, sss1, sss2, sss3, sss4
#endif
#else
      double precision param(22)
      double precision sss0, sss1, sss2, sss3, sss4
#endif
      type(nwad_dble)::PX, GX, TX, F, EUEG, Fscc, Chi
      type(nwad_dble)::Rs, D, E, W, U, Zeta, PotLC
      double precision FP, FG, FT, Tol_Rho
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
        Zeta = F1
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
        Call nwxc_c_lsda_p(tol_rho,
     R       RS,Zeta,PotLC,dLdS,dLdZ,d2LdSS,d2LdSZ,d2LdZZ)
#else
        Call nwxc_c_lsda(tol_rho,
     R       RS,Zeta,PotLC,dLdS,dLdZ,d2LdSS,d2LdSZ,d2LdZZ)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
        Call nwxc_c_lsda_d2(tol_rho,
     R       RS,Zeta,PotLC,dLdS,dLdZ,d2LdSS,d2LdSZ,d2LdZZ)
#else
        Call nwxc_c_lsda_d3(tol_rho,
     R       RS,Zeta,PotLC,dLdS,dLdZ,d2LdSS,d2LdSZ,d2LdZZ)
#endif
        EUEG = PX*PotLC
        D = TX - Pt25*GX/PX
C        DUEG = FDUEG*PX**F53
        Chi = GX/(PX**F83)
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
#include "nwxc_c_m06.F"
#endif
#ifndef SECOND_DERIV
#define SECOND_DERIV
c
c     Compile source again for the 2nd derivative case
c
#include "nwxc_c_m06.F"
#endif
#ifndef THIRD_DERIV
#define THIRD_DERIV
c
c     Compile source again for the 3rd derivative case
c
#include "nwxc_c_m06.F"
#endif
#undef NWAD_PRINT
C> @}
