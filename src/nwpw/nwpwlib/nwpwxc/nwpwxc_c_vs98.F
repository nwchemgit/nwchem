C> \ingroup nwpwxc
C> @{
C>
C> \file nwpwxc_c_vs98.F
C> Implementation of the VS98 correlation functional
C>
C> @}
C>
C> \ingroup nwpwxc_priv
C> @{
C>
C> \brief The Voorhis and Scuseria correlation functional
C>
C> The VS98 functional [1,2] is a meta-GGA. This routine implements
C> the correlation component.
C>
C> [1] T van Voorhis, GE Scuseria,
C>     "A novel form for the exchange-correlation energy functional",
C>     J.Chem.Phys. <b>109</b>, 400-410 (1998), DOI:
C>     <a href="https://doi.org/10.1063/1.476577">
C>     10.1063/1.476577</a>.
C>
C> [2] T van Voorhis, GE Scuseria,
C>     Erratum: "A novel form for the exchange-correlation energy
C>     functional",
C>     J.Chem.Phys. <b>129</b>, 219901-219901 (2008), DOI:
C>     <a href="https://doi.org/10.1063/1.3005348">
C>     10.1063/1.3005348</a>.
C>
c    VSXC correlation functional          
c           META GGA
C         utilizes ingredients:
c                              rho   -  density
c                              delrho - gradient of density
c                              tau (tauN)- K.S kinetic energy density
c                              ijzy - 1  VS98
c                              ijzy - 2  M06-L  
c                              ijzy - 3  M06-HF
c                              ijzy - 4  M06
c                              ijzy - 5  M06-2X
c
      Subroutine nwpwxc_c_vs98(param, tol_rho, ipol, nq, wght, rho,
     &                       rgamma, tau, func, Amat, Cmat, Mmat)
c
c$Id$
c
c  Reference
c   [a] T. V. Voorhis and G. E. Scuseria, J. Chem. Phys. 109, 400 (1998). 
c   [b] Y. Zhao and D. G. Truhlar, J. Chem. Phys. 125, 194101 (2006).
 
c
      implicit none
c
#include "nwpwxc_param.fh"
c
c     Input and other parameters
c
      double precision param(*) !< [Input] Parameters of functional (table 1)
                                !< - param(1): \f$ a_{\sigma\sigma'} \f$
                                !< - param(2): \f$ b_{\sigma\sigma'} \f$
                                !< - param(3): \f$ c_{\sigma\sigma'} \f$
                                !< - param(4): \f$ d_{\sigma\sigma'} \f$
                                !< - param(5): \f$ e_{\sigma\sigma'} \f$
                                !< - param(6): \f$ f_{\sigma\sigma'} \f$
                                !< - param(7): \f$ a_{\sigma\sigma} \f$
                                !< - param(8): \f$ b_{\sigma\sigma} \f$
                                !< - param(9): \f$ c_{\sigma\sigma} \f$
                                !< - param(10): \f$ d_{\sigma\sigma} \f$
                                !< - param(11): \f$ e_{\sigma\sigma} \f$
                                !< - param(12): \f$ f_{\sigma\sigma} \f$
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
      double precision DTol,F1, F2, F3, F4, gab, cf 
      Data F1/1.0d0/,F2/2.0d0/,
     & F3/3.0d0/,F4/4.0d0/,gab/0.00304966d0/,
     & cf/9.115599720d0/ 
c
c     Local
c 
      integer n
      double precision tauN

c    call to the vs98css subroutine
      double precision PA,GAA,TA,FA,FPA,FGA,FTA,EUA,EUEGA,ChiA,EUPA
     &,ChiAP,ChiAG,ZA,ZAP,ZAT
      double precision PB,GBB,TB,FB,FPB,FGB,FTB,EUB,EUEGB,ChiB,EUPB
     &,ChiBP,ChiBG,ZB,ZBP,ZBT
c
      double precision Pi, F43, F13, Pi34, F6, PotLC, 
     &RS,RSP,Zeta,dZdA,dZdB,dLdS,dLdZ,d2LdSS,d2LdSZ,d2LdZZ
      double precision P, EUEG, ZAB, XAB, kab, xk, zk
      double precision dgdx,dgdz,dgdPA,dgdGA,dgdTA,dgdPB,dgdGB,dgdTB
      double precision EUEGPA,EUEGPB,gcab
      double precision r7, r8, r9, r10, r11, r12
      double precision r13, r14, r15, r16, r17, r18

      
c
c     ======> BOTH SPIN-RESTRICETED AND UNRESTRICTED <======
c
c     DTol=1.0d-7
      dtol=tol_rho
C     Parameters for VS98 / M06-L / M06-HF / M06 / M06-2X
      r7=   param(1)
      r8=   param(2)
      r9=   param(3)
      r10=  param(4)
      r11=  param(5)
      r12=  param(6)
      r13=  param(7)
      r14=  param(8)
      r15=  param(9)
      r16=  param(10)
      r17=  param(11)
      r18=  param(12)
c
      Pi = F4*ATan(F1)
      F6=6.0d0
      F43 = F4 / F3
      Pi34 = F3 / (F4*Pi)
      F13 = F1 / F3
c
      do 20 n = 1, nq
        if (ipol.eq.1) then
          if (rho(n,R_A).lt.DTol) goto 20
        else
          if (rho(n,R_A)+rho(n,R_B).lt.DTol) goto 20
        endif
        if (ipol.eq.1) then
c
c         get the density, gradient, and tau for the alpha spin
c         from the total 
c
          PA = rho(n,R_A)/F2
          PB = 0.0d0
          GAA = rgamma(n,G_AA)/4.0d0
          if(sqrt(gaa).lt.dtol) goto 20
c         In the bc95css subroutine, we use 2*TA as the tau, so we do
c         not divide the tau by 2 here

          TA = tau(n,T_A) 
          if(ta.lt.dtol) goto 20
                  
          Call nwpwxc_vs98ss(tol_rho,PA,GAA,TA,FA,FPA,FGA,FTA,EUA,ZA,
     &                     ChiA,EUPA,ChiAP,ChiAG,ZAP,ZAT,
     &                     r13,r14,r15,r16,r17,r18)
          PB = PA
          GBB = GAA
          TB = TA
          FB = FA
          FPB = FPA
          FGB = FGA
          FTB = FTA
          EUB = EUA
          ZB = ZA
          ChiB = ChiA
          EUPB = EUPA
          ChiBP = ChiAP
          ChiBG = ChiAG
          ZBP = ZAP
          ZBT = ZAT
 
          func(n)        = func(n)        + 2.0d0*FA*wght
          Amat(n,D1_RA)  = Amat(n,D1_RA)  + FPA*wght
          Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + FGA*wght
          Mmat(n,D1_TA)  = Mmat(n,D1_TA)  + FTA*wght
 
 
c UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUnrestricted
      else  ! ipol=2
c
c        ======> SPIN-UNRESTRICTED <======
c
c
c       alpha
c
         
         PA = rho(n,R_A)
         if (PA.le.DTol) go to 25
         GAA = rgamma(n,G_AA)
c
c        In the bc95css subroutine, we use 2*TA as the tau 
c
         TA = tau(n,T_A)*2.0d0

         Call nwpwxc_vs98ss(tol_rho,PA,GAA,TA,FA,FPA,FGA,FTA,EUA,ZA,
     &                    ChiA,EUPA,ChiAP,ChiAG,ZAP,ZAT,
     &                    r13,r14,r15,r16,r17,r18)
         func(n)        = func(n)        + FA*wght
         Amat(n,D1_RA)  = Amat(n,D1_RA)  + FPA*wght
         Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + FGA*wght
c        2*0.5=1.0 for Mmat
         Mmat(n,D1_TA)  = Mmat(n,D1_TA)  + FTA*wght

c
c        In the vs98ss subroutine, we use 2*TA as the tau, 
c
c
c       Beta 
c
 25      continue
         PB = rho(n,R_B)
         if(PB.le.DTol) go to 30
         GBB = rgamma(n,G_BB)
         TB = tau(n,T_B)*2.0d0

         Call nwpwxc_vs98ss(tol_rho,PB,GBB,TB,FB,FPB,FGB,FTB,EUB,ZB,
     &                    ChiB,EUPB,ChiBP,ChiBG,ZBP,ZBT,
     &                    r13,r14,r15,r16,r17,r18)
         func(n)        = func(n)        + FB*wght
         Amat(n,D1_RB)  = Amat(n,D1_RB)  + FPB*wght
         Cmat(n,D1_GBB) = Cmat(n,D1_GBB) + FGB*wght
         Mmat(n,D1_TB)  = Mmat(n,D1_TB)  + FTB*wght
         
#if 0
      write (0,'(A,3F20.6)') "BAmat Cmat Mmat",FPB,FGB,FTB
#endif
      endif
 30   continue
      P = PA + PB
      If(PA.gt.DTol.and.PB.gt.DTol) then
          RS = (Pi34/P) ** F13 
          RSP = -RS/(F3*P)
          Zeta = (PA-PB)/P
          dZdA = (F1-Zeta)/P
          dZdB = (-F1-Zeta)/P
          Call nwpwxc_c_lsda(tol_rho,RS,Zeta,PotLC,
     &               dLdS,dLdZ,d2LdSS,d2LdSZ,d2LdZZ)
          EUEG = P*PotLC - EUA - EUB
          ZAB = ZA + ZB
          XAB = ChiA+ChiB
          kab = F1 + gab*(XAB+ZAB)
          xk = XAB/kab
          zk = ZAB/kab
          call nwpwxc_gvt4(gcab,dgdx,dgdz,xk,zk,kab,gab,gab,
     &                   r7,r8,r9,r10,r11,r12)
          func(n) = func(n) + gcab*EUEG*wght
          dgdPA = dgdx*ChiAP + dgdz*ZAP
          dgdGA = dgdx*ChiAG
          dgdTA = dgdz*ZAT
          dgdPB = dgdx*ChiBP + dgdz*ZBP
          dgdGB = dgdx*ChiBG
          dgdTB = dgdz*ZBT
          EUEGPA = PotLC + P*dLdS*RSP + P*dLdZ*dZdA - EUPA
          EUEGPB = PotLC + P*dLdS*RSP + P*dLdZ*dZdB - EUPB
          if (ipol.eq.1) then 
            Amat(n,D1_RA)  = Amat(n,D1_RA) 
     &                     + (EUEGPA*gcab   + EUEG*dgdPA)*wght
            Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + EUEG*dgdGA*wght
            Mmat(n,D1_TA)  = Mmat(n,D1_TA)  + EUEG*dgdTA*wght
          else
            Amat(n,D1_RA)  = Amat(n,D1_RA)
     &                     + (EUEGPA*gcab   + EUEG*dgdPA)*wght
            Amat(n,D1_RB)  = Amat(n,D1_RB)
     &                     + (EUEGPB*gcab   + EUEG*dgdPB)*wght
            Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + EUEG*dgdGA*wght
            Cmat(n,D1_GBB) = Cmat(n,D1_GBB) + EUEG*dgdGB*wght
            Mmat(n,D1_TA)  = Mmat(n,D1_TA)  + EUEG*dgdTA*wght
            Mmat(n,D1_TB)  = Mmat(n,D1_TB)  + EUEG*dgdTB*wght
          endif
      endIf
c      write (*,*) "Amat(n,1),Cmat(n,1),Mmat(n,1)",Amat(n,1),Cmat(n,1)
c     & ,Mmat(n,1),ipol
c      stop
20    continue
      end

      Subroutine nwpwxc_c_vs98_d2()
      call errquit(' cvs98: d2 not coded ',0,0)
      return
      end
C>
C> \brief Compute the same-spin part of the VS98 correlation functional
C>
C> This routine evaluates the same-spin part of the VS98 functional for
C> 1 grid point and 1 spin-case.
C>
      Subroutine nwpwxc_vs98ss(tol_rho,PX,GX,TX,F,FP,FG,FT,EUEG,Z,Chi,
     +                       EUEGP,ChiP,ChiG,ZP,ZT,
     +                       r13,r14,r15,r16,r17,r18)
      Implicit none
C
C     Compute the same-spin part of the vs98 correlation functional for one grid
C     point and one spin-case.
C

      double precision tol_rho
      double precision r13, r14, r15, r16, r17, r18
      double precision PX, GX, TX, F, FP, FG, FT, DTol, Z, ZP, ZT
      double precision EUEG, Chi, EUEGP, ChiP, ChiG, cf, gcc
      double precision Zero, Pt25, F1, F2, F3, F4, F5, F6, F8, F11
      double precision Pi, Pi34, F13, F23, F43, F53, F83, F113
      double precision RS, D, RSP, PotLC, DX, DZ, dgdP, dgdG, dgdT
      double precision E,DP, DG, DT, rhoo, rho43, rho53, rho83
      double precision rrho, F4o3, rho13, kc, xk, zk, gc, dgdx, dgdz
      double precision d2LdSS, d2LdSZ, d2LdZZ, dLdS, dLdZ

      Data Zero/0.0d0/, Pt25/0.25d0/, F1/1.0d0/, F2/2.0d0/, F3/3.0d0/,
     $  F4/4.0d0/, F5/5.0d0/, F6/6.0d0/, F8/8.0d0/, F11/11.0d0/,
     $  gcc/0.00515088d0/,cf/9.115599720d0/
 
 
      F4o3 = 4.0d0/3.0d0
      dtol=tol_rho
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
        Z  = Zero
        ZP = Zero
        ZT = Zero
      else
        Pi = F4*ATan(F1)
        Pi34 = F3 / (F4*Pi)
        F13 = F1 / F3
        F23 = F2 / F3
        F43 = F2 * F23
        F53 = F5 / F3
        F83 = F8 / F3
        F113 = F11 / F3
        rhoo = PX 
        rrho = 1.0d0/rhoo
        rho43 = rhoo**F4o3
        rho13 = rho43*rrho
        rho53 = rhoo**F53
        rho83 = rho53*rhoo
        
        RS = (Pi34/PX) ** F13
        Call nwpwxc_c_lsda(tol_rho,
     A       RS,F1,PotLC,dLdS,dLdZ,d2LdSS,d2LdSZ,d2LdZZ)
        EUEG = PX*PotLC
        Chi = GX/rho83
        Z = (TX/rho53) - cf
        kc = F1 + gcc*(Chi + Z)
        xk = Chi/kc
        zk = Z/kc
        D = F1 - Chi/(F4*(Z + cf)) 
        call nwpwxc_gvt4(gc,dgdx,dgdz,xk,zk,kc,gcc,gcc,
     &                 r13,r14,r15,r16,r17,r18)
        E = D*EUEG*gc
c         write (*,*) "Chi, Z, gc", CHi, Z, gc
        F = E 
c
        RSP = -RS/(F3*Px)
        ChiG = F1/PX**F83
        ChiP = -F83*Chi/PX
        ZP = -F53 * TX/rho83
        ZT =  F1/rho53
        DZ = Chi/(F4*(Z + cf)*(Z + cf)) 
        DX = -F1/(F4*(Z + cf))
        DP = DZ*ZP + DX*ChiP
        DG = DX*ChiG
        DT = DZ*ZT
        dgdP = dgdx*ChiP + dgdz*ZP
        dgdG = dgdx*ChiG 
        dgdT = dgdz*ZT
        EUEGP = PotLC + PX*dLdS*RSP
        FP = DP*EUEG*gc + D*EUEGP*gc + D*EUEG*dgdP
        FG = DG*EUEG*gc + D*EUEG*dgdG
        FT = DT*EUEG*gc + D*EUEG*dgdT
       Endif
       Return
       End
C> @}
