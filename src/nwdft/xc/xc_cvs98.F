#include "dft2drv.fh"
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
c                              ijzy - 6  revM06
c                              ijzy - 7  revM06-L
c
      Subroutine xc_cvs98(tol_rho, cfac, lcfac, nlcfac, rho, delrho,  
     &                     nq, ipol, Ec, qwght, ldew, func,
     &                     tau, Amat, Cmat, Mmat, ijzy)


c
c$Id$
c
c  Reference
c   [a] T. V. Voorhis and G. E. Scuseria, J. Chem. Phys. 109, 400 (1998). 
c   [b] Y. Zhao and D. G. Truhlar, J. Chem. Phys. 125, 194101 (2006).
 
c
      implicit none
c
#include "errquit.fh"
c
c
c     Input and other parameters
c
      integer ipol, nq

      double precision cfac
      logical lcfac, nlcfac

      logical lfac, nlfac
      double precision fac
      double precision tol_rho
c
      logical ldew
      double precision func(*)
c
c     Threshold parameters
c
      double precision DTol,F1, F2, F3, F4, gab, cf 
      Data F1/1.0d0/,F2/2.0d0/,
     & F3/3.0d0/,F4/4.0d0/,gab/0.00304966d0/,
     & cf/9.115599720d0/ 
c
c     Correlation energy
c
      double precision Ec
c
c     Charge Density 
c
      double precision rho(nq,ipol*(ipol+1)/2)
c
c     Charge Density Gradient
c
      double precision delrho(nq,3,ipol)
      
c
c     Kinetic Energy Density
c
      double precision tau(nq,ipol), tauN
 
c
c     Quadrature Weights
c
      double precision qwght(nq)
c
c     Sampling Matrices for the XC Potential
c
      double precision Amat(nq,ipol), Cmat(nq,*)
      double precision Mmat(nq,*)

      integer n, ijzy

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

      
c
c     ======> BOTH SPIN-RESTRICETED AND UNRESTRICTED <======
c
c      DTol=1.0d-7
      dtol=tol_rho
C     Parameters for VS98 
      if (ijzy.eq.1) then
              r7=   7.035010d-01
              r8=   7.694574d-03
              r9=   5.152765d-02
              r10=   3.394308d-05
              r11=  -1.269420d-03
              r12=   1.296118d-03
C     Parameters for M06-L
      elseif (ijzy.eq.2) then
              r7=      3.957626D-01
              r8=      -5.614546D-01
              r9=      1.403963D-02
              r10=     9.831442D-04
              r11=     -3.577176D-03
              r12=     0.000000D+00
C     Parameters for M06-HF
      elseif (ijzy.eq.3) then
              r7=    -6.746338D-01
              r8=    -1.534002D-01
              r9=    -9.021521D-02
              r10=   -1.292037D-03
              r11=   -2.352983D-04
              r12=   0.000000D+00

C     Parameters for M06
      elseif (ijzy.eq.4) then
               r7= -2.741539D+00
               r8= -6.720113D-01
               r9= -7.932688D-02
               r10=1.918681D-03
               r11=-2.032902D-03
               r12=0.000000D+00

C     Parameters for M06-2X
      elseif (ijzy.eq.5) then
              r7=  1.166404D-01
              r8=  -9.120847D-02
              r9=  -6.726189D-02
              r10= 6.720580D-05
              r11= 8.448011D-04
              r12= 0.000000D+00
C     Parameters for revM06-L
      elseif (ijzy.eq.6) then
              r7=  4.007146D-01
              r8=  1.5796569D-02
              r9=  -3.2680984D-02
              r10= 0.000000D+00
              r11= 0.000000D+00
              r12= 1.260132D-03
C     Parameters for revM06
      elseif (ijzy.eq.7) then
              r7=  -3.390666720D-01
              r8=  3.790156384D-03
              r9=  -2.762485975D-02
              r10= 0.000000D+00
              r11= 0.000000D+00
              r12= 4.076285162D-04
      else
         call errquit("xc_cvs98: illegal value of ijzy",ijzy,UERR)
      endif
      Pi = F4*ATan(F1)
      F6=6.0d0
      F43 = F4 / F3
      Pi34 = F3 / (F4*Pi)
      F13 = F1 / F3


      do 20 n = 1, nq
       if (rho(n,1).lt.DTol) goto 20
       if (ipol.eq.1) then
c
c    get the density, gradient, and tau for the alpha spin from the total 
c
         PA = rho(n,1)/F2
         GAA = (    delrho(n,1,1)*delrho(n,1,1) +
     &                 delrho(n,2,1)*delrho(n,2,1) +
     &                 delrho(n,3,1)*delrho(n,3,1))/4.0d0
c  In the bc95css subroutine, we use 2*TA as the tau, so we do not divide 
c  the tau by 2 here

         TA = tau(n,1) 
                  
         Call vs98ss(tol_rho,PA,GAA,TA,FA,FPA,FGA,FTA,EUA,ZA,
     &                ChiA,EUPA,ChiAP,ChiAG,ZAP,ZAT,ijzy)
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

         Ec = Ec + 2.0d0*FA*qwght(n)            !factor of 2 account for both spin
         if(ldew) func(n)=func(n)+ 2.0d0*FA
         Amat(n,1)=Amat(n,1)+ FPA
         Cmat(n,D1_GAA)=  Cmat(n,D1_GAA) + FGA
         Mmat(n,1)=  Mmat(n,1) + FTA
c if 0
c       write (0,'(A,3F20.6)') " PA,EUA",PA,EUA
c endif
 
 
c UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUnrestricted
      else  ! ipol=2
c
c        ======> SPIN-UNRESTRICTED <======
c
c
c       alpha
c
         
         PA = rho(n,2)
         if (PA.le.DTol) go to 25
         GAA =   delrho(n,1,1)*delrho(n,1,1) +
     &           delrho(n,2,1)*delrho(n,2,1) +
     &          delrho(n,3,1)*delrho(n,3,1)
c
c  In the bc95css subroutine, we use 2*TA as the tau 
c
         TA = tau(n,1)*2.0d0

         Call vs98ss(tol_rho,PA,GAA,TA,FA,FPA,FGA,FTA,EUA,ZA,
     &                ChiA,EUPA,ChiAP,ChiAG,ZAP,ZAT,ijzy)
         Ec = Ec + FA*qwght(n)    
         if(ldew) func(n)=func(n)+ FA
         Amat(n,1)=Amat(n,1)+ FPA
         Cmat(n,D1_GAA)=  Cmat(n,D1_GAA) + FGA
c      2*0.5=1.0 for Mmat
         Mmat(n,1)=  Mmat(n,1) + FTA
#if 0
      write (0,'(A,3F20.6)') "AAmat Cmat Mmat",FPA,FGA,FTA
#endif

c
c  In the vs98ss subroutine, we use 2*TA as the tau, 
c
c
c       Beta 
c
 25       continue
         PB = rho(n,3)
         GBB =   delrho(n,1,2)*delrho(n,1,2) +
     &           delrho(n,2,2)*delrho(n,2,2) +
     &          delrho(n,3,2)*delrho(n,3,2)

         TB = tau(n,2)*2.0d0

         Call vs98ss(tol_rho,PB,GBB,TB,FB,FPB,FGB,FTB,EUB,ZB,
     &                ChiB,EUPB,ChiBP,ChiBG,ZBP,ZBT,ijzy)
         Ec = Ec + FB*qwght(n)          
         if(ldew) func(n)=func(n)+ FB
         Amat(n,2)= Amat(n,2)+ FPB
         Cmat(n,D1_GBB)=  Cmat(n,D1_GBB) + FGB
         Mmat(n,2)=  Mmat(n,2) + FTB
         
#if 0
      write (0,'(A,3F20.6)') "BAmat Cmat Mmat",FPB,FGB,FTB
#endif
      endif
 30   continue
      P = rho(n,1) 
      If(PA.gt.DTol.and.PB.gt.DTol) then
          RS = (Pi34/P) ** F13 
          RSP = -RS/(F3*P)
          Zeta = (PA-PB)/P
          dZdA = (F1-Zeta)/P
          dZdB = (-F1-Zeta)/P
          Call lsdac(tol_rho,
     A         RS,Zeta,PotLC,dLdS,dLdZ,d2LdSS,d2LdSZ,
     $      d2LdZZ)
          EUEG = P*PotLC - EUA - EUB
          ZAB = ZA + ZB
          XAB = ChiA+ChiB
          kab = F1 + gab*(XAB+ZAB)
          xk = XAB/kab
          zk = ZAB/kab
       call gvt4(gcab,dgdx,dgdz,xk,zk,kab,gab,gab,r7,r8,r9,r10,r11,r12)
          Ec = Ec + gcab*EUEG*qwght(n)
          if(ldew) func(n)=func(n)+ gcab*EUEG
          dgdPA = dgdx*ChiAP + dgdz*ZAP
          dgdGA = dgdx*ChiAG
          dgdTA = dgdz*ZAT
          dgdPB = dgdx*ChiBP + dgdz*ZBP
          dgdGB = dgdx*ChiBG
          dgdTB = dgdz*ZBT
          EUEGPA = PotLC + P*dLdS*RSP + P*dLdZ*dZdA - EUPA
          EUEGPB = PotLC + P*dLdS*RSP + P*dLdZ*dZdB - EUPB
          if (ipol.eq.1) then 
           Amat(n,1) = Amat(n,1) + (EUEGPA*gcab + EUEG*dgdPA)
           Cmat(n,D1_GAA)=  Cmat(n,D1_GAA) + EUEG*dgdGA
           Mmat(n,1)=  Mmat(n,1) + EUEG*dgdTA 
          else
            Amat(n,1) = Amat(n,1) + (EUEGPA*gcab + EUEG*dgdPA)
            Amat(n,2) = Amat(n,2) + (EUEGPB*gcab + EUEG*dgdPB)
            Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + EUEG*dgdGA 
            Cmat(n,D1_GBB) = Cmat(n,D1_GBB) + EUEG*dgdGB
            Mmat(n,1)=  Mmat(n,1) + EUEG*dgdTA
            Mmat(n,2)=  Mmat(n,2) + EUEG*dgdTB 
          endif
      endIf
c      write (*,*) "Amat(n,1),Cmat(n,1),Mmat(n,1)",Amat(n,1),Cmat(n,1)
c     & ,Mmat(n,1),ipol
c      stop
20    continue
      end

      Subroutine xc_cvs98_d2()
      call errquit(' cvs98: d2 not coded ',0,0)
      return
      end

      Subroutine vs98ss(tol_rho,PX,GX,TX,F,FP,FG,FT,EUEG,Z,Chi,EUEGP,
     &                   ChiP,ChiG,ZP,ZT,ijzy)
      Implicit none
c
#include "errquit.fh"
C
C     Compute the same-spin part of the vs98 correlation functional for one grid
C     point and one spin-case.
C

      integer ijzy
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
C     Parameters for VS98 
      if (ijzy.eq.1) then
              r13=   3.270912d-01
              r14=  -3.228915d-02
              r15=  -2.942406d-02
              r16=   2.134222d-03
              r17=  -5.451559d-03
              r18=   1.577575d-02
C     Parameters for M06-L
      elseif (ijzy.eq.2) then
              r13=   4.650534D-01
              r14=   1.617589D-01
              r15=   1.833657D-01
              r16=   4.692100D-04
              r17=  -4.990573D-03
              r18=   0.000000D+00
C     Parameters for M06-HF
      elseif (ijzy.eq.3) then
              r13=   8.976746D-01
              r14=  -2.345830D-01
              r15=   2.368173D-01
              r16=  -9.913890D-04
              r17=  -1.146165D-02
              r18=   0.000000D+00
C     Parameters for M06
      elseif (ijzy.eq.4) then
               r13=  4.905945D-01
               r14= -1.437348D-01
               r15=  2.357824D-01
               r16=  1.871015D-03
               r17= -3.788963D-03
               r18=  0.000000D+00
C     Parameters for M06-2X
      elseif (ijzy.eq.5) then
              r13=  6.902145D-01
              r14=  9.847204D-02
              r15=  2.214797D-01
              r16= -1.968264D-03
              r17= -6.775479D-03
              r18=  0.000000D+00
C     Parameters for revM06-L
      elseif (ijzy.eq.6) then
              r13=  -5.38821292D-01
              r14=  -2.829603D-02
              r15=  2.3889696D-02
              r16=  0.000000D+00
              r17=  0.000000D+00
              r18=  -2.437902D-03
C     Parameters for revM06
      elseif (ijzy.eq.7) then
              r13=  -1.467095900D-01
              r14=  -1.832187007D-04
              r15=  8.484372430D-02
              r16=  0.000000D+00
              r17=  0.000000D+00
              r18=  2.280677172D-04
      else
        call errquit("vs98ss: illegal value of ijzy",ijzy,UERR)
      endif

      dtol=tol_rho
      If(PX.le.DTol.or.gx.le.dtol.or.tx.le.dtol) then
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
        Call lsdac(tol_rho,
     A       RS,F1,PotLC,dLdS,dLdZ,d2LdSS,d2LdSZ,d2LdZZ)
        EUEG = PX*PotLC
        Chi = GX/rho83
        Z = (TX/rho53) - cf
        kc = F1 + gcc*(Chi + Z)
        xk = Chi/kc
        zk = Z/kc
        D = F1 - Chi/(F4*(Z + cf)) 
        call gvt4(gc,dgdx,dgdz,xk,zk,kc,gcc,gcc,r13,r14,r15,r16,r17,r18)
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


