#include "dft2drv.fh"
c    Bc95 correlation functional          
c           META GGA
C         utilizes ingredients:
c                              rho   -  density
c                              delrho - gradient of density
c                              tau (tauN)- K.S kinetic energy density

      Subroutine xc_bc95(tol_rho, cfac, lcfac, nlcfac, rho, delrho,  
     &                     nq, ipol, Ec, qwght, ldew, func,
     &                     tau, Amat, Cmat, Mmat,ijmswitch)


c
c$Id$
c
c  Reference
c    Becke, A. D. J. Chem. Phys. 1996, 104, 1040.
c
      implicit none
c
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
      double precision DTol,F1, F2, F3, F4,COpp 
      Data F1/1.0d0/,F2/2.0d0/,
     & F3/3.0d0/,F4/4.0d0/ 
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
      double precision delrho(nq,3,ipol), gammaval, gam12
      
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

      integer n, ijmswitch

c    call to the bc95css subroutine
      double precision PA,GAA,TA,FA,FPA,FGA,FTA,EUA,EUEGA,ChiA,EUPA
     &,ChiAP,ChiAG
      double precision PB,GBB,TB,FB,FPB,FGB,FTB,EUB,EUEGB,ChiB,EUPB
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
      if (ijmswitch.eq.1) then
C     Parameters for PW6B95 Correlation
        COpp=0.00262d0 
      elseif (ijmswitch.eq.2) then
C     Parameters for PWB6K Correlation
        COpp=0.00353d0
      else
C     Parameters for B95 (data statement is unsafe as values may survive
C     between calls)
        COpp=0.0031d0
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
!         TA=0.0005d0
                  
         Call bc95ss(PA,GAA,TA,FA,FPA,FGA,FTA,EUA,
     &                ChiA,EUPA,ChiAP,ChiAG,ijmswitch)
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

         Ec = Ec + 2.0d0*FA*qwght(n)            !factor of 2 account for both spin
         if(ldew) func(n)=func(n)+ 2.0d0*FA
         Amat(n,1)=Amat(n,1)+ FPA
         Cmat(n,D1_GAA)=  Cmat(n,D1_GAA) + FGA
         Mmat(n,1)=  Mmat(n,1) + FTA
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
         
         PA = rho(n,2)
         if (PA.le.DTol) go to 25
         GAA =   delrho(n,1,1)*delrho(n,1,1) +
     &           delrho(n,2,1)*delrho(n,2,1) +
     &          delrho(n,3,1)*delrho(n,3,1)
c
c  In the bc95css subroutine, we use 2*TA as the tau 
c
         TA = tau(n,1)*2.0d0

         Call bc95ss(PA,GAA,TA,FA,FPA,FGA,FTA,EUA,
     &                ChiA,EUPA,ChiAP,ChiAG,ijmswitch)
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
c  In the bc95css subroutine, we use 2*TA as the tau, 
c
c
c       Beta 
c
 25       continue
         PB = rho(n,3)
         if(PB.le.DTol) go to 30
         GBB =   delrho(n,1,2)*delrho(n,1,2) +
     &           delrho(n,2,2)*delrho(n,2,2) +
     &          delrho(n,3,2)*delrho(n,3,2)

         TB = tau(n,2)*2.0d0

         Call bc95ss(PB,GBB,TB,FB,FPB,FGB,FTB,EUB,
     &                ChiB,EUPB,ChiBP,ChiBG,ijmswitch)
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
          Call lsdac(dtol,
     D         RS,Zeta,PotLC,dLdS,dLdZ,d2LdSS,d2LdSZ,
     $      d2LdZZ)
          EUEG = P*PotLC - EUA - EUB
          Denom = F1 + COpp*(ChiA+ChiB)
          Ec = Ec + sop*EUEG*qwght(n)/Denom
          if(ldew) func(n)=func(n)+ sop*EUEG/Denom
          DenPA = COpp*ChiAP
          DenPB = COpp*ChiBP
          DenGA = COpp*ChiAG
          DenGB = COpp*ChiBG
          EUEGPA = PotLC + P*dLdS*RSP + P*dLdZ*dZdA - EUPA
          EUEGPB = PotLC + P*dLdS*RSP + P*dLdZ*dZdB - EUPB
          if (ipol.eq.1) then 
            Amat(n,1) = Amat(n,1) + 
     &               sop*(EUEGPA/Denom - EUEG*DenPA/Denom**2)
            Cmat(n,D1_GAA)=  Cmat(n,D1_GAA) - sop*(EUEG*DenGA/Denom**2)
          else
            Amat(n,1) = Amat(n,1) + 
     &                 sop*(EUEGPA/Denom - EUEG*DenPA/Denom**2) 
            Amat(n,2) = Amat(n,2) +
     &                 sop*(EUEGPB/Denom - EUEG*DenPB/Denom**2) 
            Cmat(n,D1_GAA) = Cmat(n,D1_GAA) - sop*EUEG*DenGA/Denom**2
            Cmat(n,D1_GBB) = Cmat(n,D1_GBB) - sop*EUEG*DenGB/Denom**2
          endif
      endIf
c      write (*,*) "Amat(n,1),Cmat(n,1),Mmat(n,1)",Amat(n,1),Cmat(n,1)
c     & ,Mmat(n,1)
c      stop
20    continue
      end

      Subroutine xc_bc95_d2()
      call errquit(' bc95: d2 not coded ',0,0)
      return
      end




      Subroutine bc95ss(PX,GX,TX,F,FP,FG,FT,EUEG,Chi,EUEGP,
     &                   ChiP,ChiG,ijmswitch)
      Implicit none
C
C     Compute the same-spin part of the bc95 correlation functional for one grid
C     point and one spin-case.
C
C
      integer ijmswitch
      double precision PX, GX, TX, F, FP, FG, FT, DTol
      double precision EUEG, Chi, EUEGP, ChiP, ChiG,Css
      double precision Zero, Pt25, F1, F2, F3, F4, F5, F6, F8, F11
      double precision Pi, Pi34, F13, F23, F43, F53, F83, F113
      double precision RS, FDUEG, D,  RSP,DUEG, Denom, PotLC
      double precision E, DenomG, DenomP, DUEGP, DP, DG, DT  
      double precision d2LdSS,d2LdSZ,d2LdZZ,dLdS,dLdZ
     


      Data Zero/0.0d0/, Pt25/0.25d0/, F1/1.0d0/, F2/2.0d0/, F3/3.0d0/,
     $  F4/4.0d0/, F5/5.0d0/, F6/6.0d0/, F8/8.0d0/, F11/11.0d0/
C
      if (ijmswitch.eq.1) then
C     Parameters for PW6B95 Correlation
       Css=0.03668d0
      elseif (ijmswitch.eq.2) then
C     Parameters for PWB6K Correlation
       Css=0.04120d0
      else
C     Parameters for B95 (data statement is unsafe as values may survive
C     between calls)
       Css=0.038d0
      endif
      DTol =1.0d-6 
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
        Call lsdac(dtol,
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


