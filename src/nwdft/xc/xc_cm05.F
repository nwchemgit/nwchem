c    M05 and M05-2X correlation functional 
c           META GGA
C         utilizes ingredients:
c                              rho   -  density
c                              delrho - gradient of density
c                              tau (tauN)- K.S kinetic energy density

      Subroutine xc_cm05(tol_rho, cfac, lcfac, nlcfac, rho, delrho, 
     &                     nq, ipol, Ec, qwght, ldew, func,
     &                     tau, Amat, Cmat, Mmat,ijzy)     
     

c
c$Id$
c
c
c     [a]       Zhao, Y.; Schultz, N. E.; Truhlar, D. G. J. Chem. Phys. 2005, 123, 161103;
c      Note that in this communication we interchanged cCab,i and cCss,i in Table 1. 
c     [b]       Zhao, Y.; Schultz, N. E.; Truhlar, D. G. J. Chem. Theory Comput. 2006, in press.    

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
c     Threshold parameters
c
      double precision DTol,F1, F2, F3, F4,COpp 
      Data COpp/0.0031d0/,F1/1.0d0/,F2/2.0d0/,
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
      logical ldew
      double precision func(*)
c
c     Sampling Matrices for the XC Potential
c
      double precision Amat(nq,ipol), Cmat(nq,*)
      double precision Mmat(nq,*)

      integer n, ijzy

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
      if (ijzy.eq.1) then
C     Parameters for M05 Correlation
         sopp0= 1.00000d0
         sopp1= 3.78569d0
         sopp2= -14.15261d0
         sopp3= -7.46589d0
         sopp4= 17.94491d0
      elseif (ijzy.eq.2) then
C     Parameters for M05-2X Correlation
         sopp0= 1.00000
         sopp1= 1.09297
         sopp2= -3.79171
         sopp3= 2.82810
         sopp4= -10.58909
      else
         call errquit("xc_cm05: illegal value of ijzy",ijzy,UERR)
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
     &                 delrho(n,3,1)*delrho(n,3,1))/4
c  In the m05css subroutine, we use 2*TA as the tau, so we do not divide 
c  the tau by 2 here

         TA = tau(n,1)
                  
         Call m05css(DTol,PA,GAA,TA,FA,FPA,FGA,FTA,EUA,
     &                ChiA,EUPA,ChiAP,ChiAG,ijzy)
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

         Ec = Ec + 2.d0*FA*qwght(n)            !factor of 2 account for both spin
         if(ldew) func(n)=func(n)+ FA*2d0
         Amat(n,1)=Amat(n,1)+ FPA
         Cmat(n,1)=  Cmat(n,1) + FGA
         Mmat(n,1)=  Mmat(n,1) + FTA
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

         PA = rho(n,2)
         if (PA.le.DTol) go to 25
         GAA =   delrho(n,1,1)*delrho(n,1,1) +
     &           delrho(n,2,1)*delrho(n,2,1) +
     &          delrho(n,3,1)*delrho(n,3,1)
c
c  In the m05css subroutine, we use 2*TA as the tau 
c
         TA = 2*tau(n,1) 

         Call m05css(DTol,PA,GAA,TA,FA,FPA,FGA,FTA,EUA,
     &                ChiA,EUPA,ChiAP,ChiAG,ijzy)
         Ec = Ec + FA*qwght(n)     
         if(ldew) func(n)=func(n)+ FA
         Amat(n,1)=Amat(n,1)+ FPA
         Cmat(n,1)=  Cmat(n,1) + FGA
         Mmat(n,1)=  Mmat(n,1) + FTA
c
c  In the m05css subroutine, we use 2*TB as the tau, 
c
c
c       Beta 
c
 25      continue
         PB = rho(n,3)
         if (PB.le.DTol) go to 30
         GBB =   delrho(n,1,2)*delrho(n,1,2) +
     &           delrho(n,2,2)*delrho(n,2,2) +
     &          delrho(n,3,2)*delrho(n,3,2)

         TB = 2*tau(n,2)

         Call m05css(DTol,PB,GBB,TB,FB,FPB,FGB,FTB,EUB,
     &                ChiB,EUPB,ChiBP,ChiBG,ijzy)
         Ec = Ec + FB*qwght(n)          
         if(ldew) func(n)=func(n)+ FB
         Amat(n,2)= Amat(n,2)+ FPB
         Cmat(n,3)=  Cmat(n,3) + FGB
         Mmat(n,2)=  Mmat(n,2) + FTB
      endif

 30   continue
      P = PA + PB
      
      If((PA.gt.DTol).and.(PB.gt.DTol)) then
          RS = (Pi34/P) ** F13
          RSP = -RS/(F3*P)
          Zeta = (PA-PB)/P
          dZdA = (F1-Zeta)/P
          dZdB = (-F1-Zeta)/P
          Call lsdac(dtol,
     D         RS,Zeta,PotLC,dLdS,dLdZ,d2LdSS,d2LdSZ,
     $      d2LdZZ)
          EUEG = P*PotLC - EUA - EUB
          U = COpp*(ChiA+ChiB)/(F1 + COpp*(ChiA+ChiB))
          W = sopp0+U*(sopp1+U*(sopp2+U*(sopp3+U*sopp4)))
          Ec = Ec + sop*EUEG*W*qwght(n)
          if(ldew) func(n)=func(n)+ sop*EUEG*W
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
            Amat(n,1) = Amat(n,1) + sop*(EUEGPA*W + EUEG*dWdPA)
            Cmat(n,1)=  Cmat(n,1) + sop*(EUEG*dWdGA)
          else
            Amat(n,1) = Amat(n,1) + sop*(EUEGPA*W + EUEG*dWdPA)
            Amat(n,2) = Amat(n,2) + sop*(EUEGPB*W + EUEG*dWdPB) 
            Cmat(n,1) = Cmat(n,1) + sop*EUEG*dWdGA
            Cmat(n,3) = Cmat(n,3) + sop*(EUEG*dWdGB)
          endif
      endIf
c      write (*,*) "PA, PB, GAA, GBB,ipol",PA, PB, GAA, GBB,ipol
c      write (*,*) "FA, FB,FGA, FGB",FA, FB,FGA, FGB
c      Stop
20    continue
      end

      Subroutine xc_cm05_d2()
      implicit none
      call errquit(' cm05: d2 not coded ',0,0)
      return
      end




      Subroutine m05css(DTol,PX,GX,TX,F,FP,FG,FT,EUEG,Chi,EUEGP,
     &                   ChiP,ChiG,ijzy)
      Implicit none
c
#include "errquit.fh"
C
C     Compute the same-spin part of the m05 correlation functional for one grid
C     point and one spin-case.
C
C
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
      if (ijzy.eq.1) then
C     Parameters for M05 Correlation
         sss0=  1.00000d0
         sss1=  3.77344d0
         sss2=  -26.04463d0
         sss3=  30.69913d0
         sss4=  -9.22695d0
      elseif (ijzy.eq.2) then
C     Parameters for M05-2X Correlation
         sss0=  1.00000d0
         sss1=  -3.05430d0
         sss2=  7.61854d0
         sss3=  1.47665d0
         sss4=  -11.92365d0
      else
         call errquit("m05css: illegal value of ijzy",ijzy,UERR)
      endif
      
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
        Call lsdac(dtol,
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



c
c Perdew 91 local correlation functional at one grid point 
c

      Subroutine lsdac(tol_rho,
     R     rs,zeta,PotLC,dLdS,dLdZ,d2LdSS,d2LdSZ,d2LdZZ)

C     Evaluate the Perdew 91 local correlation functional and its
C     derivatives at one point. 
      implicit none
      double precision tol_rho
      double precision rs
      double precision PotLC,dLdS, dLdZ, d2dSS, d2dSZ, d2dZZ
      double precision eps0c(6), eps1c(6), epsc(6)
      double precision F1, F2, F3, F4, F6, F8, F9, F12, F24, F36
      double precision GammaI,Zeta,FZeta,dfZdz,d2fZdz 
      double precision EU,dEUdRS,d2UdRS
      double precision EP,dEPdRS,d2PdRS
      double precision AlphaM,dAMdRS,d2AdRS
      double precision GZ, HZ, dGZ, dHZ, d2GZ, d2HZ
      double precision d2LdSS,d2LdSZ,d2LdZZ,fzzi,
     Z     z2,z3,z4,third,pi
 
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
      
      Call EvFZet(tol_rho,GammaI,Zeta,FZeta,dfZdz,d2fZdz)
      Call EvPWLC(eps0c(1),eps0c(2),eps0c(3),eps0c(4),eps0c(5),
     &    eps0c(6), RS,EU,dEUdRS,d2UdRS)
      Call EvPWLC(eps1c(1),eps1c(2),eps1c(3),eps1c(4),eps1c(5),
     &    eps1c(6),RS,EP,dEPdRS,d2PdRS)
      Call EvPWLC(epsc(1),epsc(2),epsc(3),epsc(4),epsc(5),
     &    epsc(6),RS,AlphaM,dAMdRS,d2AdRS)
      Z2 = Zeta*Zeta
      Z3 = Zeta*Z2
      Z4 = Zeta*Z3
      GZ = FZeta*Z4
      HZ = FZZI*(FZeta-GZ)
      PotLC = EU*(F1-GZ) + EP*GZ - AlphaM*HZ
      dLdS = dEUdRS*(F1-GZ) + dEPdRS*GZ - dAMdRS*HZ
      dGZ = dfZdz*Z4 + F4*FZeta*Z3
      dHZ = FZZI*(dFZdz-dGZ)
      dLdz = (EP-EU)*dGZ - AlphaM*dHZ
      d2GZ = d2fZdz*Z4 + F8*Z3*dfZdz + F12*FZeta*Z2
      d2HZ = FZZI*(d2FZdz-d2GZ)
      d2LdSS = d2UdRS*(F1-GZ) + d2PdRS*GZ - d2AdRS*HZ
      d2LdSZ = (dEPdRS-dEUdRS)*dGZ - dAMdRS*dHZ
      d2LdZZ = (EP-EU)*d2GZ - AlphaM*d2HZ
      Return
      End

c
c   f(zeta)
c
      Subroutine EvFZet(tol_rho,S,Zeta,FZeta,dfZdz,d2fZdz)
      Implicit none
c
c     evaluate f(Zeta) and its derivatives for lsdac.
c
      double precision tol_rho
      double precision Small
      double precision S, Zeta, FZeta,dfZdz,d2fZdz
      double precision Zero, One, Two, Three, Four, Nine, F8, F27
      double precision OMZ, OPZ, OMZ2, OPZ2, OMZ3, OPZ3
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
      OMZ2 = OMZ**2
      OPZ2 = OPZ**2
      F13 = One / Three
      F43 = Four / Three
      F49 = Four / Nine
      F827 = F8 / F27
      If(OMZ.gt.Small) then
        OMZ3 = OMZ ** F13
        fZeta = fZeta + OMZ*OMZ3
        dfZdz = dfZdz - OMZ3
        d2fZdz = d2fZdz + OMZ3/OMZ
        endIf
      If(OPZ.gt.Small) then
        OPZ3 = OPZ ** F13
        fZeta = fZeta + OPZ*OPZ3
        dfZdz = dfZdz + OPZ3
        d2fZdz = d2fZdz + OPZ3/OPZ
        endIf
      fZeta = fZeta * S
      dfZdz = dfZdz * F43 * S
      d2fZdz = d2fZdz * F49 * S
      Return
      End

c
c  pw91 local correlation
c 
      Subroutine EvPWLC(A,A1,B1,B2,B3,B4,RS,V,dVdRS,d2VdRS)
      Implicit none
C
C     Evaluate the interpolation function for PW91 local correlation.
C
      double precision A,A1,B1,B2,B3,B4,RS,V,dVdRS,d2VdRS
      double precision F1,F2, F3, F4
      double precision Q0, RS12, RS32,Q1,Q2
      double precision dQ0dRS,dQ1dRS,dQ2dRS
      double precision d2Q1dS, d2Q2dS
      data F1/1.0d0/, F2/2.0d0/, F3/3.0d0/, F4/4.0d0/
C
      Q0 = -F2*A*(F1+A1*RS)
      RS12 = Sqrt(RS)
      RS32 = RS*RS12
      Q1 = F2*A*(B1*RS12+B2*RS+B3*RS32+B4*RS*RS)
      Q2 = Log(F1+F1/Q1)
      V = Q0*Q2
C
      dQ0dRS = -F2*A*A1
      dQ1dRS = A*(B1/RS12+F2*B2+F3*B3*RS12+F4*B4*RS)
      dQ2dRS = -dQ1dRS/(Q1+Q1**2)
      dVdRS = dQ0dRS*Q2 + Q0*dQ2dRS
C
      d2Q1dS = A*(F3*B3/(RS12*F2)-B1/(RS32*F2)+F4*B4)
      d2Q2dS = (F2*Q1+F1)*(dQ1dRS/(Q1+Q1**2))**2 - d2Q1dS/(Q1+Q1**2)
      d2VdRS = F2*dQ0dRS*dQ2dRS + Q0*d2Q2dS
      Return
      End


c------------------------------------------------------------------------
c     dlDF correlation functional 
c           META GGA
C         utilizes ingredients:
c                              rho   -  density
c                              delrho - gradient of density
c                              tau (tauN)- K.S kinetic energy density
c     References: 
c     [a]	Pernal,Podeszwa,Patkowski,Szalewicz, PRL 103 263201 (2009)

      Subroutine xc_cdldf(tol_rho, cfac, lcfac, nlcfac, rho, delrho, 
     &                     nq, ipol, Ec, qwght, ldew, func,
     &                     tau, Amat, Cmat, Mmat)      

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
c     Threshold parameters
c
      double precision DTol,F1, F2, F3, F4,COpp 
      Data COpp/0.0031d0/,F1/1.0d0/,F2/2.0d0/,
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
      double precision tau(nq,ipol)
 
c
c     Quadrature Weights
c
      double precision qwght(nq)
c
      logical ldew
      double precision func(*)
c
c     Sampling Matrices for the XC Potential
c
      double precision Amat(nq,ipol), Cmat(nq,*)
      double precision Mmat(nq,*)

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
       if (rho(n,1).lt.DTol) goto 20
       if (ipol.eq.1) then
c
c    get the density, gradient, and tau for the alpha spin from the total 
c
         PA = rho(n,1)/F2
         GAA = (    delrho(n,1,1)*delrho(n,1,1) +
     &                 delrho(n,2,1)*delrho(n,2,1) +
     &                 delrho(n,3,1)*delrho(n,3,1))/4
c  In the dldfcss subroutine, we use 2*TA as the tau, so we do not divide 
c  the tau by 2 here

         TA = tau(n,1)
                  
         Call dldfcss(DTol,PA,GAA,TA,FA,FPA,FGA,FTA,EUA,
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

         Ec = Ec + 2.d0*FA*qwght(n)            !factor of 2 account for both spin
         if(ldew) func(n)=func(n)+ FA*2d0
         Amat(n,1)=Amat(n,1)+ FPA
         Cmat(n,1)=  Cmat(n,1) + FGA
         Mmat(n,1)=  Mmat(n,1) + FTA
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

         PA = rho(n,2)
         if (PA.le.DTol) go to 25
         GAA =   delrho(n,1,1)*delrho(n,1,1) +
     &           delrho(n,2,1)*delrho(n,2,1) +
     &          delrho(n,3,1)*delrho(n,3,1)
c
c  In the dldfcss subroutine, we use 2*TA as the tau 
c
         TA = 2*tau(n,1)

         Call dldfcss(DTol,PA,GAA,TA,FA,FPA,FGA,FTA,EUA,
     &                ChiA,EUPA,ChiAP,ChiAG)
         Ec = Ec + FA*qwght(n)     
         if(ldew) func(n)=func(n)+ FA
         Amat(n,1)=Amat(n,1)+ FPA
         Cmat(n,1)=  Cmat(n,1) + FGA
         Mmat(n,1)=  Mmat(n,1) + FTA
c
c  In the dldfcss subroutine, we use 2*TB as the tau, 
c
c
c       Beta 
c
 25      continue
         PB = rho(n,3)
         if (PB.le.DTol) go to 30
         GBB =   delrho(n,1,2)*delrho(n,1,2) +
     &           delrho(n,2,2)*delrho(n,2,2) +
     &          delrho(n,3,2)*delrho(n,3,2)

         TB = 2*tau(n,2) 

         Call dldfcss(DTol,PB,GBB,TB,FB,FPB,FGB,FTB,EUB,
     &                ChiB,EUPB,ChiBP,ChiBG)
         Ec = Ec + FB*qwght(n)          
         if(ldew) func(n)=func(n)+ FB
         Amat(n,2)= Amat(n,2)+ FPB
         Cmat(n,3)=  Cmat(n,3) + FGB
         Mmat(n,2)=  Mmat(n,2) + FTB
      endif

 30   continue
      P = PA + PB
      
      If((PA.gt.DTol).and.(PB.gt.DTol)) then
          RS = (Pi34/P) ** F13
          RSP = -RS/(F3*P)
          Zeta = (PA-PB)/P
          dZdA = (F1-Zeta)/P
          dZdB = (-F1-Zeta)/P
          Call lsdac(dtol,
     D         RS,Zeta,PotLC,dLdS,dLdZ,d2LdSS,d2LdSZ,
     $      d2LdZZ)
          EUEG = P*PotLC - EUA - EUB
          U = COpp*(ChiA+ChiB)/(F1 + COpp*(ChiA+ChiB))
          W = sopp0+U*(sopp1+U*(sopp2+U*(sopp3+U*sopp4)))
          Ec = Ec + sop*EUEG*W*qwght(n)
          if(ldew) func(n)=func(n)+ sop*EUEG*W
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
            Amat(n,1) = Amat(n,1) + sop*(EUEGPA*W + EUEG*dWdPA)
            Cmat(n,1)=  Cmat(n,1) + sop*(EUEG*dWdGA)
          else
            Amat(n,1) = Amat(n,1) + sop*(EUEGPA*W + EUEG*dWdPA)
            Amat(n,2) = Amat(n,2) + sop*(EUEGPB*W + EUEG*dWdPB) 
            Cmat(n,1) = Cmat(n,1) + sop*EUEG*dWdGA
            Cmat(n,3) = Cmat(n,3) + sop*(EUEG*dWdGB)
          endif
      endIf
c      write (*,*) "PA, PB, GAA, GBB,ipol",PA, PB, GAA, GBB,ipol
c      write (*,*) "FA, FB,FGA, FGB",FA, FB,FGA, FGB
c      Stop
20    continue
      end


      Subroutine dldfcss(DTol,PX,GX,TX,F,FP,FG,FT,EUEG,Chi,EUEGP,
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
        Call lsdac(dtol,
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


      Subroutine xc_cdldf_d2()
      implicit none
      call errquit(' cdldf: d2 not coded ',0,0)
      return
      end


