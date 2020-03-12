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
c                              ijzy - 5 revM06-L
c                              ijzy - 6 revM06

      Subroutine xc_cm06(tol_rho, cfac, lcfac, nlcfac, rho, delrho, 
     &                     nq, ipol, Ec, qwght, ldew, func,
     &                     tau, Amat, Cmat, Mmat,ijzy)     
     

c
c$Id$
c
c
c     [a]   Zhao, Y. and  Truhlar, D. G. J. Chem. Phys. 2006, 125, 194101;
c     [b]   Zhao, Y. and  Truhlar, D. G. J. Phys. Chem. A (2006),110(49),13126-13130.
c     [c]   Grafenstein, J., Izotov, D. and Cremer, D. J. Chem. Phys. 2007, 127, 214103.
    

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
      double precision F1, F2, F3, F4,COpp 
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
      sop=1.0d0
      sopp0= 0d0
      sopp1= 0d0
      sopp2= 0d0
      sopp3= 0d0
      sopp4= 0d0
      if (ijzy.eq.1) then
C     Parameters for M06-L Correlation
         sopp0= 6.042374D-01
         sopp1= 1.776783D+02
         sopp2= -2.513252D+02
         sopp3= 7.635173D+01
         sopp4= -1.255699D+01
      elseif (ijzy.eq.2) then
C     Parameters for M06-HF Correlation
         sopp0= 1.674634D+00
         sopp1= 5.732017D+01
         sopp2= 5.955416D+01
         sopp3= -2.311007D+02
         sopp4= 1.255199D+02
      elseif (ijzy.eq.3) then
C     Parameters for M06 Correlation
         sopp0= 3.741539D+00
         sopp1= 2.187098D+02
         sopp2= -4.531252D+02
         sopp3= 2.936479D+02
         sopp4= -6.287470D+01
      elseif (ijzy.eq.4) then
C     Parameters for M06-2X Correlation
         sopp0= 8.833596D-01
         sopp1= 3.357972D+01
         sopp2= -7.043548D+01
         sopp3= 4.978271D+01
         sopp4= -1.852891D+01
      elseif (ijzy.eq.5) then
C     Parameters for revM06-L Correlation
         sopp0= 3.44360696D-01
         sopp1= -5.57080242D-01
         sopp2= -2.009821162D+00
         sopp3= -1.857641887D+00
         sopp4= -1.076639864D+00
      elseif (ijzy.eq.6) then
C     Parameters for revM06 Correlation
         sopp0= 1.222401598D+00
         sopp1= 6.613907336D-01
         sopp2= -1.884581043D+00
         sopp3= -2.780360568D+00
         sopp4= -3.068579344D+00
      else
         call errquit("xc_cm06: illegal value of ijzy",ijzy,UERR)
      endif

      call xc_cvs98(tol_rho, 1.0d0, lfac, nlcfac,
     &           rho, delrho,  nq, ipol,
     &           Ec, qwght, ldew,func,tau,Amat,Cmat,Mmat,ijzy+1)


      
      Pi = F4*ATan(F1)
      F6=6.0d0
      F43 = F4 / F3
      Pi34 = F3 / (F4*Pi)
      F13 = F1 / F3

      do 20 n = 1, nq
       if (rho(n,1).lt.Tol_Rho) goto 20
       if (ipol.eq.1) then
c
c    get the density, gradient, and tau for the alpha spin from the total 
c
         PA = rho(n,1)/F2
         GAA = (    delrho(n,1,1)*delrho(n,1,1) +
     &                 delrho(n,2,1)*delrho(n,2,1) +
     &                 delrho(n,3,1)*delrho(n,3,1))/F4
         if(sqrt(gaa).lt.tol_rho) goto 20
c  In the m05css subroutine, we use 2*TA as the tau, so we do not divide 
c  the tau by 2 here

         TA = tau(n,1)
         if(ta.lt.tol_rho) goto 20
                  
         Call m06css(Tol_Rho,PA,GAA,TA,FA,FPA,FGA,FTA,EUA,
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
         ChiA  = 0.0d0
         ChiAP = 0.0d0
         ChiAG = 0.0d0
         ChiB  = 0.0d0
         ChiBP = 0.0d0
         ChiBG = 0.0d0
         EUA   = 0.0d0
         EUB   = 0.0d0
         EUPA  = 0.0d0
         EUPB  = 0.0d0
c
c       alpha
c

         PA = rho(n,2)
         if (PA.le.Tol_Rho) go to 25
         GAA =   delrho(n,1,1)*delrho(n,1,1) +
     &           delrho(n,2,1)*delrho(n,2,1) +
     &          delrho(n,3,1)*delrho(n,3,1)
c
c  In the m05css subroutine, we use 2*TA as the tau 
c
         TA = 2*tau(n,1) 

         if(sqrt(gaa).lt.tol_rho) goto 25
         if(ta.lt.tol_rho) goto 25

         Call m06css(Tol_Rho,PA,GAA,TA,FA,FPA,FGA,FTA,EUA,
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
         if (PB.le.Tol_Rho) go to 30
         GBB =   delrho(n,1,2)*delrho(n,1,2) +
     &           delrho(n,2,2)*delrho(n,2,2) +
     &          delrho(n,3,2)*delrho(n,3,2)

         TB = 2*tau(n,2)

         if(sqrt(gbb).lt.tol_rho) goto 30
         if(tb.lt.tol_rho) goto 30

         Call m06css(Tol_Rho,PB,GBB,TB,FB,FPB,FGB,FTB,EUB,
     &                ChiB,EUPB,ChiBP,ChiBG,ijzy)
         Ec = Ec + FB*qwght(n)          
         if(ldew) func(n)=func(n)+ FB
         Amat(n,2)= Amat(n,2)+ FPB
         Cmat(n,3)=  Cmat(n,3) + FGB
         Mmat(n,2)=  Mmat(n,2) + FTB
      endif

 30   continue
      P = PA + PB
      
      If((PA.gt.Tol_Rho).and.(PB.gt.Tol_Rho)) then
          RS = (Pi34/P) ** F13
          RSP = -RS/(F3*P)
          Zeta = (PA-PB)/P
          dZdA = (F1-Zeta)/P
          dZdB = (-F1-Zeta)/P
          Call lsdac(tol_rho,
     R         RS,Zeta,PotLC,dLdS,dLdZ,d2LdSS,d2LdSZ,
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

      Subroutine xc_cm06_d2()
      implicit none
      call errquit(' xc06: d2 not coded ',0,0)
      return
      end




      Subroutine m06css(Tol_Rho,PX,GX,TX,F,FP,FG,FT,EUEG,Chi,EUEGP,
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
      sss0=  0d0
      sss1=  0d0
      sss2=  0d0
      sss3=  0d0
      sss4=  0d0
      ss=1.0
      if (ijzy.eq.1) then
C     Parameters for M06-L Correlation
         sss0=  5.349466D-01
         sss1=  5.396620D-01
         sss2=  -3.161217D+01
         sss3=  5.149592D+01
         sss4=  -2.919613D+01
      elseif (ijzy.eq.2) then
C     Parameters for M06-HF Correlation
         sss0=  1.023254D-01
         sss1=  -2.453783D+00
         sss2=  2.913180D+01
         sss3=  -3.494358D+01
         sss4=  2.315955D+01
      elseif (ijzy.eq.3) then
C     Parameters for M06 Correlation
         sss0=  5.094055D-01
         sss1=  -1.491085D+00
         sss2=  1.723922D+01
         sss3=  -3.859018D+01
         sss4=  2.845044D+01
      elseif (ijzy.eq.4) then
C     Parameters for M06-2X Correlation
         sss0=  3.097855D-01
         sss1=  -5.528642D+00
         sss2=  1.347420D+01
         sss3=  -3.213623D+01
         sss4=  2.846742D+01
      elseif (ijzy.eq.5) then
C     Parameters for revM06-L Correlation
         sss0=  1.227659748D+00
         sss1=  8.55201283D-01
         sss2=  -3.113346677D+00
         sss3=  -2.239678026D+00
         sss4=  3.54638962D-01
      elseif (ijzy.eq.6) then
C     Parameters for revM06 Correlation
         sss0=  9.017224575D-01
         sss1=  2.079991827D-01
         sss2=  -1.823747562D+00
         sss3=  -1.384430429D+00
         sss4=  -4.423253381D-01
      else
         call errquit("m06css: illegal value of ijzy",ijzy,UERR)
      endif
      
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
        Call lsdac(tol_rho,
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


