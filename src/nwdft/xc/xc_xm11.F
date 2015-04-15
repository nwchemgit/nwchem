      Subroutine xc_xm11(tol_rho, fac,lfac,nlfac, rho, delrho, 
     &                     Amat, Cmat, nq, ipol, Ex, 
     &                     qwght, ldew, func, tau, Mmat, ijzy)


c   
c$Id$
c
c
c
c**********************************************************************c
c                                                                      c
c  xc_xm11 evaluates the exchange part of the M08 and M11 suite of     c
c  functionals on the grid.                                            c
c  !!! Second derivatives are not available yet.                       c
c                                                                      c
c  Ref: (a) Zhao, Y.  and Truhlar, D. G. JCTC, 2008, 4 , 1849          c
c       (b) Peverati, R. and Truhlar, D. G. J.P.C.Lett. 2011, 2, 2810  c
c       (c) Peverati, R. and Truhlar, D. G. J.P.C.Lett. 2011, 3, 117   c
c                                                                      c
c       ijzy - 1 M08-HX (a)                                            c
c       ijzy - 2 M08-SO (a)                                            c
c       ijzy - 3 M11 (b)                                               c
c       ijzy - 4 M11-L (c)                                             c
c                                                                      c
c Coded by Roberto Peverati (12/11)                                    c
c                                                                      c
c**********************************************************************c    
c
      implicit none
c
#include "errquit.fh"
c      
      double precision fac, Ex
      integer nq, ipol
      logical lfac, nlfac,ldew,   uselc
      double precision func(*)  ! value of the functional [output]
c
c     Charge Density & Its Cube Root
c
      double precision rho(nq,ipol*(ipol+1)/2)
c
c     Charge Density Gradient
c
      double precision delrho(nq,3,ipol)
c
c     Quadrature Weights
c
      double precision qwght(nq)
c
c     Sampling Matrices for the XC Potential & Energy
c
      double precision Amat(nq,ipol), Cmat(nq,*), Mmat(nq,*)
      double precision tol_rho, pi
c
c     kinetic energy density   or  tau
c
      double precision tau(nq,ipol)
      double precision tauN,tauu,DTol
c
c      functional derivatives
c
      double precision dWdT, dTdR, dTdTau
c
c     Intermediate derivative results, etc.
c
      integer n, ijzy
c
      double precision Ax, s, s2
      double precision kapa,kapas,mu,mus
c
      double precision f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11
      double precision F1o3,F2o3,F3o5,F4o3,F5o3,F48,F81 
      double precision Fsig1, Fsig2, Fsig3, Fsig4, Fx1, Fx2
      double precision ElSR, ElLR
      double precision PDUM
      double precision GGA1, GGA2, GGA3, GGA4
      double precision Emu, X, deno
      double precision ds2drho, ds2dg, dfx1ds2
      double precision dfx2ds2, df1dw
      double precision dfx1drho,dfx1dg,dfx2drho,dfx2dg,df2dw
      double precision df3dw, df4dw, delsrdr, dellrdr
      double precision dgga1dr, dgga2dr, dgga3dr, dgga4dr
      double precision df1dr, df1dtau, df2dr, df2dtau
      double precision df3dr, df3dtau, df4dr, df4dtau
      double precision dgga1dg, dgga2dg, dgga3dg, dgga4dg
c
      double precision at00,at01,at02,at03,at04,at05,at06
      double precision at07,at08,at09,at10,at11
      double precision bt00,bt01,bt02,bt03,bt04,bt05,bt06
      double precision bt07,bt08,bt09,bt10,bt11
      double precision ct00,ct01,ct02,ct03,ct04,ct05,ct06
      double precision ct07,ct08,ct09,ct10,ct11
      double precision dt00,dt01,dt02,dt03,dt04,dt05,dt06
      double precision dt07,dt08,dt09,dt10,dt11
      double precision rho43, rho13, rhoo, rho53
      double precision Gamma2, Gamma
      double precision TauUEG, Tsig, Wsig, W1, W2, W3, W4, W5, W6
      double precision W7, W8, W9, W10, W11
c
      parameter( F0=0.0D+00,  F1=1.0D+00,  F2=2.0D+00,
     $           F3=3.0D+00,  F4=4.0D+00,  F5=5.0D+00,
     $           F6=6.0D+00,  F7=7.0D+00,  F8=8.0D+00,
     $           F9=9.0D+00,  F10=10.0D+00,F11=11.0D+00)
c
        pi=acos(-1d0)      
c
        ct00= 0D+00
        ct01= 0D+00
        ct02= 0D+00
        ct03= 0D+00
        ct04= 0D+00
        ct05= 0D+00
        ct06= 0D+00
        ct07= 0D+00
        ct08= 0D+00
        ct09= 0D+00
        ct10= 0D+00
        ct11= 0D+00
C
        dt00= 0D+00
        dt01= 0D+00
        dt02= 0D+00
        dt03= 0D+00
        dt04= 0D+00
        dt05= 0D+00
        dt06= 0D+00
        dt07= 0D+00
        dt08= 0D+00
        dt09= 0D+00
        dt10= 0D+00
        dt11= 0D+00
        at00=0.000000D+00
        at01=0.000000D+00
        at02=0.000000D+00
        at03=0.000000D+00
        at04=0.000000D+00
        at05=0.000000D+00
        at06=0.000000D+00
        at07=0.000000D+00
        at08=0.000000D+00
        at09=0.000000D+00
        at10=0.000000D+00 
        at11=0.000000D+00 
        bt00=0.000000D+00
        bt01=0.000000D+00
        bt02=0.000000D+00
        bt03=0.000000D+00
        bt04=0.000000D+00
        bt05=0.000000D+00
        bt06=0.000000D+00
        bt07=0.000000D+00
        bt08=0.000000D+00
        bt09=0.000000D+00
        bt10=0.000000D+00 
        bt11=0.000000D+00 
        UseLC=.False.
        Emu =0.00D+00
C
      if (ijzy.eq.1) then
C     Parameters for M08-HX
        at00=  1.3340172D+00
        at01= -9.4751087D+00
        at02= -1.2541893D+01
        at03=  9.1369974D+00
        at04=  3.4717204D+01
        at05=  5.8831807D+01
        at06=  7.1369574D+01
        at07=  2.3312961D+01
        at08=  4.8314679D+00
        at09= -6.5044167D+00
        at10= -1.4058265D+01
        at11=  1.2880570D+01

        bt00= -8.5631823D-01
        bt01=  9.2810354D+00
        bt02=  1.2260749D+01
        bt03= -5.5189665D+00
        bt04= -3.5534989D+01
        bt05= -8.2049996D+01
        bt06= -6.8586558D+01
        bt07=  3.6085694D+01
        bt08= -9.3740983D+00
        bt09= -5.9731688D+01
        bt10=  1.6587868D+01
        bt11=  1.3993203D+01
C
        UseLC=.False.
C
       elseif (ijzy.eq.2) then
C     Parameters for M08-SO
        at00= -3.4888428D-01
        at01= -5.8157416D+00
        at02=  3.7550810D+01
        at03=  6.3727406D+01
        at04= -5.3742313D+01
        at05= -9.8595529D+01
        at06=  1.6282216D+01
        at07=  1.7513468D+01
        at08= -6.7627553D+00
        at09=  1.1106658D+01
        at10=  1.5663545D+00
        at11=  8.7603470D+00

        bt00=  7.8098428D-01
        bt01=  5.4538178D+00
        bt02= -3.7853348D+01
        bt03= -6.2295080D+01
        bt04=  4.6713254D+01
        bt05=  8.7321376D+01
        bt06=  1.6053446D+01
        bt07=  2.0126920D+01
        bt08= -4.0343695D+01
        bt09= -5.8577565D+01
        bt10=  2.0890272D+01
        bt11=  1.0946903D+01
C
        UseLC=.False.
C
      elseif (ijzy.eq.3) then
C     Parameters for M11
        at00= -0.18399900D+00
        at01= -1.39046703D+01
        at02=  1.18206837D+01
        at03=  3.10098465D+01
        at04= -5.19625696D+01
        at05=  1.55750312D+01
        at06= -6.94775730D+00
        at07= -1.58465014D+02
        at08= -1.48447565D+00
        at09=  5.51042124D+01
        at10= -1.34714184D+01
        at11=  0.00000000D+00

        bt00=  0.75599900D+00
        bt01=  1.37137944D+01
        bt02= -1.27998304D+01
        bt03= -2.93428814D+01
        bt04=  5.91075674D+01
        bt05= -2.27604866D+01
        bt06= -1.02769340D+01
        bt07=  1.64752731D+02
        bt08=  1.85349258D+01
        bt09= -5.56825639D+01
        bt10=  7.47980859D+00
        bt11=  0.00000000D+00
C
        UseLC=.True.
        Emu =0.25D+00
C
      elseif (ijzy.eq.4) then
C     Parameters for M11-L
        at00=  8.121131D-01
        at01=  1.738124D+01
        at02=  1.154007D+00
        at03=  6.869556D+01
        at04=  1.016864D+02
        at05= -5.887467D+00
        at06=  4.517409D+01
        at07= -2.773149D+00
        at08= -2.617211D+01
        at09=  0.000000D+00
        at10=  0.000000D+00 
        at11=  0.000000D+00
C
        bt00=  1.878869D-01
        bt01= -1.653877D+01
        bt02=  6.755753D-01
        bt03= -7.567572D+01
        bt04= -1.040272D+02
        bt05=  1.831853D+01
        bt06= -5.573352D+01
        bt07= -3.520210D+00
        bt08=  3.724276D+01
        bt09=  0.000000D+00
        bt10=  0.000000D+00
        bt11=  0.000000D+00
C
        ct00= -4.386615D-01
        ct01= -1.214016D+02
        ct02= -1.393573D+02
        ct03= -2.046649D+00
        ct04=  2.804098D+01
        ct05= -1.312258D+01
        ct06= -6.361819D+00
        ct07= -8.055758D-01
        ct08=  3.736551D+00
        ct09=  0.000000D+00
        ct10=  0.000000D+00
        ct11=  0.000000D+00
C
        dt00=  1.438662D+00
        dt01=  1.209465D+02
        dt02=  1.328252D+02
        dt03=  1.296355D+01
        dt04=  5.854866D+00
        dt05= -3.378162D+00
        dt06= -4.423393D+01
        dt07=  6.844475D+00
        dt08=  1.949541D+01
        dt09=  0.000000D+00
        dt10=  0.000000D+00
        dt11=  0.000000D+00
C
        UseLC=.True.
        Emu =0.25D+00
C
      else
        call errquit("xc_xm11: illegal value of ijzy",ijzy,UERR)
      endif
      DTol=tol_rho
      F1o3 = F1/F3 
      F2o3 = F2/F3
      F3o5 = F3/F5
      F4o3 = F4/F3 
      F5o3 = F5/F3
      F48 = 48.0d0
      F81 = 81.0d0
      Ax = -(F3/F2) * (F4o3*Pi)**(-F1o3) 
C     RPBE parameters
      Mus = F10/F81
      kapas = 0.552d0
C     PBE parameters 
      Mu = 0.21951d0
      kapa = 0.804d0
c
      if (ipol.eq.1 )then
c
c        ======> SPIN-RESTRICTED <======
c                     or
c                SPIN-UNPOLARIZED
c
c
         do 10 n = 1, nq
            if (rho(n,1).lt.DTol) goto 10
            rhoo = rho(n,1)/F2
            rho43 = rhoo**F4o3  
            rho13 = rho43/rhoo
            rho53 = rhoo**F5o3
c            
            tauN = tau(n,1)
         if(taun.lt.dtol) goto 10
            tauu=tauN
            TAUUEG=F3O5*((F6*PI*PI)**F2O3)*RHO53
            Tsig =TauUEG/tauu
            Wsig =(Tsig - F1)/(Tsig + F1)
            W1=Wsig 
            W2=Wsig*W1
            W3=Wsig*W2
            W4=Wsig*W3
            W5=Wsig*W4
            W6=Wsig*W5
            W7=Wsig*W6
            W8=Wsig*W7
            W9=Wsig*W8
            W10=Wsig*W9
            W11=Wsig*W10
            Fsig1 =(at00    + at01*W1 + at02*W2 + at03*W3
     $            + at04*W4 + at05*W5 + at06*W6 + at07*W7
     $            + at08*W8 + at09*W9 + at10*W10+ at11*W11)
            Fsig2 =(bt00    + bt01*W1 + bt02*W2 + bt03*W3
     $            + bt04*W4 + bt05*W5 + bt06*W6 + bt07*W7
     $            + bt08*W8 + bt09*W9 + bt10*W10+ bt11*W11)
            Fsig3 =(ct00    + ct01*W1 + ct02*W2 + ct03*W3
     $            + ct04*W4 + ct05*W5 + ct06*W6 + ct07*W7
     $            + ct08*W8 + ct09*W9 + ct10*W10+ ct11*W11)
            Fsig4 =(dt00    + dt01*W1 + dt02*W2 + dt03*W3
     $            + dt04*W4 + dt05*W5 + dt06*W6 + dt07*W7
     $            + dt08*W8 + dt09*W9 + dt10*W10+ dt11*W11)

            Gamma2 =(delrho(n,1,1)*delrho(n,1,1) +
     &              delrho(n,2,1)*delrho(n,2,1) +
     &              delrho(n,3,1)*delrho(n,3,1))/F4
            Gamma = dsqrt(Gamma2)
         if(gamma.lt.dtol) goto 10
         
         X = GAMMA/RHO43
         S = X/(F48*PI*PI)**F1o3
         s2     = s*s
         Deno = (F1 + Mu*s2/kapa)
         fx1=F1+kapa*(F1-F1/Deno)
         fx2=F1+kapas*(F1-DExp(-Mus*s2/kapas))
         If(UseLC) then
           CALL LRCLSDA(EMU,RHOO,ElSR,PDUM)
           ElLR = Ax*Rho43-ElSR
         else
           ElSR = Ax*Rho43
           ElLR = F0
         endIf
         GGA1 = ElSR*fx1
         GGA2 = ElSR*fx2
         GGA3 = ElLR*fx1
         GGA4 = ElLR*fx2
C
          Ex = Ex +F2*(GGA1*Fsig1 + GGA2*Fsig2
     $            +    GGA3*Fsig3 + GGA4*Fsig4)*qwght(n)
          if(ldew) func(n)=func(n)+F2*(GGA1*Fsig1+GGA2*Fsig2
     $                            +    GGA3*Fsig3+GGA4*Fsig4)

c
c     functional derivatives 
c
            ds2dRho = -(F8/F3) * s2/rhoo
            ds2dG = s2/Gamma2
C
            dfx1ds2 = Mu*(F1/(Deno*Deno)) 
            dfx1dRho = dfx1ds2*ds2dRho
            dfx1dG = dfx1ds2*ds2dG
C
            dfx2ds2 = Mus*DExp(-Mus*s2/kapas)
            dfx2dRho = dfx2ds2*ds2dRho
            dfx2dG = dfx2ds2*ds2dG
c
            dF1dW = (at01 + F2*at02*W1 + F3*at03*W2
     $                    + F4*at04*W3 + F5*at05*W4
     $                    + F6*at06*W5 + F7*at07*W6
     $                    + F8*at08*W7 + F9*at09*W8
     $                    + F10*at10*W9+F11*at11*W10)
            dF2dW = (bt01 + F2*bt02*W1 + F3*bt03*W2
     $                    + F4*bt04*W3 + F5*bt05*W4
     $                    + F6*bt06*W5 + F7*bt07*W6
     $                    + F8*bt08*W7 + F9*bt09*W8
     $                    + F10*Bt10*W9+F11*Bt11*W10)
            dF3dW = (ct01 + F2*ct02*W1 + F3*ct03*W2
     $                    + F4*ct04*W3 + F5*ct05*W4
     $                    + F6*ct06*W5 + F7*ct07*W6
     $                    + F8*ct08*W7 + F9*ct09*W8
     $                    + F10*ct10*W9+F11*ct11*W10)
            dF4dW = (dt01 + F2*dt02*W1 + F3*dt03*W2
     $                    + F4*dt04*W3 + F5*dt05*W4
     $                    + F6*dt06*W5 + F7*dt07*W6
     $                    + F8*dt08*W7 + F9*dt09*W8
     $                    + F10*dt10*W9+F11*dt11*W10)
c
            dWdT = F2/((F1 + Tsig)**F2)
            dTdR = ((F6*PI*PI)**F2o3)*(rhoo**F2o3)/tauN
            dTdTau = -TauUEG/tauN**F2
C
           If(UseLC) then
             dElSRdR = PDUM
             dElLRdR = Ax*F4o3*Rho13-PDUM
           else
             dElSRdR=Ax*F4o3*Rho13
             dElLRdR=F0
           endIf  
           dGGA1dR = dElSRdR*fx1 + ElSR*dfx1dRho
           dGGA2dR = dElSRdR*fx2 + ElSR*dfx2dRho 
           dGGA3dR = dElLRdR*fx1 + ElLR*dfx1dRho
           dGGA4dR = dElLRdR*fx2 + ElLR*dfx2dRho 
c
           dF1dR = dF1dW*dWdT*dTdR
           dF1dTau=dF1dW*dWdT*dTdTau
           dF2dR = dF2dW*dWdT*dTdR
           dF2dTau=dF2dW*dWdT*dTdTau
           dF3dR = dF3dW*dWdT*dTdR
           dF3dTau=dF3dW*dWdT*dTdTau
           dF4dR = dF4dW*dWdT*dTdR
           dF4dTau=dF4dW*dWdT*dTdTau
c
           dGGA1dG = ElSR*dfx1dG
           dGGA2dG = ElSR*dfx2dG
           dGGA3dG = ElLR*dfx1dG
           dGGA4dG = ElLR*dfx2dG
c
           Amat(n,1) = Amat(n,1)   +dGGA1dR*Fsig1 + GGA1*dF1dR
     $                             +dGGA2dR*Fsig2 + GGA2*dF2dR
     $                             +dGGA3dR*Fsig3 + GGA3*dF3dR
     $                             +dGGA4dR*Fsig4 + GGA4*dF4dR
           Cmat(n,1)=  Cmat(n,1)  +dGGA1dG*Fsig1 + dGGA2dG*Fsig2
     $                            +dGGA3dG*Fsig3 + dGGA4dG*Fsig4
           Mmat(n,1)=  Mmat(n,1)   +GGA1*dF1dTau + GGA2*dF2dTau
     $                             +GGA3*dF3dTau + GGA4*dF4dTau
c    
10      continue
c
c UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUnrestricted
      else  ! ipol=2
c
c        ======> SPIN-UNRESTRICTED <======

c
c  use spin density functional theory ie n-->2n
c                               Ex=(1/2)Ex[2*alpha] + (1/2)Ex[2*beta]
c
c     Alpha            ALPHA               ALPHA
c
         do 20 n = 1, nq
           if (rho(n,1).lt.DTol) goto 20
           if (rho(n,2).lt.DTol) goto 25           
            rhoo  = rho(n,2)
            rho43 = rhoo**F4o3  
            rho13 = rho43/rhoo
            rho53 = rhoo**F5o3
c       
            tauN = tau(n,1)*F2
     
         if(taun.lt.dtol) goto 25
            tauu=tauN
            TAUUEG=F3O5*((F6*PI*PI)**F2O3)*RHO53
            Tsig =TauUEG/tauu
            Wsig =(Tsig - F1)/(Tsig + F1)
            W1=Wsig 
            W2=Wsig*W1
            W3=Wsig*W2
            W4=Wsig*W3
            W5=Wsig*W4
            W6=Wsig*W5
            W7=Wsig*W6
            W8=Wsig*W7
            W9=Wsig*W8
            W10=Wsig*W9
            W11=Wsig*W10
            Fsig1 =(at00    + at01*W1 + at02*W2 + at03*W3
     $            + at04*W4 + at05*W5 + at06*W6 + at07*W7
     $            + at08*W8 + at09*W9 + at10*W10+ at11*W11)
            Fsig2 =(bt00    + bt01*W1 + bt02*W2 + bt03*W3
     $            + bt04*W4 + bt05*W5 + bt06*W6 + bt07*W7
     $            + bt08*W8 + bt09*W9 + bt10*W10+ bt11*W11)
            Fsig3 =(ct00    + ct01*W1 + ct02*W2 + ct03*W3
     $            + ct04*W4 + ct05*W5 + ct06*W6 + ct07*W7
     $            + ct08*W8 + ct09*W9 + ct10*W10+ ct11*W11)
            Fsig4 =(dt00    + dt01*W1 + dt02*W2 + dt03*W3
     $            + dt04*W4 + dt05*W5 + dt06*W6 + dt07*W7
     $            + dt08*W8 + dt09*W9 + dt10*W10+ dt11*W11)

            Gamma2 =(delrho(n,1,1)*delrho(n,1,1) +
     &              delrho(n,2,1)*delrho(n,2,1) +
     &              delrho(n,3,1)*delrho(n,3,1))
            Gamma = dsqrt(Gamma2)
         if(gamma.lt.dtol) goto 25
         
         X = GAMMA/RHO43
         S = X/(F48*PI*PI)**F1o3
         s2     = s*s
         Deno = (F1 + Mu*s2/kapa)
         fx1=F1+kapa*(F1-F1/Deno)
         fx2=F1+kapas*(F1-DExp(-Mus*s2/kapas))
         If(UseLC) then
           CALL LRCLSDA(EMU,RHOO,ElSR,PDUM)
           ElLR = Ax*Rho43-ElSR
         else
           ElSR = Ax*Rho43
           ElLR = F0
         endIf
         GGA1 = ElSR*fx1
         GGA2 = ElSR*fx2
         GGA3 = ElLR*fx1
         GGA4 = ElLR*fx2
C
          Ex = Ex +   (GGA1*Fsig1 + GGA2*Fsig2
     $            +    GGA3*Fsig3 + GGA4*Fsig4)*qwght(n)
          if(ldew) func(n)=func(n)+   (GGA1*Fsig1+GGA2*Fsig2
     $                            +    GGA3*Fsig3+GGA4*Fsig4)
c
c     functional derivatives 
c
            ds2dRho = -(F8/F3) * s2/rhoo
            ds2dG = s2/Gamma2
C
            dfx1ds2 = Mu*(F1/(Deno*Deno)) 
            dfx1dRho = dfx1ds2*ds2dRho
            dfx1dG = dfx1ds2*ds2dG
C
            dfx2ds2 = Mus*DExp(-Mus*s2/kapas)
            dfx2dRho = dfx2ds2*ds2dRho
            dfx2dG = dfx2ds2*ds2dG
c
            dF1dW = (at01 + F2*at02*W1 + F3*at03*W2
     $                    + F4*at04*W3 + F5*at05*W4
     $                    + F6*at06*W5 + F7*at07*W6
     $                    + F8*at08*W7 + F9*at09*W8
     $                    + F10*at10*W9+F11*at11*W10)
            dF2dW = (bt01 + F2*bt02*W1 + F3*bt03*W2
     $                    + F4*bt04*W3 + F5*bt05*W4
     $                    + F6*bt06*W5 + F7*bt07*W6
     $                    + F8*bt08*W7 + F9*bt09*W8
     $                    + F10*Bt10*W9+F11*Bt11*W10)
            dF3dW = (ct01 + F2*ct02*W1 + F3*ct03*W2
     $                    + F4*ct04*W3 + F5*ct05*W4
     $                    + F6*ct06*W5 + F7*ct07*W6
     $                    + F8*ct08*W7 + F9*ct09*W8
     $                    + F10*ct10*W9+F11*ct11*W10)
            dF4dW = (dt01 + F2*dt02*W1 + F3*dt03*W2
     $                    + F4*dt04*W3 + F5*dt05*W4
     $                    + F6*dt06*W5 + F7*dt07*W6
     $                    + F8*dt08*W7 + F9*dt09*W8
     $                    + F10*dt10*W9+F11*dt11*W10)

            dWdT = F2/((F1 + Tsig)**F2)
            dTdR = ((F6*PI*PI)**F2o3)*(rhoo**F2o3)/tauN
            dTdTau = -TauUEG/tauN**F2
C
           If(UseLC) then
             dElSRdR = PDUM
             dElLRdR = Ax*F4o3*Rho13-PDUM
           else
             dElSRdR=Ax*F4o3*Rho13
             dElLRdR=F0
           endIf  
           dGGA1dR = dElSRdR*fx1 + ElSR*dfx1dRho
           dGGA2dR = dElSRdR*fx2 + ElSR*dfx2dRho 
           dGGA3dR = dElLRdR*fx1 + ElLR*dfx1dRho
           dGGA4dR = dElLRdR*fx2 + ElLR*dfx2dRho 
c
           dF1dR = dF1dW*dWdT*dTdR
           dF1dTau=dF1dW*dWdT*dTdTau
           dF2dR = dF2dW*dWdT*dTdR
           dF2dTau=dF2dW*dWdT*dTdTau
           dF3dR = dF3dW*dWdT*dTdR
           dF3dTau=dF3dW*dWdT*dTdTau
           dF4dR = dF4dW*dWdT*dTdR
           dF4dTau=dF4dW*dWdT*dTdTau
c
           dGGA1dG = ElSR*dfx1dG
           dGGA2dG = ElSR*dfx2dG
           dGGA3dG = ElLR*dfx1dG
           dGGA4dG = ElLR*dfx2dG
c
           Amat(n,1) = Amat(n,1)   +dGGA1dR*Fsig1 + GGA1*dF1dR
     $                             +dGGA2dR*Fsig2 + GGA2*dF2dR
     $                             +dGGA3dR*Fsig3 + GGA3*dF3dR
     $                             +dGGA4dR*Fsig4 + GGA4*dF4dR
           Cmat(n,1)=  Cmat(n,1)   +dGGA1dG*Fsig1 + dGGA2dG*Fsig2
     $                             +dGGA3dG*Fsig3 + dGGA4dG*Fsig4
           Mmat(n,1)=  Mmat(n,1)   +GGA1*dF1dTau  + GGA2*dF2dTau
     $                             +GGA3*dF3dTau  + GGA4*dF4dTau
c
25         continue
c
c     Beta               BETA           BETA
c
            if (rho(n,3).lt.DTol) goto 20
            rhoo  = rho(n,3)
            rho43 = rhoo**F4o3  
            rho13 = rho43/rhoo
            rho53 = rhoo**F5o3
c
             
            tauN = tau(n,2)*F2
     
         if(taun.lt.dtol) goto 20
            tauu=tauN
            TAUUEG=F3O5*((F6*PI*PI)**F2O3)*RHO53
            Tsig =TauUEG/tauu
            Wsig =(Tsig - F1)/(Tsig + F1)
            W1=Wsig 
            W2=Wsig*W1
            W3=Wsig*W2
            W4=Wsig*W3
            W5=Wsig*W4
            W6=Wsig*W5
            W7=Wsig*W6
            W8=Wsig*W7
            W9=Wsig*W8
            W10=Wsig*W9
            W11=Wsig*W10
            Fsig1 =(at00    + at01*W1 + at02*W2 + at03*W3
     $            + at04*W4 + at05*W5 + at06*W6 + at07*W7
     $            + at08*W8 + at09*W9 + at10*W10+ at11*W11)
            Fsig2 =(bt00    + bt01*W1 + bt02*W2 + bt03*W3
     $            + bt04*W4 + bt05*W5 + bt06*W6 + bt07*W7
     $            + bt08*W8 + bt09*W9 + bt10*W10+ bt11*W11)
            Fsig3 =(ct00    + ct01*W1 + ct02*W2 + ct03*W3
     $            + ct04*W4 + ct05*W5 + ct06*W6 + ct07*W7
     $            + ct08*W8 + ct09*W9 + ct10*W10+ ct11*W11)
            Fsig4 =(dt00    + dt01*W1 + dt02*W2 + dt03*W3
     $            + dt04*W4 + dt05*W5 + dt06*W6 + dt07*W7
     $            + dt08*W8 + dt09*W9 + dt10*W10+ dt11*W11)

            Gamma2 =(delrho(n,1,2)*delrho(n,1,2) +
     &              delrho(n,2,2)*delrho(n,2,2) +
     &              delrho(n,3,2)*delrho(n,3,2))
            Gamma = dsqrt(Gamma2)
         if(gamma.lt.dtol) goto 20
         
         X = GAMMA/RHO43
         S = X/(F48*PI*PI)**F1o3
         s2     = s*s
         Deno = (F1 + Mu*s2/kapa)
         fx1=F1+kapa*(F1-F1/Deno)
         fx2=F1+kapas*(F1-DExp(-Mus*s2/kapas))
         If(UseLC) then
           CALL LRCLSDA(EMU,RHOO,ElSR,PDUM)
           ElLR = Ax*Rho43-ElSR
         else
           ElSR = Ax*Rho43
           ElLR = F0
         endIf
         GGA1 = ElSR*fx1
         GGA2 = ElSR*fx2
         GGA3 = ElLR*fx1
         GGA4 = ElLR*fx2
C
          Ex = Ex +   (GGA1*Fsig1 + GGA2*Fsig2
     $            +    GGA3*Fsig3 + GGA4*Fsig4)*qwght(n)
          if(ldew) func(n)=func(n)+   (GGA1*Fsig1+GGA2*Fsig2
     $                            +    GGA3*Fsig3+GGA4*Fsig4)

c
c     functional derivatives 
c
            ds2dRho = -(F8/F3) * s2/rhoo
            ds2dG = s2/Gamma2
C
            dfx1ds2 = Mu*(F1/(Deno*Deno)) 
            dfx1dRho = dfx1ds2*ds2dRho
            dfx1dG = dfx1ds2*ds2dG
C
            dfx2ds2 = Mus*DExp(-Mus*s2/kapas)
            dfx2dRho = dfx2ds2*ds2dRho
            dfx2dG = dfx2ds2*ds2dG
c
            dF1dW = (at01 + F2*at02*W1 + F3*at03*W2
     $                    + F4*at04*W3 + F5*at05*W4
     $                    + F6*at06*W5 + F7*at07*W6
     $                    + F8*at08*W7 + F9*at09*W8
     $                    + F10*at10*W9+F11*at11*W10)
            dF2dW = (bt01 + F2*bt02*W1 + F3*bt03*W2
     $                    + F4*bt04*W3 + F5*bt05*W4
     $                    + F6*bt06*W5 + F7*bt07*W6
     $                    + F8*bt08*W7 + F9*bt09*W8
     $                    + F10*Bt10*W9+F11*Bt11*W10)
            dF3dW = (ct01 + F2*ct02*W1 + F3*ct03*W2
     $                    + F4*ct04*W3 + F5*ct05*W4
     $                    + F6*ct06*W5 + F7*ct07*W6
     $                    + F8*ct08*W7 + F9*ct09*W8
     $                    + F10*ct10*W9+F11*ct11*W10)
            dF4dW = (dt01 + F2*dt02*W1 + F3*dt03*W2
     $                    + F4*dt04*W3 + F5*dt05*W4
     $                    + F6*dt06*W5 + F7*dt07*W6
     $                    + F8*dt08*W7 + F9*dt09*W8
     $                    + F10*dt10*W9+F11*dt11*W10)

            dWdT = F2/((F1 + Tsig)**F2)
            dTdR = ((F6*PI*PI)**F2o3)*(rhoo**F2o3)/tauN
            dTdTau = -TauUEG/tauN**F2
C
           If(UseLC) then
             dElSRdR = PDUM
             dElLRdR = Ax*F4o3*Rho13-PDUM
           else
             dElSRdR=Ax*F4o3*Rho13
             dElLRdR=F0
           endIf  
           dGGA1dR = dElSRdR*fx1 + ElSR*dfx1dRho
           dGGA2dR = dElSRdR*fx2 + ElSR*dfx2dRho 
           dGGA3dR = dElLRdR*fx1 + ElLR*dfx1dRho
           dGGA4dR = dElLRdR*fx2 + ElLR*dfx2dRho 
c
           dF1dR = dF1dW*dWdT*dTdR
           dF1dTau=dF1dW*dWdT*dTdTau
           dF2dR = dF2dW*dWdT*dTdR
           dF2dTau=dF2dW*dWdT*dTdTau
           dF3dR = dF3dW*dWdT*dTdR
           dF3dTau=dF3dW*dWdT*dTdTau
           dF4dR = dF4dW*dWdT*dTdR
           dF4dTau=dF4dW*dWdT*dTdTau
c
           dGGA1dG = ElSR*dfx1dG
           dGGA2dG = ElSR*dfx2dG
           dGGA3dG = ElLR*dfx1dG
           dGGA4dG = ElLR*dfx2dG
c
           Amat(n,2) = Amat(n,2)   +dGGA1dR*Fsig1 + GGA1*dF1dR
     $                             +dGGA2dR*Fsig2 + GGA2*dF2dR
     $                             +dGGA3dR*Fsig3 + GGA3*dF3dR
     $                             +dGGA4dR*Fsig4 + GGA4*dF4dR
           Cmat(n,3)=  Cmat(n,3)   +dGGA1dG*Fsig1 + dGGA2dG*Fsig2
     $                             +dGGA3dG*Fsig3 + dGGA4dG*Fsig4
           Mmat(n,2)=  Mmat(n,2)   +GGA1*dF1dTau  + GGA2*dF2dTau
     $                             +GGA3*dF3dTau  + GGA4*dF4dTau
c
20      continue
      endif
      return
      end
c
      Subroutine xc_xm11_d2()
      call errquit(' not coded ',0,0)
      return
      end
c
      SUBROUTINE LRCLSDA(Emu,Rho,F,D1F)
c
c***********************************************
c                                               
c   INPUT:                                      
c      Emu - Value of mu (or omega)
c      Rho - Spin density                 
c                                               
c   OUTPUT:                                     
c      F      - Functional value               
c      D1F    - First derivative               
c                                               
c***********************************************
c
      IMPLICIT REAL*8 (a-h,o-z)
      Save F1, F2, F3, F4, F5, F6, F7, F8, F9
      DATA F1/1.0D+00/,F2/2.0D+00/,F3/3.0D+00/,F4/4.0D+00/,F5/5.0D+00/,
     $     F6/6.0D+00/,F7/7.0D+00/,F8/8.0D+00/,F9/9.0D+00/
C
      PARAMETER( PI = 3.1415926535897932384626433832795D+00 )
C
      F1o2 = F1 / F2
      F1o3 = F1 / F3
      F1o4 = F1 / F4
      F4o3 = F4 / F3
      F8o3 = F8 / F3
      PI12 = SQRT(Pi)
C
      AX   = -(F3/F2) * (F4o3*PI)**(-F1o3)
      Cmu  = (F6*Pi**F2)**F1o3   
C
      Rho13 = Rho**F1o3
      Rho43 = Rho**F4o3
c
      tmu  = Emu/(F2*Cmu*Rho13)
      tmu2 = tmu*tmu
      tmu3 = tmu*tmu2
c
      W    = DExp(-F1o4/tmu2)
      ERFV = DErf( F1o2/tmu)
      dtmudR = -F1o3*tmu / Rho
c
      Fsr = F1-F4o3*tmu*(-F6*tmu+F8*tmu3+W*
     $        (F4*tmu-F8*tmu3)+F2*PI12*ERFV)
      dFsrdtmu = F8o3*(F2*tmu*(F3-F8*tmu2+W*
     $          (-F1+F8*tmu2))-PI12*ERFV)
c
      F = Ax*Rho43*Fsr
      D1F = Ax*F4o3*Rho13*Fsr + Ax*Rho43*(dFsrdtmu*dtmudR)
c
      RETURN
      END

