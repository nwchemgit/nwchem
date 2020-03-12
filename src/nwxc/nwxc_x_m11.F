#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if !defined(NWAD_PRINT)
C> \ingroup nwxc
C> @{
C>
C> \file nwxc_x_m11.F
C> Implementation of the M08 and M11 exchange functionals
C>
C> @}
#endif
#endif
C>
C> \ingroup nwxc_priv
C> @{
C>
C> \brief The M08 exchange functionals
C>
C> This routine implements exchange meta-GGA functionals 
C> of the M08 [1] family of functionals.
C>
C> ### References ###
C>
C> [1] Y. Zhao, D. G. Truhlar,
C>     "Exploring the limit of accuracy of the global hybrid meta 
C>     density functional for main-group thermochemistry, kinetics,
C>     and noncovalent interactions",
C>     J. Chem. Theory Comput. <b>4</b>, 1849-1868 (2008), DOI:
C>     <a href="https://doi.org/10.1021/ct800246v">
C>     10.1021/ct800246v</a>.
C>
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
      Subroutine nwxc_x_m08_p(param, tol_rho, ipol, nq, wght, rho,
     &                        rgamma, tau, func)
#else
      Subroutine nwxc_x_m08(param, tol_rho, ipol, nq, wght, rho, rgamma,
     &                      tau, func)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      Subroutine nwxc_x_m08_d2(param, tol_rho, ipol, nq, wght, rho,
     &                         rgamma, tau, func)
#else
      Subroutine nwxc_x_m08_d3(param, tol_rho, ipol, nq, wght, rho,
     &                         rgamma, tau, func)
#endif


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
#include "nwad.fh"
c
      implicit none
c
#include "intf_nwxc_lrclsda.fh"
c
#include "nwxc_param.fh"
c
#if defined(NWAD_PRINT)
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      type(nwad_dble)::param(*) !< [Input] Parameters of the functional
      type(nwad_dble)::at00,at01,at02,at03,at04,at05,at06
      type(nwad_dble)::at07,at08,at09,at10,at11
      type(nwad_dble)::bt00,bt01,bt02,bt03,bt04,bt05,bt06
      type(nwad_dble)::bt07,bt08,bt09,bt10,bt11
      type(nwad_dble)::ct00,ct01,ct02,ct03,ct04,ct05,ct06
      type(nwad_dble)::ct07,ct08,ct09,ct10,ct11
      type(nwad_dble)::dt00,dt01,dt02,dt03,dt04,dt05,dt06
      type(nwad_dble)::dt07,dt08,dt09,dt10,dt11
#else
      double precision param(*) !< [Input] Parameters of the functional
      double precision at00,at01,at02,at03,at04,at05,at06
      double precision at07,at08,at09,at10,at11
      double precision bt00,bt01,bt02,bt03,bt04,bt05,bt06
      double precision bt07,bt08,bt09,bt10,bt11
      double precision ct00,ct01,ct02,ct03,ct04,ct05,ct06
      double precision ct07,ct08,ct09,ct10,ct11
      double precision dt00,dt01,dt02,dt03,dt04,dt05,dt06
      double precision dt07,dt08,dt09,dt10,dt11
#endif
#else
      double precision param(*) !< [Input] Parameters of the functional
                                !< (see Table 2 of [1])
                                !< - param(1): \f$ a_0 \f$
                                !< - param(2): \f$ a_1 \f$
                                !< - param(3): \f$ a_2 \f$
                                !< - param(4): \f$ a_3 \f$
                                !< - param(5): \f$ a_4 \f$
                                !< - param(6): \f$ a_5 \f$
                                !< - param(7): \f$ a_6 \f$
                                !< - param(8): \f$ a_7 \f$
                                !< - param(9): \f$ a_8 \f$
                                !< - param(10): \f$ a_9 \f$
                                !< - param(11): \f$ a_{10} \f$
                                !< - param(12): \f$ a_{11} \f$
                                !< - param(13): \f$ b_0 \f$
                                !< - param(14): \f$ b_1 \f$
                                !< - param(15): \f$ b_2 \f$
                                !< - param(16): \f$ b_3 \f$
                                !< - param(17): \f$ b_4 \f$
                                !< - param(18): \f$ b_5 \f$
                                !< - param(19): \f$ b_6 \f$
                                !< - param(20): \f$ b_7 \f$
                                !< - param(21): \f$ b_8 \f$
                                !< - param(22): \f$ b_9 \f$
                                !< - param(23): \f$ b_{10} \f$
                                !< - param(24): \f$ b_{11} \f$
                                !< - param(25): \f$ c_0 \f$
                                !< - param(26): \f$ c_1 \f$
                                !< - param(27): \f$ c_2 \f$
                                !< - param(28): \f$ c_3 \f$
                                !< - param(29): \f$ c_4 \f$
                                !< - param(30): \f$ c_5 \f$
                                !< - param(31): \f$ c_6 \f$
                                !< - param(32): \f$ c_7 \f$
                                !< - param(33): \f$ c_8 \f$
                                !< - param(34): \f$ c_9 \f$
                                !< - param(35): \f$ c_{10} \f$
                                !< - param(36): \f$ c_{11} \f$
                                !< - param(37): \f$ d_0 \f$
                                !< - param(38): \f$ d_1 \f$
                                !< - param(39): \f$ d_2 \f$
                                !< - param(40): \f$ d_3 \f$
                                !< - param(41): \f$ d_4 \f$
                                !< - param(42): \f$ d_5 \f$
                                !< - param(43): \f$ d_6 \f$
                                !< - param(44): \f$ d_7 \f$
                                !< - param(45): \f$ d_8 \f$
                                !< - param(46): \f$ d_9 \f$
                                !< - param(47): \f$ d_{10} \f$
                                !< - param(48): \f$ d_{11} \f$
      double precision at00,at01,at02,at03,at04,at05,at06
      double precision at07,at08,at09,at10,at11
      double precision bt00,bt01,bt02,bt03,bt04,bt05,bt06
      double precision bt07,bt08,bt09,bt10,bt11
      double precision ct00,ct01,ct02,ct03,ct04,ct05,ct06
      double precision ct07,ct08,ct09,ct10,ct11
      double precision dt00,dt01,dt02,dt03,dt04,dt05,dt06
      double precision dt07,dt08,dt09,dt10,dt11
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
      type(nwad_dble)::tauN,tauu
      double precision pi,DTol
c
c      functional derivatives
c
      double precision dWdT, dTdR, dTdTau
c
c     Intermediate derivative results, etc.
c
      integer n
c
      double precision Ax
      type(nwad_dble)::s, s2
      double precision kapa,kapas,mu,mus
c
      double precision f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11
      double precision F1o3,F2o3,F3o5,F4o3,F5o3,F48,F81 
      type(nwad_dble)::Fsig1, Fsig2, Fsig3, Fsig4, Fx1, Fx2
      type(nwad_dble)::ElSR, ElLR
c     double precision PDUM
      type(nwad_dble)::GGA1, GGA2, GGA3, GGA4
      type(nwad_dble):: X, deno
      double precision Emu
      double precision ds2drho, ds2dg, dfx1ds2
      double precision dfx2ds2, df1dw
      double precision dfx1drho,dfx1dg,dfx2drho,dfx2dg,df2dw
      double precision df3dw, df4dw, delsrdr, dellrdr
      double precision dgga1dr, dgga2dr, dgga3dr, dgga4dr
      double precision df1dr, df1dtau, df2dr, df2dtau
      double precision df3dr, df3dtau, df4dr, df4dtau
      double precision dgga1dg, dgga2dg, dgga3dg, dgga4dg
c
      type(nwad_dble)::rho43, rhoo, rho53
c     type(nwad_dble)::Gamma2, Gamma
      type(nwad_dble)::Gamma2
      type(nwad_dble)::TauUEG, Tsig, Wsig, W1, W2, W3, W4, W5, W6
      type(nwad_dble)::W7, W8, W9, W10, W11
      logical UseLC
c
      parameter( F0=0.0D+00,  F1=1.0D+00,  F2=2.0D+00,
     $           F3=3.0D+00,  F4=4.0D+00,  F5=5.0D+00,
     $           F6=6.0D+00,  F7=7.0D+00,  F8=8.0D+00,
     $           F9=9.0D+00,  F10=10.0D+00,F11=11.0D+00)
c
        pi=acos(-1d0)      
c
        at00=param( 1)
        at01=param( 2)
        at02=param( 3)
        at03=param( 4)
        at04=param( 5)
        at05=param( 6)
        at06=param( 7)
        at07=param( 8)
        at08=param( 9)
        at09=param(10)
        at10=param(11) 
        at11=param(12) 
c
        bt00=param(13)
        bt01=param(14)
        bt02=param(15)
        bt03=param(16)
        bt04=param(17)
        bt05=param(18)
        bt06=param(19)
        bt07=param(20)
        bt08=param(21)
        bt09=param(22)
        bt10=param(23) 
        bt11=param(24) 
c
        ct00=param(25)
        ct01=param(26)
        ct02=param(27)
        ct03=param(28)
        ct04=param(29)
        ct05=param(30)
        ct06=param(31)
        ct07=param(32)
        ct08=param(33)
        ct09=param(34)
        ct10=param(35)
        ct11=param(36)
C
        dt00=param(37)
        dt01=param(38)
        dt02=param(39)
        dt03=param(40)
        dt04=param(41)
        dt05=param(42)
        dt06=param(43)
        dt07=param(44)
        dt08=param(45)
        dt09=param(46)
        dt10=param(47)
        dt11=param(48)
c
        UseLC=.False.
        Emu  = 0.0d0
C
c     if (ijzy.eq.1) then
C     Parameters for M08-HX
c       at00=  1.3340172D+00
c       at01= -9.4751087D+00
c       at02= -1.2541893D+01
c       at03=  9.1369974D+00
c       at04=  3.4717204D+01
c       at05=  5.8831807D+01
c       at06=  7.1369574D+01
c       at07=  2.3312961D+01
c       at08=  4.8314679D+00
c       at09= -6.5044167D+00
c       at10= -1.4058265D+01
c       at11=  1.2880570D+01

c       bt00= -8.5631823D-01
c       bt01=  9.2810354D+00
c       bt02=  1.2260749D+01
c       bt03= -5.5189665D+00
c       bt04= -3.5534989D+01
c       bt05= -8.2049996D+01
c       bt06= -6.8586558D+01
c       bt07=  3.6085694D+01
c       bt08= -9.3740983D+00
c       bt09= -5.9731688D+01
c       bt10=  1.6587868D+01
c       bt11=  1.3993203D+01
c
c       UseLC=.False.
c
c      elseif (ijzy.eq.2) then
c     Parameters for M08-SO
c       at00= -3.4888428D-01
c       at01= -5.8157416D+00
c       at02=  3.7550810D+01
c       at03=  6.3727406D+01
c       at04= -5.3742313D+01
c       at05= -9.8595529D+01
c       at06=  1.6282216D+01
c       at07=  1.7513468D+01
c       at08= -6.7627553D+00
c       at09=  1.1106658D+01
c       at10=  1.5663545D+00
c       at11=  8.7603470D+00

c       bt00=  7.8098428D-01
c       bt01=  5.4538178D+00
c       bt02= -3.7853348D+01
c       bt03= -6.2295080D+01
c       bt04=  4.6713254D+01
c       bt05=  8.7321376D+01
c       bt06=  1.6053446D+01
c       bt07=  2.0126920D+01
c       bt08= -4.0343695D+01
c       bt09= -5.8577565D+01
c       bt10=  2.0890272D+01
c       bt11=  1.0946903D+01
c
c       UseLC=.False.
c
c     elseif (ijzy.eq.3) then
c     Parameters for M11
c       at00= -0.18399900D+00
c       at01= -1.39046703D+01
c       at02=  1.18206837D+01
c       at03=  3.10098465D+01
c       at04= -5.19625696D+01
c       at05=  1.55750312D+01
c       at06= -6.94775730D+00
c       at07= -1.58465014D+02
c       at08= -1.48447565D+00
c       at09=  5.51042124D+01
c       at10= -1.34714184D+01
c       at11=  0.00000000D+00

c       bt00=  0.75599900D+00
c       bt01=  1.37137944D+01
c       bt02= -1.27998304D+01
c       bt03= -2.93428814D+01
c       bt04=  5.91075674D+01
c       bt05= -2.27604866D+01
c       bt06= -1.02769340D+01
c       bt07=  1.64752731D+02
c       bt08=  1.85349258D+01
c       bt09= -5.56825639D+01
c       bt10=  7.47980859D+00
c       bt11=  0.00000000D+00
c
c       UseLC=.True.
c       Emu =0.25D+00
c
c     elseif (ijzy.eq.4) then
c     Parameters for M11-L
c       at00=  8.121131D-01
c       at01=  1.738124D+01
c       at02=  1.154007D+00
c       at03=  6.869556D+01
c       at04=  1.016864D+02
c       at05= -5.887467D+00
c       at06=  4.517409D+01
c       at07= -2.773149D+00
c       at08= -2.617211D+01
c       at09=  0.000000D+00
c       at10=  0.000000D+00 
c       at11=  0.000000D+00
c
c       bt00=  1.878869D-01
c       bt01= -1.653877D+01
c       bt02=  6.755753D-01
c       bt03= -7.567572D+01
c       bt04= -1.040272D+02
c       bt05=  1.831853D+01
c       bt06= -5.573352D+01
c       bt07= -3.520210D+00
c       bt08=  3.724276D+01
c       bt09=  0.000000D+00
c       bt10=  0.000000D+00
c       bt11=  0.000000D+00
c
c       ct00= -4.386615D-01
c       ct01= -1.214016D+02
c       ct02= -1.393573D+02
c       ct03= -2.046649D+00
c       ct04=  2.804098D+01
c       ct05= -1.312258D+01
c       ct06= -6.361819D+00
c       ct07= -8.055758D-01
c       ct08=  3.736551D+00
c       ct09=  0.000000D+00
c       ct10=  0.000000D+00
c       ct11=  0.000000D+00
c
c       dt00=  1.438662D+00
c       dt01=  1.209465D+02
c       dt02=  1.328252D+02
c       dt03=  1.296355D+01
c       dt04=  5.854866D+00
c       dt05= -3.378162D+00
c       dt06= -4.423393D+01
c       dt07=  6.844475D+00
c       dt08=  1.949541D+01
c       dt09=  0.000000D+00
c       dt10=  0.000000D+00
c       dt11=  0.000000D+00
c
c       UseLC=.True.
c       Emu =0.25D+00
c
c     endif
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
            if (rho(n,R_T).lt.DTol) goto 10
            rhoo = rho(n,R_T)/F2
            rho43 = rhoo**F4o3  
c           rho13 = rho43/rhoo
            rho53 = rhoo**F5o3
c            
            tauN = tau(n,T_T)
            tauu=tauN
            TAUUEG=F3O5*((F6*PI*PI)**F2O3)*RHO53
c           Tsig =TauUEG/tauu
c           Wsig =(Tsig - F1)/(Tsig + F1)
            Wsig = (TauUEG-tauu)/(TauUEG+tauu)
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

c           Gamma2 =(delrho(n,1,1)*delrho(n,1,1) +
c    &              delrho(n,2,1)*delrho(n,2,1) +
c    &              delrho(n,3,1)*delrho(n,3,1))/F4
            Gamma2 = rgamma(n,G_TT)/F4
c           Gamma = sqrt(Gamma2)
c        
c        X = GAMMA/RHO43
c        S = X/(F48*PI*PI)**F1o3
c        s2     = s*s
         s2 = Gamma2/(rho43*rho43)/(F48*PI*PI)**F2o3
         Deno = (F1 + Mu*s2/kapa)
         fx1=F1+kapa*(F1-F1/Deno)
         fx2=F1+kapas*(F1-Exp(-Mus*s2/kapas))
         If(UseLC) then
#if defined(THIRD_DERIV)
           CALL nwxc_LRCLSDA_d3(EMU,RHOO,ElSR)
#elif defined(SECOND_DERIV)
           CALL nwxc_LRCLSDA_d2(EMU,RHOO,ElSR)
#else
#if defined(NWAD_PRINT)
           CALL nwxc_LRCLSDA_p(EMU,RHOO,ElSR)
#else
           CALL nwxc_LRCLSDA(EMU,RHOO,ElSR)
#endif
#endif
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
c        Ex = Ex +F2*(GGA1*Fsig1 + GGA2*Fsig2
c    $           +    GGA3*Fsig3 + GGA4*Fsig4)*qwght(n)
         func(n)=func(n)+F2*(GGA1*Fsig1+GGA2*Fsig2
     $                  +    GGA3*Fsig3+GGA4*Fsig4)*wght

c
c     functional derivatives 
c
c           ds2dRho = -(F8/F3) * s2/rhoo
c           ds2dG = s2/Gamma2
C
c           dfx1ds2 = Mu*(F1/(Deno*Deno)) 
c           dfx1dRho = dfx1ds2*ds2dRho
c           dfx1dG = dfx1ds2*ds2dG
C
c           dfx2ds2 = Mus*DExp(-Mus*s2/kapas)
c           dfx2dRho = dfx2ds2*ds2dRho
c           dfx2dG = dfx2ds2*ds2dG
c
c           dF1dW = (at01 + F2*at02*W1 + F3*at03*W2
c    $                    + F4*at04*W3 + F5*at05*W4
c    $                    + F6*at06*W5 + F7*at07*W6
c    $                    + F8*at08*W7 + F9*at09*W8
c    $                    + F10*at10*W9+F11*at11*W10)
c           dF2dW = (bt01 + F2*bt02*W1 + F3*bt03*W2
c    $                    + F4*bt04*W3 + F5*bt05*W4
c    $                    + F6*bt06*W5 + F7*bt07*W6
c    $                    + F8*bt08*W7 + F9*bt09*W8
c    $                    + F10*Bt10*W9+F11*Bt11*W10)
c           dF3dW = (ct01 + F2*ct02*W1 + F3*ct03*W2
c    $                    + F4*ct04*W3 + F5*ct05*W4
c    $                    + F6*ct06*W5 + F7*ct07*W6
c    $                    + F8*ct08*W7 + F9*ct09*W8
c    $                    + F10*ct10*W9+F11*ct11*W10)
c           dF4dW = (dt01 + F2*dt02*W1 + F3*dt03*W2
c    $                    + F4*dt04*W3 + F5*dt05*W4
c    $                    + F6*dt06*W5 + F7*dt07*W6
c    $                    + F8*dt08*W7 + F9*dt09*W8
c    $                    + F10*dt10*W9+F11*dt11*W10)
c
c           dWdT = F2/((F1 + Tsig)**F2)
c           dTdR = ((F6*PI*PI)**F2o3)*(rhoo**F2o3)/tauN
c           dTdTau = -TauUEG/tauN**F2
C
c          If(UseLC) then
c            dElSRdR = PDUM
c            dElLRdR = Ax*F4o3*Rho13-PDUM
c          else
c            dElSRdR=Ax*F4o3*Rho13
c            dElLRdR=F0
c          endIf  
c          dGGA1dR = dElSRdR*fx1 + ElSR*dfx1dRho
c          dGGA2dR = dElSRdR*fx2 + ElSR*dfx2dRho 
c          dGGA3dR = dElLRdR*fx1 + ElLR*dfx1dRho
c          dGGA4dR = dElLRdR*fx2 + ElLR*dfx2dRho 
c
c          dF1dR = dF1dW*dWdT*dTdR
c          dF1dTau=dF1dW*dWdT*dTdTau
c          dF2dR = dF2dW*dWdT*dTdR
c          dF2dTau=dF2dW*dWdT*dTdTau
c          dF3dR = dF3dW*dWdT*dTdR
c          dF3dTau=dF3dW*dWdT*dTdTau
c          dF4dR = dF4dW*dWdT*dTdR
c          dF4dTau=dF4dW*dWdT*dTdTau
c
c          dGGA1dG = ElSR*dfx1dG
c          dGGA2dG = ElSR*dfx2dG
c          dGGA3dG = ElLR*dfx1dG
c          dGGA4dG = ElLR*dfx2dG
c
c          Amat(n,D1_RA) = Amat(n,D1_RA) 
c    $                   + (dGGA1dR*Fsig1 + GGA1*dF1dR
c    $                     +dGGA2dR*Fsig2 + GGA2*dF2dR
c    $                     +dGGA3dR*Fsig3 + GGA3*dF3dR
c    $                     +dGGA4dR*Fsig4 + GGA4*dF4dR)*wght
c          Cmat(n,D1_GAA)=  Cmat(n,D1_GAA)
c    $                   + (dGGA1dG*Fsig1 + dGGA2dG*Fsig2
c    $                     +dGGA3dG*Fsig3 + dGGA4dG*Fsig4)*wght
c          Mmat(n,D1_TA) =  Mmat(n,D1_TA)
c    $                   + (GGA1*dF1dTau + GGA2*dF2dTau
c    $                     +GGA3*dF3dTau + GGA4*dF4dTau)*wght
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
           if (rho(n,R_A)+rho(n,R_B).lt.DTol) goto 20
           if (rho(n,R_A).lt.0.5d0*DTol) goto 25           
            rhoo  = rho(n,R_A)
            rho43 = rhoo**F4o3  
c           rho13 = rho43/rhoo
            rho53 = rhoo**F5o3
c       
            tauN = tau(n,T_A)*F2
     
            tauu=tauN
            TAUUEG=F3O5*((F6*PI*PI)**F2O3)*RHO53
c           Tsig =TauUEG/tauu
c           Wsig =(Tsig - F1)/(Tsig + F1)
            Wsig =(TauUEG-tauu)/(TauUEG+tauu)
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

c           Gamma2 =(delrho(n,1,1)*delrho(n,1,1) +
c    &              delrho(n,2,1)*delrho(n,2,1) +
c    &              delrho(n,3,1)*delrho(n,3,1))
            Gamma2 = rgamma(n,G_AA)
c           Gamma = sqrt(Gamma2)
c        
c        X = GAMMA/RHO43
c        S = X/(F48*PI*PI)**F1o3
c        s2     = s*s
         s2 = Gamma2/(rho43*rho43)/(F48*PI*PI)**F2o3
         Deno = (F1 + Mu*s2/kapa)
         fx1=F1+kapa*(F1-F1/Deno)
         fx2=F1+kapas*(F1-Exp(-Mus*s2/kapas))
         If(UseLC) then
#if defined(THIRD_DERIV)
           CALL nwxc_LRCLSDA_d3(EMU,RHOO,ElSR)
#elif defined(SECOND_DERIV)
           CALL nwxc_LRCLSDA_d2(EMU,RHOO,ElSR)
#else
#if defined(NWAD_PRINT)
           CALL nwxc_LRCLSDA_p(EMU,RHOO,ElSR)
#else
           CALL nwxc_LRCLSDA(EMU,RHOO,ElSR)
#endif
#endif
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
c         Ex = Ex +   (GGA1*Fsig1 + GGA2*Fsig2
c    $            +    GGA3*Fsig3 + GGA4*Fsig4)*qwght(n)
          func(n)=func(n)+   (GGA1*Fsig1+GGA2*Fsig2
     $                   +    GGA3*Fsig3+GGA4*Fsig4)*wght
c
c     functional derivatives 
c
c           ds2dRho = -(F8/F3) * s2/rhoo
c           ds2dG = s2/Gamma2
c
c           dfx1ds2 = Mu*(F1/(Deno*Deno)) 
c           dfx1dRho = dfx1ds2*ds2dRho
c           dfx1dG = dfx1ds2*ds2dG
c
c           dfx2ds2 = Mus*DExp(-Mus*s2/kapas)
c           dfx2dRho = dfx2ds2*ds2dRho
c           dfx2dG = dfx2ds2*ds2dG
c
c           dF1dW = (at01 + F2*at02*W1 + F3*at03*W2
c    $                    + F4*at04*W3 + F5*at05*W4
c    $                    + F6*at06*W5 + F7*at07*W6
c    $                    + F8*at08*W7 + F9*at09*W8
c    $                    + F10*at10*W9+F11*at11*W10)
c           dF2dW = (bt01 + F2*bt02*W1 + F3*bt03*W2
c    $                    + F4*bt04*W3 + F5*bt05*W4
c    $                    + F6*bt06*W5 + F7*bt07*W6
c    $                    + F8*bt08*W7 + F9*bt09*W8
c    $                    + F10*Bt10*W9+F11*Bt11*W10)
c           dF3dW = (ct01 + F2*ct02*W1 + F3*ct03*W2
c    $                    + F4*ct04*W3 + F5*ct05*W4
c    $                    + F6*ct06*W5 + F7*ct07*W6
c    $                    + F8*ct08*W7 + F9*ct09*W8
c    $                    + F10*ct10*W9+F11*ct11*W10)
c           dF4dW = (dt01 + F2*dt02*W1 + F3*dt03*W2
c    $                    + F4*dt04*W3 + F5*dt05*W4
c    $                    + F6*dt06*W5 + F7*dt07*W6
c    $                    + F8*dt08*W7 + F9*dt09*W8
c    $                    + F10*dt10*W9+F11*dt11*W10)

c           dWdT = F2/((F1 + Tsig)**F2)
c           dTdR = ((F6*PI*PI)**F2o3)*(rhoo**F2o3)/tauN
c           dTdTau = -TauUEG/tauN**F2
c
c          If(UseLC) then
c            dElSRdR = PDUM
c            dElLRdR = Ax*F4o3*Rho13-PDUM
c          else
c            dElSRdR=Ax*F4o3*Rho13
c            dElLRdR=F0
c          endIf  
c          dGGA1dR = dElSRdR*fx1 + ElSR*dfx1dRho
c          dGGA2dR = dElSRdR*fx2 + ElSR*dfx2dRho 
c          dGGA3dR = dElLRdR*fx1 + ElLR*dfx1dRho
c          dGGA4dR = dElLRdR*fx2 + ElLR*dfx2dRho 
c
c          dF1dR = dF1dW*dWdT*dTdR
c          dF1dTau=dF1dW*dWdT*dTdTau
c          dF2dR = dF2dW*dWdT*dTdR
c          dF2dTau=dF2dW*dWdT*dTdTau
c          dF3dR = dF3dW*dWdT*dTdR
c          dF3dTau=dF3dW*dWdT*dTdTau
c          dF4dR = dF4dW*dWdT*dTdR
c          dF4dTau=dF4dW*dWdT*dTdTau
c
c          dGGA1dG = ElSR*dfx1dG
c          dGGA2dG = ElSR*dfx2dG
c          dGGA3dG = ElLR*dfx1dG
c          dGGA4dG = ElLR*dfx2dG
c
c          Amat(n,D1_RA) = Amat(n,D1_RA)
c    $                   + (dGGA1dR*Fsig1 + GGA1*dF1dR
c    $                     +dGGA2dR*Fsig2 + GGA2*dF2dR
c    $                     +dGGA3dR*Fsig3 + GGA3*dF3dR
c    $                     +dGGA4dR*Fsig4 + GGA4*dF4dR)*wght
c          Cmat(n,D1_GAA)= Cmat(n,D1_GAA)
c    $                   + (dGGA1dG*Fsig1 + dGGA2dG*Fsig2
c    $                     +dGGA3dG*Fsig3 + dGGA4dG*Fsig4)*wght
c          Mmat(n,D1_TA) = Mmat(n,D1_TA)
c    $                   + (GGA1*dF1dTau  + GGA2*dF2dTau
c    $                     +GGA3*dF3dTau  + GGA4*dF4dTau)*wght
c
25         continue
c
c     Beta               BETA           BETA
c
            if (rho(n,R_B).lt.0.5d0*DTol) goto 20
            rhoo  = rho(n,R_B)
            rho43 = rhoo**F4o3  
c           rho13 = rho43/rhoo
            rho53 = rhoo**F5o3
c
             
            tauN = tau(n,T_B)*F2
     
            tauu=tauN
            TAUUEG=F3O5*((F6*PI*PI)**F2O3)*RHO53
c           Tsig =TauUEG/tauu
c           Wsig =(Tsig - F1)/(Tsig + F1)
            Wsig =(TauUEG-tauu)/(TauUEG+tauu)
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

c           Gamma2 =(delrho(n,1,2)*delrho(n,1,2) +
c    &              delrho(n,2,2)*delrho(n,2,2) +
c    &              delrho(n,3,2)*delrho(n,3,2))
            Gamma2 = rgamma(n,G_BB)
c           Gamma = sqrt(Gamma2)
c        
c        X = GAMMA/RHO43
c        S = X/(F48*PI*PI)**F1o3
c        s2     = s*s
         s2 = Gamma2/(rho43*rho43)/(F48*PI*PI)**F2o3
         Deno = (F1 + Mu*s2/kapa)
         fx1=F1+kapa*(F1-F1/Deno)
         fx2=F1+kapas*(F1-Exp(-Mus*s2/kapas))
         If(UseLC) then
#if defined(THIRD_DERIV)
           CALL nwxc_LRCLSDA_d3(EMU,RHOO,ElSR)
#elif defined(SECOND_DERIV)
           CALL nwxc_LRCLSDA_d2(EMU,RHOO,ElSR)
#else
#if defined(NWAD_PRINT)
           CALL nwxc_LRCLSDA_p(EMU,RHOO,ElSR)
#else
           CALL nwxc_LRCLSDA(EMU,RHOO,ElSR)
#endif
#endif
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
c         Ex = Ex +   (GGA1*Fsig1 + GGA2*Fsig2
c    $            +    GGA3*Fsig3 + GGA4*Fsig4)*qwght(n)
          func(n)=func(n)+   (GGA1*Fsig1+GGA2*Fsig2
     $                   +    GGA3*Fsig3+GGA4*Fsig4)*wght
c
c     functional derivatives 
c
c           ds2dRho = -(F8/F3) * s2/rhoo
c           ds2dG = s2/Gamma2
c
c           dfx1ds2 = Mu*(F1/(Deno*Deno)) 
c           dfx1dRho = dfx1ds2*ds2dRho
c           dfx1dG = dfx1ds2*ds2dG
c
c           dfx2ds2 = Mus*DExp(-Mus*s2/kapas)
c           dfx2dRho = dfx2ds2*ds2dRho
c           dfx2dG = dfx2ds2*ds2dG
c
c           dF1dW = (at01 + F2*at02*W1 + F3*at03*W2
c    $                    + F4*at04*W3 + F5*at05*W4
c    $                    + F6*at06*W5 + F7*at07*W6
c    $                    + F8*at08*W7 + F9*at09*W8
c    $                    + F10*at10*W9+F11*at11*W10)
c           dF2dW = (bt01 + F2*bt02*W1 + F3*bt03*W2
c    $                    + F4*bt04*W3 + F5*bt05*W4
c    $                    + F6*bt06*W5 + F7*bt07*W6
c    $                    + F8*bt08*W7 + F9*bt09*W8
c    $                    + F10*Bt10*W9+F11*Bt11*W10)
c           dF3dW = (ct01 + F2*ct02*W1 + F3*ct03*W2
c    $                    + F4*ct04*W3 + F5*ct05*W4
c    $                    + F6*ct06*W5 + F7*ct07*W6
c    $                    + F8*ct08*W7 + F9*ct09*W8
c    $                    + F10*ct10*W9+F11*ct11*W10)
c           dF4dW = (dt01 + F2*dt02*W1 + F3*dt03*W2
c    $                    + F4*dt04*W3 + F5*dt05*W4
c    $                    + F6*dt06*W5 + F7*dt07*W6
c    $                    + F8*dt08*W7 + F9*dt09*W8
c    $                    + F10*dt10*W9+F11*dt11*W10)

c           dWdT = F2/((F1 + Tsig)**F2)
c           dTdR = ((F6*PI*PI)**F2o3)*(rhoo**F2o3)/tauN
c           dTdTau = -TauUEG/tauN**F2
C
c          If(UseLC) then
c            dElSRdR = PDUM
c            dElLRdR = Ax*F4o3*Rho13-PDUM
c          else
c            dElSRdR=Ax*F4o3*Rho13
c            dElLRdR=F0
c          endIf  
c          dGGA1dR = dElSRdR*fx1 + ElSR*dfx1dRho
c          dGGA2dR = dElSRdR*fx2 + ElSR*dfx2dRho 
c          dGGA3dR = dElLRdR*fx1 + ElLR*dfx1dRho
c          dGGA4dR = dElLRdR*fx2 + ElLR*dfx2dRho 
c
c          dF1dR = dF1dW*dWdT*dTdR
c          dF1dTau=dF1dW*dWdT*dTdTau
c          dF2dR = dF2dW*dWdT*dTdR
c          dF2dTau=dF2dW*dWdT*dTdTau
c          dF3dR = dF3dW*dWdT*dTdR
c          dF3dTau=dF3dW*dWdT*dTdTau
c          dF4dR = dF4dW*dWdT*dTdR
c          dF4dTau=dF4dW*dWdT*dTdTau
c
c          dGGA1dG = ElSR*dfx1dG
c          dGGA2dG = ElSR*dfx2dG
c          dGGA3dG = ElLR*dfx1dG
c          dGGA4dG = ElLR*dfx2dG
c
c          Amat(n,D1_RB) = Amat(n,D1_RB)
c    $                   + (dGGA1dR*Fsig1 + GGA1*dF1dR
c    $                     +dGGA2dR*Fsig2 + GGA2*dF2dR
c    $                     +dGGA3dR*Fsig3 + GGA3*dF3dR
c    $                     +dGGA4dR*Fsig4 + GGA4*dF4dR)*wght
c          Cmat(n,D1_GBB)= Cmat(n,D1_GBB)
c    $                   + (dGGA1dG*Fsig1 + dGGA2dG*Fsig2
c    $                     +dGGA3dG*Fsig3 + dGGA4dG*Fsig4)*wght
c          Mmat(n,D1_TB) = Mmat(n,D1_TB)
c    $                   + (GGA1*dF1dTau  + GGA2*dF2dTau
c    $                     +GGA3*dF3dTau  + GGA4*dF4dTau)*wght
c
20      continue
      endif
      return
      end
c
C>
C> \brief The M11 exchange functionals
C>
C> This routine implements exchange meta-GGA functionals 
C> of the M11 [1,2,3] family of functionals.
C>
C> ### References ###
C>
C> [1] Y. Zhao, D. G. Truhlar,
C>     "Exploring the limit of accuracy of the global hybrid meta 
C>     density functional for main-group thermochemistry, kinetics,
C>     and noncovalent interactions",
C>     J. Chem. Theory Comput. <b>4</b>, 1849-1868 (2008), DOI:
C>     <a href="https://doi.org/10.1021/ct800246v">
C>     10.1021/ct800246v</a>.
C>
C> [2] R. Peverati, D. G. Truhlar,
C>     "Improving the accuracy of hybrid meta-GGA density functionals
C>     by range separation",
C>     J. Phys. Chem. Lett. <b>2</b>, 2810-2817 (2011), DOI:
C>     <a href="https://doi.org/10.1021/jz201170d">
C>     10.1021/jz201170d</a>.
C>
C> [3] R. Peverati, D. G. Truhlar,
C>     "M11-L: A local density functional that provides improved
C>     accuracy for electronic structure calculations in chemistry
C>     and physics",
C>     J. Phys. Chem. Lett. <b>3</b>, 117-124 (2011), DOI:
C>     <a href="https://doi.org/10.1021/jz201525m">
C>     10.1021/jz201525m</a>.
C>
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
      Subroutine nwxc_x_m11_p(param, tol_rho, ipol, nq, wght, rho,
     &                        rgamma, tau, func)
#else
      Subroutine nwxc_x_m11(param, tol_rho, ipol, nq, wght, rho,
     &                      rgamma, tau, func)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      Subroutine nwxc_x_m11_d2(param, tol_rho, ipol, nq, wght, rho,
     &                         rgamma, tau, func)
#else
      Subroutine nwxc_x_m11_d3(param, tol_rho, ipol, nq, wght, rho,
     &                         rgamma, tau, func)
#endif
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
#include "nwad.fh"
c
      implicit none
c
#include "nwxc_param.fh"
c
#if defined(NWAD_PRINT)
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      type(nwad_dble)::param(*) !< [Input] Parameters of the functional
      type(nwad_dble)::at00,at01,at02,at03,at04,at05,at06
      type(nwad_dble)::at07,at08,at09,at10,at11
      type(nwad_dble)::bt00,bt01,bt02,bt03,bt04,bt05,bt06
      type(nwad_dble)::bt07,bt08,bt09,bt10,bt11
      type(nwad_dble)::ct00,ct01,ct02,ct03,ct04,ct05,ct06
      type(nwad_dble)::ct07,ct08,ct09,ct10,ct11
      type(nwad_dble)::dt00,dt01,dt02,dt03,dt04,dt05,dt06
      type(nwad_dble)::dt07,dt08,dt09,dt10,dt11
#else
      double precision param(*) !< [Input] Parameters of the functional
      double precision at00,at01,at02,at03,at04,at05,at06
      double precision at07,at08,at09,at10,at11
      double precision bt00,bt01,bt02,bt03,bt04,bt05,bt06
      double precision bt07,bt08,bt09,bt10,bt11
      double precision ct00,ct01,ct02,ct03,ct04,ct05,ct06
      double precision ct07,ct08,ct09,ct10,ct11
      double precision dt00,dt01,dt02,dt03,dt04,dt05,dt06
      double precision dt07,dt08,dt09,dt10,dt11
#endif
#else
      double precision param(*) !< [Input] Parameters of the functional
                                !< (see Table 2 of [1])
                                !< - param(1): \f$ a_0 \f$
                                !< - param(2): \f$ a_1 \f$
                                !< - param(3): \f$ a_2 \f$
                                !< - param(4): \f$ a_3 \f$
                                !< - param(5): \f$ a_4 \f$
                                !< - param(6): \f$ a_5 \f$
                                !< - param(7): \f$ a_6 \f$
                                !< - param(8): \f$ a_7 \f$
                                !< - param(9): \f$ a_8 \f$
                                !< - param(10): \f$ a_9 \f$
                                !< - param(11): \f$ a_{10} \f$
                                !< - param(12): \f$ a_{11} \f$
                                !< - param(13): \f$ b_0 \f$
                                !< - param(14): \f$ b_1 \f$
                                !< - param(15): \f$ b_2 \f$
                                !< - param(16): \f$ b_3 \f$
                                !< - param(17): \f$ b_4 \f$
                                !< - param(18): \f$ b_5 \f$
                                !< - param(19): \f$ b_6 \f$
                                !< - param(20): \f$ b_7 \f$
                                !< - param(21): \f$ b_8 \f$
                                !< - param(22): \f$ b_9 \f$
                                !< - param(23): \f$ b_{10} \f$
                                !< - param(24): \f$ b_{11} \f$
                                !< - param(25): \f$ c_0 \f$
                                !< - param(26): \f$ c_1 \f$
                                !< - param(27): \f$ c_2 \f$
                                !< - param(28): \f$ c_3 \f$
                                !< - param(29): \f$ c_4 \f$
                                !< - param(30): \f$ c_5 \f$
                                !< - param(31): \f$ c_6 \f$
                                !< - param(32): \f$ c_7 \f$
                                !< - param(33): \f$ c_8 \f$
                                !< - param(34): \f$ c_9 \f$
                                !< - param(35): \f$ c_{10} \f$
                                !< - param(36): \f$ c_{11} \f$
                                !< - param(37): \f$ d_0 \f$
                                !< - param(38): \f$ d_1 \f$
                                !< - param(39): \f$ d_2 \f$
                                !< - param(40): \f$ d_3 \f$
                                !< - param(41): \f$ d_4 \f$
                                !< - param(42): \f$ d_5 \f$
                                !< - param(43): \f$ d_6 \f$
                                !< - param(44): \f$ d_7 \f$
                                !< - param(45): \f$ d_8 \f$
                                !< - param(46): \f$ d_9 \f$
                                !< - param(47): \f$ d_{10} \f$
                                !< - param(48): \f$ d_{11} \f$
      double precision at00,at01,at02,at03,at04,at05,at06
      double precision at07,at08,at09,at10,at11
      double precision bt00,bt01,bt02,bt03,bt04,bt05,bt06
      double precision bt07,bt08,bt09,bt10,bt11
      double precision ct00,ct01,ct02,ct03,ct04,ct05,ct06
      double precision ct07,ct08,ct09,ct10,ct11
      double precision dt00,dt01,dt02,dt03,dt04,dt05,dt06
      double precision dt07,dt08,dt09,dt10,dt11
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
      double precision pi,DTol
      type(nwad_dble)::tauN,tauu
c
c      functional derivatives
c
      double precision dWdT, dTdR, dTdTau
c
c     Intermediate derivative results, etc.
c
      integer n
c
      double precision Ax
      type(nwad_dble)::s, s2
      double precision kapa,kapas,mu,mus
c
      double precision f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11
      double precision F1o3,F2o3,F3o5,F4o3,F5o3,F48,F81 
      type(nwad_dble)::Fsig1, Fsig2, Fsig3, Fsig4, Fx1, Fx2
      type(nwad_dble)::ElSR, ElLR
      double precision PDUM
      type(nwad_dble)::GGA1, GGA2, GGA3, GGA4
      double precision Emu
      type(nwad_dble)::X, deno
      double precision ds2drho, ds2dg, dfx1ds2
      double precision dfx2ds2, df1dw
      double precision dfx1drho,dfx1dg,dfx2drho,dfx2dg,df2dw
      double precision df3dw, df4dw, delsrdr, dellrdr
      double precision dgga1dr, dgga2dr, dgga3dr, dgga4dr
      double precision df1dr, df1dtau, df2dr, df2dtau
      double precision df3dr, df3dtau, df4dr, df4dtau
      double precision dgga1dg, dgga2dg, dgga3dg, dgga4dg
c
      type(nwad_dble)::rho43, rhoo, rho53
      type(nwad_dble)::Gamma2, Gamma
      type(nwad_dble)::TauUEG, Tsig, Wsig, W1, W2, W3, W4, W5, W6
      type(nwad_dble)::W7, W8, W9, W10, W11
      logical UseLC
c
      parameter( F0=0.0D+00,  F1=1.0D+00,  F2=2.0D+00,
     $           F3=3.0D+00,  F4=4.0D+00,  F5=5.0D+00,
     $           F6=6.0D+00,  F7=7.0D+00,  F8=8.0D+00,
     $           F9=9.0D+00,  F10=10.0D+00,F11=11.0D+00)
c
        pi=acos(-1d0)      
c
        at00=param( 1)
        at01=param( 2)
        at02=param( 3)
        at03=param( 4)
        at04=param( 5)
        at05=param( 6)
        at06=param( 7)
        at07=param( 8)
        at08=param( 9)
        at09=param(10)
        at10=param(11) 
        at11=param(12) 
c
        bt00=param(13)
        bt01=param(14)
        bt02=param(15)
        bt03=param(16)
        bt04=param(17)
        bt05=param(18)
        bt06=param(19)
        bt07=param(20)
        bt08=param(21)
        bt09=param(22)
        bt10=param(23) 
        bt11=param(24) 
c
        ct00=param(25)
        ct01=param(26)
        ct02=param(27)
        ct03=param(28)
        ct04=param(29)
        ct05=param(30)
        ct06=param(31)
        ct07=param(32)
        ct08=param(33)
        ct09=param(34)
        ct10=param(35)
        ct11=param(36)
C
        dt00=param(37)
        dt01=param(38)
        dt02=param(39)
        dt03=param(40)
        dt04=param(41)
        dt05=param(42)
        dt06=param(43)
        dt07=param(44)
        dt08=param(45)
        dt09=param(46)
        dt10=param(47)
        dt11=param(48)
c
        UseLC=.True.
        Emu  = 0.25d0
C
c     if (ijzy.eq.1) then
C     Parameters for M08-HX
c       at00=  1.3340172D+00
c       at01= -9.4751087D+00
c       at02= -1.2541893D+01
c       at03=  9.1369974D+00
c       at04=  3.4717204D+01
c       at05=  5.8831807D+01
c       at06=  7.1369574D+01
c       at07=  2.3312961D+01
c       at08=  4.8314679D+00
c       at09= -6.5044167D+00
c       at10= -1.4058265D+01
c       at11=  1.2880570D+01

c       bt00= -8.5631823D-01
c       bt01=  9.2810354D+00
c       bt02=  1.2260749D+01
c       bt03= -5.5189665D+00
c       bt04= -3.5534989D+01
c       bt05= -8.2049996D+01
c       bt06= -6.8586558D+01
c       bt07=  3.6085694D+01
c       bt08= -9.3740983D+00
c       bt09= -5.9731688D+01
c       bt10=  1.6587868D+01
c       bt11=  1.3993203D+01
c
c       UseLC=.False.
c
c      elseif (ijzy.eq.2) then
c     Parameters for M08-SO
c       at00= -3.4888428D-01
c       at01= -5.8157416D+00
c       at02=  3.7550810D+01
c       at03=  6.3727406D+01
c       at04= -5.3742313D+01
c       at05= -9.8595529D+01
c       at06=  1.6282216D+01
c       at07=  1.7513468D+01
c       at08= -6.7627553D+00
c       at09=  1.1106658D+01
c       at10=  1.5663545D+00
c       at11=  8.7603470D+00

c       bt00=  7.8098428D-01
c       bt01=  5.4538178D+00
c       bt02= -3.7853348D+01
c       bt03= -6.2295080D+01
c       bt04=  4.6713254D+01
c       bt05=  8.7321376D+01
c       bt06=  1.6053446D+01
c       bt07=  2.0126920D+01
c       bt08= -4.0343695D+01
c       bt09= -5.8577565D+01
c       bt10=  2.0890272D+01
c       bt11=  1.0946903D+01
c
c       UseLC=.False.
c
c     elseif (ijzy.eq.3) then
c     Parameters for M11
c       at00= -0.18399900D+00
c       at01= -1.39046703D+01
c       at02=  1.18206837D+01
c       at03=  3.10098465D+01
c       at04= -5.19625696D+01
c       at05=  1.55750312D+01
c       at06= -6.94775730D+00
c       at07= -1.58465014D+02
c       at08= -1.48447565D+00
c       at09=  5.51042124D+01
c       at10= -1.34714184D+01
c       at11=  0.00000000D+00

c       bt00=  0.75599900D+00
c       bt01=  1.37137944D+01
c       bt02= -1.27998304D+01
c       bt03= -2.93428814D+01
c       bt04=  5.91075674D+01
c       bt05= -2.27604866D+01
c       bt06= -1.02769340D+01
c       bt07=  1.64752731D+02
c       bt08=  1.85349258D+01
c       bt09= -5.56825639D+01
c       bt10=  7.47980859D+00
c       bt11=  0.00000000D+00
c
c       UseLC=.True.
c       Emu =0.25D+00
c
c     elseif (ijzy.eq.4) then
c     Parameters for M11-L
c       at00=  8.121131D-01
c       at01=  1.738124D+01
c       at02=  1.154007D+00
c       at03=  6.869556D+01
c       at04=  1.016864D+02
c       at05= -5.887467D+00
c       at06=  4.517409D+01
c       at07= -2.773149D+00
c       at08= -2.617211D+01
c       at09=  0.000000D+00
c       at10=  0.000000D+00 
c       at11=  0.000000D+00
c
c       bt00=  1.878869D-01
c       bt01= -1.653877D+01
c       bt02=  6.755753D-01
c       bt03= -7.567572D+01
c       bt04= -1.040272D+02
c       bt05=  1.831853D+01
c       bt06= -5.573352D+01
c       bt07= -3.520210D+00
c       bt08=  3.724276D+01
c       bt09=  0.000000D+00
c       bt10=  0.000000D+00
c       bt11=  0.000000D+00
c
c       ct00= -4.386615D-01
c       ct01= -1.214016D+02
c       ct02= -1.393573D+02
c       ct03= -2.046649D+00
c       ct04=  2.804098D+01
c       ct05= -1.312258D+01
c       ct06= -6.361819D+00
c       ct07= -8.055758D-01
c       ct08=  3.736551D+00
c       ct09=  0.000000D+00
c       ct10=  0.000000D+00
c       ct11=  0.000000D+00
c
c       dt00=  1.438662D+00
c       dt01=  1.209465D+02
c       dt02=  1.328252D+02
c       dt03=  1.296355D+01
c       dt04=  5.854866D+00
c       dt05= -3.378162D+00
c       dt06= -4.423393D+01
c       dt07=  6.844475D+00
c       dt08=  1.949541D+01
c       dt09=  0.000000D+00
c       dt10=  0.000000D+00
c       dt11=  0.000000D+00
c
c       UseLC=.True.
c       Emu =0.25D+00
c
c     endif
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
            if (rho(n,R_T).lt.DTol) goto 10
            rhoo = rho(n,R_T)/F2
            rho43 = rhoo**F4o3  
c           rho13 = rhoo**F1o3
            rho53 = rhoo**F5o3
c            
            tauN = tau(n,T_T)
            tauu=tauN
            TAUUEG=F3O5*((F6*PI*PI)**F2O3)*RHO53
c           Tsig =TauUEG/tauu
c           Wsig =(Tsig - F1)/(Tsig + F1)
            Wsig =(TauUEG-tauu)/(TauUEG+tauu)
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

c           Gamma2 =(delrho(n,1,1)*delrho(n,1,1) +
c    &              delrho(n,2,1)*delrho(n,2,1) +
c    &              delrho(n,3,1)*delrho(n,3,1))/F4
            Gamma2 = rgamma(n,G_TT)/F4
c           Gamma = sqrt(Gamma2)
c        
c        X = GAMMA/RHO43
c        S = X/(F48*PI*PI)**F1o3
c        s2     = s*s
         s2 = Gamma2/(rho43*rho43)/(F48*PI*PI)**F2o3
         Deno = (F1 + Mu*s2/kapa)
         fx1=F1+kapa*(F1-F1/Deno)
         fx2=F1+kapas*(F1-Exp(-Mus*s2/kapas))
         If(UseLC) then
#if defined(THIRD_DERIV)
           CALL nwxc_LRCLSDA_d3(EMU,RHOO,ElSR)
#elif defined(SECOND_DERIV)
           CALL nwxc_LRCLSDA_d2(EMU,RHOO,ElSR)
#else
#if defined(NWAD_PRINT)
           CALL nwxc_LRCLSDA_p(EMU,RHOO,ElSR)
#else
           CALL nwxc_LRCLSDA(EMU,RHOO,ElSR)
#endif
#endif
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
c        Ex = Ex +F2*(GGA1*Fsig1 + GGA2*Fsig2
c    $           +    GGA3*Fsig3 + GGA4*Fsig4)*qwght(n)
         func(n)=func(n)+F2*(GGA1*Fsig1+GGA2*Fsig2
     $                  +    GGA3*Fsig3+GGA4*Fsig4)*wght

c
c     functional derivatives 
c
c           ds2dRho = -(F8/F3) * s2/rhoo
c           ds2dG = s2/Gamma2
c
c           dfx1ds2 = Mu*(F1/(Deno*Deno)) 
c           dfx1dRho = dfx1ds2*ds2dRho
c           dfx1dG = dfx1ds2*ds2dG
c
c           dfx2ds2 = Mus*DExp(-Mus*s2/kapas)
c           dfx2dRho = dfx2ds2*ds2dRho
c           dfx2dG = dfx2ds2*ds2dG
c
c           dF1dW = (at01 + F2*at02*W1 + F3*at03*W2
c    $                    + F4*at04*W3 + F5*at05*W4
c    $                    + F6*at06*W5 + F7*at07*W6
c    $                    + F8*at08*W7 + F9*at09*W8
c    $                    + F10*at10*W9+F11*at11*W10)
c           dF2dW = (bt01 + F2*bt02*W1 + F3*bt03*W2
c    $                    + F4*bt04*W3 + F5*bt05*W4
c    $                    + F6*bt06*W5 + F7*bt07*W6
c    $                    + F8*bt08*W7 + F9*bt09*W8
c    $                    + F10*Bt10*W9+F11*Bt11*W10)
c           dF3dW = (ct01 + F2*ct02*W1 + F3*ct03*W2
c    $                    + F4*ct04*W3 + F5*ct05*W4
c    $                    + F6*ct06*W5 + F7*ct07*W6
c    $                    + F8*ct08*W7 + F9*ct09*W8
c    $                    + F10*ct10*W9+F11*ct11*W10)
c           dF4dW = (dt01 + F2*dt02*W1 + F3*dt03*W2
c    $                    + F4*dt04*W3 + F5*dt05*W4
c    $                    + F6*dt06*W5 + F7*dt07*W6
c    $                    + F8*dt08*W7 + F9*dt09*W8
c    $                    + F10*dt10*W9+F11*dt11*W10)
c
c           dWdT = F2/((F1 + Tsig)**F2)
c           dTdR = ((F6*PI*PI)**F2o3)*(rhoo**F2o3)/tauN
c           dTdTau = -TauUEG/tauN**F2
c
c          If(UseLC) then
c            dElSRdR = PDUM
c            dElLRdR = Ax*F4o3*Rho13-PDUM
c          else
c            dElSRdR=Ax*F4o3*Rho13
c            dElLRdR=F0
c          endIf  
c          dGGA1dR = dElSRdR*fx1 + ElSR*dfx1dRho
c          dGGA2dR = dElSRdR*fx2 + ElSR*dfx2dRho 
c          dGGA3dR = dElLRdR*fx1 + ElLR*dfx1dRho
c          dGGA4dR = dElLRdR*fx2 + ElLR*dfx2dRho 
c
c          dF1dR = dF1dW*dWdT*dTdR
c          dF1dTau=dF1dW*dWdT*dTdTau
c          dF2dR = dF2dW*dWdT*dTdR
c          dF2dTau=dF2dW*dWdT*dTdTau
c          dF3dR = dF3dW*dWdT*dTdR
c          dF3dTau=dF3dW*dWdT*dTdTau
c          dF4dR = dF4dW*dWdT*dTdR
c          dF4dTau=dF4dW*dWdT*dTdTau
c
c          dGGA1dG = ElSR*dfx1dG
c          dGGA2dG = ElSR*dfx2dG
c          dGGA3dG = ElLR*dfx1dG
c          dGGA4dG = ElLR*dfx2dG
c
c          Amat(n,D1_RA) = Amat(n,D1_RA) 
c    $                   + (dGGA1dR*Fsig1 + GGA1*dF1dR
c    $                     +dGGA2dR*Fsig2 + GGA2*dF2dR
c    $                     +dGGA3dR*Fsig3 + GGA3*dF3dR
c    $                     +dGGA4dR*Fsig4 + GGA4*dF4dR)*wght
c          Cmat(n,D1_GAA)=  Cmat(n,D1_GAA)
c    $                   + (dGGA1dG*Fsig1 + dGGA2dG*Fsig2
c    $                     +dGGA3dG*Fsig3 + dGGA4dG*Fsig4)*wght
c          Mmat(n,D1_TA) =  Mmat(n,D1_TA)
c    $                   + (GGA1*dF1dTau + GGA2*dF2dTau
c    $                     +GGA3*dF3dTau + GGA4*dF4dTau)*wght
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
           if (rho(n,R_A)+rho(n,R_B).lt.DTol) goto 20
           if (rho(n,R_A).lt.0.5d0*DTol) goto 25           
            rhoo  = rho(n,R_A)
            rho43 = rhoo**F4o3  
c           rho13 = rhoo**F1o3
            rho53 = rhoo**F5o3
c       
            tauN = tau(n,T_A)*F2
     
            tauu=tauN
            TAUUEG=F3O5*((F6*PI*PI)**F2O3)*RHO53
c           Tsig =TauUEG/tauu
c           Wsig =(Tsig - F1)/(Tsig + F1)
            Wsig = (TauUEG-tauu)/(TauUEG+tauu)
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

c           Gamma2 =(delrho(n,1,1)*delrho(n,1,1) +
c    &              delrho(n,2,1)*delrho(n,2,1) +
c    &              delrho(n,3,1)*delrho(n,3,1))
            Gamma2 = rgamma(n,G_AA)
c           Gamma = sqrt(Gamma2)
c        
c        X = GAMMA/RHO43
c        S = X/(F48*PI*PI)**F1o3
c        s2     = s*s
         s2 = Gamma2/(rho43*rho43)/(F48*PI*PI)**F2o3
         Deno = (F1 + Mu*s2/kapa)
         fx1=F1+kapa*(F1-F1/Deno)
         fx2=F1+kapas*(F1-Exp(-Mus*s2/kapas))
         If(UseLC) then
#if defined(THIRD_DERIV)
           CALL nwxc_LRCLSDA_d3(EMU,RHOO,ElSR)
#elif defined(SECOND_DERIV)
           CALL nwxc_LRCLSDA_d2(EMU,RHOO,ElSR)
#else
#if defined(NWAD_PRINT)
           CALL nwxc_LRCLSDA_p(EMU,RHOO,ElSR)
#else
           CALL nwxc_LRCLSDA(EMU,RHOO,ElSR)
#endif
#endif
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
c         Ex = Ex +   (GGA1*Fsig1 + GGA2*Fsig2
c    $            +    GGA3*Fsig3 + GGA4*Fsig4)*qwght(n)
          func(n)=func(n)+   (GGA1*Fsig1+GGA2*Fsig2
     $                   +    GGA3*Fsig3+GGA4*Fsig4)*wght
c
c     functional derivatives 
c
c           ds2dRho = -(F8/F3) * s2/rhoo
c           ds2dG = s2/Gamma2
c
c           dfx1ds2 = Mu*(F1/(Deno*Deno)) 
c           dfx1dRho = dfx1ds2*ds2dRho
c           dfx1dG = dfx1ds2*ds2dG
c
c           dfx2ds2 = Mus*DExp(-Mus*s2/kapas)
c           dfx2dRho = dfx2ds2*ds2dRho
c           dfx2dG = dfx2ds2*ds2dG
c
c           dF1dW = (at01 + F2*at02*W1 + F3*at03*W2
c    $                    + F4*at04*W3 + F5*at05*W4
c    $                    + F6*at06*W5 + F7*at07*W6
c    $                    + F8*at08*W7 + F9*at09*W8
c    $                    + F10*at10*W9+F11*at11*W10)
c           dF2dW = (bt01 + F2*bt02*W1 + F3*bt03*W2
c    $                    + F4*bt04*W3 + F5*bt05*W4
c    $                    + F6*bt06*W5 + F7*bt07*W6
c    $                    + F8*bt08*W7 + F9*bt09*W8
c    $                    + F10*Bt10*W9+F11*Bt11*W10)
c           dF3dW = (ct01 + F2*ct02*W1 + F3*ct03*W2
c    $                    + F4*ct04*W3 + F5*ct05*W4
c    $                    + F6*ct06*W5 + F7*ct07*W6
c    $                    + F8*ct08*W7 + F9*ct09*W8
c    $                    + F10*ct10*W9+F11*ct11*W10)
c           dF4dW = (dt01 + F2*dt02*W1 + F3*dt03*W2
c    $                    + F4*dt04*W3 + F5*dt05*W4
c    $                    + F6*dt06*W5 + F7*dt07*W6
c    $                    + F8*dt08*W7 + F9*dt09*W8
c    $                    + F10*dt10*W9+F11*dt11*W10)

c           dWdT = F2/((F1 + Tsig)**F2)
c           dTdR = ((F6*PI*PI)**F2o3)*(rhoo**F2o3)/tauN
c           dTdTau = -TauUEG/tauN**F2
c
c          If(UseLC) then
c            dElSRdR = PDUM
c            dElLRdR = Ax*F4o3*Rho13-PDUM
c          else
c            dElSRdR=Ax*F4o3*Rho13
c            dElLRdR=F0
c          endIf  
c          dGGA1dR = dElSRdR*fx1 + ElSR*dfx1dRho
c          dGGA2dR = dElSRdR*fx2 + ElSR*dfx2dRho 
c          dGGA3dR = dElLRdR*fx1 + ElLR*dfx1dRho
c          dGGA4dR = dElLRdR*fx2 + ElLR*dfx2dRho 
c
c          dF1dR = dF1dW*dWdT*dTdR
c          dF1dTau=dF1dW*dWdT*dTdTau
c          dF2dR = dF2dW*dWdT*dTdR
c          dF2dTau=dF2dW*dWdT*dTdTau
c          dF3dR = dF3dW*dWdT*dTdR
c          dF3dTau=dF3dW*dWdT*dTdTau
c          dF4dR = dF4dW*dWdT*dTdR
c          dF4dTau=dF4dW*dWdT*dTdTau
c
c          dGGA1dG = ElSR*dfx1dG
c          dGGA2dG = ElSR*dfx2dG
c          dGGA3dG = ElLR*dfx1dG
c          dGGA4dG = ElLR*dfx2dG
c
c          Amat(n,D1_RA) = Amat(n,D1_RA)
c    $                   + (dGGA1dR*Fsig1 + GGA1*dF1dR
c    $                     +dGGA2dR*Fsig2 + GGA2*dF2dR
c    $                     +dGGA3dR*Fsig3 + GGA3*dF3dR
c    $                     +dGGA4dR*Fsig4 + GGA4*dF4dR)*wght
c          Cmat(n,D1_GAA)= Cmat(n,D1_GAA)
c    $                   + (dGGA1dG*Fsig1 + dGGA2dG*Fsig2
c    $                     +dGGA3dG*Fsig3 + dGGA4dG*Fsig4)*wght
c          Mmat(n,D1_TA) = Mmat(n,D1_TA)
c    $                   + (GGA1*dF1dTau  + GGA2*dF2dTau
c    $                     +GGA3*dF3dTau  + GGA4*dF4dTau)*wght
c
25         continue
c
c     Beta               BETA           BETA
c
            if (rho(n,R_B).lt.0.5d0*DTol) goto 20
            rhoo  = rho(n,R_B)
            rho43 = rhoo**F4o3  
c           rho13 = rhoo**F1o3
            rho53 = rhoo**F5o3
c
             
            tauN = tau(n,T_B)*F2
     
            tauu=tauN
            TAUUEG=F3O5*((F6*PI*PI)**F2O3)*RHO53
c           Tsig =TauUEG/tauu
c           Wsig =(Tsig - F1)/(Tsig + F1)
            Wsig = (TauUEG-tauu)/(TauUEG+tauu)
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

c           Gamma2 =(delrho(n,1,2)*delrho(n,1,2) +
c    &              delrho(n,2,2)*delrho(n,2,2) +
c    &              delrho(n,3,2)*delrho(n,3,2))
            Gamma2 = rgamma(n,G_BB)
c           Gamma = sqrt(Gamma2)
c        
c        X = GAMMA/RHO43
c        S = X/(F48*PI*PI)**F1o3
c        s2     = s*s
         s2 = Gamma2/(rho43*rho43)/(F48*PI*PI)**F2o3
         Deno = (F1 + Mu*s2/kapa)
         fx1=F1+kapa*(F1-F1/Deno)
         fx2=F1+kapas*(F1-Exp(-Mus*s2/kapas))
         If(UseLC) then
#if defined(THIRD_DERIV)
           CALL nwxc_LRCLSDA_d3(EMU,RHOO,ElSR)
#elif defined(SECOND_DERIV)
           CALL nwxc_LRCLSDA_d2(EMU,RHOO,ElSR)
#else
#if defined(NWAD_PRINT)
           CALL nwxc_LRCLSDA_p(EMU,RHOO,ElSR)
#else
           CALL nwxc_LRCLSDA(EMU,RHOO,ElSR)
#endif
#endif
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
c         Ex = Ex +   (GGA1*Fsig1 + GGA2*Fsig2
c    $            +    GGA3*Fsig3 + GGA4*Fsig4)*qwght(n)
          func(n)=func(n)+   (GGA1*Fsig1+GGA2*Fsig2
     $                   +    GGA3*Fsig3+GGA4*Fsig4)*wght
c
c     functional derivatives 
c
c           ds2dRho = -(F8/F3) * s2/rhoo
c           ds2dG = s2/Gamma2
c
c           dfx1ds2 = Mu*(F1/(Deno*Deno)) 
c           dfx1dRho = dfx1ds2*ds2dRho
c           dfx1dG = dfx1ds2*ds2dG
C
c           dfx2ds2 = Mus*DExp(-Mus*s2/kapas)
c           dfx2dRho = dfx2ds2*ds2dRho
c           dfx2dG = dfx2ds2*ds2dG
c
c           dF1dW = (at01 + F2*at02*W1 + F3*at03*W2
c    $                    + F4*at04*W3 + F5*at05*W4
c    $                    + F6*at06*W5 + F7*at07*W6
c    $                    + F8*at08*W7 + F9*at09*W8
c    $                    + F10*at10*W9+F11*at11*W10)
c           dF2dW = (bt01 + F2*bt02*W1 + F3*bt03*W2
c    $                    + F4*bt04*W3 + F5*bt05*W4
c    $                    + F6*bt06*W5 + F7*bt07*W6
c    $                    + F8*bt08*W7 + F9*bt09*W8
c    $                    + F10*Bt10*W9+F11*Bt11*W10)
c           dF3dW = (ct01 + F2*ct02*W1 + F3*ct03*W2
c    $                    + F4*ct04*W3 + F5*ct05*W4
c    $                    + F6*ct06*W5 + F7*ct07*W6
c    $                    + F8*ct08*W7 + F9*ct09*W8
c    $                    + F10*ct10*W9+F11*ct11*W10)
c           dF4dW = (dt01 + F2*dt02*W1 + F3*dt03*W2
c    $                    + F4*dt04*W3 + F5*dt05*W4
c    $                    + F6*dt06*W5 + F7*dt07*W6
c    $                    + F8*dt08*W7 + F9*dt09*W8
c    $                    + F10*dt10*W9+F11*dt11*W10)

c           dWdT = F2/((F1 + Tsig)**F2)
c           dTdR = ((F6*PI*PI)**F2o3)*(rhoo**F2o3)/tauN
c           dTdTau = -TauUEG/tauN**F2
C
c          If(UseLC) then
c            dElSRdR = PDUM
c            dElLRdR = Ax*F4o3*Rho13-PDUM
c          else
c            dElSRdR=Ax*F4o3*Rho13
c            dElLRdR=F0
c          endIf  
c          dGGA1dR = dElSRdR*fx1 + ElSR*dfx1dRho
c          dGGA2dR = dElSRdR*fx2 + ElSR*dfx2dRho 
c          dGGA3dR = dElLRdR*fx1 + ElLR*dfx1dRho
c          dGGA4dR = dElLRdR*fx2 + ElLR*dfx2dRho 
c
c          dF1dR = dF1dW*dWdT*dTdR
c          dF1dTau=dF1dW*dWdT*dTdTau
c          dF2dR = dF2dW*dWdT*dTdR
c          dF2dTau=dF2dW*dWdT*dTdTau
c          dF3dR = dF3dW*dWdT*dTdR
c          dF3dTau=dF3dW*dWdT*dTdTau
c          dF4dR = dF4dW*dWdT*dTdR
c          dF4dTau=dF4dW*dWdT*dTdTau
c
c          dGGA1dG = ElSR*dfx1dG
c          dGGA2dG = ElSR*dfx2dG
c          dGGA3dG = ElLR*dfx1dG
c          dGGA4dG = ElLR*dfx2dG
c
c          Amat(n,D1_RB) = Amat(n,D1_RB)
c    $                   + (dGGA1dR*Fsig1 + GGA1*dF1dR
c    $                     +dGGA2dR*Fsig2 + GGA2*dF2dR
c    $                     +dGGA3dR*Fsig3 + GGA3*dF3dR
c    $                     +dGGA4dR*Fsig4 + GGA4*dF4dR)*wght
c          Cmat(n,D1_GBB)= Cmat(n,D1_GBB)
c    $                   + (dGGA1dG*Fsig1 + dGGA2dG*Fsig2
c    $                     +dGGA3dG*Fsig3 + dGGA4dG*Fsig4)*wght
c          Mmat(n,D1_TB) = Mmat(n,D1_TB)
c    $                   + (GGA1*dF1dTau  + GGA2*dF2dTau
c    $                     +GGA3*dF3dTau  + GGA4*dF4dTau)*wght
c
20      continue
      endif
      return
      end
c
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
      SUBROUTINE nwxc_LRCLSDA_p(Emu,Rho,F)
#else
      SUBROUTINE nwxc_LRCLSDA(Emu,Rho,F)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      SUBROUTINE nwxc_LRCLSDA_d2(Emu,Rho,F)
#else
      SUBROUTINE nwxc_LRCLSDA_d3(Emu,Rho,F)
#endif
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
#include "nwad.fh"
c
c     IMPLICIT REAL*8 (a-h,o-z)
      implicit none
      double precision Emu
      type(nwad_dble)::Rho,F
      type(nwad_dble)::tmu,tmu2,tmu3
      type(nwad_dble)::erfv,rho13,rho43
      type(nwad_dble)::fsr,w
      double precision f1o2,f1o3,f1o4,f4o3,f8o3,pi12
      double precision AX, Cmu
      double precision D1F, Pi
      double precision F1, F2, F3, F4, F5, F6, F7, F8, F9
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
      W    = Exp(-F1o4/tmu2)
      ERFV = Erf( F1o2/tmu)
c     dtmudR = -F1o3*tmu / Rho
c
      Fsr = F1-F4o3*tmu*(-F6*tmu+F8*tmu3+W*
     $        (F4*tmu-F8*tmu3)+F2*PI12*ERFV)
c     dFsrdtmu = F8o3*(F2*tmu*(F3-F8*tmu2+W*
c    $          (-F1+F8*tmu2))-PI12*ERFV)
c
      F = Ax*Rho43*Fsr
c     D1F = Ax*F4o3*Rho13*Fsr + Ax*Rho43*(dFsrdtmu*dtmudR)
c
      RETURN
      END
#if !defined(NWAD_PRINT)
#define NWAD_PRINT
c
c Recompile for 2nd order derivatives
c
#include "nwxc_x_m11.F"
#endif
#if !defined(SECOND_DERIV)
#define SECOND_DERIV
c
c Recompile for 2nd order derivatives
c
#include "nwxc_x_m11.F"
#endif
#if !defined(THIRD_DERIV)
#define THIRD_DERIV
c
c Recompile for 3rd order derivatives
c
#include "nwxc_x_m11.F"
#endif
#undef NWAD_PRINT
c
C>
C> @}
