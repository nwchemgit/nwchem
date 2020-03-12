#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if !defined(NWAD_PRINT)
C> \ingroup nwxc
C> @{
C>
C> \file nwxc_x_sogga.F
C> The SOGGA, SOGGA11 and SOGGA-X exchange functionals
C>
C> @}
#endif
#endif
C>
C> \ingroup nwxc_priv
C> @{
C>
C> \brief Evaluate the SOGGA, SOGGA11 and SOGGA11-X exchange functionals
C>
C> The SOGGA, SOGGA11 and SOGGA11-X functionals are GGA functionals that
C> are formulated and optimized such that the coefficient of the
C> gradient correction term matches that of the appropriate expansion ofC> the electronic energy [1,2,3].
C>
C> ### References ###
C>
C> [1] Y. Zhao, D.G. Truhlar,
C>     "Construction of a generalized gradient approximation by 
C>     restoring the density-gradient expansion and enforcing a tight
C>     Lieb-Oxford bound", J. Chem. Phys. <b>128</b> (2008) 184109,
C>     DOI:
C>     <a href="https://doi.org/10.1063/1.2912068">
C>     10.1063/1.2912068</a>.
C>
C> [2] R. Peverati, Y. Zhao, D.G. Truhlar, 
C>     "Generalized gradient approximation that recovers the
C>     second-order density-gradient expansion with optimized
C>     across-the-board performance", J. Phys. Chem. Lett. <b>2</b>
C>     (2011) 1991-1997, DOI:
C>     <a href="https://doi.org/10.1021/jz200616w">
C>     10.1021/jz200616w</a>.
C>
C> [3] R. Peverati, D.G. Truhlar, 
C>     "Communication: A global hybrid generalized gradient
C>     approximation to the exchange-correlation functional that
C>     satisfies the second-order density-gradient constraint and has
C>     broad applicability in chemistry", J. Chem. Phys. <b>135</b>
C>     (2011) 191102, DOI:
C>     <a href="https://doi.org/10.1063/1.3663871">
C>     10.1063/1.3663871</a>.
C>
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
      Subroutine nwxc_x_sogga_p(param, tol_rho, ipol, nq, wght, rho,
     &                          rgamma, ffunc)
#else
      Subroutine nwxc_x_sogga(param, tol_rho, ipol, nq, wght, rho,
     &                        rgamma, ffunc)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      Subroutine nwxc_x_sogga_d2(param, tol_rho, ipol, nq, wght, rho,
     &                           rgamma, ffunc)
#else
      Subroutine nwxc_x_sogga_d3(param, tol_rho, ipol, nq, wght, rho,
     &                           rgamma, ffunc)
#endif
c   
c$Id$
c
c**********************************************************************c
c                                                                      c
c  SOGGA11X evaluates the exchange part of the SOGGA, SOGGA11          c
c  and SOGGA11-X functionals on the grid.                              c
c                                                                      c
c     a) Zhao and Truhlar, J.Chem.Phys., 128, 184109 (2008)            c
c     b) Peverati, Zhao and Truhlar, J.Phys.Chem.Lett, 2, 1991 (2011)  c
c     c) Peverati and Truhlar, J.Chem.Phys, 135, 191102 (2011)         c
c                                                                      c
c      ijzy = 1 - SOGGA functional (a) - it requires PBE correlation   c
c      ijzy = 2 - SOGGA11 functional (b)                               c
c      ijzy = 3 - SOGGA11-X functional (c)                             c
c                                                                      c
c Coded by Roberto Peverati (12/11)                                    c
c                                                                      c
c**********************************************************************c
c
#include "nwad.fh"
c
      implicit none
#include "nwxc_param.fh"
c
c     Input and other parameters
c
#if defined(NWAD_PRINT)
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      type(nwad_dble)::param(*)!< [Input] Parameters of the functional,
      type(nwad_dble)::CxA0,CxA1,CxA2,CxA3,CxA4,CxA5
      type(nwad_dble)::CxB0,CxB1,CxB2,CxB3,CxB4,CxB5
      type(nwad_dble)::FA0
      type(nwad_dble)::FB0
#else
      double precision param(*)!< [Input] Parameters of the functional,
      double precision CxA0,CxA1,CxA2,CxA3,CxA4,CxA5
      double precision CxB0,CxB1,CxB2,CxB3,CxB4,CxB5
      double precision FA0
      double precision FB0
#endif
#else
      double precision param(*)!< [Input] Parameters of the functional,
                               !< see Table 1 of [1] and Table 1 of [2].
                               !< - param(1): \f$ a_0 \f$
                               !< - param(2): \f$ a_1 \f$
                               !< - param(3): \f$ a_2 \f$
                               !< - param(4): \f$ a_3 \f$
                               !< - param(5): \f$ a_4 \f$
                               !< - param(6): \f$ a_5 \f$
                               !< - param(7): \f$ b_0 \f$
                               !< - param(8): \f$ b_1 \f$
                               !< - param(9): \f$ b_2 \f$
                               !< - param(10): \f$ b_3 \f$
                               !< - param(11): \f$ b_4 \f$
                               !< - param(12): \f$ b_5 \f$
      double precision CxA0,CxA1,CxA2,CxA3,CxA4,CxA5
      double precision CxB0,CxB1,CxB2,CxB3,CxB4,CxB5
      double precision FA0
      double precision FB0
#endif
      double precision tol_rho !< [Input] The lower limit on the density
      integer ipol             !< [Input] The number of spin channels
      integer nq               !< [Input] The number of points
      double precision wght    !< [Input] The weight of the functional
c
c     Charge Density 
c
      type(nwad_dble)::rho(nq,*)    !< [Input] The density
c
c     Charge Density Gradient
c
      type(nwad_dble)::rgamma(nq,*) !< [Input] The norm of the density gradients
c
c     Sampling Matrices for the XC Potential
c
      type(nwad_dble)::ffunc(nq)    !< [Output] The value of the functional
c     double precision Amat(nq,*)   !< [Output] The derivative wrt rho
c     double precision Cmat(nq,*)   !< [Output] The derivative wrt rgamma
c      
      double precision pi
c
c     Intermediate derivative results, etc.
c
      integer n
c
      type(nwad_dble)::rho43, rho13, rhoo
c
      double precision AS, ASO, AX, DELOCDR
      double precision DFA1DG, DFA1DR, DFA1DY
      double precision DFA2DG, DFA2DR, DFA2DY
      double precision DFA3DG, DFA3DR, DFA3DY
      double precision DFA4DG, DFA4DR, DFA4DY
      double precision DFA5DG, DFA5DR, DFA5DY
      double precision DFB1DG, DFB1DR, DFB1DY
      double precision DFB2DG, DFB2DR, DFB2DY
      double precision DFB3DG, DFB3DR, DFB3DY
      double precision DFB4DG, DFB4DR, DFB4DY
      double precision DFB5DG, DFB5DR, DFB5DY
      double precision DFEXPDPON, DFFRACDPON, DFGGAXDG, DFGGAXDR
      double precision DYDG, DYDR, DTOL
      double precision MU
      type(nwad_dble)::FA1, FA2, FA3, FA4, FA5
      type(nwad_dble)::FB1, FB2, FB3, FB4, FB5
c     type(nwad_dble)::Gam12, Gam, X, FEXP, FFRAC, FGGAX, S, Y, ELOC
      type(nwad_dble)::Gam, X, FEXP, FFRAC, FGGAX, Y, ELOC
      type(nwad_dble)::PON
c      
      double precision f1,f2,f3,f4,f5,f8
      double precision F1o3,F4o3,F48
      parameter( F1=1.0D+00,  F2=2.0D+00,  F3=3.0D+00,  
     $           F4=4.0D+00,  F5=5.0D+00,  F8=8.0D+00,
     $           F48=48.0D+00)
c
      pi=acos(-1d0)      
c
c     if (ijzy.eq.1) then
c SOGGA
c      CxA0 = 0.5d0
c      CxA1 = 0.276d0
c      CxA2 = 0d0
c      CxA3 = 0d0
c      CxA4 = 0d0
c      CxA5 = 0d0
c      CxB0 = 0.5d0
c      CxB1 = 0.276d0
c      CxB2 = 0d0
c      CxB3 = 0d0
c      CxB4 = 0d0
c      CxB5 = 0d0
c     elseif (ijzy.eq.2) then
c SOGGA11
c      CxA0 =  0.50000d0
c      CxA1 = -2.95535d0
c      CxA2 =  15.7974d0
c      CxA3 = -91.1804d0
c      CxA4 =  96.2030d0
c      CxA5 =  0.18683d0
c      CxB0 =  0.50000d0
c      CxB1 =  3.50743d0
c      CxB2 = -12.9523d0
c      CxB3 =  49.7870d0
c      CxB4 = -33.2545d0
c      CxB5 = -11.1396d0
c     elseif (ijzy.eq.3) then
c SOGGA11-X
c      CxA0 =  2.99250d-01
c      CxA1 =  3.21638d+00
c      CxA2 = -3.55605d+00
c      CxA3 =  7.65852d+00
c      CxA4 = -1.12830d+01
c      CxA5 =  5.25813d+00
c      CxB0 =  2.99250d-01
c      CxB1 = -2.88595d+00
c      CxB2 =  3.23617d+00
c      CxB3 = -2.45393d+00
c      CxB4 = -3.75495d+00
c      CxB5 =  3.96613d+00
c     endif
      CxA0 = param(1)
      CxA1 = param(2)
      CxA2 = param(3)
      CxA3 = param(4)
      CxA4 = param(5)
      CxA5 = param(6)
      CxB0 = param(7)
      CxB1 = param(8)
      CxB2 = param(9)
      CxB3 = param(10)
      CxB4 = param(11)
      CxB5 = param(12)
c
      DTol = tol_rho
      F1o3 = F1/F3 
      F4o3 = F4/F3
      Pi   = ACos(-F1)
      AsO  = (F48*PI*PI)**F1o3
      As   = F1/AsO
      Ax   = -(F3/F2) * (F4o3*Pi)**(-F1o3)
      mu = 0.2236536053d0
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
            rho13 = rho43/rhoo
            Gam = rgamma(n,G_TT)/F4
c           Gam =(delrho(n,1,1)*delrho(n,1,1) +
c    &              delrho(n,2,1)*delrho(n,2,1) +
c    &              delrho(n,3,1)*delrho(n,3,1))/F4
c           Gam12 = sqrt(Gam)
c           if(gam12.lt.dtol) goto 10
c
            Eloc = Ax*Rho43
c           x = Gam12/Rho43
c           s = As*x
c           y = s*s
            y = As*As*Gam/(Rho43*Rho43)
            PON = mu*y
            Ffrac = F1-F1/(F1+PON)
            Fexp  = F1-exp(-PON)
            fa0 = CxA0
            fa1 = CxA1 *Ffrac
            fa2 = CxA2 *Ffrac**F2
            fa3 = CxA3 *Ffrac**F3
            fa4 = CxA4 *Ffrac**F4
            fa5 = CxA5 *Ffrac**F5
            fb0 = CxB0
            fb1 = CxB1 *Fexp
            fb2 = CxB2 *Fexp**F2
            fb3 = CxB3 *Fexp**F3
            fb4 = CxB4 *Fexp**F4
            fb5 = CxB5 *Fexp**F5
c
            Fggax = fa0+fa1+fa2+fa3+fa4+fa5 +
     $              fb0+fb1+fb2+fb3+fb4+fb5
C
C     1st derivatives.
C

c           dElocdR=Ax*F4o3*Rho13
c           dydR = -(F8/F3)*y/Rhoo
c           dydG   = y/Gam
c           dFfracdPON = F1/((F1+PON)**F2)
c           dFexpdPON  = exp(-PON)
c           dfa1dy = CxA1 *mu*dFfracdPON
c           dfa2dy = CxA2 *mu*F2*Ffrac*dFfracdPON
c           dfa3dy = CxA3 *mu*(F3 *Ffrac**F2)*dFfracdPON
c           dfa4dy = CxA4 *mu*(F4 *Ffrac**F3)*dFfracdPON
c           dfa5dy = CxA5 *mu*(F5 *Ffrac**F4)*dFfracdPON
c           dfa1dR = dfa1dy *dydR
c           dfa2dR = dfa2dy *dydR
c           dfa3dR = dfa3dy *dydR
c           dfa4dR = dfa4dy *dydR
c           dfa5dR = dfa5dy *dydR
c           dfa1dG = dfa1dy *dydG
c           dfa2dG = dfa2dy *dydG
c           dfa3dG = dfa3dy *dydG
c           dfa4dG = dfa4dy *dydG
c           dfa5dG = dfa5dy *dydG
c           dfb1dy = CxB1 *mu*dFexpdPON
c           dfb2dy = CxB2 *mu*F2*Fexp*dFexpdPON
c           dfb3dy = CxB3 *mu*(F3 *Fexp**F2)*dFexpdPON
c           dfb4dy = CxB4 *mu*(F4 *Fexp**F3)*dFexpdPON
c           dfb5dy = CxB5 *mu*(F5 *Fexp**F4)*dFexpdPON
c           dfb1dR = dfb1dy *dydR
c           dfb2dR = dfb2dy *dydR
c           dfb3dR = dfb3dy *dydR
c           dfb4dR = dfb4dy *dydR
c           dfb5dR = dfb5dy *dydR
c           dfb1dG = dfb1dy *dydG
c           dfb2dG = dfb2dy *dydG
c           dfb3dG = dfb3dy *dydG
c           dfb4dG = dfb4dy *dydG
c           dfb5dG = dfb5dy *dydG
c
c           dFggaxdR = dfa1dR+dfa2dR+dfa3dR+dfa4dR+dfa5dR +
c    $                 dfb1dR+dfb2dR+dfb3dR+dfb4dR+dfb5dR
c                     
c           dFggaxdG = dfa1dG+dfa2dG+dfa3dG+dfa4dG+dfa5dG +
c    $                 dfb1dG+dfb2dG+dfb3dG+dfb4dG+dfb5dG
c
c          Ex = Ex +F2*(Eloc*Fggax)*qwght(n)
           ffunc(n)=ffunc(n)+F2*(Eloc*Fggax)*wght
c          Amat(n,D1_RA)=Amat(n,D1_RA)
c    $                  +(dElocdR*Fggax+Eloc*dFggaxdR)*wght
c          Cmat(n,D1_GAA)=Cmat(n,D1_GAA)+Eloc*dFggaxdG*wght
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
           if (rho(n,R_A).lt.DTol) goto 25           
           rhoo  = rho(n,R_A)
           rho43 = rhoo**F4o3  
           rho13 = rho43/rhoo
           Gam = rgamma(n,G_AA)
c          Gam =(delrho(n,1,1)*delrho(n,1,1) +
c    &           delrho(n,2,1)*delrho(n,2,1) +
c    &           delrho(n,3,1)*delrho(n,3,1))
c          Gam12 = sqrt(Gam)
c          if(gam12.lt.dtol) goto 25
c
           Eloc = Ax*Rho43
c          x = Gam12/Rho43
c          s = As*x
c          y = s*s
           y = As*As*Gam/(Rho43*Rho43)
           PON = mu*y
           Ffrac = F1-F1/(F1+PON)
           Fexp  = F1-exp(-PON)
           fa0 = CxA0
           fa1 = CxA1 *Ffrac
           fa2 = CxA2 *Ffrac**F2
           fa3 = CxA3 *Ffrac**F3
           fa4 = CxA4 *Ffrac**F4
           fa5 = CxA5 *Ffrac**F5
           fb0 = CxB0
           fb1 = CxB1 *Fexp
           fb2 = CxB2 *Fexp**F2
           fb3 = CxB3 *Fexp**F3
           fb4 = CxB4 *Fexp**F4
           fb5 = CxB5 *Fexp**F5
c          
           Fggax = fa0+fa1+fa2+fa3+fa4+fa5 +
     $             fb0+fb1+fb2+fb3+fb4+fb5
C
C     1st derivatives.
C
c          dElocdR=Ax*F4o3*Rho13
c          dydR = -(F8/F3)*y/Rhoo
c          dydG   = y/Gam
c          dFfracdPON = F1/((F1+PON)**F2)
c          dFexpdPON  = exp(-PON)
c          dfa1dy = CxA1 *mu*dFfracdPON
c          dfa2dy = CxA2 *mu*F2*Ffrac*dFfracdPON
c          dfa3dy = CxA3 *mu*(F3 *Ffrac**F2)*dFfracdPON
c          dfa4dy = CxA4 *mu*(F4 *Ffrac**F3)*dFfracdPON
c          dfa5dy = CxA5 *mu*(F5 *Ffrac**F4)*dFfracdPON
c          dfa1dR = dfa1dy *dydR
c          dfa2dR = dfa2dy *dydR
c          dfa3dR = dfa3dy *dydR
c          dfa4dR = dfa4dy *dydR
c          dfa5dR = dfa5dy *dydR
c          dfa1dG = dfa1dy *dydG
c          dfa2dG = dfa2dy *dydG
c          dfa3dG = dfa3dy *dydG
c          dfa4dG = dfa4dy *dydG
c          dfa5dG = dfa5dy *dydG
c          dfb1dy = CxB1 *mu*dFexpdPON
c          dfb2dy = CxB2 *mu*F2*Fexp*dFexpdPON
c          dfb3dy = CxB3 *mu*(F3 *Fexp**F2)*dFexpdPON
c          dfb4dy = CxB4 *mu*(F4 *Fexp**F3)*dFexpdPON
c          dfb5dy = CxB5 *mu*(F5 *Fexp**F4)*dFexpdPON
c          dfb1dR = dfb1dy *dydR
c          dfb2dR = dfb2dy *dydR
c          dfb3dR = dfb3dy *dydR
c          dfb4dR = dfb4dy *dydR
c          dfb5dR = dfb5dy *dydR
c          dfb1dG = dfb1dy *dydG
c          dfb2dG = dfb2dy *dydG
c          dfb3dG = dfb3dy *dydG
c          dfb4dG = dfb4dy *dydG
c          dfb5dG = dfb5dy *dydG
c
c          dFggaxdR = dfa1dR+dfa2dR+dfa3dR+dfa4dR+dfa5dR +
c    $                dfb1dR+dfb2dR+dfb3dR+dfb4dR+dfb5dR
c                      
c          dFggaxdG = dfa1dG+dfa2dG+dfa3dG+dfa4dG+dfa5dG +
c    $                dfb1dG+dfb2dG+dfb3dG+dfb4dG+dfb5dG
c
c          Ex = Ex + (Eloc*Fggax)*qwght(n)
           ffunc(n)=ffunc(n)+ Eloc*Fggax*wght
c          Amat(n,D1_RA) = Amat(n,D1_RA) + dElocdR*Fggax*wght 
c    $                                   + Eloc*dFggaxdR*wght
c          Cmat(n,D1_GAA)=  Cmat(n,D1_GAA)   + Eloc*dFggaxdG*wght
c
25         continue
c
c     Beta               BETA           BETA
c
           if (rho(n,R_B).lt.DTol) goto 20
           rhoo  = rho(n,R_B)
           rho43 = rhoo**F4o3  
           rho13 = rho43/rhoo
c
           Gam = rgamma(n,G_BB)
c          Gam =(delrho(n,1,2)*delrho(n,1,2) +
c    &           delrho(n,2,2)*delrho(n,2,2) +
c    &           delrho(n,3,2)*delrho(n,3,2))
c          Gam12 = sqrt(Gam)
c          if(gam12.lt.dtol) goto 20
c
           Eloc = Ax*Rho43
c          x = Gam12/Rho43
c          s = As*x
c          y = s*s
           y = As*As*Gam/(Rho43*Rho43)
           PON = mu*y
           Ffrac = F1-F1/(F1+PON)
           Fexp  = F1-exp(-PON)
           fa0 = CxA0
           fa1 = CxA1 *Ffrac
           fa2 = CxA2 *Ffrac**F2
           fa3 = CxA3 *Ffrac**F3
           fa4 = CxA4 *Ffrac**F4
           fa5 = CxA5 *Ffrac**F5
           fb0 = CxB0
           fb1 = CxB1 *Fexp
           fb2 = CxB2 *Fexp**F2
           fb3 = CxB3 *Fexp**F3
           fb4 = CxB4 *Fexp**F4
           fb5 = CxB5 *Fexp**F5
c
           Fggax = fa0+fa1+fa2+fa3+fa4+fa5 +
     $             fb0+fb1+fb2+fb3+fb4+fb5
C
C     1st derivatives.
C

c          dElocdR=Ax*F4o3*Rho13
c          dydR = -(F8/F3)*y/Rhoo
c          dydG   = y/Gam
c          dFfracdPON = F1/((F1+PON)**F2)
c          dFexpdPON  = exp(-PON)
c          dfa1dy = CxA1 *mu*dFfracdPON
c          dfa2dy = CxA2 *mu*F2*Ffrac*dFfracdPON
c          dfa3dy = CxA3 *mu*(F3 *Ffrac**F2)*dFfracdPON
c          dfa4dy = CxA4 *mu*(F4 *Ffrac**F3)*dFfracdPON
c          dfa5dy = CxA5 *mu*(F5 *Ffrac**F4)*dFfracdPON
c          dfa1dR = dfa1dy *dydR
c          dfa2dR = dfa2dy *dydR
c          dfa3dR = dfa3dy *dydR
c          dfa4dR = dfa4dy *dydR
c          dfa5dR = dfa5dy *dydR
c          dfa1dG = dfa1dy *dydG
c          dfa2dG = dfa2dy *dydG
c          dfa3dG = dfa3dy *dydG
c          dfa4dG = dfa4dy *dydG
c          dfa5dG = dfa5dy *dydG
c          dfb1dy = CxB1 *mu*dFexpdPON
c          dfb2dy = CxB2 *mu*F2*Fexp*dFexpdPON
c          dfb3dy = CxB3 *mu*(F3 *Fexp**F2)*dFexpdPON
c          dfb4dy = CxB4 *mu*(F4 *Fexp**F3)*dFexpdPON
c          dfb5dy = CxB5 *mu*(F5 *Fexp**F4)*dFexpdPON
c          dfb1dR = dfb1dy *dydR
c          dfb2dR = dfb2dy *dydR
c          dfb3dR = dfb3dy *dydR
c          dfb4dR = dfb4dy *dydR
c          dfb5dR = dfb5dy *dydR
c          dfb1dG = dfb1dy *dydG
c          dfb2dG = dfb2dy *dydG
c          dfb3dG = dfb3dy *dydG
c          dfb4dG = dfb4dy *dydG
c          dfb5dG = dfb5dy *dydG
c
c          dFggaxdR = dfa1dR+dfa2dR+dfa3dR+dfa4dR+dfa5dR +
c    $                dfb1dR+dfb2dR+dfb3dR+dfb4dR+dfb5dR
c                     
c          dFggaxdG = dfa1dG+dfa2dG+dfa3dG+dfa4dG+dfa5dG +
c    $                dfb1dG+dfb2dG+dfb3dG+dfb4dG+dfb5dG
c
c          Ex = Ex + (Eloc*Fggax)*qwght(n)
           ffunc(n)=ffunc(n)+ Eloc*Fggax*wght

c          Amat(n,D1_RB) = Amat(n,D1_RB) + dElocdR*Fggax*wght
c    $                                   + Eloc*dFggaxdR*wght
c                     
c          Cmat(n,D1_GBB)=  Cmat(n,D1_GBB)  + Eloc*dFggaxdG*wght
c
20      continue
      endif
      return
      end
#ifndef NWAD_PRINT
#define NWAD_PRINT
c
c     compile again for 2nd derivatives
c
#include "nwxc_x_sogga.F"
#endif
#ifndef SECOND_DERIV
#define SECOND_DERIV
c
c     compile again for 2nd derivatives
c
#include "nwxc_x_sogga.F"
#endif
#ifndef THIRD_DERIV
#define THIRD_DERIV
c
c     compile again for 3rd derivatives
c
#include "nwxc_x_sogga.F"
#endif
#undef NWAD_PRINT
C>
C> @}
