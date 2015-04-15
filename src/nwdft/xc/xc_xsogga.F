      Subroutine xc_xsogga(tol_rho, fac,lfac,nlfac, rho, delrho, 
     &                     Amat, Cmat, nq, ipol, Ex, 
     &                     qwght, ldew, func, ijzy)
c   
c$Id$
c
c
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
c      ijzy = 1 - SOGGA functional (a) - it requres PBE correlation    c
c      ijzy = 2 - SOGGA11 functional (b)                               c
c      ijzy = 3 - SOGGA11-X functional (c)                             c
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
      logical lfac, nlfac,ldew
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
      double precision Amat(nq,ipol), Cmat(nq,*)
      double precision tol_rho, pi
c
c     Intermediate derivative results, etc.
c
      integer n, ijzy
c
      double precision CxA0,CxA1,CxA2,CxA3,CxA4,CxA5
      double precision CxB0,CxB1,CxB2,CxB3,CxB4,CxB5
      double precision rho43, rho13, rhoo
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
      double precision DYDG, DYDR, DTOL, ELOC
      double precision FA0, FA1, FA2, FA3, FA4, FA5
      double precision FB0, FB1, FB2, FB3, FB4, FB5
      double precision FEXP, FFRAC, FGGAX, MU, PON, S, X, Y       
      double precision Gam12, Gam      
c      
      double precision f1,f2,f3,f4,f5,f8
      double precision F1o3,F4o3,F48
      parameter( F1=1.0D+00,  F2=2.0D+00,  F3=3.0D+00,  
     $           F4=4.0D+00,  F5=5.0D+00,  F8=8.0D+00,
     $           F48=48.0D+00)
c
        pi=acos(-1d0)      
c
      if (ijzy.eq.1) then
c SOGGA11
       CxA0 = 0.5d0
       CxA1 = 0.276d0
       CxA2 = 0d0
       CxA3 = 0d0
       CxA4 = 0d0
       CxA5 = 0d0
       CxB0 = 0.5d0
       CxB1 = 0.276d0
       CxB2 = 0d0
       CxB3 = 0d0
       CxB4 = 0d0
       CxB5 = 0d0
      elseif (ijzy.eq.2) then
c SOGGA11
       CxA0 =  0.50000d0
       CxA1 = -2.95535d0
       CxA2 =  15.7974d0
       CxA3 = -91.1804d0
       CxA4 =  96.2030d0
       CxA5 =  0.18683d0
       CxB0 =  0.50000d0
       CxB1 =  3.50743d0
       CxB2 = -12.9523d0
       CxB3 =  49.7870d0
       CxB4 = -33.2545d0
       CxB5 = -11.1396d0
      elseif (ijzy.eq.3) then
c SOGGA11-X
       CxA0 =  2.99250d-01
       CxA1 =  3.21638d+00
       CxA2 = -3.55605d+00
       CxA3 =  7.65852d+00
       CxA4 = -1.12830d+01
       CxA5 =  5.25813d+00
       CxB0 =  2.99250d-01
       CxB1 = -2.88595d+00
       CxB2 =  3.23617d+00
       CxB3 = -2.45393d+00
       CxB4 = -3.75495d+00
       CxB5 =  3.96613d+00
      else
       call errquit("xc_xsogga: illegal value of ijzy",ijzy,UERR)
      endif
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
            if (rho(n,1).lt.DTol) goto 10
            rhoo = rho(n,1)/F2
            rho43 = rhoo**F4o3  
            rho13 = rho43/rhoo
            Gam =(delrho(n,1,1)*delrho(n,1,1) +
     &              delrho(n,2,1)*delrho(n,2,1) +
     &              delrho(n,3,1)*delrho(n,3,1))/F4
            Gam12 = dsqrt(Gam)
            if(gam12.lt.dtol) goto 10
c
          Eloc = Ax*Rho43
          x = Gam12/Rho43
          s = As*x
          y = s*s
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
     $            fb0+fb1+fb2+fb3+fb4+fb5
C
C     1st derivatives.
C

          dElocdR=Ax*F4o3*Rho13
          dydR = -(F8/F3)*y/Rhoo
          dydG   = y/Gam
          dFfracdPON = F1/((F1+PON)**F2)
          dFexpdPON  = exp(-PON)
          dfa1dy = CxA1 *mu*dFfracdPON
          dfa2dy = CxA2 *mu*F2*Ffrac*dFfracdPON
          dfa3dy = CxA3 *mu*(F3 *Ffrac**F2)*dFfracdPON
          dfa4dy = CxA4 *mu*(F4 *Ffrac**F3)*dFfracdPON
          dfa5dy = CxA5 *mu*(F5 *Ffrac**F4)*dFfracdPON
          dfa1dR = dfa1dy *dydR
          dfa2dR = dfa2dy *dydR
          dfa3dR = dfa3dy *dydR
          dfa4dR = dfa4dy *dydR
          dfa5dR = dfa5dy *dydR
          dfa1dG = dfa1dy *dydG
          dfa2dG = dfa2dy *dydG
          dfa3dG = dfa3dy *dydG
          dfa4dG = dfa4dy *dydG
          dfa5dG = dfa5dy *dydG
          dfb1dy = CxB1 *mu*dFexpdPON
          dfb2dy = CxB2 *mu*F2*Fexp*dFexpdPON
          dfb3dy = CxB3 *mu*(F3 *Fexp**F2)*dFexpdPON
          dfb4dy = CxB4 *mu*(F4 *Fexp**F3)*dFexpdPON
          dfb5dy = CxB5 *mu*(F5 *Fexp**F4)*dFexpdPON
          dfb1dR = dfb1dy *dydR
          dfb2dR = dfb2dy *dydR
          dfb3dR = dfb3dy *dydR
          dfb4dR = dfb4dy *dydR
          dfb5dR = dfb5dy *dydR
          dfb1dG = dfb1dy *dydG
          dfb2dG = dfb2dy *dydG
          dfb3dG = dfb3dy *dydG
          dfb4dG = dfb4dy *dydG
          dfb5dG = dfb5dy *dydG
c
          dFggaxdR = dfa1dR+dfa2dR+dfa3dR+dfa4dR+dfa5dR +
     $               dfb1dR+dfb2dR+dfb3dR+dfb4dR+dfb5dR
c                     
          dFggaxdG = dfa1dG+dfa2dG+dfa3dG+dfa4dG+dfa5dG +
     $               dfb1dG+dfb2dG+dfb3dG+dfb4dG+dfb5dG
c
           Ex = Ex +F2*(Eloc*Fggax)*qwght(n)
           if(ldew) func(n)=func(n)+F2*(Eloc*Fggax)
           Amat(n,1) = Amat(n,1)   +dElocdR*Fggax + Eloc*dFggaxdR        
           Cmat(n,1)=  Cmat(n,1)  + Eloc*dFggaxdG
           
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
            Gam =(delrho(n,1,1)*delrho(n,1,1) +
     &              delrho(n,2,1)*delrho(n,2,1) +
     &              delrho(n,3,1)*delrho(n,3,1))
            Gam12 = dsqrt(Gam)
         if(gam12.lt.dtol) goto 25
c
          Eloc = Ax*Rho43
          x = Gam12/Rho43
          s = As*x
          y = s*s
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
     $            fb0+fb1+fb2+fb3+fb4+fb5
C
C     1st derivatives.
C
           dElocdR=Ax*F4o3*Rho13
           dydR = -(F8/F3)*y/Rhoo
           dydG   = y/Gam
           dFfracdPON = F1/((F1+PON)**F2)
           dFexpdPON  = exp(-PON)
           dfa1dy = CxA1 *mu*dFfracdPON
           dfa2dy = CxA2 *mu*F2*Ffrac*dFfracdPON
           dfa3dy = CxA3 *mu*(F3 *Ffrac**F2)*dFfracdPON
           dfa4dy = CxA4 *mu*(F4 *Ffrac**F3)*dFfracdPON
           dfa5dy = CxA5 *mu*(F5 *Ffrac**F4)*dFfracdPON
           dfa1dR = dfa1dy *dydR
           dfa2dR = dfa2dy *dydR
           dfa3dR = dfa3dy *dydR
           dfa4dR = dfa4dy *dydR
           dfa5dR = dfa5dy *dydR
           dfa1dG = dfa1dy *dydG
           dfa2dG = dfa2dy *dydG
           dfa3dG = dfa3dy *dydG
           dfa4dG = dfa4dy *dydG
           dfa5dG = dfa5dy *dydG
           dfb1dy = CxB1 *mu*dFexpdPON
           dfb2dy = CxB2 *mu*F2*Fexp*dFexpdPON
           dfb3dy = CxB3 *mu*(F3 *Fexp**F2)*dFexpdPON
           dfb4dy = CxB4 *mu*(F4 *Fexp**F3)*dFexpdPON
           dfb5dy = CxB5 *mu*(F5 *Fexp**F4)*dFexpdPON
           dfb1dR = dfb1dy *dydR
           dfb2dR = dfb2dy *dydR
           dfb3dR = dfb3dy *dydR
           dfb4dR = dfb4dy *dydR
           dfb5dR = dfb5dy *dydR
           dfb1dG = dfb1dy *dydG
           dfb2dG = dfb2dy *dydG
           dfb3dG = dfb3dy *dydG
           dfb4dG = dfb4dy *dydG
           dfb5dG = dfb5dy *dydG
c
           dFggaxdR = dfa1dR+dfa2dR+dfa3dR+dfa4dR+dfa5dR +
     $                dfb1dR+dfb2dR+dfb3dR+dfb4dR+dfb5dR
c                      
           dFggaxdG = dfa1dG+dfa2dG+dfa3dG+dfa4dG+dfa5dG +
     $                dfb1dG+dfb2dG+dfb3dG+dfb4dG+dfb5dG
c
           Ex = Ex + (Eloc*Fggax)*qwght(n)
           if(ldew) func(n)=func(n)+ Eloc*Fggax
           Amat(n,1) = Amat(n,1)   + dElocdR*Fggax + Eloc*dFggaxdR
           Cmat(n,1)=  Cmat(n,1)   + Eloc*dFggaxdG
c
25         continue
c
c     Beta               BETA           BETA
c
            if (rho(n,3).lt.DTol) goto 20
            rhoo  = rho(n,3)
            rho43 = rhoo**F4o3  
            rho13 = rho43/rhoo
c
            Gam =(delrho(n,1,2)*delrho(n,1,2) +
     &            delrho(n,2,2)*delrho(n,2,2) +
     &            delrho(n,3,2)*delrho(n,3,2))
            Gam12 = dsqrt(Gam)
         if(gam12.lt.dtol) goto 20
c
          Eloc = Ax*Rho43
          x = Gam12/Rho43
          s = As*x
          y = s*s
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
     $            fb0+fb1+fb2+fb3+fb4+fb5
C
C     1st derivatives.
C

          dElocdR=Ax*F4o3*Rho13
          dydR = -(F8/F3)*y/Rhoo
          dydG   = y/Gam
          dFfracdPON = F1/((F1+PON)**F2)
          dFexpdPON  = exp(-PON)
          dfa1dy = CxA1 *mu*dFfracdPON
          dfa2dy = CxA2 *mu*F2*Ffrac*dFfracdPON
          dfa3dy = CxA3 *mu*(F3 *Ffrac**F2)*dFfracdPON
          dfa4dy = CxA4 *mu*(F4 *Ffrac**F3)*dFfracdPON
          dfa5dy = CxA5 *mu*(F5 *Ffrac**F4)*dFfracdPON
          dfa1dR = dfa1dy *dydR
          dfa2dR = dfa2dy *dydR
          dfa3dR = dfa3dy *dydR
          dfa4dR = dfa4dy *dydR
          dfa5dR = dfa5dy *dydR
          dfa1dG = dfa1dy *dydG
          dfa2dG = dfa2dy *dydG
          dfa3dG = dfa3dy *dydG
          dfa4dG = dfa4dy *dydG
          dfa5dG = dfa5dy *dydG
          dfb1dy = CxB1 *mu*dFexpdPON
          dfb2dy = CxB2 *mu*F2*Fexp*dFexpdPON
          dfb3dy = CxB3 *mu*(F3 *Fexp**F2)*dFexpdPON
          dfb4dy = CxB4 *mu*(F4 *Fexp**F3)*dFexpdPON
          dfb5dy = CxB5 *mu*(F5 *Fexp**F4)*dFexpdPON
          dfb1dR = dfb1dy *dydR
          dfb2dR = dfb2dy *dydR
          dfb3dR = dfb3dy *dydR
          dfb4dR = dfb4dy *dydR
          dfb5dR = dfb5dy *dydR
          dfb1dG = dfb1dy *dydG
          dfb2dG = dfb2dy *dydG
          dfb3dG = dfb3dy *dydG
          dfb4dG = dfb4dy *dydG
          dfb5dG = dfb5dy *dydG
c
          dFggaxdR = dfa1dR+dfa2dR+dfa3dR+dfa4dR+dfa5dR +
     $               dfb1dR+dfb2dR+dfb3dR+dfb4dR+dfb5dR
c                     
          dFggaxdG = dfa1dG+dfa2dG+dfa3dG+dfa4dG+dfa5dG +
     $               dfb1dG+dfb2dG+dfb3dG+dfb4dG+dfb5dG
c
          Ex = Ex + (Eloc*Fggax)*qwght(n)
          if(ldew) func(n)=func(n)+ Eloc*Fggax


           Amat(n,2) = Amat(n,2)  +dElocdR*Fggax + Eloc*dFggaxdR
                      
           Cmat(n,3)=  Cmat(n,3)  + Eloc*dFggaxdG
    
c
20      continue
      endif
      return
      end
c
      Subroutine xc_xsogga_d2()
      call errquit(' not coded ',0,0)
      return
      end

