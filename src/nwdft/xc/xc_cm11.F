      Subroutine xc_cm11(tol_rho, cfac, lcfac, nlcfac, rho, delrho, 
     &                     nq, ipol, Ec, qwght, ldew, func,
     &                     tau, Amat, Cmat, Mmat,ijzy)     

c
c$Id$
c
c
c**********************************************************************c
c                                                                      c
c  xc_cm11 evaluates the correlation part of the M08 and M11 suite of  c
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
#include "dft2drv.fh"
c
c     Input and other parameters
c
      integer ipol, nq
      double precision dummy(1)

      double precision cfac
      logical lcfac, nlcfac

      logical lfac, nlfac
      double precision fac
      double precision tol_rho
c
c     Constants in PBE functional
c
      double precision GAMMA, BETA, PI
      parameter (GAMMA = 0.03109069086965489503494086371273d0)
      parameter (BETA = 0.06672455060314922d0)
      parameter (PI = 3.1415926535897932385d0)
c
c     Threshold parameters
c
      double precision EXPTOL
      parameter (EXPTOL = 40.0d0)
c
c     Correlation energy
c
      double precision Ec
c
c     Charge Density 
c
      double precision rho(nq,ipol*(ipol+1)/2)
      double precision rho_t(3)
c
c     Charge Density Gradient
c
      double precision delrho(nq,3,ipol)
      double precision dsqgamma
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
      double precision Amat(nq,ipol), Cmat(nq,*), Mmat(nq,*)
c
c     Intermediate derivative results, etc.
c
      integer n, ijzy
      
      double precision at0, at1, at2, at3, at4, at5
      double precision at6, at7, at8, at9, at10, at11
      double precision bt0, bt1,  bt2, bt3, bt4, bt5
      double precision bt6, bt7, bt8, bt9, bt10, bt11
      double precision rhoa, rhob, rhoval, gammaval, TA, TB
      double precision E1, E2
      double precision de1drb, de2dga, de2dgb, de2dgc, de2dk, de2dra
      double precision de2drb, df1dr, df1dtau, df1dw, df2dr, df2dtau
      double precision de2dgaa, de2dgab, de2dgbb
      double precision de1dk, de1dra
      double precision df2dw, dggadga, dggadgb, dggadgc, dggadra
      double precision dggadrb, dgrhody, dhdga, dhdgb, dhdgc, dhdy
      double precision dlsdadra, dlsdadrb, dtdr, dtdtau, dwdt
      double precision fsig1, fsig2
      double precision taup, tauueg, tsig, wsig 
      
      double precision nepsc, dnepscdn(2)
      double precision epsc, depscdna, depscdnb
      double precision H0, dH0dna, dH0dnb, dH0dg

      double precision phi, dphidna, dphidnb, dphidzeta
      double precision zeta, dzetadna, dzetadnb
      double precision arglog, darglogdna, darglogdnb, darglogdg
      double precision fAt, dfAtdt, dfAtdA
      double precision fAtnum, dfAtnumdt, dfAtnumdA
      double precision fAtden, dfAtdendt, dfAtdendA
      double precision dfAtdna, dfAtdnb, dfAtdg
      double precision A, dAdna, dAdnb
      double precision t, dtdna, dtdnb, dtdg
      double precision ks, dksdna, dksdnb
      double precision argexp, dargexpdna, dargexpdnb
      double precision expinA
    
      double precision LSDA,GGA
 
      double precision f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11
      double precision f1o3, f2o3, f4o3, f5o3     
      Save F1, F2, F3, F4, F5, F6, F7, F8, F9, F10, F11
      DATA F1/1.0D+00/,  F2/2.0D+00/,  F3/3.0D+00/,
     $     F4/4.0D+00/,  F5/5.0D+00/,  F6/6.0D+00/,
     $     F7/7.0D+00/,  F8/8.0D+00/,  F9/9.0D+00/, 
     $     F10/10.0D+00/,F11/11.0D+00/

      
        at0=0.000000D+00
        at1=0.000000D+00
        at2=0.000000D+00
        at3=0.000000D+00
        at4=0.000000D+00
        at5=0.000000D+00
        at6=0.000000D+00
        at7=0.000000D+00
        at8=0.000000D+00
        at9=0.000000D+00
        at10=0.000000D+00 
        at11=0.000000D+00 
        bt0=0.000000D+00
        bt1=0.000000D+00
        bt2=0.000000D+00
        bt3=0.000000D+00
        bt4=0.000000D+00
        bt5=0.000000D+00
        bt6=0.000000D+00
        bt7=0.000000D+00
        bt8=0.000000D+00
        bt9=0.000000D+00
        bt10=0.000000D+00 
        bt11=0.000000D+00 
c
c     ======> BOTH SPIN-RESTRICETED AND UNRESTRICTED <======
c

      if (ijzy.eq.1) then
C     Parameters for M08-HX
        at0=    1.0000000D+00
        at1=   -4.0661387D-01
        at2=   -3.3232530D+00
        at3=    1.5540980D+00
        at4=    4.4248033D+01
        at5=   -8.4351930D+01
        at6=   -1.1955581D+02
        at7=    3.9147081D+02
        at8=    1.8363851D+02
        at9=   -6.3268223D+02
        at10=  -1.1297403D+02
        at11=   3.3629312D+02

        bt0=    1.3812334D+00
        bt1=   -2.4683806D+00
        bt2=   -1.1901501D+01
        bt3=   -5.4112667D+01
        bt4=    1.0055846D+01
        bt5=    1.4800687D+02
        bt6=    1.1561420D+02
        bt7=    2.5591815D+02
        bt8=    2.1320772D+02
        bt9=   -4.8412067D+02
        bt10=  -4.3430813D+02
        bt11=   5.6627964D+01
       elseif (ijzy.eq.2) then
C     Parameters for M08-SO 
        at0=    1.0000000D+00
        at1=    0.0000000D+00
        at2=   -3.9980886D+00
        at3=    1.2982340D+01
        at4=    1.0117507D+02
        at5=   -8.9541984D+01
        at6=   -3.5640242D+02
        at7=    2.0698803D+02
        at8=    4.6037780D+02
        at9=   -2.4510559D+02
        at10=  -1.9638425D+02
        at11=   1.1881459D+02

        bt0=    1.0000000D+00
        bt1=   -4.4117403D+00
        bt2=   -6.4128622D+00
        bt3=    4.7583635D+01
        bt4=    1.8630053D+02
        bt5=   -1.2800784D+02
        bt6=   -5.5385258D+02
        bt7=    1.3873727D+02
        bt8=    4.1646537D+02
        bt9=   -2.6626577D+02
        bt10=   5.6676300D+01
        bt11=   3.1673746D+02
       elseif (ijzy.eq.3) then
C     Parameters for M11
        at0=   1.0000000D+00
        at1=   0.0000000D+00
        at2=  -3.8933250D+00
        at3=  -2.1688455D+00
        at4=   9.3497200D+00
        at5=  -1.9845140D+01
        at6=   2.3455253D+00
        at7=   7.9246513D+01
        at8=   9.6042757D+00
        at9=  -6.7856719D+01
        at10= -9.1841067D+00
        at11=  0.0000000D+00

        bt0=   7.2239798D-01
        bt1=   4.3730564D-01
        bt2=  -1.6088809D+01
        bt3=  -6.5542437D+01
        bt4=   3.2057230D+01
        bt5=   1.8617888D+02
        bt6=   2.0483468D+01
        bt7=  -7.0853739D+01
        bt8=   4.4483915D+01
        bt9=  -9.4484747D+01
        bt10= -1.1459868D+02
        bt11=  0.0000000D+00
       elseif (ijzy.eq.4) then
C     Parameters for M11-L
        at0=   1.000000D+00
        at1=   0.000000D+00
        at2=   2.750880D+00
        at3=  -1.562287D+01
        at4=   9.363381D+00
        at5=   2.141024D+01
        at6=  -1.424975D+01
        at7=  -1.134712D+01
        at8=   1.022365D+01
        at9=   0.000000D+00
        at10=  0.000000D+00
        at11=  0.000000D+00
C
        bt0=   1.000000D+00
        bt1=  -9.082060D+00
        bt2=   6.134682D+00
        bt3=  -1.333216D+01
        bt4=  -1.464115D+01
        bt5=   1.713143D+01
        bt6=   2.480738D+00
        bt7=  -1.007036D+01
        bt8=  -1.117521D-01
        bt9=   0.000000D+00
        bt10=  0.000000D+00
        bt11=  0.000000D+00
       endif


      F1o3 = F1/F3
      F2o3 = F2/F3
      F4o3 = F4/F3
      F5o3 = F5/F3
      


c
c     ======> BOTH SPIN-RESTRICETED AND UNRESTRICTED <======
c
      do 20 n = 1, nq

      if (ipol.eq.1) then
        RHOA = rho(n,1)/F2
        RHOB = RHOA
        rho_t(1) = rho(n,1)        
        gammaval =(delrho(n,1,1)*delrho(n,1,1) +
     &             delrho(n,2,1)*delrho(n,2,1) +
     &             delrho(n,3,1)*delrho(n,3,1))
        TA = tau(n,1)/F2
        TB=TA
      else
        RHOA = rho(n,2)
        RHOB = rho(n,3)
        rho_t(1) = rho(n,1)
        rho_t(2) = rho(n,2)
        rho_t(3) = rho(n,3)
        gammaval = delrho(n,1,1)*delrho(n,1,1) +
     &             delrho(n,1,2)*delrho(n,1,2) +
     &             delrho(n,2,1)*delrho(n,2,1) +
     &             delrho(n,2,2)*delrho(n,2,2) +
     &             delrho(n,3,1)*delrho(n,3,1) +
     &             delrho(n,3,2)*delrho(n,3,2) +
     &       2.d0*(delrho(n,1,1)*delrho(n,1,2) +
     &             delrho(n,2,1)*delrho(n,2,2) +
     &             delrho(n,3,1)*delrho(n,3,2))
        TA =  tau(n,1)
        TB =  tau(n,2)
      endif
      RHOVAL  = RHOA + RHOB
      TauP = TA+TB
      dsqgamma = max(dsqrt(gammaval),tol_rho)
      if (dabs(rhoval).le.tol_rho) goto 20
      if(taup.lt.tol_rho) goto 20 
c
      nepsc = 0.0d0
      dnepscdn(1) = 0.0d0
      if (ipol.eq.2) dnepscdn(2) = 0.0d0
c
c        call for LDA bit
c        this implementation temporarily assigns the pw91LDA for
c        use in the metaGGA local part
c
            call xc_pw91lda(tol_rho,1d0,.true.,.false.,rho_t,
     &         dnepscdn,1,ipol,nepsc,1d0,
     &         .false.,dummy)
        
c
c        PBE non-local part
c        This part is a modified version of xc_cMpbe96.F
c
         if(abs(nepsc).lt.tol_rho*tol_rho) goto 20
c
         epsc = nepsc/rhoval
         if (ipol.eq.1) then
            depscdna = dnepscdn(1)/rhoval-nepsc/(rhoval**2)
            depscdnb = depscdna
         else
            depscdna = dnepscdn(1)/rhoval-nepsc/(rhoval**2)
            depscdnb = dnepscdn(2)/rhoval-nepsc/(rhoval**2)
         endif
c
         ks = 2.0d0*(3.0d0*PI*PI*rhoval)**(1.0d0/6.0d0)/dsqrt(PI)
         dksdna = (1.0d0/6.0d0)*ks/rhoval
         dksdnb = dksdna
c
         if (ipol.eq.1) then
            zeta = 0.0d0
         else
            zeta = (rho(n,2)-rho(n,3))/rhoval
         endif
         if(zeta.lt.-1.0d0) zeta=-1.0d0
         if(zeta.gt. 1.0d0) zeta= 1.0d0
         if (ipol.eq.1) then
            dzetadna = 1.0d0/rhoval
            dzetadnb = -1.0d0/rhoval
         else
            dzetadna =  2.0d0*rho(n,3)/(rhoval**2)
            dzetadnb = -2.0d0*rho(n,2)/(rhoval**2)
         endif
c
         phi = 0.5d0*((1.0d0+zeta)**(2.0d0/3.0d0)
     &               +(1.0d0-zeta)**(2.0d0/3.0d0))
         if (dabs(1.0d0-zeta).lt.tol_rho) then
            dphidzeta = 0.5d0*(2.0d0/3.0d0)*(
     &             (1.0d0+zeta)**(2.0d0/3.0d0)/(1.0d0+zeta))
         else if (dabs(1.0d0+zeta).lt.tol_rho) then
            dphidzeta = 0.5d0*(2.0d0/3.0d0)*(
     &            -(1.0d0-zeta)**(2.0d0/3.0d0)/(1.0d0-zeta))
         else
            dphidzeta = 0.5d0*(2.0d0/3.0d0)*(
     &         (1.0d0+zeta)**(2.0d0/3.0d0)/(1.0d0+zeta)
     &        -(1.0d0-zeta)**(2.0d0/3.0d0)/(1.0d0-zeta))
         endif
         dphidna = dphidzeta*dzetadna
         dphidnb = dphidzeta*dzetadnb
c
         t = dsqgamma/(2.0d0*phi*ks*rhoval)
         dtdna = -t/rhoval-t/phi*dphidna-t/ks*dksdna
         dtdnb = -t/rhoval-t/phi*dphidnb-t/ks*dksdnb
c
         argexp = -epsc/GAMMA/(phi**3)
         dargexpdna = -depscdna/GAMMA/(phi**3)
     &                +3.0d0*epsc/GAMMA/(phi**4)*dphidna
         dargexpdnb = -depscdnb/GAMMA/(phi**3)
     &                +3.0d0*epsc/GAMMA/(phi**4)*dphidnb
c
         if (dabs(argexp).lt.EXPTOL) then
            expinA=dexp(argexp)
         else
            expinA=0.0d0
         endif
         A = BETA/GAMMA/(expinA-1.0d0)
         dAdna = -BETA/GAMMA*dargexpdna*expinA/(expinA-1.0d0)**2
         dAdnb = -BETA/GAMMA*dargexpdnb*expinA/(expinA-1.0d0)**2
c
         fAtnum = 1.0d0+A*t**2
         fAtden = 1.0d0+A*t**2+A**2*t**4
         fAt = fAtnum/fAtden
         dfAtnumdt = 2.0d0*A*t
         dfAtnumdA = t**2
         dfAtdendt = 2.0d0*A*t+4.0d0*A**2*t**3
         dfAtdendA = t**2+2.0d0*A*t**4
         dfAtdt = (dfAtnumdt*fAtden-fAtnum*dfAtdendt)/(fAtden**2)
         dfAtdA = (dfAtnumdA*fAtden-fAtnum*dfAtdendA)/(fAtden**2)
         dfAtdna = dfAtdt * dtdna + dfAtdA * dAdna
         dfAtdnb = dfAtdt * dtdnb + dfAtdA * dAdnb
c
         arglog = 1.0d0 + BETA/GAMMA*t**2*fAt
         darglogdna = BETA/GAMMA*(2.0d0*t*dtdna*fAt
     &                            +t*t*dfAtdna)
         darglogdnb = BETA/GAMMA*(2.0d0*t*dtdnb*fAt
     &                            +t*t*dfAtdnb)
c
         H0 = GAMMA*(phi**3)*dlog(arglog)
         dH0dna = GAMMA*(3.0d0*(phi**2)*dphidna*dlog(arglog)
     &                  +(phi**3)*darglogdna/arglog)
         dH0dnb = GAMMA*(3.0d0*(phi**2)*dphidnb*dlog(arglog)
     &                  +(phi**3)*darglogdnb/arglog)       
c
        TauUEG=F3*(F3*Pi*Pi)**F2o3*rhoval**F5o3/F10
        Tsig =TauUEG/TauP
        Wsig =(Tsig - F1)/(Tsig + F1)
        Fsig1=(at0 + Wsig*(at1 + Wsig*(at2 + Wsig*(at3 + Wsig*(
     &            at4 + Wsig*(at5 + Wsig*(at6 + Wsig*(at7 + Wsig*(
     &            at8 + Wsig*(at9 + Wsig*(at10+Wsig*at11)))))))))))

        Fsig2=(bt0 + Wsig*(bt1 + Wsig*(bt2 + Wsig*(bt3 + Wsig*(
     &            bt4 + Wsig*(bt5 + Wsig*(bt6 + Wsig*(bt7 + Wsig*(
     &            bt8 + Wsig*(bt9 + Wsig*(bt10+Wsig*bt11)))))))))))
        LSDA = rhoval*epsc
        GGA  = rhoval*H0
        E1   = LSDA*Fsig1
        E2   = GGA*Fsig2
c
c        derivatives of the meta-GGA part, and total
c
         dF1dW=( at1 + Wsig*(F2  *at2 + Wsig*(F3*at3 + Wsig*(
     &            F4 *at4 + Wsig*(F5 *at5 + Wsig*(F6  *at6 + Wsig*(
     &            F7*at7 + Wsig*(F8*at8 + Wsig*(F9 *at9 + Wsig*(
     &            F10  *at10+ Wsig*F11*at11))))))))))
         dF2dW=( bt1 + Wsig*(F2  *bt2 + Wsig*(F3*bt3 + Wsig*(
     &            F4 *bt4 + Wsig*(F5 *bt5 + Wsig*(F6  *bt6 + Wsig*(
     &            F7*bt7 + Wsig*(F8*bt8 + Wsig*(F9 *bt9 + Wsig*(
     &            F10  *bt10+ Wsig*F11*bt11))))))))))
         dWdT = F2/((F1 + Tsig)**F2)
         dTdR = Tsig*F5/(F3*rhoval) 
         dTdTau = -Tsig/TauP
         dF1dR = dF1dW*dWdT*dTdR
         dF1dTau=dF1dW*dWdT*dTdTau
         dF2dR = dF2dW*dWdT*dTdR
         dF2dTau=dF2dW*dWdT*dTdTau

         dLSDAdRA=dnepscdn(1)
         if (ipol.eq.2) dLSDAdRB=dnepscdn(2)
         dGGAdRA=rhoval*dH0dna+H0
         if (ipol.eq.2) dGGAdRB=rhoval*dH0dnb+H0
C
         dE1dRA = dLSDAdRA*Fsig1 + LSDA*dF1dR
         if (ipol.eq.2) dE1dRB = dLSDAdRB*Fsig1 + LSDA*dF1dR

C
         dE2dRA = dGGAdRA*Fsig2 + GGA*dF2dR
         if (ipol.eq.2) dE2dRB = dGGAdRB*Fsig2 + GGA*dF2dR 
c
         dtdg = 0.25d0/(phi*ks*rhoval)/dsqgamma
         dfAtdg = dfAtdt*dtdg
         darglogdg = BETA/GAMMA*(2.0d0*t*dtdg*fAt+t*t*dfAtdg)
         dH0dg = GAMMA*(phi**3)*darglogdg/arglog
C
         dE2dGAA =    rhoval*dH0dg*Fsig2
         dE2dGAB = F2*rhoval*dH0dg*Fsig2
         dE2dGBB =    rhoval*dH0dg*Fsig2
c
         dE1dK = LSDA*dF1dTau
         dE2dK = GGA*dF2dTau
c
         Ec = Ec+ (E1+E2)*qwght(n)
         if(ldew) func(n) = func(n)+(E1+E2)
         Amat(n,1) = Amat(n,1) + (dE1dRA + dE2dRA)
         if (ipol.eq.2) Amat(n,2) = Amat(n,2) + (dE1dRB+dE2dRB)
         Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + dE2dGAA
         Cmat(n,D1_GAB) = Cmat(n,D1_GAB) + dE2dGAB
         if (ipol.eq.2) Cmat(n,D1_GBB) = Cmat(n,D1_GBB) + dE2dGBB         
         Mmat(n,1)=Mmat(n,1)+ (de1dk+de2dk)/F2
         if (ipol.eq.2) Mmat(n,2)=Mmat(n,2)+(de1dk+de2dk)/F2         
   20 continue
      end
c
      Subroutine xc_cm11_d2()
      implicit none
      call errquit(' not coded ',0,0)
      return
      end
