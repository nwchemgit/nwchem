      Subroutine xc_perdew91(tol_rho, cfac, lcfac, nlcfac, rho, delrho, 
     &                       Amat, Cmat, nq, ipol, Ec, qwght, 
     ,     ldew, ffunc)
c
c     $Id: xc_perdew91.F,v 1.27 2001-04-12 02:28:27 edo Exp $
c
      Implicit none
c
#include "dft2drv.fh"
c
      double precision TOLL, EXPTOL
      double precision alpha, beta, pgamma, delta, beta10
      double precision zzz, fff, CINF
      double precision CX37, aa
      double precision dummy(1)
      double precision ONE, ONE3, THREE
      double precision FOURTHIRDS, SEV6
      double precision FIVE3, TWO3, FIVE6
      double precision pi
      double precision eps
      double precision cfac(15)
      logical lcfac(15), nlcfac(15)
      logical lfac, nlfac, lfacl, nlfacl,ldew
      double precision ffunc(*)  ! value of the functional [output]
      double precision fac, facl
      integer lnq
      double precision lqwght
c
      Parameter (TOLL = 1.D-40, EXPTOL = 40.d0)
      Parameter (alpha = 0.023266D0, beta = 7.389D-6, 
     &           pgamma = 8.723d0, delta = 0.472d0,  
     &           beta10 = 10000.d0*beta)
      parameter (zzz = 0.001667d0, fff = 0.002568d0)
      parameter (CINF = zzz+fff)
      parameter (CX37 = -zzz*3.d0/7.d0)
      parameter (aa = 0.09d0)
      Parameter (ONE = 1.D0, ONE3 = 1.d0/3.d0, THREE = 3.d0)
      Parameter (FOURTHIRDS = 4.D0/3.D0, SEV6 = 7.d0/6.d0)
      parameter (FIVE3 = 5.d0/3.d0, TWO3 = 2.d0/3.d0, FIVE6 = 5.d0/6.d0)
      parameter (pi = 3.1415926535897932385d0)
      parameter (eps = 1.e-8)
c
      integer ipol, nq
      double precision tol_rho
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
c     Quadrature Weights
c
      double precision qwght(nq)
c
c     Sampling Matrices for the XC Potential & Energy
c
      double precision Amat(nq,ipol), Cmat(nq,*)
      double precision dp81drho(2), rho_t(3)
c     
c     J.P.Perdew, J.A.Chevary, S.H.Vosko, K.A.Jackson, M.R.Pederson, 
c     D.J.~Singh and C.~Fiolhais,  Phys. Rev. B  46, 6671 (1992).     
c     
      integer n
      double precision expo1, expo2, dexpo1, dexpo2
      double precision TPP, nu, bb, twoab, twoab2, atwob2
      double precision rsfact, rhoval, rho13, rs, rs2, rs3, 
     &                 drsdrho, arho, drho
      double precision gammaval, dsqgamma, dtdgamma
      double precision anum, aden, CC, CCPRIME
      double precision egas, degasda, degasdb
      double precision g, g2, g3
      double precision gzeta, gzprime, gzprime_1, gzprime_m1, x
      double precision dgzetadag, dgzetadbg, gggr
      double precision akf, aks, dksdrho, dkfdrho
      double precision t, t2, t3, t4, dtdrhoa, dtdrhob
      double precision A, A1p2At2, aden1, alognd, dAdrhoa, dAdrhob
      double precision func, H0, H1, dH1drhoa, dH0drhoa, 
     &                 dH1drhob, dH0drhob
      double precision dH1dgamma, dH0dgamma
      double precision gaa, gbb
      double precision zeta
c
      gzeta(x) = 0.5d0*((ONE + x)**TWO3 + 
     &                  (abs(ONE - x))**TWO3)
      gzprime(x) = 0.5d0*TWO3*((ONE + x)**(-ONE3) - 
     &                         (abs(ONE - x))**(-ONE3))
c
      gzprime_1  =  0.5d0*TWO3*(2.0d0**(-ONE3))
      gzprime_m1 = -0.5d0*TWO3*(2.0d0**(-ONE3))
      TPP = (3.d0*pi*pi)**ONE3
      nu = (16.d0/pi)*TPP
      bb = nu*CINF
      twoab = (aa+aa)/bb
      twoab2 = twoab/bb
      atwob2 = 1.d0/twoab2
      fac = cfac(5)
      lfac = lcfac(5)
      nlfac = nlcfac(5)
      rsfact = (0.75d0/pi)**ONE3
c
      if (ipol.eq.1)then
c
c        ======> SPIN-RESTRICTED <======
c
         do 10 n = 1, nq
            rhoval = rho(n,1)
            if (abs(rhoval).le.tol_rho)goto 10
            arho = 1.d0/rhoval
            rho13  = rhoval**ONE3
            rs = rsfact/rho13
            rs2 = rs*rs
            rs3 = rs2*rs
            drsdrho = -ONE3*rs*arho
            gammaval = delrho(n,1,1)*delrho(n,1,1) +
     &                 delrho(n,2,1)*delrho(n,2,1) +
     &                 delrho(n,3,1)*delrho(n,3,1)
            dsqgamma = dsqrt(gammaval)
c         
c           C(n)
c         
            anum = fff+alpha*rs+beta*rs2
            aden = 1.d0/(1.d0+pgamma*rs+delta*rs2+beta10*rs3)
            CC = zzz+(anum*aden)
            CCPRIME = aden*((alpha+beta*(rs+rs)) -
     &                (anum*aden)*(pgamma+(rs+rs)*delta + 
     &                beta10*rs2*THREE))
            CCPRIME = CCPRIME*drsdrho
            egas = 0.d0
            dp81drho(1) = 0.d0
c
c           call for LDA bit
c
            lnq = 1
            lqwght = 1.0d0
c
            if (abs(cfac(1)).gt.eps)then
                 facl = cfac(1)
                lfacl = lcfac(1)
               nlfacl = nlcfac(1)
               call xc_vwn_5(tol_rho,facl,lfacl,nlfacl,rhoval,dp81drho,
     &                       lnq,ipol,egas,lqwght,.false.,dummy)
            endif
c
            if (abs(cfac(3)).gt.eps)then
                 facl = cfac(3)
                lfacl = lcfac(3)
               nlfacl = nlcfac(3)
               call xc_p81(tol_rho,facl,lfacl,nlfacl,rhoval,dp81drho,
     &                     lnq,ipol,egas,lqwght,.false.,dummy)
            endif
c
            if (abs(cfac(6)).gt.eps.or.lfac)then
               facl = cfac(6)
               if (abs(facl).lt.eps)facl = 1.0d0
                lfacl = lcfac(6)
               nlfacl = nlcfac(6)
               call xc_pw91lda(tol_rho,facl,lfacl,nlfacl,rhoval,
     &                         dp81drho,lnq,ipol,egas,lqwght,
     ,              .false.,dummy)
            endif
c
            if (abs(cfac(7)).gt.eps)then
                 facl = cfac(7)
                lfacl = lcfac(7)
               nlfacl = nlcfac(7)
               call xc_vwn_1_rpa(tol_rho,facl,lfacl,nlfacl,rhoval,
     &                           dp81drho,lnq,ipol,egas,lqwght,
     &                           .false.,dummy)
            endif
c
            if (abs(cfac(8)).gt.eps)then
                 facl = cfac(8)
                lfacl = lcfac(8)
               nlfacl = nlcfac(8)
               call xc_vwn_1(tol_rho,facl,lfacl,nlfacl,rhoval,dp81drho,
     &                       lnq,ipol,egas,lqwght,.false.,dummy)
            endif
c
            if (abs(cfac(9)).gt.eps)then
                 facl = cfac(9)
                lfacl = lcfac(9)
               nlfacl = nlcfac(9)
               call xc_vwn_2(tol_rho,facl,lfacl,nlfacl,rhoval,dp81drho,
     &                       lnq,ipol,egas,lqwght,.false.,dummy)
            endif
c
            if (abs(cfac(10)).gt.eps)then
                 facl = cfac(10)
                lfacl = lcfac(10)
               nlfacl = nlcfac(10)
               call xc_vwn_3(tol_rho,facl,lfacl,nlfacl,rhoval,dp81drho,
     &                       lnq,ipol,egas,lqwght,.false.,dummy)
            endif
c
            if (abs(cfac(11)).gt.eps)then
                 facl = cfac(11)
                lfacl = lcfac(11)
               nlfacl = nlcfac(11)
               call xc_vwn_4(tol_rho,facl,lfacl,nlfacl,rhoval,dp81drho,
     &                       lnq,ipol,egas,lqwght,.false.,dummy)
            endif
c
            degasda = arho*(dp81drho(1) - egas*arho)
            akf = TPP*rho13
            aks = sqrt(4.d0*akf/pi)
            dksdrho = aks*arho/6.d0
            dkfdrho = akf*arho*ONE3
            t = 0.5d0*dsqgamma*arho/aks
            t2 = t*t
            t3 = t2*t
            t4 = t3*t
            dtdrhoa = -t*arho*SEV6
            dtdgamma = 0.5d0*t/gammaval
            dexpo1 = twoab2*(egas*arho)
            dexpo2 = (10.d0*(aks/akf)*t)**2
            if (dexpo1.gt.EXPTOL)then
               expo1 = 0.d0
            else
               expo1 = exp(-dexpo1)
            endif
            if (dexpo2.gt.EXPTOL)then
               expo2 = 0.d0
            else
               expo2 = exp(-dexpo2)
            endif
            A = twoab/(expo1-ONE)
            A1p2At2 = ONE+2.d0*A*t2
            aden1 = ONE/(ONE+A*t2+A*A*t4)
            alognd = (t2+A*t4)*aden1
            if (twoab*alognd.gt.TOLL)then
               H0 = atwob2*log(ONE+twoab*alognd)
            else
               H0 = 0.d0
            endif
            H1 = nu*(CC-CINF-CX37)*t2*expo2
c
            func = H0 + H1
c
            if (lfac)then
               Ec = Ec + egas*qwght(n)*fac
               if(ldew) ffunc(n)=ffunc(n)+egas*fac
            endif
            if (nlfac)then
               Ec = Ec + (func*rhoval)*qwght(n)*fac
               if(ldew) ffunc(n)=ffunc(n)+(func*rhoval)*fac
            endif
c
c           df/drho = func + rhoval*dfunc/drho
c
            dH1drhoa = H1*(CCPRIME/(CC-CINF-CX37) + 2.d0*dtdrhoa/t -
     &                     2.d0*dexpo2*(dtdrhoa/t + dksdrho/aks - 
     &                     dkfdrho/akf))
            dAdrhoa = A/(expo1-ONE)*expo1*twoab2*degasda
            dH0drhoa = bb*aden1*A1p2At2*(2.d0*t*dtdrhoa + 
     &                 dAdrhoa*t4/A1p2At2 - alognd*(2.d0*A*t*dtdrhoa +
     &                 dAdrhoa*t2))/(ONE+twoab*alognd)
            if (lfac)then
               Amat(n,1) = Amat(n,1) +  dp81drho(1)*fac
            endif
            if (nlfac)then
               Amat(n,1) = Amat(n,1) +  (func +
     &                     rhoval*(dH0drhoa+dH1drhoa))*fac
            endif
c
c           df/dgamma
c
            if (dsqgamma.gt.TOLL)then
               dH1dgamma = H1/gammaval*(ONE - dexpo2)
               dH0dgamma = bb/(ONE + twoab*alognd)*
     &                     aden1*2.d0*t*dtdgamma*A1p2At2*
     &                     (ONE - alognd*A)
               gaa = (dH0dgamma + dH1dgamma)*rhoval*fac*2.d0
               Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + gaa
#if 0
               Cmat(n,1,1) = Cmat(n,1,1) + gaa*delrho(n,1,1)
               Cmat(n,2,1) = Cmat(n,2,1) + gaa*delrho(n,2,1)
               Cmat(n,3,1) = Cmat(n,3,1) + gaa*delrho(n,3,1)
#endif
            endif
   10    continue
      else
c
c        ======> SPIN-UNRESTRICTED <======
c
         do 20 n = 1, nq
            rhoval = rho(n,1)
            rho_t(1) = rho(n,1)
            rho_t(2) = rho(n,2)
            rho_t(3) = rho(n,3)
            if (abs(rhoval).le.tol_rho) goto 20
            arho = 1.d0/rhoval
            rho13  = rhoval**ONE3
            rs = rsfact/rho13
            rs2 = rs*rs
            rs3 = rs2*rs
            drsdrho = -ONE3*rs*arho
            gammaval = delrho(n,1,1)*delrho(n,1,1) +
     &                 delrho(n,1,2)*delrho(n,1,2) +
     &                 delrho(n,2,1)*delrho(n,2,1) +
     &                 delrho(n,2,2)*delrho(n,2,2) +
     &                 delrho(n,3,1)*delrho(n,3,1) +
     &                 delrho(n,3,2)*delrho(n,3,2) +
     &           2.d0*(delrho(n,1,1)*delrho(n,1,2) +
     &                 delrho(n,2,1)*delrho(n,2,2) +
     &                 delrho(n,3,1)*delrho(n,3,2))
            dsqgamma = max(dsqrt(gammaval),tol_rho)
            drho = rho(n,2) - rho(n,3)
            zeta = drho*arho
c         
c           C(n)
c         
            anum = fff+alpha*rs+beta*rs2
            aden = 1.d0/(1.d0+pgamma*rs+delta*rs2+beta10*rs3)
            CC = zzz+(anum*aden)
            CCPRIME = aden*((alpha+beta*(rs+rs)) -
     &                (anum*aden)*(pgamma+(rs+rs)*delta + 
     &                beta10*rs2*THREE))
            CCPRIME = CCPRIME*drsdrho
            egas = 0.d0
            dp81drho(1) = 0.d0
            dp81drho(2) = 0.d0
c
c           call for LDA bit
c
            lnq = 1
            lqwght = 1.0d0
c
            if (abs(cfac(1)).gt.eps)then
                 facl = cfac(1)
                lfacl = lcfac(1)
               nlfacl = nlcfac(1)
               call xc_vwn_5(tol_rho,facl,lfacl,nlfacl,rho_t,
     &                       dp81drho,lnq,ipol,egas,lqwght,
     &                       .false.,dummy)
            endif
c
            if (abs(cfac(3)).gt.eps)then
                 facl = cfac(3)
                lfacl = lcfac(3)
               nlfacl = nlcfac(3)
               call xc_p81(tol_rho,facl,lfacl,nlfacl,rho_t,dp81drho,
     &                     lnq,ipol,egas,lqwght,.false.,dummy)
            endif
c
            if (abs(cfac(6)).gt.eps.or.lfac)then
                 facl = cfac(6)
                 if (abs(facl).lt.eps)facl = 1.0d0
                lfacl = lcfac(6)
               nlfacl = nlcfac(6)
               call xc_pw91lda(tol_rho,facl,lfacl,nlfacl,rho_t,
     &                         dp81drho,lnq,ipol,egas,lqwght,
     ,              .false.,dummy)
            endif
c
            if (abs(cfac(7)).gt.eps)then
                 facl = cfac(7)
                lfacl = lcfac(7)
               nlfacl = nlcfac(7)
               call xc_vwn_1_rpa(tol_rho,facl,lfacl,nlfacl,rho_t,
     &                           dp81drho,lnq,ipol,egas,lqwght,
     &                           .false.,dummy)
            endif
c
            if (abs(cfac(8)).gt.eps)then
                 facl = cfac(8)
                lfacl = lcfac(8)
               nlfacl = nlcfac(8)
               call xc_vwn_1(tol_rho,facl,lfacl,nlfacl,rho_t,
     &                       dp81drho,lnq,ipol,egas,lqwght,
     &                       .false.,dummy)
            endif
c
            if (abs(cfac(9)).gt.eps)then
                 facl = cfac(9)
                lfacl = lcfac(9)
               nlfacl = nlcfac(9)
               call xc_vwn_2(tol_rho,facl,lfacl,nlfacl,rho_t,
     &                       dp81drho,lnq,ipol,egas,lqwght,
     &                       .false.,dummy)
            endif
c
            if (abs(cfac(10)).gt.eps)then
                 facl = cfac(10)
                lfacl = lcfac(10)
               nlfacl = nlcfac(10)
               call xc_vwn_3(tol_rho,facl,lfacl,nlfacl,rho_t,
     &                       dp81drho,lnq,ipol,egas,lqwght,
     &                       .false.,dummy)
            endif
c
            if (abs(cfac(11)).gt.eps)then
                 facl = cfac(11)
                lfacl = lcfac(11)
               nlfacl = nlcfac(11)
               call xc_vwn_4(tol_rho,facl,lfacl,nlfacl,rho_t,
     &                       dp81drho,lnq,ipol,egas,lqwght,
     &                       .false.,dummy)
            endif
c
            degasda = arho*(dp81drho(1) - egas*arho)
            degasdb = arho*(dp81drho(2) - egas*arho)
            g = max(tol_rho,gzeta(zeta))
            g2 = g*g
            g3 = g2*g
            akf = TPP*rho13
            aks = sqrt(4.d0*akf/pi)
            dksdrho = aks*arho/6.d0
            dkfdrho = akf*arho*ONE3
            if (abs(1.d0-zeta).lt.tol_rho)then
               dgzetadag = 0.d0
               dgzetadbg = -2.d0*(gzprime_1/g)*arho
            elseif (abs(1.d0+zeta).lt.tol_rho)then
               dgzetadag = 0.d0
               dgzetadbg = -2.d0*(gzprime_m1/g)*arho
            else
               gggr = (gzprime(zeta)/g)*arho
               dgzetadag = (1.d0-zeta)*gggr
               dgzetadbg = -(1.d0+zeta)*gggr
            endif
c
            t = 0.5d0*dsqgamma*arho/(g*aks)
            t2 = t*t
            t3 = t2*t
            t4 = t3*t
            dtdrhoa = -t*(dgzetadag+dksdrho/aks+arho)
            dtdrhob = -t*(dgzetadbg+dksdrho/aks+arho)
            dtdgamma = 0.5d0*t/gammaval
c
            dexpo1 = twoab2*(egas*arho)/g3
            dexpo2 = (10.d0*g2*(aks/akf)*t)**2
c
            expo1 = TOLL
            if (dexpo1.lt.EXPTOL)expo1 = exp(-dexpo1)
            expo2 = TOLL
            if (dexpo2.lt.EXPTOL)expo2 = exp(-dexpo2)
c
            A = twoab/(expo1-ONE)
            A1p2At2 = ONE+2.d0*A*t2
            aden1 = ONE/(ONE+A*t2+A*A*t4)
            alognd = (t2+A*t4)*aden1
            H0 = atwob2*g3*log(ONE+twoab*alognd)
            H1 = nu*(CC-CINF-CX37)*g3*t2*expo2
c
            func = H0 + H1
c
            if (lfac)then
               Ec = Ec + egas*qwght(n)*fac
               if(ldew) ffunc(n)=ffunc(n)+egas*fac
            endif
            if (nlfac)then
               Ec = Ec + (func*rhoval)*qwght(n)*fac
               if(ldew) ffunc(n)=ffunc(n)+(func*rhoval)*fac
            endif
c
c           df/drho = func + rhoval*dfunc/drho
c
            dH1drhoa = H1*(CCPRIME/(CC-CINF-CX37) +
     &                     THREE*dgzetadag + 2.d0*dtdrhoa/t -
     &                     2.d0*dexpo2*(2.d0*dgzetadag +
     &                     dtdrhoa/t + dksdrho/aks - dkfdrho/akf))
            dH1drhob = H1*(CCPRIME/(CC-CINF-CX37) +
     &                     THREE*dgzetadbg + 2.d0*dtdrhob/t -
     &                     2.d0*dexpo2*(2.d0*dgzetadbg +
     &                     dtdrhob/t + dksdrho/aks - dkfdrho/akf))
            dAdrhoa = A/(expo1-ONE)*expo1*twoab2 /
     &                  g3*(-THREE*egas*arho*dgzetadag + degasda)
            dAdrhob = A/(expo1-ONE)*expo1*twoab2 /
     &                  g3*(-THREE*egas*arho*dgzetadbg + degasdb)
            dH0drhoa = Three*H0*dgzetadag + g3*bb*aden1*
     &                 (2.d0*t*dtdrhoa*A1p2At2 + dAdrhoa*t4 -
     &                  alognd*A1p2At2*(2.d0*A*t*dtdrhoa +
     &                  dAdrhoa*t2))/(ONE+twoab*alognd)
            dH0drhob = Three*H0*dgzetadbg + g3*bb*aden1*
     &                 (2.d0*t*dtdrhob*A1p2At2 + dAdrhob*t4 -
     &                  alognd*A1p2At2*(2.d0*A*t*dtdrhob +
     &                  dAdrhob*t2))/(ONE+twoab*alognd)
c
            if (lfac)then
               Amat(n,1) = Amat(n,1) +  dp81drho(1)*fac
               Amat(n,2) = Amat(n,2) +  dp81drho(2)*fac
            endif
            if (nlfac)then
               Amat(n,1) = Amat(n,1) +  (func +
     &                     rhoval*(dH0drhoa+dH1drhoa))*fac
               Amat(n,2) = Amat(n,2) +  (func +
     &                     rhoval*(dH0drhob+dH1drhob))*fac
            endif
c
            if (dsqgamma.gt.TOLL)then
               dH1dgamma = H1/gammaval*(ONE - dexpo2)
               dH0dgamma = g3*bb/(ONE + twoab*alognd)*
     &                     aden1*2.d0*t*dtdgamma*A1p2At2*
     &                     (ONE - alognd*A)
               gaa = 2.d0*(dH0dgamma + dH1dgamma)*rhoval*fac
               gbb = gaa
c               gbb = gaa*fac
               Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + gaa*0.5d0
               Cmat(n,D1_GAB) = Cmat(n,D1_GAB) + gaa
               Cmat(n,D1_GBB) = Cmat(n,D1_GBB) + gbb*0.5d0
#if 0
               Cmat(n,1,1) = Cmat(n,1,1) + gaa*(delrho(n,1,1) + 
     &                                          delrho(n,1,2))
               Cmat(n,2,1) = Cmat(n,2,1) + gaa*(delrho(n,2,1) + 
     &                                          delrho(n,2,2))
               Cmat(n,3,1) = Cmat(n,3,1) + gaa*(delrho(n,3,1) + 
     &                                          delrho(n,3,2))
               Cmat(n,1,2) = Cmat(n,1,2) + gbb*(delrho(n,1,2) + 
     &                                          delrho(n,1,1))
               Cmat(n,2,2) = Cmat(n,2,2) + gbb*(delrho(n,2,2) + 
     &                                          delrho(n,2,1))
               Cmat(n,3,2) = Cmat(n,3,2) + gbb*(delrho(n,3,2) + 
     &                                          delrho(n,3,1))
#endif
            endif
   20    continue
      endif
c
      return
      end

