      Subroutine xc_cpbe96(tol_rho, cfac, lcfac, nlcfac, rho, delrho, 
     &                    Amat, Cmat, nq, ipol, Ec, qwght,
     &                    ldew,ffunc)     
c
c$Id: xc_pbe96.F,v 1.10 2001-12-11 18:41:57 edo Exp $
c
      Implicit none
c
#include "dft2drv.fh"
c
      double precision TOLL, EXPTOL
      double precision dummy(1)
      double precision ONE, ONE3, THREE
      double precision SEV6
      double precision TWO3
      double precision pi
      double precision eps
      double precision cfac(15)
      logical lcfac(15), nlcfac(15),ldew
      logical lfac, nlfac, lfacl, nlfacl
      double precision ffunc(*)
      double precision fac, facl
      integer lnq
      double precision lqwght
c     pbe
      double precision bet,delt,gamma

c
      Parameter (TOLL = 1.D-40, EXPTOL = 40.d0)
      Parameter (ONE = 1.D0, ONE3 = 1.d0/3.d0, THREE = 3.d0)
      Parameter (SEV6 = 7.d0/6.d0)
      parameter (TWO3 = 2.d0/3.d0)
      parameter (pi = 3.1415926535897932385d0)
      parameter (eps = 1.e-8)
      parameter(gamma=0.03109069086965489503494086371273d0)
      parameter(bet=0.06672455060314922d0,delt=bet/gamma)
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
c References:
c [a] J.P.~Perdew, K.~Burke, and M.~Ernzerhof,
c     {\sl Generalized gradient approximation made simple}, sub.
c     to Phys. Rev.Lett. May 1996.
c [b] J. P. Perdew, K. Burke, and Y. Wang, {\sl Real-space cutoff
c     construction of a generalized gradient approximation:  The PW91
c     density functional}, submitted to Phys. Rev. B, Feb. 1996.
c [c] J. P. Perdew and Y. Wang, Phys. Rev. B {\bf 45}, 13244 (1992).
c     
      integer n
      double precision expo1,  dexpo1
      double precision TPP,  twoab, twoab2, atwob2
      double precision rhoval, rho13, 
     &                 arho, drho
      double precision gammaval, dsqgamma, dtdgamma
      double precision egas, degasda, degasdb
      double precision g, g2, g3
      double precision gzeta, gzprime, gzprime_1, gzprime_m1, x
      double precision dgzetadag, dgzetadbg, gggr
      double precision akf, aks, dksdrho
      double precision t, t2, t3, t4, dtdrhoa, dtdrhob
      double precision A, A1p2At2, aden1, alognd, dAdrhoa, dAdrhob
      double precision func, H0, dH0drhoa,  dH0drhob
      double precision dH0dgamma
      double precision gaa, gbb
      double precision zeta
c
      gzeta(x) = 0.5d0*((ONE + x)**TWO3 + 
     &     sign(1d0,one-x)*abs(ONE - x)**TWO3)
      gzprime(x) = 0.5d0*TWO3*((ONE + x)**(-ONE3) - 
     &     sign(1d0,one-x)*abs(ONE - x)**(-ONE3))
c
      gzprime_1  =  0.5d0*TWO3*(2.0d0**(-ONE3))
      gzprime_m1 = -0.5d0*TWO3*(2.0d0**(-ONE3))
      TPP = (3.d0*pi*pi)**ONE3
c
c     twoab becomes what is called DELT in pbe.asc
c
      twoab = delt
c
c     twoab2 becomes what is called 1/gamma in pbe.asc
c
      twoab2 = 1.d0/gamma
      atwob2 = gamma
      fac = cfac(12)
      lfac = lcfac(12)
      nlfac = nlcfac(12)
c
      if (ipol.eq.1)then
c
c        ======> SPIN-RESTRICTED <======
c
         do 10 n = 1, nq
            rhoval = rho(n,1)
            if (abs(rhoval).le.tol_rho)goto 10
            arho = 1.d0/rhoval
            rho13  = abs(rhoval)**ONE3
            gammaval = delrho(n,1,1)*delrho(n,1,1) +
     &                 delrho(n,2,1)*delrho(n,2,1) +
     &                 delrho(n,3,1)*delrho(n,3,1)
            dsqgamma = dsqrt(gammaval)
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
               call xc_vwn_5(tol_rho,facl,rhoval,dp81drho,
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
     ,                                  .false.,dummy)
            endif
c
            if (abs(cfac(7)).gt.eps)then
                 facl = cfac(7)
               call xc_vwn_1_rpa(tol_rho,facl,rhoval,
     &                           dp81drho,lnq,ipol,egas,lqwght,
     &                           .false.,dummy)
            endif
c
            if (abs(cfac(8)).gt.eps)then
                 facl = cfac(8)
               call xc_vwn_1(tol_rho,facl,rhoval,dp81drho,
     &                       lnq,ipol,egas,lqwght,.false.,dummy)
            endif
c
            if (abs(cfac(9)).gt.eps)then
                 facl = cfac(9)
               call xc_vwn_2(tol_rho,facl,rhoval,dp81drho,
     &                       lnq,ipol,egas,lqwght,.false.,dummy)
            endif
c
            if (abs(cfac(10)).gt.eps)then
                 facl = cfac(10)
               call xc_vwn_3(tol_rho,facl,rhoval,dp81drho,
     &                       lnq,ipol,egas,lqwght,.false.,dummy)
            endif
c
            if (abs(cfac(11)).gt.eps)then
                 facl = cfac(11)
               call xc_vwn_4(tol_rho,facl,rhoval,dp81drho,
     &                       lnq,ipol,egas,lqwght,.false.,dummy)
            endif
c
            degasda = arho*(dp81drho(1) - egas*arho)
            akf = TPP*rho13
            aks = sqrt(4.d0*akf/pi)
            dksdrho = aks*arho/6.d0
            t = 0.5d0*dsqgamma*arho/aks
            t2 = t*t
            t3 = t2*t
            t4 = t3*t
            dtdrhoa = -t*arho*SEV6
            dtdgamma = 0.5d0*t/gammaval
            dexpo1 = twoab2*(egas*arho)
            if(abs(egas).lt.tol_rho*tol_rho) goto 10
            if(abs(dexpo1).gt.tol_rho) then
            if (dexpo1.gt.EXPTOL)then
               expo1 = 0.d0
            else
               expo1 = exp(-dexpo1)
            endif
            A = twoab/(expo1-ONE)
            else
            A = 0d0
            endif
            A1p2At2 = ONE+2.d0*A*t2
            aden1 = ONE/(ONE+A*t2+A*A*t4)
            alognd = (t2+A*t4)*aden1
            if (twoab*alognd.gt.TOLL)then
               H0 = atwob2*log(ONE+twoab*alognd)
            else
               H0 = 0.d0
            endif
c
            func = H0 
c
            if (lfac)then
               Ec = Ec + egas*qwght(n)*fac
               if(ldew)ffunc(n)=ffunc(n)+egas*fac
            endif
            if (nlfac)then
               Ec = Ec + (func*rhoval)*qwght(n)*fac
               if(ldew)ffunc(n)=ffunc(n)+(func*rhoval)*fac
            endif
c
c           df/drho = func + rhoval*dfunc/drho
c
            dAdrhoa = A/(expo1-ONE)*expo1*twoab2*degasda
            dH0drhoa = bet*aden1*A1p2At2*(2.d0*t*dtdrhoa + 
     &                 dAdrhoa*t4/A1p2At2 - alognd*(2.d0*A*t*dtdrhoa +
     &                 dAdrhoa*t2))/(ONE+twoab*alognd)
            if (lfac)then
               Amat(n,1) = Amat(n,1) +  dp81drho(1)*fac
            endif
            if (nlfac)then
               Amat(n,1) = Amat(n,1) +  (func +
     &                     rhoval*dH0drhoa)*fac
            endif
c
c           df/dgamma
c
            if (dsqgamma.gt.TOLL)then
               dH0dgamma = bet/(ONE + twoab*alognd)*
     &                     aden1*2.d0*t*dtdgamma*A1p2At2*
     &                     (ONE - alognd*A)
               gaa = dH0dgamma*rhoval*fac*2.d0
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
            rho13  = abs(rhoval)**ONE3
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
               call xc_vwn_5(tol_rho,facl,rho_t,
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
               call xc_vwn_1_rpa(tol_rho,facl,rho_t,
     &                           dp81drho,lnq,ipol,egas,lqwght,
     &                           .false.,dummy)
            endif
c
            if (abs(cfac(8)).gt.eps)then
                 facl = cfac(8)
               call xc_vwn_1(tol_rho,facl,rho_t,
     &                       dp81drho,lnq,ipol,egas,lqwght,
     &                       .false.,dummy)
            endif
c
            if (abs(cfac(9)).gt.eps)then
                 facl = cfac(9)
               call xc_vwn_2(tol_rho,facl,rho_t,
     &                       dp81drho,lnq,ipol,egas,lqwght,
     &                       .false.,dummy)
            endif
c
            if (abs(cfac(10)).gt.eps)then
                 facl = cfac(10)
               call xc_vwn_3(tol_rho,facl,rho_t,
     &                       dp81drho,lnq,ipol,egas,lqwght,
     &                       .false.,dummy)
            endif
c
            if (abs(cfac(11)).gt.eps)then
                 facl = cfac(11)
               call xc_vwn_4(tol_rho,facl,rho_t,
     &                       dp81drho,lnq,ipol,egas,lqwght,
     &                       .false.,dummy)
            endif
            if(abs(egas).lt.tol_rho*tol_rho) goto 20
c
            degasda = arho*(dp81drho(1) - egas*arho)
            degasdb = arho*(dp81drho(2) - egas*arho)
            g = max(tol_rho,gzeta(zeta))
            g2 = g*g
            g3 = g2*g
            akf = TPP*rho13
            aks = sqrt(4.d0*akf/pi)
            dksdrho = aks*arho/6.d0
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
c
            expo1 = TOLL
            if (abs(dexpo1).lt.EXPTOL)expo1 = exp(-dexpo1)
c
            A = twoab/(expo1-ONE)
            A1p2At2 = ONE+2.d0*A*t2
            aden1 = ONE/(ONE+A*t2+A*A*t4)
            alognd = (t2+A*t4)*aden1
            H0 = atwob2*g3*log(ONE+twoab*alognd)
c
            func = H0 
c
            if (lfac)then
               Ec = Ec + egas*qwght(n)*fac
               if(ldew)ffunc(n)=ffunc(n)+egas*fac               
            endif
            if (nlfac)then
               Ec = Ec + (func*rhoval)*qwght(n)*fac
               if(ldew)ffunc(n)=ffunc(n)+(func*rhoval)*fac
            endif
c
c           df/drho = func + rhoval*dfunc/drho
c
            dAdrhoa = A/(expo1-ONE)*expo1*twoab2 /
     &                  g3*(-THREE*egas*arho*dgzetadag + degasda)
            dAdrhob = A/(expo1-ONE)*expo1*twoab2 /
     &                  g3*(-THREE*egas*arho*dgzetadbg + degasdb)
            dH0drhoa = Three*H0*dgzetadag + g3*bet*aden1*
     &                 (2.d0*t*dtdrhoa*A1p2At2 + dAdrhoa*t4 -
     &                  alognd*A1p2At2*(2.d0*A*t*dtdrhoa +
     &                  dAdrhoa*t2))/(ONE+twoab*alognd)
            dH0drhob = Three*H0*dgzetadbg + g3*bet*aden1*
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
     &                     rhoval*dH0drhoa)*fac
               Amat(n,2) = Amat(n,2) +  (func +
     &                     rhoval*dH0drhob)*fac
            endif
c
            if (dsqgamma.gt.TOLL)then
               dH0dgamma = g3*bet/(ONE + twoab*alognd)*
     &                     aden1*2.d0*t*dtdgamma*A1p2At2*
     &                     (ONE - alognd*A)
               gaa = 2.d0*dH0dgamma*rhoval*fac
               gbb = gaa
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

