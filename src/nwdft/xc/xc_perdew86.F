#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      Subroutine xc_perdew86(tol_rho, fac, lfac, nlfac, rho, delrho, 
     &                       Amat, Cmat, nq, ipol, Ec, qwght,
     ,     ldew, ffunc)
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      Subroutine xc_perdew86_d2(tol_rho, fac, lfac, nlfac, rho, delrho, 
     &                       Amat, Amat2, Cmat, Cmat2, nq, ipol, Ec,
     ,     qwght, ldew, ffunc)
#else
      Subroutine xc_perdew86_d3(tol_rho, fac, lfac, nlfac, rho, delrho, 
     1     Amat, Amat2, Amat3, Cmat, Cmat2, Cmat3, nq, ipol, Ec,
     2     qwght, ldew, ffunc)
#endif
c
c$Id$
c
      implicit none
c
#include "dft2drv.fh"
#include "dft3drv.fh"
c
      double precision tol_rho, fac ! [input]
      integer nq, ipol              ! [input]
      double precision Ec           ! [input/output]
      logical lfac, nlfac, ldew
      double precision ffunc(*)  ! value of the functional [output]
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
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
      double precision Amat2(nq,NCOL_AMAT2), Cmat2(nq,NCOL_CMAT2)
#endif
#ifdef THIRD_DERIV
      double precision Amat3(nq,NCOL_AMAT3), Cmat3(nq,NCOL_CMAT3)
#endif
      double precision TOLL, EXPTOL, alpha, beta, pgamma, delta, 
     &                 beta10, ftilde, zzz, fff, pfff, CINF, ONE, 
     &                 ONE3, THREE, FOUR3, SEV6, FIVE3, 
     &                 TWO3, FIVE6, pi
      double precision SEVEN3, EIGHT3
      Parameter (TOLL = 1.D-40, EXPTOL = 80.d0)
      Parameter (alpha = 0.023266D0, beta  =  7.389D-6, 
     &   pgamma = 8.723d0, delta = 0.472d0,  beta10 = 10000.d0*beta)
      parameter (ftilde = 0.11d0, zzz = 0.001667d0, fff = 0.002568d0)
      parameter(pfff = 1.745d0, CINF = zzz+fff)
      Parameter (ONE = 1.D0, ONE3 = 1.d0/3.d0, THREE = 3.d0)
      Parameter (FOUR3 = 4.D0/3.D0, SEV6 = 7.d0/6.d0)
      parameter (FIVE3 = 5.d0/3.d0, TWO3 = 2.d0/3.d0, FIVE6 = 5.d0/6.d0)
      parameter (pi = 3.1415926535897932385d0)
      parameter (SEVEN3 = 7.0d0/3.0d0, EIGHT3 = 8.0d0/3.0d0)
c     
c     Mlynarski Salahub PRB 43, 1399 (1991)
c     
      integer n
      double precision rsfact, rs, rs2, rs3
      double precision rhoval, rho13, rho43, rho76, arho
      double precision d1rs
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
      double precision d2rs
#endif
#ifdef THIRD_DERIV
      double precision d3rs
#endif
      double precision gamma, gam12
      double precision anum, aden, d1anum, d1aden, Cn, d1Cn,
     &     expfac, phi, d1phi(2), dlnphi, func, d1f(3),
     &     dlnfrho(2), dlnfgam
      double precision zeta, d1z(2), d, dm1, adp, d1d(2), t,
     &     dt12, d1dt12
      double precision aden2
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
      double precision d2anum, d2aden, rrho2, d2z(3), dpp, d2d(3),
     &     d2phi(3), d2dt12, d2Cn
      double precision aden3
      double precision arho2
      double precision d2lnphi
      double precision d2f(3)
      double precision d2lnfrho(3), d2lnfrg(2), d2lnfgam
#endif
#ifdef THIRD_DERIV
      double precision d3lnphi
      double precision d3anum, d3aden, d3Cn, d3phi(4)
      double precision d3lnfrho(4), d3lnfgam
      double precision d3f(3)
      double precision aden4
      double precision arho3
#endif
c
      rsfact = (0.75d0/pi)**ONE3
c
      if (ipol.eq.1 )then
c
c        ======> SPIN-RESTRICTED <======
c
         do 10 n = 1, nq
            rhoval = rho(n,1)
            if (rhoval.lt.tol_rho) goto 10
            arho=1.d0/rhoval
            rho13 = abs(rhoval)**ONE3
            rho43 = rhoval*rho13
            rho76 = abs(rhoval)**SEV6
            rs = rsfact/rho13
            rs2 = rs*rs
            rs3 = rs2*rs
            d1rs = -ONE3*rs*arho
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
            d2rs = -FOUR3*d1rs*arho
#endif
#ifdef THIRD_DERIV
            d3rs = -SEVEN3*d2rs*arho 
#endif
            gamma = delrho(n,1,1)*delrho(n,1,1) +
     &              delrho(n,2,1)*delrho(n,2,1) +
     &              delrho(n,3,1)*delrho(n,3,1)
            gam12 = sqrt(abs(gamma))
c         
c           C(n)
c         
            anum = fff+alpha*rs+beta*rs2
            aden = 1.d0+pgamma*rs+delta*rs2+beta10*rs3
            Cn = zzz + anum/aden
            d1anum = alpha + 2d0*beta*rs
            d1aden = pgamma + 2d0*delta*rs + 3d0*beta10*rs2
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
            d2anum = 2d0*beta
            d2aden = 2d0*delta + 6d0*beta10*rs
#endif
#ifdef THIRD_DERIV
            d3anum = 0.0d0
            d3aden = 6.0d0*beta10 
#endif
c     First compute rs derivative
            aden2 = aden*aden
c
            d1Cn = d1anum/aden - anum*d1aden/aden2
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
            aden3 = aden2*aden
c
c            d2Cn = d2anum/aden - (2d0*d1anum*d1aden+anum*d2aden)/aden**2
c     &           + 2d0*anum*d1aden**2/aden**3
            d2Cn = d2anum/aden - (2d0*d1anum*d1aden+anum*d2aden)/aden2
     &           + 2d0*anum*d1aden**2/aden3
#endif
#ifdef THIRD_DERIV
            aden4 = aden3*aden
c
            d3Cn = -( 3.0d0*d2anum*d1aden + 3.0d0*d1anum*d2aden
     1              + anum*d3aden )/aden2
     2           + 6.0d0*( d1anum*d1aden**2
     3                   + anum*d2aden*d1aden )/aden3
     4           - 6.0d0*anum*d1aden**3/aden4
#endif
c     Convert to rho derivative
#ifdef THIRD_DERIV
            d3Cn = d3Cn*d1rs*d1rs*d1rs
     1           + 3.0d0*d2Cn*d2rs*d1rs
     2           + d1Cn*d3rs
#endif
c
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
            d2Cn = d2Cn*d1rs*d1rs + d1Cn*d2rs
#endif
            d1Cn = d1Cn*d1rs
c         
c           phi(n,gradn)
c         
            expfac = 0.d0
            phi = (pfff*ftilde)*(CINF/Cn)*gam12/rho76
            if (phi.lt.EXPTOL) expfac = exp(-phi)
            dlnphi = -(d1Cn/Cn + SEV6/rhoval)
            d1phi(1) = phi*dlnphi
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
            arho2 = arho*arho
c
            d2lnphi = (d1Cn/Cn)**2 - d2Cn/Cn + SEV6*arho2
c
c            d2phi(1) = d1phi(1)*dlnphi
c     &               + phi*((d1Cn/Cn)**2 - d2Cn/Cn + SEV6/rhoval**2)
            d2phi(1) = d1phi(1)*dlnphi + phi*d2lnphi
#endif
c
#ifdef THIRD_DERIV
            arho3 = arho2*arho
c
            d3lnphi = -2.0d0*(d1Cn/Cn)**3
     1              + 3.0d0*(d2Cn/Cn)*(d1Cn/Cn)
     2              - d3Cn/Cn
     3              - SEVEN3*arho3
            d3phi(1) = d2phi(1)*dlnphi
     1               + 2.0d0*d1phi(1)*d2lnphi
     2               + phi*d3lnphi
#endif
c         
c           functional
c         
            func = expfac*Cn*gamma/rho43
            dlnfrho(1) = d1Cn/Cn - (d1phi(1) + FOUR3/rhoval)
            d1f(1) = dlnfrho(1)*func
            Amat(n,1) = Amat(n,1) + d1f(1)*fac
            if (gam12.gt.TOLL)then
               d1phi(2) = phi / (2d0*gamma)
               dlnfgam = 1d0/gamma - d1phi(2)
               d1f(3) = func*dlnfgam
               Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + d1f(3)*fac
               Cmat(n,D1_GAB) = Cmat(n,D1_GAB) + d1f(3)*2D0*fac
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
               d2phi(2) = d1phi(2)*dlnphi
               d2phi(3) =-d1phi(2)/(2d0*gamma)
c rr terms
c!!! Which of the following are actually needed for restricted?
c!!! Should treat derivatives of d as zero? d is a constant?
c Daniel (11-19-12): d is a constant (it equals 1) for a restricted 
c calculation, since there is no spin-polarization.  Thus, the
c derivatives are zero.
               d2lnfrho(1) = -d2phi(1) - (d1Cn/Cn)**2 + d2Cn/Cn
     1                     + FOUR3*arho2
c
               d2f(1) = d1f(1)*dlnfrho(1)
     1                + func*d2lnfrho(1)
c
               t = d2f(1)*fac
c
               Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + t
               Amat2(n,D2_RA_RB) = Amat2(n,D2_RA_RB) + t
c OLD CODE
c               t = -d2phi(1) - (d1Cn/Cn)**2 + d2Cn/Cn + FOUR3/rhoval**2
c               Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA)
c     &              + (d1f(1)*dlnfrho(1)
c     &              + func*t)*fac
c               Amat2(n,D2_RA_RB) = Amat2(n,D2_RA_RB)
c     &              + (d1f(1)*dlnfrho(1)
c     &              + func*t)*fac
c OLD CODE
#if 0
               Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA)
     &              + (d1f(1)*dlnfrho(1)
     &              + func*(d1d(1)*d1d(1)*dm1**2-d2d(1)*dm1+t))*fac
               Amat2(n,D2_RA_RB) = Amat2(n,D2_RA_RB)
     &              + (d1f(1)*dlnfrho(2)
     &              + func*(d1d(1)*d1d(2)*dm1**2-d2d(2)*dm1+t))*fac
#endif
c rg terms
               d2lnfrg(1) = -d2phi(2)
               d2f(2) = (d1f(1)*dlnfgam + func*d2lnfrg(1))
               t = d2f(2)*fac
c
               Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA) + t
               Cmat2(n,D2_RA_GAB) = Cmat2(n,D2_RA_GAB) + t*2d0
               Cmat2(n,D2_RA_GBB) = Cmat2(n,D2_RA_GBB) + t
c OLD CODE
c               t = (d1f(1)*dlnfgam - func*d2phi(2))*fac
c               Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA) + t
c               Cmat2(n,D2_RA_GAB) = Cmat2(n,D2_RA_GAB) + t*2d0
c               Cmat2(n,D2_RA_GBB) = Cmat2(n,D2_RA_GBB) + t
c OLD CODE
c gg terms
               d2lnfgam = -1.0d0/gamma**2 - d2phi(3)
               d2f(3) = d1f(3)*dlnfgam + func*d2lnfgam
               t = d2f(3)*fac
c
               Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA) + t
               Cmat2(n,D2_GAA_GBB) = Cmat2(n,D2_GAA_GBB) + t
               Cmat2(n,D2_GAA_GAB) = Cmat2(n,D2_GAA_GAB) + t*2d0
               Cmat2(n,D2_GAB_GAB) = Cmat2(n,D2_GAB_GAB) + t*4d0
c OLD CODE
c               t = (d1f(3)*dlnfgam - func*(1d0/gamma**2+d2phi(3)))*fac
c               Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA) + t
c               Cmat2(n,D2_GAA_GBB) = Cmat2(n,D2_GAA_GBB) + t
c               Cmat2(n,D2_GAA_GAB) = Cmat2(n,D2_GAA_GAB) + t*2d0
c               Cmat2(n,D2_GAB_GAB) = Cmat2(n,D2_GAB_GAB) + t*4d0
c OLD CODE
#endif
c
#ifdef THIRD_DERIV
c rrr terms
               d3lnfrho(1) = -d3phi(1)
     1                     + 2.0d0*(d1Cn/Cn)**3
     2                     - 3.0d0*(d2Cn/Cn)*(d1Cn/Cn)
     3                     + d3Cn/Cn
     4                     - EIGHT3*arho3
c
               d3f(1) = d2f(1)*dlnfrho(1)
     1                + 2.0d0*d1f(1)*d2lnfrho(1)
     2                + func*d3lnfrho(1)
c
               t = d3f(1)*fac
c
               Amat3(n,D3_RA_RA_RA) = Amat3(n,D3_RA_RA_RA) + t
               Amat3(n,D3_RA_RA_RB) = Amat3(n,D3_RA_RA_RB) + t
               Amat3(n,D3_RA_RB_RB) = Amat3(n,D3_RA_RB_RB) + t
c rrg terms
               d3phi(2) = d2phi(2)*dlnphi + d1phi(2)*d2lnphi
c
               t = ( d2f(2)*dlnfrho(1)
     1             - d1f(1)*d2phi(2)
     2             + d1f(3)*d2lnfrho(1)
     3             - func*d3phi(2) )*fac
c
               Cmat3(n,D3_RA_RA_GAA) = Cmat3(n,D3_RA_RA_GAA) + t
               Cmat3(n,D3_RA_RA_GAB) = Cmat3(n,D3_RA_RA_GAB) + t*2.0d0
               Cmat3(n,D3_RA_RA_GBB) = Cmat3(n,D3_RA_RA_GBB) + t
               Cmat3(n,D3_RA_RB_GAA) = Cmat3(n,D3_RA_RB_GAA) + t
               Cmat3(n,D3_RA_RB_GAB) = Cmat3(n,D3_RA_RB_GAB) + t*2.0d0
               Cmat3(n,D3_RA_RB_GBB) = Cmat3(n,D3_RA_RB_GBB) + t
c rgg terms
               d3phi(3) = -d2phi(3)*dlnphi
c
               t = ( d2f(2)*dlnfgam
     1             + d1f(1)*d2lnfgam
     2             + d1f(3)*d2lnfrg(1)
     3             + func*d3phi(3) )*fac
c
               Cmat3(n,D3_RA_GAA_GAA) = Cmat3(n,D3_RA_GAA_GAA) + t
               Cmat3(n,D3_RA_GAA_GAB) = Cmat3(n,D3_RA_GAA_GAB) + t*2.0d0
               Cmat3(n,D3_RA_GAA_GBB) = Cmat3(n,D3_RA_GAA_GBB) + t
               Cmat3(n,D3_RA_GAB_GAB) = Cmat3(n,D3_RA_GAB_GAB) + t*4.0d0
               Cmat3(n,D3_RA_GAB_GBB) = Cmat3(n,D3_RA_GAB_GBB) + t*2.0d0
               Cmat3(n,D3_RA_GBB_GBB) = Cmat3(n,D3_RA_GBB_GBB) + t
c ggg terms              
               d3phi(4) = -3.0d0*d2phi(3)/(2.0d0*gamma)
               d3lnfgam = 2.0d0/gamma**3 - d3phi(4)
c
               t = ( d2f(3)*dlnfgam
     1             + 2.0d0*d1f(3)*d2lnfgam
     2             + func*d3lnfgam )*fac
c
               Cmat3(n,D3_GAA_GAA_GAA) = Cmat3(n,D3_GAA_GAA_GAA) + t
               Cmat3(n,D3_GAA_GAA_GAB) = Cmat3(n,D3_GAA_GAA_GAB) 
     1                                 + t*2.0d0
               Cmat3(n,D3_GAA_GAA_GBB) = Cmat3(n,D3_GAA_GAA_GBB) + t
               Cmat3(n,D3_GAA_GAB_GAB) = Cmat3(n,D3_GAA_GAB_GAB) 
     1                                 + t*4.0d0
               Cmat3(n,D3_GAA_GAB_GBB) = Cmat3(n,D3_GAA_GAB_GBB) 
     1                                 + t*2.0d0
               Cmat3(n,D3_GAA_GBB_GBB) = Cmat3(n,D3_GAA_GBB_GBB) + t
               Cmat3(n,D3_GAB_GAB_GAB) = Cmat3(n,D3_GAB_GAB_GAB) 
     1                                 + t*8.0d0
#endif
            endif
            Ec = Ec + func*qwght(n)*fac
            if (ldew) ffunc(n)=ffunc(n)+func*fac
   10    continue
      else
c
c        ======> SPIN-UNRESTRICTED <======
c
         do 20 n = 1, nq
            rhoval = rho(n,1)
            if (rhoval.lt.tol_rho) goto 20
            arho=1.d0/rhoval
            rho13  = abs(rhoval)**ONE3
            rho43  = rhoval*rho13
            rho76  = abs(rhoval)**SEV6
            rs = rsfact/rho13
            rs2 = rs*rs
            rs3 = rs2*rs
            d1rs = -ONE3*rs*arho
#ifdef SECOND_DERIV
            d2rs = -FOUR3*d1rs*arho
#endif
            gamma = delrho(n,1,1)*delrho(n,1,1) +
     &              delrho(n,2,1)*delrho(n,2,1) +
     &              delrho(n,3,1)*delrho(n,3,1) +
     &              delrho(n,1,2)*delrho(n,1,2) +
     &              delrho(n,2,2)*delrho(n,2,2) +
     &              delrho(n,3,2)*delrho(n,3,2) +
     &        2.d0*(delrho(n,1,1)*delrho(n,1,2) +
     &              delrho(n,2,1)*delrho(n,2,2) +
     &              delrho(n,3,1)*delrho(n,3,2))
            gam12 = sqrt(abs(gamma))
            zeta = (rho(n,2) - rho(n,3))*arho
            if(zeta.lt.-1d0) zeta=-1d0
            if(zeta.gt.1d0) zeta=1d0
            d1z(1) =  (1.d0 - zeta)*arho
            d1z(2) = -(1.d0 + zeta)*arho
#ifdef SECOND_DERIV
            rrho2 = 2.d0*arho*arho
c           1 = aa, 2 = ab, 3 = bb
            d2z(1) =-rrho2*(1.d0-zeta)
            d2z(2) = rrho2*zeta
            d2z(3) = rrho2*(1.d0+zeta)
#endif
c
c           d(zeta)
c
            dt12 = ((ONE+zeta)*.5d0)**FIVE3 + ((ONE-zeta)*.5d0)**FIVE3
            d1dt12 = FIVE3*0.5d0*(
     &           ((ONE+zeta)*.5d0)**TWO3 - ((ONE-zeta)*.5d0)**TWO3 )
            d = 2.d0**ONE3*dsqrt(dt12)
            dm1 = 1.d0/d
            adp = 0.5d0*d/dt12*d1dt12
            d1d(1) = adp*d1z(1)
            d1d(2) = adp*d1z(2)
#ifdef SECOND_DERIV
            if ((1.d0-zeta).lt.tol_rho) then
              d2dt12 = FIVE3*TWO3*0.25d0*(((ONE+zeta)*.5d0)**(-ONE3))
            else if ((1.d0+zeta).lt.tol_rho) then
              d2dt12 = FIVE3*TWO3*0.25d0*(((ONE-zeta)*.5d0)**(-ONE3))
            else
              d2dt12 = FIVE3*TWO3*0.25d0*(
     &         ((ONE+zeta)*.5d0)**(-ONE3) + ((ONE-zeta)*.5d0)**(-ONE3) )
            end if
c
            dpp =-0.5d0*adp/dt12*d1dt12
     &        + 2.d0**(-TWO3)*d2dt12/dsqrt(dt12)
            d2d(1) = dpp*d1z(1)*d1z(1) + adp*d2z(1)
            d2d(2) = dpp*d1z(1)*d1z(2) + adp*d2z(2)
            d2d(3) = dpp*d1z(2)*d1z(2) + adp*d2z(3)
#endif
c         
c           C(n)
c         
            anum = fff+alpha*rs+beta*rs2
            aden = 1.d0+pgamma*rs+delta*rs2+beta10*rs3
            Cn = zzz + anum/aden
            d1anum = alpha + 2d0*beta*rs
            d1aden = pgamma + 2d0*delta*rs + 3d0*beta10*rs2
#ifdef SECOND_DERIV
            d2anum = 2d0*beta
            d2aden = 2d0*delta + 6d0*beta10*rs
#endif
c     First compute rs derivative
            d1Cn = d1anum/aden - anum*d1aden/aden**2
#ifdef SECOND_DERIV
            d2Cn = d2anum/aden - (2d0*d1anum*d1aden+anum*d2aden)/aden**2
     &           + 2d0*anum*d1aden**2/aden**3
#endif
c     Convert to rho derivative
#ifdef SECOND_DERIV
            d2Cn = d2Cn*d1rs*d1rs + d1Cn*d2rs
#endif
            d1Cn = d1Cn*d1rs
c         
c           phi(n,gradn)
c         
            expfac = 0.d0
            phi = (pfff*ftilde)*(CINF/Cn)*gam12/rho76
            if (phi.lt.EXPTOL) expfac = exp(-phi)
            dlnphi = -(d1Cn/Cn + SEV6/rhoval)
            d1phi(1) = phi*dlnphi
#ifdef SECOND_DERIV
            d2phi(1) = d1phi(1)*dlnphi
     &               + phi*((d1Cn/Cn)**2 - d2Cn/Cn + SEV6/rhoval**2)
#endif
c         
c           functional
c         
            func = expfac*Cn*gamma/rho43*dm1
            t = d1Cn/Cn - (d1phi(1) + FOUR3/rhoval)
            dlnfrho(1) = t - dm1*d1d(1)
            dlnfrho(2) = t - dm1*d1d(2)
            d1f(1) = dlnfrho(1)*func
            d1f(2) = dlnfrho(2)*func
            Amat(n,1) = Amat(n,1) + d1f(1)*fac
            Amat(n,2) = Amat(n,2) + d1f(2)*fac
            if (gam12.gt.TOLL)then
               d1phi(2) = phi / (2d0*gamma)
               dlnfgam = 1d0/gamma - d1phi(2)
               d1f(3) = func*dlnfgam
               Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + d1f(3)*fac
               Cmat(n,D1_GAB) = Cmat(n,D1_GAB) + d1f(3)*2D0*fac
               Cmat(n,D1_GBB) = Cmat(n,D1_GBB) + d1f(3)*fac
#ifdef SECOND_DERIV
               d2phi(2) = d1phi(2)*dlnphi
               d2phi(3) =-d1phi(2)/(2d0*gamma)
c
               t = -d2phi(1) - (d1Cn/Cn)**2 + d2Cn/Cn + FOUR3/rhoval**2
               Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA)
     &              + (d1f(1)*dlnfrho(1)
     &              + func*(d1d(1)*d1d(1)*dm1**2-d2d(1)*dm1+t))*fac
               Amat2(n,D2_RA_RB) = Amat2(n,D2_RA_RB)
     &              + (d1f(1)*dlnfrho(2)
     &              + func*(d1d(1)*d1d(2)*dm1**2-d2d(2)*dm1+t))*fac
               Amat2(n,D2_RB_RB) = Amat2(n,D2_RB_RB)
     &              + (d1f(2)*dlnfrho(2)
     &              + func*(d1d(2)*d1d(2)*dm1**2-d2d(3)*dm1+t))*fac
c
               t = (d1f(1)*dlnfgam - func*d2phi(2))*fac
               Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA) + t
               Cmat2(n,D2_RA_GAB) = Cmat2(n,D2_RA_GAB) + t*2d0
               Cmat2(n,D2_RA_GBB) = Cmat2(n,D2_RA_GBB) + t
               t = (d1f(2)*dlnfgam - func*d2phi(2))*fac
               Cmat2(n,D2_RB_GAA) = Cmat2(n,D2_RB_GAA) + t
               Cmat2(n,D2_RB_GAB) = Cmat2(n,D2_RB_GAB) + t*2d0
               Cmat2(n,D2_RB_GBB) = Cmat2(n,D2_RB_GBB) + t
c
               t = (d1f(3)*dlnfgam - func*(1d0/gamma**2+d2phi(3)))*fac
               Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA) + t
               Cmat2(n,D2_GAA_GBB) = Cmat2(n,D2_GAA_GBB) + t
               Cmat2(n,D2_GBB_GBB) = Cmat2(n,D2_GBB_GBB) + t
               Cmat2(n,D2_GAA_GAB) = Cmat2(n,D2_GAA_GAB) + t*2d0
               Cmat2(n,D2_GAB_GBB) = Cmat2(n,D2_GAB_GBB) + t*2d0
               Cmat2(n,D2_GAB_GAB) = Cmat2(n,D2_GAB_GAB) + t*4d0
#endif
            endif
            Ec = Ec + func*qwght(n)*fac
            if (ldew) ffunc(n)=ffunc(n)+func*fac
   20    continue
      endif
      return
      end

#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      Subroutine xc_p81(tol_rho, fac, lfac, nlfac, rho, Amat, nq, ipol, 
     &                  Ec, qwght, ldew, func)
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#include "dft2drv.fh"
      Subroutine xc_p81_d2(tol_rho, fac, lfac, nlfac, rho, Amat, Amat2,
     &                     nq, ipol, Ec, qwght, ldew, func)
#else 
#include "dft3drv.fh"
      Subroutine xc_p81_d3(tol_rho, fac, lfac, nlfac, rho, Amat, Amat2,
     &                     Amat3, nq, ipol, Ec, qwght, ldew, func)
#endif
c Daniel (4-2-13): Third derivatives aren't implemented for Perdew 81
c LDA yet.  The preprocessor stuff is here to keep the compiler from
c complaining.  We need this functional for Perdew 86 to work.
c
c     Ceperley Alder LDA from Perdew Zunger PRB 23, 5048 (1981)
c
      implicit none
c
      integer nq, ipol
      logical lfac, nlfac, ldew
      double precision func(*)  ! value of the functional [output]
      double precision Ec, fac
c
c     Charge Density 
c
      double precision rho(nq,ipol*(ipol+1)/2)
c
c     Quadrature Weights
c
      double precision qwght(nq)
c
c     Sampling Matrices for the XC Potential & Energy
c
      double precision Amat(nq,ipol)
#ifdef SECOND_DERIV
c      double precision Amat2(nq,*)
      double precision Amat2(nq,NCOL_AMAT2)
#endif
c
#ifdef THIRD_DERIV
      double precision Amat3(nq,NCOL_AMAT3)
#endif
      double precision A(2), B(2), C(2), D(2), G(2), B1(2), B2(2),
     &                 pi, tol_rho, ONE3, FOUR3, TWO3
      double precision FIVE3, SEVEN3
      save A, B, C, D, G, B1, B2
      parameter (pi = 3.1415926535897932385d0)
      Parameter (ONE3 = 1.d0/3.d0, FOUR3 = 4.D0/3.D0)
      Parameter (TWO3 = 2.d0/3.d0)
      Parameter (FIVE3 = 5.0d0/3.0d0, SEVEN3 = 7.0d0/3.0d0)
      integer n, i
      double precision rhoval, rs, alnrs, d1rs, e(2), d1e(2), rden(2),
     &                 d1den(2), d1zeta(2), d1ersz(2), d1edrho(2), eps,
     &                 sqrtrs, fz, d1fz, zeta
#ifdef SECOND_DERIV
      double precision d2rs, d2e(2), d2den(2), d2zeta(3), d2ersz(3),
     &                 d2edrho(3), d2fzeta, d2fz, rrho2
#endif
#ifdef THIRD_DERIV
      double precision d3rs, d3fz, rrho3, d3zeta(4), d3den(2), d3e(2),
     1                 d3ersz(4), d3edrho(4)
#endif
      double precision x, fzeta, d1fzeta, rsfact
      fzeta(x) = ((1.d0+x)**FOUR3 + 
     &            (1.d0-x)**FOUR3 - 2.d0) / (2.d0**FOUR3-2.d0)
      d1fzeta(x) = FOUR3*((1.d0+x)**ONE3 - 
     &                    (1.d0-x)**ONE3) / (2.d0**FOUR3-2.d0)
#ifdef SECOND_DERIV
      d2fzeta(x) = ONE3*FOUR3*((1.d0+x)**(-TWO3) + 
     &                    (1.d0-x)**(-TWO3)) / (2.d0**FOUR3-2.d0)
#endif
      data A / 0.0311d0, 0.01555d0 /
      data B / -0.048d0, -0.0269d0 /
      data C / 0.0020d0, 0.0007d0 /
      data D / -0.0116d0, -0.0048d0 /
      data G / -.1423d0, -.0843d0 /
      data B1 / 1.0529d0, 1.3981d0 /
      data B2 / 0.3334d0, 0.2611d0 /
c
      rsfact = (0.75d0/pi)**ONE3
c
c     ======> BOTH SPIN-RESTRICTED AND UNRESTRICTED <======
c
      do n = 1, nq
         if (rho(n,1).gt.tol_rho)then
            rhoval = rho(n,1)
            if (ipol.eq.1) then
               zeta = 0.0d0
               d1zeta(1) = 1.d0/rhoval
               d1zeta(2) =-1.d0/rhoval
               fz = 0d0
               d1fz = 0d0
            else
               zeta = (rho(n,2)-rho(n,3))/rhoval
               if(zeta.lt.-1d0) zeta=-1d0
               if(zeta.gt.1d0) zeta=1d0
               fz = fzeta(zeta)
               d1fz = d1fzeta(zeta)
               d1zeta(1) = (1.d0-zeta)/rhoval
               d1zeta(2) =-(1.d0+zeta)/rhoval
            endif
            rs = rsfact/abs(rhoval)**ONE3
            d1rs = -ONE3*rs/rhoval
#ifdef SECOND_DERIV
            d2rs = -FOUR3*d1rs/rhoval
            if ((1.d0-zeta).lt.tol_rho) then
              d2fz = (1.d0+zeta)**(-TWO3)
            else if ((1.d0+zeta).lt.tol_rho) then
              d2fz = (1.d0-zeta)**(-TWO3)
            else
              d2fz = (1.d0+zeta)**(-TWO3) + (1.d0-zeta)**(-TWO3)
            end if
            d2fz = d2fz*ONE3*FOUR3/(2.d0**FOUR3-2.d0)
c
            rrho2 = 2.d0/(rhoval*rhoval)
c           1 = aa, 2 = ab, 3 = bb
            d2zeta(1) =-rrho2*(1.d0-zeta)
            d2zeta(2) = rrho2*zeta
            d2zeta(3) = rrho2*(1.d0+zeta)
#endif
c
#ifdef THIRD_DERIV
            d3rs = -SEVEN3*d2rs/rhoval 
            if ((1.d0-zeta).lt.tol_rho) then
              d3fz = (1.d0+zeta)**(-FIVE3)
            else if ((1.d0+zeta).lt.tol_rho) then
              d3fz = (1.d0-zeta)**(-FIVE3)
            else
              d3fz = (1.d0+zeta)**(-FIVE3) + (1.d0-zeta)**(-FIVE3)
            end if
            d3fz = -d3fz*TWO3*ONE3*FOUR3/(2.d0**FOUR3-2.d0)
c
            rrho3 = 2.0d0/(rhoval*rhoval*rhoval)
c
c 1 = aaa, 2 = aab, 3 = abb, 4 = bbb
            d3zeta(1) = 3.0d0*rrho3*(1.0d0 - zeta)
            d3zeta(2) = rrho3*(1.0d0 - 3.0d0*zeta)
            d3zeta(3) = -rrho3*(1.0d0 + 3.0d0*zeta)
            d3zeta(4) = -3.0d0*rrho3*(1.0d0 + zeta)
#endif
            if (rs.lt.1.d0)then
               alnrs = log(rs)
               do i = 1, 2
                  e(i) = A(i)*alnrs+B(i)+C(i)*rs*alnrs+D(i)*rs
                  d1e(i) = A(i)/rs+C(i)*(alnrs+1d0)+D(i)
#ifdef SECOND_DERIV
                  d2e(i) = (C(i)-A(i)/rs)/rs
#endif
#ifdef THIRD_DERIV
                  d3e(i) = 2.0d0*A(i)/(rs*rs*rs)
     1                   - C(i)/(rs*rs) 
#endif
               enddo
            else
               sqrtrs = sqrt(rs)
               do i = 1, 2
                  rden(i) = 1.d0/(1.d0+B1(i)*sqrtrs+B2(i)*rs)
                  d1den(i) = B1(i)/(2.d0*sqrtrs)+B2(i)
                  e(i) = G(i)*rden(i)
                  d1e(i) = -G(i)*d1den(i)*rden(i)**2
#ifdef SECOND_DERIV
                  d2den(i) = -B1(i)/(4.d0*rs*sqrtrs)
                  d2e(i) = G(i)*rden(i)**2
     &                 *(2.d0*d1den(i)**2*rden(i)-d2den(i))
#endif
#ifdef THIRD_DERIV
                  d3den(i) = 3.0d0*B1(i)/(8.0d0*rs*rs*sqrtrs)
                  d3e(i) = G(i)*rden(i)*rden(i)*
     1                   ( 6.0d0*( d1den(i)*d2den(i)*rden(i)
     2                           - d1den(i)*d1den(i)*d1den(i)*
     3                             rden(i)*rden(i) ) 
     4                   - d3den(i) )
#endif
               enddo
            endif
            eps = e(1) + fz*(e(2)-e(1))
            d1ersz(1) = d1e(1) + fz*(d1e(2)-d1e(1))
            d1ersz(2) = d1fz*(e(2)-e(1))
            d1edrho(1) = d1ersz(1)*d1rs + d1ersz(2)*d1zeta(1)
            d1edrho(2) = d1ersz(1)*d1rs + d1ersz(2)*d1zeta(2)
            Ec = Ec + eps*qwght(n)*rhoval*fac
            if (ldew) func(n) = func(n) + eps*rhoval*fac
            Amat(n,1) = Amat(n,1) + (eps + rhoval*d1edrho(1))*fac
            if (ipol.eq.2)
     &      Amat(n,2) = Amat(n,2) + (eps + rhoval*d1edrho(2))*fac
#ifdef SECOND_DERIV
c           1 = rsrs, 2 = rsz, 3 = zz
            d2ersz(1) = d2e(1) + fz*(d2e(2)-d2e(1))
            d2ersz(2) = d1fz*(d1e(2)-d1e(1))
            d2ersz(3) = d2fz*(e(2)-e(1))
c           1 = aa, 2 = ab, 3 = bb
            d2edrho(1) = d2ersz(1)*d1rs*d1rs
     &                 + d2ersz(2)*d1rs*d1zeta(1)*2.d0
     &                 + d2ersz(3)*d1zeta(1)*d1zeta(1)
     &                 + d1ersz(1)*d2rs
     &                 + d1ersz(2)*d2zeta(1)
            d2edrho(2) = d2ersz(1)*d1rs*d1rs
     &                 + d2ersz(2)*d1rs*(d1zeta(1)+d1zeta(2))
     &                 + d2ersz(3)*d1zeta(1)*d1zeta(2)
     &                 + d1ersz(1)*d2rs
     &                 + d1ersz(2)*d2zeta(2)
            d2edrho(3) = d2ersz(1)*d1rs*d1rs
     &                 + d2ersz(2)*d1rs*d1zeta(2)*2.d0
     &                 + d2ersz(3)*d1zeta(2)*d1zeta(2)
     &                 + d1ersz(1)*d2rs
     &                 + d1ersz(2)*d2zeta(3)
            Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA)
     &           + (2.d0*d1edrho(1) + rhoval*d2edrho(1))*fac
            Amat2(n,D2_RA_RB) = Amat2(n,D2_RA_RB)
     &           + (d1edrho(1) + d1edrho(2) + rhoval*d2edrho(2))*fac
            if (ipol.eq.2)
     &      Amat2(n,D2_RB_RB) = Amat2(n,D2_RB_RB)
     &           + (2.d0*d1edrho(2) + rhoval*d2edrho(3))*fac
#endif
#ifdef THIRD_DERIV
c 1 = rsrsrs, 2 = rsrsz, 3 = rszz, 4 = zzz
            d3ersz(1) = d3e(1) + fz*(d3e(2)-d3e(1))
            d3ersz(2) = d1fz*(d2e(2)-d2e(1))
            d3ersz(3) = d2fz*(d1e(2)-d1e(1))
            d3ersz(4) = d3fz*(e(2)-e(1))
c
c 1 = aaa, 2 = aab, 3 = abb, 4 = bbb
            d3edrho(1) = d3ersz(1)*d1rs*d1rs*d1rs
     1                 + d2ersz(1)*d1rs*d2rs*3.0d0
     2                 + d3ersz(3)*d1rs*d1zeta(1)*d1zeta(1)*3.0d0
     3                 + d2ersz(2)*d1rs*d2zeta(1)*3.0d0
     4                 + d1ersz(1)*d3rs
     5                 + d2ersz(2)*d1zeta(1)*d2rs*3.0d0
     6                 + d3ersz(2)*d1zeta(1)*d1rs*d1rs*3.0d0
     7                 + d3ersz(4)*d1zeta(1)*d1zeta(1)*d1zeta(1)
     8                 + d2ersz(3)*d1zeta(1)*d2zeta(1)*3.0d0
     9                 + d1ersz(2)*d3zeta(1)
            d3edrho(2) = d3ersz(1)*d1rs*d1rs*d1rs
     1                 + d2ersz(1)*d1rs*d2rs*3.0d0
     2                 + d3ersz(3)*d1rs*(d1zeta(1)*d1zeta(1)
     3                                 + d1zeta(1)*d1zeta(2)*2.0d0)
     4                 + d2ersz(2)*d1rs*(d2zeta(2)*2.0d0
     5                                 + d2zeta(1))
     6                 + d1ersz(1)*d3rs
     7                 + d2ersz(2)*d2rs*(d1zeta(1)*2.0d0
     8                                 + d1zeta(2))
     9                 + d3ersz(2)*d1rs*d1rs*(d1zeta(2)
     A                                      + d1zeta(1)*2.0d0)
     B                 + d3ersz(4)*d1zeta(2)*d1zeta(1)*d1zeta(1)
     C                 + d2ersz(3)*(d1zeta(1)*d2zeta(2)*2.0d0
     D                            + d1zeta(2)*d2zeta(1))
     E                 + d1ersz(2)*d3zeta(2)
            d3edrho(3) = d3ersz(1)*d1rs*d1rs*d1rs
     1                 + d2ersz(1)*d1rs*d2rs*3.0d0
     2                 + d3ersz(3)*d1rs*(d1zeta(2)*d1zeta(2)
     3                                 + d1zeta(2)*d1zeta(1)*2.0d0)
     4                 + d2ersz(2)*d1rs*(d2zeta(2)*2.0d0
     5                                 + d2zeta(3))
     6                 + d1ersz(1)*d3rs
     7                 + d2ersz(2)*d2rs*(d1zeta(2)*2.0d0
     8                                 + d1zeta(1))
     9                 + d3ersz(2)*d1rs*d1rs*(d1zeta(1)
     A                                      + d1zeta(2)*2.0d0)
     B                 + d3ersz(4)*d1zeta(1)*d1zeta(2)*d1zeta(2)
     C                 + d2ersz(3)*(d1zeta(2)*d2zeta(2)*2.0d0
     D                            + d1zeta(1)*d2zeta(3))
     E                 + d1ersz(2)*d3zeta(3)
            d3edrho(4) = d3ersz(1)*d1rs*d1rs*d1rs
     1                 + d2ersz(1)*d1rs*d2rs*3.0d0
     2                 + d3ersz(3)*d1rs*d1zeta(2)*d1zeta(2)*3.0d0
     3                 + d2ersz(2)*d1rs*d2zeta(3)*3.0d0
     4                 + d1ersz(1)*d3rs
     5                 + d2ersz(2)*d1zeta(2)*d2rs*3.0d0
     6                 + d3ersz(2)*d1zeta(2)*d1rs*d1rs*3.0d0
     7                 + d3ersz(4)*d1zeta(2)*d1zeta(2)*d1zeta(2)
     8                 + d2ersz(3)*d1zeta(2)*d2zeta(3)*3.0d0
     9                 + d1ersz(2)*d3zeta(4)
c
            Amat3(n,D3_RA_RA_RA) = Amat3(n,D3_RA_RA_RA)
     1           + ( 3.0d0*d2edrho(1) + rhoval*d3edrho(1) )*fac
            Amat3(n,D3_RA_RA_RB) = Amat3(n,D3_RA_RA_RB)
     1           + ( d2edrho(1) + 2.0d0*d2edrho(2) 
     2             + rhoval*d3edrho(2) )*fac
            Amat3(n,D3_RA_RB_RB) = Amat3(n,D3_RA_RB_RB)
     1           + ( 2.0d0*d2edrho(2) + d2edrho(3) 
     2             + rhoval*d3edrho(3) )*fac
            if (ipol.eq.2)
     1      Amat3(n,D3_RB_RB_RB) = Amat3(n,D3_RB_RB_RB)
     2           + ( 3.0d0*d2edrho(3) + rhoval*d3edrho(4) )*fac
#endif
         endif
      enddo
      return
      end
c
#ifndef SECOND_DERIV
#define SECOND_DERIV
c
c     Compile source again for the 2nd derivative case
c
#include "xc_perdew86.F"
#endif
c
#ifndef THIRD_DERIV
#define THIRD_DERIV
c
c     Compile source again for the 3rd derivative case
c
#include "xc_perdew86.F"
#endif
