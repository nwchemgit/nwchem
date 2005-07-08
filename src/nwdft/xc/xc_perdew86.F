#ifndef SECOND_DERIV
      Subroutine xc_perdew86(tol_rho, fac, lfac, nlfac, rho, delrho, 
     &                       Amat, Cmat, nq, ipol, Ec, qwght,
     ,     ldew, ffunc)
#else
      Subroutine xc_perdew86_d2(tol_rho, fac, lfac, nlfac, rho, delrho, 
     &                       Amat, Amat2, Cmat, Cmat2, nq, ipol, Ec,
     ,     qwght, ldew, ffunc)
#endif
c
c$Id: xc_perdew86.F,v 1.23 2005-07-08 00:25:21 edo Exp $
c
      implicit none
c
#include "dft2drv.fh"
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
#ifdef SECOND_DERIV
      double precision Amat2(nq,NCOL_AMAT2), Cmat2(nq,NCOL_CMAT2)
#endif
      double precision TOLL, EXPTOL, alpha, beta, pgamma, delta, 
     &                 beta10, ftilde, zzz, fff, pfff, CINF, ONE, 
     &                 ONE3, THREE, FOUR3, SEV6, FIVE3, 
     &                 TWO3, FIVE6, pi
      Parameter (TOLL = 1.D-40, EXPTOL = 80.d0)
      Parameter (alpha = 0.023266D0, beta  =  7.389D-6, 
     &   pgamma = 8.723d0, delta = 0.472d0,  beta10 = 10000.d0*beta)
      parameter (ftilde = 0.11d0, zzz = 0.001667d0, fff = 0.002568d0)
      parameter(pfff = 1.745d0, CINF = zzz+fff)
      Parameter (ONE = 1.D0, ONE3 = 1.d0/3.d0, THREE = 3.d0)
      Parameter (FOUR3 = 4.D0/3.D0, SEV6 = 7.d0/6.d0)
      parameter (FIVE3 = 5.d0/3.d0, TWO3 = 2.d0/3.d0, FIVE6 = 5.d0/6.d0)
      parameter (pi = 3.1415926535897932385d0)
c     
c     Mlynarski Salahub PRB 43, 1399 (1991)
c     
      integer n
      double precision rsfact, rs, rs2, rs3
      double precision rhoval, rho13, rho43, rho76, arho
      double precision d1rs
#ifdef SECOND_DERIV
      double precision d2rs
#endif
      double precision gamma, gam12
      double precision anum, aden, d1anum, d1aden, Cn, d1Cn,
     &     expfac, phi, d1phi(2), dlnphi, func, d1f(3),
     &     dlnfrho(2), dlnfgam
      double precision zeta, d1z(2), d, dm1, adp, d1d(2), t,
     &     dt12, d1dt12
#ifdef SECOND_DERIV
      double precision d2anum, d2aden, rrho2, d2z(3), dpp, d2d(3),
     &     d2phi(3), d2dt12, d2Cn
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
            if (abs(rhoval).lt.tol_rho) goto 10
            arho=1.d0/rhoval
            rho13 = abs(rhoval)**ONE3
            rho43 = rhoval*rho13
            rho76 = abs(rhoval)**SEV6
            rs = rsfact/rho13
            rs2 = rs*rs
            rs3 = rs2*rs
            d1rs = -ONE3*rs*arho
#ifdef SECOND_DERIV
            d2rs = -FOUR3*d1rs*arho
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
#ifdef SECOND_DERIV
               d2phi(2) = d1phi(2)*dlnphi
               d2phi(3) =-d1phi(2)/(2d0*gamma)
c!!! Which of the following are actually needed for restricted?
c!!! Should treat derivatives of d as zero? d is a constant?
               t = -d2phi(1) - (d1Cn/Cn)**2 + d2Cn/Cn + FOUR3/rhoval**2
               Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA)
     &              + (d1f(1)*dlnfrho(1)
     &              + func*t)*fac
               Amat2(n,D2_RA_RB) = Amat2(n,D2_RA_RB)
     &              + (d1f(1)*dlnfrho(1)
     &              + func*t)*fac
#if 0
               Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA)
     &              + (d1f(1)*dlnfrho(1)
     &              + func*(d1d(1)*d1d(1)*dm1**2-d2d(1)*dm1+t))*fac
               Amat2(n,D2_RA_RB) = Amat2(n,D2_RA_RB)
     &              + (d1f(1)*dlnfrho(2)
     &              + func*(d1d(1)*d1d(2)*dm1**2-d2d(2)*dm1+t))*fac
#endif
c
               t = (d1f(1)*dlnfgam - func*d2phi(2))*fac
               Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA) + t
               Cmat2(n,D2_RA_GAB) = Cmat2(n,D2_RA_GAB) + t*2d0
               Cmat2(n,D2_RA_GBB) = Cmat2(n,D2_RA_GBB) + t
c
               t = (d1f(3)*dlnfgam - func*(1d0/gamma**2+d2phi(3)))*fac
               Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA) + t
               Cmat2(n,D2_GAA_GBB) = Cmat2(n,D2_GAA_GBB) + t
               Cmat2(n,D2_GAA_GAB) = Cmat2(n,D2_GAA_GAB) + t*2d0
               Cmat2(n,D2_GAB_GAB) = Cmat2(n,D2_GAB_GAB) + t*4d0
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
            if (abs(rhoval).lt.tol_rho) goto 20
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
            d2dt12 = FIVE3*TWO3*0.25d0*(
     &         ((ONE+zeta)*.5d0)**(-ONE3) + ((ONE-zeta)*.5d0)**(-ONE3) )
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

#ifndef SECOND_DERIV
      Subroutine xc_p81(tol_rho, fac, lfac, nlfac, rho, Amat, nq, ipol, 
     &                  Ec, qwght, ldew, func)
#else
#include "dft2drv.fh"
      Subroutine xc_p81_d2(tol_rho, fac, lfac, nlfac, rho, Amat, Amat2,
     &                     nq, ipol, Ec, qwght, ldew, func)
#endif
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
      double precision Amat2(nq,*)
#endif
      double precision A(2), B(2), C(2), D(2), G(2), B1(2), B2(2),
     &                 pi, tol_rho, ONE3, FOUR3, TWO3
      save A, B, C, D, G, B1, B2
      parameter (pi = 3.1415926535897932385d0)
      Parameter (ONE3 = 1.d0/3.d0, FOUR3 = 4.D0/3.D0)
      Parameter (TWO3 = 2.d0/3.d0)
      integer n, i
      double precision rhoval, rs, alnrs, d1rs, e(2), d1e(2), rden(2),
     &                 d1den(2), d1zeta(2), d1ersz(2), d1edrho(2), eps,
     &                 sqrtrs, fz, d1fz, zeta
#ifdef SECOND_DERIV
      double precision d2rs, d2e(2), d2den(2), d2zeta(3), d2ersz(3),
     &                 d2edrho(3), d2fzeta, d2fz, rrho2
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
         if (abs(rho(n,1)).gt.tol_rho)then
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
            d2fz = d2fzeta(zeta)
            rrho2 = 2.d0/(rhoval*rhoval)
c           1 = aa, 2 = ab, 3 = bb
            d2zeta(1) =-rrho2*(1.d0-zeta)
            d2zeta(2) = rrho2*zeta
            d2zeta(3) = rrho2*(1.d0+zeta)
#endif
            if (rs.lt.1.d0)then
               alnrs = log(rs)
               do i = 1, 2
                  e(i) = A(i)*alnrs+B(i)+C(i)*rs*alnrs+D(i)*rs
                  d1e(i) = A(i)/rs+C(i)*(alnrs+1d0)+D(i)
#ifdef SECOND_DERIV
                  d2e(i) = (C(i)-A(i)/rs)/rs
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
