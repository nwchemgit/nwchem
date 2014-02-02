#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      Subroutine xc_lyp88(tol_rho, fac,  rho, delrho,
     &                    Amat, Cmat, nq, ipol, Ec, qwght, ldew, func)
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      Subroutine xc_lyp88_d2(tol_rho, fac,  rho, delrho,
     &                       Amat, Amat2, Cmat, Cmat2, nq, ipol, Ec,
     &                       qwght, ldew, func)
#else
      Subroutine xc_lyp88_d3(tol_rho, fac,  rho, delrho,
     &                       Amat, Amat2, Amat3, Cmat, Cmat2, Cmat3, 
     &                       nq, ipol, Ec, qwght, ldew, func)
#endif
c
C$Id$
c
      implicit none 
c
#include "dft2drv.fh"
#include "dft3drv.fh"
c
      double precision fac ! [input]
      integer nq
      integer ipol
      double precision Ec
      logical ldew
      double precision func(*)  ! value of the functional [output]
c
c     Charge Density & Its Cube Root
c
      double precision rho(nq,(ipol*(ipol+1))/2)
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
      double precision F13, F43, F113, F83, F53, F19, F79, P1, tol_rho,
     &                 A, B, C, D
c
      Parameter (F13 = 1.D0/3.D0, F43 = 4.D0*F13, F113 = 11.D0*F13,
     &           F83 = 8.D0*F13, F53 = 5.D0*F13, F19 = 1.D0/9.D0,
     &           F79 = 7.D0*F19)
#ifdef THIRD_DERIV
      double precision F23, F73, F223
      Parameter (F23 = 2.0D0*F13, F73 = 7.0D0*F13, F223 = 22.0*F13)
#endif
c
c     P1 = 2**(11/3)*(3/10)*(3*PI**2)**(2/3)
c
      Parameter (P1 = 0.3646239897876487D+02)
c
c     Colle-Salvetti Empirical Parameters
c
      Parameter (A = 0.04918D0)
      Parameter (B = 0.13200D0)
      Parameter (C = 0.25330D0)
      Parameter (D = 0.34900D0)
c
c     Compute the partial derivatives of the correlation functional of 
c     Lee, Yang and Parr.
c
c     References:
c
c     Colle & Salvetti, Theor. Chim. Acta 37, 329 (1975)
c     Lee, Yang & Parr, Phys. Rev. B 37, 785 (1988)
c     Miehlich, Savin, Stoll & Preuss, Chem. Phys. Lett. 157, 200 (1989)
c     Johnson, Gill & Pople, J. Chem. Phys. 98, 5612 (1993)
c
      integer n
      double precision c1, c2, ab, ratcd
      double precision rrho, rhoa, rhob, rrho2, rhoa2, rhob2, rhoab,
     &     rhoa53, rhob53, rhoa83, rhob83, rho2,
     &     h1, h2, h3, d1h1, d1h2, d1h3, om, d1om, de, d1de, de11, de47,
     &     t, d1t(2)
c
      double precision xrarb, d1xrarb(2)
      double precision tm_in, d1tm_in(2)
c
      double precision gaa, gab, gbb
      double precision f1, f2, d1f1(2), d1f2(2), f, d1f(5),
     &     d2fgaa(2), d2fgab(2), d2fgbb(2)
c
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
      double precision d2h1, d2h2, d2h3, d2om, d2de, rrho3, d2f1(3),
     &     d2f2(3), rhoa113, rhob113, d3fgaa(3), d3fgab(3), d3fgbb(3),
     &     d2t(3), d2f(3)
      double precision dt
      double precision d2xrarb(3)
      double precision d2tm_in(3)
#endif
#ifdef THIRD_DERIV
      double precision rrho4, rhoa23, d3h1, d3h2, d3h3, d3om, d3de,
     1     d3f1(4), d3f2(4), d4fgaa(4), d4fgab(4), d4fgbb(4), d3t(4),
     2     d3f(4)
      double precision rhob23
      double precision ddt
      double precision d3xrarb(4)
      double precision d3tm_in(4)
#endif
c
c     Coefficients of first two terms in LYP functional and other
c     commonly occurring factors
c
      c1 = -4.0d0*a
      c2 = -P1*a*b
      ab = a*b
      ratcd = c/d
c
      if (ipol.eq.1)then
c
c        ======> SPIN-RESTRICTED <======
c
         do 10 n = 1, nq
            if (rho(n,1).lt.tol_rho)goto 10
c            rrho = 1d0/rho(n,1)
            rrho = 1.0d0/rho(n,1)
            rhoa = 0.5d0*rho(n,1)
            rrho2 = rrho*rrho
c            rho2 = 1d0/rrho2
            rho2 = 1.0d0/rrho2
            rhoa2 = rhoa*rhoa
            rhoab = rhoa2
c            rhoa53 = abs(rhoa)**F53*sign(1d0,rhoa)
            rhoa53 = abs(rhoa)**F53*sign(1.0d0,rhoa)
c            rhoa83 = abs(rhoa)**F83*sign(1d0,rhoa)
            rhoa83 = abs(rhoa)**F83*sign(1.0d0,rhoa)
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
            rrho3 = rrho*rrho2
            rhoa113 = rhoa*rhoa83
#endif
c
#ifdef THIRD_DERIV
            rrho4 = rrho*rrho3
            rhoa23 = abs(rhoa)**F23*sign(1.0d0,rhoa) 
#endif
c
            h2 = d*abs(rrho)**F13
            d1h2 = -F13*h2*rrho
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
            d2h2 = -F43*d1h2*rrho
#endif
c
#ifdef THIRD_DERIV
            d3h2 = -F73*d2h2*rrho
#endif 
c
            h3 = ratcd*h2
            d1h3 = ratcd*d1h2
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
            d2h3 = ratcd*d2h2
#endif
c
#ifdef THIRD_DERIV
            d3h3 = ratcd*d3h2
#endif
c
c            h1 = 1d0/(1d0+h2)
            h1 = 1.0d0/(1.0d0+h2)
            d1h1 = -h1*h1*d1h2
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c            d2h1 = -(2d0*h1*d1h1*d1h2 + h1*h1*d2h2)
            d2h1 = -(2.0d0*h1*d1h1*d1h2 + h1*h1*d2h2)
#endif
c
#ifdef THIRD_DERIV
            d3h1 = -6.0d0*d1h1*d1h1*d1h2
     1           - 6.0d0*h1*d2h2*d1h1
     2           - h1*h1*d3h2
#endif
c
!            om = exp(-h3)*h1*rrho**F113
            om = exp(-h3)*h1*abs(rrho)**F113
            t = d1h3+h1*d1h2+F113*rrho
            d1om = -om*t
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
            dt = d2h3+d1h1*d1h2+h1*d2h2-F113*rrho2
c            d2om = -(d1om*t+om*(d2h3+d1h1*d1h2+h1*d2h2-F113*rrho2))
            d2om = -(d1om*t+om*dt)
#endif
c
#ifdef THIRD_DERIV
            ddt = d3h3 + d2h1*d1h2 + 2.0d0*d1h1*d2h2
     1          + h1*d3h2 + F223*rrho3
            d3om = -(ddt*om + 2.0d0*d1om*dt + d2om*t)
#endif
c
            de = h3+h1*h2
            d1de = d1h3 + d1h1*h2 + h1*d1h2
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c            d2de = d2h3 + d2h1*h2 + h1*d2h2 + 2d0*d1h1*d1h2
            d2de = d2h3 + d2h1*h2 + h1*d2h2 + 2.0d0*d1h1*d1h2
#endif
c
#ifdef THIRD_DERIV
            d3de = d3h3 + d3h1*h2 + 3.0d0*d2h1*d1h2
     1           + 3.0d0*d1h1*d2h2 + h1*d3h2
#endif
c
            f1 = h1*rhoab*rrho
            d1f1(1) = d1h1*rhoab*rrho - h1*rhoab*rrho2
            d1f1(1) = d1f1(1) + h1*rhoa*rrho
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c Daniel (10-30-12): I'd suggest rewriting the second derivatives of
c term 1 in a much simpler format.  This relies on the fact that 
c 2*rhoa (or 2*rhob) is rho for a restricted calculation.
c            d2f1(1) = d2h1*rhoab*rrho + 2d0*d1h1*(rhoa*rrho-rhoab*rrho2)
c     &           + 2d0*h1*(-rhoa*rrho2+rhoab*rrho3)
c            d2f1(2) = d2h1*rhoab*rrho + d1h1*(1d0-2d0*rhoab*rrho2)
c     &           + 2d0*h1*rhoab*rrho3
            d2f1(1) = d2h1*rhoab*rrho
     1              + d1h1*(1.0d0 - rhoa*rrho)
     2              + h1*(-rrho + rhoa*rrho2)
            d2f1(2) = d2h1*rhoab*rrho
     1              + d1h1*(1.0d0-rhoa*rrho)
     2              + h1*rhoa*rrho2
#endif
c
#ifdef THIRD_DERIV
            d3f1(1) = d3h1*rhoab*rrho
     1              + d2h1*( 1.0d0 - rhoab*rrho2)
     2              + d1h1*(-1.5d0*rrho)
     3              + h1*(1.5d0*rrho2)
            d3f1(2) = d3h1*rhoab*rrho
     1              + d2h1*(1.0d0 - rhoab*rrho2)
     2              + d1h1*(0.50d0*rrho)
     3              + h1*(-0.50d0*rrho2)
            d3f1(3) = d3f1(2)
            d3f1(4) = d3f1(1)
#endif
c
            xrarb = rhoab*(rhoa83+rhoa83)
            f2 = om*xrarb
c
            d1xrarb(1) = rhoa*(F113*rhoa83+rhoa83)
c
c            d1f2(1) = d1om*rhoab*(rhoa83+rhoa83)
c            d1f2(1) = d1f2(1) + om*rhoa*(F113*rhoa83+rhoa83)
            d1f2(1) = d1om*xrarb + om*d1xrarb(1)
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c            d2f2(1) = d2om*rhoab*(rhoa83+rhoa83)
c     &          + 2d0*d1om*rhoa*(F113*rhoa83+rhoa83)
c     &          +       om*rhoa*F113*F83*rhoa53
c            d2f2(2) = d2om*rhoab*(rhoa83+rhoa83)
c     &           + d1om*(rhoa113+rhoa113+F113*(rhoa*rhoa83+rhoa*rhoa83))
c     &           +   om*F113*(rhoa83+rhoa83)
c
            d2xrarb(1) = rhoa*F113*F83*rhoa53
            d2xrarb(2) = F113*(rhoa83+rhoa83)
c
            d2f2(1) = d2om*xrarb
     &              + 2.0d0*d1om*d1xrarb(1)
     &              + om*d2xrarb(1)
            d2f2(2) = d2om*xrarb
     &              + 2.0d0*d1om*d1xrarb(1)
     &              + om*d2xrarb(2)
#endif
c
#ifdef THIRD_DERIV
            d3xrarb(1) = rhoa*F113*F83*F53*rhoa23
            d3xrarb(2) = F113*F83*rhoa53
            d3xrarb(3) = F113*F83*rhoa53
c
            d3f2(1) = d3om*xrarb
     1              + 3.0d0*d2om*d1xrarb(1)
     2              + 3.0d0*d1om*d2xrarb(1)
     3              + om*d3xrarb(1)
            d3f2(2) = d3om*xrarb
     1              + 3.0d0*d2om*d1xrarb(1)
     2              + d1om*(d2xrarb(1) + 2.0d0*d2xrarb(2))
     3              + om*d3xrarb(2)
            d3f2(3) = d3f2(2)
            d3f2(4) = d3f2(1)
#endif
c
            gaa =(delrho(n,1,1)*delrho(n,1,1) +
     &            delrho(n,2,1)*delrho(n,2,1) +
     &            delrho(n,3,1)*delrho(n,3,1))*0.25d0
c
c            de11 = de - 11d0
c            de47 = 47d0 - 7d0*de
            de11 = de - 11.0d0
            de47 = 47.0d0 - 7.0d0*de
c
c Daniel (10-30-12): tm_in is what I call Qi (the inside term)
            tm_in = 1.0d0 - 3.0d0*de - de11*rhoa*rrho
c Daniel (10-23-12): "t" is what I call Q or S.
c            t = F19*rhoab*(1d0-3d0*de-de11*rhoa*rrho)-rhoa2
            t = F19*rhoab*tm_in - rhoa2
c Daniel (10-30-12): Derivatives of the inside term
            d1tm_in(1) = -(3.0d0+rhoa*rrho)*d1de-de11*rhoa*rrho2
            d1tm_in(2) = -(3.0d0+rhoa*rrho)*d1de+de11*rhoa*rrho2
c Daniel (10-23-12): d1t(1) is the derivative with respect to rhoa,
c and d1t(2) is the derivative with respect to rhob.
c            d1t(1) = F19*(rhoa*(1d0-3d0*de-de11*rhoa*rrho)
c     &             - rhoab*((3d0+rhoa*rrho)*d1de+de11*rhoa*rrho2))
c            d1t(2) = F19*(rhoa*(1d0-3d0*de-de11*rhoa*rrho)
c     &             + rhoab*(-(3d0+rhoa*rrho)*d1de+de11*rhoa*rrho2))
c     &             - 2d0*rhoa
            d1t(1) = F19*( rhoa*tm_in + rhoab*d1tm_in(1) )
            d1t(2) = F19*( rhoa*tm_in + rhoab*d1tm_in(2) )
     &             - 2.0d0*rhoa
c
            d1f(3) = -ab*om*t
c
            d2fgaa(1) = -ab*( d1om*t + om*d1t(1) )
            d2fgaa(2) = -ab*( d1om*t + om*d1t(2) )
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c Daniel (10-30-12): Derivatives of the inside term, Qi.
            d2tm_in(1) = -(3.0d0+rhoa*rrho)*d2de
     1                 - 2.0d0*d1de*rhoa*rrho2
     2                 + 2.0d0*de11*rhoa*rrho3
            d2tm_in(2) = -(3.0d0+rhoa*rrho) ! Written without d2de
            d2tm_in(3) = -(3.0d0+rhoa*rrho)*d2de
     1                 + 2.0d0*d1de*rhoa*rrho2
     2                 - 2.0d0*de11*rhoa*rrho3
c
c            d2t(1) = -F19*(
c     &           2d0*rhoa*((3d0+rhoa*rrho)*d1de+de11*rhoa*rrho2)
c     &           +  rhoab*((3d0+rhoa*rrho)*d2de+2d0*d1de*rhoa*rrho2
c     &                                         -2d0*de11*rhoa*rrho3))
c            d2t(2) = F19*(1d0-3d0*de-de11*rhoa*rrho
c     &           - rho(n,1)*(3d0+rhoa*rrho)*d1de
c     &           - rhoab*((3d0+rhoa*rrho)*d2de))
c            d2t(3) = -F19*(
c     &           2d0*rhoa*((3d0+rhoa*rrho)*d1de-de11*rhoa*rrho2)
c     &           +  rhoab*((3d0+rhoa*rrho)*d2de-2d0*d1de*rhoa*rrho2
c     &                                         +2d0*de11*rhoa*rrho3))
c     &           - 2d0
            d2t(1) = F19*( 2.0d0*rhoa*d1tm_in(1)
     1                   + rhoab*d2tm_in(1) )
            d2t(2) = F19*( tm_in
     1                   + rho(n,1)*d2tm_in(2)*d1de
     2                   + rhoab*d2tm_in(2)*d2de )
            d2t(3) = F19*( 2.0d0*rhoa*d1tm_in(2)
     1                   + rhoab*d2tm_in(3) )
     2             - 2.0d0
c
c            d3fgaa(1) = -ab*(d2om*t+2d0*d1om*d1t(1)+om*d2t(1))
c            d3fgaa(2) = -ab*(d2om*t+d1om*(d1t(1)+d1t(2))+om*d2t(2))
c            d3fgaa(3) = -ab*(d2om*t+2d0*d1om*d1t(2)+om*d2t(3))
            d3fgaa(1) = -ab*( d2om*t + 2.0d0*d1om*d1t(1) + om*d2t(1) )
            d3fgaa(2) = -ab*( d2om*t + d1om*( d1t(1) + d1t(2) ) 
     1                      + om*d2t(2) )
            d3fgaa(3) = -ab*( d2om*t + 2.0d0*d1om*d1t(2) + om*d2t(3) )
#endif
c
#ifdef THIRD_DERIV
            d3tm_in(1) = -( 3.0d0 + rhoa*rrho )*d3de
     1                 - 3.0d0*d2de*rhoa*rrho2
     2                 + 6.0d0*d1de*rhoa*rrho3
     3                 - 6.0d0*de11*rhoa*rrho4
            d3tm_in(2) = -( 3.0d0 + rhoa*rrho )*d3de
     1                 - d2de*rhoa*rrho2
     2                 + 2.0d0*d1de*rhoa*rrho3
     3                 - 2.0d0*de11*rhoa*rrho4
            d3tm_in(3) = -( 3.0d0 + rhoa*rrho )*d3de
     1                 + d2de*rhoa*rrho2
     2                 - 2.0d0*d1de*rhoa*rrho3
     3                 + 2.0d0*de11*rhoa*rrho4
            d3tm_in(4) = -( 3.0d0 + rhoa*rrho )*d3de
     1                 + 3.0d0*d2de*rhoa*rrho2
     2                 - 6.0d0*d1de*rhoa*rrho3
     3                 + 6.0d0*de11*rhoa*rrho4
c
            d3t(1) = F19*( 3.0d0*rhoa*d2tm_in(1)
     1                   + rhoab*d3tm_in(1) )
            d3t(2) = F19*( 2.0d0*d1tm_in(1)
     1                   + rhoa*d2tm_in(1)
     2                   + 2.0d0*rhoa*d2tm_in(2)*d2de
     3                   + rhoab*d3tm_in(2) )
            d3t(3) = F19*( 2.0d0*d1tm_in(2)
     1                   + 2.0d0*rhoa*d2tm_in(2)*d2de
     2                   + rhoa*d2tm_in(3)
     3                   + rhoab*d3tm_in(3) )
            d3t(4) = F19*( 3.0d0*rhoa*d2tm_in(3)
     1                   + rhoab*d3tm_in(4) )
c
            d4fgaa(1) = -ab*( d3om*t + 3.0d0*d2om*d1t(1)
     1                      + 3.0d0*d1om*d2t(1) + om*d3t(1) )
            d4fgaa(2) = -ab*( d3om*t + d2om*(2.0d0*d1t(1) + d1t(2))
     1                      + d1om*(d2t(1) + 2.0d0*d2t(2)) + om*d3t(2) )
            d4fgaa(3) = -ab*( d3om*t + d2om*(d1t(1) + 2.0d0*d1t(2))
     1                      + d1om*(2.0d0*d2t(2) + d2t(3)) + om*d3t(3) )
            d4fgaa(4) = -ab*( d3om*t + 3.0d0*d2om*d1t(2)
     1                      + 3.0d0*d1om*d2t(3) + om*d3t(4) )
#endif
c
c Daniel (10-23-12): "t" is what I call R.
            t = F19*rhoab*de47-F43*rho2
c
            d1t(1) = F19*rhoa*de47 - F79*rhoab*d1de - F83*rho(n,1)
c
            d1f(4) = -ab*om*t
c
            d2fgab(1) = -ab*( d1om*t + om*d1t(1) )
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c            d2t(1) = -F79*(2d0*rhoa*d1de+rhoab*d2de) - F83
c            d2t(2) = F19*de47 - F79*(rho(n,1)*d1de+rhoab*d2de) - F83
            d2t(1) = -F79*( 2.0d0*rhoa*d1de + rhoab*d2de ) - F83
            d2t(2) = F19*de47 - F79*( rho(n,1)*d1de + rhoab*d2de ) 
     1             - F83
c
c            d3fgab(1) = -ab*(d2om*t+2d0*d1om*d1t(1)+om*d2t(1))
c            d3fgab(2) = -ab*(d2om*t+2d0*d1om*d1t(1)+om*d2t(2))
            d3fgab(1) = -ab*( d2om*t + 2.0d0*d1om*d1t(1) + om*d2t(1) )
            d3fgab(2) = -ab*( d2om*t + 2.0d0*d1om*d1t(1) + om*d2t(2) )
c To keep compilers quiet (WE MAY NEED TO FIX THIS)
            d3fgab(3) = d3fgab(1) 
#endif
c
#ifdef THIRD_DERIV
            d3t(1) = -F79*( 3.0d0*rhoa*d2de + rhoab*d3de )
            d3t(2) = -F79*( 2.0d0*d1de + 3.0d0*rhoa*d2de
     1                    + rhoab*d3de )
            d3t(3) = d3t(2)
            d3t(4) = d3t(1)
c
            d4fgab(1) = -ab*( d3om*t + 3.0d0*d2om*d1t(1)
     1                      + 3.0d0*d1om*d2t(1) + om*d3t(1) )
            d4fgab(2) = -ab*( d3om*t + 3.0d0*d2om*d1t(1)
     1                      + d1om*(d2t(1) + 2.0d0*d2t(2))
     2                      + om*d3t(2) )
            d4fgab(3) = -ab*( d3om*t + 3.0d0*d2om*d1t(1)
     1                      + d1om*(d2t(1) + 2.0d0*d2t(2))
     2                      + om*d3t(3) )
            d4fgab(4) = -ab*( d3om*t + 3.0d0*d2om*d1t(1)
     1                      + 3.0d0*d1om*d2t(1) + om*d3t(4) )
#endif
c
            d2fgbb(1) = d2fgaa(2)
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
            d3fgbb(1) = d3fgaa(3)
            d3fgbb(2) = d3fgaa(2)
            d3fgbb(3) = d3fgaa(1)
#endif
c
#ifdef THIRD_DERIV
            d4fgbb(1) = d4fgaa(4)
            d4fgbb(2) = d4fgaa(3)
            d4fgbb(3) = d4fgaa(2)
            d4fgbb(4) = d4fgaa(1)
#endif
c
            f = c1*f1 + c2*f2 + gaa*(2d0*d1f(3) + d1f(4))
            d1f(1) = c1*d1f1(1) + c2*d1f2(1)
     &             + gaa*(d2fgaa(1) + d2fgab(1) + d2fgbb(1))
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
            d2f(1) = c1*d2f1(1) + c2*d2f2(1)
     &             + gaa*(d3fgaa(1) + d3fgab(1) + d3fgbb(1))
            d2f(2) = c1*d2f1(2) + c2*d2f2(2)
     &             + gaa*(d3fgaa(2) + d3fgab(2) + d3fgbb(2))
#endif
c
#ifdef THIRD_DERIV
            d3f(1) = c1*d3f1(1) + c2*d3f2(1)
     1             + gaa*(d4fgaa(1) + d4fgab(1) + d4fgbb(1))
            d3f(2) = c1*d3f1(2) + c2*d3f2(2)
     1             + gaa*(d4fgaa(2) + d4fgab(2) + d4fgbb(2))
            d3f(3) = c1*d3f1(3) + c2*d3f2(3)
     1             + gaa*(d4fgaa(3) + d4fgab(3) + d4fgbb(3))
            d3f(4) = c1*d3f1(4) + c2*d3f2(4)
     1             + gaa*(d4fgaa(4) + d4fgab(4) + d4fgbb(4))
#endif
c
            Ec = Ec + f*fac*qwght(n)
            if (ldew) func(n) = func(n) + f*fac
            Amat(n,1) = Amat(n,1) + d1f(1)*fac
            Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + d1f(3)*fac
            Cmat(n,D1_GAB) = Cmat(n,D1_GAB) + d1f(4)*fac
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
            Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + d2f(1)*fac
            Amat2(n,D2_RA_RB) = Amat2(n,D2_RA_RB) + d2f(2)*fac
            Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA) + d2fgaa(1)*fac
            Cmat2(n,D2_RA_GAB) = Cmat2(n,D2_RA_GAB) + d2fgab(1)*fac
            Cmat2(n,D2_RA_GBB) = Cmat2(n,D2_RA_GBB) + d2fgbb(1)*fac
#endif
c
#ifdef THIRD_DERIV
            Amat3(n,D3_RA_RA_RA) = Amat3(n,D3_RA_RA_RA) + d3f(1)*fac
            Amat3(n,D3_RA_RA_RB) = Amat3(n,D3_RA_RA_RB) + d3f(2)*fac
            Amat3(n,D3_RA_RB_RB) = Amat3(n,D3_RA_RB_RB) + d3f(3)*fac
c
            Cmat3(n,D3_RA_RA_GAA) = Cmat3(n,D3_RA_RA_GAA) 
     1                            + d3fgaa(1)*fac
            Cmat3(n,D3_RA_RA_GAB) = Cmat3(n,D3_RA_RA_GAB) 
     1                            + d3fgab(1)*fac
            Cmat3(n,D3_RA_RA_GBB) = Cmat3(n,D3_RA_RA_GBB) 
     1                            + d3fgbb(1)*fac
c
            Cmat3(n,D3_RA_RB_GAA) = Cmat3(n,D3_RA_RB_GAA) 
     1                            + d3fgaa(2)*fac
            Cmat3(n,D3_RA_RB_GAB) = Cmat3(n,D3_RA_RB_GAB) 
     1                            + d3fgab(2)*fac
            Cmat3(n,D3_RA_RB_GBB) = Cmat3(n,D3_RA_RB_GBB) 
     1                            + d3fgbb(2)*fac
#endif
c
   10    continue
      else
c
c        ======> SPIN-UNRESTRICTED <======
c
         do 20 n = 1,nq
            if (rho(n,1).lt.tol_rho)goto 20
            rrho = 1d0/rho(n,1)
            rhoa = max(0.0d0,rho(n,2))
            rhob = max(0.0d0,rho(n,3))
            rrho2 = rrho*rrho
            rho2 = 1d0/rrho2
            rhoa2 = rhoa*rhoa
            rhob2 = rhob*rhob
            rhoab = rhoa*rhob
            rhoa53 = rhoa**F53
            rhob53 = rhob**F53
csedo            rhoa53 = abs(rhoa)**F53
csedo            rhob53 = abs(rhob)**F53
            rhoa83 = rhoa*rhoa53
            rhob83 = rhob*rhob53
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
            rrho3 = rrho*rrho2
            rhoa113 = rhoa*rhoa83
            rhob113 = rhob*rhob83
#endif
c
#ifdef THIRD_DERIV
            rrho4 = rrho*rrho3
            rhoa23 = rhoa**F23 
            rhob23 = rhob**F23 
#endif
c
cedo            h2 = d*abs(rrho)**F13*sign(1d0,rrho)
            h2 = d*rrho**F13
            d1h2 = -F13*h2*rrho
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
            d2h2 = -F43*d1h2*rrho
#endif
c
#ifdef THIRD_DERIV
            d3h2 = -F73*d2h2*rrho
#endif
c
            h3 = ratcd*h2
            d1h3 = ratcd*d1h2
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
            d2h3 = ratcd*d2h2
#endif
c
#ifdef THIRD_DERIV
            d3h3 = ratcd*d3h2
#endif
c
            h1 = 1d0/(1d0+h2)
            d1h1 = -h1*h1*d1h2
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
            d2h1 = -(2d0*h1*d1h1*d1h2 + h1*h1*d2h2)
#endif
c
#ifdef THIRD_DERIV
            d3h1 = -6.0d0*d1h1*d1h1*d1h2
     1           - 6.0d0*h1*d2h2*d1h1
     2           - h1*h1*d3h2
#endif
c
            om = exp(-h3)*h1*rrho**F113
cedo            om = exp(-h3)*h1*abs(rrho)**F113*sign(1d0,rrho)
            t = d1h3+h1*d1h2+F113*rrho
            d1om = -om*t
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
            dt = d2h3+d1h1*d1h2+h1*d2h2-F113*rrho2
c            d2om = -(d1om*t+om*(d2h3+d1h1*d1h2+h1*d2h2-F113*rrho2))
            d2om = -(d1om*t+om*dt)
#endif
c
#ifdef THIRD_DERIV
            ddt = d3h3 + d2h1*d1h2 + 2.0d0*d1h1*d2h2
     1          + h1*d3h2 + F223*rrho3
            d3om = -(ddt*om + 2.0d0*d1om*dt + d2om*t)
#endif
c
            de = h3+h1*h2
            d1de = d1h3 + d1h1*h2 + h1*d1h2
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
            d2de = d2h3 + d2h1*h2 + h1*d2h2 + 2d0*d1h1*d1h2
#endif
c
#ifdef THIRD_DERIV
            d3de = d3h3 + d3h1*h2 + 3.0d0*d2h1*d1h2
     1           + 3.0d0*d1h1*d2h2 + h1*d3h2
#endif
c
c Daniel (3-21-13): f1 is J/(-4*a) in my notes.
            f1 = h1*rhoab*rrho
            d1f1(1) = d1h1*rhoab*rrho - h1*rhoab*rrho2
            d1f1(2) = d1f1(1)
            d1f1(1) = d1f1(1) + h1*rhob*rrho
            d1f1(2) = d1f1(2) + h1*rhoa*rrho
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
            d2f1(1) = d2h1*rhoab*rrho + 2d0*d1h1*(rhob*rrho-rhoab*rrho2)
     &           + 2d0*h1*(-rhob*rrho2+rhoab*rrho3)
            d2f1(2) = d2h1*rhoab*rrho + d1h1*(1d0-2d0*rhoab*rrho2)
     &           + 2d0*h1*rhoab*rrho3
            d2f1(3) = d2h1*rhoab*rrho + 2d0*d1h1*(rhoa*rrho-rhoab*rrho2)
     &           + 2d0*h1*(-rhoa*rrho2+rhoab*rrho3)
#endif
c
#ifdef THIRD_DERIV
            d3f1(1) = d3h1*rhoab*rrho 
     1              + 3.0d0*d2h1*(rhob*rrho-rhoab*rrho2)
     2              + 6.0d0*d1h1*(-rhob*rrho2+rhoab*rrho3)
     3              + 6.0d0*h1*(rhob*rrho3-rhoab*rrho4)
            d3f1(2) = d3h1*rhoab*rrho
     1              + d2h1*( (2.0d0*rhob + rhoa)*rrho
     2                     - 3.0d0*rhoab*rrho2 )
     3              + 2.0d0*d1h1*( rrho 
     4                           - (2.0d0*rhob + rhoa)*rrho2
     5                           + 3.0d0*rhoab*rrho3 )
     6              + 2.0d0*h1*( -rrho2
     7                         + (2.0d0*rhob + rhoa)*rrho3
     8                         - 3.0d0*rhoab*rrho4 )
            d3f1(3) = d3h1*rhoab*rrho
     1              + d2h1*( (rhob + 2.0d0*rhoa)*rrho
     2                     - 3.0d0*rhoab*rrho2 )
     3              + 2.0d0*d1h1*( rrho
     4                           - (rhob + 2.0d0*rhoa)*rrho2
     5                           + 3.0d0*rhoab*rrho3 )
     6              + 2.0d0*h1*( -rrho2
     7                         + (rhob + 2.0d0*rhoa)*rrho3
     8                         - 3.0d0*rhoab*rrho4 )
            d3f1(4) = d3h1*rhoab*rrho
     1              + 3.0d0*d2h1*(rhoa*rrho-rhoab*rrho2)
     2              + 6.0d0*d1h1*(-rhoa*rrho2+rhoab*rrho3)
     3              + 6.0d0*h1*(rhoa*rrho3-rhoab*rrho4)
#endif
c
c Daniel (10-30-12): Define xrarb here
            xrarb = rhoab*(rhoa83+rhob83)
c
c            f2 =om*rhoab*(rhoa83+rhob83)
            f2 =om*xrarb
c
            d1xrarb(1) = rhob*(F113*rhoa83+rhob83)
            d1xrarb(2) = rhoa*(F113*rhob83+rhoa83)
c
c            d1f2(1) = d1om*rhoab*(rhoa83+rhob83)
            d1f2(1) = d1om*xrarb
            d1f2(2) = d1f2(1)
c            d1f2(1) = d1f2(1) + om*rhob*(F113*rhoa83+rhob83)
            d1f2(1) = d1f2(1) + om*d1xrarb(1)
c            d1f2(2) = d1f2(2) + om*rhoa*(F113*rhob83+rhoa83)
            d1f2(2) = d1f2(2) + om*d1xrarb(2)
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
            d2xrarb(1) = rhob*F113*F83*rhoa53
            d2xrarb(2) = F113*(rhoa83+rhob83)
            d2xrarb(3) = rhoa*F113*F83*rhob53
c
c            d2f2(1) = d2om*rhoab*(rhoa83+rhob83)
c     &          + 2d0*d1om*rhob*(F113*rhoa83+rhob83)
c     &          +       om*rhob*F113*F83*rhoa53
c            d2f2(2) = d2om*rhoab*(rhoa83+rhob83)
c     &           + d1om*(rhoa113+rhob113+F113*(rhob*rhoa83+rhoa*rhob83))
c     &           +   om*F113*(rhoa83+rhob83)
c            d2f2(3) = d2om*rhoab*(rhoa83+rhob83)
c     &          + 2d0*d1om*rhoa*(F113*rhob83+rhoa83)
c     &          +       om*rhoa*F113*F83*rhob53
            d2f2(1) = d2om*xrarb
     &          + 2d0*d1om*d1xrarb(1)
     &          +       om*d2xrarb(1)
            d2f2(2) = d2om*xrarb
     &           + d1om*(d1xrarb(1) + d1xrarb(2))
     &           +   om*d2xrarb(2)
            d2f2(3) = d2om*xrarb
     &          + 2d0*d1om*d1xrarb(2)
     &          +       om*d2xrarb(3)
#endif
c
#ifdef THIRD_DERIV
            d3xrarb(1) = rhob*F113*F83*F53*rhoa23
            d3xrarb(2) = F113*F83*rhoa53
            d3xrarb(3) = F113*F83*rhob53
            d3xrarb(4) = rhoa*F113*F83*F53*rhob23
c
            d3f2(1) = d3om*xrarb
     1              + 3.0d0*d2om*d1xrarb(1)
     2              + 3.0d0*d1om*d2xrarb(1)
     3              + om*d3xrarb(1)
            d3f2(2) = d3om*xrarb
     1              + d2om*( d1xrarb(2) + 2.0d0*d1xrarb(1) )
     2              + d1om*( 2.0d0*d2xrarb(2) + d2xrarb(1) )
     3              + om*d3xrarb(2) 
            d3f2(3) = d3om*xrarb
     1              + d2om*( 2.0d0*d1xrarb(2) + d1xrarb(1) )
     2              + d1om*( d2xrarb(3) + 2.0d0*d2xrarb(2) )
     3              + om*d3xrarb(3) 
            d3f2(4) = d3om*xrarb
     1              + 3.0d0*d2om*d1xrarb(2)
     2              + 3.0d0*d1om*d2xrarb(3)
     3              + om*d3xrarb(4) 
#endif
c
            gaa = delrho(n,1,1)*delrho(n,1,1) +
     &            delrho(n,2,1)*delrho(n,2,1) +
     &            delrho(n,3,1)*delrho(n,3,1)
            gab = delrho(n,1,1)*delrho(n,1,2) +
     &            delrho(n,2,1)*delrho(n,2,2) +
     &            delrho(n,3,1)*delrho(n,3,2)
            gbb = delrho(n,1,2)*delrho(n,1,2) +
     &            delrho(n,2,2)*delrho(n,2,2) +
     &            delrho(n,3,2)*delrho(n,3,2)
c
            de11 = de - 11d0
            de47 = 47d0 - 7d0*de
c Daniel (10-30-12): tm_in is what I call Qi (the inside term)
            tm_in = 1.0d0 - 3.0d0*de - de11*rhoa*rrho
c
c            t = F19*rhoab*(1d0-3d0*de-de11*rhoa*rrho)-rhob2
            t = F19*rhoab*tm_in-rhob2
c Daniel (10-30-12): Derivatives of the inside term
            d1tm_in(1) = -(3.0d0+rhoa*rrho)*d1de-de11*rhob*rrho2
            d1tm_in(2) = -(3.0d0+rhoa*rrho)*d1de+de11*rhoa*rrho2
c
c            d1t(1) = F19*(rhob*(1d0-3d0*de-de11*rhoa*rrho)
c     &             - rhoab*((3d0+rhoa*rrho)*d1de+de11*rhob*rrho2))
c            d1t(2) = F19*(rhoa*(1d0-3d0*de-de11*rhoa*rrho)
c     &             + rhoab*(-(3d0+rhoa*rrho)*d1de+de11*rhoa*rrho2))
c     &             - 2d0*rhob
            d1t(1) = F19*(rhob*tm_in + rhoab*d1tm_in(1))
            d1t(2) = F19*(rhoa*tm_in + rhoab*d1tm_in(2)) - 2.0d0*rhob
c
            d1f(3) = -ab*om*t
c
            d2fgaa(1) = -ab*( d1om*t + om*d1t(1) )
            d2fgaa(2) = -ab*( d1om*t + om*d1t(2) )
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c Daniel (10-30-12): Derivatives of the inside term, Qi.
            d2tm_in(1) = -(3.0d0+rhoa*rrho)*d2de
     1                 - 2.0d0*d1de*rhob*rrho2
     2                 + 2.0d0*de11*rhob*rrho3
            d2tm_in(2) = -(3.0d0+rhoa*rrho)*d2de
     1                 + ( rhoa - rhob)*d1de*rrho2
     2                 + ( rhob - rhoa)*de11*rrho3
            d2tm_in(3) = -(3.0d0+rhoa*rrho)*d2de
     1                 + 2.0d0*d1de*rhoa*rrho2
     2                 - 2.0d0*de11*rhoa*rrho3
c
c            d2t(1) = -F19*(
c     &           2d0*rhob*((3d0+rhoa*rrho)*d1de+de11*rhob*rrho2)
c     &           +  rhoab*((3d0+rhoa*rrho)*d2de+2d0*d1de*rhob*rrho2
c     &                                         -2d0*de11*rhob*rrho3))
c            d2t(2) = F19*(1d0-3d0*de-de11*rhoa*rrho
c     &           - rho(n,1)*(3d0+rhoa*rrho)*d1de
c     &           - rhoab*((3d0+rhoa*rrho)*d2de-d1de*(rhoa-rhob)*rrho2
c     &                                        +de11*(rhoa-rhob)*rrho3))
c            d2t(3) = -F19*(
c     &           2d0*rhoa*((3d0+rhoa*rrho)*d1de-de11*rhoa*rrho2)
c     &           +  rhoab*((3d0+rhoa*rrho)*d2de-2d0*d1de*rhoa*rrho2
c     &                                         +2d0*de11*rhoa*rrho3))
c     &           - 2d0
            d2t(1) = F19*( 2.0d0*rhob*d1tm_in(1)
     &                   + rhoab*d2tm_in(1) )
            d2t(2) = F19*( tm_in
     &                   - rho(n,1)*(3.0d0+rhoa*rrho)*d1de
     &                   + rhoab*d2tm_in(2) )
            d2t(3) = F19*( 2.0d0*rhoa*d1tm_in(2)
     &                   + rhoab*d2tm_in(3) )
     &             - 2.0d0
c
            d3fgaa(1) = -ab*( d2om*t + 2d0*d1om*d1t(1) + om*d2t(1) )
            d3fgaa(2) = -ab*( d2om*t + d1om*( d1t(1) + d1t(2) ) 
     1                      + om*d2t(2) )
            d3fgaa(3) = -ab*( d2om*t + 2d0*d1om*d1t(2) + om*d2t(3) )
#endif
c
#ifdef THIRD_DERIV
            d3tm_in(1) = -( 3.0d0 + rhoa*rrho )*d3de
     1                 - 3.0d0*d2de*rhob*rrho2
     2                 + 6.0d0*d1de*rhob*rrho3
     3                 - 6.0d0*de11*rhob*rrho4
            d3tm_in(2) = -( 3.0d0 + rhoa*rrho )*d3de
     1                 + d2de*rhoa*rrho2
     2                 - 2.0d0*d2de*rhob*rrho2
     3                 - 2.0d0*d1de*rhoa*rrho3
     4                 + 4.0d0*d1de*rhob*rrho3
     5                 - 2.0d0*de11*(2.0d0*rhob-rhoa)*rrho4
            d3tm_in(3) = -( 3.0d0 + rhoa*rrho )*d3de
     1                 + 2.0d0*d2de*rhoa*rrho2
     2                 - d2de*rhob*rrho2
     3                 - 4.0d0*d1de*rhoa*rrho3
     4                 + 2.0d0*d1de*rhob*rrho3
     5                 - 2.0d0*de11*(rhob-2.0d0*rhoa)*rrho4 
            d3tm_in(4) = -( 3.0d0 + rhoa*rrho )*d3de
     1                 + 3.0d0*d2de*rhoa*rrho2
     2                 - 6.0d0*d1de*rhoa*rrho3
     3                 + 6.0d0*de11*rhoa*rrho4
c
            d3t(1) = F19*( 3.0d0*rhob*d2tm_in(1)
     1                   + rhoab*d3tm_in(1) )
            d3t(2) = F19*( 2.0d0*d1tm_in(1)
     1                   + rhoa*d2tm_in(1)
     2                   + 2.0d0*rhob*d2tm_in(2)
     3                   + rhoab*d3tm_in(2) )
            d3t(3) = F19*( 2.0d0*d1tm_in(2)
     1                   + 2.0d0*rhoa*d2tm_in(2)
     2                   + rhob*d2tm_in(3)
     3                   + rhoab*d3tm_in(3) )
            d3t(4) = F19*( 3.0d0*rhoa*d2tm_in(3)
     1                   + rhoab*d3tm_in(4) )
c
            d4fgaa(1) = -ab*( d3om*t + 3.0d0*d2om*d1t(1)
     1                      + 3.0d0*d1om*d2t(1) + om*d3t(1) )
            d4fgaa(2) = -ab*( d3om*t + d2om*(2.0d0*d1t(1) + d1t(2))
     1                      + d1om*(d2t(1) + 2.0d0*d2t(2)) + om*d3t(2) )
            d4fgaa(3) = -ab*( d3om*t + d2om*(d1t(1) + 2.0d0*d1t(2))
     1                      + d1om*(2.0d0*d2t(2) + d2t(3)) + om*d3t(3) )
            d4fgaa(4) = -ab*( d3om*t + 3.0d0*d2om*d1t(2)
     1                      + 3.0d0*d1om*d2t(3) + om*d3t(4) )
#endif
c
            t = F19*rhoab*de47 - F43*rho2
            d1t(1) = F19*rhob*de47 - F79*rhoab*d1de - F83*rho(n,1)
            d1t(2) = F19*rhoa*de47 - F79*rhoab*d1de - F83*rho(n,1)
            d1f(4) = -ab*om*t
            d2fgab(1) = -ab*( d1om*t + om*d1t(1) )
            d2fgab(2) = -ab*( d1om*t + om*d1t(2) )
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
            d2t(1) = -F79*( 2d0*rhob*d1de + rhoab*d2de ) - F83
            d2t(2) = F19*de47 - F79*( rho(n,1)*d1de + rhoab*d2de ) 
     1             - F83
            d2t(3) = -F79*( 2d0*rhoa*d1de + rhoab*d2de ) - F83
            d3fgab(1) = -ab*( d2om*t + 2d0*d1om*d1t(1) + om*d2t(1) )
            d3fgab(2) = -ab*( d2om*t + d1om*( d1t(1) + d1t(2) ) 
     1                      + om*d2t(2) )
            d3fgab(3) = -ab*( d2om*t + 2d0*d1om*d1t(2) + om*d2t(3) )
#endif
c
#ifdef THIRD_DERIV
            d3t(1) = -F79*( 3.0d0*rhob*d2de + rhoab*d3de )
            d3t(2) = -F79*( 2.0d0*d1de + 2.0d0*rhob*d2de
     1                    + rhoa*d2de + rhoab*d3de )
            d3t(3) = -F79*( 2.0d0*d1de + rhob*d2de
     1                    + 2.0d0*rhoa*d2de + rhoab*d3de ) 
            d3t(4) = -F79*( 3.0d0*rhoa*d2de + rhoab*d3de ) 
c
            d4fgab(1) = -ab*( d3om*t + 3.0d0*d2om*d1t(1)
     1                      + 3.0d0*d1om*d2t(1) + om*d3t(1) )
            d4fgab(2) = -ab*( d3om*t + 3.0d0*d2om*d1t(1)
     1                      + d1om*(d2t(1) + 2.0d0*d2t(2))
     2                      + om*d3t(2) )
            d4fgab(3) = -ab*( d3om*t + 3.0d0*d2om*d1t(1)
     1                      + d1om*(d2t(1) + 2.0d0*d2t(2))
     2                      + om*d3t(3) )
            d4fgab(4) = -ab*( d3om*t + 3.0d0*d2om*d1t(1)
     1                      + 3.0d0*d1om*d2t(1) + om*d3t(4) )            
#endif
c Daniel (3-21-13): tm_in is what I call Si (the inside term)
            tm_in = 1.0d0 - 3.0d0*de - de11*rhob*rrho
c
c            t = F19*rhoab*( 1d0 - 3d0*de - de11*rhob*rrho ) - rhoa2
            t = F19*rhoab*tm_in - rhoa2
c Daniel (10-30-12): Derivatives of the inside term
            d1tm_in(1) = -(3.0d0+rhob*rrho)*d1de + de11*rhob*rrho2
            d1tm_in(2) = -(3.0d0+rhob*rrho)*d1de - de11*rhoa*rrho2
c
c            d1t(1) = F19*( rhob*( 1d0 - 3d0*de - de11*rhob*rrho )
c     1             + rhoab*( -( 3d0 + rhob*rrho )*d1de 
c     2                     + de11*rhob*rrho2 ) )
c     3             - 2d0*rhoa
c            d1t(2) = F19*( rhoa*( 1d0 - 3d0*de - de11*rhob*rrho )
c     1             - rhoab*( ( 3d0 + rhob*rrho )*d1de 
c     2                     + de11*rhoa*rrho2 ) )
            d1t(1) = F19*( rhob*tm_in + rhoab*d1tm_in(1) )
     1             - 2.0d0*rhoa
            d1t(2) = F19*( rhoa*tm_in + rhoab*d1tm_in(2) )
            d1f(5) = -ab*om*t
            d2fgbb(1) = -ab*( d1om*t + om*d1t(1) )
            d2fgbb(2) = -ab*( d1om*t + om*d1t(2) )
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c Daniel (10-30-12): Derivatives of the inside term, Ri.
            d2tm_in(1) = -(3.0d0+rhob*rrho)*d2de
     1                 + 2.0d0*d1de*rhob*rrho2
     2                 - 2.0d0*de11*rhob*rrho3
            d2tm_in(2) = -(3.0d0+rhoa*rrho)*d2de
     1                 + (rhob - rhoa)*d1de*rrho2
     2                 + (rhoa - rhob)*de11*rrho3
            d2tm_in(3) = -(3.0d0+rhoa*rrho)*d2de
     1                 - 2.0d0*d1de*rhoa*rrho2
     2                 + 2.0d0*de11*rhoa*rrho3
c
c            d2t(1) = -F19*(
c     &           2d0*rhob*((3d0+rhob*rrho)*d1de-de11*rhob*rrho2)
c     &           +  rhoab*((3d0+rhob*rrho)*d2de-2d0*d1de*rhob*rrho2
c     &                                         +2d0*de11*rhob*rrho3))
c     &           - 2d0
c            d2t(2) = F19*(1d0-3d0*de-de11*rhob*rrho
c     &           - rho(n,1)*(3d0+rhob*rrho)*d1de
c     &           - rhoab*((3d0+rhob*rrho)*d2de+d1de*(rhoa-rhob)*rrho2
c     &                                        -de11*(rhoa-rhob)*rrho3))
c            d2t(3) = -F19*(
c     &           2d0*rhoa*((3d0+rhob*rrho)*d1de+de11*rhoa*rrho2)
c     &           +  rhoab*((3d0+rhob*rrho)*d2de+2d0*d1de*rhoa*rrho2
c     &                                         -2d0*de11*rhoa*rrho3))
            d2t(1) = F19*( 2.0d0*rhob*d1tm_in(1) + rhoab*d2tm_in(1) )
     1             - 2.0d0
            d2t(2) = F19*( tm_in
     1                   - rho(n,1)*( 3.0d0 + rhob*rrho )*d1de
     2                   + rhoab*d2tm_in(2) )
            d2t(3) = F19*( 2.0d0*rhoa*d1tm_in(2) + rhoab*d2tm_in(3) )
c
            d3fgbb(1) = -ab*( d2om*t + 2.0d0*d1om*d1t(1) + om*d2t(1))
            d3fgbb(2) = -ab*( d2om*t + d1om*( d1t(1) + d1t(2) ) 
     1                      +om*d2t(2) )
            d3fgbb(3) = -ab*( d2om*t + 2.0d0*d1om*d1t(2) + om*d2t(3) )
#endif
c 
#ifdef THIRD_DERIV
            d3tm_in(1) = -( 3.0d0 + rhoa*rrho )*d3de
     1                 + 3.0d0*d2de*rhob*rrho2
     2                 - 6.0d0*d1de*rhob*rrho3
     3                 + 6.0d0*de11*rhob*rrho4
            d3tm_in(2) = -( 3.0d0 + rhoa*rrho )*d3de
     1                 - d2de*rhoa*rrho2
     2                 + 2.0d0*d2de*rhob*rrho2
     3                 + 2.0d0*d1de*rhoa*rrho3
     4                 - 4.0d0*d1de*rhob*rrho3
     5                 - 2.0d0*de11*(rhoa-2.0d0*rhob)*rrho4
            d3tm_in(3) = -( 3.0d0 + rhoa*rrho )*d3de
     1                 - 2.0d0*d2de*rhoa*rrho2
     2                 + d2de*rhob*rrho2
     3                 + 4.0d0*d1de*rhoa*rrho3
     4                 - 2.0d0*d1de*rhob*rrho3
     5                 - 2.0d0*de11*(2.0d0*rhoa-rhoa)*rrho4
            d3tm_in(4) = -( 3.0d0 + rhoa*rrho )*d3de
     1                 - 3.0d0*d2de*rhoa*rrho2
     2                 + 6.0d0*d1de*rhoa*rrho3
     3                 - 6.0d0*de11*rhoa*rrho4
c
            d3t(1) = F19*( 3.0d0*rhob*d2tm_in(1)
     1                   + rhoab*d3tm_in(1) )
            d3t(2) = F19*( 2.0d0*d1tm_in(1)
     1                   + rhoa*d2tm_in(1)
     2                   + 2.0d0*rhob*d2tm_in(2)
     3                   + rhoab*d3tm_in(2) )
            d3t(3) = F19*( 2.0d0*d1tm_in(2)
     1                   + 2.0d0*rhoa*d2tm_in(2)
     2                   + rhob*d2tm_in(3)
     3                   + rhoab*d3tm_in(3) )
            d3t(4) = F19*( 3.0d0*rhoa*d2tm_in(3)
     1                   + rhoab*d3tm_in(4) )
c
            d4fgbb(1) = -ab*( d3om*t + 3.0d0*d2om*d1t(1)
     1                      + 3.0d0*d1om*d2t(1) + om*d3t(1) )
            d4fgbb(2) = -ab*( d3om*t + d2om*(2.0d0*d1t(1) + d1t(2))
     1                      + d1om*(d2t(1) + 2.0d0*d2t(2)) + om*d3t(2) )
            d4fgbb(3) = -ab*( d3om*t + d2om*(d1t(1) + 2.0d0*d1t(2))
     1                      + d1om*(2.0d0*d2t(2) + d2t(3)) + om*d3t(3) )
            d4fgbb(4) = -ab*( d3om*t + 3.0d0*d2om*d1t(2)
     1                      + 3.0d0*d1om*d2t(3) + om*d3t(4) )
#endif
c
            f = c1*f1 + c2*f2 + gaa*d1f(3) + gab*d1f(4) + gbb*d1f(5)
            d1f(1) = c1*d1f1(1) + c2*d1f2(1)
     &             + gaa*d2fgaa(1) + gab*d2fgab(1) + gbb*d2fgbb(1)
            d1f(2) = c1*d1f1(2) + c2*d1f2(2)
     &             + gaa*d2fgaa(2) + gab*d2fgab(2) + gbb*d2fgbb(2)
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
            d2f(1) = c1*d2f1(1) + c2*d2f2(1)
     &             + gaa*d3fgaa(1) + gab*d3fgab(1) + gbb*d3fgbb(1)
            d2f(2) = c1*d2f1(2) + c2*d2f2(2)
     &             + gaa*d3fgaa(2) + gab*d3fgab(2) + gbb*d3fgbb(2)
            d2f(3) = c1*d2f1(3) + c2*d2f2(3)
     &             + gaa*d3fgaa(3) + gab*d3fgab(3) + gbb*d3fgbb(3)
#endif
c
#ifdef THIRD_DERIV
            d3f(1) = c1*d3f1(1) + c2*d3f2(1)
     1             + gaa*d4fgaa(1) + gab*d4fgab(1) + gbb*d4fgbb(1)
            d3f(2) = c1*d3f1(2) + c2*d3f2(2)
     1             + gaa*d4fgaa(2) + gab*d4fgab(2) + gbb*d4fgbb(2)
            d3f(3) = c1*d3f1(3) + c2*d3f2(3)
     1             + gaa*d4fgaa(3) + gab*d4fgab(3) + gbb*d4fgbb(3)
            d3f(4) = c1*d3f1(4) + c2*d3f2(4)
     1             + gaa*d4fgaa(4) + gab*d4fgab(4) + gbb*d4fgbb(4)
#endif
c
            Ec = Ec + f*fac*qwght(n)
            if (ldew) func(n) = func(n) + f*fac
            Amat(n,1) = Amat(n,1) + d1f(1)*fac
            Amat(n,2) = Amat(n,2) + d1f(2)*fac
            Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + d1f(3)*fac
            Cmat(n,D1_GAB) = Cmat(n,D1_GAB) + d1f(4)*fac
            Cmat(n,D1_GBB) = Cmat(n,D1_GBB) + d1f(5)*fac
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
            Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + d2f(1)*fac
            Amat2(n,D2_RA_RB) = Amat2(n,D2_RA_RB) + d2f(2)*fac
            Amat2(n,D2_RB_RB) = Amat2(n,D2_RB_RB) + d2f(3)*fac
            Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA) + d2fgaa(1)*fac
            Cmat2(n,D2_RA_GAB) = Cmat2(n,D2_RA_GAB) + d2fgab(1)*fac
            Cmat2(n,D2_RA_GBB) = Cmat2(n,D2_RA_GBB) + d2fgbb(1)*fac
            Cmat2(n,D2_RB_GAA) = Cmat2(n,D2_RB_GAA) + d2fgaa(2)*fac
            Cmat2(n,D2_RB_GAB) = Cmat2(n,D2_RB_GAB) + d2fgab(2)*fac
            Cmat2(n,D2_RB_GBB) = Cmat2(n,D2_RB_GBB) + d2fgbb(2)*fac
#endif
c
#ifdef THIRD_DERIV
            Amat3(n,D3_RA_RA_RA) = Amat3(n,D3_RA_RA_RA) + d3f(1)*fac
            Amat3(n,D3_RA_RA_RB) = Amat3(n,D3_RA_RA_RB) + d3f(2)*fac
            Amat3(n,D3_RA_RB_RB) = Amat3(n,D3_RA_RB_RB) + d3f(3)*fac
            Amat3(n,D3_RB_RB_RB) = Amat3(n,D3_RB_RB_RB) + d3f(4)*fac
c
            Cmat3(n,D3_RA_RA_GAA) = Cmat3(n,D3_RA_RA_GAA)
     1                            + d3fgaa(1)*fac
            Cmat3(n,D3_RA_RA_GAB) = Cmat3(n,D3_RA_RA_GAB)
     1                            + d3fgab(1)*fac
            Cmat3(n,D3_RA_RA_GBB) = Cmat3(n,D3_RA_RA_GBB)
     1                            + d3fgbb(1)*fac
c
            Cmat3(n,D3_RA_RB_GAA) = Cmat3(n,D3_RA_RB_GAA)
     1                            + d3fgaa(2)*fac
            Cmat3(n,D3_RA_RB_GAB) = Cmat3(n,D3_RA_RB_GAB)
     1                            + d3fgab(2)*fac
            Cmat3(n,D3_RA_RB_GBB) = Cmat3(n,D3_RA_RB_GBB)
     1                            + d3fgbb(2)*fac
c
            Cmat3(n,D3_RB_RB_GAA) = Cmat3(n,D3_RB_RB_GAA)
     1                            + d3fgaa(3)*fac
            Cmat3(n,D3_RB_RB_GAB) = Cmat3(n,D3_RB_RB_GAB)
     1                            + d3fgab(3)*fac
            Cmat3(n,D3_RB_RB_GBB) = Cmat3(n,D3_RB_RB_GBB)
     1                            + d3fgbb(3)*fac
#endif
c
   20    continue
      endif
      return
      end
#ifndef SECOND_DERIV
#define SECOND_DERIV
c
c     Compile source again for the 2nd derivative case
c
#include "xc_lyp88.F"
#endif
c
#ifndef THIRD_DERIV
#define THIRD_DERIV
c
c     Compile source again for the 3rd derivative case
c
#include "xc_lyp88.F"
#endif
