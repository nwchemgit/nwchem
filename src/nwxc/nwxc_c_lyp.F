#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if !defined(NWAD_PRINT)
C> \ingroup nwxc
C> @{
C>
C> \file nwxc_c_lyp.F
C> The LYP correlation functional
C>
C> @}
#endif
#endif
C>
C> \ingroup nwxc_priv
C> @{
C>
C> \brief Evaluate the LYP correlation functional
C>
C> Evaluate the  LYP correlation functional [1-4].
C>
C> ### References ###
C>
C> [1] C. Lee, W. Yang, R.G. Parr, 
C> "Development of the Colle-Salvetti correlation-energy formula into
C> a functional of the electron density", Phys. Rev. B <b>37</b>,
C> 785-789 (1988), DOI: 
C> <a href="https://doi.org/10.1103/PhysRevB.37.785">
C> 10.1103/PhysRevB.37.785</a>.
C>
C> [2]  R. Colle, O. Salvetti,
C> "Approximate calculation of the correlation energy for the closed
C> shells", Theor. Chim. Acta <b>37</b>, 329-334 (1975), DOI:
C> <a href="https://doi.org/10.1007/BF01028401">
C> 10.1007/BF01028401</a>.
C>
C> [3] B. Miehlich, A. Savin, H. Stoll, H. Preuss,
C> "Results obtained with the correlation energy density functionals of
C> Becke, and Lee, Yang and Parr", Chem. Phys. Lett. <b>157</b>, 200-206
C> (1989), DOI: <a href="https://doi.org/10.1016/0009-2614(89)87234-3">
C> 10.1016/0009-2614(89)87234-3</a>.
C>
C> [4] B.G. Johnson, P.M.W. Gill, J.A. Pople,
C> "The performance of a family of density functional methods",
C> J. Chem. Phys. <b>98</b>, 5612-5626 (1993), DOI:
C> <a href="https://doi.org/10.1063/1.464906">
C> 10.1063/1.464906</a>.
C>
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
      Subroutine nwxc_c_lyp_p(tol_rho, ipol, nq, wght, rho, rgamma, 
     &                      func)
#else
      Subroutine nwxc_c_lyp(tol_rho, ipol, nq, wght, rho, rgamma, 
     &                      func)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      Subroutine nwxc_c_lyp_d2(tol_rho, ipol, nq, wght, rho,
     &                         rgamma, func)
#else
      Subroutine nwxc_c_lyp_d3(tol_rho, ipol, nq, wght, rho,
     &                         rgamma, func)
#endif
c
C$Id$
c
#include "nwad.fh"
      implicit none 
c
#include "nwxc_param.fh"
c
c     Input and other parameters
c
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
      type(nwad_dble)::func(nq)    !< [Output] The value of the functional
c     double precision Amat(nq,*)   !< [Output] The derivative wrt rho
c     double precision Cmat(nq,*)   !< [Output] The derivative wrt rgamma
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c
c     Sampling Matrices for the XC Kernel
c
c     double precision Amat2(nq,*)  !< [Output] The 2nd derivative wrt rho
c     double precision Cmat2(nq,*)  !< [Output] The 2nd derivative wrt rgamma
                                    !< and possibly rho
#endif
#if defined(THIRD_DERIV)
c
c     Sampling Matrices for the XC Kernel
c
c     double precision Amat3(nq,*)  !< [Output] The 3rd derivative wrt rho
c     double precision Cmat3(nq,*)  !< [Output] The 3rd derivative wrt rgamma
                                    !< and possibly rho
#endif
      double precision F13, F43, F113, F83, F53, F19, F79, P1, 
     &                 A, B, C, D
c
      Parameter (F13 = 1.D0/3.D0, F43 = 4.D0*F13, F113 = 11.D0*F13,
     &           F83 = 8.D0*F13, F53 = 5.D0*F13, F19 = 1.D0/9.D0,
     &           F79 = 7.D0*F19)
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
      type(nwad_dble)::rrho, rhoa, rhob, d1t(2)
      type(nwad_dble)::rhoa83, rhob83
      type(nwad_dble)::rrho2, rhoa2, rhob2, rhoab, rho2
      type(nwad_dble)::h1, h2, h3, om, t, tm_in
      type(nwad_dble)::de, de11, de47
      type(nwad_dble)::zero
      type(nwad_dble)::xrarb
c
      double precision d1xrarb(2)
      double precision d1tm_in(2)
c
      type(nwad_dble)::gaa, gab, gbb
      type(nwad_dble)::f1, f2, d1f1(2), d1f2(2), f, d1f(5),
     &     d2fgaa(2), d2fgab(2), d2fgbb(2)
c
c     Coefficients of first two terms in LYP functional and other
c     commonly occurring factors
c
      zero = 0.0d0
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
            if (rho(n,R_T).lt.tol_rho)goto 10
            rrho = 1.0d0/rho(n,R_T)
            rhoa = 0.5d0*rho(n,R_T)
            rrho2 = rrho*rrho
            rho2 = 1.0d0/rrho2
            rhoa2 = rhoa*rhoa
            rhoab = rhoa2
c           rhoa53 = rhoa**F53
            rhoa83 = rhoa**F83
c
            h2 = d*rrho**F13
c           d1h2 = -F13*h2*rrho
c
            h3 = ratcd*h2
c           d1h3 = ratcd*d1h2
c
            h1 = 1d0/(1d0+h2)
c           d1h1 = -h1*h1*d1h2
c
            om = exp(-h3)*h1*rrho**F113
c           t = d1h3+h1*d1h2+F113*rrho
c           d1om = -om*t
c
            de = h3+h1*h2
c           d1de = d1h3 + d1h1*h2 + h1*d1h2
c
            f1 = h1*rhoab*rrho
c           d1f1(1) = d1h1*rhoab*rrho - h1*rhoab*rrho2
c           d1f1(1) = d1f1(1) + h1*rhoa*rrho
c
            xrarb = rhoab*(rhoa83+rhoa83)
            f2 = om*xrarb
c           d1xrarb(1) = rhoa*(F113*rhoa83+rhoa83)
c           d1f2(1) = d1om*xrarb + om*d1xrarb(1)
c
c           gaa =(delrho(n,1,1)*delrho(n,1,1) +
c    &            delrho(n,2,1)*delrho(n,2,1) +
c    &            delrho(n,3,1)*delrho(n,3,1))*0.25d0
            gaa = rgamma(n,G_TT)*0.25d0
c
            de11 = de - 11.0d0
            de47 = 47.0d0 - 7.0d0*de
c
c Daniel (10-30-12): tm_in is what I call Qi (the inside term)
            tm_in = 1.0d0 - 3.0d0*de - de11*rhoa*rrho
c Daniel (10-23-12): "t" is what I call Q or S.
            t = F19*rhoab*tm_in - rhoa2
c Daniel (10-30-12): Derivatives of the inside term
c           d1tm_in(1) = -(3.0d0+rhoa*rrho)*d1de-de11*rhoa*rrho2
c           d1tm_in(2) = -(3.0d0+rhoa*rrho)*d1de+de11*rhoa*rrho2
c Daniel (10-23-12): d1t(1) is the derivative with respect to rhoa,
c and d1t(2) is the derivative with respect to rhob.
c           d1t(1) = F19*( rhoa*tm_in + rhoab*d1tm_in(1) )
c           d1t(2) = F19*( rhoa*tm_in + rhoab*d1tm_in(2) )
c    &             - 2.0d0*rhoa
c
            d1f(3) = -ab*om*t
c
c           d2fgaa(1) = -ab*( d1om*t + om*d1t(1) )
c           d2fgaa(2) = -ab*( d1om*t + om*d1t(2) )
c
c Daniel (10-23-12): "t" is what I call R.
            t = F19*rhoab*de47-F43*rho2
c           d1t(1) = F19*rhoa*de47 - F79*rhoab*d1de - F83*rho(n,R_T)
            d1f(4) = -ab*om*t
c           d2fgab(1) = -ab*(d1om*t+om*d1t(1))
c
c           d2fgbb(1) = d2fgaa(2)
c
            f = c1*f1 + c2*f2 + gaa*(2d0*d1f(3) + d1f(4))
c           d1f(1) = c1*d1f1(1) + c2*d1f2(1)
c    &             + gaa*(d2fgaa(1) + d2fgab(1) + d2fgbb(1))
c
            func(n) = func(n) + f*wght
c           Amat(n,D1_RA) = Amat(n,D1_RA) + d1f(1)*wght
c           Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + d1f(3)*wght
c           Cmat(n,D1_GAB) = Cmat(n,D1_GAB) + d1f(4)*wght
c
   10    continue
      else
c
c        ======> SPIN-UNRESTRICTED <======
c
         do 20 n = 1,nq
            if (rho(n,R_A)+rho(n,R_B).lt.tol_rho)goto 20
            rrho = 1.0d0/(rho(n,R_A)+rho(n,R_B))
c           rhoa = max(zero,rho(n,R_A))
c           rhob = max(zero,rho(n,R_B))
            if (rho(n,R_A).ge.0.5d0*tol_rho) then
              rhoa = rho(n,R_A)
            else
              rhoa = 0.0d0
              rrho = 1.0d0/rho(n,R_B)
            endif
            if (rho(n,R_B).ge.0.5d0*tol_rho) then
              rhob = rho(n,R_B)
            else
              rhob = 0.0d0
              rrho = 1.0d0/rho(n,R_A)
            endif
            rrho2 = rrho*rrho
            rho2 = 1d0/rrho2
            rhoa2 = rhoa*rhoa
            rhob2 = rhob*rhob
            rhoab = rhoa*rhob
c           rhoa53 = rhoa**F53
c           rhob53 = rhob**F53
            rhoa83 = rhoa**F83
            rhob83 = rhob**F83
c
            h2 = d*rrho**F13
c           d1h2 = -F13*h2*rrho
c
            h3 = ratcd*h2
c           d1h3 = ratcd*d1h2
c
            h1 = 1d0/(1d0+h2)
c           d1h1 = -h1*h1*d1h2
c
            om = exp(-h3)*h1*rrho**F113
c           t = d1h3+h1*d1h2+F113*rrho
c           d1om = -om*t
c
            de = h3+h1*h2
c           d1de = d1h3 + d1h1*h2 + h1*d1h2
c
c Daniel (3-21-13): f1 is J/(-4*a) in my notes.
            f1 = h1*rhoab*rrho
c           d1f1(1) = d1h1*rhoab*rrho - h1*rhoab*rrho2
c           d1f1(2) = d1f1(1)
c           d1f1(1) = d1f1(1) + h1*rhob*rrho
c           d1f1(2) = d1f1(2) + h1*rhoa*rrho
c
c Daniel (10-30-12): Define xrarb here
            xrarb = rhoab*(rhoa83+rhob83)
c
            f2 =om*xrarb
c
c           d1xrarb(1) = rhob*(F113*rhoa83+rhob83)
c           d1xrarb(2) = rhoa*(F113*rhob83+rhoa83)
c
c           d1f2(1) = d1om*xrarb
c           d1f2(2) = d1f2(1)
c           d1f2(1) = d1f2(1) + om*d1xrarb(1)
c           d1f2(2) = d1f2(2) + om*d1xrarb(2)
c
            gaa = rgamma(n,G_AA)
            gbb = rgamma(n,G_BB)
            gab = rgamma(n,G_AB)
            if (rho(n,R_A).lt.0.5d0*tol_rho) then
              gaa = 0.0d0
              gab = 0.0d0
            endif
            if (rho(n,R_B).lt.0.5d0*tol_rho) then
              gbb = 0.0d0
              gab = 0.0d0
            endif
c
            de11 = de - 11d0
            de47 = 47d0 - 7d0*de
c Daniel (10-30-12): tm_in is what I call Qi (the inside term)
            tm_in = 1.0d0 - 3.0d0*de - de11*rhoa*rrho
            t = F19*rhoab*tm_in-rhob2
c Daniel (10-30-12): Derivatives of the inside term
c           d1tm_in(1) = -(3.0d0+rhoa*rrho)*d1de-de11*rhob*rrho2
c           d1tm_in(2) = -(3.0d0+rhoa*rrho)*d1de+de11*rhoa*rrho2
c           d1t(1) = F19*(rhob*tm_in + rhoab*d1tm_in(1))
c           d1t(2) = F19*(rhoa*tm_in + rhoab*d1tm_in(2)) - 2.0d0*rhob
c
            d1f(3) = -ab*om*t
c           d2fgaa(1) = -ab*(d1om*t+om*d1t(1))
c           d2fgaa(2) = -ab*(d1om*t+om*d1t(2))
c
            t = F19*rhoab*de47-F43*rho2
c           d1t(1) = F19*rhob*de47 - F79*rhoab*d1de 
c    &             - F83*(rho(n,R_A)+rho(n,R_B))
c           d1t(2) = F19*rhoa*de47 - F79*rhoab*d1de
c    &             - F83*(rho(n,R_A)+rho(n,R_B))
            d1f(4) = -ab*om*t
c           d2fgab(1) = -ab*(d1om*t+om*d1t(1))
c           d2fgab(2) = -ab*(d1om*t+om*d1t(2))
c Daniel (3-21-13): tm_in is what I call Si (the inside term)
            tm_in = 1.0d0 - 3.0d0*de - de11*rhob*rrho
            t = F19*rhoab*tm_in - rhoa2
c Daniel (10-30-12): Derivatives of the inside term
c           d1tm_in(1) = -(3.0d0+rhob*rrho)*d1de + de11*rhob*rrho2
c           d1tm_in(2) = -(3.0d0+rhob*rrho)*d1de - de11*rhoa*rrho2
c           d1t(1) = F19*( rhob*tm_in + rhoab*d1tm_in(1) )
c    1             - 2.0d0*rhoa
c           d1t(2) = F19*( rhoa*tm_in + rhoab*d1tm_in(2) )
            d1f(5) = -ab*om*t
c           d2fgbb(1) = -ab*( d1om*t + om*d1t(1) )
c           d2fgbb(2) = -ab*( d1om*t + om*d1t(2) )
c
            f = c1*f1 + c2*f2 + gaa*d1f(3) + gab*d1f(4) + gbb*d1f(5)
c           d1f(1) = c1*d1f1(1) + c2*d1f2(1)
c    &             + gaa*d2fgaa(1) + gab*d2fgab(1) + gbb*d2fgbb(1)
c           d1f(2) = c1*d1f1(2) + c2*d1f2(2)
c    &             + gaa*d2fgaa(2) + gab*d2fgab(2) + gbb*d2fgbb(2)
c
            func(n) = func(n) + f*wght
c           Amat(n,D1_RA) = Amat(n,D1_RA) + d1f(1)*wght
c           Amat(n,D1_RB) = Amat(n,D1_RB) + d1f(2)*wght
c           Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + d1f(3)*wght
c           Cmat(n,D1_GAB) = Cmat(n,D1_GAB) + d1f(4)*wght
c           Cmat(n,D1_GBB) = Cmat(n,D1_GBB) + d1f(5)*wght
c
   20    continue
      endif
      return
      end
#ifndef NWAD_PRINT
#define NWAD_PRINT
c
c     Compile source again for Maxima
c
#include "nwxc_c_lyp.F"
#endif
#ifndef SECOND_DERIV
#define SECOND_DERIV
c
c     Compile source again for the 2nd derivative case
c
#include "nwxc_c_lyp.F"
#endif
#ifndef THIRD_DERIV
#define THIRD_DERIV
c
c     Compile source again for the 3rd derivative case
c
#include "nwxc_c_lyp.F"
#endif
#undef NWAD_PRINT
C> @}
