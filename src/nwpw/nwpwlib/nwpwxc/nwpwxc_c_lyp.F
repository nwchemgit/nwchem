#ifndef SECOND_DERIV
C> \ingroup nwpwxc
C> @{
C>
C> \file nwpwxc_c_lyp.F
C> The LYP correlation functional
C>
C> @}
#endif
C>
C> \ingroup nwpwxc_priv
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
#ifndef SECOND_DERIV
      Subroutine nwpwxc_c_lyp(tol_rho, ipol, nq, wght, rho, rgamma, 
     &                      func, Amat, Cmat) 
#else
      Subroutine nwpwxc_c_lyp_d2(tol_rho, ipol, nq, wght, rho,
     &                         rgamma, func, Amat, Amat2, Cmat, Cmat2)
#endif
c
C$Id$
c
      implicit none 
c
#include "nwpwxc_param.fh"
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
      double precision rho(nq,*)    !< [Input] The density
c
c     Charge Density Gradient
c
      double precision rgamma(nq,*) !< [Input] The norm of the density gradients
c
c     Sampling Matrices for the XC Potential
c
      double precision func(nq)    !< [Output] The value of the functional
      double precision Amat(nq,*)   !< [Output] The derivative wrt rho
      double precision Cmat(nq,*)   !< [Output] The derivative wrt rgamma
#ifdef SECOND_DERIV
c
c     Sampling Matrices for the XC Kernel
c
      double precision Amat2(nq,*)  !< [Output] The 2nd derivative wrt rho
      double precision Cmat2(nq,*)  !< [Output] The 2nd derivative wrt rgamma
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
      double precision rrho, rhoa, rhob, rrho2, rhoa2, rhob2, rhoab,
     &     rhoa53, rhob53, rhoa83, rhob83, rho2,
     &     h1, h2, h3, d1h1, d1h2, d1h3, om, d1om, de, d1de, de11, de47,
     &     t, d1t(2)
      double precision gaa, gab, gbb
      double precision f1, f2, d1f1(2), d1f2(2), f, d1f(5),
     &     d2fgaa(2), d2fgab(2), d2fgbb(2)
#ifdef SECOND_DERIV
      double precision d2h1, d2h2, d2h3, d2om, d2de, rrho3, d2f1(3),
     &     d2f2(3), rhoa113, rhob113, d3fgaa(3), d3fgab(3), d3fgbb(3),
     &     d2t(3), d2f(3)
#endif
c
c     Coefficients of first two terms in LYP functional and other
c     commonly occurring factors
c
      c1 = -4d0*a
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
            rrho = 1d0/rho(n,R_T)
            rhoa = 0.5d0*rho(n,R_T)
            rrho2 = rrho*rrho
            rho2 = 1d0/rrho2
            rhoa2 = rhoa*rhoa
            rhoab = rhoa2
            rhoa53 = rhoa**F53
            rhoa83 = rhoa**F83
#ifdef SECOND_DERIV
            rrho3 = rrho*rrho2
            rhoa113 = rhoa*rhoa83
#endif
c
            h2 = d*rrho**F13
            d1h2 = -F13*h2*rrho
#ifdef SECOND_DERIV
            d2h2 = -F43*d1h2*rrho
#endif
c
            h3 = ratcd*h2
            d1h3 = ratcd*d1h2
#ifdef SECOND_DERIV
            d2h3 = ratcd*d2h2
#endif
c
            h1 = 1d0/(1d0+h2)
            d1h1 = -h1*h1*d1h2
#ifdef SECOND_DERIV
            d2h1 = -(2d0*h1*d1h1*d1h2 + h1*h1*d2h2)
#endif
c
            om = exp(-h3)*h1*rrho**F113
            t = d1h3+h1*d1h2+F113*rrho
            d1om = -om*t
#ifdef SECOND_DERIV
            d2om = -(d1om*t+om*(d2h3+d1h1*d1h2+h1*d2h2-F113*rrho2))
#endif
c
            de = h3+h1*h2
            d1de = d1h3 + d1h1*h2 + h1*d1h2
#ifdef SECOND_DERIV
            d2de = d2h3 + d2h1*h2 + h1*d2h2 + 2d0*d1h1*d1h2
#endif
c
            f1 = h1*rhoab*rrho
            d1f1(1) = d1h1*rhoab*rrho - h1*rhoab*rrho2
            d1f1(1) = d1f1(1) + h1*rhoa*rrho
#ifdef SECOND_DERIV
            d2f1(1) = d2h1*rhoab*rrho + 2d0*d1h1*(rhoa*rrho-rhoab*rrho2)
     &           + 2d0*h1*(-rhoa*rrho2+rhoab*rrho3)
            d2f1(2) = d2h1*rhoab*rrho + d1h1*(1d0-2d0*rhoab*rrho2)
     &           + 2d0*h1*rhoab*rrho3
#endif
c
            f2 = om*rhoab*(rhoa83+rhoa83)
            d1f2(1) = d1om*rhoab*(rhoa83+rhoa83)
            d1f2(1) = d1f2(1) + om*rhoa*(F113*rhoa83+rhoa83)
#ifdef SECOND_DERIV
            d2f2(1) = d2om*rhoab*(rhoa83+rhoa83)
     &          + 2d0*d1om*rhoa*(F113*rhoa83+rhoa83)
     &          +       om*rhoa*F113*F83*rhoa53
            d2f2(2) = d2om*rhoab*(rhoa83+rhoa83)
     &           + d1om*(rhoa113+rhoa113+F113*(rhoa*rhoa83+rhoa*rhoa83))
     &           +   om*F113*(rhoa83+rhoa83)
#endif
c
c           gaa =(delrho(n,1,1)*delrho(n,1,1) +
c    &            delrho(n,2,1)*delrho(n,2,1) +
c    &            delrho(n,3,1)*delrho(n,3,1))*0.25d0
            gaa = rgamma(n,G_TT)*0.25d0
c
            de11 = de - 11d0
            de47 = 47d0 - 7d0*de
c
            t = F19*rhoab*(1d0-3d0*de-de11*rhoa*rrho)-rhoa2
            d1t(1) = F19*(rhoa*(1d0-3d0*de-de11*rhoa*rrho)
     &             - rhoab*((3d0+rhoa*rrho)*d1de+de11*rhoa*rrho2))
            d1t(2) = F19*(rhoa*(1d0-3d0*de-de11*rhoa*rrho)
     &             + rhoab*(-(3d0+rhoa*rrho)*d1de+de11*rhoa*rrho2))
     &             - 2d0*rhoa
            d1f(3) = -ab*om*t
            d2fgaa(1) = -ab*(d1om*t+om*d1t(1))
            d2fgaa(2) = -ab*(d1om*t+om*d1t(2))
#ifdef SECOND_DERIV
            d2t(1) = -F19*(
     &           2d0*rhoa*((3d0+rhoa*rrho)*d1de+de11*rhoa*rrho2)
     &           +  rhoab*((3d0+rhoa*rrho)*d2de+2d0*d1de*rhoa*rrho2
     &                                         -2d0*de11*rhoa*rrho3))
            d2t(2) = F19*(1d0-3d0*de-de11*rhoa*rrho
     &           - rho(n,R_T)*(3d0+rhoa*rrho)*d1de
     &           - rhoab*((3d0+rhoa*rrho)*d2de))
            d2t(3) = -F19*(
     &           2d0*rhoa*((3d0+rhoa*rrho)*d1de-de11*rhoa*rrho2)
     &           +  rhoab*((3d0+rhoa*rrho)*d2de-2d0*d1de*rhoa*rrho2
     &                                         +2d0*de11*rhoa*rrho3))
     &           - 2d0
            d3fgaa(1) = -ab*(d2om*t+2d0*d1om*d1t(1)+om*d2t(1))
            d3fgaa(2) = -ab*(d2om*t+d1om*(d1t(1)+d1t(2))+om*d2t(2))
            d3fgaa(3) = -ab*(d2om*t+2d0*d1om*d1t(2)+om*d2t(3))
#endif
c
            t = F19*rhoab*de47-F43*rho2
            d1t(1) = F19*rhoa*de47 - F79*rhoab*d1de - F83*rho(n,R_T)
            d1f(4) = -ab*om*t
            d2fgab(1) = -ab*(d1om*t+om*d1t(1))
#ifdef SECOND_DERIV
            d2t(1) = -F79*(2d0*rhoa*d1de+rhoab*d2de) - F83
            d2t(2) = F19*de47 - F79*(rho(n,R_T)*d1de+rhoab*d2de) - F83
            d3fgab(1) = -ab*(d2om*t+2d0*d1om*d1t(1)+om*d2t(1))
            d3fgab(2) = -ab*(d2om*t+2d0*d1om*d1t(1)+om*d2t(2))
#endif
c
            d2fgbb(1) = d2fgaa(2)
#ifdef SECOND_DERIV
            d3fgbb(1) = d3fgaa(3)
            d3fgbb(2) = d3fgaa(2)
#endif
c
            f = c1*f1 + c2*f2 + gaa*(2d0*d1f(3) + d1f(4))
            d1f(1) = c1*d1f1(1) + c2*d1f2(1)
     &             + gaa*(d2fgaa(1) + d2fgab(1) + d2fgbb(1))
#ifdef SECOND_DERIV
            d2f(1) = c1*d2f1(1) + c2*d2f2(1)
     &             + gaa*(d3fgaa(1) + d3fgab(1) + d3fgbb(1))
            d2f(2) = c1*d2f1(2) + c2*d2f2(2)
     &             + gaa*(d3fgaa(2) + d3fgab(2) + d3fgbb(2))
#endif
c
            func(n) = func(n) + f*wght
            Amat(n,D1_RA) = Amat(n,D1_RA) + d1f(1)*wght
            Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + d1f(3)*wght
            Cmat(n,D1_GAB) = Cmat(n,D1_GAB) + d1f(4)*wght
#ifdef SECOND_DERIV
            Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + d2f(1)*wght
            Amat2(n,D2_RA_RB) = Amat2(n,D2_RA_RB) + d2f(2)*wght
            Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA) + d2fgaa(1)*wght
            Cmat2(n,D2_RA_GAB) = Cmat2(n,D2_RA_GAB) + d2fgab(1)*wght
            Cmat2(n,D2_RA_GBB) = Cmat2(n,D2_RA_GBB) + d2fgbb(1)*wght
#endif
c
   10    continue
      else
c
c        ======> SPIN-UNRESTRICTED <======
c
         do 20 n = 1,nq
            if (rho(n,R_A)+rho(n,R_B).lt.tol_rho)goto 20
            rrho = 1d0/(rho(n,R_A)+rho(n,R_B))
            rhoa = max(0.0d0,rho(n,R_A))
            rhob = max(0.0d0,rho(n,R_B))
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
#ifdef SECOND_DERIV
            rrho3 = rrho*rrho2
            rhoa113 = rhoa*rhoa83
            rhob113 = rhob*rhob83
#endif
c
cedo            h2 = d*abs(rrho)**F13*sign(1d0,rrho)
            h2 = d*rrho**F13
            d1h2 = -F13*h2*rrho
#ifdef SECOND_DERIV
            d2h2 = -F43*d1h2*rrho
#endif
c
            h3 = ratcd*h2
            d1h3 = ratcd*d1h2
#ifdef SECOND_DERIV
            d2h3 = ratcd*d2h2
#endif
c
            h1 = 1d0/(1d0+h2)
            d1h1 = -h1*h1*d1h2
#ifdef SECOND_DERIV
            d2h1 = -(2d0*h1*d1h1*d1h2 + h1*h1*d2h2)
#endif
c
            om = exp(-h3)*h1*rrho**F113
cedo            om = exp(-h3)*h1*abs(rrho)**F113*sign(1d0,rrho)
            t = d1h3+h1*d1h2+F113*rrho
            d1om = -om*t
#ifdef SECOND_DERIV
            d2om = -(d1om*t+om*(d2h3+d1h1*d1h2+h1*d2h2-F113*rrho2))
#endif
c
            de = h3+h1*h2
            d1de = d1h3 + d1h1*h2 + h1*d1h2
#ifdef SECOND_DERIV
            d2de = d2h3 + d2h1*h2 + h1*d2h2 + 2d0*d1h1*d1h2
#endif
c
            f1 = h1*rhoab*rrho
            d1f1(1) = d1h1*rhoab*rrho - h1*rhoab*rrho2
            d1f1(2) = d1f1(1)
            d1f1(1) = d1f1(1) + h1*rhob*rrho
            d1f1(2) = d1f1(2) + h1*rhoa*rrho
#ifdef SECOND_DERIV
            d2f1(1) = d2h1*rhoab*rrho + 2d0*d1h1*(rhob*rrho-rhoab*rrho2)
     &           + 2d0*h1*(-rhob*rrho2+rhoab*rrho3)
            d2f1(2) = d2h1*rhoab*rrho + d1h1*(1d0-2d0*rhoab*rrho2)
     &           + 2d0*h1*rhoab*rrho3
            d2f1(3) = d2h1*rhoab*rrho + 2d0*d1h1*(rhoa*rrho-rhoab*rrho2)
     &           + 2d0*h1*(-rhoa*rrho2+rhoab*rrho3)
#endif
c
            f2 = om*rhoab*(rhoa83+rhob83)
            d1f2(1) = d1om*rhoab*(rhoa83+rhob83)
            d1f2(2) = d1f2(1)
            d1f2(1) = d1f2(1) + om*rhob*(F113*rhoa83+rhob83)
            d1f2(2) = d1f2(2) + om*rhoa*(F113*rhob83+rhoa83)
#ifdef SECOND_DERIV
            d2f2(1) = d2om*rhoab*(rhoa83+rhob83)
     &          + 2d0*d1om*rhob*(F113*rhoa83+rhob83)
     &          +       om*rhob*F113*F83*rhoa53
            d2f2(2) = d2om*rhoab*(rhoa83+rhob83)
     &           + d1om*(rhoa113+rhob113+F113*(rhob*rhoa83+rhoa*rhob83))
     &           +   om*F113*(rhoa83+rhob83)
            d2f2(3) = d2om*rhoab*(rhoa83+rhob83)
     &          + 2d0*d1om*rhoa*(F113*rhob83+rhoa83)
     &          +       om*rhoa*F113*F83*rhob53
#endif
c
c           gaa = delrho(n,1,1)*delrho(n,1,1) +
c    &            delrho(n,2,1)*delrho(n,2,1) +
c    &            delrho(n,3,1)*delrho(n,3,1)
c           gab = delrho(n,1,1)*delrho(n,1,2) +
c    &            delrho(n,2,1)*delrho(n,2,2) +
c    &            delrho(n,3,1)*delrho(n,3,2)
c           gbb = delrho(n,1,2)*delrho(n,1,2) +
c    &            delrho(n,2,2)*delrho(n,2,2) +
c    &            delrho(n,3,2)*delrho(n,3,2)
c
            gaa = rgamma(n,G_AA)
            gbb = rgamma(n,G_BB)
            gab = rgamma(n,G_AB)
c
            de11 = de - 11d0
            de47 = 47d0 - 7d0*de
c
            t = F19*rhoab*(1d0-3d0*de-de11*rhoa*rrho)-rhob2
            d1t(1) = F19*(rhob*(1d0-3d0*de-de11*rhoa*rrho)
     &             - rhoab*((3d0+rhoa*rrho)*d1de+de11*rhob*rrho2))
            d1t(2) = F19*(rhoa*(1d0-3d0*de-de11*rhoa*rrho)
     &             + rhoab*(-(3d0+rhoa*rrho)*d1de+de11*rhoa*rrho2))
     &             - 2d0*rhob
            d1f(3) = -ab*om*t
            d2fgaa(1) = -ab*(d1om*t+om*d1t(1))
            d2fgaa(2) = -ab*(d1om*t+om*d1t(2))
#ifdef SECOND_DERIV
            d2t(1) = -F19*(
     &           2d0*rhob*((3d0+rhoa*rrho)*d1de+de11*rhob*rrho2)
     &           +  rhoab*((3d0+rhoa*rrho)*d2de+2d0*d1de*rhob*rrho2
     &                                         -2d0*de11*rhob*rrho3))
            d2t(2) = F19*(1d0-3d0*de-de11*rhoa*rrho
     &           - rho(n,1)*(3d0+rhoa*rrho)*d1de
     &           - rhoab*((3d0+rhoa*rrho)*d2de-d1de*(rhoa-rhob)*rrho2
     &                                        +de11*(rhoa-rhob)*rrho3))
            d2t(3) = -F19*(
     &           2d0*rhoa*((3d0+rhoa*rrho)*d1de-de11*rhoa*rrho2)
     &           +  rhoab*((3d0+rhoa*rrho)*d2de-2d0*d1de*rhoa*rrho2
     &                                         +2d0*de11*rhoa*rrho3))
     &           - 2d0
            d3fgaa(1) = -ab*(d2om*t+2d0*d1om*d1t(1)+om*d2t(1))
            d3fgaa(2) = -ab*(d2om*t+d1om*(d1t(1)+d1t(2))+om*d2t(2))
            d3fgaa(3) = -ab*(d2om*t+2d0*d1om*d1t(2)+om*d2t(3))
#endif
c
            t = F19*rhoab*de47-F43*rho2
            d1t(1) = F19*rhob*de47 - F79*rhoab*d1de 
     &             - F83*(rho(n,R_A)+rho(n,R_B))
            d1t(2) = F19*rhoa*de47 - F79*rhoab*d1de
     &             - F83*(rho(n,R_A)+rho(n,R_B))
            d1f(4) = -ab*om*t
            d2fgab(1) = -ab*(d1om*t+om*d1t(1))
            d2fgab(2) = -ab*(d1om*t+om*d1t(2))
#ifdef SECOND_DERIV
            d2t(1) = -F79*(2d0*rhob*d1de+rhoab*d2de) - F83
            d2t(2) = F19*de47 - F79*(rho(n,1)*d1de+rhoab*d2de) - F83
            d2t(3) = -F79*(2d0*rhoa*d1de+rhoab*d2de) - F83
            d3fgab(1) = -ab*(d2om*t+2d0*d1om*d1t(1)+om*d2t(1))
            d3fgab(2) = -ab*(d2om*t+d1om*(d1t(1)+d1t(2))+om*d2t(2))
            d3fgab(3) = -ab*(d2om*t+2d0*d1om*d1t(2)+om*d2t(3))
#endif
c
            t = F19*rhoab*(1d0-3d0*de-de11*rhob*rrho)-rhoa2
            d1t(1) = F19*(rhob*(1d0-3d0*de-de11*rhob*rrho)
     &             + rhoab*(-(3d0+rhob*rrho)*d1de+de11*rhob*rrho2))
     &             - 2d0*rhoa
            d1t(2) = F19*(rhoa*(1d0-3d0*de-de11*rhob*rrho)
     &             - rhoab*((3d0+rhob*rrho)*d1de+de11*rhoa*rrho2))
            d1f(5) = -ab*om*t
            d2fgbb(1) = -ab*(d1om*t+om*d1t(1))
            d2fgbb(2) = -ab*(d1om*t+om*d1t(2))
#ifdef SECOND_DERIV
            d2t(1) = -F19*(
     &           2d0*rhob*((3d0+rhob*rrho)*d1de-de11*rhob*rrho2)
     &           +  rhoab*((3d0+rhob*rrho)*d2de-2d0*d1de*rhob*rrho2
     &                                         +2d0*de11*rhob*rrho3))
     &           - 2d0
            d2t(2) = F19*(1d0-3d0*de-de11*rhob*rrho
     &           - rho(n,1)*(3d0+rhob*rrho)*d1de
     &           - rhoab*((3d0+rhob*rrho)*d2de+d1de*(rhoa-rhob)*rrho2
     &                                        -de11*(rhoa-rhob)*rrho3))
            d2t(3) = -F19*(
     &           2d0*rhoa*((3d0+rhob*rrho)*d1de+de11*rhoa*rrho2)
     &           +  rhoab*((3d0+rhob*rrho)*d2de+2d0*d1de*rhoa*rrho2
     &                                         -2d0*de11*rhoa*rrho3))
            d3fgbb(1) = -ab*(d2om*t+2d0*d1om*d1t(1)+om*d2t(1))
            d3fgbb(2) = -ab*(d2om*t+d1om*(d1t(1)+d1t(2))+om*d2t(2))
            d3fgbb(3) = -ab*(d2om*t+2d0*d1om*d1t(2)+om*d2t(3))
#endif
c
            f = c1*f1 + c2*f2 + gaa*d1f(3) + gab*d1f(4) + gbb*d1f(5)
            d1f(1) = c1*d1f1(1) + c2*d1f2(1)
     &             + gaa*d2fgaa(1) + gab*d2fgab(1) + gbb*d2fgbb(1)
            d1f(2) = c1*d1f1(2) + c2*d1f2(2)
     &             + gaa*d2fgaa(2) + gab*d2fgab(2) + gbb*d2fgbb(2)
#ifdef SECOND_DERIV
            d2f(1) = c1*d2f1(1) + c2*d2f2(1)
     &             + gaa*d3fgaa(1) + gab*d3fgab(1) + gbb*d3fgbb(1)
            d2f(2) = c1*d2f1(2) + c2*d2f2(2)
     &             + gaa*d3fgaa(2) + gab*d3fgab(2) + gbb*d3fgbb(2)
            d2f(3) = c1*d2f1(3) + c2*d2f2(3)
     &             + gaa*d3fgaa(3) + gab*d3fgab(3) + gbb*d3fgbb(3)
#endif
c
            func(n) = func(n) + f*wght
            Amat(n,D1_RA) = Amat(n,D1_RA) + d1f(1)*wght
            Amat(n,D1_RB) = Amat(n,D1_RB) + d1f(2)*wght
            Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + d1f(3)*wght
            Cmat(n,D1_GAB) = Cmat(n,D1_GAB) + d1f(4)*wght
            Cmat(n,D1_GBB) = Cmat(n,D1_GBB) + d1f(5)*wght
#ifdef SECOND_DERIV
            Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + d2f(1)*wght
            Amat2(n,D2_RA_RB) = Amat2(n,D2_RA_RB) + d2f(2)*wght
            Amat2(n,D2_RB_RB) = Amat2(n,D2_RB_RB) + d2f(3)*wght
            Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA) + d2fgaa(1)*wght
            Cmat2(n,D2_RA_GAB) = Cmat2(n,D2_RA_GAB) + d2fgab(1)*wght
            Cmat2(n,D2_RA_GBB) = Cmat2(n,D2_RA_GBB) + d2fgbb(1)*wght
            Cmat2(n,D2_RB_GAA) = Cmat2(n,D2_RB_GAA) + d2fgaa(2)*wght
            Cmat2(n,D2_RB_GAB) = Cmat2(n,D2_RB_GAB) + d2fgab(2)*wght
            Cmat2(n,D2_RB_GBB) = Cmat2(n,D2_RB_GBB) + d2fgbb(2)*wght
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
#include "nwpwxc_c_lyp.F"
#endif
C> @}
