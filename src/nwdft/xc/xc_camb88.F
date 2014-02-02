#if defined(FUJITSU_VPP)
!ocl scalar
#endif
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      Subroutine xc_camb88(tol_rho, fac, lfac, nlfac, rho, delrho, 
     &                      Amat, Cmat, nq, ipol, Ex, qwght,ldew,func)
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      Subroutine xc_camb88_d2(tol_rho, fac, lfac, nlfac, rho, delrho, 
     &                         Amat, Amat2, Cmat, Cmat2, nq, ipol, Ex,
     &                         qwght,ldew,func)
#else
      Subroutine xc_camb88_d3(tol_rho, fac, lfac, nlfac, rho, delrho,
     1                         Amat, Amat2, Amat3, Cmat, Cmat2, Cmat3, 
     2                         nq, ipol, Ex, qwght,ldew,func)
#endif
c
C$Id$
c
c     Coulomb attenuated Becke88 functional
c
      implicit none
c      
#include "dft2drv.fh"
#include "dft3drv.fh"
c
      double precision tol_rho, fac, Ex
      integer nq, ipol
      logical lfac, nlfac,ldew
      double precision func(*)  ! value of the functional [output]
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
c     Sampling Matrices for the XC Potential
c
      double precision Amat(nq,ipol), Cmat(nq,*)
c
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c
c     Second Derivatives of the Exchange Energy Functional
c
      double precision Amat2(nq,NCOL_AMAT2), Cmat2(nq,NCOL_CMAT2)
      double precision A2tmp, C2tmp, C3tmp
#endif
#ifdef THIRD_DERIV
c
c     Third Derivatives of the Exchange Energy Functional
c
      double precision Amat3(nq,NCOL_AMAT3), Cmat3(nq,NCOL_CMAT3)
      double precision A3tmp, C4tmp, C5tmp, C6tmp
#endif
c
      double precision BETA
      Parameter (BETA = 0.0042D0)
c
c References:
c
c    Becke, Phys. Rev. A 38, 3098 (1988)
c    Johnson, Gill & Pople, J. Chem. Phys. 98, 5612 (1993)
c
c***************************************************************************
c
      integer n
      double precision arcsinh, darcsinh, d2arcsinh
      double precision C, rho13, rho43, gamma, x, g, gdenom, dg,
     &     dgdenom, t, Etmp, Atmp, Ctmp
      double precision gdenom2
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
      double precision rhom23, d2g, d2gdenom
      double precision gdenom3
#endif
c
#ifdef THIRD_DERIV
      double precision rhom53, d3g, d3gdenom
      double precision gdenom4
#endif
c
      arcsinh(x)=log(x+dsqrt(1d0+x*x))
      darcsinh(x)=1d0/dsqrt(1d0+x*x)
      d2arcsinh(x) = -x/dsqrt(1d0+x*x)**3
c
c     Uniform electron gas constant
c
      C = -(1.5d0)*(0.75d0/acos(-1d0))**(1d0/3d0)
c
      if (ipol.eq.1) then
c
c        ======> SPIN-RESTRICTED <======
c
         do 10 n = 1, nq
            if (rho(n,1).lt.tol_rho) goto 10
c
c           Spin alpha:
c
            rho13 = (0.5d0*rho(n,1))**(1.d0/3.d0)
            rho43 = rho13**4 
            gamma = delrho(n,1,1)*delrho(n,1,1) +
     &              delrho(n,2,1)*delrho(n,2,1) +
     &              delrho(n,3,1)*delrho(n,3,1)
            if (dsqrt(gamma).gt.tol_rho)then
               gamma = 0.25d0 * gamma
               x = dsqrt(gamma) / rho43
            else
               x = 0d0
            endif
c
            gdenom = 1d0 + 6d0*BETA*x*arcsinh(x)
            gdenom2 = gdenom*gdenom
            g = -BETA*x*x / gdenom
            dgdenom = 6d0*BETA*(arcsinh(x) + x*darcsinh(x))
c            dg = BETA*x*(x*dgdenom - 2d0*gdenom) / gdenom**2
            dg = BETA*x*(x*dgdenom - 2d0*gdenom) / gdenom2
c
            Etmp = 0.d0
            Atmp = 0.d0
            Ctmp = 0.d0
            if (lfac) then
               Etmp = 2d0*rho43*C*fac
               Atmp = (4d0/3d0)*rho13*C*fac
            endif
c
            if (nlfac) then
               Etmp = Etmp + 2d0*rho43*g*fac
               Atmp = Atmp + (4d0/3d0)*rho13*(g-x*dg)*fac
            endif
c
            if (x.gt.tol_rho) then
               Ctmp = 0.5d0 * dg / sqrt(gamma) * fac
            endif
c
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
            A2tmp = 0d0
            C2tmp = 0d0
            C3tmp = 0d0
            if(lfac) g = g + C           ! Add local contribution back to g
            rhom23 = rho13 / (0.5d0*rho(n,1))
            d2gdenom = 6d0*BETA*darcsinh(x)*(2d0 - x*x/(x*x+1d0))
            gdenom3 = gdenom2*gdenom
            d2g = -2d0*BETA/gdenom + 4d0*BETA*x*dgdenom/gdenom2
     &           + BETA*x*x*d2gdenom/gdenom2
     &           - 2d0*BETA*x*x*(dgdenom)**2/gdenom3
c
c rr
            A2tmp = (4d0/9d0)*rhom23*(g-x*dg+4d0*x*x*d2g)*fac
c rg
            C2tmp = - (4d0/3d0)*(rhom23**2/rho(n,1))*d2g*fac
c gg
            if (x.gt.tol_rho) then
               C3tmp = - 0.25d0*gamma**(-1.5d0)*(dg-x*d2g)*fac
            endif
#endif
#ifdef THIRD_DERIV
            A3tmp = 0.0d0
            C4tmp = 0.0d0
            C5tmp = 0.0d0
            C6tmp = 0.0d0
c
            rhom53 = rhom23 / (0.5d0*rho(n,1))
c
            d3gdenom = 6.0d0*BETA*
     1          d2arcsinh(x)*( 3.0d0
     2                       - (2.0d0*x*x - 1.0d0)/(1.0d0 + x*x))
c
            gdenom4 = gdenom3*gdenom
c
            d3g = 6.0d0*BETA*dgdenom/gdenom2
     1          - 12.0d0*BETA*x*dgdenom*dgdenom/gdenom3
     2          + 6.0d0*BETA*x*d2gdenom/gdenom2
     3          + 6.0d0*BETA*x*x*dgdenom*dgdenom*dgdenom/gdenom4
     4          - 6.0d0*BETA*x*x*dgdenom*d2gdenom/gdenom3
     5          + BETA*x*x*d3gdenom/gdenom2
c
c rrr
            A3tmp = (8.0d0/27.0d0)*rhom53*(-g + x*dg
     1                                    - 18.0d0*x*x*d2g
     2                                    - 8.0d0*x*x*x*d3g)*fac
c
c rrg
            C4tmp = (2.0d0/9.0d0)*(rhom23/gamma)*( 7.0d0*x*x*d2g
     1                                           + 4.0d0*x*x*x*d3g)*fac
c
c rgg
            C5tmp = -(8.0d0/3.0d0)*(rhom23/rho(n,1)**3)/dsqrt(gamma)
     1               *d3g*fac
c
c ggg
            if (x.gt.tol_rho) then
              C6tmp = (1.0d0/8.0d0)*gamma**(-2.5d0)*( 3.0d0*dg
     1                                              - 3.0d0*x*d2g
     2                                              + x*x*d3g)*fac
            endif
#endif
c
#ifdef THIRD_DERIV
            call xc_att_xc_d3(rho(n,1),ipol,Etmp,Atmp,Ctmp,A2tmp,
     &           C2tmp,C3tmp,A3tmp,C4tmp,C5tmp,C6tmp)
c
            Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + A2tmp
            Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA) + C2tmp
            Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA) + C3tmp
c
            Amat3(n,D3_RA_RA_RA) = Amat3(n,D3_RA_RA_RA) + A3tmp
            Cmat3(n,D3_RA_RA_GAA) = Cmat3(n,D3_RA_RA_GAA) + C4tmp
            Cmat3(n,D3_RA_GAA_GAA) = Cmat3(n,D3_RA_GAA_GAA) + C5tmp
            Cmat3(n,D3_GAA_GAA_GAA) = Cmat3(n,D3_GAA_GAA_GAA) + C6tmp
#elif defined(SECOND_DERIV)
            call xc_att_xc_d2(rho(n,1),ipol,Etmp,Atmp,Ctmp,A2tmp,
     &           C2tmp,C3tmp)
c
            Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + A2tmp
            Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA) + C2tmp
            Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA) + C3tmp
#else
            call xc_att_xc(rho(n,1),ipol,Etmp,Atmp,Ctmp)
#endif
            Ex = Ex + qwght(n)*Etmp
            if(ldew) func(n) = func(n) + Etmp
            Amat(n,1) = Amat(n,1) + Atmp
            Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + Ctmp
 10      continue
c
      else
c
c        ======> SPIN-UNRESTRICTED <======
c
         do 20 n = 1, nq
            if (rho(n,1).lt.tol_rho) goto 20
            if (rho(n,2).lt.tol_rho) goto 25
c
c           Spin alpha:
c
            rho13 = rho(n,2)**(1.d0/3.d0)
            rho43 = rho13*rho(n,2)
            gamma = delrho(n,1,1)*delrho(n,1,1) +
     &              delrho(n,2,1)*delrho(n,2,1) +
     &              delrho(n,3,1)*delrho(n,3,1)
            if (dsqrt(gamma).gt.tol_rho)then
               x = dsqrt(gamma) / rho43
            else
               x = 0d0
            endif
c
            gdenom = 1d0 + 6d0*BETA*x*arcsinh(x)
            g = -BETA*x*x / gdenom
            dgdenom = 6d0*BETA*(arcsinh(x) + x*darcsinh(x))
            gdenom2 = gdenom*gdenom
            dg = BETA*x*(x*dgdenom - 2d0*gdenom) / gdenom2
c
            Etmp = 0.d0
            Atmp = 0.d0
            Ctmp = 0.d0
            if (lfac) then
               Etmp = rho43*C*fac
               Atmp = (4d0/3d0)*rho13*C*fac
            endif
c
            if (nlfac) then
               Etmp = Etmp + rho43*g*fac
               Atmp = Atmp + (4d0/3d0)*rho13*(g-x*dg)*fac
            endif
c
            if (x.gt.tol_rho) then
               Ctmp = 0.5d0*dg / sqrt(gamma) * fac
            endif
c
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
            if (lfac) g = g + C           ! Add local contribution back to g
            rhom23 = rho13 / rho(n,2)
            d2gdenom = 6d0*BETA*darcsinh(x)*(2d0 - x*x/(x*x+1d0))
            gdenom3 = gdenom2*gdenom
            d2g = -2d0*BETA/gdenom + 4d0*BETA*x*dgdenom/gdenom2
     &           + BETA*x*x*d2gdenom/gdenom2
     &           - 2d0*BETA*x*x*(dgdenom)**2/gdenom3
c            d2g = -2d0*BETA/gdenom + 4d0*BETA*x*dgdenom/gdenom**2
c     &           + BETA*x*x*d2gdenom/gdenom**2
c     &           - 2d0*BETA*x*x*(dgdenom)**2/gdenom**3
c
            A2tmp = (4d0/9d0)*rhom23*(g-x*dg+4d0*x*x*d2g)*fac
c            C2tmp =  (2d0/3d0)*(rhom23**2/rho(n,2))*d2g*fac
            C2tmp = -(2d0/3d0)*(rhom23**2/rho(n,2))*d2g*fac
            if (x.gt.tol_rho) then
               C3tmp = -0.25d0*gamma**(-1.5d0)*(dg-x*d2g)*fac
            endif
#endif
#ifdef THIRD_DERIV
            rhom53 = rhom23 / rho(n,2)
c
            d3gdenom = 6.0d0*BETA*
     1          d2arcsinh(x)*( 3.0d0
     2                       - (2.0d0*x*x - 1.0d0)/(1.0d0 + x*x))
c
            gdenom4 = gdenom3*gdenom
c
            d3g = 6.0d0*BETA*dgdenom/gdenom2
     1          - 12.0d0*BETA*x*dgdenom*dgdenom/gdenom3
     2          + 6.0d0*BETA*x*d2gdenom/gdenom2
     3          + 6.0d0*BETA*x*x*dgdenom*dgdenom*dgdenom/gdenom4
     4          - 6.0d0*BETA*x*x*dgdenom*d2gdenom/gdenom3
     5          + BETA*x*x*d3gdenom/gdenom2
c
            A3tmp = (8.0d0/27.0d0)*rhom53*(-g + x*dg
     1                                    - 18.0d0*x*x*d2g
     2                                    - 8.0d0*x*x*x*d3g)*fac
c
            C4tmp = (2.0d0/9.0d0)*(rhom23/gamma)*( 7.0d0*x*x*d2g
     1                                           + 4.0d0*x*x*x*d3g)*fac
c
            C5tmp = -(1.0d0/3.0d0)*(rhom23/rho(n,2)**3)/dsqrt(gamma)
     1               *d3g*fac
c
            if (x.gt.tol_rho) then
              C6tmp = (1.0d0/8.0d0)*gamma**(-2.5d0)*( 3.0d0*dg
     1                                              - 3.0d0*x*d2g
     2                                              + x*x*d3g)*fac
            endif
#endif
c
#ifdef THIRD_DERIV
            call xc_att_xc_d3(rho(n,2),ipol,Etmp,Atmp,Ctmp,A2tmp,
     &           C2tmp,C3tmp,A3tmp,C4tmp,C5tmp,C6tmp)
c
            Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + A2tmp
            Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA) + C2tmp
            Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA) + C3tmp
c
            Amat3(n,D3_RA_RA_RA) = Amat3(n,D3_RA_RA_RA) + A3tmp
            Cmat3(n,D3_RA_RA_GAA) = Cmat3(n,D3_RA_RA_GAA) + C4tmp
            Cmat3(n,D3_RA_GAA_GAA) = Cmat3(n,D3_RA_GAA_GAA) + C5tmp
            Cmat3(n,D3_GAA_GAA_GAA) = Cmat3(n,D3_GAA_GAA_GAA) + C6tmp
#elif defined(SECOND_DERIV)
            call xc_att_xc_d2(rho(n,2),ipol,Etmp,Atmp,Ctmp,A2tmp,C2tmp,
     &           C3tmp)
            Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + A2tmp
            Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA) + C2tmp
            Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA) + C3tmp
#else
           call xc_att_xc(rho(n,2),ipol,Etmp,Atmp,Ctmp)
#endif
            Ex = Ex + qwght(n)*Etmp
            if(ldew) func(n) = func(n) + Etmp
            Amat(n,1) = Amat(n,1) + Atmp
            Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + Ctmp
c
 25         continue
c
c           Spin beta:
c
            if (rho(n,3).lt.tol_rho) goto 20
c
            rho13 = rho(n,3)**(1.d0/3.d0)
            rho43 = rho13*rho(n,3)
            gamma = delrho(n,1,2)*delrho(n,1,2) +
     &              delrho(n,2,2)*delrho(n,2,2) +
     &              delrho(n,3,2)*delrho(n,3,2)
            if (dsqrt(gamma).gt.tol_rho)then
               x = dsqrt(gamma) / rho43
            else
               x = 0d0
            endif
c
            gdenom = 1d0 + 6d0*BETA*x*arcsinh(x)
            g = -BETA*x*x / gdenom
            dgdenom = 6d0*BETA*(arcsinh(x) + x*darcsinh(x))
            gdenom2 = gdenom*gdenom
            dg = BETA*x*(x*dgdenom - 2d0*gdenom) / gdenom2
c
            Etmp = 0.d0
            Atmp = 0.d0
            Ctmp = 0.d0
            if (lfac) then
               Etmp = rho43*C*fac
               Atmp = (4d0/3d0)*rho13*C*fac
            endif
c
            if (nlfac) then
               Etmp = Etmp + rho43*g*fac
               Atmp = Atmp + (4d0/3d0)*rho13*(g-x*dg)*fac
            endif
c
            if (x.gt.tol_rho) then
               Ctmp = 0.5d0*dg / sqrt(gamma) * fac
            endif
c
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
            A2tmp = 0d0
            C2tmp = 0d0
            C3tmp = 0d0
            if(lfac) g = g + C           ! Add local contribution back to g
            rhom23 = rho13 / rho(n,3)
            d2gdenom = 6d0*BETA*darcsinh(x)*(2d0 - x*x/(x*x+1d0))
            gdenom3 = gdenom2*gdenom
            d2g = -2d0*BETA/gdenom + 4d0*BETA*x*dgdenom/gdenom2
     &           + BETA*x*x*d2gdenom/gdenom2
     &           - 2d0*BETA*x*x*(dgdenom)**2/gdenom3
            A2tmp = (4d0/9d0)*rhom23*(g-x*dg+4d0*x*x*d2g)*fac
            C2tmp = -(2d0/3d0)*(rhom23**2/rho(n,3))*d2g*fac
            if (x.gt.tol_rho) then
               C3tmp = - 0.25d0*gamma**(-1.5d0)*(dg-x*d2g)*fac
            endif
#endif
#ifdef THIRD_DERIV
            rhom53 = rhom23 / rho(n,3)
c
            d3gdenom = 6.0d0*BETA*
     1          d2arcsinh(x)*( 3.0d0
     2                       - (2.0d0*x*x - 1.0d0)/(1.0d0 + x*x))
c
            gdenom4 = gdenom3*gdenom
c
            d3g = 6.0d0*BETA*dgdenom/gdenom2
     1          - 12.0d0*BETA*x*dgdenom*dgdenom/gdenom3
     2          + 6.0d0*BETA*x*d2gdenom/gdenom2
     3          + 6.0d0*BETA*x*x*dgdenom*dgdenom*dgdenom/gdenom4
     4          - 6.0d0*BETA*x*x*dgdenom*d2gdenom/gdenom3
     5          + BETA*x*x*d3gdenom/gdenom2
c
            A3tmp = (8.0d0/27.0d0)*rhom53*(-g + x*dg
     1                                    - 18.0d0*x*x*d2g
     2                                    - 8.0d0*x*x*x*d3g)*fac
c
            C4tmp = (2.0d0/9.0d0)*(rhom23/gamma)*( 7.0d0*x*x*d2g
     1                                           + 4.0d0*x*x*x*d3g)*fac
c
            C5tmp = -(1.0d0/3.0d0)*(rhom23/rho(n,3)**3)/dsqrt(gamma)
     1               *d3g*fac
c
            if (x.gt.tol_rho) then
              C6tmp = (1.0d0/8.0d0)*gamma**(-2.5d0)*( 3.0d0*dg
     1                                              - 3.0d0*x*d2g
     2                                              + x*x*d3g)*fac
            endif
#endif
c
#ifdef THIRD_DERIV
            call xc_att_xc_d3(rho(n,3),ipol,Etmp,Atmp,Ctmp,A2tmp,
     &           C2tmp,C3tmp,A3tmp,C4tmp,C5tmp,C6tmp)
c
            Amat2(n,D2_RB_RB) = Amat2(n,D2_RB_RB) + A2tmp
            Cmat2(n,D2_RB_GBB) = Cmat2(n,D2_RB_GBB) + C2tmp
            Cmat2(n,D2_GBB_GBB) = Cmat2(n,D2_GBB_GBB) + C3tmp
c
            Amat3(n,D3_RB_RB_RB) = Amat3(n,D3_RB_RB_RB) + A3tmp
            Cmat3(n,D3_RB_RB_GBB) = Cmat3(n,D3_RB_RB_GBB) + C4tmp
            Cmat3(n,D3_RB_GBB_GBB) = Cmat3(n,D3_RB_GBB_GBB) + C5tmp
            Cmat3(n,D3_GBB_GBB_GBB) = Cmat3(n,D3_GBB_GBB_GBB) + C6tmp
#elif defined(SECOND_DERIV)
            call xc_att_xc_d2(rho(n,3),ipol,Etmp,Atmp,Ctmp,A2tmp,C2tmp,
     &           C3tmp)
            Amat2(n,D2_RB_RB) = Amat2(n,D2_RB_RB) + A2tmp
            Cmat2(n,D2_RB_GBB) = Cmat2(n,D2_RB_GBB) + C2tmp
            Cmat2(n,D2_GBB_GBB) = Cmat2(n,D2_GBB_GBB) + C3tmp
#else
            call xc_att_xc(rho(n,3),ipol,Etmp,Atmp,Ctmp)
#endif
            Ex = Ex + qwght(n)*Etmp
            if(ldew) func(n) = func(n) + Etmp
            Amat(n,2) = Amat(n,2) + Atmp
            Cmat(n,D1_GBB) = Cmat(n,D1_GBB) + Ctmp
 20      continue
c
      endif
c
      return
      end
#ifndef SECOND_DERIV
#define SECOND_DERIV
c
c     Compile source again for the 2nd derivative case
c
#include "xc_camb88.F"
#endif
c
#ifndef THIRD_DERIV
#define THIRD_DERIV
c
c     Compile source again for the 3rd derivative case
c
#include "xc_camb88.F"
#endif
