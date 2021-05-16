#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if !defined(NWAD_PRINT)
C> \ingroup nwxc
C> @{
C>
C> \file nwxc_x_camb88.F
C> The CAM-B88 exchange functional
C>
C> @}
#endif
#endif
C> \ingroup nwxc_priv
C> @{
C>
C> \brief Evaluate the CAM-B88 exchange functional
C>
C> Evaluate the non-local part of CAM-B88 hybrid-GGA [1,2]. This routine is
C> also used to implement CAM-B3LYP.
C>
C> ### References ###
C>
C> [1] T. Yanai, D.P. Tew, N.C. Handy,
C> "A new hybrid exchange-correlation functional using the Coulomb-attenuating
C> method (CAM-B3LYP)",
C> Chem. Phys. Lett. <b>393</b>, 51-57 (2004), DOI:
C> <a href="https://doi.org/10.1016/j.cplett.2004.06.011">
C> 10.1016/j.cplett.2004.06.011</a>.
C>
C> [2] A.D. Becke,
C> "Density-functional exchange-energy approximation with correct
C> asymptotic behavior",
C> Phys. Rev. A <b>38</b>, 3098-3100 (1998), DOI:
C> <a href="https://doi.org/10.1103/PhysRevA.38.3098">
C> 10.1103/PhysRevA.38.3098</a>.
C>
#if defined(FUJITSU_VPP)
!ocl scalar
#endif
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
      Subroutine nwxc_x_camb88_p(param,tol_rho, ipol, nq, wght, 
     &                         rho, rgamma, func)
#else
      Subroutine nwxc_x_camb88(param,tol_rho, ipol, nq, wght, 
     &                         rho, rgamma, func)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      Subroutine nwxc_x_camb88_d2(param, tol_rho, ipol, nq, wght, rho,
     &                            rgamma, func)
#else
      Subroutine nwxc_x_camb88_d3(param, tol_rho, ipol, nq, wght, rho,
     &                            rgamma, func)
#endif
c
C$Id$
c
c     Coulomb attenuated Becke88 functional
c
#include "nwad.fh"
c
      implicit none
c      
#include "nwxc_param.fh"
c
c     Input and other parameters
c
#if defined(NWAD_PRINT)
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      type(nwad_dble)::param(*) !< [Input] Parameters of functional
#else
      double precision param(*) !< [Input] Parameters of functional
#endif
#else
      double precision param(*) !< [Input] Parameters of functional
                                !< - param(1): \f$ \alpha_{CAM} \f$
                                !< - param(2): \f$ \beta_{CAM} \f$
                                !< - param(3): \f$ \omega_{CAM} \f$
#endif
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
c     type(nwad_dble)::arcsinh, darcsinh, d2arcsinh
      double precision C, dg,
     &     dgdenom, t, Atmp, Ctmp
      type(nwad_dble)::rho43, gamma, x
      type(nwad_dble)::Etmp, g, gdenom
      double precision gdenom2
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c     double precision rhom23, d2g, d2gdenom
c     double precision gdenom3
c     double precision A2tmp, C2tmp, C3tmp
#endif
#if defined(THIRD_DERIV)
c     double precision rhom53, d3g, d3gdenom
c     double precision gdenom4
c     double precision A3tmp, C4tmp, C5tmp, C6tmp
#endif
c
c     arcsinh(x)=log(x+sqrt(1d0+x*x))
c     darcsinh(x)=1d0/sqrt(1d0+x*x)
c     d2arcsinh(x) = -x/sqrt(1d0+x*x)**3.0d0
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
            if (rho(n,R_T).lt.tol_rho) goto 10
c
c           Spin alpha:
c
c           rho13 = (0.5d0*rho(n,R_T))**(1.d0/3.d0)
c           rho43 = rho13**4.0d0
            rho43 = (0.5d0*rho(n,R_T))**(4.d0/3.d0)
c           gamma = delrho(n,1,1)*delrho(n,1,1) +
c    &              delrho(n,2,1)*delrho(n,2,1) +
c    &              delrho(n,3,1)*delrho(n,3,1)
            gamma = rgamma(n,G_TT)
            if (gamma.gt.tol_rho*tol_rho)then
               gamma = 0.25d0 * gamma
               x = sqrt(gamma) / rho43
            else
               x = 0.0d0
            endif
c
            gdenom = 1d0 + 6d0*BETA*x*asinh(x)
            g = -BETA*x*x / gdenom
c           gdenom2 = gdenom*gdenom
c           dgdenom = 6d0*BETA*(arcsinh(x) + x*darcsinh(x))
c           dg = BETA*x*(x*dgdenom - 2d0*gdenom) / gdenom2
c
            Etmp = 0.d0
c           Atmp = 0.d0
c           Ctmp = 0.d0
c           if (lfac) then
               Etmp = 2d0*rho43*C*wght
c              Atmp = (4d0/3d0)*rho13*C*wght
c           endif
c
c           if (nlfac) then
               Etmp = Etmp + 2d0*rho43*g*wght
c              Atmp = Atmp + (4d0/3d0)*rho13*(g-x*dg)*wght
c           endif
c
c           if (x.gt.tol_rho) then
c              Ctmp = 0.5d0 * dg / sqrt(gamma) * wght
c           endif
c
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c           A2tmp = 0d0
c           C2tmp = 0d0
c           C3tmp = 0d0
c           if(lfac) g = g + C           ! Add local contribution back to g
c           g = g + C           ! Add local contribution back to g
c           rhom23 = rho13 / (0.5d0*rho(n,R_T))
c           d2gdenom = 6d0*BETA*darcsinh(x)*(2d0 - x*x/(x*x+1d0))
c           gdenom3 = gdenom2*gdenom
c           d2g = -2d0*BETA/gdenom + 4d0*BETA*x*dgdenom/gdenom2
c    &           + BETA*x*x*d2gdenom/gdenom2
c    &           - 2d0*BETA*x*x*(dgdenom)**2/gdenom3
c
c           A2tmp = (4d0/9d0)*rhom23*(g-x*dg+4d0*x*x*d2g)*wght
c           C2tmp = - (4d0/3d0)*(rhom23**2/rho(n,R_T))*d2g*wght
c           if (x.gt.tol_rho) then
c              C3tmp = - 0.25d0*gamma**(-1.5d0)*(dg-x*d2g)*wght
c           endif
#endif
#if defined(THIRD_DERIV)
c           A3tmp = 0.0d0
c           C4tmp = 0.0d0
c           C5tmp = 0.0d0
c           C6tmp = 0.0d0
c
c           rhom53 = rhom23 / (0.5d0*rho(n,R_T))
c
c           d3gdenom = 6.0d0*BETA*
c    1          d2arcsinh(x)*( 3.0d0
c    2                       - (2.0d0*x*x - 1.0d0)/(1.0d0 + x*x))
c
c           gdenom4 = gdenom3*gdenom
c
c           d3g = 6.0d0*BETA*dgdenom/gdenom2
c    1          - 12.0d0*BETA*x*dgdenom*dgdenom/gdenom3
c    2          + 6.0d0*BETA*x*d2gdenom/gdenom2
c    3          + 6.0d0*BETA*x*x*dgdenom*dgdenom*dgdenom/gdenom4
c    4          - 6.0d0*BETA*x*x*dgdenom*d2gdenom/gdenom3
c    5          + BETA*x*x*d3gdenom/gdenom2
c
c rrr
c           A3tmp = (8.0d0/27.0d0)*rhom53*(-g + x*dg
c    1                                    - 18.0d0*x*x*d2g
c    2                                    - 8.0d0*x*x*x*d3g)*wght
c
c rrg
c           C4tmp = (2.0d0/9.0d0)*(rhom23/gamma)*( 7.0d0*x*x*d2g
c    1                                           + 4.0d0*x*x*x*d3g)*wght
c
c rgg
c           C5tmp = -(8.0d0/3.0d0)*(rhom23/rho(n,1)**3)/dsqrt(gamma)
c    1               *d3g*wght
c
c ggg

c           if (x.gt.tol_rho) then
c             C6tmp = (1.0d0/8.0d0)*gamma**(-2.5d0)*( 3.0d0*dg
c    1                                              - 3.0d0*x*d2g
c    2                                              + x*x*d3g)*wght
c           endif
#endif
#if defined(THIRD_DERIV)
            call nwxc_x_att_d3(param,tol_rho,rho(n,R_T),ipol,Etmp)
c
c           Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + A2tmp
c           Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA) + C2tmp
c           Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA) + C3tmp
c
c           Amat3(n,D3_RA_RA_RA) = Amat3(n,D3_RA_RA_RA) + A3tmp
c           Cmat3(n,D3_RA_RA_GAA) = Cmat3(n,D3_RA_RA_GAA) + C4tmp
c           Cmat3(n,D3_RA_GAA_GAA) = Cmat3(n,D3_RA_GAA_GAA) + C5tmp
c           Cmat3(n,D3_GAA_GAA_GAA) = Cmat3(n,D3_GAA_GAA_GAA) + C6tmp
#elif defined(SECOND_DERIV)
            call nwxc_x_att_d2(param,tol_rho,rho(n,R_T),ipol,Etmp)
c           Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + A2tmp
c           Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA) + C2tmp
c           Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA) + C3tmp
#else
#if defined(NWAD_PRINT)
            call nwxc_x_att_p(param,tol_rho,rho(n,R_T),ipol,Etmp)
#else
            call nwxc_x_att(param,tol_rho,rho(n,R_T),ipol,Etmp)
#endif
#endif
            func(n) = func(n) + Etmp
c           Amat(n,D1_RA) = Amat(n,D1_RA) + Atmp
c           Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + Ctmp
 10      continue
c
      else
c
c        ======> SPIN-UNRESTRICTED <======
c
         do 20 n = 1, nq
            if (rho(n,R_A).lt.0.5d0*tol_rho) goto 25
c
c           Spin alpha:
c
c           rho13 = rho(n,R_A)**(1.d0/3.d0)
c           rho43 = rho13*rho(n,R_A)
            rho43 = rho(n,R_A)**(4.d0/3.d0)
c           gamma = delrho(n,1,1)*delrho(n,1,1) +
c    &              delrho(n,2,1)*delrho(n,2,1) +
c    &              delrho(n,3,1)*delrho(n,3,1)
            gamma = rgamma(n,G_AA)
            if (gamma.gt.tol_rho*tol_rho)then
               x = sqrt(gamma) / rho43
            else
               x = 0d0
            endif
c
            gdenom = 1d0 + 6d0*BETA*x*asinh(x)
            g = -BETA*x*x / gdenom
c           dgdenom = 6d0*BETA*(arcsinh(x) + x*darcsinh(x))
c           gdenom2 = gdenom*gdenom
c           dg = BETA*x*(x*dgdenom - 2d0*gdenom) / gdenom2
c
            Etmp = 0.d0
c           Atmp = 0.d0
c           Ctmp = 0.d0
c           if (lfac) then
               Etmp = rho43*C*wght
c              Atmp = (4d0/3d0)*rho13*C*wght
c           endif
c
c           if (nlfac) then
               Etmp = Etmp + rho43*g*wght
c              Atmp = Atmp + (4d0/3d0)*rho13*(g-x*dg)*wght
c           endif
c
c           if (x.gt.tol_rho) then
c              Ctmp = 0.5d0*dg / sqrt(gamma) * wght
c           endif
c
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c           if (lfac) g = g + C           ! Add local contribution back to g
c           g = g + C           ! Add local contribution back to g
c           rhom23 = rho13 / rho(n,R_A)
c           d2gdenom = 6d0*BETA*darcsinh(x)*(2d0 - x*x/(x*x+1d0))
c           gdenom3 = gdenom2*gdenom
c           d2g = -2d0*BETA/gdenom + 4d0*BETA*x*dgdenom/gdenom2
c    &           + BETA*x*x*d2gdenom/gdenom2
c    &           - 2d0*BETA*x*x*(dgdenom)**2/gdenom3
c
c           A2tmp = (4d0/9d0)*rhom23*(g-x*dg+4d0*x*x*d2g)*wght
c           C2tmp =  (2d0/3d0)*(rhom23**2/rho(n,R_A))*d2g*wght
c           if (x.gt.tol_rho) then
c              C3tmp = -0.25d0*gamma**(-1.5d0)*(dg-x*d2g)*wght
c           endif
#endif
#if defined(THIRD_DERIV)
c           rhom53 = rhom23 / rho(n,R_A)
c
c           d3gdenom = 6.0d0*BETA*
c    1          d2arcsinh(x)*( 3.0d0
c    2                       - (2.0d0*x*x - 1.0d0)/(1.0d0 + x*x))
c
c           gdenom4 = gdenom3*gdenom
c
c           d3g = 6.0d0*BETA*dgdenom/gdenom2
c    1          - 12.0d0*BETA*x*dgdenom*dgdenom/gdenom3
c    2          + 6.0d0*BETA*x*d2gdenom/gdenom2
c    3          + 6.0d0*BETA*x*x*dgdenom*dgdenom*dgdenom/gdenom4
c    4          - 6.0d0*BETA*x*x*dgdenom*d2gdenom/gdenom3
c    5          + BETA*x*x*d3gdenom/gdenom2
c
c           A3tmp = (8.0d0/27.0d0)*rhom53*(-g + x*dg
c    1                                    - 18.0d0*x*x*d2g
c    2                                    - 8.0d0*x*x*x*d3g)*wght
c
c           C4tmp = (2.0d0/9.0d0)*(rhom23/gamma)*( 7.0d0*x*x*d2g
c    1                                           + 4.0d0*x*x*x*d3g)*wght
c
c           C5tmp = -(1.0d0/3.0d0)*(rhom23/rho(n,2)**3)/dsqrt(gamma)
c    1               *d3g*wght
c
c           if (x.gt.tol_rho) then
c             C6tmp = (1.0d0/8.0d0)*gamma**(-2.5d0)*( 3.0d0*dg
c    1                                              - 3.0d0*x*d2g
c    2                                              + x*x*d3g)*wght
c           endif
#endif
#if defined(THIRD_DERIV)
            call nwxc_x_att_d3(param,tol_rho,rho(n,R_A),ipol,Etmp)
c
c           Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + A2tmp
c           Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA) + C2tmp
c           Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA) + C3tmp
c
c           Amat3(n,D3_RA_RA_RA) = Amat3(n,D3_RA_RA_RA) + A3tmp
c           Cmat3(n,D3_RA_RA_GAA) = Cmat3(n,D3_RA_RA_GAA) + C4tmp
c           Cmat3(n,D3_RA_GAA_GAA) = Cmat3(n,D3_RA_GAA_GAA) + C5tmp
c           Cmat3(n,D3_GAA_GAA_GAA) = Cmat3(n,D3_GAA_GAA_GAA) + C6tmp
#elif defined(SECOND_DERIV)
            call nwxc_x_att_d2(param,tol_rho,rho(n,R_A),ipol,Etmp)
c           Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + A2tmp
c           Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA) + C2tmp
c           Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA) + C3tmp
#else
#if defined(NWAD_PRINT)
            call nwxc_x_att_p(param,tol_rho,rho(n,R_A),ipol,Etmp)
#else
            call nwxc_x_att(param,tol_rho,rho(n,R_A),ipol,Etmp)
#endif
#endif
            func(n) = func(n) + Etmp
c           Amat(n,D1_RA) = Amat(n,D1_RA) + Atmp
c           Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + Ctmp
c
 25         continue
c
c           Spin beta:
c
            if (rho(n,R_B).lt.0.5d0*tol_rho) goto 20
c
c           rho13 = rho(n,R_B)**(1.d0/3.d0)
c           rho43 = rho13*rho(n,R_B)
            rho43 = rho(n,R_B)**(4.d0/3.d0)
c           gamma = delrho(n,1,2)*delrho(n,1,2) +
c    &              delrho(n,2,2)*delrho(n,2,2) +
c    &              delrho(n,3,2)*delrho(n,3,2)
            gamma = rgamma(n,G_BB)
            if (gamma.gt.tol_rho*tol_rho)then
               x = sqrt(gamma) / rho43
            else
               x = 0d0
            endif
c
            gdenom = 1d0 + 6d0*BETA*x*asinh(x)
            g = -BETA*x*x / gdenom
c           dgdenom = 6d0*BETA*(arcsinh(x) + x*darcsinh(x))
c           gdenom2 = gdenom*gdenom
c           dg = BETA*x*(x*dgdenom - 2d0*gdenom) / gdenom2
c
            Etmp = 0.d0
c           Atmp = 0.d0
c           Ctmp = 0.d0
c           if (lfac) then
               Etmp = rho43*C*wght
c              Atmp = (4d0/3d0)*rho13*C*wght
c           endif
c
c           if (nlfac) then
               Etmp = Etmp + rho43*g*wght
c              Atmp = Atmp + (4d0/3d0)*rho13*(g-x*dg)*wght
c           endif
c
c           if (x.gt.tol_rho) then
c              Ctmp = 0.5d0*dg / sqrt(gamma) * wght
c           endif
c
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c           A2tmp = 0d0
c           C2tmp = 0d0
c           C3tmp = 0d0
c           if(lfac) g = g + C           ! Add local contribution back to g
c           g = g + C           ! Add local contribution back to g
c           rhom23 = rho13 / rho(n,R_B)
c           d2gdenom = 6d0*BETA*darcsinh(x)*(2d0 - x*x/(x*x+1d0))
c           gdenom3 = gdenom2*gdenom
c           d2g = -2d0*BETA/gdenom + 4d0*BETA*x*dgdenom/gdenom2
c    &           + BETA*x*x*d2gdenom/gdenom2
c    &           - 2d0*BETA*x*x*(dgdenom)**2/gdenom3
c
c           A2tmp = (4d0/9d0)*rhom23*(g-x*dg+4d0*x*x*d2g)*wght
c           C2tmp = -(2d0/3d0)*(rhom23**2/rho(n,R_B))*d2g*wght
c           if (x.gt.tol_rho) then
c              C3tmp = - 0.25d0*gamma**(-1.5d0)*(dg-x*d2g)*wght
c           endif
#endif
#if defined(THIRD_DERIV)
c           rhom53 = rhom23 / rho(n,R_B)
c
c           d3gdenom = 6.0d0*BETA*
c    1          d2arcsinh(x)*( 3.0d0
c    2                       - (2.0d0*x*x - 1.0d0)/(1.0d0 + x*x))
c
c           gdenom4 = gdenom3*gdenom
c
c           d3g = 6.0d0*BETA*dgdenom/gdenom2
c    1          - 12.0d0*BETA*x*dgdenom*dgdenom/gdenom3
c    2          + 6.0d0*BETA*x*d2gdenom/gdenom2
c    3          + 6.0d0*BETA*x*x*dgdenom*dgdenom*dgdenom/gdenom4
c    4          - 6.0d0*BETA*x*x*dgdenom*d2gdenom/gdenom3
c    5          + BETA*x*x*d3gdenom/gdenom2
c
c           A3tmp = (8.0d0/27.0d0)*rhom53*(-g + x*dg
c    1                                    - 18.0d0*x*x*d2g
c    2                                    - 8.0d0*x*x*x*d3g)*wght
c
c           C4tmp = (2.0d0/9.0d0)*(rhom23/gamma)*( 7.0d0*x*x*d2g
c    1                                           + 4.0d0*x*x*x*d3g)*wght
c
c           C5tmp = -(1.0d0/3.0d0)*(rhom23/rho(n,3)**3)/dsqrt(gamma)
c    1               *d3g*wght
c
c           if (x.gt.tol_rho) then
c             C6tmp = (1.0d0/8.0d0)*gamma**(-2.5d0)*( 3.0d0*dg
c    1                                              - 3.0d0*x*d2g
c    2                                              + x*x*d3g)*wght
c           endif

#endif
#if defined(THIRD_DERIV)
            call nwxc_x_att_d3(param,tol_rho,rho(n,R_B),ipol,Etmp)
c
c           Amat2(n,D2_RB_RB) = Amat2(n,D2_RB_RB) + A2tmp
c           Cmat2(n,D2_RB_GBB) = Cmat2(n,D2_RB_GBB) + C2tmp
c           Cmat2(n,D2_GBB_GBB) = Cmat2(n,D2_GBB_GBB) + C3tmp
c
c           Amat3(n,D3_RB_RB_RB) = Amat3(n,D3_RB_RB_RB) + A3tmp
c           Cmat3(n,D3_RB_RB_GBB) = Cmat3(n,D3_RB_RB_GBB) + C4tmp
c           Cmat3(n,D3_RB_GBB_GBB) = Cmat3(n,D3_RB_GBB_GBB) + C5tmp
c           Cmat3(n,D3_GBB_GBB_GBB) = Cmat3(n,D3_GBB_GBB_GBB) + C6tmp
#elif defined(SECOND_DERIV)
            call nwxc_x_att_d2(param,tol_rho,rho(n,R_B),ipol,Etmp)
c           Amat2(n,D2_RB_RB) = Amat2(n,D2_RB_RB) + A2tmp
c           Cmat2(n,D2_RB_GBB) = Cmat2(n,D2_RB_GBB) + C2tmp
c           Cmat2(n,D2_GBB_GBB) = Cmat2(n,D2_GBB_GBB) + C3tmp
#else
#if defined(NWAD_PRINT)
            call nwxc_x_att_p(param,tol_rho,rho(n,R_B),ipol,Etmp)
#else
            call nwxc_x_att(param,tol_rho,rho(n,R_B),ipol,Etmp)
#endif
#endif
            func(n) = func(n) + Etmp
c           Amat(n,D1_RB) = Amat(n,D1_RB) + Atmp
c           Cmat(n,D1_GBB) = Cmat(n,D1_GBB) + Ctmp
 20      continue
c
      endif
c
      return
      end
#ifndef NWAD_PRINT
#define NWAD_PRINT
c
c     Compile source again for Maxima
c
#include "nwxc_x_camb88.F"
#endif
#ifndef SECOND_DERIV
#define SECOND_DERIV
c
c     Compile source again for the 2nd derivative case
c
#include "nwxc_x_camb88.F"
#endif
#ifndef THIRD_DERIV
#define THIRD_DERIV
c
c     Compile source again for the 3rd derivative case
c
#include "nwxc_x_camb88.F"
#endif
#undef NWAD_PRINT
C>
C> @}
