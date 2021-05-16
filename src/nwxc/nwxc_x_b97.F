#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if !defined(NWAD_PRINT)
C> \ingroup nwxc
C> @{
C>
C> \file nwxc_x_b97.F
C> The B97 family of exchange functionals
C>
C> @}
#include "nwxcP_xc_b97.fh"
#endif
#endif
C>
C> \ingroup nwxc_priv
C> @{
C>
C> \brief Evaluate the B97 family of exchange functionals
C>
C> This code evaluates exchange functionals from the 
C> B97 family of functionals [1,2].
C>
C> ### References ###
C>
C> [1] A.D. Becke, "Density-functional thermochemistry. V. Systematic
C>     optimization of exchange-correlation functionals", J. Chem. Phys.
C>     107 (1997) 8554-8560, DOI:
C>     <a href="https://doi.org/10.1063/1.475007">
C>     10.1063/1.475007</a>.
C>
C> [2] S. Grimme, "Semiempirical GGA-type density functional constructed
C>     with a long-range dispersion correction", J. Comput. Chem. 27 (2006)
C>     1787-1799, DOI:
C>     <a href="https://doi.org/10.1002/jcc.20495">
C>     10.1002/jcc.20495</a>.
C>
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
      Subroutine nwxc_x_b97_p(param,tol_rho,ipol,nq,wght,rho,rgamma,
     &                      func)
#else
      Subroutine nwxc_x_b97(param,tol_rho,ipol,nq,wght,rho,rgamma,
     &                      func)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      Subroutine nwxc_x_b97_d2(param,tol_rho,ipol,nq,wght,rho,rgamma,
     &                      func)
#else
      Subroutine nwxc_x_b97_d3(param,tol_rho,ipol,nq,wght,rho,rgamma,
     &                      func)
#endif
c
c     $Id$
c
#include "nwad.fh"
      implicit none
c
#include "nwxc_param.fh"
c
c     Input and other parameters
c
#if defined(NWAD_PRINT)
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      type(nwad_dble)::param(*)
#else
      double precision param(*)
#endif
#else
      double precision param(*)!< [Input] Parameters of functional as
                               !< defined in [1]:
                               !< - param(1): \f$m\f$ of Eqs.(20).
                               !< - param(2): \f$C_{X\sigma,0}\f$
                               !< - param(3): \f$C_{X\sigma,1}\f$
                               !< - param(4): \f$C_{X\sigma,2}\f$
                               !< - param(5): \f$C_{X\sigma,3}\f$
                               !< - param(6): \f$C_{X\sigma,4}\f$
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
      type(nwad_dble)::func(nq)     !< [Output] The value of the functional
c     double precision Amat(nq,*)   !< [Output] The derivative wrt rho
c     double precision Cmat(nq,*)   !< [Output] The derivative wrt rgamma
c
c     Local variables
c
      integer i
      type(nwad_dble)::rho_a(0:1)
      type(nwad_dble)::rho_b(0:1)
      type(nwad_dble)::FX(0:_FXC_NUMDERI)
c
c     Code
c
      if (ipol.eq.1) then
        do i = 1, nq
          rho_a(0) = rho(i,R_T)*0.5d0
          rho_b(0) = rho_a(0)
          rho_a(1) = rgamma(i,G_TT)*0.25d0
          rho_b(1) = rho_a(1)
          if (rho_a(0).gt.tol_rho) then
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
            call nwxcp_x_b97_p(rho_a,rho_b,1.0d0,tol_rho,FX,param)
#else
            call nwxcp_x_b97(rho_a,rho_b,1.0d0,tol_rho,FX,param)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
            call nwxcp_x_b97_d2(rho_a,rho_b,1.0d0,tol_rho,FX,param)
#else
            call nwxcp_x_b97_d3(rho_a,rho_b,1.0d0,tol_rho,FX,param)
#endif
            func(i)        = func(i)        + FX(_FXC_E)*wght
!           Amat(i,D1_RA)  = Amat(i,D1_RA)  + FX(_FXC_RA)*wght
!           Cmat(i,D1_GAA) = Cmat(i,D1_GAA) + FX(_FXC_GAA)*wght
          endif
        enddo
      else
        do i = 1, nq
          rho_a(0) = rho(i,R_A)
          rho_b(0) = rho(i,R_B)
          rho_a(1) = rgamma(i,G_AA)
          rho_b(1) = rgamma(i,G_BB)
          if (rho_a(0).gt.tol_rho.or.rho_b(0).gt.tol_rho) then
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
            call nwxcp_x_b97_p(rho_a,rho_b,1.0d0,tol_rho,FX,param)
#else
            call nwxcp_x_b97(rho_a,rho_b,1.0d0,tol_rho,FX,param)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
            call nwxcp_x_b97_d2(rho_a,rho_b,1.0d0,tol_rho,FX,param)
#else
            call nwxcp_x_b97_d3(rho_a,rho_b,1.0d0,tol_rho,FX,param)
#endif
            func(i)        = func(i)        + FX(_FXC_E)*wght
!           Amat(i,D1_RA)  = Amat(i,D1_RA)  + FX(_FXC_RA)*wght
!           Amat(i,D1_RB)  = Amat(i,D1_RB)  + FX(_FXC_RB)*wght
!           Cmat(i,D1_GAA) = Cmat(i,D1_GAA) + FX(_FXC_GAA)*wght
!           Cmat(i,D1_GBB) = Cmat(i,D1_GBB) + FX(_FXC_GBB)*wght
          endif
        enddo
      endif
c
      end
c
#ifndef NWAD_PRINT
#define NWAD_PRINT
c
c     Compile source again for Maxima
c
#include "nwxc_x_b97.F"
#endif
#ifndef SECOND_DERIV
#define SECOND_DERIV
c
c     Compile source again for the 2nd derivative case
c
#include "nwxc_x_b97.F"
#endif
#ifndef THIRD_DERIV
#define THIRD_DERIV
c
c     Compile source again for the 3rd derivative case
c
#include "nwxc_x_b97.F"
#endif
#undef NWAD_PRINT
C>
C> @}
