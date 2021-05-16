#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
C> \ingroup nwxc
C> @{
C>
C> \file nwxc_c_b97.F
C> The B97 family of correlation functionals
C>
C> @}
#include "nwxcP_xc_b97.fh"
#endif
C>
C> \ingroup nwxc_priv
C> @{
C>
C> \brief Evaluate the B97 family of correlation functionals
C>
C> This code evaluates correlation functionals from the 
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
C>     with a long-range dispersion correction", J. Comput. Chem. 27
C>     (2006) 1787-1799, DOI:
C>     <a href="https://doi.org/10.1002/jcc.20495">
C>     10.1002/jcc.20495</a>.
C>
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#ifdef NWAD_PRINT
      Subroutine nwxc_c_b97_p(param,tol_rho,ipol,nq,wght,rho,rgamma,
     &                      func)
#else
      Subroutine nwxc_c_b97(param,tol_rho,ipol,nq,wght,rho,rgamma,
     &                      func)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      Subroutine nwxc_c_b97_d2(param,tol_rho,ipol,nq,wght,rho,rgamma,
     &                         func)
#else
      Subroutine nwxc_c_b97_d3(param,tol_rho,ipol,nq,wght,rho,rgamma,
     &                         func)
#endif
c
c     $Id$
c
#include "nwad.fh"
      implicit none
#include "intf_nwxcP_c_b97.fh"
c
#include "nwxc_param.fh"
c
c     Input and other parameters
c
#ifdef NWAD_PRINT
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      type(nwad_dble)::param(*)
#else
      double precision param(*)
#endif
#else
      double precision param(*)!< [Input] Parameters of functional as
                               !< defined in [1]:
                               !< - param(1): \f$m\f$ of Eqs.(20).
                               !< - param(2): \f$C_{C\sigma\sigma,0}\f$
                               !< - param(3): \f$C_{C\alpha\beta,0}\f$
                               !< - param(4): \f$C_{C\sigma\sigma,1}\f$
                               !< - param(5): \f$C_{C\alpha\beta,1}\f$
                               !< - param(6): \f$C_{C\sigma\sigma,2}\f$
                               !< - param(7): \f$C_{C\alpha\beta,2}\f$
                               !< - param(8): \f$C_{C\sigma\sigma,3}\f$
                               !< - param(9): \f$C_{C\alpha\beta,3}\f$
                               !< - param(10): \f$C_{C\sigma\sigma,4}\f$
                               !< - param(11): \f$C_{C\alpha\beta,4}\f$
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
      type(nwad_dble)::FC(0:_FXC_NUMDERI)
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
#ifdef NWAD_PRINT
            call nwxcp_c_b97_p(rho_a,rho_b,1.0d0,tol_rho,FC,param)
#else
            call nwxcp_c_b97(rho_a,rho_b,1.0d0,tol_rho,FC,param)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
            call nwxcp_c_b97_d2(rho_a,rho_b,1.0d0,tol_rho,FC,param)
#else
            call nwxcp_c_b97_d3(rho_a,rho_b,1.0d0,tol_rho,FC,param)
#endif
            func(i)        = func(i)        + FC(_FXC_E)*wght
!           Amat(i,D1_RA)  = Amat(i,D1_RA)  + FC(_FXC_RA)*wght
!           Cmat(i,D1_GAA) = Cmat(i,D1_GAA) + FC(_FXC_GAA)*wght
          endif
        enddo
      else
        do i = 1, nq
          rho_a(0) = 0.0d0
          rho_a(1) = 0.0d0
          rho_b(0) = 0.0d0
          rho_b(1) = 0.0d0
          if (rho(i,R_A).gt.tol_rho) then
            rho_a(0) = rho(i,R_A)
            rho_a(1) = rgamma(i,G_AA)
          endif
          if (rho(i,R_B).gt.tol_rho) then
            rho_b(0) = rho(i,R_B)
            rho_b(1) = rgamma(i,G_BB)
          endif
          if (rho_a(0).gt.tol_rho.or.rho_b(0).gt.tol_rho) then
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#ifdef NWAD_PRINT
            call nwxcp_c_b97_p(rho_a,rho_b,1.0d0,tol_rho,FC,param)
#else
            call nwxcp_c_b97(rho_a,rho_b,1.0d0,tol_rho,FC,param)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
            call nwxcp_c_b97_d2(rho_a,rho_b,1.0d0,tol_rho,FC,param)
#else
            call nwxcp_c_b97_d3(rho_a,rho_b,1.0d0,tol_rho,FC,param)
#endif
            func(i)        = func(i)        + FC(_FXC_E)*wght
!           Amat(i,D1_RA)  = Amat(i,D1_RA)  + FC(_FXC_RA)*wght
!           Amat(i,D1_RB)  = Amat(i,D1_RB)  + FC(_FXC_RB)*wght
!           Cmat(i,D1_GAA) = Cmat(i,D1_GAA) + FC(_FXC_GAA)*wght
!           Cmat(i,D1_GBB) = Cmat(i,D1_GBB) + FC(_FXC_GBB)*wght
          endif
        enddo
      endif
c
      end
c
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#ifdef NWAD_PRINT
      subroutine nwxcp_c_pwlda_p(ra,rb,FCLDA)
#else
      subroutine nwxcp_c_pwlda(ra,rb,FCLDA)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      subroutine nwxcp_c_pwlda_d2(ra,rb,FCLDA)
#else
      subroutine nwxcp_c_pwlda_d3(ra,rb,FCLDA)
#endif
#include "nwad.fh"
      implicit none
#include "intf_nwxc_c_pw91lda.fh"
#include "nwxc_param.fh"
c
      type(nwad_dble)::ra
      type(nwad_dble)::rb
c     type(nwad_dble)::FCLDA(0:_FCLDA_ELEMENTS)
      type(nwad_dble)::FCLDA(0:_FXC_E)
c
      type(nwad_dble)::ec(1)
      type(nwad_dble)::rho(2)
c     double precision Amat(2)
c     double precision Amat2(3)
c
      rho(R_A)        = ra
      rho(R_B)        = rb
      ec(1)           = 0.0d0
!     Amat(D1_RA)     = 0.0d0
!     Amat(D1_RB)     = 0.0d0
!     Amat2(D2_RA_RA) = 0.0d0
!     Amat2(D2_RA_RB) = 0.0d0
!     Amat2(D2_RB_RB) = 0.0d0
c
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#ifdef NWAD_PRINT
      call nwxc_c_pw91lda_p(1.0d-20,2,1,1.0d0,rho,ec)
#else
      call nwxc_c_pw91lda(1.0d-20,2,1,1.0d0,rho,ec)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      call nwxc_c_pw91lda_d2(1.0d-20,2,1,1.0d0,rho,ec)
#else
      call nwxc_c_pw91lda_d3(1.0d-20,2,1,1.0d0,rho,ec)
#endif
c
      FCLDA(_FXC_E)    = ec(1)
!     FCLDA(_FXC_RA)   = Amat(D1_RA)
!     FCLDA(_FXC_RB)   = Amat(D1_RB)
!     FCLDA(_FXC_RARA) = Amat2(D2_RA_RA)
!     FCLDA(_FXC_RARB) = Amat2(D2_RA_RB)
!     FCLDA(_FXC_RBRB) = Amat2(D2_RB_RB)
c
      end
#ifndef NWAD_PRINT
#define NWAD_PRINT
c
c     Compile source again for Maxima
c
#include "nwxc_c_b97.F"
#endif
#ifndef SECOND_DERIV
#define SECOND_DERIV
c
c     Compile source again for the 2nd derivative case
c
#include "nwxc_c_b97.F"
#endif
#ifndef THIRD_DERIV
#define THIRD_DERIV
c
c     Compile source again for the 3rd derivative case
c
#include "nwxc_c_b97.F"
#endif
#undef NWAD_PRINT
C>
C> @}
