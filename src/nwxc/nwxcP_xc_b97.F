#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if !defined(NWAD_PRINT)
!> 
!> \ingroup nwxc
!> @{
!> 
!> \file nwxcP_xc_b97.c
!> The code for the B97 family of functionals
!> 
!> \brief Evaluate functionals from the B97 family of functionals
!> 
!> These functions evaluate the exchange and correlation terms of the B97 
!> family of functionals [1]. The code was provided by Tobias Risthaus [2]
!> and based on [3].
!> The code was initially modified to integrate it into the NWDFT
!> framework. Later it was modified to integrate it into the NWXC 
!> framework. These latter changes included separating out the correlation
!> and exchange terms as well as the usual changes to the top-level 
!> interface. 
!>
!> With the introduction of the NWAD framework for automatic differentiation
!> the routines were translated into Fortran.
!> 
!> Obviously we maintained the structure that allows 3 and 5 term functionals.
!> 
!> ### References ###
!> 
!> [1] A.D. Becke, "Density-functional thermochemistry. V. Systematic
!>     optimization of exchange-correlation functionals", J. Chem. Phys.
!>     107 (1997) 8554-8560, DOI:
!>     <a href="https://doi.org/10.1063/1.475007">
!>     10.1063/1.475007</a>.
!> 
!> [2] Tobias Risthaus (tobias.risthaus@thch.uni-bonn.de),
!>     Mulliken Center for Theoretical Chemistry,
!>     Institut für Physikalische und Theoretische Chemie,
!>     Universität Bonn, Beringstr. 4, D-53115 Bonn, Germany.
!> 
!> [3] S. Grimme, "Semiempirical GGA-type density functional constructed
!>     with a long-range dispersion correction", J. Comput. Chem. 27 (2006)
!>     1787-1799, DOI:
!>     <a href="https://doi.org/10.1002/jcc.20495">
!>     10.1002/jcc.20495</a>.
!> @}
!>
#endif
#endif
!>
!> \ingroup nwxc_priv
!> \@{
#include "nwxcP_xc_b97.fh"

#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#ifdef NWAD_PRINT
      subroutine nwxcP_xc_gss3tar_p(r, g, alpha, a, f)
#else
      subroutine nwxcP_xc_gss3tar(r, g, alpha, a, f)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      subroutine nwxcP_xc_gss3tar_d2(r, g, alpha, a, f)
#else
      subroutine nwxcP_xc_gss3tar_d3(r, g, alpha, a, f)
#endif
#include "nwad.fh"
      implicit none
      type(nwad_dble):: r, g, f
      double precision alpha
#ifdef NWAD_PRINT
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      type(nwad_dble)::a(3)
#else
      double precision a(3)
#endif
#else
      double precision a(3)
#endif
!     B97 like 
!     S. Grimme, J. Comput. Chem., 2006, 27, 1787-1799 
!     doi: 10.1002/jcc.20495
!     for now, this is only valid for 3 terms. 

      type(nwad_dble):: r83,r113,r143,s2,s2p1,u

      r83=r**(8.0d0/3.0d0)
!     r113=r83*r
!     r143=r113*r
      s2=g/r83
      s2p1=1.0d0+alpha*s2

      u=alpha*s2/s2p1

      f=a(1) + a(2)*u + a(3)*u**2

      end

#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#ifdef NWAD_PRINT
      subroutine nwxcP_xc_gss5tar_p(r, g, alpha, a, f)
#else
      subroutine nwxcP_xc_gss5tar(r, g, alpha, a, f)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      subroutine nwxcP_xc_gss5tar_d2(r, g, alpha, a, f)
#else
      subroutine nwxcP_xc_gss5tar_d3(r, g, alpha, a, f)
#endif
#include "nwad.fh"
      implicit none
      type(nwad_dble):: r, g, f
      double precision alpha
#ifdef NWAD_PRINT
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      type(nwad_dble)::a(5)
#else
      double precision a(5)
#endif
#else
      double precision a(5)
#endif

!     B97 like 
!     S. Grimme, J. Comput. Chem., 2006, 27, 1787-1799 
!     doi: 10.1002/jcc.20495
!     for now, this is only valid for 5 terms. 

      type(nwad_dble):: r83,r113,r143,s2,s2p1,u

      r83=r**(8.0d0/3.0d0)
!     r113=r83*r
!     r143=r113*r
      s2=g/r83
      s2p1=1.0d0+alpha*s2

      u=alpha*s2/s2p1

      f=a(1) + a(2)*u + a(3)*(u**2) + a(4)*(u**3) + a(5)*(u**4)

      end

#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#ifdef NWAD_PRINT
      subroutine nwxcP_xc_gab3tar_p(ra, rb, ga, gb, alpha, a, f)
#else
      subroutine nwxcP_xc_gab3tar(ra, rb, ga, gb, alpha, a, f)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      subroutine nwxcP_xc_gab3tar_d2(ra, rb, ga, gb, alpha, a, f)
#else
      subroutine nwxcP_xc_gab3tar_d3(ra, rb, ga, gb, alpha, a, f)
#endif
#include "nwad.fh"
      implicit none
      double precision alpha
#ifdef NWAD_PRINT
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      type(nwad_dble)::a(3)
#else
      double precision a(3)
#endif
#else
      double precision a(3)
#endif
      type(nwad_dble):: ra, rb, ga, gb, f

!     B97 like 
!     S. Grimme, J. Comput. Chem., 2006, 27, 1787-1799 
!     doi: 10.1002/jcc.20495
!     for now, this is only valid for 3 terms. 

      type(nwad_dble):: ra83,rb83,ra113,rb113,ra143,rb143,sa2,sb2,
     +                  sv2,sv2p1,u

      ra83=ra**(8.0d0/3.0d0)
!     ra113=ra83*ra
!     ra143=ra113*ra
      rb83=rb**(8.0d0/3.0d0)
!     rb113=rb83*rb
!     rb143=rb113*rb
      sa2=ga/ra83
      sb2=gb/rb83
      sv2=0.5d0*(sa2 + sb2)                  ! called s_av in some papers 
      sv2p1=1.0d0+alpha*sv2

      u=alpha*sv2/sv2p1
      f=a(1) + a(2)*u + a(3)*u**2

      end

#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#ifdef NWAD_PRINT
      subroutine nwxcP_xc_gab5tar_p(ra, rb, ga, gb, alpha, a, f)
#else
      subroutine nwxcP_xc_gab5tar(ra, rb, ga, gb, alpha, a, f)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      subroutine nwxcP_xc_gab5tar_d2(ra, rb, ga, gb, alpha, a, f)
#else
      subroutine nwxcP_xc_gab5tar_d3(ra, rb, ga, gb, alpha, a, f)
#endif
#include "nwad.fh"
      implicit none
      double precision alpha
#ifdef NWAD_PRINT
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      type(nwad_dble)::a(5)
#else
      double precision a(5)
#endif
#else
      double precision a(5)
#endif
      type(nwad_dble):: ra, rb, ga, gb, f

!     B97 like 
!     S. Grimme, J. Comput. Chem., 2006, 27, 1787-1799 
!     doi: 10.1002/jcc.20495
!     for now, this is only valid for 5 terms. 

      type(nwad_dble):: ra83,rb83,ra113,rb113,ra143,rb143,sa2,sb2,
     +                  sv2,sv2p1,u

      ra83=ra**(8.0d0/3.0d0)
!     ra113=ra83*ra
!     ra143=ra113*ra
      rb83=rb**(8.0d0/3.0d0)
!     rb113=rb83*rb
!     rb143=rb113*rb
      sa2=ga/ra83
      sb2=gb/rb83
      sv2=0.5d0*(sa2 + sb2)                  ! called s_av in some papers 
      sv2p1=1.0d0+alpha*sv2

      u=alpha*sv2/sv2p1
      f=a(1) + a(2)*u + a(3)*u**2 + a(4)*u**3 + a(5)*u**4

      end


!/ =============================================================================
!/
!/                      B97 exchange
!/
!/ =============================================================================


#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#ifdef NWAD_PRINT
      subroutine nwxcp_x_b97_p(rho_a,rho_b,ScalGGAX,tol_rho,FX,sol)
#else
      subroutine nwxcp_x_b97(rho_a,rho_b,ScalGGAX,tol_rho,FX,sol)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      subroutine nwxcp_x_b97_d2(rho_a,rho_b,ScalGGAX,tol_rho,FX,sol)
#else
      subroutine nwxcp_x_b97_d3(rho_a,rho_b,ScalGGAX,tol_rho,FX,sol)
#endif
#include "nwad.fh"
      implicit none
      type(nwad_dble):: rho_a(2), rho_b(2), FX(0:0)
      double precision ScalGGAX, tol_rho
#ifdef NWAD_PRINT
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      type(nwad_dble)::sol(6), xss(5)
#else
      double precision sol(6), xss(5)
#endif
#else
      double precision sol(6), xss(5)
#endif
      integer max_pow_u, i
      double precision gamma_xss
      parameter(gamma_xss = 0.004d0)

      type(nwad_dble):: ra, rb, gaa, gbb
      type(nwad_dble):: ra43, rb43
      ! enhancement factors and derivatives thereof 
      type(nwad_dble):: exaa, exbb
      double precision CX,CX43,CX49

      ra   = rho_a(1)
      rb   = rho_b(1)
      gaa  = rho_a(2) ! gamma_aa = G_AA
      gbb  = rho_b(2) ! gamma_bb = G_BB
#ifdef NWAD_PRINT
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      max_pow_u = value(sol(1))
#else
      max_pow_u = sol(1)
#endif
#else
      max_pow_u = sol(1)
#endif


      do i = 1, max_pow_u+1
        xss(i)=sol(i+1)
      enddo

      !  CX    = 0.930525736349
      CX = 0.9305257363490993d0
      ! CX43  = (4.0d0/3.0d0)*CX
      ! CX49  = (4.0d0/9.0d0)*CX

      if(ra > tol_rho ) then
        if(max_pow_u == 2 ) then
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#ifdef NWAD_PRINT
          call nwxcP_xc_gss3tar_p(ra,gaa,gamma_xss,xss,exaa)
#else
          call nwxcP_xc_gss3tar(ra,gaa,gamma_xss,xss,exaa)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
          call nwxcP_xc_gss3tar_d2(ra,gaa,gamma_xss,xss,exaa)
#else
          call nwxcP_xc_gss3tar_d3(ra,gaa,gamma_xss,xss,exaa)
#endif
        else
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#ifdef NWAD_PRINT
          call nwxcP_xc_gss5tar_p(ra,gaa,gamma_xss,xss,exaa)
#else
          call nwxcP_xc_gss5tar(ra,gaa,gamma_xss,xss,exaa)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
          call nwxcP_xc_gss5tar_d2(ra,gaa,gamma_xss,xss,exaa)
#else
          call nwxcP_xc_gss5tar_d3(ra,gaa,gamma_xss,xss,exaa)
#endif
        endif
        ! Slater-X and derivatives 
        ra43 = -CX  *ra**(4.0d0/3.0d0)
        ! ra13 = -CX43*ra**(1.0d0/3.0d0)
        ! raM23= -CX49*ra**(-2.0d0/3.0d0)
      else
        ra43 = 0.0d0
        exaa = 0.0d0
        ! ra13 = 0.0
        ! raM23= 0.0
      endif

      if(rb > tol_rho ) then
        if(max_pow_u == 2 )then
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#ifdef NWAD_PRINT
          call nwxcP_xc_gss3tar_p(rb,gbb,gamma_xss,xss,exbb)
#else
          call nwxcP_xc_gss3tar(rb,gbb,gamma_xss,xss,exbb)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
          call nwxcP_xc_gss3tar_d2(rb,gbb,gamma_xss,xss,exbb)
#else
          call nwxcP_xc_gss3tar_d3(rb,gbb,gamma_xss,xss,exbb)
#endif
        else
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#ifdef NWAD_PRINT
          call nwxcP_xc_gss5tar_p(rb,gbb,gamma_xss,xss,exbb)
#else
          call nwxcP_xc_gss5tar(rb,gbb,gamma_xss,xss,exbb)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
          call nwxcP_xc_gss5tar_d2(rb,gbb,gamma_xss,xss,exbb)
#else
          call nwxcP_xc_gss5tar_d3(rb,gbb,gamma_xss,xss,exbb)
#endif
        endif
        ! Slater-X, more or less
        rb43 = -CX  *rb**(4.0d0/3.0d0)
        ! rb13 = -CX43*rb**(1.0d0/3.0d0)
        ! rbM23= -CX49*rb**(-2.0d0/3.0)
      else
        rb43 = 0.0d0
        exbb = 0.0d0
        ! rb13 = 0.0
        ! rbM23 = 0.0
      endif

      ! ScalGGAX is the ACM parameter for the hybrid
      ! recall that B97 does not obey the UEG limit 

      FX( _FXC_E ) = ScalGGAX*(ra43*exaa + rb43*exbb)

      end

! =============================================================================
!
!                      B97 correlation 
!
! =============================================================================

#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#ifdef NWAD_PRINT
      subroutine nwxcp_c_b97_p(rho_a,rho_b,ScalGGAC,tol_rho,FC,sol)
#else
      subroutine nwxcp_c_b97(rho_a,rho_b,ScalGGAC,tol_rho,FC,sol)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      subroutine nwxcp_c_b97_d2(rho_a,rho_b,ScalGGAC,tol_rho,FC,sol)
#else
      subroutine nwxcp_c_b97_d3(rho_a,rho_b,ScalGGAC,tol_rho,FC,sol)
#endif
#include "nwad.fh"
      implicit none
#include "intf_nwxcP_c_pwlda.fh"
#include "intf_nwxcP_xc_gab3tar.fh"
#include "intf_nwxcP_xc_gss3tar.fh"
#include "intf_nwxcP_xc_gab5tar.fh"
#include "intf_nwxcP_xc_gss5tar.fh"
      type(nwad_dble):: rho_a(2), rho_b(2), FC(0:0)
      double precision ScalGGAC, tol_rho
#ifdef NWAD_PRINT
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      type(nwad_dble)::sol(11), css(5), cab(5)
#else
      double precision sol(11), css(5), cab(5)
#endif
#else
      double precision sol(11), css(5), cab(5)
#endif

      type(nwad_dble):: ra,rb,gaa,gbb, ec_a0, ec_b0, ec_ab
      type(nwad_dble):: gcaa, gcbb, gcab
      integer max_pow_u, i
      type(nwad_dble)::  FCLDA(0:0), t

      double precision gamma_css, gamma_cab
      parameter(gamma_css = 0.2d0, gamma_cab = 0.006d0)

      ra   = rho_a(1)
      rb   = rho_b(1)
      gaa  = rho_a(2)
      gbb  = rho_b(2)
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#ifdef NWAD_PRINT
      max_pow_u = value(sol(1))
#else
      max_pow_u = sol(1)
#endif
#else
      max_pow_u = sol(1)
#endif


      do i = 1, max_pow_u+1
        css(i)=sol((i-1)*2+2)
        cab(i)=sol((i-1)*2+3)
      enddo

      ! normal
      if((ra > tol_rho) .or. (rb > tol_rho)) then
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#ifdef NWAD_PRINT
        call nwxcP_c_pwlda_p(ra, rb, FCLDA)
#else
        call nwxcP_c_pwlda(ra, rb, FCLDA)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
        call nwxcP_c_pwlda_d2(ra, rb, FCLDA)
#else
        call nwxcP_c_pwlda_d3(ra, rb, FCLDA)
#endif
        ! LDA derivatives
        ec_ab = FCLDA(_FXC_E) 
      else
        ec_ab = 0.0d0
      endif

      ! alpha density only
      if(ra > tol_rho) then
        t = 0.0d0
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#ifdef NWAD_PRINT
        call nwxcP_c_pwlda_p(ra, t, FCLDA)
#else
        call nwxcP_c_pwlda(ra, t, FCLDA)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
        call nwxcP_c_pwlda_d2(ra, t, FCLDA)
#else
        call nwxcP_c_pwlda_d3(ra, t, FCLDA)
#endif
        ec_a0 = FCLDA(_FXC_E)
      else 
        ec_a0 = 0.0d0
      endif

      ! beta density only
      if(rb > tol_rho) then
        t = 0.0d0
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#ifdef NWAD_PRINT
        call nwxcP_c_pwlda_p(t, rb, FCLDA)
#else
        call nwxcP_c_pwlda(t, rb, FCLDA)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
        call nwxcP_c_pwlda_d2(t, rb, FCLDA)
#else
        call nwxcP_c_pwlda_d3(t, rb, FCLDA)
#endif
        ec_b0 = FCLDA(_FXC_E)
      else
        ec_b0 = 0.0d0
      endif

      ec_ab = ec_ab - ec_a0 - ec_b0

      if(ra > tol_rho) then
        if(max_pow_u == 2 ) then
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#ifdef NWAD_PRINT
          call nwxcP_xc_gss3tar_p(ra,gaa,gamma_css,css,gcaa)
#else
          call nwxcP_xc_gss3tar(ra,gaa,gamma_css,css,gcaa)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
          call nwxcP_xc_gss3tar_d2(ra,gaa,gamma_css,css,gcaa)
#else
          call nwxcP_xc_gss3tar_d3(ra,gaa,gamma_css,css,gcaa)
#endif
        else
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#ifdef NWAD_PRINT
          call nwxcP_xc_gss5tar_p(ra,gaa,gamma_css,css,gcaa)
#else
          call nwxcP_xc_gss5tar(ra,gaa,gamma_css,css,gcaa)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
          call nwxcP_xc_gss5tar_d2(ra,gaa,gamma_css,css,gcaa)
#else
          call nwxcP_xc_gss5tar_d3(ra,gaa,gamma_css,css,gcaa)
#endif
        endif
      else
        ! gcaa = css(1) + css(2) + css(3) !  u -> 1 if s -> infinity
        gcaa = 0.0d0 !  Elda(r)*gcaa -> 0 if r -> 0 hence OK to zero gcaa
      endif

      if(rb > tol_rho) then
        if(max_pow_u == 2 ) then
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#ifdef NWAD_PRINT
          call nwxcP_xc_gss3tar_p(rb,gbb,gamma_css,css,gcbb)
#else
          call nwxcP_xc_gss3tar(rb,gbb,gamma_css,css,gcbb)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
          call nwxcP_xc_gss3tar_d2(rb,gbb,gamma_css,css,gcbb)
#else
          call nwxcP_xc_gss3tar_d3(rb,gbb,gamma_css,css,gcbb)
#endif
        else
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#ifdef NWAD_PRINT
          call nwxcP_xc_gss5tar_p(rb,gbb,gamma_css,css,gcbb)
#else
          call nwxcP_xc_gss5tar(rb,gbb,gamma_css,css,gcbb)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
          call nwxcP_xc_gss5tar_d2(rb,gbb,gamma_css,css,gcbb)
#else
          call nwxcP_xc_gss5tar_d3(rb,gbb,gamma_css,css,gcbb)
#endif
        endif
      else
        ! gcbb = css(1) + css(2) + css(3) !  u -> 1 if s -> infinity
        gcbb = 0.0d0 !  Elda(r)*gcbb -> 0 if r -> 0 hence OK to zero gcbb
      endif

      if(ra > tol_rho .and. rb > tol_rho) then
        if(max_pow_u == 2 ) then
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#ifdef NWAD_PRINT
          call nwxcP_xc_gab3tar_p(ra,rb,gaa,gbb,gamma_cab,cab,gcab)
#else
          call nwxcP_xc_gab3tar(ra,rb,gaa,gbb,gamma_cab,cab,gcab)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
          call nwxcP_xc_gab3tar_d2(ra,rb,gaa,gbb,gamma_cab,cab,gcab)
#else
          call nwxcP_xc_gab3tar_d3(ra,rb,gaa,gbb,gamma_cab,cab,gcab)
#endif
        else
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#ifdef NWAD_PRINT
          call nwxcP_xc_gab5tar_p(ra,rb,gaa,gbb,gamma_cab,cab,gcab)
#else
          call nwxcP_xc_gab5tar(ra,rb,gaa,gbb,gamma_cab,cab,gcab)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
          call nwxcP_xc_gab5tar_d2(ra,rb,gaa,gbb,gamma_cab,cab,gcab)
#else
          call nwxcP_xc_gab5tar_d3(ra,rb,gaa,gbb,gamma_cab,cab,gcab)
#endif
        endif
      else
        ! This term is not necessary for exchange as there is no cross
        ! term between the spin-channels in that case.
        ! gcab = cab(1) + cab(2) + cab(3) !  u == 1 if s->inf. 
        gcab = 0.0d0 !  Elda(r)*gcab -> 0 if r -> 0 hence OK to zero gcab
      endif

      ! ScalGGAC is the ACM D parameter as used in double hybrids 

      FC( _FXC_E ) = ScalGGAC*(ec_a0*gcaa + ec_b0*gcbb + ec_ab*gcab)

      end
#ifndef NWAD_PRINT
#define NWAD_PRINT
c
c     Compile source again for Maxima
c
#include "nwxcP_xc_b97.F"
#endif
#ifndef SECOND_DERIV
#define SECOND_DERIV
c
c     Compile source again for the 2nd derivative case
c
#include "nwxcP_xc_b97.F"
#endif
#ifndef THIRD_DERIV
#define THIRD_DERIV
c
c     Compile source again for the 3rd derivative case
c
#include "nwxcP_xc_b97.F"
#endif
#undef NWAD_PRINT
!>
!> @}
!> 
!> $Id$
