#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if !defined(NWAD_PRINT)
C> \ingroup nwxc
C> @{
C>
C> \file nwxc_x_opt.F
C> The OPT exchange functional
C>
C> @}
#endif
#endif
C>
C> \ingroup nwxc_priv
C> @{
C>
C> \brief The OPT exchange functional
C>
C> The OPTimized exchange functional [1,2].
C>
C> [1] N.C. Handy, A.J. Cohen, "Left-right correlation energy",
C>     Mol. Phys. <b>99</b>, 403-412 (2001), DOI:
C>     <A HREF="https://doi.org/10.1080/00268970010018431">
C>     10.1080/00268970010018431</A>.
C>
C> [2] N.C. Handy, A.J. Cohen, "Dynamic correlation",
C>     Mol. Phys. <b>99</b>, 607-615 (2001), DOI:
C>     <A HREF="https://doi.org/10.1080/00268970010023435">
C>     10.1080/00268970010023435</A>.
c
C$Id$
c
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
      subroutine nwxc_x_opt_p(tol_rho, ipol, nq, wght, rho, rgamma,
     &                        func)
#else
      subroutine nwxc_x_opt(tol_rho, ipol, nq, wght, rho, rgamma, func)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      subroutine nwxc_x_opt_d2(tol_rho, ipol, nq, wght, rho, rgamma,
     &                         func)
#else
      subroutine nwxc_x_opt_d3(tol_rho, ipol, nq, wght, rho, rgamma,
     &                         func)
#endif
c
#include "nwad.fh"
c
      implicit none
c      
#include "nwxc_param.fh"
c
      double precision tol_rho !< [Input] The lower limit on the density
      integer nq               !< [Input] The number of points
      integer ipol             !< [Input] The number of spin channels
      double precision wght    !< [Input] The weight of the functional
c
c     Charge Density
c
      type(nwad_dble)::rho(nq,*) !< [Input] The density
c
c     Charge Density Gradient
c
      type(nwad_dble)::rgamma(nq,*) !< [Input] The norm of the density gradients
c
c     Sampling Matrices for the XC Potential & Energy
c
      type(nwad_dble)::func(nq)   !< [Output] The value of the functional
c     double precision Amat(nq,*) !< [Output] The derivative wrt rho
c     double precision Cmat(nq,*) !< [Output] The derivative wrt rgamma
c
c
c References:
c
c    Becke,  (1986)
c    Handy NC, Cohen AJ, Mol Phys 99 (5): 403-412 MAR 2001
c    idem, Mol Phys 99 (7); 607-615 2001
c
c***************************************************************************
c
      integer n
      type(nwad_dble)::rho13, rho43, gamma, x, g,  dg,
     &     t, hrho
      double precision gamma86,gx
      Parameter (gamma86=0.006d0)
c     type(nwad_dble)::ux,uxp
      double precision p
      parameter (p = 2.0d0)
c     ux(x,gx)=gx*x*x/(1d0+gx*x*x)
c     uxp(x,gx)=gx*x*2d0/(1d0+gx*x*x)**2
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
            hrho  = 0.5d0*rho(n,R_T)
            rho13 = hrho**(1.d0/3.d0)
            rho43 = rho13*hrho
c           gamma = delrho(n,1,1)*delrho(n,1,1) +
c    &              delrho(n,2,1)*delrho(n,2,1) +
c    &              delrho(n,3,1)*delrho(n,3,1)
            gamma = rgamma(n,G_TT)
c           if (sqrt(gamma).gt.tol_rho)then
c              gamma = 0.25d0 * gamma
c              x = sqrt(gamma) / rho43
            x = gamma/(4.0d0*rho43*rho43)
c           else
c              x = 0d0
c           endif
c
            g = -ux(x,gamma86)**p
c           dg = -p*ux(x,gamma86)**(p-1)*uxp(x,gamma86)
c
c
            func(n) = func(n) + 2.d0*rho43*g*wght
c           Amat(n,D1_RA) = Amat(n,D1_RA)+(4d0/3d0)*rho13*(g-x*dg)*wght
c
c           if (x.gt.tol_rho) then
c              t = 0.5d0 * dg / sqrt(gamma) * wght
c              Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + t
c           endif
c
c
 10      continue
c
      else
c
c        ======> SPIN-UNRESTRICTED <======
c
         do 20 n = 1, nq
            if (rho(n,R_A)+rho(n,R_B).lt.tol_rho) goto 20
            if (rho(n,R_A).lt.tol_rho) goto 25
c
c           Spin alpha:
c
            rho13 = rho(n,R_A)**(1.d0/3.d0)
            rho43 = rho13*rho(n,R_A)
c           gamma = delrho(n,1,1)*delrho(n,1,1) +
c    &              delrho(n,2,1)*delrho(n,2,1) +
c    &              delrho(n,3,1)*delrho(n,3,1)
            gamma = rgamma(n,G_AA)
c           if (sqrt(gamma).gt.tol_rho)then
c              x = sqrt(gamma) / rho43
            x = gamma/(rho43*rho43)
c           else
c              x = 0d0
c           endif
c
            g = -ux(x,gamma86)**p
c           dg = -p*ux(x,gamma86)**(p-1)*uxp(x,gamma86)
c
c
            func(n) = func(n) + rho43*g*wght
c           Amat(n,D1_RA) = Amat(n,D1_RA)+(4d0/3d0)*rho13*(g-x*dg)*wght
c     
c           if (x.gt.tol_rho) then
c              t = dg / sqrt(gamma) * wght
c              Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + t * 0.5d0
c           endif
c
c
 25         continue
c
c           Spin beta:
c
            if (rho(n,R_B).lt.tol_rho) goto 20
c
            rho13 = rho(n,R_B)**(1.d0/3.d0)
            rho43 = rho13*rho(n,R_B)
c           gamma = delrho(n,1,2)*delrho(n,1,2) +
c    &              delrho(n,2,2)*delrho(n,2,2) +
c    &              delrho(n,3,2)*delrho(n,3,2)
            gamma = rgamma(n,G_BB)
c           if (sqrt(gamma).gt.tol_rho)then
c              x = sqrt(gamma) / rho43
            x = gamma/(rho43*rho43)
c           else
c              x = 0d0
c           endif
            g = -ux(x,gamma86)**p
c           dg = -p*ux(x,gamma86)**(p-1)*uxp(x,gamma86)
c
c
            func(n) = func(n) +rho43*g*wght
c           Amat(n,D1_RB) = Amat(n,D1_RB)+(4d0/3d0)*rho13*(g-x*dg)*wght
c     
c           if (x.gt.tol_rho) then
c              t = dg / sqrt(gamma) * wght
c              Cmat(n,D1_GBB) = Cmat(n,D1_GBB) + t * 0.5d0
c           endif
c
c
 20      continue
c
      endif
c
      return
c
      contains
c
c     The combination of statement functions and derived types with
c     overloaded operators is not properly supported in the Fortran 
c     standard (apparently). Therefore the statement functions from
c     the original subroutine had to be transformed into contained 
c     functions.
c
c     WARNING: This code is EXTREMELY EVIL! Although there appears to be
c     only one function there are actually three with the same name,
c     each one returning results of a different data type. The trick is
c     that depending on the data type the the subroutine that contains
c     these functions changes its name and hence these different
c     functions of the same name do not lead to conflicts. The 
c     alternative would have been to add a forest of conditional 
c     compilation constructs (#ifdef's) to change the function names
c     in the declarations and all places where they are used. That 
c     would have been extremely ugly, so we are between a rock and a
c     hard place on this one.
c
      function ux(x,gx) result(s)
#include "nwad.fh"
        implicit none
        type(nwad_dble),  intent(in) :: x
        double precision, intent(in) :: gx
        type(nwad_dble)              :: s
        s = gx*x/(1d0+gx*x)
      end function
c
      end
#if !defined(NWAD_PRINT)
#define NWAD_PRINT
c
c     Compile again for Maxima
c
#include "nwxc_x_opt.F"
#endif
#if !defined(SECOND_DERIV)
#define SECOND_DERIV
c
c     Compile again for 2nd derivatives
c
#include "nwxc_x_opt.F"
#endif
#if !defined(THIRD_DERIV)
#define THIRD_DERIV
c
c     Compile again for 3rd derivatives
c
#include "nwxc_x_opt.F"
#endif
#undef NWAD_PRINT
C> @}
