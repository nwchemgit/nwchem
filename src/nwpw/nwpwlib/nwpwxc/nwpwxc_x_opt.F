C> \ingroup nwpwxc
C> @{
C>
C> \file nwpwxc_x_opt.F
C> The OPT exchange functional
C>
C> @}
C>
C> \ingroup nwpwxc_priv
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
      subroutine nwpwxc_x_opt(tol_rho,ipol,nq,wght,rho,rgamma,func,
     &                      Amat, Cmat)
      implicit none
c      
#include "nwpwxc_param.fh"
c
      double precision tol_rho !< [Input] The lower limit on the density
      integer nq               !< [Input] The number of points
      integer ipol             !< [Input] The number of spin channels
      double precision wght    !< [Input] The weight of the functional
c
c     Charge Density
c
      double precision rho(nq,*) !< [Input] The density
c
c     Charge Density Gradient
c
      double precision rgamma(nq,*) !< [Input] The norm of the density gradients
c
c     Sampling Matrices for the XC Potential & Energy
c
      double precision func(nq)   !< [Output] The value of the functional
      double precision Amat(nq,*) !< [Output] The derivative wrt rho
      double precision Cmat(nq,*) !< [Output] The derivative wrt rgamma
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
      double precision rho13, rho43, gamma, x, g,  dg,
     &     t, hrho
      double precision gamma86
      Parameter (gamma86=0.006d0)
      double precision ux,uxp,gx
      integer p
      parameter (p = 2)
      ux(x,gx)=gx*x*x/(1d0+gx*x*x)
      uxp(x,gx)=gx*x*2d0/(1d0+gx*x*x)**2
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
            if (dsqrt(gamma).gt.tol_rho)then
               gamma = 0.25d0 * gamma
               x = sqrt(gamma) / rho43
            else
               x = 0d0
            endif
c
            g = -ux(x,gamma86)**p
            dg = -p*ux(x,gamma86)**(p-1)*uxp(x,gamma86)
c
c
            func(n) = func(n) + 2.d0*rho43*g*wght
            Amat(n,D1_RA) = Amat(n,D1_RA)+(4d0/3d0)*rho13*(g-x*dg)*wght
c
            if (x.gt.tol_rho) then
               t = 0.5d0 * dg / sqrt(gamma) * wght
               Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + t
            endif
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
            if (dsqrt(gamma).gt.tol_rho)then
               x = dsqrt(gamma) / rho43
            else
               x = 0d0
            endif
c
            g = -ux(x,gamma86)**p
            dg = -p*ux(x,gamma86)**(p-1)*uxp(x,gamma86)
c
c
            func(n) = func(n) + rho43*g*wght
            Amat(n,D1_RA) = Amat(n,D1_RA)+(4d0/3d0)*rho13*(g-x*dg)*wght
c     
            if (x.gt.tol_rho) then
               t = dg / sqrt(gamma) * wght
               Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + t * 0.5d0
            endif
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
            if (dsqrt(gamma).gt.tol_rho)then
               x = dsqrt(gamma) / rho43
            else
               x = 0d0
            endif
            g = -ux(x,gamma86)**p
            dg = -p*ux(x,gamma86)**(p-1)*uxp(x,gamma86)
c
c
            func(n) = func(n) +rho43*g*wght
            Amat(n,D1_RB) = Amat(n,D1_RB)+(4d0/3d0)*rho13*(g-x*dg)*wght
c     
            if (x.gt.tol_rho) then
               t = dg / sqrt(gamma) * wght
               Cmat(n,D1_GBB) = Cmat(n,D1_GBB) + t * 0.5d0
            endif
c
c
 20      continue
c
      endif
c
      return
      end
C> @}
