!$Id$
!
!$Log:

      subroutine get_spher_grid(ntheta,nphi,angle_phi,
     >                          cos_theta,w_theta,w_phi)
      implicit none

      integer ntheta
      integer nphi
      double precision cos_theta(ntheta)
      double precision angle_phi(nphi)
      double precision w_theta(ntheta)
      double precision w_phi(nphi)

#include "paw_params.fh"

*     *** local variables ***
      integer i

*     *** gaussian quadrature angular grid for cos_theta ***
      call gauss_weights(-1.0d0,1.0d0,cos_theta,w_theta,ntheta)

      if (nphi.gt.1) then
*       *** linear angular grid for angle_phi***
        do i=1,nphi
         angle_phi(i) = 2.0d0*pi*(i-1)/(nphi-1)
         w_phi(i) = 2.0d0*pi/(nphi-1)
        end do
        w_phi(1)    = 0.5d0*w_phi(1)
        w_phi(nphi) = w_phi(1)
      else
        angle_phi(1) = 0.0d0
        w_phi(1)     = 2.0d0*pi
      end if

      end
