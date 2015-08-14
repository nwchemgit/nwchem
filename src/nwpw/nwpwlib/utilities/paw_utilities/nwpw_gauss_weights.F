c $Id: nwpw_xc.F 26658 2015-01-06 23:07:46Z bylaska $

*     *****************************************************
*     *                                                   *
*     *          nwpw_gauss_weights                       *
*     *                                                   *
*     *****************************************************
      subroutine nwpw_gauss_weights(x1, x2, x, w, n)
      implicit none
      integer n
      double precision x1, x2
      double precision x(*), w(*)

!    *** local variables ***
      integer i, j, m, niter
      double precision eps
      parameter (eps = 3.0d-14)
      double precision p1, p2, p3, pp, xl, xm, z, z1,pi

      pi = 4.0d0*datan(1.0d0)
      m = (n + 1)/2  
      xm = 0.5d0*(x2 + x1)
      xl = 0.5d0*(x2 - x1)

      do i = 1, m
         z = cos(pi*(i - 0.25d0)/(n + 0.5d0))
         niter = 0
    1    continue

         niter = niter + 1
         if (niter .ge. 1000000)
     >     call errquit('cannot converge in gauss_weights',0,1)

         p1 = 1.0d0
         p2 = 0.0d0
         do j = 1, n
            p3 = p2
            p2 = p1
            p1 = ((2.0d0*j - 1.0d0)*z*p2 - (j - 1.0d0)*p3)/dble(j)
         end do

         pp = n*(z*p1 - p2)/(z*z - 1.0d0)
         z1 = z
         z = z1 - p1/pp

         if (abs(z - z1) .gt. eps) go to 1
         x(i) = xm - xl*z
         x(n+1-i) = xm + xl*z
         w(i) = 2.0d0*xl/((1.0d0 - z*z)*pp*pp)
         w(n+1-i) = w(i)
      end do
      return
      end

