c $Id: paw_gauss_weights.F,v 1.3 2003-10-21 02:05:18 marat Exp $
c
c $Log: not supported by cvs2svn $
c Revision 1.2  2003/02/06 06:15:44  marat
c ...
c
c

!*************************************************
!  
!     Name    : gauss_weights
!  
!  
!     Purpose : generates the weights and coordinates for
!               Gauss-Legendre integration.
!  
!  
!     Notes : double is presision is required
!  
!*************************************************
      subroutine gauss_weights(x1, x2, x, w, n)
      
      integer n
      double precision x1, x2
      double precision x(*), w(*)
      
!    *** LOCAL VARIABLES**
      integer i, j, m, niter
      double precision eps
      parameter (eps = 3.e-14)
      double precision p1, p2, p3, pp, xl, xm, z, z1
      
      m = (n + 1)/2
      xm = 0.5*(x2 + x1)
      xl = 0.5*(x2 - x1)
      
      do i = 1, m
      
         z = cos(3.141592654*(i - .25)/(n + 0.5))
         
         niter = 0
    1    continue
    
         niter = niter + 1
         if (niter .ge. 1000000)
     >     call errquit('cannot converge in gauss_weights',0,1)
     
         p1 = 1.0
         p2 = 0.0

         do j = 1, n

            p3 = p2
            p2 = p1
            p1 = ((2.0*j - 1.0)*z*p2 - (j - 1.0)*p3)/j

         end do

         pp = n*(z*p1 - p2)/(z*z - 1.0)
         z1 = z
         z = z1 - p1/pp

         if (abs(z - z1) .gt. eps) go to 1
         x(i) = xm - xl*z
         x(n+1-i) = xm + xl*z
         w(i) = 2.0*xl/((1.0 - z*z)*pp*pp)
         w(n+1-i) = w(i)

      end do

      return

      end

