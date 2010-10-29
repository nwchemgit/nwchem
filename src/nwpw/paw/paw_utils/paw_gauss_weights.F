c $Id$
c
c $Log: not supported by cvs2svn $
c Revision 1.3  2003/10/21 02:05:18  marat
c switched to new errquit by running global replace operation
c see the script below (note it will not work on multiline errquit calls)
c *********************************************************
c #!/bin/sh
c
c e=`find . -name "*F" -print`
c
c for f in $e
c do
c cp $f $f.bak
c sed  's|\(^[ ].*call[ ]*errquit([^,]*\)\(,[^,]*\)\()\)|\1,0\2\3|' $f.bak > $f
c #rm $f.bak
c done
c **********************************************************
c
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
      parameter (eps = 3.d-14)
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
            p1 = ((2.0d0*j - 1.0d0)*z*p2 - (j - 1.0d0)*p3)/j

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



!*************************************************
!
!     Name    : gauss_chebyshev_weights
!
!
!     Purpose : generates the weights and coordinates for
!               Gauss-Chebyshev integration.
!
!
!     Notes : double is presision is required
!
!*************************************************
      subroutine gauss_chebyshev_weights(x, w, n)
      implicit none
      double precision x(*), w(*)
      integer n

      integer i
      double precision pi

      pi = 4.0d0*datan(1.0d0)
      do i=1,n
         w(i) = pi/dble(n)
         x(i) = dcos( ((2*i-1)/dble(2*n))*pi)
      end do

      return
      end


