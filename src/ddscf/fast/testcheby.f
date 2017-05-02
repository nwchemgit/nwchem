      program tester
*
* $Id$
*
      implicit none
c
      integer n, order
      parameter (n=20, order=13)
      double precision g(n,n,n), f(n)
      integer i, j, k
      double precision x, y, z, value, tn_interp_3d_point, test
      double precision util_random, tn_cube_eval, xx, yy, zz, maxerr
c
      do i = 1, n
         f(i) = 1 + 0.5d0*(i-1 + 0.5d0*(i-1))
      enddo
      call tn_lsq_fit(n,1,order,f)
c
      do k = 1, n
         do j = 1, n
            do i = 1, n
               g(i,j,k) = exp(-dble(i+j+k)**2 / dble(n)**2)
            enddo
         enddo
      enddo
c
      call tn_lsq_fit_cube(n,order,g)
c
      maxerr = 0.0d0
      do i = 1, 2000
         x = util_random(0)*n + 1
         y = util_random(0)*n + 1
         z = util_random(0)*n + 1
*         value = tn_interp_3d_point(g,n,n,n,x,y,z,order)
         test = exp(-dble(x+y+z)**2 / dble(n)**2)

c
         xx = 2.0d0*(x - 1.0d0 - dble(n-1)*0.5d0)/dble(n-1)
         yy = 2.0d0*(y - 1.0d0 - dble(n-1)*0.5d0)/dble(n-1)
         zz = 2.0d0*(z - 1.0d0 - dble(n-1)*0.5d0)/dble(n-1)
c
         value = tn_cube_eval(g,n,n,order,xx,yy,zz)
*         write(6,*) xx, yy, zz


         write(6,1) x, y, z, value, test, value-test
 1       format(3f8.4,2x,2f8.4,2x,1p,d9.2)
c
         maxerr= max(maxerr,abs(value-test))

      enddo
      write(6,*) ' maxerr ', maxerr
c
      end
