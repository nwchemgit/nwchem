      subroutine cosmo_fibonacci(n,xyz,dsurf,dvol)
      implicit none
      integer n
      double precision xyz(3,n)
      double precision dsurf, dvol

      integer i,ii
      double precision eps, goldenratio, twopigr, phi, theta, sint
      parameter(goldenratio = 0.5d0 + 0.5d0*sqrt(5.0d0) )
      twopigr = 2d0*acos(-1d0)/goldenratio

      dsurf = 4d0*acos(-1d0)/dble(n)
      dvol  = dsurf/3d0

      do i=0,n-1
        ii = i+1
        xyz(3,ii) = 1d0 - (2d0*i+1d0)/dble(n)
        phi = twopigr * i
        theta = acos(xyz(3,ii))
        xyz(1,ii) = cos(phi)*sin(theta)
        xyz(2,ii) = sin(phi)*sin(theta)
      enddo

      end subroutine
