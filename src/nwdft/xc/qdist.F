      Subroutine qdist(rq0,rq,qxyz,xyz,nq,ncenters)

      Implicit real*8 (a-h,o-z)
      Implicit integer (i-n)

c--> Distance Squared Between Quadrature Points & Integration Centers

      Dimension rq(nq,ncenters)
      Dimension rq0(ncenters)

c--> Cartesian Coordinates of Quadrature Points

      Dimension qxyz(3,nq)

c--> Cartesian Coordinates of Integration Centers

      Dimension xyz(3,ncenters)
c
c Evaluate the distances (squared) between each center and the sampling points.
c
c***************************************************************************

      do 20 m = 1,ncenters

       r2min = 1.D+06

       do 10 n = 1,nq

        x = qxyz(1,n) - xyz(1,m)
        y = qxyz(2,n) - xyz(2,m)
        z = qxyz(3,n) - xyz(3,m)

        r2 = x*x + y*y + z*z

        r2min = min(r2,r2min)

        rq(n,m) = r2

   10  continue

       rq0(m) = r2min

   20 continue

      end
