      Subroutine qdist(rq0, rq, qxyz, xyz, nq, ncenters)
c
C$Id$
c
      implicit none
      integer nq, ncenters
c
c     Distance Squared Between Quadrature Points & Integration Centers
c
      double precision rq(nq,ncenters)
      double precision rq0(ncenters)
c
c     Cartesian Coordinates of Quadrature Points
c
      double precision qxyz(3,nq)
c
c     Cartesian Coordinates of Integration Centers
c
      double precision xyz(3,ncenters)
c
c     Evaluate the distances (squared) between each center and the 
c     sampling points.
c
      integer m, n
      double precision r2min, x, y, z, r2
c
      do 20 m = 1, ncenters
         r2min = 1.D+06
         do 10 n = 1, nq
            x = qxyz(1,n) - xyz(1,m)
            y = qxyz(2,n) - xyz(2,m)
            z = qxyz(3,n) - xyz(3,m)
            r2 = x*x + y*y + z*z
            r2min = min(r2,r2min)
            rq(n,m) = r2
   10    continue
         rq0(m) = r2min
   20 continue
      return
      end
