      subroutine grid_setspac(rtdb, xyz, xyzw, Rij, rq, p, ictr, nq)
c
c$Id: grid_setspac.F,v 1.2 1999-04-29 16:36:58 d3h449 Exp $
c
      implicit none
      integer rtdb
c
#include "cdft.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "mafdecls.fh"
#include "cgridfile.fh"
c
      integer ictr
      integer nq
      double precision xyzw(4,nq), xyz(3,ncenters)
      double precision Rij((ncenters*(ncenters+1))/2)
      double precision rq(nq,ncenters)
      double precision p(nq,ncenters)
c
c     should have all possible spatial weight possibilities here ... for now 
c     try Becke weights.
c
      call grid_beckew(xyzw, xyz, rij, rq, p, ictr, nq)
c
      return
      end 
