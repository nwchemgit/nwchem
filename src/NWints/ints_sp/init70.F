      subroutine init70
c $Id: init70.F,v 1.1 1994-04-06 23:51:55 d3e129 Exp $
      implicit none
c
c     called once by genr70 to initialize internal data
c     structures for the rotated axis integrals
c
#include "maxc.fh"
#include "auxvar.fh"
#include "picon.fh"
      double precision pi
c
      pidiv4 = atan(1.0d0)
      pi = 4.0d0*pidiv4
      pito52 = 2.0d0*pi**2.5d0
      root3 = sqrt(3.0d0)
      root5 = sqrt(5.0d0)
c
      call tabgen
c
      error2 = 1.0d-11
      error1 = sqrt(error2)
      var(1) = 15.0d0
      var(2) = 5.0d0
c
      end
