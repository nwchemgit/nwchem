      subroutine init70
c $Id$
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
c     Generate tables for evaluation of incomplete gamma function
c
      call tabgen
c
c     Error2 and error1 control the accuracy of 2-e integrals.
c     For use in energy expressions 10^-11 is appropriate.
c     For use in the Schwarz() screening 10^-22 is required.
c     
c     The routine int_sp_acc has been provided so that this
c     value can be changed according to circumstances.
c     
      error2 = 1.0d-22
      error1 = sqrt(error2)
      var(1) = 15.0d0
      var(2) = 5.0d0
c
      end
      subroutine int_sp_acc(tol)
c $Id$
      implicit none
#include "maxc.fh"
      double precision tol
c
      error2 = tol
      error1 = sqrt(error2)
c
      end
