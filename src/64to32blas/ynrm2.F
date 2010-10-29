      DOUBLE PRECISION FUNCTION YNRM2 ( N, X, INCX )
*
* $Id$
*
      implicit none
#include "y64.fh"
      INTEGER                           INCX, N
      DOUBLE PRECISION                  X( * )
c
      INTGR4 incx4,n4
      double precision dnrm2
      external dnrm2
c
      n4=n
      incx4=incx
c
      ynrm2=DNRM2 ( N4, X, INCX4 )
c
      return
      end
