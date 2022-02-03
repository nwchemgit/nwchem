      subroutine ypftri( transr, uplo, n, a, info)
*
* $Id$
*
      implicit none
#include "y64.fh"
      integer            info, n
      character          transr, uplo
      double precision   a( * )
c
      INTGR4  n4,info4
c
      n4=n
c
      call dpftri( transr, uplo, n4, a, info4)
      
      info=info4
      return
      end
