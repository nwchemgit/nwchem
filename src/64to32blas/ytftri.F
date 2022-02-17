      subroutine ytftri( transr, uplo, diag, n, a, info)
*
* $Id$
*
      implicit none
#include "y64.fh"
      integer            info, n
      character          transr, uplo, diag
      double precision   a( * )
c
      INTGR4  n4,info4
c
      n4=n
c
      call dtftri( transr, uplo, diag, n4, a, info4)
      
      info=info4
      return
      end
