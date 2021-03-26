      subroutine ytrtri( uplo, diag, n, a, lda, info )
*
* $id$
*
      implicit none
#include "y64.fh"
      character          diag, uplo
      integer            info, lda, n
*     ..
*     .. array arguments ..
      double precision   a( lda, * )
c
      INTGR4 lda4,info4,n4
c
      lda4=lda
      n4=n
c      
      call dtrtri( uplo, diag, n4, a, lda4, info4 )
c
      info=info4
c      
      return
      end
