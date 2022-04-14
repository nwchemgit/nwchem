      SUBROUTINE jgetrf( M, N,  A, LDA, ipiv, INFO )
*
* $id$
*
      implicit none
#include "y64.fh"
      integer ipiv(*)
      integer            info, lda, m, n
      complex*16   a( lda, * )
c
      INTGR4 m4,n4,info4,lda4
c
      n4=n
      m4=n
      lda4=lda
c
      call zgetrf( m4, n4,  a, lda4, ipiv, info4 )
      info=info4
      return
      end
