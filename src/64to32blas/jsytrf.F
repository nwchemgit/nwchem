      subroutine jsytrf( uplo, n, a, lda, ipiv, work, lwork, info)
*
* $Id$
*
      implicit none
#include "y64.fh"
      character          uplo
      integer            info, lda, lwork, n
      integer            ipiv( * )
      double complex   a( lda, * ), work( * )
c
      INTGR4 n4,lda4,lwork4,info4
c
      n4=n
      lda4=lda
      lwork4=lwork
      call zsytrf( uplo, n4, a, lda4, ipiv, work,
     $                  lwork4, info4 )
      info=info4
      return
      end

