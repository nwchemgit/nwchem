      SUBROUTINE ygetri( N,  A, LDA, ipiv, WORK, LWORK, INFO )
*
* $id$
*
      implicit none
#include "y64.fh"
      integer ipiv(*)
      integer            info, lda, n, lwork
      double precision   a( lda, * )
      double precision work(*)
c
      INTGR4 n4,info4,lda4,lwork4
c
      n4=n
      lda4=lda
      lwork4=lwork
c
      call dgetri( n4,  a, lda4, ipiv, work, lwork4, info4 )
      info=info4
      return
      end
