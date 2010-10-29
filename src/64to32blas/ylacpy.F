      SUBROUTINE YLACPY( UPLO, M, N, A, LDA, B, LDB )
c
* $Id$
c
      implicit none
#include "y64.fh"
      CHARACTER          UPLO
      INTEGER            LDA, LDB, M, N
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
c
      INTGR4            LDA4, LDB4, M4, N4
c
      lda4=lda
      ldb4=ldb
      m4=m
      n4=n
c
      call DLACPY( UPLO, M4, N4, A, LDA4, B, LDB4 )
c      
      return
      end
