      SUBROUTINE YLASET( UPLO, M, N, ALPHA, BETA, A, LDA )
c
* $Id: ylaset.F,v 1.1 2005-11-02 17:14:01 edo Exp $
c
      implicit none
#include "y64.fh"
      CHARACTER          UPLO
      INTEGER            LDA, M, N
      DOUBLE PRECISION   ALPHA, BETA
      DOUBLE PRECISION   A( LDA, * )
c
      INTGR4            LDA4, M4, N4
c
      lda4=lda
      m4=m
      n4=n
c
      call DLASET( UPLO, M4, N4, ALPHA, BETA, A, LDA4 )
*
      RETURN
      END
