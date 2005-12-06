      SUBROUTINE YPOTRF( UPLO, N, A, LDA, INFO )
c
* $Id: ypotrf.f,v 1.1 2005-12-06 16:59:50 edo Exp $
c
      CHARACTER          UPLO
      INTEGER            INFO, LDA, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
c
      integer*4 n4,lda4,info4
c
      n4=n
      lda4=lda
c
      call DPOTRF( UPLO, N4, A, LDA4, INFO4 )
      info=info4
      return
      end
