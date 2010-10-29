      SUBROUTINE YORGHR( N, ILO, IHI, A, LDA, TAU, WORK, LWORK, INFO )
c $Id$
*
      implicit none
#include "y64.fh"
      INTEGER            IHI, ILO, INFO, LDA, LWORK, N
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( LWORK )
c
      INTGR4            IHI4, ILO4, INFO4, LDA4, LWORK4, N4
c
      IHI4=ihi
      ILO4=ilo
      LDA4=lda
      LWORK4=lwork
      N4=n
c
      call DORGHR( N4, ILO4, IHI4, A, LDA4, TAU, WORK, LWORK4, INFO4 )
c
      info=info4
c
      return
      END
