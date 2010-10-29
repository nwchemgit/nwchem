      SUBROUTINE YLASCL( TYPE, KL, KU, CFROM, CTO, M, N, A, LDA, INFO )
c
* $Id$
c
      implicit none
#include "y64.fh"
      CHARACTER          TYPE
      INTEGER            INFO, KL, KU, LDA, M, N
      DOUBLE PRECISION   CFROM, CTO
      DOUBLE PRECISION   A( LDA, * )
c
      INTGR4            INFO4, KL4, KU4, LDA4, M4, N4
c
      KL4=kl
      KU4=ku
      LDA4=lda 
      M4=m 
      N4=n
c
      call DLASCL( TYPE, KL4, KU4, CFROM, CTO, M4, N4, A, LDA4, INFO4 )
      info=info4
      return
      end
