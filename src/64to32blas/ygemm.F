      SUBROUTINE YGEMM (TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB,
     $                   BETA, C, LDC )
C     $Id$
      implicit none
#include "y64.fh"
      CHARACTER*1        TRANSA, TRANSB
      INTEGER            M, N, K, LDA, LDB, LDC
      DOUBLE PRECISION   ALPHA, BETA
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( LDC, * )
c
      INTGR4            M4, N4, K4, LDA4, LDB4, LDC4
c
      M4=M
      N4=N
      K4=K
      LDA4=LDA
      LDB4=LDB
      LDC4=ldc
#ifdef USE_YFLOP
      flops_ycount = flops_ycount + m*n*k*2
#endif
c
      CALL DGEMM(TRANSA, TRANSB, M4, N4, K4, ALPHA, A, LDA4, B, LDB4,
     $                   BETA, C, LDC4 )
      RETURN
      END
