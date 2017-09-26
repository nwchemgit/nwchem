      SUBROUTINE jGEMM (TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB,
     $                   BETA, C, LDC )
C     $Id: xgemm.F 19697 2010-10-29 16:57:34Z d3y133 $
      implicit none
#include "y64.fh"
*     .. Scalar Arguments ..
      CHARACTER*1        TRANSA, TRANSB
      INTEGER            M, N, K, LDA, LDB, LDC
      DOUBLE COMPLEX   ALPHA, BETA
*     .. Array Arguments ..
      DOUBLE COMPLEX   A( LDA, * ), B( LDB, * ), C( LDC, * )
c
      INTGR4            M4, N4, K4, LDA4, LDB4, LDC4
c
      M4=M
      N4=N
      K4=K
      LDA4=LDA
      LDB4=LDB
      LDC4=ldc
c
      CALL ZGEMM(TRANSA, TRANSB, M4, N4, K4, ALPHA, A, LDA4, B, LDB4,
     $                   BETA, C, LDC4 )
      RETURN
      END
