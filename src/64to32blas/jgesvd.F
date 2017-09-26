      SUBROUTINE jGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, 
     $                   VT, LDVT,  WORK, LWORK, RWORK, INFO )
* $Id: ygesvd.F 19697 2010-10-29 16:57:34Z d3y133 $
      implicit none
#include "y64.fh"
      CHARACTER          JOBU, JOBVT
      INTEGER            INFO, LDA, LDU, LDVT, LWORK, M, N
      DOUBLE PRECISION   A( LDA, * ), S( * ), U( LDU, * ),
     $                   VT( LDVT, * ), WORK( * ), RWORK(*)
c
      INTGR4            INFO4, LDA4, LDU4, LDVT4, LWORK4, M4, N4
c
      lda4=lda
      ldu4=ldu
      ldvt4=ldvt
      m4=m
      n4=n
      lwork4=lwork
c
      call ZGESVD( JOBU, JOBVT, M4, N4, A, LDA4, S, U, LDU4, 
     $                   VT, LDVT4,  WORK, LWORK4, RWORK, INFO4 )
      info=info4

      RETURN
      END
