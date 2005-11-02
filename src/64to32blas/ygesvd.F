      SUBROUTINE YGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT,
     $                   WORK, LWORK, INFO )
* $Id: ygesvd.F,v 1.1 2005-11-02 17:14:00 edo Exp $
      implicit none
#include "y64.fh"
      CHARACTER          JOBU, JOBVT
      INTEGER            INFO, LDA, LDU, LDVT, LWORK, M, N
      DOUBLE PRECISION   A( LDA, * ), S( * ), U( LDU, * ),
     $                   VT( LDVT, * ), WORK( * )
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
      call DGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT,
     $                   WORK, LWORK, INFO )
      info=info4

      RETURN
      END
