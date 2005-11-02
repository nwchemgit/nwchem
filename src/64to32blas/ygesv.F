      SUBROUTINE YGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
*
* $Id: ygesv.F,v 1.1 2005-11-02 17:14:00 edo Exp $
*     .. Scalar Arguments ..
      implicit none
#include "y64.fh"
      INTEGER            INFO, LDA, LDB, N, NRHS
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
c
      INTGR4            INFO4, LDA4, LDB4, N4, NRHS4
C      integer maxn
C      parameter (maxn=10000)
C      INTGR4            IPIV4(maxn)
      integer i
c
cedo      if(n.gt.maxn) then
cedo         write(6,*) ' ygesv: increase maxn to: ',n
cedo         call errquit(' ygesv: maxn value ',0,0)
cedo      endif
      LDA4=lda
      ldb4=ldb
      n4=n
      nrhs4=nrhs
c
      call DGESV( N4, NRHS4, A, LDA4, IPIV, B, LDB4, INFO4 )
c
cedo      do i=1,N
cedo         ipiv(i)=ipiv4(i)
cedo      enddo
      info=info4
      RETURN
      END
