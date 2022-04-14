      subroutine ytfsm( transr, side, uplo, trans, diag, m, n, alpha, a,
     $                  b, ldb )
*
* $Id$
*
      implicit none
#include "y64.fh"
      double precision   alpha
      integer            m, n, ldb
      character          trans, transr, uplo, side, diag
      double precision   a( * ), b(ldb, * )
c
      INTGR4  m4,n4,ldb4
c
      m4=m
      n4=n
      ldb4=ldb
c
      call dtfsm( transr, side, uplo, trans, diag, m4, n4, alpha, a, 
     $                  b, ldb4 )
      return
      end
