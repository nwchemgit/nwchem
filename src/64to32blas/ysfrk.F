      subroutine ysfrk( transr, uplo, trans, n, k, alpha, a, lda, beta,
     $                  c )
*
* $Id$
*
      implicit none
#include "y64.fh"
      double precision   alpha, beta
      integer            k, lda, n
      character          trans, transr, uplo
      double precision   a( lda, * ), c( * )
c
      INTGR4  n4,k4,lda4
c
      n4=n
      k4=k
      lda4=lda
c
      call dsfrk( transr, uplo, trans, n4, k4, alpha, a, lda4, beta,
     $                  c )
      return
      end
