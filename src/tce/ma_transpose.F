      subroutine ma_transpose(a,n)
c
c $Id: ma_transpose.F,v 1.1 2002-10-18 23:42:07 sohirata Exp $
c
      implicit none
#include "mafdecls.fh"
#include "tce.fh"
      integer i,j,n
      double precision a(n,n)
      double precision swap
c
      do i = 1,n
        do j = i,n
          swap = a(i,j)
          a(i,j) = a(j,i)      
          a(j,i) = swap
        enddo
      enddo
      return
      end
