      subroutine ma_scale(a,n,scalar)
c
c $Id: ma_scale.F,v 1.1 2002-10-18 23:42:06 sohirata Exp $
c
      implicit none
#include "mafdecls.fh"
#include "tce.fh"
      integer i,n
      double precision a(n)
      double precision scalar
c
      do i = 1,n
        a(i) = a(i) * scalar     
      enddo
      return
      end
