      subroutine ma_print_compact(vector,ndim,mdim,title)
c
c $Id: ma_print_compact.F,v 1.2 2004-10-28 18:16:40 edo Exp $
c

      implicit none
#include "stdio.fh"
      integer ndim,mdim
      double precision vector(ndim,mdim)
      character*(*) title
      integer i,j

      write(LuOut,*)
      write(LuOut,*) title
      write(LuOut,9000)
      do i = 1,mdim
        do j = 1,ndim
          if (dabs(vector(j,i)).gt.1.0d-14)
     1      write(LuOut,9010) j,i,vector(j,i)
        enddo
      enddo
 9000 format('-------------------------')
 9010 format(2i5,1X,20e17.10)
      return
      end
