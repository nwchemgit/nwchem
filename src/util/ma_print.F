      subroutine ma_print(vector,ndim,mdim,title)
c
c $Id$
c
c Printout an MA or other vector
c Written by So Hirata, Feb 2002.
c (c) Battelle, PNNL, 2002.

      implicit none
#include "stdio.fh"
      integer ndim,mdim
      double precision vector(ndim,mdim)
      character*(*) title
      integer i,j,k

      write(LuOut,*)
      write(LuOut,*) title
      write(LuOut,9000)
      do k = 0,mdim/6
        write(LuOut,*)
        write(LuOut,9010) (j,j=k*6+1,min(mdim,k*6+6))
        do i=1,ndim
          write(LuOut,9020) i,(vector(i,j),j=k*6+1,min(mdim,k*6+6))
        enddo
      enddo
      write(LuOut,9000)
      write(LuOut,*)
 9000 format('-------------------------')
 9010 format(5x,20i10)
 9020 format(i5,20f10.5)
      return
      end
