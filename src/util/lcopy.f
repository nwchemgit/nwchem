      Subroutine lcopy(n,src,ix,dest,iy)
*
* $Id#
*
      implicit none
      integer n, ix, iy
      logical src(*), dest(*)
      integer i,ii,jj

      ii = 1
      jj = 1
      do i=1,n
        dest(jj) = src(ii)
        ii = ii + ix
        jj = jj + iy
      enddo
      return
      end

c $Id$
