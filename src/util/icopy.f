      Subroutine icopy(n,src,ix,dest,iy)
*
* $Id: icopy.f,v 1.3 1997-10-31 20:45:33 d3e129 Exp $
*
      implicit none
      integer n, src(*), ix, dest(*), iy
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

