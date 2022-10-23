      subroutine  cdcopy(n,dx,incx,dy,incy)
c
c     copies a vector, x, to a vector, y.
c     uses rolled loops for increments equal to one.
c     written by ron shepard, based on cdcopy written by
c     jack dongarra, linpack, 3/11/78.

c     need to avoid blas version for integer rep incompatibility
c     with columbus

      real*8 dx(*),dy(*)
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if ( n.le.0 ) return
      if ( incx.eq.1 .and. incy.eq.1 ) go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if ( incx.lt.0 ) ix = (-n+1)*incx + 1
      if ( incy.lt.0 ) iy = (-n+1)*incy + 1
      do 10 i = 1, n
         dy(iy) = dx(ix)
         ix = ix + incx
         iy = iy + incy
10    continue
      return
c
c        code for both increments equal to 1
c
20    continue
      do 50 i = 1, n
         dy(i) = dx(i)
50    continue
      return
      end
