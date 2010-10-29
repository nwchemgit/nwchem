      integer function idamin(n,dx,incx)
*
* $Id$
*
c
c     finds the index of element having min. absolute value.
C     Trivial modification of BLAS idamax routine!
c
      double precision dx(*),dmin
      integer i,incx,ix,n
c
      idamin = 0
      if( n.lt.1 .or. incx.le.0 ) return
      idamin = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      dmin = dabs(dx(1))
      ix = ix + incx
      do 10 i = 2,n
         if(dabs(dx(ix)).ge.dmin) go to 5
         idamin = i
         dmin = dabs(dx(ix))
    5    ix = ix + incx
   10 continue
      return
c
c        code for increment equal to 1
c
   20 dmin = dabs(dx(1))
      do 30 i = 2,n
         if(dabs(dx(i)).ge.dmin) go to 30
         idamin = i
         dmin = dabs(dx(i))
   30 continue
      return
      end
