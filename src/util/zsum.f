      double complex function zsum(n,sx,incx)
*
* $Id: zsum.f 19707 2010-10-29 17:59:36Z d3y133 $
*
c
c     takes the sum of the array elements.
c     uses unrolled loops for increment equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c
      double complex  sx(*),stemp
      integer i,incx,m,mp1,n,nincx
c
      zsum = cmplx(0.0d0,0.0d0)
      stemp = zsum
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      nincx = n*incx
      do 10 i = 1,nincx,incx
        stemp = stemp + sx(i)
   10 continue
      zsum = stemp
      return
c
c        code for increment equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,6)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        stemp = stemp + sx(i)
   30 continue
      if( n .lt. 6 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,6
        stemp = stemp + sx(i) + sx(i + 1) + sx(i + 2)
     *  + sx(i + 3) + sx(i + 4) + sx(i + 5)
   50 continue
   60 zsum = stemp
      return
      end
