*
* $Id$
*
      subroutine screen(thresh,x,ix,n)
      implicit none
      double precision thresh, x(*)
      integer ix, n, i, itop
c     
c     zero elements of x whose absolute value is less than thresh
c     
      if (ix.eq.1) then
         do 10 i = 1,n
            if (abs(x(i)).lt.thresh) x(i) = 0.0d0
 10      continue
      else
         itop = (n-1)*ix + 1
         do 20 i = 1,itop,ix
            if (abs(x(i)).lt.thresh) x(i) = 0.0d0
 20      continue
      endif
c     
      return
      end
