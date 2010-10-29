*
* $Id$
*
      subroutine dfill(val,x,iskip,n)
      implicit none
      double precision val, x(*)
      integer iskip, n, ipt, i
c      
      ipt=1
      do 10 i=1,n
          x(ipt)=val
10        ipt=ipt+iskip
      return
      end
