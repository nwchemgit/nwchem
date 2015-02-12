      subroutine dfill(n,val,a,ia)
C$Id$
      implicit none
      integer n, ia
      double precision val, a(*)
      integer i
c
c     initialise double precision array to scalar value
c
#if 1
         call dcopy(n,val,0,a,ia)
#else
      if (ia.eq.1) then
         do 10 i = 1, n
            a(i) = val
 10      continue
      else
         do 20 i = 1,(n-1)*ia+1,ia
            a(i) = val
 20      continue
      endif
#endif
c
      end
