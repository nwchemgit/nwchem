      subroutine ifill(n,val,a,ia)
C$Id$
      implicit none
      integer n, val, a(*), ia, i
c
c     initialise integer precision array to scalar value
c
      if (ia.eq.1) then
         do 10 i = 1, n
            a(i) = val
 10      continue
      else
         do 20 i = 1,(n-1)*ia+1,ia
            a(i) = val
 20      continue
      endif
c
      end
      subroutine ifill2(n,val,a,ia)
C$Id$
      implicit none
      integer n, val, ia, i
      integer*2 a(*)
c
c     initialise integer precision array to scalar value
c
      if (ia.eq.1) then
         do 10 i = 1, n
            a(i) = val
 10      continue
      else
         do 20 i = 1,(n-1)*ia+1,ia
            a(i) = val
 20      continue
      endif
c
      end
      subroutine lfill(n,val,a,ia)
C$Id$
      implicit none
      integer n, ia, i
      logical a(*),val
c
c     initialise logical array to false
c
      if (ia.eq.1) then
         do 10 i = 1, n
            a(i) = val
 10      continue
      else
         do 20 i = 1,(n-1)*ia+1,ia
            a(i) = val
 20      continue
      endif
c
      end
