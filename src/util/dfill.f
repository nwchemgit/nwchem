      subroutine dfill(n,val,a,ia)
C$Id$
      implicit none
      integer n, ia
      double precision val, a(*)
      integer i
c
c     initialise double precision array to scalar value
c
         call dcopy(n,val,0,a,ia)
c
      end
