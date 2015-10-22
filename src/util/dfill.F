C>
C> \brief Routine to initialize an array (Deprecated)
C>
C> The DCOPY call that implements this routine is the preferred
C> initialization mechanism for performance reasons.
C>
      subroutine dfill(n,val,a,ia)
C$Id: dfill.f 26932 2015-03-09 23:07:47Z d3y133 $
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
