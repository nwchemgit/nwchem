C> \ingroup wfn1
C> @{
C>
C> \brief Find three points closest to a particular point
C>
C> This routine takes an ordered table of points and an index into the table.
C> Given this information it returns the indeces of the three points that are
C> closest together and where one of the points is the given point. Also
C> \f$xx(ix1) < xx(ix2) < xx(ix3)\f$.
C>
      subroutine wfn1_closest_three(npts,xx,ix,ix1,ix2,ix3)
      implicit none
c
      integer npts              !< [Input] The number of points
      double precision xx(npts) !< [Input] The coordinates of the points
      integer ix                !< [Input] The index of the given point
      integer ix1               !< [Output] The index of point 1
      integer ix2               !< [Output] The index of point 2
      integer ix3               !< [Output] The index of point 3
c
      double precision dist1, dist2, dist3
c
      if (npts.eq.3) then
c
c       there is only one possible set of three points
c
        ix1 = 1
        ix2 = 2
        ix3 = 3
c
      else if (ix.eq.1) then
c
c       there are no previous points
c
        ix1 = ix
        ix2 = ix1+1
        ix3 = ix2+1
c
      else if (ix.eq.npts) then
c
c       there are no subsequent points
c
        ix3 = ix
        ix2 = ix3-1
        ix1 = ix2-1
c
      else
c
c       just take the points on either side
c
        ix1 = ix-1
        ix2 = ix
        ix3 = ix+1
c
c     else if (ix.eq.2) then
c
c       ix is either the middle point or the lowest point
c
c       dist1 = abs(xx(ix-1)-xx(ix))+abs(xx(ix+1)-xx(ix))
c       dist2 = abs(xx(ix+1)-xx(ix))+abs(xx(ix+2)-xx(ix))
c       if (dist1.le.dist2) then
c         ix1 = ix - 1
c         ix2 = ix
c         ix3 = ix + 1
c       else
c         ix1 = ix
c         ix2 = ix + 1
c         ix3 = ix + 2
c       endif
c
c     else if (ix.eq.npts-1) then
c
c       ix is either the middle point or the highest point
c
c       dist1 = abs(xx(ix-1)-xx(ix))+abs(xx(ix-2)-xx(ix))
c       dist2 = abs(xx(ix-1)-xx(ix))+abs(xx(ix+1)-xx(ix))
c       if (dist2.le.dist1) then
c         ix1 = ix - 1
c         ix2 = ix
c         ix3 = ix + 1
c       else
c         ix1 = ix - 2
c         ix2 = ix - 1
c         ix3 = ix
c       endif
c
c     else
c
c       ix is either the lowest, the middle or the highest point
c
c       dist1 = abs(xx(ix-2)-xx(ix))+abs(xx(ix-1)-xx(ix))
c       dist2 = abs(xx(ix+1)-xx(ix))+abs(xx(ix-1)-xx(ix))
c       dist3 = abs(xx(ix+2)-xx(ix))+abs(xx(ix+1)-xx(ix))
c       if (dist2.le.dist1.and.dist2.le.dist3) then
c         ix1 = ix - 1
c         ix2 = ix
c         ix3 = ix + 1
c       else if (dist1.le.dist2.and.dist1.le.dist3) then
c         ix1 = ix - 2
c         ix2 = ix - 1
c         ix3 = ix
c       else
c         ix1 = ix
c         ix2 = ix + 1
c         ix3 = ix + 2
c       endif
      endif
c
      end
C>
C> @}
