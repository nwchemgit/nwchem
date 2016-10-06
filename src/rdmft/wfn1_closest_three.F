C> \ingroup wfn1
C> @{
C>
C> \brief Find three points closest to a particular point
C>
C> This routine takes an ordered table of points and an index into the
C> table. Given this information it returns the indeces of the three
C> points that are closest together and where one of the points is the
C> given point. Also \f$xx(ix1) < xx(ix2) < xx(ix3)\f$.
C>
      subroutine wfn1_closest_three(npts,xx,ix,ix1,ix2,ix3)
      implicit none
c
      integer,          intent(in)  :: npts     !< The number of points
      double precision, intent(in)  :: xx(npts) !< The coordinates of
                                                !< the points
      integer,          intent(in)  :: ix       !< The index of the
                                                !< given point
      integer,          intent(out) :: ix1      !< The index of point 1
      integer,          intent(out) :: ix2      !< The index of point 2
      integer,          intent(out) :: ix3      !< The index of point 3
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
      endif
c
      end
C>
C> @}
