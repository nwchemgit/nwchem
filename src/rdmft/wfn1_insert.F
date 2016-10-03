C> \ingroup wfn1
C> @{
C>
C> \brief Insert a new data point in an ordered table
C>
C> The table is sorted in increasing order. This routine adds another
C> data point to the table.
C>
      subroutine wfn1_insert(mxpts,npts,xx,fx,df,xnew,fnew,dfnew)
      implicit none
c
#include "errquit.fh"
c
      integer, intent(in)    :: mxpts !< The dimension of the table
      integer, intent(inout) :: npts  !< The number of points already in
                                      !< the table, must be greater
                                      !< than 0
      double precision, intent(inout) :: xx(mxpts) !< The coordinates
      double precision, intent(inout) :: fx(mxpts) !< The function values
      double precision, intent(inout) :: df(mxpts) !< The function values
      double precision, intent(in)    :: xnew  !< The coordinate of the
                                               !< new point
      double precision, intent(in)    :: fnew  !< The function value of
                                               !< the new point
      double precision, intent(in)    :: dfnew !< The function value of
                                               !< the new point
c
      integer ii, jj !< Counters
c
      if (npts.ge.mxpts) call errquit("wfn1_insert: table full",0,UERR)
c
      ii = npts
      do while (xnew.lt.xx(ii).and.ii.gt.1)
        ii = ii - 1
      enddo
      if (xnew.lt.xx(ii)) then
        do jj = npts, ii, -1
          xx(jj+1) = xx(jj)
          fx(jj+1) = fx(jj)
          df(jj+1) = df(jj)
        enddo
        xx(ii) = xnew
        fx(ii) = fnew
        df(ii) = dfnew
      else
        do jj = npts, ii+1, -1
          xx(jj+1) = xx(jj)
          fx(jj+1) = fx(jj)
          df(jj+1) = df(jj)
        enddo
        xx(ii+1) = xnew
        fx(ii+1) = fnew
        df(ii+1) = dfnew
      endif
      npts = npts + 1
c
      end
C>
C> @}
