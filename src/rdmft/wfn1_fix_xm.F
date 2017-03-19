C> \ingroup wfn1
C> @{
C>
C> \brief Make sure the new x-position is not known already
C>
C> Sometimes a x-position is suggested that has been tried
C> before already. Recomputing this result does not provide any
C> new information. Hence if such a value is proposed it is 
C> replaced by something else.
C>
      subroutine wfn1_fix_xm(npts,xx,xm)
      implicit none
      integer,          intent(in)    :: npts
      double precision, intent(in)    :: xx(npts)
      double precision, intent(inout) :: xm
c
      double precision, parameter :: tol = 1.0d-8
c
      integer ii
c
      do ii = 1, npts
        if (abs(xx(ii)-xm).lt.tol) then
          if (ii.eq.1) then
            xm = xm + xx(ii)-xx(ii+1)
          else if (ii.eq.npts) then
            xm = xm + xx(ii)-xx(ii-1)
          else
            if (xx(ii).gt.xm) then
              xm = 0.50d0*(xx(ii-1)+xx(ii))
            else
              xm = 0.50d0*(xx(ii)+xx(ii+1))
            endif
          endif
        endif
      enddo
      end subroutine wfn1_fix_xm
C> @}
