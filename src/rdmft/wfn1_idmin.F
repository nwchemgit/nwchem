C> \ingroup wfn1
C> @{
C>
C> \brief Find minimal value in a vector
C>
C> This routine is similar to idamax but for some reason does not seem
C> to be part of BLAS. It is different in that it returns the index
C> of the element with the smallest value.
C>
      integer function wfn1_idmin(n,d,inc)
      implicit none
      integer,          intent(in) :: n    !< The number of elements in
                                           !< the array
      double precision, intent(in) :: d(*) !< The data array
      integer,          intent(in) :: inc  !< The increment between
                                           !< values in the array
c
      integer indx
      integer i
c
      indx = 1
      do i = 2, n
        if (d(1+(i-1)*inc).le.d(1+(indx-1)*inc)) then
          indx = i
        endif
      enddo
c
      wfn1_idmin = indx
c
      return
      end
C>
C> @}
