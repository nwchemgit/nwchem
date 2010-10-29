c $Id$

      subroutine util_set_ma_char_arrays(nt,tin,tout)
      implicit none
      integer nt
      character*1 tin(nt)
      character*1 tout(nt)

      integer i

      do i=1,nt
       tout(i) = tin(i) 
      end do

      end

