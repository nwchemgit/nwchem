      subroutine util_print_centered(unit, string, center, ounder)
C$Id: print_center.F,v 1.3 1995-02-02 18:09:46 d3g681 Exp $
      implicit none
#include "inp.fh"
      integer unit              ! [input] fortran unit for output
      character*(*) string      ! [input] string to print
      integer center            ! [input] column to center on
      logical ounder            ! [input] if true then underline
c
c     write the string to speficied fortran unit centered about
c     the given column.  optionally underline.
c
      character*80 blank, underline
      integer empty, length, i
c
      blank = ' '
      length = inp_strlen(string)
      empty = center - length/2
      if (empty .lt. 0) empty = 0
      if (empty .gt. 80) empty = 80
      write(unit, 1) blank(1:empty), string(1:length)
 1    format(1x,a,a)
      if (ounder) then
         length = min(80,length)
         do i = 1, length
            underline(i:i) = '-'
         enddo
         write(unit, 1) blank(1:empty), underline(1:length)
      endif
c
      end
