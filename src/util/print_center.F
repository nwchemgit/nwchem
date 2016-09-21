      subroutine util_print_centered(unit, string, center, ounder)
C$Id$
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
#define BLNKLGTH 80
      character*BLNKLGTH blank, underline
      integer empty, length, i
c
      do i=1,BLNKLGTH
         blank(i:i) = ' '
      enddo
      length = inp_strlen(string)
      empty = center - length/2
      if (empty .lt. 0) empty = 0
      if (empty .gt. BLNKLGTH) empty = BLNKLGTH
      write(unit, 1) blank(1:empty), string(1:length)
    1 format(1x,a,a)
      if (ounder) then
         length = min(80,length)
         do i = 1, length
            underline(i:i) = '-'
         enddo
         write(unit, 1) blank(1:empty), underline(1:length)
      endif
c
      end
