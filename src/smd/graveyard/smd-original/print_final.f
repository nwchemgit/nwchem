      SUBROUTINE print_final()

      implicit none

      include 'p_input.inc'
      include 'p_array.inc'
      include 'cm_atom.inc'

      integer  i

      write(output,'(/,a4,1x,i12)')'COOR',natms

      do i=1,natms

       write(output,'(i3,3(1x,g20.10))')atmtype(i),
     $  ccc(i,1),ccc(i,2),ccc(i,3)

      enddo

      return

      END
c $Id$
