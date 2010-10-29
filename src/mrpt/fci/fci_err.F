*
* $Id$
*
      subroutine fci_err(string, code)
      implicit none
      character*(*) string
      integer code
c
      write(6,1) string, code
 1    format(/
     $     ' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'/
     $     1x, a, 1x, i8/
     $     ' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'/
     $     )
c
      stop 99
c
      end
