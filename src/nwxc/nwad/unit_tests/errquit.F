C
C A dummy implementation of ERRQUIT needed because the multivariate
C implementation may call this. Of course when this gets called you are
C in trouble, but making sure that the routine at least exists avoids
C linking problems for the unit_tests.
C
      subroutine errquit(c,i1,i2)
      implicit none
      character*(*) c
      integer i1,i2
      character*6 c2
      write(*,*)"ERROR: ",c,i1,i2
      stop
      end
