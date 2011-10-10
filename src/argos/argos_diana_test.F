      subroutine argos_diana_test(x,nsa)
      integer nsa
      real*8 x(nsa,3)
c
      integer i,j
c
      do 1 i=1,10
      write(*,'(i5,3f12.6)') i,(x(i,j),j=1,3)
    1 continue
c
      return
      end
c $Id$
