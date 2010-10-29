c***********************************************************************
c
c     subroutine header
c
c     writes out a header to new section
c
c***********************************************************************
c
      subroutine dft_header(message)
C$Id$
      implicit none

      character*(*) message
      character*80 border
#include "stdio.fh"
c
      integer lmsg,i
c
c-->  Determine length of message.
c
      lmsg = len(message)
c
c-->  Define character for border.
c
      do 1000 i=1,lmsg
         border(i:i) = '-'
1000  continue
c
c-->  Write header.
c
      write(LuOut,'(/,5x,a)') border(1:lmsg)
      write(LuOut,'(5x,a)') message
      write(LuOut,'(5x,a)') border(1:lmsg)
c
      return
      end
