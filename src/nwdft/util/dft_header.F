c***********************************************************************
c
c     subroutine header
c
c     writes out a header to new section
c
c***********************************************************************
c
      subroutine dft_header(message)
C$Id: dft_header.F,v 1.1 1995-04-27 03:16:43 og845 Exp $
      implicit real*8 (a-h,o-z)
      character*(*) message
      character*80 border
#include "stdio.fh"
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
