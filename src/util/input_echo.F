      subroutine input_echo(ir,iw)
C$Id: input_echo.F,v 1.1 1996-07-26 16:43:39 d3g681 Exp $
      implicit none
#include "inp.fh"
      character*256 line
      integer ir, iw
c
      integer lenline
c
      rewind ir
      write(iw,*)
      write(iw,10000)
10000 format(/,/,30('='),' echo of input deck ',30('='))
00001 continue
      read(ir,'(a)',err=00002,end=00002)line
      lenline = inp_strlen(line)
      write(iw,'(a)')line(1:lenline)
      goto 00001
00002 continue
      write(iw,10001)
10001 format(80('='),/,/)
      rewind ir
      call inp_init(ir,iw)
      end
