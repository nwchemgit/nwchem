C> \ingroup selci
C> @{
      subroutine restk
*
* $Id$
*
      implicit real*8 (a-h,o-z), integer(i-n)
#include "errquit.fh"
      parameter(irest1 = 69,irest2=70, istde=0,istdo=6)
      integer idum(7)
c     
#include "cselcifiles.fh"
c     
      do 00100 i = 1,7
         idum(i) = -i
00100 continue
      open(unit=irest1,
     &     file=names(1),
     &     access='sequential',
     &     form='unformatted',
     &     status='unknown',
     &     err=90911)
      open(unit=irest2,
     &     file=names(2),
     &     access='sequential',
     &     form='unformatted',
     &     status='unknown',
     &     err=90911)
      write(irest1) idum
      write(irest2) idum
      call selci_restkc(irest1,irest2)
      return
90911 call errquit(' error opening restart file ',0, DISK_ERR)
c
      end
      subroutine selci_restkc(irest1,irest2)
      implicit integer(i-z)
c
      close(unit=irest1,status='delete')
      close(unit=irest2,status='delete')
c
      return
      end
C> @}
