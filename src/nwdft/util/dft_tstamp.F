      subroutine dft_tstamp (name)
C$Id: dft_tstamp.F,v 1.2 1995-06-20 17:46:45 og845 Exp $
      implicit real*8 (a-h,o-z)
#include "global.fh"
#include "tcgmsg.fh"
      character*8 istime,isdate
      character*(*) name
#if defined(CRAY)
      character*24 datime
      call date(datime)
      isdate=datime(1:9)
      call clock(istime)
#else
      character*24 datime
      call fdate(datime)
      isdate=datime(4:11)
      istime=datime(12:19)
#endif
      if( ga_nodeid().eq.0 )then
        write (6,100) name,isdate,istime
        call flush(6)
      endif
      call ga_sync
  100 format (/,/,1x,'    <<<<<<   ',a25,'   >>>>>>',/,
     &            5x,'<<<<<< date: ',a8,' time: ',a8,' >>>>>>',/)
      return
      end
