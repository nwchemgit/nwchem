      subroutine dft_tstamp (name)
C$Id: dft_tstamp.F,v 1.5 1995-11-22 07:32:06 d3e129 Exp $
      implicit real*8 (a-h,o-z)
#include "global.fh"
#include "tcgmsg.fh"
#include "stdio.fh"
      character*8 istime,isdate
      character*(*) name
#if defined(CRAY)
      character*24 datime
      call date(datime)
      isdate=datime(1:9)
      call clock(istime)
#else
#if defined(USE_FDATE)
      character*24 datime
      call fdate(datime)
      isdate=datime(4:11)
      istime=datime(12:19)
#else
c                1         2    
c:     :123456789012345678901234
c:ctime:Mon Aug 31 16:02:05 1992
c:fdate:Mon Nov 9 10:48:45 1992
c
      character*26 datime
      call util_date(datime)
      isdate=datime(4:11)
      istime=datime(12:19)
#endif
#endif
      if( ga_nodeid().eq.0 )then
        write (LuOut,100) name,isdate,istime
        call flush(LuOut)
      endif
      call ga_sync
  100 format (/,/,1x,'    <<<<<<   ',a25,'   >>>>>>',/,
     &            5x,'<<<<<< date: ',a8,' time: ',a8,' >>>>>>',/)
      return
      end
