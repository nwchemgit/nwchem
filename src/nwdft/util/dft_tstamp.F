      subroutine dft_tstamp (name)
C$Id$
      implicit none


#include "global.fh"
c$$$#include "tcgmsg.fh"
#include "stdio.fh"
      character*8 istime,isdate
      character*(*) name
c
c     RJH. Commented out unnecessary machine dependent code.
c 
c$$$#if defined(CRAY)
c$$$      character*24 datime
c$$$      call date(datime)
c$$$      isdate=datime(1:9)
c$$$      call clock(istime)
c$$$#else
c$$$#if defined(USE_FDATE)
c$$$      character*24 datime
c$$$      call fdate(datime)
c$$$      isdate=datime(4:11)
c$$$      istime=datime(12:19)
c$$$#else
c$$$c                1         2    
c$$$c:     :123456789012345678901234
c$$$c:ctime:Mon Aug 31 16:02:05 1992
c$$$c:fdate:Mon Nov 9 10:48:45 1992
c$$$c
      character*26 datime
      call util_date(datime)
      isdate=datime(4:11)
      istime=datime(12:19)
c$$$#endif
c$$$#endif
      if( ga_nodeid().eq.0 )then
        write (LuOut,100) name,isdate,istime
        call util_flush(LuOut)
      endif
      call ga_sync
  100 format (/,/,1x,'    <<<<<<   ',a25,'   >>>>>>',/,
     &            5x,'<<<<<< date: ',a8,' time: ',a8,' >>>>>>',/)
      return
      end
