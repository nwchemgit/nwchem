      subroutine util_rtdb_speak(rtdb)
*
* $Id$
*
      implicit none
#include "rtdb.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "inp.fh"
#include "util.fh"
      integer rtdb
c
      character*255 hostname
      integer port
c
      if (.not. rtdb_cget(rtdb,'speech:hostname',1,hostname))
     $     return
      if (.not. rtdb_get(rtdb,'speech:port',mt_int,1,port))
     $     return
c
      if (ga_nodeid() .gt. 0) return
c
      if (util_print('speech',print_low)) then
         write(6,1) hostname(1:inp_strlen(hostname)), port
 1       format(/'  Speech output at ',a,' port number',i6/)
         call util_flush(6)
      endif
c
      call util_speak_init(hostname, port)
c
      end

