      logical function task_property(rtdb)
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"
#include "global.fh"
      integer rtdb
c
      logical status
      logical property
      external property
c
      status = property(rtdb)
c
      if (.not. rtdb_put(rtdb, 'task:status', mt_log, 1, status))
     $     call errquit('task_property: failed to store status',0, 
     $                  RTDB_ERR)
c
      task_property = status
c
      end
