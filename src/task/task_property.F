C> \ingroup task
C> @{
C>
C> \brief Generic interface to property calculation capabilities
C>
C> \return Return .true. if successfull, and .false. otherwise
C>
      logical function task_property(rtdb)
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"
#include "global.fh"
      integer rtdb !< [Input] The RTDB handle
c
      logical status
      logical property
      external property
      logical tce_prop_start
      external tce_prop_start
      character*32 theory
c
      if (.not. rtdb_cget(rtdb, 'task:theory', 1, theory))
     $     call errquit('task:energy: theory not specified',0, RTDB_ERR)
c
c      if (theory .eq. 'tce') then
c        status = tce_prop_start(rtdb)
c      else 
        status = property(rtdb)
c      endif
c
      if (.not. rtdb_put(rtdb, 'task:status', mt_log, 1, status))
     $     call errquit('task_property: failed to store status',0, 
     $                  RTDB_ERR)
c
      task_property = status
c
      end
C> @}
c $Id$
