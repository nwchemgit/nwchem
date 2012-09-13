C Removed task_cholesky and task_imaginary from the TASK API as they
C are not called anywhere.
C  \ingroup task
C  @{
      logical function task_cholesky(rtdb)
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
      logical tce_cholesky
      external tce_cholesky
      character*32 theory
c
      if (.not. rtdb_cget(rtdb, 'task:theory', 1, theory))
     $ call errquit('task_cholesky: theory not specified',0, RTDB_ERR)
c
      if (theory .eq. 'tce') then
        status = tce_cholesky(rtdb)
      else 
       call errquit('task_cholesky: only valid for tce',0,RTDB_ERR)
      endif
c
      if (.not. rtdb_put(rtdb, 'task:status', mt_log, 1, status))
     $     call errquit('task_cholesky: failed to store status',0, 
     $                  RTDB_ERR)
c
      task_cholesky = status
c
      end
C  
      logical function task_imaginary(rtdb)
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
      logical tce_imaginary
      external tce_imaginary
      character*32 theory
c
      if (.not. rtdb_cget(rtdb, 'task:theory', 1, theory))
     $ call errquit('task_imaginary: theory not specified',0, RTDB_ERR)
c
      if (theory .eq. 'tce') then
        status = tce_imaginary(rtdb)
      else 
       call errquit('task_imaginary: only valid for tce',0,RTDB_ERR)
      endif
c
      if (.not. rtdb_put(rtdb, 'task:status', mt_log, 1, status))
     $     call errquit('task_imaginary: failed to store status',0, 
     $                  RTDB_ERR)
c
      task_imaginary = status
c
      end
C  @}
c $Id$
