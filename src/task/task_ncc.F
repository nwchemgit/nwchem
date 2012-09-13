C Taken task_ncc out of the TASK module as it seems not to be used
C anywhere.
C  \ingroup task
C  @{
      logical function task_ncc(rtdb)
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"
#include "global.fh"
      integer rtdb
c
      logical stat
      logical ncc_driver
      external ncc_driver
      character*32 theory
c
c      print*,'top of task_ncc'
c
      task_ncc = .false.
c
      if (.not. rtdb_cget(rtdb, 'task:theory', 1, theory))
     $ call errquit('task_ncc: theory not specified',0, RTDB_ERR)
c
      if (theory .eq. 'ncc') then
        if (.not. ncc_driver(rtdb) ) then
          call errquit('task_ncc: ncc_driver failed',0,RTDB_ERR)
        else
          task_ncc = .true.
        endif
      else 
        call errquit('task_ncc: only valid for ncc',0,RTDB_ERR)
      endif
c
      if (.not. rtdb_put(rtdb, 'task:status', mt_log, 1, stat)) then
        call errquit('task_ncc: failed to store stat',0,RTDB_ERR)
      endif
c
c      print*,'end of task_ncc'
c
      end
C  @}
c $Id$
