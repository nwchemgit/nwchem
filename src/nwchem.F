      program nwchem
      implicit none
#include "mafdecls.h"
#include "global.fh"
#include "rtdb.fh"
#include "tcgmsg.fh"
      integer rtdb
      integer stack
      integer heap
      logical status
#ifdef PARALLEL_DIAG
      integer junk(4)
#endif
c
      character*255 input_filename
c
c     Create parallel processes and initialize IPC layer
c
      call pbeginf()
#ifdef PARALLEL_DIAG
      call mxpara(junk)
#endif
c
c     Only process 0 opens the input file 
c     (note that ga_nodeid() does not work yet!)
c
      if (nodeid() .eq. 0) then   
c
c     Get the name of the input file from the command line
c
         call get_input_filename(input_filename)
c
         open(5, file=input_filename, status='old', form='formatted',
     $        err=10000)
         goto 10001
10000    call errquit('nwchem: failed to open the input file', 0)
10001    continue
      else
#ifndef KSR
         close(5)
#endif
      endif
c
c     Look for memory directive in the input ... must eventually
c     open the rtdb first so that can get memory directive out of that
c     if it is not in the input
c
      call input_mem_size(stack, heap)
c
c     Initialize local memory allocator
c
      if (.not. ma_init(mt_dbl, stack, heap))
     $     call errquit('nwchem: ma_init failed', -1)
      status = ma_set_auto_verify(.true.)
c
c     Intialize global memory allocator
c
      call ga_initialize()
c
c     Read the input and open the RTDB
c
      call start(rtdb)
c
c     Call driver routine to do something useful
c
      call control(rtdb)
c
c     Close the RTDB
c
      if (.not. rtdb_close(rtdb, 'keep'))
     $     call errquit('nwchem: rtdb_close failed', rtdb)
c
c     Tidy up the global memory
c
      call ga_terminate()
c
c     Tidy up the parallel environment
c
      call pend()
c
c     Done
c
      end
      subroutine start(rtdb)
      implicit none
#include "mafdecls.h"
#include "rtdb.fh"
      integer rtdb              ! [output]
      character*255 rtdb_name
      character*255 file_prefix
      logical startup
c
c     Figure out if this is a start/restart and what the database
c     and file prefix are
c
      call input_start_opt(rtdb_name, file_prefix, startup)
c
c     Actually open the database and insert the file_prefix
c
c     Note that only process 0 has the database name ... that's OK.
c
      if (startup) then
         if (.not. rtdb_open(rtdb_name, 'empty', rtdb))
     $        call errquit('start: rtdb_open empty failed ', 0)
      else
         if (.not. rtdb_open(rtdb_name, 'old', rtdb))
     $        call errquit('start: rtdb_open old failed ', 0)
      endif
      if (.not. rtdb_cput(rtdb, 'file_prefix', 1, file_prefix))
     $     call errquit('start: rtdb_put of file_prefix failed', 0)
c
c     Parse input data and shove into the database
c
      call input_parse(rtdb)
c
c     Done
c
      end
      subroutine control(rtdb)
      implicit none
#include "rtdb.fh"
#include "inp.fh"
      integer rtdb              ! [input] data base handle
c     
c     Drive the calculation.  Eventually this routine will hook up with the
c     yet-to-be-written task module.  Right now task_input makes an rtdb
c     entry task that comprises just the name of the routine to be called.
c     If no task is present then perform some default action (exit).
c     
      character*255 task_name
      integer ntasks, ind
      logical status 
      logical scf, gradients, stpr_walk
      external scf, gradients, stpr_walk
      parameter (ntasks = 8)
      character*10 task_list(ntasks)
      data task_list /'rjhtest', 'raktest', 'rtdbprint', 'scf',
     $   'debtest', 'gradients', 'optscf','atwtest'/
c     
      if (.not. rtdb_cget(rtdb, 'task', 1, task_name)) then
c
c     Default action ... do nothing
c
         return
      endif
c     
      if (inp_match(ntasks, .false., task_name, task_list, ind)) then
         goto (10, 20, 30, 40, 50, 60, 70, 80) ind
c
         write(6,*) ' Unimplemented task ', 
     $        task_name(1:inp_strlen(task_name))
         goto 1000
c     
c     rjhtest
c     
 10      call rjhtest(rtdb)
         goto 1000
c     
c     raktest
c     
 20      call raktest(rtdb)
         goto 1000
c     
c     rtdbprint
c     
 30      if (.not. rtdb_print(rtdb, .true.))
     $        call errquit('control: rtdb_print failed', 0)
         goto 1001
c     
c     scf ... eventually there will be no such task name
c     
 40      status = scf(rtdb)
         goto 1000
C
C     debtest ... I want my own task!
C
 50      call debtest(rtdb)
         goto 1000
c
c    gradients (scf has to be done first)
c
 60      call scf_gradient(rtdb)
         goto 1000
c
c     optscf ... scf optimization
c
 70      call scf_optimize(rtdb)
         goto 1000
c
c     atwtest (4-index)
c
 80      call atwtest(rtdb)
         goto 1000
      else
c
         write(6,*) ' task_name = ', 
     $        task_name(1:inp_strlen(task_name))
         call errquit('control: unknown task name', 0)
      endif
c     
c     Eventually this will be where the next task is identified
c     
 1000 if (.not. rtdb_print(rtdb, .false.))
     $     call errquit('control: rtdb_print failed', 0)
 1001 continue
c
c     Finished ... eventually need to track status of tasks so that
c     restarts happen correctly.  For now assume that the task completed
c     and delete task from the database so that a restart will have no
c     effect.
c
      status = rtdb_delete(rtdb, 'task') ! quietly discard the status
c     
      end
