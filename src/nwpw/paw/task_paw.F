*
* $Id: task_paw.F,v 1.8 2003-10-17 22:56:54 carlfahl Exp $
*

      logical function task_paw(rtdb)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
c
      integer rtdb
      logical value,dummy
      integer MASTER
      parameter (MASTER=0)
c
      integer ind               ! Index of matched directive
      integer num_dirs          ! No. of known directives
      parameter (num_dirs = 2)
      character*30 dirs(num_dirs)
      character*255 test
      data dirs / 
     >            'steepest_descent',
     >            'Car-Parrinello' /

c     *** define external functions ****
      logical  paw_sd,paw_cpmd
      external paw_sd,paw_cpmd
c
c
      value = .false.
      task_paw = .false. 
 10   value = rtdb_cget(rtdb,'task:operation',1,test)

      if (.not. inp_match(num_dirs, .false., test, dirs, ind))
     >     call errquit('task paw: unknown directive', 0, INPUT_ERR)
c
c
c
      goto ( 100,  200) ind
      call errquit('task paw: unimplemented directive', ind, INPUT_ERR)


c
c  steepest_descent
c
 100  if (ga_nodeid().eq.MASTER) then
        write(*,*)  ">>>> PAW Parallel Module - steepest_descent <<<<"
      end if

      value = .true.
      value = rtdb_parallel(.true.)
      value = paw_sd(rtdb)
      call ga_sync()
      goto 9999

c
c  Car-Parrinello
c
 200  if (ga_nodeid().eq.MASTER) then
        write(*,*)  ">>>> PAW Parallel Module - Car-Parrinello <<<<"
      end if
      value = .true.
      value = rtdb_parallel(.true.)
      value = paw_cpmd(rtdb)
      call ga_sync()
      goto 9999


 9999 continue
      dummy = rtdb_parallel(.true.)
      task_paw = value
      

      return
      end
