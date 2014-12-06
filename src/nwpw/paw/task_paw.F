*
* $Id$
*
C> \ingroup task
C> @{
      logical function task_paw(rtdb)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "global.fh"
#include "bafdecls.fh"
#include "btdb.fh"
c
      integer rtdb
      logical value,dummy
      integer MASTER
      parameter (MASTER=0)
c
      integer ind               ! Index of matched directive
      integer num_dirs          ! No. of known directives
      parameter (num_dirs = 4)
      character*30 dirs(num_dirs)
      character*255 test
      data dirs / 
     >            'steepest_descent',
     >            'Car-Parrinello',
     >            'noit_energy',
     >            'stress' /

c     *** define external functions ****
      logical  paw_sd,paw_cpmd,paw_noit_energy,paw_stress
      external paw_sd,paw_cpmd,paw_noit_energy,paw_stress
c
c
      value = .false.
      task_paw = .false. 
 10   value = btdb_cget(rtdb,'task:operation',1,test)

      if (.not. inp_match(num_dirs, .false., test, dirs, ind))
     >     call errquit('task paw: unknown directive', 0, INPUT_ERR)
c
c
c
      goto ( 100,  200, 300, 400) ind
      call errquit('task paw: unimplemented directive', ind, INPUT_ERR)


c
c  steepest_descent
c
 100  if (ga_nodeid().eq.MASTER) then
        write(*,*)  ">>>> PAW Parallel Module - steepest_descent <<<<"
      end if

      value = .true.
      value = btdb_parallel(.true.)
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
      value = btdb_parallel(.true.)
      value = paw_cpmd(rtdb)
      call ga_sync()
      goto 9999

c
c  noit_energy
c
 300  if (ga_nodeid().eq.MASTER) then
        write(*,*)  ">>>> PAW Parallel Module - noit_energy <<<<"
      end if
      value = .true.
      value = btdb_parallel(.true.)
      value = paw_noit_energy(rtdb)
      call ga_sync()
      goto 9999

c
c  stress
c
 400  if (ga_nodeid().eq.MASTER) then
        write(*,*)  ">>>> PAW Parallel Module - stress <<<<"
      end if
      value = .true.
      value = btdb_parallel(.true.)
      value = paw_stress(rtdb)
      call ga_sync()
      goto 9999




 9999 continue
      dummy = btdb_parallel(.true.)
      task_paw = value
      

      return
      end
C> @}
