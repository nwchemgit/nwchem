*
* $Id$
*
C> \ingroup task
C> @{
      logical function task_ofpw(rtdb)
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
      parameter (num_dirs = 5)
      character*30 dirs(num_dirs)
      character*255 test
      data dirs / 'steepest_descent',
     >            'Car-Parrinello',
     >            'stress',
     >            'noit_energy',
     >            'Born-Oppenheimer'/

c     *** define external functions ****
      logical  ofpwsd
      external ofpwsd
c
c
      value = .false.
      task_ofpw = .false. 
 10   value = btdb_cget(rtdb,'task:operation',1,test)

      if (.not. inp_match(num_dirs, .false., test, dirs, ind))
     >     call errquit('task ofpw: unknown directive', 0, INPUT_ERR)
c
c
c
      goto ( 100,  200, 300, 400, 500 ) ind
      call errquit('task ofpw: unimplemented directive', ind, INPUT_ERR)

c
c  steepest_descent
c
 100  if (ga_nodeid().eq.MASTER) then
        write(*,*)  ">>>> OFPW Parallel Module - steepest_descent <<<<"
      end if

      value = .true.
      value = btdb_parallel(.true.)
      value = ofpwsd(rtdb)
      call ga_sync()
      goto 9999


c
c  Car-Parrinello
c
 200  if (ga_nodeid().eq.MASTER) then
        write(*,*)  ">>>> OFPW Parallel Module - Car-Parrinello <<<<"
      end if
      value = .true.
      value = btdb_parallel(.true.)
      !value = ofcpmd(rtdb)
      call ga_sync()
      goto 9999


c
c  stress - debug
c
 300  if (ga_nodeid().eq.MASTER) then
        write(*,*)  ">>>> OFPW Parallel Module - stress <<<<"
      end if
      value = .true.
      value = btdb_parallel(.true.)
      !value = ofpw_stress(rtdb)
      call ga_sync()
      goto 9999

c
c  noit_energy - debug
c
 400  if (ga_nodeid().eq.MASTER) then
        write(*,*)  ">>>> OFPW Parallel Module - noit_energy <<<<"
      end if
      value = .true.
      value = btdb_parallel(.true.)
      !value = ofpw_noit_energy(rtdb)
      call ga_sync()
      goto 9999

c
c  Born-Oppenheimer
c
 500  if (ga_nodeid().eq.MASTER) then
        write(*,*)  ">>>> OFPW Parallel Module - Born-Oppenheimer <<<<"
      end if
      value = .true.
      value = btdb_parallel(.true.)
      !value = ofcgmd(rtdb,0)
      call ga_sync()
      goto 9999


 9999 continue
      dummy = btdb_parallel(.true.)
      task_ofpw = value

      return
      end
C> @}
