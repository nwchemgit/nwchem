*
* $Id: task_paw.F,v 1.5 2002-08-01 01:18:27 bylaska Exp $
*

      logical function task_paw(rtdb)
      implicit none
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
      logical  paw_sd
      external paw_sd
c
c
      value = .false.
      task_paw = .false. 
 10   value = rtdb_cget(rtdb,'task:operation',1,test)

      if (.not. inp_match(num_dirs, .false., test, dirs, ind))
     >     call errquit('task paw: unknown directive', 0)
c
c
c
      goto ( 100,  200) ind
      call errquit('task paw: unimplemented directive', ind)


c
c  steepest_descent
c
 100  if (ga_nodeid().eq.MASTER) then
!        write(*,*) "generating PAW basis"
!       call init_paw_basis(2,(/"h","o"/))
!       call init_paw_kin_matrix()
!       call init_paw_ion_matrix()
!       call init_paw_vloc_matrix()
!       call init_paw_core_matrix()

!      call dealloc_paw_basis_data()
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
        write(*,*) "not finished!"
      end if
      call ga_sync()
      goto 9999


 9999 continue
      dummy = rtdb_parallel(.true.)
      task_paw = value
      

      return
      end
