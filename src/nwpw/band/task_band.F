*
* $Id: task_band.F,v 1.3 2004-11-29 16:05:16 bylaska Exp $
*

      logical function task_band(rtdb)
      implicit none

#include "inp.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "errquit.fh"
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
      data dirs / 'structure',
     >            'noit_energy',
     >            'stress',
     >            'cell_optimize' /

c     *** define external functions ****
      logical  band_structure,band_noit_energy,band_stress
      logical  band_cell_optimize
      external band_structure,band_noit_energy,band_stress
      external band_cell_optimize
c
c
      value = .false.
      task_band = .false. 
 10   value = rtdb_cget(rtdb,'task:operation',1,test)

      if (.not. inp_match(num_dirs, .false., test, dirs, ind))
     >     call errquit('task band: unknown directive', 0, INPUT_ERR)
c
c
c
      goto ( 100, 200, 300, 400 ) ind
      call errquit('task band: unimplemented directive', ind, INPUT_ERR)


c
c  band_structure
c
 100  value = .true.
      
      if (ga_nodeid().eq.MASTER) then 
        write(*,*)  ">>>> BAND parallel Module - band_structure <<<<"  
      end if
      value = rtdb_parallel(.true.)
      value = band_structure(rtdb)
      call ga_sync()
      goto 9999
c
c  noit_energy - debug
c
 200  if (ga_nodeid().eq.MASTER) then
        write(*,*)  ">>>> BAND Parallel Module - noit_energy <<<<"
      end if
      value = .true.
      value = rtdb_parallel(.true.)
      value = band_noit_energy(rtdb)
      call ga_sync()
      goto 9999

c
c  stress
c
 300  if (ga_nodeid().eq.MASTER) then
        write(*,*)  ">>>> BAND Parallel Module - stress <<<<"
      end if
      value = .true.
      value = rtdb_parallel(.true.)
      value = band_stress(rtdb)
      call ga_sync()
      goto 9999
c
c  cell_optimize
c
 400  if (ga_nodeid().eq.MASTER) then
        write(*,*)  ">>>> BAND Parallel Module - cell_optimize <<<<"
      end if
      value = .true.
      value = rtdb_parallel(.true.)
      value = band_cell_optimize(rtdb)
      call ga_sync()
      goto 9999





 9999 continue
      dummy = rtdb_parallel(.true.)
      task_band = value

      return
      end
