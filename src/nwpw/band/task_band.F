*
* $Id$
*
C> \ingroup task
C> @{
      logical function task_band(rtdb)
      implicit none

#include "inp.fh"
#include "global.fh"
#include "bafdecls.fh"
#include "btdb.fh"
#include "errquit.fh"
c
      integer rtdb
      logical value,dummy
      integer MASTER
      parameter (MASTER=0)
c
      integer ind               ! Index of matched directive
      integer num_dirs          ! No. of known directives
      parameter (num_dirs = 9)
      character*30 dirs(num_dirs)
      character*255 test
      data dirs / 'structure',
     >            'dos',
     >            'noit_energy',
     >            'stress',
     >            'cell_optimize',
     >            'band_dplot',
     >            'dos_dplot',
     >            'steepest_descent',
     >            'car-parrinello' /

c     *** define external functions ****
      logical  band_structure,band_noit_energy,band_stress
      logical  band_cell_optimize,band_dplot,band_sd,band_cpmd
      external band_structure,band_noit_energy,band_stress
      external band_cell_optimize,band_dplot,band_sd,band_cpmd
c
c
      value = .false.
      task_band = .false. 
 10   value = btdb_cget(rtdb,'task:operation',1,test)

      if (.not. inp_match(num_dirs, .false., test, dirs, ind))
     >     call errquit('task band: unknown directive', 0, INPUT_ERR)
c
c
c
      goto ( 100, 110, 200, 300, 400, 500, 600, 700, 800 ) ind
      call errquit('task band: unimplemented directive', ind, INPUT_ERR)


c
c  band_structure
c
 100  value = .true.
      
      if (ga_nodeid().eq.MASTER) then 
        write(*,*)  ">>>> BAND parallel Module - band_structure <<<<"  
      end if
      value = btdb_parallel(.true.)
      value = band_structure(rtdb,0)
      call ga_sync()
      goto 9999

c
c  DOS calculation
c
 110  value = .true.

      if (ga_nodeid().eq.MASTER) then
        write(*,*)  ">>>> BAND parallel Module - dos <<<<"
      end if
      value = btdb_parallel(.true.)
      value = band_structure(rtdb,1)
      call ga_sync()
      goto 9999

c
c  noit_energy - debug
c
 200  if (ga_nodeid().eq.MASTER) then
        write(*,*)  ">>>> BAND Parallel Module - noit_energy <<<<"
      end if
      value = .true.
      value = btdb_parallel(.true.)
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
      value = btdb_parallel(.true.)
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
      value = btdb_parallel(.true.)
      value = band_cell_optimize(rtdb)
      call ga_sync()
      goto 9999

c
c  band_dplot
c
 500  if (ga_nodeid().eq.MASTER) then
        write(*,*)  ">>>> BAND Parallel Module - band_dplot <<<<"
      end if
      value = .true.
      value = btdb_parallel(.true.)
      value = band_dplot(rtdb)
      call ga_sync()
      goto 9999

c
c  DOS_dplot calculation
c
 600  value = .true.

      if (ga_nodeid().eq.MASTER) then
        write(*,*)  ">>>> BAND parallel Module - dos_dplot <<<<"
      end if
      value = btdb_parallel(.true.)
      value = band_structure(rtdb,2)
      call ga_sync()
      goto 9999

c
c  steepest_descent
c
 700  if (ga_nodeid().eq.MASTER) then
        write(*,*)  ">>>> BAND Parallel Module - steepest_descent <<<<"
      end if

      value = .true.
      value = btdb_parallel(.true.)

      value = band_sd(rtdb)

      call ga_sync()
      goto 9999


c
c  Car-Parrinello
c
 800  if (ga_nodeid().eq.MASTER) then
        write(*,*)  ">>>> BAND Parallel Module - Car-Parrinello <<<<"
      end if

      value = .true.
      value = btdb_parallel(.true.)

      value = band_cpmd(rtdb)

      call ga_sync()
      goto 9999




 9999 continue
      dummy = btdb_parallel(.true.)
      task_band = value

      return
      end
C> @}
