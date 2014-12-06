      logical function band_cell_optimize(rtdb)
      implicit none
#include "errquit.fh"
      integer rtdb

#include "bafdecls.fh"
#include "global.fh"
#include "msgtypesf.h"
#include "btdb.fh"

*     **** local variables ***
      integer taskid,MASTER
      parameter (MASTER=0)
      logical value,value2

      character*50     cell_name
      integer          optimize_strategy,cycles
      integer          optimize_lattice(6)
      double precision lattice_tolerances(6)

      integer i

*     **** external functions ****
      logical  band_stress,band_update_unitcell
      external band_stress,band_update_unitcell

      taskid = ga_nodeid()

*     **** read in initial data ***
      value = btdb_cget(rtdb,'cell_optimize:cell_name',1,cell_name)

      value = value.and.
     >        btdb_get(rtdb,'cell_optimize:optimize_strategy',
     >                 mt_int,1,optimize_strategy)

      value = value.and.
     >        btdb_get(rtdb,'cell_optimize:cycles',
     >                 mt_int,1,cycles)

      value = value.and.
     >        btdb_get(rtdb,'cell_optimize:optimize_lattice',
     >                 mt_int,6,optimize_lattice)

      value = value.and.
     >        btdb_get(rtdb,'cell_optimize:lattice_tolerances',
     >                 mt_dbl,6,lattice_tolerances)


*     *** exit if error found ***
      if(.not. value) 
     >    call errquit('cell_optimize: error reading rtdb',0, RTDB_ERR)


      if (taskid.eq.MASTER) then
        write(*,*) '  ******************************************** '
        write(*,*) '  *                                          * '
        write(*,*) '  *    BAND Cell Optimize - version 0.1      * '
        write(*,*) '  *                                          * '
        write(*,*) '  *      version last updated 5/23/01        * '
        write(*,*) '  *                                          * '
        write(*,*) '  ******************************************** '
        call nwpw_message(1)

        write(*,*) 
        write(*,*) ' Optimizing cell_name       :',cell_name
        write(*,*) ' Optimizing strategy        :',optimize_strategy
        write(*,*) ' optimize_lattice           :',optimize_lattice
        write(*,*) ' Tolerances                 :',lattice_tolerances
        write(*,*) ' Number optimization cycles :',cycles
        write(*,*) 
      endif

*     **** read in initial data ***
      value = .true.
      value2 = btdb_parallel(.true.)


      call ga_sync()


*     **** start optimization loop ***
      do i=1,cycles

*         **** calculate stress tensor ***
          value = band_stress(rtdb)
          if(.not. value) 
     >      call errquit('pspw_stress did not converge: error',0,
     &       CALC_ERR)


*         **** reformat/adjust unit cell ****
          value2 = btdb_parallel(.false.)
          if (taskid.eq.MASTER) value = band_update_unitcell(rtdb)
          value2 = btdb_parallel(.true.)
          if(.not.value) call errquit(
     >       'cell_optimize: error reformatting unit cell',0,
     &       CALC_ERR)

          call ga_sync()
      end do
      value = .false.
      goto 999

  998 value = .true.

*     **** set to true for debugging ****
  999 value = .true.
      band_cell_optimize = value
      return
      end
 

c $Id$
