      logical function cell_optimize(rtdb)
      implicit none
      integer rtdb

#include "mafdecls.fh"
#include "global.fh"
#include "msgtypesf.h"
#include "rtdb.fh"

*     **** local variables ***
      integer taskid,MASTER
      parameter (MASTER=0)
      logical value,value2

      character*30     cell_name
      integer          optimize_strategy,cycles
      integer          optimize_lattice(6)
      double precision lattice_tolerances(6)

      integer i

*     **** external functions ****
      logical  psp_stress_formatter
      logical  pspw_stress,psp_formatter,wvfnc_reformat,update_unitcell
      external psp_stress_formatter
      external pspw_stress,psp_formatter,wvfnc_reformat,update_unitcell

      taskid = ga_nodeid()

*     **** read in initial data ***
      value = rtdb_cget(rtdb,'cell_optimize:cell_name',1,cell_name)

      value = value.and.
     >        rtdb_get(rtdb,'cell_optimize:optimize_strategy',
     >                 mt_int,1,optimize_strategy)

      value = value.and.
     >        rtdb_get(rtdb,'cell_optimize:cycles',
     >                 mt_int,1,cycles)

      value = value.and.
     >        rtdb_get(rtdb,'cell_optimize:optimize_lattice',
     >                 mt_int,6,optimize_lattice)

      value = value.and.
     >        rtdb_get(rtdb,'cell_optimize:lattice_tolerances',
     >                 mt_dbl,6,lattice_tolerances)


*     *** exit if error found ***
      if(.not. value) 
     >    call errquit('cell_optimize: error reading rtdb',0)


      if (taskid.eq.MASTER) then
        write(*,*) '  ******************************************** '
        write(*,*) '  *                                          * '
        write(*,*) '  *    PSPW Cell Optimize - version 0.1      * '
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
      value2 = rtdb_parallel(.true.)


      call ga_sync()

*      **** reformat pseudopotentials ****
       value = psp_formatter(rtdb)
       if(.not. value) call errquit(
     >    'cell_optimize: error reformatting pseudopotentials',0)

       value = psp_stress_formatter(rtdb)
       if(.not. value) call errquit(
     >   'cell_optimize: error reformatting stress pseudopotentials',0)

*     **** start optimization loop ***
      do i=1,cycles

*         **** calculate stress tensor ***
          value = pspw_stress(rtdb)
          if(.not. value) 
     >      call errquit('pspw_stress did not converge: error',0)


*         **** reformat/adjust unit cell ****
          value2 = rtdb_parallel(.false.)
          if (taskid.eq.MASTER) value = update_unitcell(rtdb)
          value2 = rtdb_parallel(.true.)
          if(.not.value) call errquit(
     >       'cell_optimize: error reformatting unit cell',0)


*         **** reformat pseudopotentials ****
          value = psp_formatter(rtdb)
          if(.not. value) call errquit(
     >       'cell_optimize: error reformatting pseudopotentials',0)

          value = psp_stress_formatter(rtdb)
          if(.not. value) call errquit(
     >  'cell_optimize: error reformatting stress pseudopotentials',0)

*         **** reformat wavefunctions ****
          value2 = rtdb_parallel(.false.)
          if (taskid.eq.MASTER) value = wvfnc_reformat(rtdb)
          value2 = rtdb_parallel(.true.)
          if(.not. value)  call errquit(
     >       'cell_optimize: error reformatting wavefunctions',0)

          call ga_sync()
      end do
      value = .false.
      goto 999

  998 value = .true.

*     **** set to true for debugging ****
  999 value = .true.
      cell_optimize = value
      return
      end
 

