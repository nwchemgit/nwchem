*
* $Id: pspw_stress_numerical.F,v 1.3 2003-10-17 22:57:08 carlfahl Exp $
*

      logical function pspw_stress_numerical(rtdb)
      implicit none
#include "errquit.fh"
      integer rtdb

#include "mafdecls.fh"
#include "global.fh"
#include "msgtypesf.h"
#include "rtdb.fh"
#include "geom.fh"

*     **** local variables ***
      integer taskid,MASTER
      parameter (MASTER=0)
      logical value,value2

      real*8 stress(3,3)
      real*8 unita(3,3)

*     **** external functions ****
      logical  pspw_stress,stress_fd
      external pspw_stress,stress_fd

      taskid = ga_nodeid()

*     **** calculate stress tensor ***
      value = pspw_stress(rtdb)
      if(.not. value) 
     >   call errquit('pspw_stress did not converge: error',0, CALC_ERR)


      if (taskid.eq.MASTER) then
        write(*,*) '  ******************************************** '
        write(*,*) '  *                                          * '
        write(*,*) '  *    PSPW Stress Numerical - version 1.0   * '
        write(*,*) '  *                                          * '
        write(*,*) '  *      version last updated 4/24/02        * '
        write(*,*) '  *                                          * '
        write(*,*) '  ******************************************** '
        call nwpw_message(1)
      end if

*     **** read in initial data ***
      value = .true.
      value2 = rtdb_parallel(.true.)
 
      value = stress_fd(rtdb)
      if(.not. value) 
     >   call errquit('pspw_stress_numerical did not converge: error',0,
     &       CALC_ERR)
 

      call ga_sync()

      pspw_stress_numerical = value
      return
      end
 

