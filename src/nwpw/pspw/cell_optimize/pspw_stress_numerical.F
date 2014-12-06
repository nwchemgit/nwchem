*
* $Id$
*

      logical function pspw_stress_numerical(rtdb)
      implicit none
#include "errquit.fh"
      integer rtdb

#include "bafdecls.fh"
#include "global.fh"
#include "msgtypesf.h"
#include "btdb.fh"
#include "geom.fh"

*     **** local variables ***
      integer taskid,MASTER
      parameter (MASTER=0)
      logical value,value2

      real*8 stress(3,3)
      real*8 unita(3,3)

*     **** external functions ****
      logical  cgsdv5,stress_fd
      external cgsdv5,stress_fd

      taskid = ga_nodeid()


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


*     **** minimize psi ***
      value = cgsdv5(rtdb,1)
      if(.not. value) 
     > call errquit(
     > 'pspw_stress_numerical:pspw_energy did not converge: error',
     >  0,CALC_ERR)

*     **** calculate stress tensor ***
      value = .true.
      value = value.and.btdb_parallel(.true.)
      value = value.and.stress_fd(rtdb)
      if(.not. value) 
     >   call errquit('pspw_stress_numerical did not converge: error',0,
     &       CALC_ERR)
 

      call ga_sync()

      pspw_stress_numerical = value
      return
      end
 

