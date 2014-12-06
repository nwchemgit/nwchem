*
* $Id$
*

      logical function paw_stress_numerical(rtdb)
      implicit none
      integer rtdb

#include "bafdecls.fh"
#include "global.fh"
#include "msgtypesf.h"
#include "btdb.fh"
#include "geom.fh"
#include "errquit.fh"

*     **** local variables ***
      integer taskid,MASTER
      parameter (MASTER=0)
      logical value,value2

      real*8 stress(3,3)
      real*8 unita(3,3)

*     **** external functions ****
      logical  paw_minimizer,paw_stress_fd
      external paw_minimizer,paw_stress_fd

      taskid = ga_nodeid()


      if (taskid.eq.MASTER) then
        write(*,*) '  ******************************************** '
        write(*,*) '  *                                          * '
        write(*,*) '  *     PAW Stress Numerical - version 1.0   * '
        write(*,*) '  *                                          * '
        write(*,*) '  *      version last updated 4/24/02        * '
        write(*,*) '  *                                          * '
        write(*,*) '  ******************************************** '
        call nwpw_message(1)
      end if


*     **** minimize psi ***
      value = paw_minimizer(rtdb,1)
      if(.not. value) 
     > call errquit(
     > 'paw_stress_numerical:paw_gradient did not converge: error',
     >  0,CALC_ERR)

*     **** calculate stress tensor ***
      value = .true.
      value = value.and.btdb_parallel(.true.)
      value = value.and.paw_stress_fd(rtdb)
      if(.not. value) 
     >   call errquit('paw_stress_numerical did not converge: error',0,
     &       CALC_ERR)
 

      call ga_sync()

      paw_stress_numerical = value
      return
      end
 

