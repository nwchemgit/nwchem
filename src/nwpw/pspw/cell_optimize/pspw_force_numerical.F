*
* $Id$
*

      logical function pspw_force_numerical(rtdb)
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


*     **** external functions ****
      logical  cgsdv5,force_fd
      external cgsdv5,force_fd

      taskid = ga_nodeid()


      if (taskid.eq.MASTER) then
        write(*,*) '  ******************************************** '
        write(*,*) '  *                                          * '
        write(*,*) '  *    PSPW Force Numerical - version 1.0    * '
        write(*,*) '  *                                          * '
        write(*,*) '  *      version last updated 3/25/07        * '
        write(*,*) '  *                                          * '
        write(*,*) '  ******************************************** '
        call nwpw_message(1)
      end if


*     **** minimize psi ***
      value = cgsdv5(rtdb,1)
      if(.not. value) 
     > call errquit(
     > 'pspw_force_numerical:pspw_energy did not converge: error',
     >  0,CALC_ERR)

*     **** calculate forces ***
      value = .true.
      value = value.and.btdb_parallel(.true.)
      value = value.and.force_fd(rtdb)
      if(.not. value) 
     >   call errquit('pspw_force_numerical did not converge: error',0,
     &       CALC_ERR)
 

      call ga_sync()

      pspw_force_numerical = value
      return
      end
 

