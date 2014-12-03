c
c     $Id$                       
c

      subroutine paw_sdminimize(flag)
      implicit none
      integer flag

#include "bafdecls.fh"
#include "stdio.fh"
#include "util.fh"

*     **** local variables ****
      integer MASTER,taskid
      parameter (MASTER=0)

      logical oprint
      integer it,it_in
      real*8 dt,dte

*     **** external functions ****
      logical  control_print
      integer  control_it_in
      real*8   control_time_step,control_fake_mass

      external control_print
      external control_it_in
      external control_time_step,control_fake_mass

      dt = control_time_step()
      dte = dt/dsqrt(control_fake_mass())
      it_in = control_it_in()
      if (flag.gt.0) it_in = flag

      do it=1,it_in
         call paw_psi1to2_sd_update(dte)
         call paw_psi_2to1()
      end do

      call Parallel_taskid(taskid)
      oprint = (taskid.eq.MASTER).and.control_print(print_medium)
      if (oprint) then
        write(luout,'(A,I3,A)') "     - ",
     >             it_in," steepest descent iterations performed"
      end if
  
      return
      end
 

