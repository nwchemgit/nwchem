c     $Id: c_sdminimize.F,v 1.1 2001-12-17 23:31:13 bylaska Exp $                       

      subroutine c_sdminimize(flag)
      implicit none
      integer flag

#include "mafdecls.fh"

*     **** local variables ****
      integer MASTER,taskid
      parameter (MASTER=0)

      integer it,it_in
      real*8 dt,dte

*     **** external functions ****
      integer  control_it_in
      real*8   control_time_step,control_fake_mass

      external control_it_in
      external control_time_step,control_fake_mass

      dt = control_time_step()
      dte = dt/dsqrt(control_fake_mass())
      it_in = control_it_in()
      if (flag.gt.0) it_in = flag

      do it=1,it_in
         call cpsi1to2_sd_update(dte)
         call cpsi_2to1()
      end do

      call Parallel_taskid(taskid)
      if (taskid.eq.MASTER) then
        write(*,*) "     - ",
     >             it_in," steepest descent iterations performed"
      end if
  
      return
      end
 

