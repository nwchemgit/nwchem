c     $Id: sdminimize.F,v 1.4 2002-01-09 19:02:18 bylaska Exp $                       

      subroutine sdminimize(flag)
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
         call psi1to2_sd_update(dte)
         call psi_2to1()
      end do

      call Parallel_taskid(taskid)
      if (taskid.eq.MASTER) then
        write(*,'(A,I3,A)') "     - ",
     >             it_in," steepest descent iterations performed"
      end if
  
      return
      end
 

