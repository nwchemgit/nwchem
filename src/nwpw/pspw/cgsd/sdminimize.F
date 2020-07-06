c
c     $Id$                       
c

      subroutine sdminimize(flag)
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

      call nwpw_cosmo_reset_Qt()
      dt = control_time_step()
      dte = dt/dsqrt(control_fake_mass())
      it_in = control_it_in()
      if (flag.gt.0) it_in = flag

      do it=1,it_in
         call psi1to2_sd_update(dte)
         call psi_2to1()
      end do

      call Parallel_taskid(taskid)
      oprint = (taskid.eq.MASTER).and.control_print(print_medium)
      if (oprint) then
        write(luout,'(A,I3,A)') "     - ",
     >             it_in," steepest descent iterations performed"
      end if
  
      return
      end
 

      subroutine sdminimize_virtual(flag)
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

      write(*,*) "Warning - sdminimize_virutal not finished"
      dt = control_time_step()
      dte = dt/dsqrt(control_fake_mass())
      it_in = control_it_in()
      if (flag.gt.0) it_in = flag

      do it=1,it_in
c         call epsi1to2_sd_update_virtual(dte)
         call epsi_2to1()
      end do

      call Parallel_taskid(taskid)
      oprint = (taskid.eq.MASTER).and.control_print(print_medium)
      if (oprint) then
        write(luout,'(A,I3,A)') "     - ",
     >             it_in," steepest descent iterations performed"
      end if

      return
      end


