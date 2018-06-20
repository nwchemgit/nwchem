*     *****************************************
*     *                                       *
*     *             vdw_DF_init               *
*     *                                       *
*     *****************************************
*
*
      subroutine vdw_DF_init()
      implicit none

#include "inp.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "util.fh"
#include "bafdecls.fh"

#include "vdw-DF.fh"

*     **** local variables ****
      integer taskid,MASTER
      parameter (MASTER=0)
      logical mprint,hprint,debug,does_it_exist
      logical from_environment,from_compile,from_nwchemrc
      integer iop,lgth,unitf,print_level,i,j
      character*255 datafile
      logical found_datafile
      integer ifound

*     **** external functions ***
      logical  control_has_vdw,control_is_vdw2
      external control_has_vdw,control_is_vdw2


*     **** return if not has_vdw ****
      has_vdw = control_has_vdw()
      is_vdw2 = control_is_vdw2()
      if (.not.has_vdw) return

      write(*,*) "I have vdw!"
      call Parallel_taskid(taskid)

      call util_file_name_noprefix('vdw_kernels.dat',
     >                                .false.,
     >                                .false.,
     >                                datafile)
      if (taskid.eq.MASTER) then
         ifound = 0
         inquire(file=datafile,exist=found_datafile)
         if (found_datafile) ifound = 1
      end if
      call Parallel_Brdcst_ivalue(MASTER,ifound)
      found_datafile = .false.
      if (ifound.eq.1) found_datafile=.true.

      !**** generate vdw_kernels.dat file ****
      if (.not.found_datafile) then
         write(*,*) "-- generate vdw_kernels.dat --"
         call vdw_DF_kernel_gen_data(datafile)
      end if


      write(*,*) "-- reading vdw_kernels.dat --"
      return

 999  continue
      call errquit('vdw_DF_init:error reading qmesh,Nqs=',Nqs,DISK_ERR)
      return

      end


*     **********************************************
*     *                                            *
*     *                vdw_DF_end                  *
*     *                                            *
*     **********************************************

      subroutine vdw_DF_end()
      implicit none

#include "vdw-DF.fh"

*     **** local variables ****
      logical value

      if (has_vdw) then
         !value = BA_free_heap(ufunc(2))
         !value = value.BA_free_heap(qmesh(2))
         !value = value.and.BA_free_heap(ya(2))
         !value = value.and.BA_free_heap(y2a(2))
c         if (.not.value) 
c     >      call errquit('vdw_DF_end:free heap failed',0,MA_ERR)
      end if

      return
      end

*     **********************************************
*     *                                            *
*     *                vdw_DF_init_poly            *
*     *                                            *
*     **********************************************
      subroutine vdw_DF_init_poly(n,x,ya,y2a,utmp)
      implicit none
      integer n
      real*8 x(*),ya(n,*),y2a(n,*),utmp(*)

*     **** local variables ****
      integer i,j
      real*8 yp1,ypn

      yp1 = 0.0d0
      ypn = 0.0d0

      do i=1,n
         do j=1,n
            if (i.eq.j) then
               ya(i,j) = 1.0d0
            else
               ya(i,j) = 0.0d0
            end if
         end do
      end do
      do j=1,n
         call nwpw_spline(x,ya(1,j),n,yp1,ypn,y2a(1,j),utmp)
      end do

      return
      end

