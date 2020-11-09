*
* $Id$
*

*     ***********************************
*     *             			*
*     *           rho_write		*
*     *             			*
*     ***********************************

      subroutine rho_write(ispin,rho2)
      implicit none 
      integer ispin
      real*8  rho2(*)

#include "bafdecls.fh"
#include "util.fh"
#include "stdio.fh"
#include "errquit.fh"

*    *** local variables ***
      integer version,l
      integer nfft3d,n2ft3d
      integer nfft(3)
      real*8  unita(3,3)
      character*50 filename
      character*255 full_filename

      integer MASTER,taskid
      parameter(MASTER=0)
      integer n,q

c     complex*16 tmp(*)
      integer tmp(2),tmp2(2)
      logical value,lprint

*     ***** local functions ****
      character*50 control_output_rho
      external     control_output_rho
      double precision control_unita
      external         control_unita
      integer  control_ngrid,control_version
      external control_ngrid,control_version
      logical  control_print
      external control_print


      call nwpw_timing_start(50)

      call ga_sync()
      call Parallel_taskid(taskid)
      call D3dB_nfft3d(1,nfft3d)
      call D3dB_n2ft3d(2,n2ft3d)

      lprint= ((taskid.eq.MASTER).and.control_print(print_medium))

      value = BA_push_get(mt_dbl,2*nfft3d,'tmp',tmp(2),tmp(1))
      if (.not. value) 
     > call errquit('rho_write:out of stack memory',0,MA_ERR)

      version = control_version()
      nfft(1) = control_ngrid(1)
      nfft(2) = control_ngrid(2)
      nfft(3) = control_ngrid(3)

      unita(1,1) = control_unita(1,1)
      unita(2,1) = control_unita(2,1)
      unita(3,1) = control_unita(3,1)
      unita(1,2) = control_unita(1,2)
      unita(2,2) = control_unita(2,2)
      unita(3,2) = control_unita(3,2)
      unita(1,3) = control_unita(1,3)
      unita(2,3) = control_unita(2,3)
      unita(3,3) = control_unita(3,3)

*     **** open RHOIN binary file ****
      filename = control_output_rho()
      if (taskid.eq.MASTER) then
         full_filename = filename
         call util_file_name_resolve(full_filename, .false.)
         l = index(full_filename,' ') -1
         if (lprint) write(LuOut,1210) full_filename(1:l)
 1210    FORMAT('    - Writing fractional density.',
     >         /'      Filename:',A)

         call openfile(6,full_filename,l,'w',l)

         if (taskid.eq.MASTER) then
            call iwrite(6,version,1)
            call iwrite(6,nfft,3)
            call dwrite(6,unita,9)
            call iwrite(6,ispin,1)
         end if
      end if
      call D3dB_r_write(1,6,rho2,dbl_mb(tmp(1)),0)
      if (ispin.gt.1) then
         call D3dB_r_write(1,6,rho2(1+n2ft3d),dbl_mb(tmp(1)),0)
      end if

*     *** close RHOIN binary file ***
      call ga_sync()
      if (taskid.eq.MASTER) then
        call closefile(6)
      end if
      call ga_sync()

      value = BA_pop_stack(tmp(2))
      if (.not. value) 
     > call errquit('rho_write:error popping stack',0,MA_ERR)


      call nwpw_timing_end(50)
      return
      end

