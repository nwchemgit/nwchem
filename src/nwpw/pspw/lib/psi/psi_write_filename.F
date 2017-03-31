*
* $Id$
*

*     ***********************************
*     *             			*
*     *           psi_write_filename	*
*     *             			*
*     ***********************************

      subroutine psi_write_filename(filename,ispin,ne,psi2)
      implicit none 
      character*(*) filename
      integer ispin,ne(2)
      double complex psi2(*)
      integer occupation

#include "bafdecls.fh"
#include "errquit.fh"

*    *** local variables ***
      integer version,l
      integer nfft3d,npack1
      integer nfft(3)
      real*8  unita(3,3)
      character*255 full_filename

      integer MASTER,taskid,taskid_i,taskid_j,taskid_p,com_p
      parameter(MASTER=0)
      integer n,q,pj

c     complex*16 tmp(*)
      integer tmp(2),tmp2(2)
      logical value,pio,doflush

*     ***** local functions ****
      double precision control_unita
      external         control_unita
      integer  control_ngrid,control_version
      external control_ngrid,control_version
      logical  control_parallel_io
      external control_parallel_io

      call nwpw_timing_start(50)
      call Parallel_taskid(taskid)
      call Parallel2d_taskid_i(taskid_i)
      call Parallel2d_taskid_j(taskid_j)
      call D3dB_nfft3d(1,nfft3d)
      call Pack_npack(1,npack1)

      doflush = .false.
      pio = control_parallel_io()
      if (pio) then
         taskid_p = taskid_i
         com_p = 1
      else
         taskid_p = taskid
         com_p = 0
      end if

      value = BA_push_get(mt_dcpl,nfft3d,'tmp',tmp(2),tmp(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,nfft3d,'tmp2',tmp2(2),tmp2(1))
      if (.not. value) 
     > call errquit('psi_write_filename:out of stack',0,MA_ERR)

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

*     **** open ELCIN binary file ****
      if (taskid_p.eq.MASTER) then
         call util_file_name_noprefix(filename,.false.,
     >                                .false.,
     >                       full_filename)
         l = index(full_filename,' ') -1
         doflush = .true.
         call openfile(6,full_filename,l,'w',l)
         if (taskid.eq.MASTER) then
            call iwrite(6,version,1)
            call iwrite(6,nfft,3)
            call dwrite(6,unita,9)
            call iwrite(6,ispin,1)
            call iwrite(6,ne,2)
            occupation = -1
            call iwrite(6,occupation,1)
         else
            if (pio) then
               call ishift_fileptr(6,4)
               call dshift_fileptr(6,9)
               call ishift_fileptr(6,4)
            end if
         end if
      end if

*     *** read in 3d blocks ***
      do n=1,(ne(1)+ne(2))
         call Dneall_ntoqp(n,q,pj)
         if (pj.eq.taskid_j) then
           call Pack_c_Copy(1,psi2(1+(q-1)*npack1),dcpl_mb(tmp2(1)))
           call Pack_c_unpack(1,dcpl_mb(tmp2(1)))          
         end if
         if (pio) then
            call D3dB_c_write_pio(1,6,dcpl_mb(tmp2(1)),
     >                                dcpl_mb(tmp(1)),pj)
         else
            call D3dB_c_write(1,6,dcpl_mb(tmp2(1)),
     >                            dcpl_mb(tmp(1)),pj)
         end if
      end do

*     **** flush the filepointers ****
      if (pio.and.doflush) call flush_fileptr(6)

*     *** close ELCIN binary file ***
      call ga_sync()
      if (taskid_p.eq.MASTER) then
        call closefile(6)
      end if
      call ga_sync()

      value =           BA_pop_stack(tmp2(2))
      value = value.and.BA_pop_stack(tmp(2))
      if (.not. value) 
     > call errquit('psi_write_filename:error popping stack',0,MA_ERR)

      call nwpw_timing_end(50)
      return
      end


