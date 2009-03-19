*
* $Id: psi_write.F,v 1.17 2009-03-19 20:42:17 bylaska Exp $
*

*     ***********************************
*     *             			*
*     *           psi_write		*
*     *             			*
*     ***********************************

      subroutine psi_write(ispin,ne,psi2,occupation,occ)
      implicit none 
      integer ispin,ne(2)
      double complex psi2(*)
      integer occupation
      double precision occ(*)

#include "mafdecls.fh"
#include "util.fh"
#include "errquit.fh"

*    *** local variables ***
      integer version,l
      integer nfft3d,npack1
      integer nfft(3)
      real*8  unita(3,3)
      character*30 filename
      character*255 full_filename

      integer MASTER,taskid,taskid_j
      parameter(MASTER=0)
      integer n,q,pj

c     complex*16 tmp(*)
      integer tmp(2),tmp2(2)
      logical value,lprint

*     ***** local functions ****
      character*30 control_output_psi
      external     control_output_psi
      double precision control_unita
      external         control_unita
      integer  control_ngrid,control_version
      external control_ngrid,control_version
      logical  control_print
      external control_print

      call Parallel_taskid(taskid)
      call Parallel2d_taskid_j(taskid_j)
      call D3dB_nfft3d(1,nfft3d)
      call Pack_npack(1,npack1)
      lprint= ((taskid.eq.MASTER).and.control_print(print_low))

      value = MA_push_get(mt_dcpl,nfft3d,'tmp',tmp(2),tmp(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,nfft3d,'tmp2',tmp2(2),tmp2(1))
      if (.not. value) 
     > call errquit('psi_write:out of stack memory',0,MA_ERR)

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
      if (taskid.eq.MASTER) then
         filename = control_output_psi()
         
         full_filename = filename
         call util_file_name_resolve(full_filename, .false.)
c        call util_file_name_noprefix(filename,.false.,
c    >                                .false.,
c    >                       full_filename)
         l = index(full_filename,' ') -1
         if (lprint) write(*,1210) full_filename(1:l)
 1210    FORMAT(/' output psi filename:',A)

         call openfile(6,full_filename,l,'w',l)
         call iwrite(6,version,1)
         call iwrite(6,nfft,3)
         call dwrite(6,unita,9)
         call iwrite(6,ispin,1)
         call iwrite(6,ne,2)
         call iwrite(6,occupation,1)
      end if

*     *** write out 3d blocks ***
      do n=1,(ne(1)+ne(2))
         call Dneall_ntoqp(n,q,pj)
         if (pj.eq.taskid_j) then
           call Pack_c_Copy(1,psi2(1+(q-1)*npack1),dcpl_mb(tmp2(1)))
           call Pack_c_unpack(1,dcpl_mb(tmp2(1)))          
         end if
         call D3dB_c_write(1,6,dcpl_mb(tmp2(1)),
     >                         dcpl_mb(tmp(1)),pj)
      end do

*     **** write the occupations - ****
      if (occupation.gt.0) then
         if (taskid.eq.MASTER) then
           call dwrite(6,occ,(ne(1)+ne(2)))
         end if
      end if


*     *** close ELCIN binary file ***
      if (taskid.eq.MASTER) then
        call closefile(6)
      end if

      value =           MA_pop_stack(tmp2(2))
      value = value.and.MA_pop_stack(tmp(2))
      if (.not. value) 
     > call errquit('psi_write:error popping stack',0,MA_ERR)
      return
      end




*     ***********************************
*     *             			*
*     *           epsi_write		*
*     *             			*
*     ***********************************

      subroutine epsi_write(ispin,ne,psi2)
      implicit none 
      integer ispin,ne(2)
      double complex psi2(*)

#include "mafdecls.fh"
#include "util.fh"
#include "errquit.fh"


*    *** local variables ***
      integer occupation
      integer version,l
      integer nfft3d,npack1
      integer nfft(3)
      real*8  unita(3,3)
      character*30 filename
      character*255 full_filename

      integer MASTER,taskid,taskid_j
      parameter(MASTER=0)
      integer n,q,pj

c     complex*16 tmp(*)
      integer tmp(2),tmp2(2)
      logical value,lprint

*     ***** local functions ****
      character*30 control_output_epsi
      external     control_output_epsi
      double precision control_unita
      external         control_unita
      integer  control_ngrid,control_version
      external control_ngrid,control_version
      logical  control_print
      external control_print

      call Parallel_taskid(taskid)
      call Parallel2d_taskid_j(taskid_j)
      call D3dB_nfft3d(1,nfft3d)
      call Pack_npack(1,npack1)

      lprint= ((taskid.eq.MASTER).and.control_print(print_low))


      value = MA_push_get(mt_dcpl,nfft3d,'tmp',tmp(2),tmp(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,nfft3d,'tmp2',tmp2(2),tmp2(1))
      if (.not. value) 
     > call errquit('epsi_write:out of stack memory',0,MA_ERR)

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
      if (taskid.eq.MASTER) then
         filename = control_output_epsi()
         
         full_filename = filename
         call util_file_name_resolve(full_filename, .false.)

         l = index(full_filename,' ') -1
         if (lprint) write(*,1220) full_filename(1:l)
 1220    FORMAT(/' output epsi filename:',A)
         call openfile(6,full_filename,l,'w',l)
         call iwrite(6,version,1)
         call iwrite(6,nfft,3)
         call dwrite(6,unita,9)
         call iwrite(6,ispin,1)
         call iwrite(6,ne,2)
         occupation = -1
         call iwrite(6,occupation,1)
      end if

*     *** read in 3d blocks ***
c      do n=1,(ne(1)+ne(2))
c         call Dnexall_ntoqp(n,q,pj)
c         if (pj.eq.taskid_j) then
c           call Pack_c_Copy(1,psi2(1+(q-1)*npack1),dcpl_mb(tmp2(1)))
c           call Pack_c_unpack(1,dcpl_mb(tmp2(1)))          
c         end if
c         call D3dB_c_write(1,6,dcpl_mb(tmp2(1)),
c     >                         dcpl_mb(tmp(1)),pj)
c      end do
      if (taskid_j.eq.0) then
      do n=1,(ne(1)+ne(2))
         call Pack_c_Copy(1,psi2(1+(n-1)*npack1),dcpl_mb(tmp2(1)))
         call Pack_c_unpack(1,dcpl_mb(tmp2(1)))          
         call D3dB_c_write(1,6,dcpl_mb(tmp2(1)),
     >                         dcpl_mb(tmp(1)),taskid_j)
      end do
      end if

*     *** close ELCIN binary file ***
      if (taskid.eq.MASTER) then
        call closefile(6)
      end if

      value =           MA_pop_stack(tmp2(2))
      value = value.and.MA_pop_stack(tmp(2))
      if (.not. value) 
     > call errquit('epsi_write:error popping stack',0,MA_ERR)
      return
      end



*     ***********************************
*     *             			*
*     *           psi_write_noocc	*
*     *             			*
*     ***********************************

      subroutine psi_write_noocc(ispin,ne,psi2)
      implicit none 
      integer ispin,ne(2)
      double complex psi2(*)
      integer occupation

#include "mafdecls.fh"
#include "util.fh"
#include "errquit.fh"

*    *** local variables ***
      integer version,l
      integer nfft3d,npack1
      integer nfft(3)
      real*8  unita(3,3)
      character*30 filename
      character*255 full_filename

      integer MASTER,taskid,taskid_j
      parameter(MASTER=0)
      integer n,q,pj

c     complex*16 tmp(*)
      integer tmp(2),tmp2(2)
      logical value,lprint

*     ***** local functions ****
      character*30 control_output_psi
      external     control_output_psi
      double precision control_unita
      external         control_unita
      integer  control_ngrid,control_version
      external control_ngrid,control_version
      logical  control_print
      external control_print

      call Parallel_taskid(taskid)
      call Parallel2d_taskid_j(taskid_j)
      call D3dB_nfft3d(1,nfft3d)
      call Pack_npack(1,npack1)
      lprint= ((taskid.eq.MASTER).and.control_print(print_low))

      value = MA_push_get(mt_dcpl,nfft3d,'tmp',tmp(2),tmp(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,nfft3d,'tmp2',tmp2(2),tmp2(1))
      if (.not. value) 
     > call errquit('psi_write:out of stack memory',0,MA_ERR)

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
      if (taskid.eq.MASTER) then
         filename = control_output_psi()
         
         full_filename = filename
         call util_file_name_resolve(full_filename, .false.)
c        call util_file_name_noprefix(filename,.false.,
c    >                                .false.,
c    >                       full_filename)
         l = index(full_filename,' ') -1
         if (lprint) write(*,1210) full_filename(1:l)
 1210    FORMAT(/' output psi filename:',A)

         call openfile(6,full_filename,l,'w',l)
         call iwrite(6,version,1)
         call iwrite(6,nfft,3)
         call dwrite(6,unita,9)
         call iwrite(6,ispin,1)
         call iwrite(6,ne,2)
         occupation = -1
         call iwrite(6,occupation,1)
      end if

*     *** write out 3d blocks ***
      do n=1,(ne(1)+ne(2))
         call Dneall_ntoqp(n,q,pj)
         if (pj.eq.taskid_j) then
           call Pack_c_Copy(1,psi2(1+(q-1)*npack1),dcpl_mb(tmp2(1)))
           call Pack_c_unpack(1,dcpl_mb(tmp2(1)))          
         end if
         call D3dB_c_write(1,6,dcpl_mb(tmp2(1)),
     >                         dcpl_mb(tmp(1)),pj)
      end do

c*     **** write the occupations - ****
c      if (occupation.gt.0) then
c         if (taskid.eq.MASTER) then
c           call dwrite(6,occ,(ne(1)+ne(2)))
c         end if
c      end if


*     *** close ELCIN binary file ***
      if (taskid.eq.MASTER) then
        call closefile(6)
      end if

      value =           MA_pop_stack(tmp2(2))
      value = value.and.MA_pop_stack(tmp(2))
      if (.not. value) 
     > call errquit('psi_write:error popping stack',0,MA_ERR)
      return
      end

