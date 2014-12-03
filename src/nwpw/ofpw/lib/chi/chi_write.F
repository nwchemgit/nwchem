*
* $Id$
*

*     ***********************************
*     *             			*
*     *           chi_write		*
*     *             			*
*     ***********************************

      subroutine chi_write(ispin,ne,psi2)
      implicit none 
      integer ispin,ne(2)
      double complex psi2(*)

#include "bafdecls.fh"
#include "errquit.fh"


*    *** local variables ***
      integer version,l
      integer nfft3d,npack1,occupation
      integer nfft(3)
      real*8  unita(3,3)
      character*50 filename
      character*255 full_filename

      integer MASTER,taskid,taskid_j
      parameter(MASTER=0)
      integer n,q,pj

c     complex*16 tmp(*)
      integer tmp(2),tmp2(2)
      logical value

*     ***** local functions ****
      character*50 control_output_psi
      external     control_output_psi
      double precision control_unita
      external         control_unita
      integer  control_ngrid,control_version
      external control_ngrid,control_version

      call Parallel_taskid(taskid)
      call Parallel2d_taskid_j(taskid_j)
      call D3dB_nfft3d(1,nfft3d)
      call Pack_npack(1,npack1)

      value = BA_push_get(mt_dcpl,nfft3d,'tmp',tmp(2),tmp(1))
      if (.not. value) 
     > call errquit('chi_write:out of stack memory',0, MA_ERR)
      value = BA_push_get(mt_dcpl,nfft3d,'tmp2',tmp2(2),tmp2(1))
      if (.not. value) 
     > call errquit('chi_write:out of stack memory',0, MA_ERR)

      version = 9  ! version 9 - ofpw, version 5 - band, version 3 and 4 - pspw
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
         write(*,1210) full_filename(1:l)
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
      do n=1,(ispin)
         call Pack_c_Copy(1,psi2(1+(n-1)*npack1),dcpl_mb(tmp2(1)))
         call Pack_c_unpack(1,dcpl_mb(tmp2(1)))          
         call D3dB_c_write(1,6,dcpl_mb(tmp2(1)),
     >                         dcpl_mb(tmp(1)),0)
      end do

*     *** close ELCIN binary file ***
      if (taskid.eq.MASTER) then
        call closefile(6)
      end if


      value = BA_pop_stack(tmp2(2))
      value = BA_pop_stack(tmp(2))
      return
      end


