*
* $Id$
*

*     ***********************************
*     *                                 *
*     *         psi_marat_write         *
*     *                                 *
*     ***********************************

      subroutine psi_marat_write(filename,ispin,ne,psi2)
      implicit none 
      integer ispin,ne(2)
      double complex psi2(*)
      character*(*) filename

#include "bafdecls.fh"


*    *** local variables ***
      integer version,l
      integer nfft3d,npack1
      integer nfft(3)
      real*8  unita(3,3),scale,occupation
      character*255 full_filename

      integer MASTER,taskid
      parameter(MASTER=0)
      integer n

c     complex*16 tmp(*)
      integer tmp(2),tmp2(2),psi_spin(2),psi_occup(2)
      logical value

*     ***** local functions ****
      character*50 control_output_psi
      external     control_output_psi
      double precision control_unita,lattice_omega
      external         control_unita,lattice_omega
      integer  control_ngrid,control_version
      external control_ngrid,control_version

      call Parallel_taskid(taskid)
      call D3dB_nfft3d(1,nfft3d)
      call Pack_npack(1,npack1)

      value = BA_push_get(mt_dcpl,nfft3d,'tmp',tmp(2),tmp(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,nfft3d,'tmp2',tmp2(2),tmp2(1))
      value = value.and.
     >        BA_push_get(mt_dbl,(ne(1)+ne(2)),
     >                    'psi_occup',psi_occup(2),psi_occup(1))
      value = value.and.
     >        BA_push_get(mt_int,(ne(1)+ne(2)),
     >                    'psi_spin',psi_spin(2),psi_spin(1))
      if (.not. value) call errquit('psi_marat_write: push stack',0,0)

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
      scale = 1.0d0/dsqrt(lattice_omega())

      occupation = 1.0d0
      if (ispin.eq.1) occupation = 2.0d0
      do n=1,ne(2)
        dbl_mb(psi_occup(1)+n-1) = occupation
        int_mb(psi_spin(1)+n-1)  = -1
      end do
      do n=ne(2)+1,ne(1)+ne(2)
        dbl_mb(psi_occup(1)+n-1) = occupation
        int_mb(psi_spin(1)+n-1)  = 1
      end do

*     **** open ELCIN binary file ****
      if (taskid.eq.MASTER) then
c         filename = "marat.psi"
         call util_file_name_noprefix(filename,.false.,
     >                                .false.,
     >                       full_filename)
         l = index(full_filename,' ') -1
         OPEN(87,file=full_filename,form="FORMATTED")
         write(87,*) nfft
         write(87,*) (dbl_mb(psi_occup(1)+n-1),n=1,(ne(1)+ne(2)))
         write(87,*) (int_mb(psi_spin(1)+n-1), n=1,(ne(1)+ne(2)))
      end if

*     *** write out 3d blocks ***
      !*** down spin ***
      if (ispin.eq.2) then
        do n=ne(1)+1,(ne(1)+ne(2))
         call Pack_c_Copy(1,psi2(1+(n-1)*npack1),dcpl_mb(tmp2(1)))
         call Pack_c_unpack(1,dcpl_mb(tmp2(1)))         
         call D3dB_cr_fft3b(1,dcpl_mb(tmp2(1))) 
c         call D3dB_r_SMul(1,scale,
c     >                      dcpl_mb(tmp2(1)),
c     >                      dcpl_mb(tmp2(1)))
         call D3dB_r_SMul1(1,scale,dcpl_mb(tmp2(1)))
         call D3dB_r_FormatWrite_paw(1,87,
     >                           dcpl_mb(tmp2(1)),
     >                           dcpl_mb(tmp(1)))
  
        end do
      end if
      !*** up spin ***
      do n=1,ne(1)
         call Pack_c_Copy(1,psi2(1+(n-1)*npack1),dcpl_mb(tmp2(1)))
         call Pack_c_unpack(1,dcpl_mb(tmp2(1)))         
         call D3dB_cr_fft3b(1,dcpl_mb(tmp2(1))) 
c         call D3dB_r_SMul(1,scale,
c     >                      dcpl_mb(tmp2(1)),
c     >                      dcpl_mb(tmp2(1)))
         call D3dB_r_SMul1(1,scale,dcpl_mb(tmp2(1)))
         call D3dB_r_FormatWrite_paw(1,87,
     >                           dcpl_mb(tmp2(1)),
     >                           dcpl_mb(tmp(1)))
      end do

*     *** close ELCIN binary file ***
      if (taskid.eq.MASTER) then
        close(87)
      end if

      value =           BA_pop_stack(psi_spin(2))
      value = value.and.BA_pop_stack(psi_occup(2))
      value = value.and.BA_pop_stack(tmp2(2))
      value = value.and.BA_pop_stack(tmp(2))
      return
      end


