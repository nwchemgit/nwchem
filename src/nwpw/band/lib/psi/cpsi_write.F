*
* $Id: cpsi_write.F,v 1.4 2002-05-03 19:44:02 bylaska Exp $
*

*     ***************************
*     *             		*
*     *        cpsi_write	*
*     *             		*
*     ***************************

      subroutine cpsi_write(ispin,ne,nbrillioun,psi2)
      implicit none 
      integer ispin,ne(2),nbrillioun
      double complex psi2(*)

#include "mafdecls.fh"


*    *** local variables ***
      integer version,l
      integer nfft3d,npack1
      integer nfft(3)
      real*8  unita(3,3)
      character*30 filename
      character*255 full_filename

      integer MASTER,taskid
      parameter(MASTER=0)
      integer n,nb,indx

c     complex*16 tmp(*)
      integer tmp(2),tmp2(2)
      logical value

*     ***** local functions ****
      character*30 control_output_psi
      external     control_output_psi
      double precision control_unita
      external         control_unita
      integer  control_ngrid
      external control_ngrid

      call Parallel_taskid(taskid)
      call C3dB_nfft3d(1,nfft3d)
      call Cram_max_npack(npack1)

      value = MA_push_get(mt_dcpl,nfft3d,'tmp',tmp(2),tmp(1))
      if (.not. value) call errquit('out of stack memory',0)
      value = MA_push_get(mt_dcpl,nfft3d,'tmp2',tmp2(2),tmp2(1))
      if (.not. value) call errquit('out of stack memory',0)

      version = 5 ! version 5 - band, version 3 and 4 - pspw
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
         call util_file_name_noprefix(filename,.false.,
     >                                .false.,
     >                       full_filename)
         l = index(full_filename,' ') -1
         call openfile(6,full_filename,l,'w',l)
         call iwrite(6,version,1)
         call iwrite(6,nfft,3)
         call dwrite(6,unita,9)
         call iwrite(6,ispin,1)
         call iwrite(6,ne,2)
         call iwrite(6,nbrillioun,1)
      end if

*     *** write out 3d blocks ***
      do nb=1,nbrillioun
      do n=1,(ne(1)+ne(2))
         indx = 1+(n-1)*npack1 + (nb-1)*npack1*(ne(1)+ne(2))
         call Cram_c_Copy(nb,psi2(indx),dcpl_mb(tmp2(1)))
         call Cram_c_unpack(nb,dcpl_mb(tmp2(1)))          
         call C3dB_c_write(1,6,dcpl_mb(tmp2(1)),
     >                         dcpl_mb(tmp(1)))
      end do
      end do

*     *** close ELCIN binary file ***
      if (taskid.eq.MASTER) then
        call closefile(6)
      end if

      value = MA_pop_stack(tmp2(2))
      value = MA_pop_stack(tmp(2))
      return
      end


