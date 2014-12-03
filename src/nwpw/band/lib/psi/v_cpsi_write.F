*
* $Id$
*

*     ***************************
*     *             		*
*     *        v_cpsi_write	*
*     *             		*
*     ***************************

      subroutine v_cpsi_write(spin_orbit,ispin,ne,nbrillioun,psi0_tag)
      implicit none 
      logical spin_orbit
      integer ispin,ne(2),nbrillioun
      integer psi0_tag

#include "bafdecls.fh"
#include "errquit.fh"
#include "stdio.fh"
#include "util.fh"

*    *** local variables ***
      integer version,l,ispin0
      integer nfft3d
      integer nfft(3)
      real*8  unita(3,3)
      character*50 filename
      character*255 full_filename

      integer MASTER,taskid,taskid_j,taskid_k
      parameter(MASTER=0)
      integer n,nb,indx,occupation
      integer q,pj,nbq,pk

c     complex*16 tmp(*)
      integer tmp(2),tmp2(2),tmpe(2)
      logical value

*     ***** local functions ****
      logical      control_print
      external     control_print
      character*50 control_output_v_psi
      external     control_output_v_psi
      double precision control_unita
      external         control_unita
      integer  control_ngrid,cpsi_data_get_chnk,cpsi_data_get_next
      integer  cpsi_data_get_ptr
      external control_ngrid,cpsi_data_get_chnk,cpsi_data_get_next
      external cpsi_data_get_ptr

      call Parallel_taskid(taskid)
      call Parallel3d_taskid_j(taskid_j)
      call Parallel3d_taskid_k(taskid_k)
      call C3dB_nfft3d(1,nfft3d)
      occupation = cpsi_data_get_next(psi0_tag)

      value = BA_push_get(mt_dcpl,nfft3d,'tmp',tmp(2),tmp(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)
      value = BA_push_get(mt_dcpl,nfft3d,'tmp2',tmp2(2),tmp2(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)
      value = BA_push_get(mt_dbl,(ne(1)+ne(2)),'tmpe',tmpe(2),tmpe(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)

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
         filename = control_output_v_psi()
         call util_file_name_noprefix(filename,.false.,
     >                                .false.,
     >                       full_filename)
         l = index(full_filename,' ') -1
         if (control_print(print_medium))
     >    write(luout,1210) full_filename(1:l)
 1210    FORMAT(/' output vpsi filename:',A)

         call openfile(6,full_filename,l,'w',l)
         call iwrite(6,version,1)
         call iwrite(6,nfft,3)
         call dwrite(6,unita,9)
         
         if (spin_orbit) then        !** if spin_orbit --> writeout ispin=3 **
            ispin0 = 3
         else
            ispin0 = ispin
         end if
         call iwrite(6,ispin0,1)
         call iwrite(6,ne,2)
         call iwrite(6,nbrillioun,1)
         call iwrite(6,occupation,1)
      end if

*     *** write out 3d blocks ***
      do nb=1,nbrillioun
         call K1dB_ktoqp(nb,nbq,pk)
         do n=1,(ne(1)+ne(2))
         pj = 0
         q = n
         if ((pj.eq.taskid_j).and.(pk.eq.taskid_k)) then
            indx = cpsi_data_get_ptr(psi0_tag,nbq,q)
            call Cram_c_Copy(nbq,dbl_mb(indx),dcpl_mb(tmp2(1)))
            call Cram_c_unpack(nbq,dcpl_mb(tmp2(1)))          
         end if
         call C3dB_c_write(1,6,dcpl_mb(tmp2(1)),
     >                         dcpl_mb(tmp(1)),pj,pk)
         end do
      end do

*     **** write the occupations - needs to be reworked for parallel nbrill,ne ****
      if (occupation.gt.0) then
         l = cpsi_data_get_next(psi0_tag)
         do nb=1,nbrillioun
            call K1dB_ktoqp(nb,nbq,pk)
            call dcopy((ne(1)+ne(2)),0.0d0,0,dbl_mb(tmpe(1)),1)
            if (pk.eq.taskid_k) then
               call dcopy((ne(1)+ne(2)),
     >                 dbl_mb(cpsi_data_get_chnk(l,nbq)),1,
     >                 dbl_mb(tmpe(1)),1)
            end if
            call K1dB_Vector_SumAll((ne(1)+ne(2)),dbl_mb(tmpe(1)))
            if (taskid.eq.MASTER) then
              call dwrite(6,dbl_mb(tmpe(1)),(ne(1)+ne(2)))
            end if
         end do
      end if

*     *** close ELCIN binary file ***
      if (taskid.eq.MASTER) then
        call closefile(6)
      end if

      value =           BA_pop_stack(tmpe(2))
      value = value.and.BA_pop_stack(tmp2(2))
      value = value.and.BA_pop_stack(tmp(2))
      if (.not. value) call errquit('error popping stack',0, MA_ERR)
      return
      end
