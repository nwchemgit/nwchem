*
* $Id: psi_read_filename.F,v 1.4 2007-09-24 16:58:14 bylaska Exp $
*

*     ***********************************
*     *             			*
*     *           psi_read_filename		*
*     *             			*
*     ***********************************

      subroutine psi_read_filename(filename,ispin,ne,psi2)
      implicit none 
      character*(*) filename
      integer ispin,ne(2)
      double complex psi2(*)
      integer occupation

#include "mafdecls.fh"
#include "errquit.fh"

*    *** local variables ***
      integer version,l
      integer nfft3d,npack1
      integer nfft(3)
      real*8  unita(3,3)
      character*255 full_filename

      integer MASTER,taskid,taskid_j
      parameter(MASTER=0)
      integer n,q,pj
      integer msglen

c     complex*16 tmp(*)
      integer tmp(2),tmp2(2)
      logical value

*     ***** local functions ****
      double precision control_unita
      external         control_unita
      integer  control_ngrid
      external control_ngrid

      call Parallel_taskid(taskid)
      call Parallel2d_taskid_j(taskid_j)
      call D3dB_nfft3d(1,nfft3d)
      call Pack_npack(1,npack1)


      value = MA_push_get(mt_dcpl,nfft3d,'tmp',tmp(2),tmp(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,nfft3d,'tmp2',tmp2(2),tmp2(1))
      if (.not. value) 
     >    call errquit('psi_read_filename:out of stack',0,MA_ERR)

*     **** open ELCIN binary file ****
      if (taskid.eq.MASTER) then
         call util_file_name_noprefix(filename,.false.,
     >                                .false.,
     >                        full_filename)
         l = index(full_filename,' ') -1
         call openfile(5,full_filename,l,'r',l)
         call iread(5,version,1)
         call iread(5,nfft,3)
         call dread(5,unita,9)
         call iread(5,ispin,1)
         call iread(5,ne,2)
         call iread(5,occupation,1)
      end if

c     **** send header to all nodes ****
      msglen = 1
      call Parallel_Brdcst_ivalues(MASTER,msglen,version)
      msglen = 3
      call Parallel_Brdcst_ivalues(MASTER,msglen,nfft)
      msglen = 9
      call Parallel_Brdcst_values(MASTER,msglen,unita)
      msglen = 1
      call Parallel_Brdcst_ivalues(MASTER,msglen,ispin)
      msglen = 2
      call Parallel_Brdcst_ivalues(MASTER,msglen,ne)


*     ***** Error checking ****
      if ( (nfft(1).ne.control_ngrid(1)) .or.
     >     (nfft(2).ne.control_ngrid(2)) .or.
     >     (nfft(3).ne.control_ngrid(3)) ) then
        call errquit('Error reading wavefunctions', 0, DISK_ERR)
      end if

      if ( (unita(1,1).ne.control_unita(1,1)) .or.
     >     (unita(2,1).ne.control_unita(2,1)) .or.
     >     (unita(3,1).ne.control_unita(3,1)) .or.
     >     (unita(1,2).ne.control_unita(1,2)) .or.
     >     (unita(2,2).ne.control_unita(2,2)) .or.
     >     (unita(3,2).ne.control_unita(3,2)) .or.
     >     (unita(1,3).ne.control_unita(1,3)) .or.
     >     (unita(2,3).ne.control_unita(2,3)) .or.
     >     (unita(3,3).ne.control_unita(3,3)) ) then
        call errquit('Error reading wavefunctions', 1, DISK_ERR)
      end if


*     *** read in 3d blocks ***
      do n=1,(ne(1)+ne(2))
         call Dneall_ntoqp(n,q,pj)
          
         call D3dB_c_read(1,5,dcpl_mb(tmp2(1)),
     >                        dcpl_mb(tmp(1)),pj)
         if (pj.eq.taskid_j) then
           call Pack_c_pack(1,dcpl_mb(tmp2(1)))
           call Pack_c_Copy(1,dcpl_mb(tmp2(1)),psi2(1+(q-1)*npack1))
         end if
      end do

*     *** close ELCIN binary file ***
      if (taskid.eq.MASTER) then
        call closefile(5)
      end if
 
      value =           MA_pop_stack(tmp2(2))
      value = value.and.MA_pop_stack(tmp(2))
      if (.not. value) 
     >  call errquit('psi_read_filename:error popping stack',0,MA_ERR)

*     end if

      return
      end
