*
* $Id: v_psi_read.F,v 1.7 2008-05-13 02:10:36 bylaska Exp $
*

*     ***********************************
*     *             			*
*     *           v_psi_read		*
*     *             			*
*     ***********************************

      subroutine v_psi_read(ispin0,ne0,psi2)
      implicit none 
      integer ispin0,ne0(2)
      double complex psi2(*)
      integer occupation

#include "mafdecls.fh"
#include "errquit.fh"

*    *** local variables ***
      integer ispin,ne(2)
      integer version,l
      integer nfft3d,npack1
      integer nfft(3)
      real*8  unita(3,3)
      character*30 filename
      character*255 full_filename

      integer MASTER,taskid,taskid_j
      parameter(MASTER=0)
      integer n,q,pj
      integer msglen

c     complex*16 tmp(*)
      integer tmp(2),tmp2(2)
      logical value,readblock

*     ***** local functions ****
      character*30 control_input_v_psi
      external     control_input_v_psi
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
     > call errquit('v_psi_read:out of stack memory',0,MA_ERR)

*     **** open ELCIN binary file ****
      if (taskid.eq.MASTER) then
         filename = control_input_v_psi()
         full_filename = filename
         call util_file_name_resolve(full_filename, .false.)
c        call util_file_name_noprefix(filename,.false.,
c    >                                .false.,
c    >                        full_filename)
         l = index(full_filename,' ') -1
         write(*,1210) full_filename(1:l)
 1210    FORMAT(/' input vpsi filename:',A)


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
      call Parallel_Brdcst_ivalues(MASTER,msglen,occupation)
      msglen = 2
      call Parallel_Brdcst_ivalues(MASTER,msglen,ne)


*     ***** Error checking ****
      readblock=.true.
      if ( (nfft(1).ne.control_ngrid(1)) .or.
     >     (nfft(2).ne.control_ngrid(2)) .or.
     >     (nfft(3).ne.control_ngrid(3)) ) then
         readblock = .false.
         if (taskid.eq.MASTER) then
            write(*,*) "Error reading v_psi, v_psi set to zero"
            write(*,*) " - bad ngrid "
            write(*,*) " - vpsi_ngrid(1)=",nfft(1),
     >                 " ngrid(1)=",control_ngrid(1)
            write(*,*) " - vpsi_ngrid(2)=",nfft(2),
     >                 " ngrid(2)=",control_ngrid(2)
            write(*,*) " - vpsi_ngrid(3)=",nfft(3),
     >                 " ngrid(3)=",control_ngrid(3)
         end if
      end if

      if ( (unita(1,1).ne.control_unita(1,1)) .or.
     >     (unita(2,1).ne.control_unita(2,1)) .or.
     >     (unita(3,1).ne.control_unita(3,1)) .or.
     >     (unita(1,2).ne.control_unita(1,2)) .or.
     >     (unita(2,2).ne.control_unita(2,2)) .or.
     >     (unita(3,2).ne.control_unita(3,2)) .or.
     >     (unita(1,3).ne.control_unita(1,3)) .or.
     >     (unita(2,3).ne.control_unita(2,3)) .or.
     >     (unita(3,3).ne.control_unita(3,3)) .or.
     >     (ispin     .ne.ispin0)             .or.
     >     (ne(1)     .ne.ne0(1))             .or.
     >     (ne(2)     .ne.ne0(2))  ) then
         readblock = .false.
         if (taskid.eq.MASTER) then
            write(6,*) "Error reading v_psi, v_psi set to zero"
            write(6,*) " - incorrect unitcell,ispin, or ne"
         end if
      end if


*     *** read in 3d blocks ***
      if (readblock) then
      do n=1,(ne(1)+ne(2))
         call Dneall_ntoqp(n,q,pj) 
         call D3dB_c_read(1,5,dcpl_mb(tmp2(1)),
     >                        dcpl_mb(tmp(1)),pj)
         if (pj.eq.taskid_j) then
           call Pack_c_pack(1,dcpl_mb(tmp2(1)))
           call Pack_c_Copy(1,dcpl_mb(tmp2(1)),psi2(1+(q-1)*npack1))
         end if
      end do
      end if

*     *** close ELCIN binary file ***
      if (taskid.eq.MASTER) then
        call closefile(5)
      end if

      value =           MA_pop_stack(tmp2(2))
      value = value.and.MA_pop_stack(tmp(2))
      if (.not. value) 
     >  call errquit('v_psi_read:error popping stack',0,MA_ERR)
      return
      end
