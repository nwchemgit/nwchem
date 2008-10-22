*
* $Id: psi_read.F,v 1.20 2008-10-22 23:56:44 bylaska Exp $
*

*     ***********************************
*     *             			*
*     *           psi_read		*
*     *             			*
*     ***********************************

      subroutine psi_read(ispin,ne,psi2,occupation,occ)
      implicit none 
      integer ispin,ne(2)
      double complex psi2(*)
      integer occupation
      real*8 occ(*)

#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "errquit.fh"

*    *** local variables ***
      integer version,l,rtdb
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
      logical value,psi_nolattice,lprint

*     ***** local functions ****
      character*30 control_input_psi
      external     control_input_psi
      double precision control_unita
      external         control_unita
      integer  control_ngrid,control_rtdb
      external control_ngrid,control_rtdb
      logical  control_print
      external control_print
      

      call Parallel_taskid(taskid)
      call Parallel2d_taskid_j(taskid_j)
      call D3dB_nfft3d(1,nfft3d)
      call Pack_npack(1,npack1)

      lprint= ((taskid.eq.MASTER).and.control_print(print_low))


      value = MA_push_get(mt_dcpl,nfft3d,'tmp',tmp(2),tmp(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)

      value = MA_push_get(mt_dcpl,nfft3d,'tmp2',tmp2(2),tmp2(1))
       if (.not. value) call errquit('psi_read:pushing stack',0, MA_ERR)

*     **** open ELCIN binary file ****
      if (taskid.eq.MASTER) then
         filename = control_input_psi()
         full_filename = filename
         call util_file_name_resolve(full_filename, .false.)
c        call util_file_name_noprefix(filename,.false.,
c    >                                .false.,
c    >                        full_filename)

         l = index(full_filename,' ') -1
         if (lprint) write(*,1210) full_filename(1:l)
 1210    FORMAT(/' input psi filename:',A)

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
      msglen = 1
      call Parallel_Brdcst_ivalues(MASTER,msglen,occupation)


*     ***** Error checking ****
      if (version.eq.5) then
       call errquit('Error cannot read band structure wavefunctions',1,
     &       INPUT_ERR)
      end if

      if ( (nfft(1).ne.control_ngrid(1)) .or.
     >     (nfft(2).ne.control_ngrid(2)) .or.
     >     (nfft(3).ne.control_ngrid(3)) ) then
        if (taskid.eq.MASTER) then
        write(*,*) "nfft :",nfft
        write(*,*) "ngrid:",control_ngrid(1),
     >                      control_ngrid(2),
     >                      control_ngrid(3)
        end if
        call errquit('Error reading wavefunctions - bad grid', 2,
     &       INPUT_ERR)
      end if

      rtdb = control_rtdb()
      if (.not.rtdb_get(rtdb,'nwpw:psi_nolattice',
     >                  mt_log,1,psi_nolattice))
     >   psi_nolattice = .true.
  
      if (.not.psi_nolattice) then
      if ( (unita(1,1).ne.control_unita(1,1)) .or.
     >     (unita(2,1).ne.control_unita(2,1)) .or.
     >     (unita(3,1).ne.control_unita(3,1)) .or.
     >     (unita(1,2).ne.control_unita(1,2)) .or.
     >     (unita(2,2).ne.control_unita(2,2)) .or.
     >     (unita(3,2).ne.control_unita(3,2)) .or.
     >     (unita(1,3).ne.control_unita(1,3)) .or.
     >     (unita(2,3).ne.control_unita(2,3)) .or.
     >     (unita(3,3).ne.control_unita(3,3)) ) then
        call errquit('Error reading wavefunctions - bad lattice', 3,
     &       INPUT_ERR)
      end if
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

      if (occupation.gt.0) then
         if (taskid.eq.MASTER) then
            call dread(5,occ,(ne(1)+ne(2)))
         end if
         call Parallel_Brdcst_values(MASTER,(ne(1)+ne(2)),occ)
      end if
      

*     *** close ELCIN binary file ***
      if (taskid.eq.MASTER) then
        call closefile(5)
      end if
 
      value = MA_pop_stack(tmp2(2))
      value = value.and.MA_pop_stack(tmp(2))
      if (.not. value) call errquit('psi_read:popping stack',4, MA_ERR)

*     end if

      return
      end


*     ***********************************
*     *             			*
*     *           epsi_read		*
*     *             			*
*     ***********************************

      subroutine epsi_read(ispin,ne,psi2)
      implicit none 
      integer ispin,ne(2)
      double complex psi2(*)
      integer occupation

#include "mafdecls.fh"
#include "rtdb.fh"
#include "errquit.fh"


*    *** local variables ***
      integer version,l,rtdb
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
      logical value,psi_nolattice

*     ***** local functions ****
      character*30 control_input_epsi
      external     control_input_epsi
      double precision control_unita
      external         control_unita
      integer  control_ngrid,control_rtdb
      external control_ngrid,control_rtdb
      

      call Parallel_taskid(taskid)
      call Parallel2d_taskid_j(taskid_j)
      call D3dB_nfft3d(1,nfft3d)
      call Pack_npack(1,npack1)


      value = MA_push_get(mt_dcpl,nfft3d,'tmp',tmp(2),tmp(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)

      value = MA_push_get(mt_dcpl,nfft3d,'tmp2',tmp2(2),tmp2(1))
       if (.not. value) call errquit('psi_read:pushing stack',0, MA_ERR)

*     **** open ELCIN binary file ****
      if (taskid.eq.MASTER) then
         filename = control_input_epsi()
         full_filename = filename
         call util_file_name_resolve(full_filename, .false.)

         l = index(full_filename,' ') -1
         write(*,1220) full_filename(1:l)
 1220    FORMAT(/' input epsi filename:',A)

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
      if (version.eq.5) then
       call errquit('Error cannot read band structure wavefunctions',1,
     &       INPUT_ERR)
      end if

      if ( (nfft(1).ne.control_ngrid(1)) .or.
     >     (nfft(2).ne.control_ngrid(2)) .or.
     >     (nfft(3).ne.control_ngrid(3)) ) then
        if (taskid.eq.MASTER) then
        write(*,*) "nfft :",nfft
        write(*,*) "ngrid:",control_ngrid(1),
     >                      control_ngrid(2),
     >                      control_ngrid(3)
        end if
        call errquit('Error reading wavefunctions - bad grid', 2,
     &       INPUT_ERR)
      end if

      rtdb = control_rtdb()
      if (.not.rtdb_get(rtdb,'nwpw:psi_nolattice',
     >                  mt_log,1,psi_nolattice))
     >   psi_nolattice = .true.
  
      if (.not.psi_nolattice) then
      if ( (unita(1,1).ne.control_unita(1,1)) .or.
     >     (unita(2,1).ne.control_unita(2,1)) .or.
     >     (unita(3,1).ne.control_unita(3,1)) .or.
     >     (unita(1,2).ne.control_unita(1,2)) .or.
     >     (unita(2,2).ne.control_unita(2,2)) .or.
     >     (unita(3,2).ne.control_unita(3,2)) .or.
     >     (unita(1,3).ne.control_unita(1,3)) .or.
     >     (unita(2,3).ne.control_unita(2,3)) .or.
     >     (unita(3,3).ne.control_unita(3,3)) ) then
        call errquit('Error reading wavefunctions - bad lattice', 3,
     &       INPUT_ERR)
      end if
      end if


*     *** read in 3d blocks - data replicated ***
      do n=1,(ne(1)+ne(2))
c         call Dnexall_ntoqp(n,q,pj)
          
         call D3dB_c_read(1,5,dcpl_mb(tmp2(1)),
     >                        dcpl_mb(tmp(1)),-1)
         call Pack_c_pack(1,dcpl_mb(tmp2(1)))
         call Pack_c_Copy(1,dcpl_mb(tmp2(1)),psi2(1+(n-1)*npack1))

c         if (pj.eq.taskid_j) then
c           call Pack_c_pack(1,dcpl_mb(tmp2(1)))
c           call Pack_c_Copy(1,dcpl_mb(tmp2(1)),psi2(1+(q-1)*npack1))
c         end if
      end do

*     *** close ELCIN binary file ***
      if (taskid.eq.MASTER) then
        call closefile(5)
      end if
 
      value = MA_pop_stack(tmp2(2))
      value = value.and.MA_pop_stack(tmp(2))
      if (.not. value) call errquit('psi_read:popping stack',4, MA_ERR)

*     end if

      return
      end


*     ***********************************
*     *             			*
*     *           psi_read_noocc	*
*     *             			*
*     ***********************************

      subroutine psi_read_noocc(ispin,ne,psi2)
      implicit none 
      integer ispin,ne(2)
      double complex psi2(*)
      integer occupation

#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "errquit.fh"

*    *** local variables ***
      integer version,l,rtdb
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
      logical value,psi_nolattice,lprint

*     ***** local functions ****
      character*30 control_input_psi
      external     control_input_psi
      double precision control_unita
      external         control_unita
      integer  control_ngrid,control_rtdb
      external control_ngrid,control_rtdb
      logical  control_print
      external control_print
      

      call Parallel_taskid(taskid)
      call Parallel2d_taskid_j(taskid_j)
      call D3dB_nfft3d(1,nfft3d)
      call Pack_npack(1,npack1)

      lprint= ((taskid.eq.MASTER).and.control_print(print_low))


      value = MA_push_get(mt_dcpl,nfft3d,'tmp',tmp(2),tmp(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)

      value = MA_push_get(mt_dcpl,nfft3d,'tmp2',tmp2(2),tmp2(1))
       if (.not. value) call errquit('psi_read:pushing stack',0, MA_ERR)

*     **** open ELCIN binary file ****
      if (taskid.eq.MASTER) then
         filename = control_input_psi()
         full_filename = filename
         call util_file_name_resolve(full_filename, .false.)
c        call util_file_name_noprefix(filename,.false.,
c    >                                .false.,
c    >                        full_filename)

         l = index(full_filename,' ') -1
         if (lprint) write(*,1210) full_filename(1:l)
 1210    FORMAT(/' input psi filename:',A)

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
      if (version.eq.5) then
       call errquit('Error cannot read band structure wavefunctions',1,
     &       INPUT_ERR)
      end if

      if ( (nfft(1).ne.control_ngrid(1)) .or.
     >     (nfft(2).ne.control_ngrid(2)) .or.
     >     (nfft(3).ne.control_ngrid(3)) ) then
        if (taskid.eq.MASTER) then
        write(*,*) "nfft :",nfft
        write(*,*) "ngrid:",control_ngrid(1),
     >                      control_ngrid(2),
     >                      control_ngrid(3)
        end if
        call errquit('Error reading wavefunctions - bad grid', 2,
     &       INPUT_ERR)
      end if

      rtdb = control_rtdb()
      if (.not.rtdb_get(rtdb,'nwpw:psi_nolattice',
     >                  mt_log,1,psi_nolattice))
     >   psi_nolattice = .true.
  
      if (.not.psi_nolattice) then
      if ( (unita(1,1).ne.control_unita(1,1)) .or.
     >     (unita(2,1).ne.control_unita(2,1)) .or.
     >     (unita(3,1).ne.control_unita(3,1)) .or.
     >     (unita(1,2).ne.control_unita(1,2)) .or.
     >     (unita(2,2).ne.control_unita(2,2)) .or.
     >     (unita(3,2).ne.control_unita(3,2)) .or.
     >     (unita(1,3).ne.control_unita(1,3)) .or.
     >     (unita(2,3).ne.control_unita(2,3)) .or.
     >     (unita(3,3).ne.control_unita(3,3)) ) then
        call errquit('Error reading wavefunctions - bad lattice', 3,
     &       INPUT_ERR)
      end if
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

c      if (occupation.gt.0) then
c         if (taskid.eq.MASTER) then
c            call dread(5,occ,(ne(1)+ne(2)))
c         end if
c         call Parallel_Brdcst_values(MASTER,(ne(1)+ne(2)),occ)
c      end if
      

*     *** close ELCIN binary file ***
      if (taskid.eq.MASTER) then
        call closefile(5)
      end if
 
      value = MA_pop_stack(tmp2(2))
      value = value.and.MA_pop_stack(tmp(2))
      if (.not. value) call errquit('psi_read:popping stack',4, MA_ERR)

*     end if

      return
      end

