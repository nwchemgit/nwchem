*
* $Id$
*

*     ***********************************
*     *             			*
*     *           rho_read		*
*     *             			*
*     ***********************************

      subroutine rho_read(ispin,rho2)
      implicit none 
      integer ispin
      real*8  rho2(*)

#include "bafdecls.fh"
#include "btdb.fh"
#include "util.fh"
#include "errquit.fh"
#include "stdio.fh"

*    *** local variables ***
      integer version,l,rtdb
      integer nfft3d,n2ft3d
      integer nfft(3)
      real*8  unita(3,3)
      character*50 filename
      character*255 full_filename

      integer MASTER,taskid,taskid_i,taskid_j,taskid_p,com_p
      parameter(MASTER=0)
      integer n,q,pj
      integer msglen,ispin0

c     complex*16 tmp(*)
      integer tmp(2)
      logical value,rho_nolattice,lprint,pio,file_exist

*     ***** local functions ****
      character*50 control_input_rho
      external     control_input_rho
      double precision control_unita
      external         control_unita
      integer  control_ngrid,control_rtdb
      external control_ngrid,control_rtdb
      logical  control_print,control_parallel_io
      external control_print,control_parallel_io
      

      call nwpw_timing_start(50)
      call Parallel_taskid(taskid)
      call Parallel2d_taskid_i(taskid_i)
      call Parallel2d_taskid_j(taskid_j)
      call D3dB_nfft3d(1,nfft3d)
      call D3dB_n2ft3d(1,n2ft3d)

      taskid_p = taskid
      com_p = 0

      lprint= ((taskid.eq.MASTER).and.control_print(print_medium))

      value = BA_push_get(mt_dbl,2*nfft3d,'tmp',tmp(2),tmp(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)

*     **** open RHOIN binary file ****
      filename = control_input_rho()
      if (taskid_p.eq.MASTER) then
         full_filename = filename
         call util_file_name_resolve(full_filename, .false.)

         inquire(file=full_filename,exist=file_exist)
         if (file_exist) then
            n = 1
            l = index(full_filename,' ') -1
            if (lprint) write(LuOut,1210) full_filename(1:l)
 1210       FORMAT('    - Reading fractional density ',
     >            '(To turn off, set nwpw:use_fractional_rho .false.).',
     >            /'      Filename:',A)

            call openfile(5,full_filename,l,'r',l)
            call iread(5,version,1)
            call iread(5,nfft,3)
            call dread(5,unita,9)
            call iread(5,ispin0,1)
         
         else
            n = 0
            write(luout,*) 
     >      "Warning - filename does not exist: "//full_filename
            call util_flush(luout)
         end if
      end if

      msglen = 1
      call Parallela_Brdcst_ivalues(com_p,MASTER,msglen,n)
      if (n.eq.0) then
         call errquit('filename does not exist',1,INPUT_ERR)
      end if

c     **** send header to all nodes ****
      msglen = 1
      call Parallela_Brdcst_ivalues(com_p,MASTER,msglen,version)
      msglen = 3
      call Parallela_Brdcst_ivalues(com_p,MASTER,msglen,nfft)
      msglen = 9
      call Parallela_Brdcst_values(com_p,MASTER,msglen,unita)
      msglen = 1
      call Parallela_Brdcst_ivalues(com_p,MASTER,msglen,ispin0)


*     ***** Error checking ****
      if (version.eq.5) then
       call errquit('Error cannot read band structure wavefunctions',1,
     &       INPUT_ERR)
      end if

      if ( (nfft(1).ne.control_ngrid(1)) .or.
     >     (nfft(2).ne.control_ngrid(2)) .or.
     >     (nfft(3).ne.control_ngrid(3)) .or.
     >     (ispin.ne.ispin0) ) then
        if (taskid_p.eq.MASTER) then
        write(LuOut,*) "taskid_i,taskid_j=",taskid_i,taskid_j
        write(LuOut,*) "ispin,ispin0 :",ispin,ispin0
        write(LuOut,*) "nfft :",nfft
        write(LuOut,*) "ngrid:",control_ngrid(1),
     >                          control_ngrid(2),
     >                          control_ngrid(3)
        end if
        call errquit('Error reading rho - bad grid', 2,
     &       INPUT_ERR)
      end if

      rtdb = control_rtdb()
c      if (.not.btdb_get(rtdb,'nwpw:rho_nolattice',
c     >                  mt_log,1,rho_nolattice))
c     >   rho_nolattice = .true.
      if (.not.btdb_get(rtdb,'nwpw:rho_nolattice',
     >                  mt_log,1,rho_nolattice))
     >   rho_nolattice = .true.
  
      if (.not.rho_nolattice) then
      if ( (unita(1,1).ne.control_unita(1,1)) .or.
     >     (unita(2,1).ne.control_unita(2,1)) .or.
     >     (unita(3,1).ne.control_unita(3,1)) .or.
     >     (unita(1,2).ne.control_unita(1,2)) .or.
     >     (unita(2,2).ne.control_unita(2,2)) .or.
     >     (unita(3,2).ne.control_unita(3,2)) .or.
     >     (unita(1,3).ne.control_unita(1,3)) .or.
     >     (unita(2,3).ne.control_unita(2,3)) .or.
     >     (unita(3,3).ne.control_unita(3,3)) ) then
        call errquit('Error reading rho - bad lattice', 3,
     &       INPUT_ERR)
      end if
      end if

      call Parallel_shared_vector_zero(.true.,ispin*n2ft3d,rho2)
      call D3dB_r_read(1,5,rho2,dbl_mb(tmp(1)),-1)
      if (ispin.gt.1) then
         call D3dB_r_read(1,5,rho2(1+n2ft3d),dbl_mb(tmp(1)),-1)
      end if


*     *** close ELCIN binary file ***
      if (taskid_p.eq.MASTER) then
        call closefile(5)
      end if

      value = BA_pop_stack(tmp(2))
      if (.not. value) call errquit('rho_read:popping stack',4, MA_ERR)

*     end if

      call nwpw_timing_end(50)
      return
      end

