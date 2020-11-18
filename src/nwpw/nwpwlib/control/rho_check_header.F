*
* $Id$
*


*     ***********************************
*     *             			*
*     *           rho_check_header	*
*     *             			*
*     ***********************************

      logical function rho_check_header(ispin0,band)
      implicit none 
      integer ispin0
      logical band

#include "stdio.fh"
#include "errquit.fh"

*    *** local variables ***
      logical filename_exist,ok
      integer l,n
      character*50 filename
      character*255 full_filename

      integer MASTER,taskid
      parameter(MASTER=0)
      integer msglen
      integer version,nfft(3),ispin,nbrill
      real*8 unita(3,3)


*     ***** local functions ****
      character*50 control_input_rho
      external     control_input_rho
      integer  control_ngrid,control_version,control_num_kvectors
      external control_ngrid,control_version,control_num_kvectors
      real*8   control_unita
      external control_unita

      call Parallel_taskid(taskid)

*     **** open RHOIN binary file ****
      filename = control_input_rho()

      ok = .true.
      if (taskid.eq.MASTER) then
         call util_file_name_noprefix(filename,.false.,
     >                                .false.,
     >                       full_filename)
         inquire(file=full_filename,exist=filename_exist)
         if (filename_exist) then
            n = 1
            l = index(full_filename,' ') -1
           
            call openfile(4,full_filename,l,'r',l)
            call iread(4,version,1)
            call iread(4,nfft,3)
            call dread(4,unita,9)
            call iread(4,ispin,1)
            if (band) call iread(4,nbrill,1)
            call closefile(4)
         else
            n = 0
            !write(luout,*) "filename does not exist: "//full_filename
         end if
      end if

      msglen = 1
      call Parallel_Brdcst_ivalues(MASTER,msglen,n)
      if (n.eq.0) then
         ok = .false.
      end if

      if (ok) then
c        **** send header to all nodes ****
         msglen = 1
         call Parallel_Brdcst_ivalues(MASTER,msglen,version)
         msglen = 3
         call Parallel_Brdcst_ivalues(MASTER,msglen,nfft)
         msglen = 9
         call Parallel_Brdcst_values(MASTER,msglen,unita)
         msglen = 1
         call Parallel_Brdcst_ivalues(MASTER,msglen,ispin)
         if (band) then
            msglen = 1
            call Parallel_Brdcst_ivalues(MASTER,msglen,nbrill)
            if (nbrill.ne.control_num_kvectors()) ok = .false.
         end if
         if ((control_version().ne.version).or.
     >       (ispin.ne.ispin0).or.
     >       (control_ngrid(1).ne.nfft(1)).or.
     >       (control_ngrid(2).ne.nfft(2)).or.
     >       (control_ngrid(3).ne.nfft(3)).or.
     >       (dabs(control_unita(1,1)-unita(1,1)).gt.1.0d-6).or.
     >       (dabs(control_unita(2,1)-unita(2,1)).gt.1.0d-6).or.
     >       (dabs(control_unita(3,1)-unita(3,1)).gt.1.0d-6).or.
     >       (dabs(control_unita(1,2)-unita(1,2)).gt.1.0d-6).or.
     >       (dabs(control_unita(2,2)-unita(2,2)).gt.1.0d-6).or.
     >       (dabs(control_unita(3,2)-unita(3,2)).gt.1.0d-6).or.
     >       (dabs(control_unita(1,3)-unita(1,3)).gt.1.0d-6).or.
     >       (dabs(control_unita(2,3)-unita(2,3)).gt.1.0d-6).or.
     >       (dabs(control_unita(3,3)-unita(3,3)).gt.1.0d-6)) then
            ok = .false.
         end if
       
      end if

      rho_check_header = ok
      return
      end


