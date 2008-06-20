*
* $Id: psi_get_header.F,v 1.5 2008-06-20 16:20:24 bylaska Exp $
*


*     ***********************************
*     *             			*
*     *           psi_get_header	*
*     *             			*
*     ***********************************

      subroutine psi_get_header(version,nfft,unita,ispin,ne)
      implicit none 
      integer version
      integer nfft(3)
      real*8  unita(3,3)
      integer ispin,ne(2)

*    *** local variables ***
      integer l
      character*30 filename
      character*255 full_filename

      integer MASTER,taskid
      parameter(MASTER=0)
      integer msglen


*     ***** local functions ****
      character*30 control_input_psi
      external     control_input_psi

      call Parallel_taskid(taskid)


*     **** open ELCIN binary file ****
      if (taskid.eq.MASTER) then
         filename = control_input_psi()
         call util_file_name_noprefix(filename,.false.,
     >                                .false.,
     >                       full_filename)
         l = index(full_filename,' ') -1
         call openfile(4,full_filename,l,'r',l)
         call iread(4,version,1)
         call iread(4,nfft,3)
         call dread(4,unita,9)
         call iread(4,ispin,1)
         call iread(4,ne,2)
         call closefile(4)
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

      return
      end


*     ***********************************
*     *             			*
*     *     psi_get_header_filename	*
*     *             			*
*     ***********************************

      subroutine psi_get_header_filename(filename,
     >                          version,nfft,unita,ispin,ne)
      implicit none 
      character*(*) filename
      integer version
      integer nfft(3)
      real*8  unita(3,3)
      integer ispin,ne(2)


*    *** local variables ***
      integer l
      character*255 full_filename

      integer MASTER,taskid
      parameter(MASTER=0)
      integer msglen


*     ***** local functions ****

      call Parallel_taskid(taskid)


*     **** open ELCIN binary file ****
      if (taskid.eq.MASTER) then
         call util_file_name_noprefix(filename,.false.,
     >                                .false.,
     >                       full_filename)
         l = index(full_filename,' ') -1
         call openfile(4,full_filename,l,'r',l)
         call iread(4,version,1)
         call iread(4,nfft,3)
         call dread(4,unita,9)
         call iread(4,ispin,1)
         call iread(4,ne,2)
         call closefile(4)
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

      return
      end

