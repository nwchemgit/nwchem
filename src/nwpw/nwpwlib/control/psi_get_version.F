*
* $Id: psi_get_version.F,v 1.2 2008-06-20 16:20:24 bylaska Exp $
*

*     ***********************************
*     *             			*
*     *           psi_get_version		*
*     *             			*
*     ***********************************

      integer function psi_get_version()
      implicit none 

*    *** local variables ***
      integer version,l
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
         call closefile(4)
      end if

c     **** send header to all nodes ****
      msglen = 1
      call Parallel_Brdcst_ivalues(MASTER,msglen,version)


      psi_get_version = version
      return
      end
