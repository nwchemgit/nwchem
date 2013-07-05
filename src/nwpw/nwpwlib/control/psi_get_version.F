*
* $Id$
*

*     ***********************************
*     *             			*
*     *           psi_get_version	*
*     *             			*
*     ***********************************
      integer function psi_get_version()
      implicit none 

*    *** local variables ***
      integer version,l
      character*50 filename
      character*255 full_filename

      integer MASTER,taskid
      parameter(MASTER=0)

*     ***** local functions ****
      character*50 control_input_psi
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
      call Parallel_Brdcst_ivalue(MASTER,version)

      psi_get_version = version
      return
      end



*     ***********************************
*     *                                 *
*     *         epsi_get_version        *
*     *                                 *
*     ***********************************

      integer function epsi_get_version()
      implicit none

*    *** local variables ***
      integer version,l
      character*50 filename
      character*255 full_filename
         
      integer MASTER,taskid
      parameter(MASTER=0)
         
           
*     ***** local functions ****
      character*50 control_input_epsi
      external     control_input_epsi
           
      call Parallel_taskid(taskid)
         
*     **** open ELCIN binary file ****
      if (taskid.eq.MASTER) then
         filename = control_input_epsi()
         call util_file_name_noprefix(filename,.false.,
     >                                .false.,
     >                       full_filename)
         l = index(full_filename,' ') -1
         call openfile(4,full_filename,l,'r',l)
         call iread(4,version,1)
         call closefile(4)
      end if

c     **** send header to all nodes ****
      call Parallel_Brdcst_ivalue(MASTER,version)
        
        
      epsi_get_version = version
      return 
      end


*     ***********************************
*     *                                 *
*     *           v_psi_get_version     *
*     *                                 *
*     ***********************************
      integer function v_psi_get_version()
      implicit none

*    *** local variables ***
      integer version,l
      character*50 filename
      character*255 full_filename

      integer MASTER,taskid
      parameter(MASTER=0)

*     ***** local functions ****
      character*50 control_input_v_psi
      external     control_input_v_psi

      call Parallel_taskid(taskid)

*     **** open ELCIN binary file ****
      if (taskid.eq.MASTER) then
         filename = control_input_v_psi()
         call util_file_name_noprefix(filename,.false.,
     >                                .false.,
     >                       full_filename)
         l = index(full_filename,' ') -1
         call openfile(4,full_filename,l,'r',l)
         call iread(4,version,1)
         call closefile(4)
      end if

c     **** send header to all nodes ****
      call Parallel_Brdcst_ivalue(MASTER,version)

      v_psi_get_version = version
      return
      end



