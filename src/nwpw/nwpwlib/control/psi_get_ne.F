*
* $Id$
*

*     ***********************************
*     *             			*
*     *           psi_get_ne		*
*     *             			*
*     ***********************************

      subroutine psi_get_ne(ispin,ne)
      implicit none 
      integer ispin,ne(2)

#include	"stdio.fh"
#include	"util.fh"

*    *** local variables ***
      logical oprint
      integer version,l
      integer nfft(3)
      real*8  unita(3,3)
      character*50 filename
      character*255 full_filename

      integer MASTER,taskid
      parameter(MASTER=0)
      integer msglen


*     ***** local functions ****
      logical      control_print
      external     control_print
      character*50 control_input_psi
      external     control_input_psi
      double precision control_unita
      external         control_unita
      integer  control_ngrid
      external control_ngrid

      call Parallel_taskid(taskid)
      oprint = ((taskid.eq.MASTER).and.control_print(print_medium))

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


*     ***** Error checking ****
      if ( (nfft(1).ne.control_ngrid(1)) .or.
     >     (nfft(2).ne.control_ngrid(2)) .or.
     >     (nfft(3).ne.control_ngrid(3)) ) then
         if (oprint) then
            write(luout,*) "Error reading psi - bad grid"
            write(luout,*) "nfft :",nfft
            write(luout,*) "ngrid:",control_ngrid(1),
     >                              control_ngrid(2),
     >                              control_ngrid(3)

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
     >     (unita(3,3).ne.control_unita(3,3)) ) then
         if (oprint) then
          write(luout,*) "Error reading psi - bad unitcell"
          write(luout,*) " - Ignored if the following  parameter is set"
          write(luout,*) " - set nwpw:psi_nolattice .true."
         end if

      end if


      return
      end

*     ***********************************
*     *             			*
*     *        psi_get_ne_excited	*
*     *             			*
*     ***********************************

      subroutine psi_get_ne_excited(ispin,ne)
      implicit none 
      integer ispin,ne(2)

#include	"stdio.fh"
#include	"util.fh"

*    *** local variables ***
      logical oprint
      integer version,l
      integer nfft(3)
      real*8  unita(3,3)
      character*50 filename
      character*255 full_filename

      integer MASTER,taskid
      parameter(MASTER=0)
      integer msglen


*     ***** local functions ****
      logical      control_print
      external     control_print
      character*50 control_input_epsi
      external     control_input_epsi
      double precision control_unita
      external         control_unita
      integer  control_ngrid
      external control_ngrid

      call Parallel_taskid(taskid)
      oprint = ((taskid.eq.MASTER).and.control_print(print_medium))


*     **** open ELCIN binary file ****
      if (taskid.eq.MASTER) then
         filename = control_input_epsi()
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


*     ***** Error checking ****
      if ( (nfft(1).ne.control_ngrid(1)) .or.
     >     (nfft(2).ne.control_ngrid(2)) .or.
     >     (nfft(3).ne.control_ngrid(3)) ) then
         if (oprint) then
            write(luout,*) "Error reading psi - bad grid"
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
     >     (unita(3,3).ne.control_unita(3,3)) ) then
         if (oprint) then
          write(luout,*) "Error reading psi - bad unitcell"
          write(luout,*) " - Ignored if the following  parameter is set"
          write(luout,*) " - set nwpw:psi_nolattice .true."
         end if

      end if


      return
      end


*     ***********************************
*     *             			*
*     *       psi_get_ne_occupation	*
*     *             			*
*     ***********************************

      subroutine psi_get_ne_occupation(ispin,ne,occupation)
      implicit none 
      integer ispin,ne(2),occupation

#include	"stdio.fh"
#include	"util.fh"

*    *** local variables ***
      logical oprint
      integer version,l
      integer nfft(3),nbrill
      real*8  unita(3,3)
      character*50 filename
      character*255 full_filename

      integer MASTER,taskid
      parameter(MASTER=0)
      integer msglen


*     ***** local functions ****
      logical      control_print
      external     control_print
      character*50 control_input_psi
      external     control_input_psi
      double precision control_unita
      external         control_unita
      integer  control_ngrid
      external control_ngrid

      call Parallel_taskid(taskid)
      oprint = ((taskid.eq.MASTER).and.control_print(print_medium))

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
         nbrill = 1
         if (version.eq.5) call iread(4,nbrill,1)
         call iread(4,occupation,1)
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
      msglen = 1
      call Parallel_Brdcst_ivalues(MASTER,msglen,occupation)


*     ***** Error checking ****
      if ( (nfft(1).ne.control_ngrid(1)) .or.
     >     (nfft(2).ne.control_ngrid(2)) .or.
     >     (nfft(3).ne.control_ngrid(3)) ) then
         if (oprint) then
            write(luout,*) "Error reading psi - bad grid"
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
     >     (unita(3,3).ne.control_unita(3,3)) ) then
         if (oprint) then
          write(luout,*) "Error reading psi - bad unitcell"
          write(luout,*) " - Ignored if the following  parameter is set"
          write(luout,*) " - set nwpw:psi_nolattice .true."
         end if

      end if


      return
      end

