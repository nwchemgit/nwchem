*
* $Id$
*

*     ***********************************
*     *                                 *
*     *       psi_filefind              *
*     *                                 *
*     ***********************************

      logical function psi_filefind()
      implicit none

#include "tcgmsg.fh"
#include "msgtypesf.h"

*     *** local variables ***
      integer MASTER,taskid
      parameter(MASTER=0)

      logical found
      integer idum,msglen
      character*50 filename
      character*255 full_filename

*     **** external functions ****
      character*50 control_input_psi
      external     control_input_psi


      call Parallel_taskid(taskid)
      if (taskid.eq.MASTER) then
         filename = control_input_psi()
         call util_file_name_noprefix(filename,.false.,
     >                             .false.,
     >                       full_filename)
         inquire(file=full_filename,exist=found)

         idum = 0
         if (found) idum = 1
      end if
      msglen = 1
      call BRDCST(9+MSGINT,idum,mitob(msglen),MASTER)

      if (idum.eq.1) then
         found = .true.
      else
         found = .false.
      end if

      psi_filefind = found
      return
      end


*     ***********************************
*     *                                 *
*     *       epsi_filefind              *
*     *                                 *
*     ***********************************

      logical function epsi_filefind()
      implicit none

#include "tcgmsg.fh"
#include "msgtypesf.h"

*     *** local variables ***
      integer MASTER,taskid
      parameter(MASTER=0)

      logical found
      integer idum,msglen
      character*50 filename
      character*255 full_filename

*     **** external functions ****
      character*50 control_input_epsi
      external     control_input_epsi


      call Parallel_taskid(taskid)
      if (taskid.eq.MASTER) then
         filename = control_input_epsi()
         call util_file_name_noprefix(filename,.false.,
     >                             .false.,
     >                       full_filename)
         inquire(file=full_filename,exist=found)

         idum = 0
         if (found) idum = 1
      end if
      msglen = 1
      call BRDCST(9+MSGINT,idum,mitob(msglen),MASTER)

      if (idum.eq.1) then
         found = .true.
      else
         found = .false.
      end if

      epsi_filefind = found
      return
      end

