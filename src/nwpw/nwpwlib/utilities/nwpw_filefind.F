*
* $Id$
*

*     ***********************************
*     *					*
*     *	 	  nwpw_filefind 	*
*     *					*
*     ***********************************

      logical function nwpw_filefind(fname)
      implicit none
      character*50 fname

#include "tcgmsg.fh"
#include "msgtypesf.h"

*     *** local variables ***
      integer MASTER,taskid
      parameter(MASTER=0)

      logical found
      integer idum
      character*255 full_filename

      call Parallel_taskid(taskid)
      idum = 0
      if (taskid.eq.MASTER) then
         call util_file_name_noprefix(fname,.false.,
     >                             .false.,
     >                       full_filename)
         inquire(file=full_filename,exist=found)

         idum = 0
         if (found) idum = 1
      end if
      call Parallel_Brdcst_ivalue(MASTER,idum)
      

      if (idum.eq.1) then
         found = .true.
      else
         found = .false.
      end if

      nwpw_filefind = found
      return
      end

