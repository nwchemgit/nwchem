*
* $Id: psi_get_version.F,v 1.1 2001-12-20 00:30:52 bylaska Exp $
*

#define TCGMSG
*     ***********************************
*     *             			*
*     *           psi_get_version		*
*     *             			*
*     ***********************************

      integer function psi_get_version()
      implicit none 


#ifdef MPI
      include 'mpif.h'
      integer mpierr
#endif
#ifdef TCGMSG
#include "tcgmsg.fh"
#include "msgtypesf.h"
#endif

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
#ifdef MPI
      msglen = 1
      call MPI_BCAST(version,msglen,MPI_INTEGER,
     >               MASTER,MPI_COMM_WORLD,mpierr)
#endif
#ifdef TCGMSG
      msglen = 1
      call BRDCST(9+MSGINT,version,mitob(msglen),MASTER)
#endif

      psi_get_version = version
      return
      end
