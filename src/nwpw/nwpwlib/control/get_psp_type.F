*
* $Id: get_psp_type.F,v 1.1 2003-12-02 19:46:09 bylaska Exp $
*

#define TCGMSG
*     ***********************************
*     *             			*
*     *           get_psp_type		*
*     *             			*
*     ***********************************

      subroutine get_psp_type(psp_filename,psp_type)
      implicit none 
      character*(*) psp_filename
      integer psp_type


#ifdef MPI
      include 'mpif.h'
      integer mpierr
#endif
#ifdef TCGMSG
#include "tcgmsg.fh"
#include "msgtypesf.h"
#endif

*    *** local variables ***
      integer MASTER,taskid
      parameter(MASTER=0)

      integer l,msglen
      character*255 full_filename


      call Parallel_taskid(taskid)

*     **** open psp_filename file ****
      if (taskid.eq.MASTER) then
        call util_file_name_noprefix(psp_filename,.false.,.false.,
     >                    full_filename)
        l = index(full_filename,' ') - 1
        open(unit=11,file=full_filename(1:l),
     >             status='old',form='formatted')
        read(11,'(I)',ERR=9998,END=9998) psp_type
        go to 9999
 9998   psp_type = 0
 9999   continue

        close(11)
      end if


c     **** send header to all nodes ****
#ifdef MPI
      msglen = 1
      call MPI_BCAST(psp_type,msglen,MPI_INTEGER,
     >               MASTER,MPI_COMM_WORLD,mpierr)
#endif

#ifdef TCGMSG
      msglen = 1
      call BRDCST(9+MSGINT,psp_type,mitob(msglen),MASTER)
#endif

      return
      end

