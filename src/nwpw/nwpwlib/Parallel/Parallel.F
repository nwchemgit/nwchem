#define TCGMSG
* Parallel.f
* Author - Eric Bylaska
*
*   These routines are to be used to keep track of the parallel message
* passing variables, as well as iniitialize and deinitialize the
* message passing routines.
*

      subroutine Parallel_Init()
*
* $Id: Parallel.F,v 1.1 2001-08-30 16:58:34 bylaska Exp $
*
      implicit none

      integer taskid,np
      common / parallel / taskid,np

#ifdef MPI
      include 'mpif.h'
      integer mpierr
#endif
#ifdef TCGMSG
#include "tcgmsg.fh"
#endif

#ifdef MPI
*     **** MPI initiializer *****
c     call MPI_INIT(mpierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,taskid,mpierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,np,mpierr)
#endif
#ifdef TCGMSG
      np     = nnodes()
      taskid = nodeid()
#endif



      return 
      end

      subroutine Parallel_Finalize()
      implicit none

      integer taskid,np
      common / parallel / taskid,np

#ifdef MPI
      include 'mpif.h'
      integer mpierr
#endif

c     call MPI_FINALIZE(mpierr) 
      return
      end


      subroutine Parallel_np(np_out)
      implicit none
      integer np_out

      integer taskid,np
      common / parallel / taskid,np

      np_out = np
      return
      end



      subroutine Parallel_taskid(task_out)
      implicit none
      integer task_out
      
      integer taskid,np
      common / parallel / taskid,np

      task_out = taskid
      return 
      end
