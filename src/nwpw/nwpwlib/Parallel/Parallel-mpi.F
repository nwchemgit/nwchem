*
* $Id: Parallel-mpi.F,v 1.3 2006-08-16 23:53:27 bylaska Exp $
*

* Parallel.f
* Author - Eric Bylaska
*
*   These routines are to be used to keep track of the parallel message
* passing variables, as well as iniitialize and deinitialize the
* message passing routines.
*


*     *************************************
*     *                                   *
*     *        Parallel_Init              *
*     *                                   *
*     *************************************

      subroutine Parallel_Init()
      implicit none

#include "Parallel.fh"
#include "mafdecls.fh"
#include "errquit.fh"

#include "tcgmsg.fh"
#include "global.fh"

#include "mpif.h"


*     **** local variables ****
      integer i

c*     **** MPI initiializer *****
cc     call MPI_INIT(mpierr)
c      call MPI_COMM_RANK(MPI_COMM_WORLD,taskid,mpierr)
c      call MPI_COMM_SIZE(MPI_COMM_WORLD,np,mpierr)

      np     = nnodes()
      taskid = nodeid()


*     **** set up 2d processor grid = np x 1****
      if (.not.MA_alloc_get(mt_int,np,'proc2d',proc2d(2),proc2d(1)))
     >  call errquit('Parallel_init:out of heap memory',0, MA_ERR)

      np_i = np
      np_j = 1
      do i=0,np-1
        int_mb(proc2d(1)+i) = i
      end do
      taskid_i = taskid
      taskid_j = 0
      comm_i   = MPI_COMM_WORLD
      comm_j   = -99 

      return 
      end



*     *************************************
*     *                                   *
*     *        Parallel2d_Init            *
*     *                                   *
*     *************************************

*     Sset up the 2d processor grid = np_i x np_j, 
*     where np_i = nrows, and np_j = np/np_i
*
      subroutine Parallel2d_Init(nrows)
      implicit none
      integer nrows

#include "Parallel.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "mpif.h"

*     *** local variables ***
      integer i,j,icount,ierr
      integer tmp(2),n,mpi_group

      np_i = nrows
      np_j = np/np_i


      icount = 0
      do j=0,np_j-1
      do i=0,np_i-1
        if (icount.eq.taskid) then
           taskid_i = i
           taskid_j = j
        end if
        int_mb(proc2d(1) + i + j*np_i) = icount
        icount = mod((icount+1),np)
      end do
      end do

      if (.not.MA_push_get(mt_int,np,'tmppp2',tmp(2),tmp(1)))
     >  call errquit('Parallel2d_init:out of stack memory',0, MA_ERR)

*     **** set global processor group ****
      call MPI_COMM_group(MPI_COMM_WORLD,mpi_group,ierr)

      do i=0,np_i-1
        int_mb(tmp(1)+i) = int_mb(proc2d(1) + i + taskid_j*np_i) 
      end do
      call MPI_Group_incl(mpi_group,np_i,int_mb(tmp(1)),group_i,ierr)
      call MPI_Comm_create(MPI_COMM_WORLD,group_i,comm_i,  ierr)


      do j=0,np_j-1
        int_mb(tmp(1)+j) = int_mb(proc2d(1) + taskid_i + j*np_i) 
      end do
      call MPI_Group_incl(mpi_group,np_j,int_mb(tmp(1)),group_j,ierr)
      call MPI_Comm_create(MPI_COMM_WORLD,group_j,comm_j,  ierr)

      if (.not.MA_pop_stack(tmp(2)))
     >  call errquit('Parallel2d_init:popping stack memory',0, MA_ERR)
      return
      end


*     *************************************
*     *                                   *
*     *        Parallel2d_Finalize        *
*     *                                   *
*     *************************************

      subroutine Parallel2d_Finalize()
      implicit none

#include "Parallel.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "mpif.h"


*     *** local variable ***
      integer mpierr

*      **** free comm_i and comm_j communicators ****
      call MPI_Comm_free(comm_i,  mpierr)
      call MPI_Group_free(group_i,mpierr)
      call MPI_Comm_free(comm_j,  mpierr)
      call MPI_Group_free(group_j,mpierr)

      return
      end


