*
* $Id: Parallel-tcgmsg.F,v 1.1 2006-08-13 01:03:25 bylaska Exp $
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

c      include 'mpif.h'
c      integer mpierr

*     **** local variables ****
      integer i

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
      comm_i   = ga_pgroup_get_world()
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
#include "global.fh"

*     *** local variables ***
      integer i,j,icount
      integer tmp(2)

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
        do i=0,np_i-1
          int_mb(tmp(1)+i) = int_mb(proc2d(1) + i + taskid_j*np_i) 
        end do
        comm_i = ga_pgroup_create(int_mb(tmp(1)),np_i)
        do j=0,np_j-1
          int_mb(tmp(1)+j) = int_mb(proc2d(1) + taskid_i + j*np_i) 
        end do
        comm_j = ga_pgroup_create(int_mb(tmp(1)),np_j)
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
#include "global.fh"

*      **** free comm_i and comm_j communicators ****
       if (.not.ga_pgroup_destroy(comm_i))
     >  call errquit('Parallel2d_Finalize: error destoying comm_i',0,0)
       if (.not.ga_pgroup_destroy(comm_j))
     >  call errquit('Parallel2d_Finalize: error destoying comm_j',0,1)

      return
      end


