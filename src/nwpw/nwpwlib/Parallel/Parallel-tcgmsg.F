*
* $Id: Parallel-tcgmsg.F,v 1.10 2008-09-22 22:53:47 bylaska Exp $
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


*     **** set up 3d processor grid = np x 1 x 1****
      if (.not.MA_alloc_get(mt_int,np,'procNd',procNd(2),procNd(1)))
     >  call errquit('Parallel_init:out of heap memory',0, MA_ERR)

      np_i = np
      np_j = 1
      np_k = 1
      do i=0,np-1
        int_mb(procNd(1)+i) = i
      end do
      taskid_i = taskid
      taskid_j = 0
      taskid_k = 0
      comm_i   = ga_pgroup_get_world()
      comm_j   = -99 
      comm_k   = -99 

      return 
      end


*     *************************************
*     *                                   *
*     *        Parallel2d_Init            *
*     *                                   *
*     *************************************

*     Sset up the 2d processor grid = np_i x np_j, 
*     where np_i = ncolumns, and np_j = np/np_i
*
      subroutine Parallel2d_Init(ncolumns)
      implicit none
      integer ncolumns

#include "Parallel.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "global.fh"

*     *** local variables ***
      integer i,j,icount
      integer tmp(2)

      np_i = np/ncolumns
      np_j = ncolumns

      if (np_j.gt.1) then

      icount = 0
      do j=0,np_j-1
      do i=0,np_i-1
        if (icount.eq.taskid) then
           taskid_i = i
           taskid_j = j
        end if
        int_mb(procNd(1) + i + j*np_i) = icount
        icount = mod((icount+1),np)
      end do
      end do

      if (.not.MA_push_get(mt_int,np,'tmppp2',tmp(2),tmp(1)))
     >  call errquit('Parallel2d_init:out of stack memory',0, MA_ERR)
        do i=0,np_i-1
          int_mb(tmp(1)+i) = int_mb(procNd(1) + i + taskid_j*np_i) 
        end do
        comm_i = ga_pgroup_create(int_mb(tmp(1)),np_i)
        do j=0,np_j-1
          int_mb(tmp(1)+j) = int_mb(procNd(1) + taskid_i + j*np_i) 
        end do
        comm_j = ga_pgroup_create(int_mb(tmp(1)),np_j)
      if (.not.MA_pop_stack(tmp(2)))
     >  call errquit('Parallel2d_init:popping stack memory',0, MA_ERR)

      end if

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

       if (np_j.gt.1) then

*      **** free comm_i and comm_j communicators ****
       if (.not.ga_pgroup_destroy(comm_i))
     >  call errquit('Parallel2d_Finalize: error destoying comm_i',0,0)
       if (.not.ga_pgroup_destroy(comm_j))
     >  call errquit('Parallel2d_Finalize: error destoying comm_j',0,1)

      end if
      return
      end


*     *************************************
*     *                                   *
*     *        Parallel3d_Init            *
*     *                                   *
*     *************************************

*     Sset up the 3d processor grid = np_i x np_j x np_k, 
*     where np_i = np/(np_j*np_k), and np_j = ncolumns and np_k=nzones
*
      subroutine Parallel3d_Init(ncolumns,nzones)
      implicit none
      integer ncolumns,nzones

#include "Parallel.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "global.fh"

*     *** local variables ***
      integer i,j,k,icount
      integer tmp(2)

      np_i = np/(ncolumns*nzones)
      np_j = ncolumns
      np_k = nzones

      if ((np_j.gt.1).or.(np_k.gt.1)) then

      icount = 0
      do k=0,np_k-1
      do j=0,np_j-1
      do i=0,np_i-1
        if (icount.eq.taskid) then
           taskid_i = i
           taskid_j = j
           taskid_k = k
        end if
        int_mb(procNd(1) + i + j*np_i + k*np_i*np_j) = icount
        icount = mod((icount+1),np)
      end do
      end do
      end do

      if (.not.MA_push_get(mt_int,np,'tmppp2',tmp(2),tmp(1)))
     >  call errquit('Parallel3d_init:out of stack memory',0, MA_ERR)
        do i=0,np_i-1
          int_mb(tmp(1)+i) = int_mb(procNd(1) 
     >                             + i 
     >                             + taskid_j*np_i
     >                             + taskid_k*np_i*np_j) 
        end do
        comm_i = ga_pgroup_create(int_mb(tmp(1)),np_i)
        do j=0,np_j-1
          int_mb(tmp(1)+j) = int_mb(procNd(1)
     >                              + taskid_i 
     >                              + j*np_i
     >                              + taskid_k*np_i*np_j) 
        end do
        comm_j = ga_pgroup_create(int_mb(tmp(1)),np_j)
        do k=0,np_k-1
          int_mb(tmp(1)+k) = int_mb(procNd(1)
     >                              + taskid_i 
     >                              + taskid_j*np_i
     >                              + k*np_i*np_j) 
        end do
        comm_k = ga_pgroup_create(int_mb(tmp(1)),np_k)
      if (.not.MA_pop_stack(tmp(2)))
     >  call errquit('Parallel3d_init:popping stack memory',0, MA_ERR)

      end if

      return
      end


*     *************************************
*     *                                   *
*     *        Parallel3d_Finalize        *
*     *                                   *
*     *************************************

      subroutine Parallel3d_Finalize()
      implicit none

#include "Parallel.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "global.fh"

      if ((np_j.gt.1).or.(np_k.gt.1)) then

*      **** free comm_i and comm_j communicators ****
       if (.not.ga_pgroup_destroy(comm_i))
     >  call errquit('Parallel3d_Finalize: error destoying comm_i',0,0)
       if (.not.ga_pgroup_destroy(comm_j))
     >  call errquit('Parallel3d_Finalize: error destoying comm_j',0,1)
       if (.not.ga_pgroup_destroy(comm_k))
     >  call errquit('Parallel3d_Finalize: error destoying comm_k',0,2)

      end if
      return
      end



*     ***********************************
*     *                                 *
*     *         Parallel_MaxAll         *
*     *                                 *
*     ***********************************

      subroutine Parallel_MaxAll(sum)
c     implicit none
      real*8  sum

#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "Parallel.fh"

      if (np.gt.1) then
         call GA_DGOP(9+MSGDBL,sum,1,'max')
      end if

      return
      end




*     ***********************************
*     *                                 *
*     *         Parallel_SumAll         *
*     *                                 *
*     ***********************************

      subroutine Parallel_SumAll(sum)
c     implicit none
      real*8  sum

#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "Parallel.fh"

      if (np.gt.1) then
         call GA_DGOP(9+MSGDBL,sum,1,'+')
      end if

      return
      end



*     ***********************************
*     *                                 *
*     *         Parallel_ISumAll        *
*     *                                 *
*     ***********************************

      subroutine Parallel_ISumAll(sum)
c     implicit none
      integer  sum

#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "Parallel.fh"

      if (np.gt.1) then
         call GA_IGOP(9+MSGINT,sum,1,'+')
      end if
      
      return
      end


*     ***********************************
*     *                                 *
*     *      Parallel_Vector_SumAll     *
*     *                                 *
*     ***********************************
      subroutine Parallel_Vector_SumAll(n,sum)
c     implicit none
      integer n
      real*8  sum(*)

#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "errquit.fh"
#include "Parallel.fh"


      call nwpw_timing_start(2)
      if (np.gt.1) then
         call GA_DGOP(9+MSGDBL,sum,n,'+')
      end if
      call nwpw_timing_end(2)
      return
      end



*     ***********************************
*     *                                 *
*     *     Parallel_Vector_ISumAll     *
*     *                                 *
*     ***********************************

      subroutine Parallel_Vector_ISumAll(n,sum)
c     implicit none
      integer n
      integer  sum(*)

#include "mafdecls.fh"
#include "errquit.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "Parallel.fh"


      call nwpw_timing_start(2)
      if (np.gt.1) then
         call GA_IGOP(9+MSGINT,sum,n,'+')
      end if
      call nwpw_timing_end(2)

      return
      end



*     ***********************************
*     *                                 *
*     *      Parallel_Brdcst_value      *
*     *                                 *
*     ***********************************

      subroutine Parallel_Brdcst_value(psend,sum)
      implicit none
      integer psend
      real*8  sum

#include "mafdecls.fh"
#include "errquit.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "Parallel.fh"

*     **** local variables ****
      integer msglen

      if (np.gt.1) then
         msglen = 1
         call BRDCST(9+MSGDBL,sum,mdtob(msglen),psend)
      end if

      return
      end


*     ***********************************
*     *                                 *
*     *      Parallel_Brdcst_values     *
*     *                                 *
*     ***********************************

      subroutine Parallel_Brdcst_values(psend,nsize,sum)
      implicit none
      integer psend,nsize
      real*8  sum(*)

#include "mafdecls.fh"
#include "errquit.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "Parallel.fh"


      if (np.gt.1) then
         call BRDCST(9+MSGDBL,sum,mdtob(nsize),psend)
      end if

      return
      end


*     ***********************************
*     *                                 *
*     *      Parallel_Brdcst_ivalues     *
*     *                                 *
*     ***********************************

      subroutine Parallel_Brdcst_ivalues(psend,nsize,isum)
      implicit none
      integer psend,nsize
      integer  isum(*)

#include "mafdecls.fh"
#include "errquit.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "Parallel.fh"

      if (np.gt.1) then
         call BRDCST(9+MSGINT,isum,mitob(nsize),psend)
      end if
      return
      end


