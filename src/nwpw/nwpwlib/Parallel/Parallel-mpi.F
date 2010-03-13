*
* $Id: Parallel-mpi.F,v 1.14 2008-11-17 17:25:42 bylaska Exp $
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

#ifdef MPI4
#include "stupid_mpi4.fh"
      stupid_world   = MPI_COMM_WORLD
      stupid_sum     = MPI_SUM
      stupid_max     = MPI_MAX
      stupid_min     = MPI_MIN
      stupid_integer = MPI_INTEGER
      stupid_double  = MPI_DOUBLE_PRECISION
      stupid_complex = MPI_DOUBLE_COMPLEX
      stupid_comm_i   = MPI_COMM_WORLD
      stupid_comm_j   = -99
      stupid_comm_k   = -99
#endif

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
      comm_i   = MPI_COMM_WORLD
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
*     where np_i = nrows, and np_j = np/np_i
*
      subroutine Parallel2d_Init(ncolumns)
      implicit none
      integer ncolumns

#include "Parallel.fh"
#include "mafdecls.fh"
#include "errquit.fh"

#ifdef MPI4
#include "stupid_mpi4.fh"
#else
#include "mpif.h"
#endif

*     *** local variables ***
      integer i,j,icount,ierr
      integer tmp(2),tmp2(2),mpi_group

      if (ncolumns.gt.1) then

      np_i = np/ncolumns
      np_j = ncolumns


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

#ifdef MPI4
      if (.not.MA_push_get(mt_int,np,'tmppp3',tmp2(2),tmp2(1)))
     >  call errquit('Parallel2d_init:out of stack memory',0, MA_ERR)
      stupid_np_i = np_i
      stupid_np_j = np_j
*     **** set global processor group ****
      call MPI_COMM_group(stupid_world,stupid_mpi_group,stupid_ierr)

      do i=0,np_i-1
        int_mb(tmp(1)+i) = int_mb(procNd(1) + i + taskid_j*np_i) 
      end do
      call stupid_icopy8to4(np_i,int_mb(tmp(1)),int_mb(tmp2(1)))
      call MPI_Group_incl(stupid_mpi_group,
     >                    stupid_np_i,
     >                    int_mb(tmp2(1)),
     >                    stupid_group_i,stupid_ierr)
      call MPI_Comm_create(stupid_world,stupid_group_i,
     >                     stupid_comm_i,stupid_ierr)

      do j=0,np_j-1
        int_mb(tmp(1)+j) = int_mb(procNd(1) + taskid_i + j*np_i) 
      end do
      call stupid_icopy8to4(np_j,int_mb(tmp(1)),int_mb(tmp2(1)))
      call MPI_Group_incl(stupid_mpi_group,stupid_np_j,
     >                    int_mb(tmp2(1)),
     >                    stupid_group_j,stupid_ierr)
      call MPI_Comm_create(stupid_world,stupid_group_j,stupid_comm_j,  
     >                     stupid_ierr)


      group_i = stupid_group_i
      comm_i  = stupid_comm_i
      group_j = stupid_group_j
      comm_j  = stupid_comm_j
      if (.not.MA_pop_stack(tmp2(2)))
     >  call errquit('Parallel2d_init:popping stack memory',0, MA_ERR)

#else
*     **** set global processor group ****
      call MPI_COMM_group(MPI_COMM_WORLD,mpi_group,ierr)

      do i=0,np_i-1
        int_mb(tmp(1)+i) = int_mb(procNd(1) + i + taskid_j*np_i) 
      end do
      call MPI_Group_incl(mpi_group,np_i,int_mb(tmp(1)),group_i,ierr)
      call MPI_Comm_create(MPI_COMM_WORLD,group_i,comm_i,  ierr)

      do j=0,np_j-1
        int_mb(tmp(1)+j) = int_mb(procNd(1) + taskid_i + j*np_i) 
      end do
      call MPI_Group_incl(mpi_group,np_j,int_mb(tmp(1)),group_j,ierr)
      call MPI_Comm_create(MPI_COMM_WORLD,group_j,comm_j,  ierr)
#endif

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

#ifdef MPI4
#include "stupid_mpi4.fh"
#else
#include "mpif.h"
#endif

#ifdef MPI4
      if (np_j.gt.1) then
*      **** free comm_i and comm_j communicators ****
      call MPI_Comm_free(stupid_comm_i,  stupid_ierr)
      call MPI_Group_free(stupid_group_i,stupid_ierr)
      call MPI_Comm_free(stupid_comm_j,  stupid_ierr)
      call MPI_Group_free(stupid_group_j,stupid_ierr)
      end if
#else
*     *** local variable ***
      integer mpierr

      if (np_j.gt.1) then
*      **** free comm_i and comm_j communicators ****
      call MPI_Comm_free(comm_i,  mpierr)
      call MPI_Group_free(group_i,mpierr)
      call MPI_Comm_free(comm_j,  mpierr)
      call MPI_Group_free(group_j,mpierr)
      end if
#endif

      return
      end



*     *************************************
*     *                                   *
*     *        Parallel3d_Init            *
*     *                                   *
*     *************************************

*     Sset up the 3d processor grid = np_i x np_j x np_k, 
*     where np_i = nrows=np/(np_j*np_k), 
*     np_j = ncolumns, and np_k = nzones
*
      subroutine Parallel3d_Init(ncolumns,nzones)
      implicit none
      integer ncolumns,nzones

#include "Parallel.fh"
#include "mafdecls.fh"
#include "errquit.fh"

#ifdef MPI4
#include "stupid_mpi4.fh"
#else
#include "mpif.h"
#endif

*     *** local variables ***
      integer i,j,k,icount,ierr
      integer tmp(2),tmp2(2),mpi_group

      np_i = np/(ncolumns*nzones)
      np_j = ncolumns
      np_k = nzones

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

#ifdef MPI4
      if (.not.MA_push_get(mt_int,np,'tmppp3',tmp2(2),tmp2(1)))
     >  call errquit('Parallel3d_init:out of stack memory',0, MA_ERR)
      stupid_np_i = np_i
      stupid_np_j = np_j
      stupid_np_k = np_k
*     **** set global processor group ****
      call MPI_COMM_group(stupid_world,stupid_mpi_group,stupid_ierr)

      do i=0,np_i-1
        int_mb(tmp(1)+i) = int_mb(procNd(1) 
     >                           + i 
     >                           + taskid_j*np_i
     >                           + taskid_k*np_i*np_j) 
      end do
      call stupid_icopy8to4(np_i,int_mb(tmp(1)),int_mb(tmp2(1)))
      call MPI_Group_incl(stupid_mpi_group,
     >                    stupid_np_i,
     >                    int_mb(tmp2(1)),
     >                    stupid_group_i,stupid_ierr)
      call MPI_Comm_create(stupid_world,stupid_group_i,
     >                     stupid_comm_i,stupid_ierr)

      do j=0,np_j-1
        int_mb(tmp(1)+j) = int_mb(procNd(1) 
     >                           + taskid_i 
     >                           + j*np_i
     >                           + taskid_k*np_i*np_j) 
      end do
      call stupid_icopy8to4(np_j,int_mb(tmp(1)),int_mb(tmp2(1)))
      call MPI_Group_incl(stupid_mpi_group,stupid_np_j,
     >                    int_mb(tmp2(1)),
     >                    stupid_group_j,stupid_ierr)
      call MPI_Comm_create(stupid_world,stupid_group_j,stupid_comm_j,  
     >                     stupid_ierr)

      do k=0,np_k-1
        int_mb(tmp(1)+k) = int_mb(procNd(1)
     >                           + taskid_i
     >                           + taskid_j*np_i
     >                           + k*np_i*np_j)
      end do
      call stupid_icopy8to4(np_k,int_mb(tmp(1)),int_mb(tmp2(1)))
      call MPI_Group_incl(stupid_mpi_group,
     >                    stupid_np_k,
     >                    int_mb(tmp2(1)),
     >                    stupid_group_k,stupid_ierr)
      call MPI_Comm_create(stupid_world,stupid_group_k,
     >                     stupid_comm_k,stupid_ierr)

      group_i = stupid_group_i
      comm_i  = stupid_comm_i
      group_j = stupid_group_j
      comm_j  = stupid_comm_j
      group_k = stupid_group_k
      comm_k  = stupid_comm_k
      if (.not.MA_pop_stack(tmp2(2)))
     >  call errquit('Parallel3d_init:popping stack memory',0, MA_ERR)

#else
*     **** set global processor group ****
      call MPI_COMM_group(MPI_COMM_WORLD,mpi_group,ierr)

      do i=0,np_i-1
        int_mb(tmp(1)+i) = int_mb(procNd(1) 
     >                           + i 
     >                           + taskid_j*np_i 
     >                           + taskid_k*np_i*np_j) 
      end do
      call MPI_Group_incl(mpi_group,np_i,int_mb(tmp(1)),group_i,ierr)
      call MPI_Comm_create(MPI_COMM_WORLD,group_i,comm_i,ierr)

      do j=0,np_j-1
        int_mb(tmp(1)+j) = int_mb(procNd(1) 
     >                           + taskid_i 
     >                           + j*np_i
     >                           + taskid_k*np_i*np_j) 
      end do
      call MPI_Group_incl(mpi_group,np_j,int_mb(tmp(1)),group_j,ierr)
      call MPI_Comm_create(MPI_COMM_WORLD,group_j,comm_j,ierr)

      do k=0,np_k-1
        int_mb(tmp(1)+k) = int_mb(procNd(1) 
     >                           + taskid_i 
     >                           + taskid_j*np_i
     >                           + k*np_i*np_j) 
      end do
      call MPI_Group_incl(mpi_group,np_k,int_mb(tmp(1)),group_k,ierr)
      call MPI_Comm_create(MPI_COMM_WORLD,group_k,comm_k,ierr)
#endif

      if (.not.MA_pop_stack(tmp(2)))
     >  call errquit('Parallel3d_init:popping stack memory',0, MA_ERR)
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

#ifdef MPI4
#include "stupid_mpi4.fh"
#else
#include "mpif.h"
#endif

#ifdef MPI4
*      **** free comm_i and comm_j communicators ****
      call MPI_Comm_free(stupid_comm_i,  stupid_ierr)
      call MPI_Group_free(stupid_group_i,stupid_ierr)
      call MPI_Comm_free(stupid_comm_j,  stupid_ierr)
      call MPI_Group_free(stupid_group_j,stupid_ierr)
      call MPI_Comm_free(stupid_comm_k,  stupid_ierr)
      call MPI_Group_free(stupid_group_k,stupid_ierr)
#else

*     *** local variable ***
      integer mpierr

*      **** free comm_i and comm_j communicators ****
      call MPI_Comm_free(comm_i,  mpierr)
      call MPI_Group_free(group_i,mpierr)
      call MPI_Comm_free(comm_j,  mpierr)
      call MPI_Group_free(group_j,mpierr)
      call MPI_Comm_free(comm_k,  mpierr)
      call MPI_Group_free(group_k,mpierr)
#endif

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

#include "Parallel.fh"

#ifdef MPI4
#include "stupid_mpi4.fh"
#else
#include "mpif.h"
#endif
      real*8 sumall
#ifdef MPI4
      if (np.gt.1) then
         stupid_msglen = 1
         call MPI_Allreduce(sum,sumall,stupid_msglen,stupid_double,
     >                      stupid_max,stupid_world,stupid_ierr)
         sum = sumall
      end if
#else
      integer msglen,mpierr
      if (np.gt.1) then
         msglen = 1
         call MPI_Allreduce(sum,sumall,msglen,MPI_DOUBLE_PRECISION,
     >                       MPI_MAX,MPI_COMM_WORLD,mpierr)
         sum = sumall
      end if
#endif
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


#include "Parallel.fh"

#ifdef MPI4
#include "stupid_mpi4.fh"
#else
#include "mpif.h"
#endif

      real*8 sumall

#ifdef MPI4
      if (np.gt.1) then
         stupid_msglen = 1
         call MPI_Allreduce(sum,sumall,stupid_msglen,stupid_double,
     >                      stupid_sum,stupid_world,stupid_ierr)
         sum = sumall
      end if
#else
      integer msglen,mpierr

      if (np.gt.1) then
         msglen = 1
         call MPI_Allreduce(sum,sumall,msglen,MPI_DOUBLE_PRECISION,
     >                       MPI_SUM,MPI_COMM_WORLD,mpierr)
         sum = sumall
      end if
#endif

      return
      end


*     ***********************************
*     *                                 *
*     *         Parallel_ISumAll        *
*     *                                 *
*     ***********************************

      subroutine Parallel_ISumAll(sum)
c     implicit none
      integer sum


#include "Parallel.fh"

#ifdef MPI4
#include "stupid_mpi4.fh"
#else
#include "mpif.h"
#endif

#ifdef MPI4
      integer*4 tsum,sumall

      if (np.gt.1) then
         stupid_msglen = 1
         tsum = sum
         call MPI_Allreduce(tsum,sumall,stupid_msglen,stupid_integer,
     >                       stupid_sum,stupid_world,stupid_ierr)
         sum = sumall
      end if
#else

      integer msglen,mpierr
      integer sumall

      if (np.gt.1) then
         msglen = 1

         call MPI_Allreduce(sum,sumall,msglen,MPI_INTEGER,
     >                       MPI_SUM,MPI_COMM_WORLD,mpierr)
         sum = sumall
      end if
#endif

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
#include "errquit.fh"
#include "Parallel.fh"

#ifdef MPI4
#include "stupid_mpi4.fh"
#else
#include "mpif.h"
#endif

*     **** local variable ****
      logical value
      integer msglen,mpierr

*     **** temporary workspace ****
      integer sumall(2)


      call nwpw_timing_start(2)
      if (np.gt.1) then

*     ***** allocate temporary space ****
      value = MA_push_get(mt_dbl,n,'sumall',sumall(2),sumall(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)

#ifdef MPI4
      stupid_msglen = n
      call MPI_Allreduce(sum,dbl_mb(sumall(1)),stupid_msglen,
     >                stupid_double,
     >                stupid_sum,stupid_world,stupid_ierr)
#else
      msglen = n
      call MPI_Allreduce(sum,dbl_mb(sumall(1)),msglen,
     >                MPI_DOUBLE_PRECISION,
     >                MPI_SUM,MPI_COMM_WORLD,mpierr)
#endif

      call dcopy(n,dbl_mb(sumall(1)),1,sum,1)
      value = MA_pop_stack(sumall(2))

      end if
      call nwpw_timing_end(2)
      return
      end


*     ***********************************
*     *                                 *
*     *      Parallel_Vector_ISumAll    *
*     *                                 *
*     ***********************************

      subroutine Parallel_Vector_ISumAll(n,sum)
c     implicit none
      integer n
      integer  sum(*)

#include "mafdecls.fh"
#include "errquit.fh"
#include "Parallel.fh"

#ifdef MPI4
#include "stupid_mpi4.fh"
#else
#include "mpif.h"
#endif

      logical value
      integer msglen,mpierr

*     **** temporary workspace ****
      integer sumall(2),sumall1(2)


      call nwpw_timing_start(2)

      if (np.gt.1) then

*     ***** allocate temporary space ****
      value = MA_push_get(mt_int,n,'sumall',sumall(2),sumall(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)

#ifdef MPI4
      if (.not.MA_push_get(mt_int,n,'sumall1',sumall1(2),sumall1(1)))
     > call errquit('out of stack memory',0, MA_ERR)

      stupid_msglen = n
      call stupid_icopy8to4(n,sum,int_mb(sumall1(1)))
      call MPI_Allreduce(int_mb(sumall1(1)),
     >                   int_mb(sumall(1)),stupid_msglen,
     >                stupid_integer,
     >                stupid_sum,stupid_world,stupid_ierr)
      call stupid_icopy4to8(n,int_mb(sumall(1)),sum)
      if (.not.MA_pop_stack(sumall1(2)))
     >   call errquit('error popping stack',0,MA_ERR)
#else
      msglen = n
      call MPI_Allreduce(sum,int_mb(sumall(1)),msglen,
     >                MPI_INTEGER,
     >                MPI_SUM,MPI_COMM_WORLD,mpierr)
      call icopy(n,int_mb(sumall(1)),1,sum,1)
#endif
      value = MA_pop_stack(sumall(2))
      if (.not. value) call errquit('error popping stack',0, MA_ERR)

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

#include "Parallel.fh"
#ifdef MPI4
#include "stupid_mpi4.fh"
#else
#include "mpif.h"
#endif


#ifdef MPI4
      integer*4 tpsend

      if (np.gt.1) then
         stupid_msglen = 1
         tpsend        = psend
         call MPI_Bcast(sum,stupid_msglen,stupid_double,
     >                  tpsend,stupid_world,stupid_ierr)
      end if
#else
      integer ierr

      if (np.gt.1) then
         call MPI_Bcast(sum,1,MPI_DOUBLE_PRECISION,
     >                  psend,MPI_COMM_WORLD,ierr)
      end if
#endif

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

#include "Parallel.fh"
#ifdef MPI4
#include "stupid_mpi4.fh"
#else
#include "mpif.h"
#endif


#ifdef MPI4
      integer*4 tpsend

      if (np.gt.1) then
         stupid_msglen = nsize
         tpsend        = psend
         call MPI_Bcast(sum,stupid_msglen,stupid_double,
     >                  tpsend,stupid_world,stupid_ierr)
      end if
#else
      integer ierr

      if (np.gt.1) then
         call MPI_Bcast(sum,nsize,MPI_DOUBLE_PRECISION,
     >                  psend,MPI_COMM_WORLD,ierr)
      end if
#endif

      return
      end




*     ***********************************
*     *                                 *
*     *      Parallel_Brdcst_ivalues    *
*     *                                 *
*     ***********************************

      subroutine Parallel_Brdcst_ivalues(psend,nsize,isum)
      implicit none
      integer psend,nsize
      integer isum(*)

#include "mafdecls.fh"
#include "errquit.fh"
#include "Parallel.fh"

#ifdef MPI4
#include "stupid_mpi4.fh"
#else
#include "mpif.h"
#endif


#ifdef MPI4
      logical   value
      integer*4 tpsend
      integer   sumall(2)

      if (np.gt.1) then
         stupid_msglen = nsize
         tpsend        = psend
         value = MA_push_get(mt_int,nsize,'sumall',sumall(2),sumall(1))
         if (.not. value) call errquit('out of stack memory',0, MA_ERR)

         call stupid_icopy8to4(nsize,isum,int_mb(sumall(1)))
         call MPI_Bcast(int_mb(sumall(1)),stupid_msglen,stupid_integer,
     >                  tpsend,stupid_world,stupid_ierr)
         call stupid_icopy4to8(nsize,int_mb(sumall(1)),isum)

         value = MA_pop_stack(sumall(2))
         if (.not. value) call errquit('error popping stack',0, MA_ERR)
      end if
#else
      integer ierr

      if (np.gt.1) then
         call MPI_Bcast(isum,nsize,MPI_INTEGER,
     >                  psend,MPI_COMM_WORLD,ierr)
      end if
#endif

      return
      end




#ifdef MPI4
      subroutine stupid_icopy8to4(n,array8,array4)
      implicit none
      integer n
      integer*8 array8(*)
      integer*4 array4(*)
      integer i
      do i=1,n
          array4(i) = array8(i)
      end do
      return
      end

      subroutine stupid_icopy4to8(n,array4,array8)
      implicit none
      integer n
      integer*4 array4(*)
      integer*8 array8(*)
      integer i
      do i=1,n
          array8(i) = array4(i)
      end do
      return
      end
#endif


*     ***********************************
*     *                                 *
*     *      Parallel_mpiWaitAll        *
*     *                                 *
*     ***********************************

      subroutine Parallel_mpiWaitAll(nreq,req)
      implicit none
      integer nreq,req(*)

#include "mafdecls.fh"
#include "errquit.fh"
#include "mpif.h"
#include "Parallel.fh"

*     *** local variables ***
#ifdef MPI4
#include "stupid_mpi4.fh"
      integer status(2),request2(2)
#else
      integer status(2),mpierr
#endif
     
      if (nreq.gt.0) then

*     **** allocate status memory ****
      if (.not.MA_push_get(mt_int,MPI_STATUS_SIZE*nreq*2,
     >                     'status',status(2),status(1)))
     > call errquit('Parallel_mpiWaitAll:out of stack',0,MA_ERR)



*     **** wait for completion of mp_send, also do a sync ****

#ifdef MPI4

      if (.not.MA_push_get(mt_int,nreq,'rqst2',request2(2),request2(1)))
     > call errquit(' Parallel_mpiWaitAll:out of stack',1,MA_ERR)
      call stupid_icopy8to4(nreq,req,int_mb(request2(1)))
      stupid_msglen = nreq
      call MPI_WAITALL(stupid_msglen,
     >                 int_mb(request2(1)),
     >                 int_mb(status(1)),stupid_ierr)
      if (.not.MA_pop_stack(request2(2)))
     > call errquit(' Parallel_mpiWaitAll:popping stack',1,MA_ERR)

c     *** danger! ***
      if (stupid_ierr.ne.0) then
         write(*,*) "MPI_WAITALL Failed on taskid=",taskid,
     >              "  mpierr=",stupid_ierr
         call util_flush(6)
      end if

#else
      call MPI_WAITALL(nreq,req,int_mb(status(1)),mpierr)
#endif

*     *** may need to check status here??? ***

      !*** deallocate status memory ***
      if (.not.MA_pop_stack(status(2)))
     > call errquit(' Parallel_mpiWaitAll:popping stack',0,MA_ERR)

      end if
      return
      end

