#define GA_DEV 1
#define USE_SUBGROUPS 1
*
* $Id$
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
#include "bafdecls.fh"
#include "errquit.fh"

#include "tcgmsg.fh"
#include "global.fh"
#include "mpif.h"

#ifdef GA_DEV


ccc#include "ga_mpi.fh"
#ifndef GA_MPI_H_
#define GA_MPI_H_
!     usage: call ga_mpi_comm(OUT integer comm)
!     usage: call ga_mpi_comm_pgroup(OUT integer comm, IN integer grp)
!     usage: call ga_mpi_comm_pgroup_default(OUT integer comm)
      external ga_mpi_comm
      external ga_mpi_comm_pgroup
      external ga_mpi_comm_pgroup_default

#endif /* GA_MPI_H_ */

#endif /* GA_DEV */


#ifdef MPI4
#include "stupid_mpi4.fh"
#endif

*     **** local variables ****
      integer MASTER
      parameter (MASTER=0)
      integer i,mpierr,mygroup,myio,mydepth
      character*50 nwfilename
      character*7 c_index_name
      external    c_index_name
      integer  util_sgroup_set_ioname,util_sgroup_mygroup
      integer  util_sgroup_depth
      external util_sgroup_set_ioname,util_sgroup_mygroup
      external util_sgroup_depth


c*     **** MPI initiializer *****
cc     call MPI_INIT(mpierr)
c      call MPI_COMM_RANK(MPI_COMM_WORLD,taskid,mpierr)
c      call MPI_COMM_SIZE(MPI_COMM_WORLD,np,mpierr)

#ifdef USE_SUBGROUPS

#ifdef GA_DEV
      call ga_mpi_comm_pgroup_default(comm_world)
#else
      comm_world = MPI_COMM_WORLD
#endif

#else
      comm_world = MPI_COMM_WORLD
#endif

#ifdef MPI4
      stupid_world   = comm_world
      stupid_sum     = MPI_SUM
      stupid_max     = MPI_MAX
      stupid_min     = MPI_MIN
      stupid_character = MPI_CHARACTER
      stupid_integer = MPI_INTEGER
      stupid_double  = MPI_DOUBLE_PRECISION
      stupid_complex = MPI_DOUBLE_COMPLEX
      stupid_comm_i   = comm_world
      stupid_comm_j   = -99
      stupid_comm_k   = -99
#endif

#ifdef USE_SUBGROUPS
#ifdef GA_DEV
      call MPI_COMM_RANK(comm_world,taskid,mpierr)
      call MPI_COMM_SIZE(comm_world,np,mpierr)
      mydepth = util_sgroup_depth()
      mygroup = util_sgroup_mygroup()
      if ((taskid.eq.MASTER).and.(mydepth.gt.0)) then
         myio = util_sgroup_set_ioname(
     >              "nwpw"//c_index_name(mygroup)//".output")
      end if
#else
      np     = nnodes()
      taskid = nodeid()
#endif

#else
      np     = nnodes()
      taskid = nodeid()
#endif


*     **** set up 3d processor grid = np x 1 x 1****
      if (.not.BA_alloc_get(mt_int,np,'procNd',procNd(2),procNd(1)))
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
      comm_i   = comm_world
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
#include "bafdecls.fh"
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

      if (.not.BA_push_get(mt_int,np,'tmppp2',tmp(2),tmp(1)))
     >  call errquit('Parallel2d_init:out of stack memory',0, MA_ERR)

#ifdef MPI4
      if (.not.BA_push_get(mt_int,np,'tmppp3',tmp2(2),tmp2(1)))
     >  call errquit('Parallel2d_init:out of stack memory',0, MA_ERR)
      stupid_np_i = np_i
      stupid_np_j = np_j
*     **** set global processor group ****
!$OMP BARRIER
!$OMP MASTER
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

!$OMP END MASTER
!$OMP BARRIER

      group_i = stupid_group_i
      comm_i  = stupid_comm_i
      group_j = stupid_group_j
      comm_j  = stupid_comm_j
      if (.not.BA_pop_stack(tmp2(2)))
     >  call errquit('Parallel2d_init:popping stack memory',0, MA_ERR)

#else
*     **** set global processor group ****
!$OMP BARRIER
!$OMP MASTER
      call MPI_COMM_group(comm_world,mpi_group,ierr)

      do i=0,np_i-1
        int_mb(tmp(1)+i) = int_mb(procNd(1) + i + taskid_j*np_i) 
      end do
      call MPI_Group_incl(mpi_group,np_i,int_mb(tmp(1)),group_i,ierr)
      call MPI_Comm_create(comm_world,group_i,comm_i,  ierr)

      do j=0,np_j-1
        int_mb(tmp(1)+j) = int_mb(procNd(1) + taskid_i + j*np_i) 
      end do
      call MPI_Group_incl(mpi_group,np_j,int_mb(tmp(1)),group_j,ierr)
      call MPI_Comm_create(comm_world,group_j,comm_j,  ierr)
!$OMP END MASTER
!$OMP BARRIER
#endif

      if (.not.BA_pop_stack(tmp(2)))
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
#include "bafdecls.fh"
#include "errquit.fh"

#ifdef MPI4
#include "stupid_mpi4.fh"
#else
#include "mpif.h"
#endif

#ifdef MPI4
      if (np_j.gt.1) then
*      **** free comm_i and comm_j communicators ****
!$OMP BARRIER
!$OMP MASTER
      call MPI_Comm_free(stupid_comm_i,  stupid_ierr)
      call MPI_Group_free(stupid_group_i,stupid_ierr)
      call MPI_Comm_free(stupid_comm_j,  stupid_ierr)
      call MPI_Group_free(stupid_group_j,stupid_ierr)
!$OMP END MASTER
!$OMP BARRIER
      end if
#else
*     *** local variable ***
      integer mpierr

      if (np_j.gt.1) then
*      **** free comm_i and comm_j communicators ****
!$OMP BARRIER
!$OMP MASTER
      call MPI_Comm_free(comm_i,  mpierr)
      call MPI_Group_free(group_i,mpierr)
      call MPI_Comm_free(comm_j,  mpierr)
      call MPI_Group_free(group_j,mpierr)
!$OMP END MASTER
!$OMP BARRIER
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
#include "bafdecls.fh"
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

      if (.not.BA_push_get(mt_int,np,'tmppp2',tmp(2),tmp(1)))
     >  call errquit('Parallel3d_init:out of stack memory',0, MA_ERR)

#ifdef MPI4
      if (.not.BA_push_get(mt_int,np,'tmppp3',tmp2(2),tmp2(1)))
     >  call errquit('Parallel3d_init:out of stack memory',0, MA_ERR)
      stupid_np_i = np_i
      stupid_np_j = np_j
      stupid_np_k = np_k
*     **** set global processor group ****
!$OMP BARRIER
!$OMP MASTER
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

!$OMP END MASTER
!$OMP BARRIER
      group_i = stupid_group_i
      comm_i  = stupid_comm_i
      group_j = stupid_group_j
      comm_j  = stupid_comm_j
      group_k = stupid_group_k
      comm_k  = stupid_comm_k
      if (.not.BA_pop_stack(tmp2(2)))
     >  call errquit('Parallel3d_init:popping stack memory',0, MA_ERR)

#else
*     **** set global processor group ****
      call MPI_COMM_group(comm_world,mpi_group,ierr)

      do i=0,np_i-1
        int_mb(tmp(1)+i) = int_mb(procNd(1) 
     >                           + i 
     >                           + taskid_j*np_i 
     >                           + taskid_k*np_i*np_j) 
      end do
      call MPI_Group_incl(mpi_group,np_i,int_mb(tmp(1)),group_i,ierr)
      call MPI_Comm_create(comm_world,group_i,comm_i,ierr)

      do j=0,np_j-1
        int_mb(tmp(1)+j) = int_mb(procNd(1) 
     >                           + taskid_i 
     >                           + j*np_i
     >                           + taskid_k*np_i*np_j) 
      end do
      call MPI_Group_incl(mpi_group,np_j,int_mb(tmp(1)),group_j,ierr)
      call MPI_Comm_create(comm_world,group_j,comm_j,ierr)

      do k=0,np_k-1
        int_mb(tmp(1)+k) = int_mb(procNd(1) 
     >                           + taskid_i 
     >                           + taskid_j*np_i
     >                           + k*np_i*np_j) 
      end do
      call MPI_Group_incl(mpi_group,np_k,int_mb(tmp(1)),group_k,ierr)
      call MPI_Comm_create(comm_world,group_k,comm_k,ierr)
#endif

      if (.not.BA_pop_stack(tmp(2)))
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
#include "bafdecls.fh"
#include "errquit.fh"

#ifdef MPI4
#include "stupid_mpi4.fh"
#else
#include "mpif.h"
#endif

#ifdef MPI4
*      **** free comm_i and comm_j communicators ****
!$OMP BARRIER
!$OMP MASTER
      call MPI_Comm_free(stupid_comm_i,  stupid_ierr)
      call MPI_Group_free(stupid_group_i,stupid_ierr)
      call MPI_Comm_free(stupid_comm_j,  stupid_ierr)
      call MPI_Group_free(stupid_group_j,stupid_ierr)
      call MPI_Comm_free(stupid_comm_k,  stupid_ierr)
      call MPI_Group_free(stupid_group_k,stupid_ierr)
!$OMP END MASTER
!$OMP BARRIER
#else

*     *** local variable ***
      integer mpierr

*      **** free comm_i and comm_j communicators ****
!$OMP BARRIER
!$OMP MASTER
      call MPI_Comm_free(comm_i,  mpierr)
      call MPI_Group_free(group_i,mpierr)
      call MPI_Comm_free(comm_j,  mpierr)
      call MPI_Group_free(group_j,mpierr)
      call MPI_Comm_free(comm_k,  mpierr)
      call MPI_Group_free(group_k,mpierr)
!$OMP END MASTER
!$OMP BARRIER
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
     >                       MPI_MAX,comm_world,mpierr)
         sum = sumall
      end if
#endif
      return
      end


*     ***********************************
*     *                                 *
*     *         Parallel_IMaxAll        *
*     *                                 *
*     ***********************************

      subroutine Parallel_IMaxAll(isum)
c     implicit none
      integer  isum

#include "Parallel.fh"

#ifdef MPI4
#include "stupid_mpi4.fh"
#else
#include "mpif.h"
#endif
#ifdef MPI4
      integer*4 isumall,isum2
      if (np.gt.1) then
         isum2 = isum
         stupid_msglen = 1
         call MPI_Allreduce(isum2,isumall,stupid_msglen,stupid_integer,
     >                      stupid_max,stupid_world,stupid_ierr)
         isum = isumall
      end if
#else
      integer msglen,mpierr
      integer isumall
      if (np.gt.1) then
         msglen = 1
         call MPI_Allreduce(isum,isumall,msglen,MPI_INTEGER,
     >                       MPI_MAX,comm_world,mpierr)
         isum = isumall
      end if
#endif
      return
      end




*     ***********************************
*     *                                 *
*     *     Parallel_SumAll_private     *
*     *                                 *
*     ***********************************

* sum is assumed to be private an replicated

      subroutine Parallel_SumAll_private(sum)
c     implicit none
      real*8  sum

#include "Parallel.fh"
#ifdef MPI4
#include "stupid_mpi4.fh"
#else
#include "mpif.h"
#endif
c     *** make sumall shared ***
      real*8 sumall
      common /psumall_block/ sumall

#ifdef MPI4
      if (np.gt.1) then
!$OMP MASTER
         stupid_msglen = 1
         call MPI_Allreduce(sum,sumall,stupid_msglen,stupid_double,
     >                      stupid_sum,stupid_world,stupid_ierr)
!$OMP END MASTER
!$OMP BARRIER
         sum = sumall
      end if
#else
      integer msglen,mpierr

      if (np.gt.1) then
!$OMP MASTER
         msglen = 1
         call MPI_Allreduce(sum,sumall,msglen,MPI_DOUBLE_PRECISION,
     >                       MPI_SUM,comm_world,mpierr)
!$OMP END MASTER
!$OMP BARRIER
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
!$OMP MASTER
         stupid_msglen = 1
         call MPI_Allreduce(sum,sumall,stupid_msglen,stupid_double,
     >                      stupid_sum,stupid_world,stupid_ierr)
         sum = sumall
!$OMP END MASTER
!$OMP BARRIER
      end if
#else
      integer msglen,mpierr

      if (np.gt.1) then
!$OMP MASTER
         msglen = 1
         call MPI_Allreduce(sum,sumall,msglen,MPI_DOUBLE_PRECISION,
     >                       MPI_SUM,comm_world,mpierr)
         sum = sumall
!$OMP END MASTER
!$OMP BARRIER
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
     >                       MPI_SUM,comm_world,mpierr)
         sum = sumall
      end if
#endif

      return
      end


*     ***********************************
*     *                                 *
*     *  Parallel_Vector_SumAll_master  *
*     *                                 *
*     ***********************************

      subroutine Parallel_Vector_SumAll_master(n,sum)
c     implicit none
      integer n
      real*8  sum(*)

#include "mafdecls.fh"
#include "Parallel.fh"
#include "mpif.h"

*     **** temporary workspace ****
      integer sumall(2)

*     **** local variable ****
      logical value
      integer msglen,mpierr


      call nwpw_timing_start(2)
      if (np.gt.1) then

*     ***** allocate temporary space ****
      value = MA_push_get(mt_dbl,n,'sumall',sumall(2),sumall(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)

      msglen = n
      call MPI_Allreduce(sum,sumall(1),msglen,
     >                MPI_DOUBLE_PRECISION,
     >                MPI_SUM,comm_world,mpierr)
      call dcopy(n,sumall(1),1,sum,1)


      value = MA_pop_stack(sumall(2))
      if (.not.value) call errquit('dealloc stack memory',0,MA_ERR)

      end if
      call nwpw_timing_end(2)
      return
      end

*     ***********************************
*     *                                 *
*     *  Parallel_Vector_SumAll_master_tmp  *
*     *                                 *
*     ***********************************

      subroutine Parallel_Vector_SumAll_master_tmp(n,sum0,sumall)
c     implicit none
      integer n
      real*8  sum0(*),sumall(*)

#include "mafdecls.fh"
#include "Parallel.fh"
#include "mpif.h"

*     **** local variable ****
      logical value
      integer msglen,mpierr


      call nwpw_timing_start(2)
      if (np.gt.1) then
      msglen = n
      call MPI_Allreduce(sum0,sumall,msglen,
     >                MPI_DOUBLE_PRECISION,
     >                MPI_SUM,comm_world,mpierr)
      call dcopy(n,sumall,1,sum0,1)
      end if
      call nwpw_timing_end(2)
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

#include "bafdecls.fh"
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
      value = BA_push_get(mt_dbl,n,'sumall',sumall(2),sumall(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)

!$OMP MASTER
#ifdef MPI4
      stupid_msglen = n
      call MPI_Allreduce(sum,dbl_mb(sumall(1)),stupid_msglen,
     >                stupid_double,
     >                stupid_sum,stupid_world,stupid_ierr)
#else
      msglen = n
      call MPI_Allreduce(sum,dbl_mb(sumall(1)),msglen,
     >                MPI_DOUBLE_PRECISION,
     >                MPI_SUM,comm_world,mpierr)
#endif
      call dcopy(n,dbl_mb(sumall(1)),1,sum,1)
!$OMP END MASTER
!$OMP BARRIER

      value = BA_pop_stack(sumall(2))
      if (.not.value) call errquit('dealloc stack memory',0,MA_ERR)

      end if
      call nwpw_timing_end(2)
      return
      end

*     ***********************************
*     *                                 *
*     *      Parallel_Vector_SumAll2    *
*     *                                 *
*     ***********************************

      subroutine Parallel_Vector_SumAll2(n,sum0,sumall)
c     implicit none
      integer n
      real*8  sum0(*),sumall(*)

#include "Parallel.fh"
#include "mpif.h"

*     **** local variable ****
      logical value
      integer msglen,mpierr

      call nwpw_timing_start(2)

!$OMP BARRIER
!$OMP MASTER
      if (np.gt.1) then
      msglen = n
      call MPI_Allreduce(sum0,sumall,msglen,
     >                MPI_DOUBLE_PRECISION,
     >                MPI_SUM,comm_world,mpierr)
      else
         call dcopy(n,sum0,1,sumall,1)
      end if
!$OMP END MASTER
!$OMP BARRIER

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

#include "bafdecls.fh"
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
      value = BA_push_get(mt_int,n,'sumall',sumall(2),sumall(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)

#ifdef MPI4
      if (.not.BA_push_get(mt_int,n,'sumall1',sumall1(2),sumall1(1)))
     > call errquit('out of stack memory',0, MA_ERR)

      stupid_msglen = n
      call stupid_icopy8to4(n,sum,int_mb(sumall1(1)))
      call MPI_Allreduce(int_mb(sumall1(1)),
     >                   int_mb(sumall(1)),stupid_msglen,
     >                stupid_integer,
     >                stupid_sum,stupid_world,stupid_ierr)
      call stupid_icopy4to8(n,int_mb(sumall(1)),sum)
      if (.not.BA_pop_stack(sumall1(2)))
     >   call errquit('error popping stack',0,MA_ERR)
#else
      msglen = n
      call MPI_Allreduce(sum,int_mb(sumall(1)),msglen,
     >                MPI_INTEGER,
     >                MPI_SUM,comm_world,mpierr)
      call icopy(n,int_mb(sumall(1)),1,sum,1)
#endif
      value = BA_pop_stack(sumall(2))
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
     >                  psend,comm_world,ierr)
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
!$OMP MASTER
         stupid_msglen = nsize
         tpsend        = psend
         call MPI_Bcast(sum,stupid_msglen,stupid_double,
     >                  tpsend,stupid_world,stupid_ierr)
!$OMP END MASTER
!$OMP BARRIER
      end if
#else
      integer ierr

      if (np.gt.1) then
!$OMP MASTER
         call MPI_Bcast(sum,nsize,MPI_DOUBLE_PRECISION,
     >                  psend,comm_world,ierr)
!$OMP END MASTER
!$OMP BARRIER
      end if
#endif

      return
      end

*     ***********************************
*     *                                 *
*     *      Parallel_Brdcst_ivalue     *
*     *                                 *
*     ***********************************

      subroutine Parallel_Brdcst_ivalue(psend,isum)
      implicit none
      integer psend
      integer isum

#include "Parallel.fh"
#ifdef MPI4
#include "stupid_mpi4.fh"
#else
#include "mpif.h"
#endif


#ifdef MPI4
      integer*4 tpsend,isum4

      if (np.gt.1) then
         stupid_msglen = 1
         tpsend        = psend
         isum4 = isum
         call MPI_Bcast(isum4,stupid_msglen,stupid_integer,
     >                  tpsend,stupid_world,stupid_ierr)
         isum = isum4
      end if
#else
      integer ierr

      if (np.gt.1) then
         call MPI_Bcast(isum,1,MPI_DOUBLE_PRECISION,
     >                  psend,comm_world,ierr)
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

#include "bafdecls.fh"
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
         value = BA_push_get(mt_int,nsize,'sumall',sumall(2),sumall(1))
         if (.not. value) call errquit('out of stack memory',0, MA_ERR)

         call stupid_icopy8to4(nsize,isum,int_mb(sumall(1)))
         call MPI_Bcast(int_mb(sumall(1)),stupid_msglen,stupid_integer,
     >                  tpsend,stupid_world,stupid_ierr)
         call stupid_icopy4to8(nsize,int_mb(sumall(1)),isum)

         value = BA_pop_stack(sumall(2))
         if (.not. value) call errquit('error popping stack',0, MA_ERR)
      end if
#else
      integer ierr

      if (np.gt.1) then
         call MPI_Bcast(isum,nsize,MPI_INTEGER,
     >                  psend,comm_world,ierr)
      end if
#endif

      return
      end




#ifdef MPI4
      subroutine stupid_icopy8to4(n,array8,array4)
      implicit none
      integer n
#ifdef EXT_INT
      integer*8 array8(*)
#else
      integer*4 array8(*)
#endif
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
#ifdef EXT_INT
      integer*8 array8(*)
#else
      integer*4 array8(*)
#endif
      integer i
      do i=1,n
          array8(i) = array4(i)
      end do
      return
      end
#endif



*     ***********************************
*     *                                 *
*     *      Parallel_mpiWaitAll0       *
*     *                                 *
*     ***********************************

      subroutine Parallel_mpiWaitAll0(nreq,req)
      implicit none
      integer nreq,req(*)

#include "bafdecls.fh"
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
      if (.not.MA_push_get(mt_int,MPI_STATUS_SIZE*nreq*4,
     >                     'status',status(2),status(1)))
     > call errquit('Parallel_mpiWaitAll0:out of stack',0,MA_ERR)



*     **** wait for completion of mp_send, also do a sync ****

#ifdef MPI4

      if (.not.MA_push_get(mt_int,nreq,'rqst2',request2(2),request2(1)))
     > call errquit(' Parallel_mpiWaitAll0:out of stack',1,MA_ERR)
      call stupid_icopy8to4(nreq,req,int_mb(request2(1)))
      stupid_msglen = nreq
      call MPI_WAITALL(stupid_msglen,
     >                 int_mb(request2(1)),
     >                 int_mb(status(1)),stupid_ierr)
c       call MPI_WAITALL(stupid_msglen,
c     >                 int_mb(request2(1)),
c     >                 MPI_STATUSES_IGNORE,stupid_ierr)
      if (.not.MA_pop_stack(request2(2)))
     > call errquit(' Parallel_mpiWaitAll0:popping stack',1,MA_ERR)

c     *** danger! ***
      if (stupid_ierr.ne.0) then
         write(*,*) "MPI_WAITALL Failed on taskid=",taskid,
     >              "  mpierr=",stupid_ierr
         call util_flush(6)
      end if

#else
      call MPI_WAITALL(nreq,req,int_mb(status(1)),mpierr)
c      call MPI_WAITALL(nreq,req,MPI_STATUSES_IGNORE,mpierr)
#endif

*     *** may need to check status here??? ***

      !*** deallocate status memory ***
      if (.not.MA_pop_stack(status(2)))
     > call errquit(' Parallel_mpiWaitAll0:popping stack',0,MA_ERR)

      end if
      return
      end



*     ***********************************
*     *                                 *
*     *      Parallel_mpiWaitAll        *
*     *                                 *
*     ***********************************

      subroutine Parallel_mpiWaitAll(nreq,req)
      implicit none
      integer nreq,req(*)

#include "bafdecls.fh"
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

!$OMP MASTER
     
      if (nreq.gt.0) then

*     **** allocate status memory ****
      if (.not.MA_push_get(mt_int,MPI_STATUS_SIZE*nreq*4,
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
c       call MPI_WAITALL(stupid_msglen,
c     >                 int_mb(request2(1)),
c     >                 MPI_STATUSES_IGNORE,stupid_ierr)
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
c      call MPI_WAITALL(nreq,req,MPI_STATUSES_IGNORE,mpierr)
#endif

*     *** may need to check status here??? ***

      !*** deallocate status memory ***
      if (.not.MA_pop_stack(status(2)))
     > call errquit(' Parallel_mpiWaitAll:popping stack',0,MA_ERR)

      end if
!$OMP END MASTER
!$OMP BARRIER
      return
      end




#ifdef MPI4
      integer*4 function getcomm(ic)
      implicit none
      integer   ic

#include "stupid_mpi4.fh"

      integer*4 icc
      if (ic.eq.1) then
        icc = stupid_comm_i
      else if (ic.eq.2) then
        icc = stupid_comm_j
      else if (ic.eq.3) then
        icc = stupid_comm_k
      else 
        icc = stupid_world
      end if

      getcomm = icc
      return
      end
#else
      integer function getcomm(ic)
      implicit none
      integer ic

#include "Parallel.fh"
#include "mpif.h"

      integer icc
      if (ic.eq.1) then
        icc = comm_i
      else if (ic.eq.2) then
        icc = comm_j
      else if (ic.eq.3) then
        icc = comm_k
      else 
        icc = comm_world
      end if

      getcomm = icc
      return
      end
#endif

*     ***********************************
*     *                                 *
*     *         Parallela_MaxAll        *
*     *                                 *
*     ***********************************

      subroutine Parallela_MaxAll(ic,sum)
c     implicit none
      real*8  sum

#include "Parallel.fh"

#ifdef MPI4
#include "stupid_mpi4.fh"
#else
#include "mpif.h"
#endif
      integer ic
      real*8 sumall
#ifdef MPI4
      integer*4 getcomm
      external  getcomm

      if (np.gt.1) then
         stupid_msglen = 1
         call MPI_Allreduce(sum,sumall,stupid_msglen,stupid_double,
     >                      stupid_max,getcomm(ic),stupid_ierr)
         sum = sumall
      end if
#else
      integer msglen,mpierr
      integer  getcomm
      external getcomm

      if (np.gt.1) then
         msglen = 1
         call MPI_Allreduce(sum,sumall,msglen,MPI_DOUBLE_PRECISION,
     >                       MPI_MAX,getcomm(ic),mpierr)
         sum = sumall
      end if
#endif
      return
      end





*     ***********************************
*     *                                 *
*     *         Parallela_SumAll        *
*     *                                 *
*     ***********************************

      subroutine Parallela_SumAll(ic,sum)
c     implicit none
      integer ic
      real*8  sum

#include "Parallel.fh"

#ifdef MPI4
#include "stupid_mpi4.fh"
#else
#include "mpif.h"
#endif

      real*8 sumall

#ifdef MPI4
      integer*4 getcomm
      external  getcomm

      if (np.gt.1) then
         stupid_msglen = 1
         call MPI_Allreduce(sum,sumall,stupid_msglen,stupid_double,
     >                      stupid_sum,getcomm(ic),stupid_ierr)
         sum = sumall
      end if
#else
      integer msglen,mpierr
      integer  getcomm
      external getcomm

      if (np.gt.1) then
         msglen = 1
         call MPI_Allreduce(sum,sumall,msglen,MPI_DOUBLE_PRECISION,
     >                       MPI_SUM,getcomm(ic),mpierr)
         sum = sumall
      end if
#endif

      return
      end


*     ***********************************
*     *                                 *
*     *         Parallela_ISumAll       *
*     *                                 *
*     ***********************************

      subroutine Parallela_ISumAll(ic,sum)
c     implicit none
      integer ic,sum


#include "Parallel.fh"

#ifdef MPI4
#include "stupid_mpi4.fh"
#else
#include "mpif.h"
#endif

#ifdef MPI4
      integer*4 tsum,sumall
      integer*4 getcomm
      external  getcomm

      if (np.gt.1) then
         stupid_msglen = 1
         tsum = sum
         call MPI_Allreduce(tsum,sumall,stupid_msglen,stupid_integer,
     >                       stupid_sum,getcomm(ic),stupid_ierr)
         sum = sumall
      end if
#else
      integer msglen,mpierr
      integer sumall
      integer  getcomm
      external getcomm

      if (np.gt.1) then
         msglen = 1

         call MPI_Allreduce(sum,sumall,msglen,MPI_INTEGER,
     >                       MPI_SUM,getcomm(ic),mpierr)
         sum = sumall
      end if
#endif

      return
      end

*     ***********************************
*     *                                 *
*     *      Parallela_Vector_SumAll    *
*     *                                 *
*     ***********************************

      subroutine Parallela_Vector_SumAll(ic,n,sum)
c     implicit none
      integer ic,n
      real*8  sum(*)

#include "bafdecls.fh"
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

#ifdef MPI4
      integer*4 getcomm
      external  getcomm
#else
      integer  getcomm
      external getcomm
#endif


      call nwpw_timing_start(2)
      if (np.gt.1) then

*     ***** allocate temporary space ****
      value = BA_push_get(mt_dbl,n,'sumall',sumall(2),sumall(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)

#ifdef MPI4
      stupid_msglen = n
      call MPI_Allreduce(sum,dbl_mb(sumall(1)),stupid_msglen,
     >                stupid_double,
     >                stupid_sum,getcomm(ic),stupid_ierr)
#else
      msglen = n
      call MPI_Allreduce(sum,dbl_mb(sumall(1)),msglen,
     >                MPI_DOUBLE_PRECISION,
     >                MPI_SUM,getcomm(ic),mpierr)
#endif

      call dcopy(n,dbl_mb(sumall(1)),1,sum,1)
      value = BA_pop_stack(sumall(2))

      end if
      call nwpw_timing_end(2)
      return
      end


*     ***********************************
*     *                                 *
*     *      Parallela_Vector_ISumAll   *
*     *                                 *
*     ***********************************

      subroutine Parallela_Vector_ISumAll(ic,n,sum)
c     implicit none
      integer ic,n
      integer  sum(*)

#include "bafdecls.fh"
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

#ifdef MPI4
      integer*4 getcomm
      external  getcomm
#else
      integer  getcomm
      external getcomm
#endif

      call nwpw_timing_start(2)

      if (np.gt.1) then

*     ***** allocate temporary space ****
      value = BA_push_get(mt_int,n,'sumall',sumall(2),sumall(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)

#ifdef MPI4
      if (.not.BA_push_get(mt_int,n,'sumall1',sumall1(2),sumall1(1)))
     > call errquit('out of stack memory',0, MA_ERR)

      stupid_msglen = n
      call stupid_icopy8to4(n,sum,int_mb(sumall1(1)))
      call MPI_Allreduce(int_mb(sumall1(1)),
     >                   int_mb(sumall(1)),stupid_msglen,
     >                stupid_integer,
     >                stupid_sum,getcomm(ic),stupid_ierr)
      call stupid_icopy4to8(n,int_mb(sumall(1)),sum)
      if (.not.BA_pop_stack(sumall1(2)))
     >   call errquit('error popping stack',0,MA_ERR)
#else
      msglen = n
      call MPI_Allreduce(sum,int_mb(sumall(1)),msglen,
     >                MPI_INTEGER,
     >                MPI_SUM,getcomm(ic),mpierr)
      call icopy(n,int_mb(sumall(1)),1,sum,1)
#endif
      value = BA_pop_stack(sumall(2))
      if (.not. value) call errquit('error popping stack',0, MA_ERR)

      end if

      call nwpw_timing_end(2)
      return
      end





*     ***********************************
*     *                                 *
*     *      Parallela_Brdcst_value     *
*     *                                 *
*     ***********************************

      subroutine Parallela_Brdcst_value(ic,psend,sum)
      implicit none
      integer ic,psend
      real*8  sum

#include "Parallel.fh"
#ifdef MPI4
#include "stupid_mpi4.fh"
#else
#include "mpif.h"
#endif


#ifdef MPI4
      integer*4 tpsend
      integer*4 getcomm
      external  getcomm

      if (np.gt.1) then
         stupid_msglen = 1
         tpsend        = psend
         call MPI_Bcast(sum,stupid_msglen,stupid_double,
     >                  tpsend,getcomm(ic),stupid_ierr)
      end if
#else
      integer ierr
      integer  getcomm
      external getcomm

      if (np.gt.1) then
         call MPI_Bcast(sum,1,MPI_DOUBLE_PRECISION,
     >                  psend,getcomm(ic),ierr)
      end if
#endif

      return
      end




*     ***********************************
*     *                                 *
*     *      Parallela_Brdcst_values    *
*     *                                 *
*     ***********************************

      subroutine Parallela_Brdcst_values(ic,psend,nsize,sum)
      implicit none
      integer ic,psend,nsize
      real*8  sum(*)

#include "Parallel.fh"
#ifdef MPI4
#include "stupid_mpi4.fh"
#else
#include "mpif.h"
#endif


#ifdef MPI4
      integer*4 tpsend
      integer*4 getcomm
      external  getcomm

      if (np.gt.1) then
         stupid_msglen = nsize
         tpsend        = psend
         call MPI_Bcast(sum,stupid_msglen,stupid_double,
     >                  tpsend,getcomm(ic),stupid_ierr)
      end if
#else
      integer ierr
      integer  getcomm
      external getcomm

      if (np.gt.1) then
         call MPI_Bcast(sum,nsize,MPI_DOUBLE_PRECISION,
     >                  psend,getcomm(ic),ierr)
      end if
#endif

      return
      end



*     ***********************************
*     *                                 *
*     *      Parallela_Brdcst_ivalue    *
*     *                                 *
*     ***********************************

      subroutine Parallela_Brdcst_ivalue(ic,psend,isum)
      implicit none
      integer ic,psend
      integer isum

#include "Parallel.fh"
#ifdef MPI4
#include "stupid_mpi4.fh"
#else
#include "mpif.h"
#endif


#ifdef MPI4
      integer*4 tpsend,isum4
      integer*4 getcomm
      external  getcomm

      if (np.gt.1) then
         stupid_msglen = 1
         tpsend        = psend
         isum4         = isum
         call MPI_Bcast(isum4,stupid_msglen,stupid_double,
     >                  tpsend,getcomm(ic),stupid_ierr)
         isum         = isum4
      end if
#else
      integer ierr
      integer  getcomm
      external getcomm

      if (np.gt.1) then
         call MPI_Bcast(isum,1,MPI_DOUBLE_PRECISION,
     >                  psend,getcomm(ic),ierr)
      end if
#endif

      return
      end





*     ***********************************
*     *                                 *
*     *      Parallela_Brdcst_ivalues   *
*     *                                 *
*     ***********************************

      subroutine Parallela_Brdcst_ivalues(ic,psend,nsize,isum)
      implicit none
      integer ic,psend,nsize
      integer isum(*)

#include "bafdecls.fh"
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
      integer*4 getcomm
      external  getcomm

      if (np.gt.1) then
         stupid_msglen = nsize
         tpsend        = psend
         value = BA_push_get(mt_int,nsize,'sumall',sumall(2),sumall(1))
         if (.not. value) call errquit('out of stack memory',0, MA_ERR)

         call stupid_icopy8to4(nsize,isum,int_mb(sumall(1)))
         call MPI_Bcast(int_mb(sumall(1)),stupid_msglen,stupid_integer,
     >                  tpsend,getcomm(ic),stupid_ierr)
         call stupid_icopy4to8(nsize,int_mb(sumall(1)),isum)

         value = BA_pop_stack(sumall(2))
         if (.not. value) call errquit('error popping stack',0, MA_ERR)
      end if
#else
      integer ierr
      integer  getcomm
      external getcomm

      if (np.gt.1) then
         call MPI_Bcast(isum,nsize,MPI_INTEGER,
     >                  psend,getcomm(ic),ierr)
      end if
#endif

      return
      end





*     ***********************************
*     *                                 *
*     *      Parallela_start_rotate     *
*     *                                 *
*     ***********************************

      subroutine Parallela_start_rotate(ic,shift,msgtype0,
     >                                  A1,nsize1,
     >                                  A2,nsize2,request)
      implicit none
      integer ic,shift,msgtype0
      real*8  A1(*)
      integer nsize1
      real*8  A2(*)
      integer nsize2
      integer request(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "Parallel.fh"

#ifdef MPI4
#include "stupid_mpi4.fh"
#else
#include "mpif.h"
#endif

*     **** local variables ****
      integer mynp,mytaskid,msglen,msgtype,mpierr,proc_to,proc_from

*     ***** external functions ****
#ifdef MPI4
      integer*4 getcomm
      external  getcomm
#else
      integer  getcomm
      external getcomm
#endif

      if (ic.eq.1) then
         mynp     = np_i
         mytaskid = taskid_i
      else if (ic.eq.2) then
         mynp     = np_j
         mytaskid = taskid_j
      else if (ic.eq.3) then
         mynp     = np_k
         mytaskid = taskid_k
      else
         mynp     = np
         mytaskid = taskid
      end if
      proc_to   = mod(mytaskid+shift+mynp,mynp)
      proc_from = mod(mytaskid-shift+mynp,mynp)

#ifdef MPI4
         if (nsize2.gt.0) then
            stupid_msglen = nsize2
            stupid_type   = msgtype0
            stupid_taskid = proc_from
            call MPI_IRECV(A2,
     >                  stupid_msglen,stupid_double,
     >                  stupid_taskid,
     >                  stupid_type,getcomm(ic),
     >                  stupid_request,stupid_ierr)
            request(1) = 0
            request(3) = 1
            request(1) = stupid_request
         else
            request(3) = 0
         end if

         if (nsize1.gt.0) then
            stupid_msglen = nsize1
            stupid_type   = msgtype0
            stupid_taskid = proc_to
            call MPI_ISEND(A1,
     >                  stupid_msglen,stupid_double,
     >                  stupid_taskid,
     >                  stupid_type,getcomm(ic),
     >                  stupid_request,stupid_ierr)
            request(2) = 0
            request(4) = 1
            request(2) = stupid_request
         else
            request(4) = 0
         end if

#else
         if (nsize2.gt.0) then
         msglen  = nsize2
         msgtype = msgtype0
         call MPI_IRECV(A2,msglen,MPI_DOUBLE_PRECISION,
     >                  proc_from,
     >                  msgtype,getcomm(ic),
     >                  request(1),mpierr)
            request(3) = 1
         else
            request(3) = 0
         end if

         if (nsize1.gt.0) then
         msglen  = nsize1
         msgtype = msgtype0
         call MPI_ISEND(A1,msglen,MPI_DOUBLE_PRECISION,
     >                  proc_to,
     >                  msgtype,getcomm(ic),
     >                  request(2),mpierr)
            request(4) = 1
         else
            request(4) = 0
         end if
#endif

      if ((request(3).eq.1).and.(request(4).eq.1)) then
         request(3) = 1
      else if (request(3).eq.1) then
         request(3) = 2
      else if (request(4).eq.1) then
         request(3) = 3
      else
         request(3) = 4
      end if

      return
      end



*     ***********************************
*     *                                 *
*     *      Parallela_start_Irotate     *
*     *                                 *
*     ***********************************

      subroutine Parallela_start_Irotate(ic,shift,msgtype0,
     >                                  IA1,nsize1,
     >                                  IA2,nsize2,request)
      implicit none
      integer ic,shift,msgtype0
      integer IA1(*)
      integer nsize1
      integer IA2(*)
      integer nsize2
      integer request(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "Parallel.fh"

#ifdef MPI4
#include "stupid_mpi4.fh"
#else
#include "mpif.h"
#endif

*     **** local variables ****
      integer mynp,mytaskid,msglen,msgtype,mpierr,proc_to,proc_from

*     ***** external functions ****
#ifdef MPI4
      integer*4 getcomm
      external  getcomm
#else
      integer  getcomm
      external getcomm
#endif

      if (ic.eq.1) then
         mynp     = np_i
         mytaskid = taskid_i
      else if (ic.eq.2) then
         mynp     = np_j
         mytaskid = taskid_j
      else if (ic.eq.3) then
         mynp     = np_k
         mytaskid = taskid_k
      else
         mynp     = np
         mytaskid = taskid
      end if
      proc_to   = mod(mytaskid+shift+mynp,mynp)
      proc_from = mod(mytaskid-shift+mynp,mynp)

#ifdef MPI4
         if (nsize2.gt.0) then
            stupid_msglen = nsize2
            stupid_type   = msgtype0
            stupid_taskid = proc_from
      write(*,*) "irecv, mytaskid,stupid_taskid,msglen=",mytaskid,
     >           stupid_taskid,stupid_msglen,IA2(1)
            call MPI_IRECV(IA2,
     >                  stupid_msglen,stupid_integer,
     >                  stupid_taskid,
     >                  stupid_type,getcomm(ic),
     >                  stupid_request,stupid_ierr)
            request(1) = 0
            request(1) = stupid_request
            request(3) = 1
         else
            request(3) = 0
         end if

         if (nsize1.gt.0) then
            stupid_msglen = nsize1
            stupid_type   = msgtype0
            stupid_taskid = proc_to
      write(*,*) "isend, mytaskid,stupid_taskid,msglen=",mytaskid,
     >           stupid_taskid,stupid_msglen,IA1(1)
            call MPI_ISEND(IA1,
     >                  stupid_msglen,stupid_integer,
     >                  stupid_taskid,
     >                  stupid_type,getcomm(ic),
     >                  stupid_request,stupid_ierr)
            request(2) = 0
            request(2) = stupid_request
            request(4) = 1
         else
            request(4) = 0
         end if

#else
         if (nsize2.gt.0) then
         msglen  = nsize2
         msgtype = msgtype0
         call MPI_IRECV(IA2,msglen,MPI_INTEGER,
     >                  proc_from,
     >                  msgtype,getcomm(ic),
     >                  request(1),mpierr)
            request(3) = 1
         else
            request(3) = 0
         end if

         if (nsize1.gt.0) then
         msglen  = nsize1
         msgtype = msgtype0
         call MPI_ISEND(IA1,msglen,MPI_INTEGER,
     >                  proc_to,
     >                  msgtype,getcomm(ic),
     >                  request(2),mpierr)
            request(4) = 1
         else
            request(4) = 0
         end if
#endif

      if ((request(3).eq.1).and.(request(4).eq.1)) then
         request(3) = 1
      else if (request(3).eq.1) then
         request(3) = 2
      else if (request(4).eq.1) then
         request(3) = 3
      else
         request(3) = 4
      end if

      return
      end



*     ***********************************
*     *                                 *
*     *      Parallela_end_rotate       *
*     *                                 *
*     ***********************************

      subroutine Parallela_end_rotate(request)
      implicit none
      integer request(*)

*     **** wait for completion of mp_send, also do a sync ****
      if (request(3).eq.1) then
         call Parallel_mpiWaitAll(2,request)
      else if (request(3).eq.2) then
         call Parallel_mpiWaitAll(1,request)
      else if (request(3).eq.3) then
         call Parallel_mpiWaitAll(1,request(2))
      endif

      return
      end





*     ***********************************
*     *                                 *
*     *      Parallel_send_values       *
*     *                                 *
*     ***********************************
      subroutine Parallel_send_values(pto,msgtype,nsize,rval)
      implicit none
      integer pto,msgtype,nsize
      real*8 rval(*)

#include "Parallel.fh"
#ifdef MPI4
#include "stupid_mpi4.fh"
#else
#include "mpif.h"
      integer mpierr
#endif

#ifdef MPI4
      stupid_msglen = nsize
      stupid_taskid = pto
      stupid_type = msgtype
      call MPI_SEND(rval,stupid_msglen,stupid_double,
     >              stupid_taskid,
     >              stupid_type,stupid_world,stupid_ierr)
#else
      call MPI_SEND(rval,nsize,MPI_DOUBLE_PRECISION,
     >              pto,msgtype,comm_world,mpierr)
#endif
      return
      end 

*     ***********************************
*     *                                 *
*     *      Parallel_send_ivalues      *
*     *                                 *
*     ***********************************
      subroutine Parallel_send_ivalues(pto,msgtype,nsize,ival)
      implicit none
      integer pto,msgtype,nsize
      integer ival(*)

#include "Parallel.fh"
#ifdef MPI4
#include "stupid_mpi4.fh"
#else
#include "mpif.h"
      integer mpierr
#endif

#ifdef MPI4
#include "bafdecls.fh"
#include "errquit.fh"

      logical value
      integer ival4(2)

*     ***** allocate temporary space ****
      value = BA_push_get(mt_int,nsize,'ival4',ival4(2),ival4(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)

      call stupid_icopy8to4(nsize,ival,int_mb(ival4(1)))
      stupid_msglen = nsize
      stupid_taskid = pto
      stupid_type = msgtype
      call MPI_SEND(int_mb(ival4(1)),stupid_msglen,stupid_integer,
     >              stupid_taskid,
     >              stupid_type,stupid_world,stupid_ierr)
      value = BA_pop_stack(ival4(2))
      if (.not.value) call errquit('popping stack memory',0,MA_ERR)
#else
      call MPI_SEND(ival,nsize,MPI_INTEGER,
     >              pto,msgtype,comm_world,mpierr)
#endif
      return
      end




*     ***********************************
*     *                                 *
*     *      Parallel_send_characters   *
*     *                                 *
*     ***********************************
      subroutine Parallel_send_characters(pto,msgtype,nsize,cval)
      implicit none
      integer pto,msgtype,nsize
      character cval(*)

#include "Parallel.fh"
#ifdef MPI4
#include "stupid_mpi4.fh"
#else
#include "mpif.h"
      integer mpierr
#endif

#ifdef MPI4
      stupid_msglen = nsize
      stupid_taskid = pto
      stupid_type   = msgtype
      call MPI_SEND(cval,stupid_msglen,stupid_character,
     >              stupid_taskid,
     >              stupid_type,stupid_world,stupid_ierr)
#else
      call MPI_SEND(cval,nsize,MPI_CHARACTER,
     >              pto,msgtype,comm_world,mpierr)
#endif
      return
      end



*     ***********************************
*     *                                 *
*     *      Parallel_recv_values       *
*     *                                 *
*     ***********************************
      subroutine Parallel_recv_values(pfrom,msgtype,nsize,rval)
      implicit none
      integer pfrom,msgtype,nsize
      real*8 rval(*)

#include "Parallel.fh"
#include "mpif.h"
#ifdef MPI4
#include "stupid_mpi4.fh"
#endif

#include "bafdecls.fh"
#include "errquit.fh"

      logical value
      integer status(2),mpierr

*     **** allocate status memory ****
      if (.not.BA_push_get(mt_int,MPI_STATUS_SIZE*4,
     >                     'status',status(2),status(1)))
     > call errquit('out of stack',0,MA_ERR)

#ifdef MPI4
      stupid_msglen = nsize
      stupid_taskid = pfrom
      stupid_type = msgtype
      call MPI_RECV(rval,stupid_msglen,stupid_double,
     >              stupid_taskid,
     >              stupid_type,stupid_world,int_mb(status(1)),
     >              stupid_ierr)
#else
      call MPI_RECV(rval,nsize,MPI_DOUBLE_PRECISION,
     >              pfrom,msgtype,comm_world,int_mb(status(1)),
     >              mpierr)
#endif

      value = BA_pop_stack(status(2))
      if (.not.value) call errquit('popping stack memory',0,MA_ERR)

      return
      end

*     ***********************************
*     *                                 *
*     *      Parallel_recv_ivalues       *
*     *                                 *
*     ***********************************
      subroutine Parallel_recv_ivalues(pfrom,msgtype,nsize,ival)
      implicit none
      integer pfrom,msgtype,nsize
      integer ival(*)

#include "Parallel.fh"
#include "mpif.h"
#ifdef MPI4
#include "stupid_mpi4.fh"
#endif

#include "bafdecls.fh"
#include "errquit.fh"

      logical value
      integer ival4(2),status(2),mpierr

*     **** allocate status memory ****
      if (.not.BA_push_get(mt_int,MPI_STATUS_SIZE*4,
     >                     'status',status(2),status(1)))
     > call errquit('out of stack',0,MA_ERR)

#ifdef MPI4
*     ***** allocate temporary space ****
      value = BA_push_get(mt_int,nsize,'ival4',ival4(2),ival4(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)

      stupid_msglen = nsize
      stupid_taskid = pfrom
      stupid_type = msgtype
      call MPI_RECV(int_mb(ival4(1)),stupid_msglen,stupid_integer,
     >              stupid_taskid,
     >              stupid_type,stupid_world,int_mb(status(1)),
     >              stupid_ierr)
      call stupid_icopy4to8(nsize,int_mb(ival4(1)),ival)
      value = BA_pop_stack(ival4(2))
      if (.not. value) call errquit('popping stack memory',0,MA_ERR)
      
#else
      call MPI_RECV(ival,nsize,MPI_INTEGER,
     >              pfrom,msgtype,comm_world,int_mb(status(1)),mpierr)
#endif

      value = BA_pop_stack(status(2))
      if (.not.value) call errquit('popping stack memory',0,MA_ERR)
      return
      end




*     ***********************************
*     *                                 *
*     *      Parallel_recv_characters   *
*     *                                 *
*     ***********************************
      subroutine Parallel_recv_characters(pfrom,msgtype,nsize,cval)
      implicit none
      integer  pfrom,msgtype,nsize
      character cval(*)

#include "Parallel.fh"

#include "mpif.h"

#ifdef MPI4
#include "stupid_mpi4.fh"
#endif

#include "bafdecls.fh"
#include "errquit.fh"

      logical value
      integer status(2),mpierr

*     **** allocate status memory ****
      if (.not.BA_push_get(mt_int,MPI_STATUS_SIZE*4,
     >                     'status',status(2),status(1)))
     > call errquit('out of stack',0,MA_ERR)

#ifdef MPI4
      stupid_msglen = nsize
      stupid_taskid = pfrom
      stupid_type   = msgtype
      call MPI_RECV(cval,stupid_msglen,stupid_character,
     >              stupid_taskid,
     >              stupid_type,stupid_world,int_mb(status(1)),
     >              stupid_ierr)
#else
      call MPI_RECV(cval,nsize,MPI_CHARACTER,
     >              pfrom,msgtype,comm_world,int_mb(status(1)),
     >              mpierr)
#endif

      value = BA_pop_stack(status(2))
      if (.not.value) call errquit('popping stack memory',0,MA_ERR)
      return
      end



*     *************************************
*     *                                   *
*     *        Parallel_index_1dblock     *
*     *                                   *
*     *************************************

      integer function Parallel_index_1dblock(m,mb,i)
      implicit none
      integer m,mb,i
      integer ms,r

      if (i.ge.mb) then
         ms = m
      else
         ms = (m/mb)*i
         r = mod(m,mb)
         if (i.lt.r) then
            ms = ms + i
         else
            ms = ms + r
         end if
      end if
      Parallel_index_1dblock = ms
      return
      end


*     *******************************************
*     *                                         *
*     *        Parallel_matrixblocking          *
*     *                                         *
*     *******************************************
*
*   This routine computes mp,np such that
*   mp*np = nthr and min |(mp/np)-(m/n)|
*
*   The justification for trying to keep the ratio of mp/np ~= m/n is that having
*  the subblocks mimimic the overall matrix will produce distribution with blocks are 
*  relatively equal in size.
*
      subroutine Parallel_matrixblocking(nthr,m,n,mb,nb)
      implicit none
      integer nthr,m,n,mb,nb

*     **** local variables ****
      integer ii,jj,mm,nn
      real*8 ratio

      ratio = dble(m)/dble(n)
      mb = nthr
      nb = 1
      do nn =1,nthr
         ii = (nthr-1)/nn - 1
         jj = (nthr+1)/nn + 1
         if (ii.lt.1) ii = 1
         if (jj.gt.nthr) jj=nthr
         do mm=ii,jj
            if ((nn*mm).eq.nthr) then
               if (dabs(dble(mm)/dble(nn) - ratio) .lt.
     >             dabs(dble(mb)/dble(nb) - ratio)) then
                  mb = mm
                  nb = nn
               end if
            end if
         end do
      end do
      return
      end




