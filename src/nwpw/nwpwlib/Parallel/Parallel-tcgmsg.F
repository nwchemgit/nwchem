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

c      include 'mpif.h'
c      integer mpierr

*     **** local variables ****
      integer i

      np     = nnodes()
      taskid = nodeid()


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
#include "bafdecls.fh"
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

      if (.not.BA_push_get(mt_int,np,'tmppp2',tmp(2),tmp(1)))
     >  call errquit('Parallel2d_init:out of stack memory',0, MA_ERR)
        do i=0,np_i-1
          int_mb(tmp(1)+i) = int_mb(procNd(1) + i + taskid_j*np_i) 
        end do
        comm_i = ga_pgroup_create(int_mb(tmp(1)),np_i)
        do j=0,np_j-1
          int_mb(tmp(1)+j) = int_mb(procNd(1) + taskid_i + j*np_i) 
        end do
        comm_j = ga_pgroup_create(int_mb(tmp(1)),np_j)
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
#include "bafdecls.fh"
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

      if (.not.BA_push_get(mt_int,np,'tmppp2',tmp(2),tmp(1)))
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
      if (.not.BA_pop_stack(tmp(2)))
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
#include "bafdecls.fh"
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
*     *         Parallel_IMaxAll        *
*     *                                 *
*     ***********************************
      subroutine Parallel_IMaxAll(isum)
c     implicit none
      integer isum

#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "Parallel.fh"

      if (np.gt.1) then
         call GA_IGOP(9+MSGINT,isum,1,'max')
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

#include "bafdecls.fh"
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

#include "bafdecls.fh"
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

#include "bafdecls.fh"
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

#include "bafdecls.fh"
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
*     *      Parallel_Brdcst_ivalue     *
*     *                                 *
*     ***********************************

      subroutine Parallel_Brdcst_ivalue(psend,isum)
      implicit none
      integer psend
      integer isum

#include "bafdecls.fh"
#include "errquit.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "Parallel.fh"

*     **** local variables ****
      integer msglen

      if (np.gt.1) then
         msglen = 1
         call BRDCST(9+MSGINT,isum,mitob(msglen),psend)
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

#include "bafdecls.fh"
#include "errquit.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "Parallel.fh"

      if (np.gt.1) then
         call BRDCST(9+MSGINT,isum,mitob(nsize),psend)
      end if
      return
      end




*     ***********************************
*     *                                 *
*     *         Parallela_MaxAll        *
*     *                                 *
*     ***********************************

      subroutine Parallela_MaxAll(ic,sum)
c     implicit none
      integer ic
      real*8  sum

#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "Parallel.fh"

      if (np.gt.1) then
         if (ic.eq.1) then
            call GA_PGROUP_DGOP(comm_i,9+MSGDBL,sum,1,'max')
         else if (ic.eq.2) then
            call GA_PGROUP_DGOP(comm_j,9+MSGDBL,sum,1,'max')
         else if (ic.eq.3) then
            call GA_PGROUP_DGOP(comm_k,9+MSGDBL,sum,1,'max')
         else
            call GA_DGOP(9+MSGDBL,sum,1,'max')
         end if
      end if

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

#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "Parallel.fh"

      if (np.gt.1) then
         if (ic.eq.1) then
            call GA_PGROUP_DGOP(comm_i,9+MSGDBL,sum,1,'+')
         else if (ic.eq.1) then
            call GA_PGROUP_DGOP(comm_j,9+MSGDBL,sum,1,'+')
         else if (ic.eq.3) then
            call GA_PGROUP_DGOP(comm_k,9+MSGDBL,sum,1,'+')
         else 
            call GA_DGOP(9+MSGDBL,sum,1,'+')
         end if
      end if

      return
      end



*     ***********************************
*     *                                 *
*     *         Parallela_ISumAll       *
*     *                                 *
*     ***********************************

      subroutine Parallela_ISumAll(ic,sum)
c     implicit none
      integer ic
      integer  sum

#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "Parallel.fh"

      if (np.gt.1) then
         if (ic.eq.1) then
            call GA_PGROUP_IGOP(comm_i,9+MSGINT,sum,1,'+')
         else if (ic.eq.2) then
            call GA_PGROUP_IGOP(comm_j,9+MSGINT,sum,1,'+')
         else if (ic.eq.3) then
            call GA_PGROUP_IGOP(comm_k,9+MSGINT,sum,1,'+')
         else 
            call GA_IGOP(9+MSGINT,sum,1,'+')
         end if
      end if
      
      return
      end


*     ***********************************
*     *                                 *
*     *      Parallela_Vector_SumAll    *
*     *                                 *
*     ***********************************
      subroutine Parallela_Vector_SumAll(ic,n,sum)
c     implicit none
      integer ic
      integer n
      real*8  sum(*)

#include "bafdecls.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "errquit.fh"
#include "Parallel.fh"


      call nwpw_timing_start(2)
      if (np.gt.1) then
         if (ic.eq.1) then
            call GA_PGROUP_DGOP(comm_i,9+MSGDBL,sum,n,'+')
         else if (ic.eq.2) then
            call GA_PGROUP_DGOP(comm_j,9+MSGDBL,sum,n,'+')
         else if (ic.eq.3) then
            call GA_PGROUP_DGOP(comm_k,9+MSGDBL,sum,n,'+')
         else
            call GA_DGOP(9+MSGDBL,sum,n,'+')
         end if
      end if
      call nwpw_timing_end(2)
      return
      end



*     ***********************************
*     *                                 *
*     *     Parallela_Vector_ISumAll    *
*     *                                 *
*     ***********************************

      subroutine Parallela_Vector_ISumAll(ic,n,sum)
c     implicit none
      integer ic
      integer n
      integer  sum(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "Parallel.fh"


      call nwpw_timing_start(2)
      if (np.gt.1) then
         if (ic.eq.1) then
            call GA_PGROUP_IGOP(comm_i,9+MSGINT,sum,n,'+')
         else if (ic.eq.2) then
            call GA_PGROUP_IGOP(comm_j,9+MSGINT,sum,n,'+')
         else if (ic.eq.3) then
            call GA_PGROUP_IGOP(comm_k,9+MSGINT,sum,n,'+')
         else 
            call GA_IGOP(9+MSGINT,sum,n,'+')
         end if
      end if
      call nwpw_timing_end(2)

      return
      end



*     ***********************************
*     *                                 *
*     *      Parallela_Brdcst_value      *
*     *                                 *
*     ***********************************

      subroutine Parallela_Brdcst_value(ic,psend,sum)
      implicit none
      integer ic
      integer psend
      real*8  sum

#include "bafdecls.fh"
#include "errquit.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "Parallel.fh"

*     **** local variables ****
      integer msglen

      if (np.gt.1) then
         msglen = 1
         if (ic.eq.1) then
          call GA_PGROUP_BRDCST(comm_i,9+MSGDBL,sum,mdtob(msglen),psend)
         else if (ic.eq.2) then
          call GA_PGROUP_BRDCST(comm_j,9+MSGDBL,sum,mdtob(msglen),psend)
         else if (ic.eq.3) then
          call GA_PGROUP_BRDCST(comm_k,9+MSGDBL,sum,mdtob(msglen),psend)
         else
            call GA_BRDCST(9+MSGDBL,sum,mdtob(msglen),psend)
         end if
      end if

      return
      end


*     ***********************************
*     *                                 *
*     *      Parallela_Brdcst_values    *
*     *                                 *
*     ***********************************

      subroutine Parallela_Brdcst_values(ic,psend,nsize,sum)
      implicit none
      integer ic
      integer psend,nsize
      real*8  sum(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "Parallel.fh"


      if (np.gt.1) then
         if (ic.eq.1) then
           call GA_PGROUP_BRDCST(comm_i,9+MSGDBL,sum,mdtob(nsize),psend)
         else if (ic.eq.2) then
           call GA_PGROUP_BRDCST(comm_j,9+MSGDBL,sum,mdtob(nsize),psend)
         else if (ic.eq.3) then
           call GA_PGROUP_BRDCST(comm_k,9+MSGDBL,sum,mdtob(nsize),psend)
         else
            call GA_BRDCST(9+MSGDBL,sum,mdtob(nsize),psend)
         end if
      end if

      return
      end



*     ***********************************
*     *                                 *
*     *      Parallela_Brdcst_ivalue    *
*     *                                 *
*     ***********************************
      subroutine Parallela_Brdcst_ivalue(ic,psend,isum)
      implicit none
      integer ic
      integer psend
      integer isum

#include "bafdecls.fh"
#include "errquit.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "Parallel.fh"

*     **** local variables ****
      integer msglen

      if (np.gt.1) then
         msglen = 1
         if (ic.eq.1) then
         call GA_PGROUP_BRDCST(comm_i,9+MSGINT,isum,mitob(msglen),psend)
         else if (ic.eq.2) then
         call GA_PGROUP_BRDCST(comm_j,9+MSGINT,isum,mitob(msglen),psend)
         else if (ic.eq.3) then
         call GA_PGROUP_BRDCST(comm_k,9+MSGINT,isum,mitob(msglen),psend)
         else
         call GA_BRDCST(9+MSGINT,isum,mitob(msglen),psend)
         end if
      end if

      return
      end




*     ***********************************
*     *                                 *
*     *      Parallela_Brdcst_ivalues   *
*     *                                 *
*     ***********************************

      subroutine Parallela_Brdcst_ivalues(ic,psend,nsize,isum)
      implicit none
      integer ic
      integer psend,nsize
      integer  isum(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "Parallel.fh"

      if (np.gt.1) then
         if (ic.eq.1) then
          call GA_PGROUP_BRDCST(comm_i,9+MSGINT,isum,mitob(nsize),psend)
         else if (ic.eq.2) then
          call GA_PGROUP_BRDCST(comm_j,9+MSGINT,isum,mitob(nsize),psend)
         else if (ic.eq.3) then
          call GA_PGROUP_BRDCST(comm_k,9+MSGINT,isum,mitob(nsize),psend)
         else 
            call GA_BRDCST(9+MSGINT,isum,mitob(nsize),psend)
         end if
      end if
      return
      end



*     ***********************************
*     *                                 *
*     *      Parallela_start_rotate     *
*     *                                 *
*     ***********************************

      subroutine Parallela_start_rotate(ic,shift,
     >                                  A1,nsize1,
     >                                  A2,nsize2,request)
      implicit none
      integer ic,shift
      real*8  A1(*)
      integer nsize1
      real*8  A2(*)
      integer nsize2
      integer request(*)

#include "bafdecls.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "errquit.fh"

#include "Parallel.fh"


*     **** local variables ****
      integer i,mynp,mytaskid,msglen,msgype,mpierr,proc_to,proc_from
      integer rcv_proc,rcv_len,psr,tmp(2),pto

*     ***** external functions ****
      integer  Parallel3d_convert_taskid_i
      integer  Parallel3d_convert_taskid_j
      integer  Parallel3d_convert_taskid_k
      external Parallel3d_convert_taskid_i
      external Parallel3d_convert_taskid_j
      external Parallel3d_convert_taskid_k

      if (ic.eq.1) then
         mynp     = np_i
         mytaskid = taskid_i
         proc_to   = mod(mytaskid+shift+mynp,mynp)
         proc_from = mod(mytaskid-shift+mynp,mynp)

         proc_to   = Parallel3d_convert_taskid_i(proc_to)
         proc_from = Parallel3d_convert_taskid_i(proc_from)
      else if (ic.eq.2) then
         mynp     = np_j
         mytaskid = taskid_j
         proc_to   = mod(mytaskid+shift+mynp,mynp)
         proc_from = mod(mytaskid-shift+mynp,mynp)

         proc_to   = Parallel3d_convert_taskid_j(proc_to)
         proc_from = Parallel3d_convert_taskid_j(proc_from)
      else if (ic.eq.3) then
         mynp     = np_k
         mytaskid = taskid_k
         proc_to   = mod(mytaskid+shift+mynp,mynp)
         proc_from = mod(mytaskid-shift+mynp,mynp)

         proc_to   = Parallel3d_convert_taskid_k(proc_to)
         proc_from = Parallel3d_convert_taskid_k(proc_from)
      else
         mynp     = np
         mytaskid = taskid
         proc_to   = mod(mytaskid+shift+mynp,mynp)
         proc_from = mod(mytaskid-shift+mynp,mynp)
      end if

*      /* determine psr - should be made w/o using tmp array! */
      if (.not.BA_push_get(mt_int,mynp,'tmppp2',tmp(2),tmp(1)))
     >  call errquit('Parallela_start_rotate:out of stack',0,MA_ERR)
      do i=0,np-1
         int_mb(tmp(1)+i) = 0
      end do
      do i=0,np-1
         pto = mod(i+shift+mynp,mynp)
         if ((int_mb(tmp(1)+i).eq.0).and.(int_mb(tmp(1)+pto).eq.0)) then
            int_mb(tmp(1)+i)   = 1
            int_mb(tmp(1)+pto) = 2
         end if
      end do
      psr = int_mb(tmp(1)+mytaskid)
      if (psr.eq.0) psr = 2
      if (.not.BA_pop_stack(tmp(2)))
     >  call errquit('Parallela_start_rotate:popping stack',0,MA_ERR)

*     **** send then receive ****
      if (psr.eq.1) then

         if (nsize1.gt.0) then
            msglen  = nsize1
            call SND(9+MSGDBL,A1,mdtob(msglen),proc_to,1)
            request(4) = 1
         else
            request(4) = 0
         end if

         if (nsize2.gt.0) then
            msglen  = nsize2
            call RCV(9+MSGDBL,A2,mdtob(msglen),rcv_len,
     >            proc_from,rcv_proc,1)
            request(3) = 1
         else
            request(3) = 0
         end if

*     **** receive then receive ****
      else
         if (nsize2.gt.0) then
            msglen  = nsize2
            call RCV(9+MSGDBL,A2,mdtob(msglen),rcv_len,
     >            proc_from,rcv_proc,1)
            request(3) = 1
         else
            request(3) = 0
         end if

         if (nsize1.gt.0) then
            msglen  = nsize1
            call SND(9+MSGDBL,A1,mdtob(msglen),proc_to,1)
            request(4) = 1
         else
            request(4) = 0
         end if
      end if

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
*     *      Parallela_start_Irotate    *
*     *                                 *
*     ***********************************
      subroutine Parallela_start_Irotate(ic,shift,
     >                                  A1,nsize1,
     >                                  A2,nsize2,request)
      implicit none
      integer ic,shift
      integer A1(*)
      integer nsize1
      integer A2(*)
      integer nsize2
      integer request(*)

#include "bafdecls.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "errquit.fh"

#include "Parallel.fh"


*     **** local variables ****
      integer i,mynp,mytaskid,msglen,msgype,mpierr,proc_to,proc_from
      integer rcv_proc,rcv_len,psr,tmp(2),pto

*     ***** external functions ****
      integer  Parallel3d_convert_taskid_i
      integer  Parallel3d_convert_taskid_j
      integer  Parallel3d_convert_taskid_k
      external Parallel3d_convert_taskid_i
      external Parallel3d_convert_taskid_j
      external Parallel3d_convert_taskid_k

      if (ic.eq.1) then
         mynp     = np_i
         mytaskid = taskid_i
         proc_to   = mod(mytaskid+shift+mynp,mynp)
         proc_from = mod(mytaskid-shift+mynp,mynp)

         proc_to   = Parallel3d_convert_taskid_i(proc_to)
         proc_from = Parallel3d_convert_taskid_i(proc_from)
      else if (ic.eq.2) then
         mynp     = np_j
         mytaskid = taskid_j
         proc_to   = mod(mytaskid+shift+mynp,mynp)
         proc_from = mod(mytaskid-shift+mynp,mynp)

         proc_to   = Parallel3d_convert_taskid_j(proc_to)
         proc_from = Parallel3d_convert_taskid_j(proc_from)
      else if (ic.eq.3) then
         mynp     = np_k
         mytaskid = taskid_k
         proc_to   = mod(mytaskid+shift+mynp,mynp)
         proc_from = mod(mytaskid-shift+mynp,mynp)

         proc_to   = Parallel3d_convert_taskid_k(proc_to)
         proc_from = Parallel3d_convert_taskid_k(proc_from)
      else
         mynp     = np
         mytaskid = taskid
         proc_to   = mod(mytaskid+shift+mynp,mynp)
         proc_from = mod(mytaskid-shift+mynp,mynp)
      end if

*      /* determine psr - should be made w/o using tmp array! */
      if (.not.BA_push_get(mt_int,mynp,'tmppp2',tmp(2),tmp(1)))
     >  call errquit('Parallela_start_rotate:out of stack',0,MA_ERR)
      do i=0,np-1
         int_mb(tmp(1)+i) = 0
      end do
      do i=0,np-1
         pto = mod(i+shift+mynp,mynp)
         if ((int_mb(tmp(1)+i).eq.0).and.(int_mb(tmp(1)+pto).eq.0)) then
            int_mb(tmp(1)+i)   = 1
            int_mb(tmp(1)+pto) = 2
         end if
      end do
      psr = int_mb(tmp(1)+mytaskid)
      if (psr.eq.0) psr = 2
      if (.not.BA_pop_stack(tmp(2)))
     >  call errquit('Parallela_start_rotate:popping stack',0,MA_ERR)

*     **** send then receive ****
      if (psr.eq.1) then

         if (nsize1.gt.0) then
            msglen  = nsize1
            call SND(9+MSGINT,A1,mitob(msglen),proc_to,1)
            request(4) = 1
         else
            request(4) = 0
         end if

         if (nsize2.gt.0) then
            msglen  = nsize2
            call RCV(9+MSGINT,A2,mitob(msglen),rcv_len,
     >            proc_from,rcv_proc,1)
            request(3) = 1
         else
            request(3) = 0
         end if

*     **** receive then receive ****
      else
         if (nsize2.gt.0) then
            msglen  = nsize2
            call RCV(9+MSGINT,A2,mitob(msglen),rcv_len,
     >            proc_from,rcv_proc,1)
            request(3) = 1
         else
            request(3) = 0
         end if

         if (nsize1.gt.0) then
            msglen  = nsize1
            call SND(9+MSGINT,A1,mitob(msglen),proc_to,1)
            request(4) = 1
         else
            request(4) = 0
         end if
      end if

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
      !*** do nothing ***

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

#include "tcgmsg.fh"
#include "msgtypesf.h"

      integer  rcv_len,rcv_proc

      call SND(9+MSGCHR,cval,nsize,pto,1)
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
      character rval(*)                 

#include "tcgmsg.fh"
#include "msgtypesf.h"

      call SND(9+MSGDBL,rval,mdtob(nsize),pto,1)
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
      integer   ival(*)

#include "tcgmsg.fh"
#include "msgtypesf.h"

      call SND(9+MSGINT,ival,mitob(nsize),pto,1)
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

#include "tcgmsg.fh"
#include "msgtypesf.h"

      integer  rcv_len,rcv_proc

      call RCV(9+MSGCHR,cval,nsize,rcv_len,pfrom,rcv_proc,1)
      return
      end


*     ***********************************
*     *                                 *
*     *      Parallel_recv_values       *
*     *                                 *
*     ***********************************
      subroutine Parallel_recv_values(pfrom,msgtype,nsize,rval)
      implicit none
      integer  pfrom,msgtype,nsize
      real*8   rval(*)

#include "tcgmsg.fh"
#include "msgtypesf.h"

      integer  rcv_len,rcv_proc

      call RCV(9+MSGDBL,rval,mdtob(nsize),rcv_len,pfrom,rcv_proc,1)
      return
      end


*     ***********************************
*     *                                 *
*     *      Parallel_recv_ivalues      *
*     *                                 *
*     ***********************************
      subroutine Parallel_recv_ivalues(pfrom,msgtype,nsize,ival)
      implicit none
      integer  pfrom,msgtype,nsize
      integer  ival(*)

#include "tcgmsg.fh"
#include "msgtypesf.h"

      integer  rcv_len,rcv_proc

      call RCV(9+MSGINT,ival,mitob(nsize),rcv_len,pfrom,rcv_proc,1)
      return
      end

