
#define NBLOCKS 4


*
*     ***********************************
*     *					*
*     *	        D1dB_SumAll		*	
*     *					*
*     ***********************************

      subroutine D1dB_SumAll(sum)
c     implicit none
      real*8  sum

#include "D1dB.fh"

#ifdef MPI4
#include "stupid_mpi4.fh"
#else
#include "mpif.h"
#endif


      integer msglen,mpierr,np
      real*8  sumall


*     **** external functions ****
      integer  Parallel2d_comm_j
      external Parallel2d_comm_j
    
      call Parallel2d_np_j(np)
      if (np.gt.1) then
#ifdef MPI4
        stupid_msglen = 1
        call MPI_Allreduce(sum,sumall,stupid_msglen,stupid_double,
     >                     stupid_sum,stupid_comm_j,stupid_ierr)
#else
        msglen = 1
        call MPI_Allreduce(sum,sumall,msglen,MPI_DOUBLE_PRECISION,
     >                      MPI_SUM,Parallel2d_comm_j(),mpierr)
#endif
        sum = sumall
      end if

      return
      end


*     ***********************************
*     *					*
*     *	        D1dB_ISumAll		*	
*     *					*
*     ***********************************

      subroutine D1dB_ISumAll(sum)
c     implicit none
      integer  sum


#include "D1dB.fh"

#ifdef MPI4
#include "stupid_mpi4.fh"
#else
#include "mpif.h"
#endif


      integer msglen,mpierr,np
      integer sumall

*     **** external functions ****
      integer  Parallel2d_comm_j
      external Parallel2d_comm_j

      call Parallel2d_np_j(np)
      if (np.gt.1) then
#ifdef MPI4
        stupid_msglen = 1
        stupid_tsum = sum
        call MPI_Allreduce(stupid_tsum,stupid_tsum2,stupid_msglen,
     >                     stupid_integer,
     >                     stupid_sum,stupid_comm_j,stupid_ierr)
        sum = stupid_tsum2
#else
        msglen = 1
        call MPI_Allreduce(sum,sumall,msglen,MPI_INTEGER,
     >                       MPI_SUM,Parallel2d_comm_j(),mpierr)
        sum = sumall
#endif
      end if

      return
      end


*     ***********************************
*     *                                 *
*     *         D1dB_MaxAll             *
*     *                                 *
*     ***********************************
      subroutine D1dB_MaxAll(sum)
c     implicit none
      real*8  sum

#include "D1dB.fh"

#ifdef MPI4
#include "stupid_mpi4.fh"
#else
#include "mpif.h"
#endif


      integer msglen,mpierr,np
      real*8  sumall

*     **** external functions ****
      integer  Parallel2d_comm_j
      external Parallel2d_comm_j

      call Parallel2d_np_j(np)
      if (np.gt.1) then
#ifdef MPI4
        stupid_msglen = 1
        call MPI_Allreduce(sum,sumall,stupid_msglen,stupid_double,
     >                      stupid_max,stupid_comm_j,stupid_ierr)
#else
        msglen = 1
        call MPI_Allreduce(sum,sumall,msglen,MPI_DOUBLE_PRECISION,
     >                      MPI_MAX,Parallel2d_comm_j(),mpierr)
#endif
        sum = sumall
      end if

      return
      end


*     ***********************************
*     *                                 *
*     *         D1dB_iMaxAll            *
*     *                                 *
*     ***********************************
      subroutine D1dB_iMaxAll(isum)
c     implicit none
      integer  isum

#include "D1dB.fh"

#ifdef MPI4
#include "stupid_mpi4.fh"
#else
#include "mpif.h"
#endif


      integer msglen,mpierr,np
      integer isumall,iisum

*     **** external functions ****
      integer  Parallel2d_comm_j
      external Parallel2d_comm_j

      call Parallel2d_np_j(np)
      if (np.gt.1) then
#ifdef MPI4
        stupid_msglen = 1
        call stupid_icopy8to4(1,isum,iisum)
        call MPI_Allreduce(iisum,isumall,stupid_msglen,stupid_integer,
     >                      stupid_max,stupid_comm_j,stupid_ierr)
        call stupid_icopy4to8(1,isumall,isum)
#else
        msglen = 1
        call MPI_Allreduce(isum,isumall,msglen,MPI_INTEGER,
     >                      MPI_MAX,Parallel2d_comm_j(),mpierr)
        isum = isumall
#endif
      end if

      return
      end





*     ***********************************
*     *					*
*     *	        D1dB_Vector_SumAll	*	
*     *					*
*     ***********************************

      subroutine D1dB_Vector_SumAll(n,sum)
c     implicit none
      integer n
      real*8  sum(*)

#include "mafdecls.fh"
#include "errquit.fh"
#include "D1dB.fh"

#ifdef MPI4
#include "stupid_mpi4.fh"
#else
#include "mpif.h"
#endif

      logical value
      integer msglen
      integer sumall(2),np,mpierr

*     **** external functions ****
      integer  Parallel2d_comm_j
      external Parallel2d_comm_j

      call Parallel2d_np_j(np)
      call nwpw_timing_start(2)
      if (np.gt.1) then

*     ***** allocate temporary space ****
      value = MA_push_get(mt_dbl,n,'sumall',sumall(2),sumall(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)

#ifdef MPI4
      stupid_msglen = n
      call MPI_Allreduce(sum,dbl_mb(sumall(1)),stupid_msglen,
     >                stupid_double,
     >                stupid_sum,stupid_comm_j,stupid_ierr)
#else
      msglen = n
      call MPI_Allreduce(sum,dbl_mb(sumall(1)),msglen,
     >                MPI_DOUBLE_PRECISION,
     >                MPI_SUM,Parallel2d_comm_j(),mpierr)
#endif

      call dcopy(n,dbl_mb(sumall(1)),1,sum,1)
      value = MA_pop_stack(sumall(2))

      end if
      call nwpw_timing_end(2)
      return
      end


*     ***********************************
*     *					*
*     *	        D1dB_Vector_ISumAll	*	
*     *					*
*     ***********************************

      subroutine D1dB_Vector_ISumAll(n,sum)
c     implicit none
      integer n
      integer  sum(*)


#include "mafdecls.fh"
#include "errquit.fh"
#include "D1dB.fh"

#ifdef MPI4
#include "stupid_mpi4.fh"
#else
#include "mpif.h"
#endif


      logical value
      integer msglen,sumall(2),mpierr,np

*     **** external functions ****
      integer  Parallel2d_comm_j
      external Parallel2d_comm_j

      call Parallel2d_np_j(np)
      call nwpw_timing_start(2)
      if (np.gt.1) then

*       ***** allocate temporary space ****
        value = MA_push_get(mt_int,n,'sumall',sumall(2),sumall(1))
        if (.not. value) call errquit('out of stack memory',0, MA_ERR)


#ifdef MPI4
        stupid_msglen = n
        call stupid_icopy8to4(n,sum,sum)
        call MPI_Allreduce(sum,int_mb(sumall(1)),stupid_msglen,
     >                stupid_integer,
     >                stupid_sum,stupid_comm_j,stupid_ierr)
      call stupid_icopy4to8(n,int_mb(sumall(1)),sum)
#else
        msglen = n
        call MPI_Allreduce(sum,int_mb(sumall(1)),msglen,
     >                MPI_INTEGER,
     >                MPI_SUM,Parallel2d_comm_j(),mpierr)
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
*     *         D1dB_start_rot          *
*     *                                 *
*     ***********************************
*
* This routine start sends A buffer to proc_to = mod(taskid_j+j,np_j)
* and receives W buffer from proc_from = mod(taskid_j-j+np_j,np_j)
* where taskid_j and np_j is the taskid and number of processors for
* the Parallel2d_comm_j communicator
*
      subroutine D1dB_start_rot(j,
     >                          A,W,lda,na,
     >                          request)
      implicit none
      integer j
      real*8 A(*),W(*)
      integer lda,na(*)
      integer request(*)

#include "D1dB.fh"

#include "mpif.h"
#ifdef MPI4
#include "stupid_mpi4.fh"
#endif

*     **** local variables ****
      integer amsglen,wmsglen
      integer proc_to,proc_from,msgtype,mpierr
      integer taskid_j

*     **** external functions ****
      integer  Parallel2d_comm_j
      external Parallel2d_comm_j

      call Parallel2d_taskid_j(taskid_j)

      proc_to   = mod(taskid_j+j,np_j)
      proc_from = mod(taskid_j-j+np_j,np_j)
      msgtype   = j
      amsglen = lda*na(taskid_j+1)
      wmsglen = lda*na(proc_from+1)

#ifdef MPI4
            if (wmsglen.gt.0) then
               stupid_msglen = wmsglen
               stupid_type   = msgtype
               stupid_taskid = proc_from
               call MPI_IRECV(W,
     >                    stupid_msglen,stupid_double,
     >                    stupid_taskid,
     >                    stupid_type,stupid_comm_j,
     >                    stupid_request,stupid_ierr)
               request(1) = stupid_request
               request(3) = 1
            else
               request(3) = 0
            end if

            if (amsglen.gt.0) then
               stupid_msglen = amsglen
               stupid_type   = msgtype
               stupid_taskid = proc_to
               call MPI_ISEND(A,
     >                     stupid_msglen,stupid_double,
     >                     stupid_taskid,
     >                     stupid_type,stupid_comm_j,
     >                     stupid_request,stupid_ierr)
               request(2) = stupid_request
               request(4) = 1
            else
               request(4) = 0
            end if
#else
            if (wmsglen.gt.0) then
               call MPI_IRECV(W,wmsglen,MPI_DOUBLE_PRECISION,
     >                    proc_from,
     >                    msgtype,Parallel2d_comm_j(),
     >                    request(1),mpierr)
               request(3) = 1
            else
               request(3) = 0
            end if
            if (amsglen.gt.0) then
               call MPI_ISEND(A,amsglen,MPI_DOUBLE_PRECISION,
     >                     proc_to,
     >                     msgtype,Parallel2d_comm_j(),
     >                     request(2),mpierr)
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
*     *         D1dB_end_rot            *
*     *                                 *
*     ***********************************
*
*  This routine waits for the sends and receives to 
* finish that were started with D1dB_start_rot routine.
* 
      subroutine D1dB_end_rot(request)
      implicit none
      integer request(*)

#include "mpif.h"
#ifdef MPI4
#include "stupid_mpi4.fh"
#endif

*     **** local variables ****
      integer mpierr
      integer status(2*MPI_STATUS_SIZE)

*     **** wait for completion of mp_send, also do a sync ****
      if (request(3).eq.1) then
#ifdef MPI4
      call stupid_icopy8to4(2,request,request)
      stupid_msglen = 2
      call MPI_WAITALL(stupid_msglen,request,status,stupid_ierr)
#else
      call MPI_WAITALL(2,request,status,mpierr)
#endif


      else if (request(3).eq.2) then
#ifdef MPI4
      call stupid_icopy8to4(1,request,request)
      stupid_msglen = 1
      call MPI_WAITALL(stupid_msglen,request,status,stupid_ierr)
#else
      call MPI_WAITALL(1,request,status,mpierr)
#endif


      else if (request(3).eq.3) then
#ifdef MPI4
      call stupid_icopy8to4(1,request(2),request(2))
      stupid_msglen = 1
      call MPI_WAITALL(stupid_msglen,request(2),status,stupid_ierr)
#else
      call MPI_WAITALL(1,request(2),status,mpierr)
#endif
      endif

      return
      end

