
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


