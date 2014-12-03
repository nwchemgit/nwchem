*
*     ***********************************
*     *					*
*     *	        K1dB_SumAll		*	
*     *					*
*     ***********************************
      subroutine K1dB_SumAll(sum)
      implicit none
      real*8  sum

#include "K1dB.fh"

#ifdef MPI4
#include "stupid_mpi4.fh"
      real*8  sumall
#else
#include "mpif.h"
      integer msglen,mpierr
      real*8  sumall

*     **** external functions ****
      integer  Parallel3d_comm_k
      external Parallel3d_comm_k
#endif

      if (np_k.gt.1) then

#ifdef MPI4
        stupid_msglen = 1
        call MPI_Allreduce(sum,sumall,stupid_msglen,stupid_double,
     >                     stupid_sum,stupid_comm_k,stupid_ierr)
#else
        msglen = 1
        call MPI_Allreduce(sum,sumall,msglen,
     >                      MPI_DOUBLE_PRECISION,
     >                      MPI_SUM,Parallel3d_comm_k(),mpierr)
#endif

        sum = sumall
      end if
      return
      end


*     ***********************************
*     *					*
*     *	        K1dB_ISumAll		*	
*     *					*
*     ***********************************

      subroutine K1dB_ISumAll(sum)
      implicit none
      integer  sum

#include "K1dB.fh"

#ifdef MPI4
#include "stupid_mpi4.fh"
#else
#include "mpif.h"
      integer msglen,mpierr
      integer sumall

*     **** external functions ****
      integer  Parallel3d_comm_k
      external Parallel3d_comm_k
#endif

      if (np_k.gt.1) then
#ifdef MPI4
        stupid_msglen = 1
        stupid_tsum   = sum
        call MPI_Allreduce(stupid_tsum,stupid_tsum2,stupid_msglen,
     >                     stupid_integer,
     >                     stupid_sum,stupid_comm_k,stupid_ierr)
        sum = stupid_tsum2
#else
        msglen = 1
        call MPI_Allreduce(sum,sumall,msglen,MPI_INTEGER,
     >                       MPI_SUM,Parallel3d_comm_k(),mpierr)
        sum = sumall
#endif
      end if
      return
      end


*     ***********************************
*     *                                 *
*     *         K1dB_MaxAll             *
*     *                                 *
*     ***********************************
      subroutine K1dB_MaxAll(sum)
c     implicit none
      real*8  sum

#include "K1dB.fh"

#ifdef MPI4
#include "stupid_mpi4.fh"
      real*8  sumall
#else
#include "mpif.h"
      integer msglen,mpierr,np
      real*8  sumall

*     **** external functions ****
      integer  Parallel3d_comm_k
      external Parallel3d_comm_k
#endif

      if (np_k.gt.1) then
#ifdef MPI4
        stupid_msglen = 1
        call MPI_Allreduce(sum,sumall,stupid_msglen,stupid_double,
     >                     stupid_max,stupid_comm_k,stupid_ierr)
#else
        msglen = 1
        call MPI_Allreduce(sum,sumall,msglen,MPI_DOUBLE_PRECISION,
     >                      MPI_MAX,Parallel3d_comm_k(),mpierr)
#endif
        sum = sumall
      end if
      return
      end


*     ***********************************
*     *                                 *
*     *         K1dB_MinAll             *
*     *                                 *
*     ***********************************
      subroutine K1dB_MinAll(sum)
c     implicit none
      real*8  sum

#include "K1dB.fh"

#ifdef MPI4
#include "stupid_mpi4.fh"
      real*8  sumall
#else
#include "mpif.h"
      integer msglen,mpierr,np
      real*8  sumall

*     **** external functions ****
      integer  Parallel3d_comm_k
      external Parallel3d_comm_k
#endif

      if (np_k.gt.1) then
#ifdef MPI4
        stupid_msglen = 1
        call MPI_Allreduce(sum,sumall,stupid_msglen,stupid_double,
     >                     stupid_min,stupid_comm_k,stupid_ierr)
#else
        msglen = 1
        call MPI_Allreduce(sum,sumall,msglen,MPI_DOUBLE_PRECISION,
     >                      MPI_MIN,Parallel3d_comm_k(),mpierr)
#endif
        sum = sumall
      end if
      return
      end





*     ***********************************
*     *					*
*     *	        K1dB_Vector_SumAll	*	
*     *					*
*     ***********************************

      subroutine K1dB_Vector_SumAll(n,sum)
c     implicit none
      integer n
      real*8  sum(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "K1dB.fh"

#ifdef MPI4
#include "stupid_mpi4.fh"
      logical value
      integer sumall(2)
#else
#include "mpif.h"
      logical value
      integer msglen
      integer sumall(2),mpierr

*     **** external functions ****
      integer  Parallel3d_comm_k
      external Parallel3d_comm_k

#endif
      call nwpw_timing_start(2)
      if (np_k.gt.1) then

*        ***** allocate temporary space ****
         value = BA_push_get(mt_dbl,n,'sumall',sumall(2),sumall(1))
         if (.not. value) call errquit('out of stack memory',0, MA_ERR)
#ifdef MPI4
         stupid_msglen = n
         call MPI_Allreduce(sum,dbl_mb(sumall(1)),
     >                stupid_msglen,
     >                stupid_double,
     >                stupid_sum,stupid_comm_k,stupid_ierr)
#else
         call MPI_Allreduce(sum,dbl_mb(sumall(1)),n,
     >                MPI_DOUBLE_PRECISION,
     >                MPI_SUM,Parallel3d_comm_k(),mpierr)
#endif
         call dcopy(n,dbl_mb(sumall(1)),1,sum,1)
         value = BA_pop_stack(sumall(2))
         if (.not. value) call errquit('pop stack',0,MA_ERR)
      end if
      call nwpw_timing_end(2)
      return
      end


*     ***********************************
*     *					*
*     *	        K1dB_Vector_ISumAll	*	
*     *					*
*     ***********************************

      subroutine K1dB_Vector_ISumAll(n,sum)
c     implicit none
      integer n
      integer  sum(*)


#include "bafdecls.fh"
#include "errquit.fh"
#include "K1dB.fh"

#ifdef MPI4
#include "stupid_mpi4.fh"

      logical value
      integer sumall(2),sumall1(2)
#else
#include "mpif.h"
      logical value
      integer sumall(2),mpierr

*     **** external functions ****
      integer  Parallel3d_comm_k
      external Parallel3d_comm_k
#endif

      call nwpw_timing_start(2)
      if (np_k.gt.1) then

*       ***** allocate temporary space ****
        value = BA_push_get(mt_int,n,'sumall',sumall(2),sumall(1))
        if (.not. value) call errquit('out of stack memory',0, MA_ERR)

#ifdef MPI4
        if (.not.BA_push_get(mt_int,n,'sumall1',sumall1(2),sumall1(1)))
     >     call errquit('out of stack memory',0, MA_ERR)
        stupid_msglen = n
        call stupid_icopy8to4(n,sum,int_mb(sumall1(1)))
        call MPI_Allreduce(int_mb(sumall1(1)),
     >                     int_mb(sumall(1)),stupid_msglen,
     >                stupid_integer,
     >                stupid_sum,stupid_comm_k,stupid_ierr)
        call stupid_icopy4to8(n,int_mb(sumall(1)),sum)
        if (.not.BA_pop_stack(sumall1(2)))
     >     call errquit('error popping stack',0, MA_ERR)
#else
        call MPI_Allreduce(sum,int_mb(sumall(1)),n,
     >                MPI_INTEGER,
     >                MPI_SUM,Parallel3d_comm_k(),mpierr)
        call icopy(n,int_mb(sumall(1)),1,sum,1)
#endif
        value = BA_pop_stack(sumall(2))
        if (.not. value) call errquit('error popping stack',0, MA_ERR)
      end if
      call nwpw_timing_end(2)

      return
      end


c $Id$
