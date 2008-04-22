
*
*     ***********************************
*     *					*
*     *	        K1dB_SumAll		*	
*     *					*
*     ***********************************

      subroutine K1dB_SumAll(sum)
c     implicit none
      real*8  sum

#include "mpif.h"
#include "K1dB.fh"

      integer msglen,mpierr,np
      real*8  sumall


*     **** external functions ****
      integer  Parallel3d_comm_k
      external Parallel3d_comm_k
    
      call Parallel3d_np_k(np)
      if (np.gt.1) then
        msglen = 1
        call MPI_Allreduce(sum,sumall,msglen,MPI_DOUBLE_PRECISION,
     >                      MPI_SUM,Parallel3d_comm_k(),mpierr)
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
c     implicit none
      integer  sum


#include "mpif.h"
#include "K1dB.fh"

      integer msglen,mpierr,np
      integer sumall

*     **** external functions ****
      integer  Parallel3d_comm_k
      external Parallel3d_comm_k

      call Parallel3d_np_k(np)
      if (np.gt.1) then
        msglen = 1

        call MPI_Allreduce(sum,sumall,msglen,MPI_INTEGER,
     >                       MPI_SUM,Parallel3d_comm_k(),mpierr)
        sum = sumall
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

#include "mpif.h"
#include "K1dB.fh"

      integer msglen,mpierr,np
      real*8  sumall

*     **** external functions ****
      integer  Parallel3d_comm_k
      external Parallel3d_comm_k

      call Parallel3d_np_k(np)
      if (np.gt.1) then
        msglen = 1
        call MPI_Allreduce(sum,sumall,msglen,MPI_DOUBLE_PRECISION,
     >                      MPI_MAX,Parallel3d_comm_k(),mpierr)
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

#include "mafdecls.fh"
#include "mpif.h"
#include "errquit.fh"
#include "K1dB.fh"


      logical value
      integer msglen
      integer sumall(2),np,mpierr

*     **** external functions ****
      integer  Parallel3d_comm_k
      external Parallel3d_comm_k

      call Parallel3d_np_k(np)
      call nwpw_timing_start(2)
      if (np.gt.1) then

*     ***** allocate temporary space ****
      value = MA_push_get(mt_dbl,n,'sumall',sumall(2),sumall(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)

      msglen = n


      call MPI_Allreduce(sum,dbl_mb(sumall(1)),msglen,
     >                MPI_DOUBLE_PRECISION,
     >                MPI_SUM,Parallel3d_comm_k(),mpierr)


      call dcopy(n,dbl_mb(sumall(1)),1,sum,1)
      value = MA_pop_stack(sumall(2))

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


#include "mafdecls.fh"
#include "mpif.h"
#include "errquit.fh"
#include "K1dB.fh"

      logical value
      integer msglen,sumall(2),mpierr,np

*     **** external functions ****
      integer  Parallel3d_comm_k
      external Parallel3d_comm_k

      call Parallel3d_np_k(np)
      call nwpw_timing_start(2)
      if (np.gt.1) then

*       ***** allocate temporary space ****
        value = MA_push_get(mt_int,n,'sumall',sumall(2),sumall(1))
        if (.not. value) call errquit('out of stack memory',0, MA_ERR)

        msglen = n

        call MPI_Allreduce(sum,int_mb(sumall(1)),msglen,
     >                MPI_INTEGER,
     >                MPI_SUM,Parallel3d_comm_k(),mpierr)


        call icopy(n,int_mb(sumall(1)),1,sum,1)
        value = MA_pop_stack(sumall(2))
        if (.not. value) call errquit('error popping stack',0, MA_ERR)

      end if
      call nwpw_timing_end(2)

      return
      end


