
#define NBLOCKS 4


*
*     ***********************************
*     *					*
*     *	        D1dBs_SumAll		*	
*     *					*
*     ***********************************

      subroutine D1dBs_SumAll(sum)
c     implicit none
      real  sum

#include "D1dB.fh"

#ifdef MPI4
#include "stupid_mpi4.fh"
#else
#include "mpif.h"
#endif


      integer msglen,mpierr,np
      real  sumall


*     **** external functions ****
      integer  Parallel2d_comm_j
      external Parallel2d_comm_j
    
      call Parallel2d_np_j(np)
      if (np.gt.1) then
!$OMP MASTER
#ifdef MPI4
        stupid_msglen = 1
        call MPI_Allreduce(sum,sumall,stupid_msglen,stupid_real,
     >                     stupid_sum,stupid_comm_j,stupid_ierr)
#else
        msglen = 1
        call MPI_Allreduce(sum,sumall,msglen,MPI_REAL,
     >                      MPI_SUM,Parallel2d_comm_j(),mpierr)
#endif
        sum = sumall
!$OMP END MASTER
!$OMP BARRIER
      end if

      return
      end



*     ***********************************
*     *                                 *
*     *         D1dBs_MaxAll            *
*     *                                 *
*     ***********************************
      subroutine D1dBs_MaxAll(sum)
c     implicit none
      real  sum

#include "D1dB.fh"

#ifdef MPI4
#include "stupid_mpi4.fh"
#else
#include "mpif.h"
#endif


      integer msglen,mpierr,np
      real  sumall

*     **** external functions ****
      integer  Parallel2d_comm_j
      external Parallel2d_comm_j

      call Parallel2d_np_j(np)
      if (np.gt.1) then
#ifdef MPI4
        stupid_msglen = 1
        call MPI_Allreduce(sum,sumall,stupid_msglen,stupid_real,
     >                      stupid_max,stupid_comm_j,stupid_ierr)
#else
        msglen = 1
        call MPI_Allreduce(sum,sumall,msglen,MPI_REAL,
     >                      MPI_MAX,Parallel2d_comm_j(),mpierr)
#endif
        sum = sumall
      end if

      return
      end





*     ***********************************
*     *					*
*     *	       D1dBs_Vector_SumAll	*	
*     *					*
*     ***********************************

      subroutine D1dBs_Vector_SumAll(n,sum)
c     implicit none
      integer n
      real  sum(*)

#include "bafdecls.fh"
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
      integer tid,nthr,nn,nnmod,istart

*     **** external functions ****
      integer  Parallel2d_comm_j
      external Parallel2d_comm_j

      call Parallel2d_np_j(np)
      call nwpw_timing_start(2)
      if (np.gt.1) then

*     ***** allocate temporary space ****
      value = BA_push_get(mt_real,n,'sumall',sumall(2),sumall(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)

!$OMP MASTER
#ifdef MPI4
      stupid_msglen = n
      call MPI_Allreduce(sum,real_mb(sumall(1)),stupid_msglen,
     >                stupid_real,
     >                stupid_sum,stupid_comm_j,stupid_ierr)
#else
      msglen = n
      call MPI_Allreduce(sum,real_mb(sumall(1)),msglen,
     >                MPI_REAL,
     >                MPI_SUM,Parallel2d_comm_j(),mpierr)
#endif
      call Parallel_shared_vector_scopy(.false.,n,
     >                                  real_mb(sumall(1)),sum)
!$OMP END MASTER
!$OMP BARRIER

      value = BA_pop_stack(sumall(2))

      end if
      call nwpw_timing_end(2)
      return
      end


*     ***********************************
*     *                                 *
*     *      D1dBs_Brdcst_values        *
*     *                                 *
*     ***********************************

      subroutine D1dBs_Brdcst_values(psend,nsize,sum)
      implicit none
      integer psend,nsize
      real    sum(*)

#ifdef MPI4
#include "stupid_mpi4.fh"
#else
#include "mpif.h"
#endif


#ifdef MPI4
      integer*4 tpsend
      integer np

      call Parallel2d_np_j(np)
      if (np.gt.1) then
!$OMP MASTER
         stupid_msglen = nsize
         tpsend        = psend
         call MPI_Bcast(sum,stupid_msglen,stupid_real,
     >                  tpsend,stupid_comm_j,stupid_ierr)
!$OMP END MASTER
!$OMP BARRIER
      end if
#else
*     **** external functions ****
      integer  Parallel2d_comm_j
      external Parallel2d_comm_j

      integer ierr,np

      call Parallel2d_np_j(np)
      if (np.gt.1) then
!$OMP MASTER
         call MPI_Bcast(sum,nsize,MPI_REAL,
     >                  psend,Parallel2d_comm_j(),ierr)
!$OMP END MASTER
!$OMP BARRIER
      end if
#endif

      return
      end




*     ***********************************
*     *                                 *
*     *         D1dBs_start_rot         *
*     *                                 *
*     ***********************************
*
* This routine start sends A buffer to proc_to = mod(taskid_j+j,np_j)
* and receives W buffer from proc_from = mod(taskid_j-j+np_j,np_j)
* where taskid_j and np_j is the taskid and number of processors for
* the Parallel2d_comm_j communicator
*
      subroutine D1dBs_start_rot(j,
     >                           A,W,lda,na,
     >                           request)
      implicit none
      integer j
      real     A(*),W(*)
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
     >                    stupid_msglen,stupid_real,
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
     >                     stupid_msglen,stupid_real,
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
               call MPI_IRECV(W,wmsglen,MPI_REAL,
     >                    proc_from,
     >                    msgtype,Parallel2d_comm_j(),
     >                    request(1),mpierr)
               request(3) = 1
            else
               request(3) = 0
            end if
            if (amsglen.gt.0) then
               call MPI_ISEND(A,amsglen,MPI_REAL,
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
*     *         D1dBs_end_rot           *
*     *                                 *
*     ***********************************
*
*  This routine waits for the sends and receives to 
* finish that were started with D1dBs_start_rot routine.
* 
      subroutine D1dBs_end_rot(request)
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
*     *         D1dBs_isendrecv         *
*     *                                 *
*     ***********************************
      subroutine D1dBs_isendrecv(pto,  ssize,sdata,
     >                           pfrom,rsize,rdata,
     >                           request,reqcnt)
      implicit none
      integer pto,ssize
      real    sdata(*)
      integer pfrom,rsize
      real    rdata(*)
      integer request(*)
      integer reqcnt

#include "D1dB.fh"

#ifdef MPI4
#include "stupid_mpi4.fh"
#else
#include "mpif.h"
#endif

*     **** local variables ****
      integer msgtype,mpierr

*     **** external functions ****
      integer  Parallel2d_comm_j
      external Parallel2d_comm_j


      msgtype   = 7

#ifdef MPI4
            if (rsize.gt.0) then
               stupid_msglen = rsize
               stupid_type   = msgtype
               stupid_taskid = pfrom
               call MPI_IRECV(rdata,
     >                    stupid_msglen,stupid_real,
     >                    stupid_taskid,
     >                    stupid_type,stupid_comm_j,
     >                    stupid_request,stupid_ierr)
               reqcnt          = reqcnt + 1
               request(reqcnt) = stupid_request
            end if

            if (ssize.gt.0) then
               stupid_msglen = ssize
               stupid_type   = msgtype
               stupid_taskid = pto
               call MPI_ISEND(sdata,
     >                     stupid_msglen,stupid_real,
     >                     stupid_taskid,
     >                     stupid_type,stupid_comm_j,
     >                     stupid_request,stupid_ierr)
               reqcnt          = reqcnt + 1
               request(reqcnt) = stupid_request
            end if
#else
            if (rsize.gt.0) then
               reqcnt = reqcnt + 1
               call MPI_IRECV(rdata,rsize,MPI_REAL,
     >                    pfrom,
     >                    msgtype,Parallel2d_comm_j(),
     >                    request(reqcnt),mpierr)
            end if
            if (ssize.gt.0) then
               reqcnt = reqcnt + 1
               call MPI_ISEND(sdata,ssize,MPI_REAL,
     >                     pto,
     >                     msgtype,Parallel2d_comm_j(),
     >                     request(reqcnt),mpierr)
            end if
#endif

      return
      end


c     ****************************************
c     *                                      *
c     *          D1dBs_Brdcst_step           *
c     *                                      *
c     ****************************************
c
c  This routine performs step l of a butterfly Broadcast all algorithm. The step
c  l spans from 0..(Level-1) where the number of levels is Level = Log(np_j)/Log(2).
c
c   Entry - l: Butterfly step 0...(Level-1)
c           na: an array of length np_j containing the number of orbitals per taskid_j
c           blocks0: number of blocks to send size=blocks0, 
c                    the exceptions are:
c                    if blocks0==0: the block size is size=2**l.
c                    if blocks0==-1: block size is size=(np_j-2**Level)/2 + 1 for l==(Level-1),
c                                    blocksize is size=2**l otherwise
c           n2ft3d: leading size of psi_rep
c           psi_rep: data array
c   Exit -
c           psi_rep: modified data array 
c           requests,reqcnt: tags for asychronous message passing
c
      subroutine D1dBs_Brdcst_step(l,na,blocks0,
     >                             n2ft3d,psi_rep,
     >                             requests,reqcnt)
      implicit none
      integer l,na(*),blocks0
      integer n2ft3d
      real    psi_rep(n2ft3d,*)
      integer requests(*),reqcnt

*     *** local variables ***
      integer taskid_j,np_j
      integer i,pr,ps,shift,size,Level
      integer pto,pfrom,rsize,ssize,rindx,sindx

*     *** local variables ***
      integer  Butter_levels
      external Butter_levels

      call Parallel2d_taskid_j(taskid_j)
      call Parallel2d_np_j(np_j)

cccc      Level = (log(dble(np_j))/log(2.0d0))
      Level = Butter_levels(np_j)

      shift = 2**l
      pfrom = mod(taskid_j     +shift,np_j)
      pto   = mod(taskid_j+np_j-shift,np_j)

*     *** hfx exception ***
      if (blocks0.lt.0) then
         if (l.eq.(Level-1)) then
            size = ((np_j-2**Level)/2) + 1
         else
            size =  shift
         end if
*     *** standard butterfly exception ***
      elseif  (blocks0.eq.0) then
         size = shift
*     *** user blocksize ***
      else
         size = blocks0
      end if

*     *** determine message sizes ***
      rsize = 0
      ssize = 0
      do i=0,(size-1)
        pr = mod(pfrom    + i,np_j)
        ps = mod(taskid_j + i,np_j)
        rsize = rsize + na(pr+1)
        ssize = ssize + na(ps+1)
      end do

*     *** determine message indexes ***
      rindx = 1
      do  i=0,(shift-1)
        ps = mod(taskid_j + i,np_j)
        rindx = rindx + na(ps+1)
      end do
      sindx = 1


*     *** will be much more complicated for synchronous ***
      reqcnt = 0
      call D1dBs_isendrecv(pto,  ssize*n2ft3d,psi_rep(1,sindx),
     >                     pfrom,rsize*n2ft3d,psi_rep(1,rindx),
     >                     requests,reqcnt)
      return
      end


c     ****************************************
c     *                                      *
c     *           D1dBs_Reduce_step          *
c     *                                      *
c     ****************************************
c
c  This routine performs step l of a butterfly Reduceall algorithm. The step
c  l spans from 0..(Level-1) where the number of levels is Level = Log(np_j)/Log(2).
c
c   Entry - l: Butterfly step 0...(Level-1)
c           na: an array of length np_j containing the number of orbitals per taskid_j
c           blocks0: number of blocks to send size=blocks0, 
c                    the exceptions are:
c                    if blocks0==0: the block size is size=2**l.
c                    if blocks0==-1: block size is size=(np_j-2**Level)/2 + 1 for l==(Level-1),
c                                    blocksize is size=2**l otherwise
c           n2ft3d: leading size of psi_rep
c           hpsi_rep: data array
c           tmp: tempory data array. Needs to be at least  n2ft3d*size
c
c   Exit - hpsi_rep: modified data array
c          requests,reqcnt: tags for asychronous message passing
c
      subroutine D1dBs_Reduce_step(l,na,blocks0,
     >                             n2ft3d,hpsi_rep,tmp)
      implicit none
      integer l,na(*),blocks0
      integer n2ft3d
      real    hpsi_rep(n2ft3d,*)
      real    tmp(*)

*     *** local variables ***
      integer taskid_j,np_j
      integer i,pr,ps,size,shift,Level,pfrom,pto
      integer rsize,ssize,rindx,sindx
      integer requests(10),reqcnt

*     *** local variables ***
      integer  Butter_levels
      external Butter_levels

      call Parallel2d_taskid_j(taskid_j)
      call Parallel2d_np_j(np_j)

cccccc      !Level = (log(dble(np_j))/log(2.0d0))
      Level = Butter_levels(np_j)

      shift = 2**l
      pfrom = mod(taskid_j+np_j-shift,np_j)
      pto   = mod(taskid_j     +shift,np_j)

*     *** hfx exception ***
      if (blocks0.lt.0) then
         if (l.eq.(Level-1)) then
            size = ((np_j-2**Level)/2) + 1
         else
            size =  shift
         end if
*     *** standard butterfly exception ***
      elseif  (blocks0.eq.0) then
         size = shift
*     *** user blocksize ***
      else
         size = blocks0
      end if

*     *** determine message sizes ***
      rsize = 0
      ssize = 0
      do i=0,(size-1)
        pr = mod(taskid_j + i,np_j)
        ps = mod(pto      + i,np_j)
        rsize = rsize + na(pr+1)
        ssize = ssize + na(ps+1)
      end do

*     *** determine message indexes ***
      rindx = 1
      sindx = 1
      do i=0,(shift-1)
         ps = mod(taskid_j  + i,np_j)
         sindx = sindx + na(ps+1)
      end do


*     *** will be much more complicated for synchronous ***
      reqcnt = 0
      call D1dBs_isendrecv(pto,  ssize*n2ft3d,hpsi_rep(1,sindx),
     >                     pfrom,rsize*n2ft3d,tmp,
     >                    requests,reqcnt)

      call D1dB_WaitAll(requests,reqcnt)
      call SAXPY_OMP(rsize*n2ft3d,1.0d0,tmp,1,hpsi_rep(1,rindx),1)

      return
      end
c $Id$
