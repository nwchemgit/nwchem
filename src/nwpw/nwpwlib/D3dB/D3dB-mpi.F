#define NBLOCKS 3

*
* $Id$
*

*     ***********************************
*     *					*
*     *	   D3dB_c_transpose_jk		*
*     *					*
*     ***********************************

      subroutine D3dB_c_transpose_jk(nb,A,tmp1,tmp2)

*****************************************************
*                                                   *
*      This routine performs the operation          *
*               A(i,k,j) <- A(i,j,k)                * 
*                                                   *
*      np = the number of worker nodes              *
*      proc#=0...(np-1)
*                                                   *
*       this transpose uses more buffer space       *
*       then transpose2                             *
*****************************************************
      implicit none
      integer nb
      complex*16  A(*)
      complex*16  tmp1(*),tmp2(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "mpif.h"
#include "D3dB.fh"

#ifdef MPI4
#include "stupid_mpi4.fh"
#endif


*     **** indexing variables ****
c     integer iq_to_i1((NFFT1/2+1)*NFFT2*NSLABS)
c     integer iq_to_i2((NFFT1/2+1)*NFFT2*NSLABS)
c     integer i1_start(NFFT3+1)
c     integer i2_start(NFFT3+1)
      integer iq_to_i1(2,NBLOCKS)
      integer iq_to_i2(2,NBLOCKS)
      integer i1_start(2,NBLOCKS)
      integer i2_start(2,NBLOCKS)
      common / trans_blk / iq_to_i1,iq_to_i2,i1_start,i2_start



*     **** local variables ***
      logical value
      integer i
      integer proc_to,proc_from,it
      integer msgtype,source
      integer msglen
      integer request(2),reqcnt
      integer taskid,np
      integer mpierr

*     **** external functions ****
      integer  Parallel2d_convert_taskid_i,Parallel_comm_world
      external Parallel2d_convert_taskid_i,Parallel_comm_world
      
      call Parallel2d_taskid_i(taskid)
      call Parallel2d_np_i(np)

*     **** allocate memory ****
      value = BA_push_get(mt_int,np,
     >                     'request',request(2),request(1))
      if (.not. value)
     > call errquit(' D3dB_c_transpose_jk:out of stack',0,MA_ERR)


*     **** pack A(i) array ****
!$OMP DO
      do i=1,nfft3d(nb) 
         tmp1(int_mb(iq_to_i1(1,nb)+i-1)) = A(i)
      end do
!$OMP END DO

*     **** it = 0, transpose data on same thread ****  
      msglen = int_mb(i2_start(1,nb)+2-1) - int_mb(i2_start(1,nb)+1-1)

c      do i=1,msglen
c         tmp2(int_mb(i2_start(1,nb)+1-1)+i-1) 
c     > = tmp1(int_mb(i1_start(1,nb)+1-1)+i-1)
c      end do
      call Parallel_shared_vector_copy(.false.,2*msglen,
     >               tmp1(int_mb(i1_start(1,nb)+1-1)),
     >               tmp2(int_mb(i2_start(1,nb)+1-1)))
         

*     **** receive packed array data ****
      reqcnt = 0
!$OMP MASTER
      do it = 1,np-1

*        **** synchronous receive of tmp ****
         proc_from = mod((taskid)-it+np, np) 

         msgtype = 9
         source=proc_from
         msglen = (int_mb(i2_start(1,nb)+it+2-1) 
     >          -  int_mb(i2_start(1,nb)+it+1-1))
         

         if (msglen.gt.0) then
            reqcnt = reqcnt + 1
#ifdef MPI4
            stupid_msglen = msglen
            stupid_type   = msgtype
            stupid_taskid = Parallel2d_convert_taskid_i(source)
            call MPI_IRECV(tmp2(int_mb(i2_start(1,nb)+it+1-1)),
     >                    stupid_msglen,stupid_complex,
     >                    stupid_taskid,
     >                    stupid_type,stupid_world,
     >                    stupid_request,stupid_ierr) 
            int_mb(request(1)+reqcnt-1) = 0
            int_mb(request(1)+reqcnt-1) = stupid_request
#else
            call MPI_IRECV(tmp2(int_mb(i2_start(1,nb)+it+1-1)),
     >                    msglen,MPI_DOUBLE_COMPLEX,
     >                    Parallel2d_convert_taskid_i(source),
     >                    msgtype,Parallel_comm_world(),
     >                    int_mb(request(1)+reqcnt-1),mpierr) 
#endif
         end if

      end do

*     **** send packed array to other processors ****
      do it = 1,np-1

*        **** synchronous send of tmp ****
         proc_to   = mod((taskid)+it, np) 
         msglen    = (int_mb(i1_start(1,nb)+it+2-1)
     >              - int_mb(i1_start(1,nb)+it+1-1))
         msgtype   = 9
      
         if (msglen.gt.0) then
#ifdef MPI4
            stupid_msglen = msglen
            stupid_type   = msgtype
            stupid_taskid = Parallel2d_convert_taskid_i(proc_to)
            call MPI_SEND(tmp1(int_mb(i1_start(1,nb)+it+1-1)),
     >                     stupid_msglen,stupid_complex,
     >                     stupid_taskid,
     >                     stupid_type,stupid_world,stupid_ierr)
#else
            call MPI_SEND(tmp1(int_mb(i1_start(1,nb)+it+1-1)),
     >                     msglen,MPI_DOUBLE_COMPLEX,
     >                     Parallel2d_convert_taskid_i(proc_to),
     >                     msgtype,Parallel_comm_world(),mpierr)
#endif
         end if

      end do

   
*     **** wait for completion of mp_send, also do a sync ****
      if (np.gt.1) call Parallel_mpiWaitAll0(reqcnt,int_mb(request(1)))

!$OMP END MASTER
!$OMP BARRIER

*     **** unpack A(i) array ****
!$OMP DO
      do i=1,nfft3d(nb) !***(nx(nb)/2+1)*ny(nb)*nq(nb)
         A(i) = tmp2(int_mb(iq_to_i2(1,nb)+i-1))
      end do
!$OMP END DO
      

      !*** deallocate memory ***
      value =           BA_pop_stack(request(2))
      if (.not. value)
     > call errquit(' D3dB_c_transpose_jk:popping stack',0,MA_ERR)

      return
      end




*     ***********************************
*     *                                 *
*     *    D3dB_c_timereverse           *
*     *                                 *
*     ***********************************

      subroutine D3dB_c_timereverse(nb,A,tmp1,tmp2)

*****************************************************
*                                                   *
*      This routine performs the operation          *
*            A(i,j,k) <- conjugate(A(i,-j,-k))      *
*                                                   *
*      np = the number of worker nodes              *
*      proc#=0...(np-1)                             *
*                                                   *
*****************************************************
      implicit none
      integer     nb
      complex*16  A(*)
      complex*16  tmp1(*)
      complex*16  tmp2(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "mpif.h"
#include "D3dB.fh"

#ifdef MPI4
#include "stupid_mpi4.fh"
#endif


*     **** indexing variables ****
c     integer iq_to_i1(2**NFFT2*NSLABS)
c     integer iq_to_i2(2**NFFT2*NSLABS)
c     integer i1_start(NFFT3+1)
c     integer i2_start(NFFT3+1)
      integer iq_to_i1(2,NBLOCKS)
      integer iq_to_i2(2,NBLOCKS)
      integer i1_start(2,NBLOCKS)
      integer i2_start(2,NBLOCKS)
      common / timereverse_blk / iq_to_i1,iq_to_i2,i1_start,i2_start


*     **** local variables ***
      logical value
      integer it
      integer source
      integer msglen,msgtype
      integer proc_from,proc_to
      integer request(2),reqcnt
      integer taskid,np
      integer index1,index2
      integer mpierr

*     **** external functions ****
      integer  Parallel2d_convert_taskid_i,Parallel_comm_world
      external Parallel2d_convert_taskid_i,Parallel_comm_world

      call Parallel2d_taskid_i(taskid)
      call Parallel2d_np_i(np)

*     **** allocate memory ****
      value = BA_push_get(mt_int,np,
     >                     'request',request(2),request(1))
      if (.not. value)
     > call errquit(' D3dB_c_timereverse:out of stack',0,MA_ERR)


*     **** pack A(i) array ****
#ifndef CRAY
!DIR$ ivdep
#endif
!$OMP DO
      do index1=int_mb(i1_start(1,nb)+1-1),
     >         (int_mb(i1_start(1,nb)+np+1-1)-1)
         tmp1(index1) = A(int_mb(iq_to_i1(1,nb)+index1-1))
      end do
!$OMP END DO

*     **** it = 0, transpose data on same thread ****
      msglen = int_mb(i2_start(1,nb)+2-1) - int_mb(i2_start(1,nb)+1-1)
c      call dcopy(2*msglen,tmp1(int_mb(i1_start(1,nb))),1,
c     >                    tmp2(int_mb(i2_start(1,nb))),1)
      call Parallel_shared_vector_copy(.false.,2*msglen,
     >                                 tmp1(int_mb(i1_start(1,nb))),
     >                                 tmp2(int_mb(i2_start(1,nb))))


*     **** receive packed array data ****
      reqcnt = 0
!$OMP MASTER
      do it = 1,np-1

*        **** synchronous receive of tmp ****
         proc_from = mod((taskid)-it+np, np)

         msgtype = 19
         source=proc_from
         msglen = (int_mb(i2_start(1,nb)+it+2-1)
     >          -  int_mb(i2_start(1,nb)+it+1-1))


         if (msglen.gt.0) then
            reqcnt = reqcnt + 1
#ifdef MPI4
            stupid_msglen = msglen
            stupid_type   = msgtype
            stupid_taskid = Parallel2d_convert_taskid_i(source)
            call MPI_IRECV(tmp2(int_mb(i2_start(1,nb)+it+1-1)),
     >                    stupid_msglen,stupid_complex,
     >                    stupid_taskid,
     >                    stupid_type,stupid_world,
     >                    stupid_request,stupid_ierr)
            int_mb(request(1)+reqcnt-1) = 0
            int_mb(request(1)+reqcnt-1) = stupid_request
#else
            call MPI_IRECV(tmp2(int_mb(i2_start(1,nb)+it+1-1)),
     >                    msglen,MPI_DOUBLE_COMPLEX,
     >                    Parallel2d_convert_taskid_i(source),
     >                    msgtype,Parallel_comm_world(),
     >                    int_mb(request(1)+reqcnt-1),mpierr)
#endif
         end if

      end do


*     **** send packed array to other processors ****
      do it = 1,np-1

*        **** synchronous send of tmp ****
         proc_to   = mod((taskid)+it, np)
         msglen    = (int_mb(i1_start(1,nb)+it+2-1)
     >              - int_mb(i1_start(1,nb)+it+1-1))
         msgtype   = 19

         if (msglen.gt.0) then
#ifdef MPI4
            stupid_msglen = msglen
            stupid_type   = msgtype
            stupid_taskid = Parallel2d_convert_taskid_i(proc_to)
            call MPI_SEND(tmp1(int_mb(i1_start(1,nb)+it+1-1)),
     >                     stupid_msglen,stupid_complex,
     >                     stupid_taskid,
     >                     stupid_type,stupid_world,stupid_ierr)
#else
            call MPI_SEND(tmp1(int_mb(i1_start(1,nb)+it+1-1)),
     >                     msglen,MPI_DOUBLE_COMPLEX,
     >                     Parallel2d_convert_taskid_i(proc_to),
     >                     msgtype,Parallel_comm_world(),mpierr)
#endif
         end if

      end do


*     **** wait for completion of mp_send, also do a sync ****
      if (np.gt.1) call Parallel_mpiWaitAll0(reqcnt,int_mb(request(1)))
!$OMP END MASTER
!$OMP BARRIER


*     **** unpack A(i) array ****
#ifndef CRAY
!DIR$ ivdep
#endif
!$OMP DO
      do index2=int_mb(i2_start(1,nb)+1-1),
     >         (int_mb(i2_start(1,nb)+np+1-1)-1)
        A(int_mb(iq_to_i2(1,nb)+index2-1))=dconjg(tmp2(index2))
      end do
!$OMP END DO


      !*** deallocate memory ***
      value =           BA_pop_stack(request(2))
      if (.not. value)
     > call errquit(' D3dB_c_timereverse:popping stack',0,MA_ERR)


      return
      end







*     ***********************************
*     *					*
*     *	   D3dB_c_transpose_ijk		*
*     *					*
*     ***********************************

      subroutine D3dB_c_transpose_ijk(nb,op,A,tmp1,tmp2)

*****************************************************
*                                                   *
*      This routine performs the operation          *
*               A(i,k,j) <- A(i,j,k)                * 
*                                                   *
*      np = the number of worker nodes              *
*      proc#=0...(np-1)
*                                                   *
*       this transpose uses more buffer space       *
*       then transpose2                             *
*****************************************************
      implicit none
      integer nb,op
      complex*16  A(*)
      complex*16  tmp1(*),tmp2(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "D3dB.fh"
#include "mpif.h"

#ifdef MPI4
#include "stupid_mpi4.fh"
#endif


*     **** indexing variables ****
      integer h_iq_to_i1(2,6,NBLOCKS)
      integer h_iq_to_i2(2,6,NBLOCKS)
      integer h_i1_start(2,6,NBLOCKS)
      integer h_i2_start(2,6,NBLOCKS)
      common / trans_blk_ijk / h_iq_to_i1,
     >                         h_iq_to_i2,
     >                         h_i1_start,
     >                         h_i2_start

*     **** local variables ***
      logical value
      integer i,nnfft3d
      integer proc_to,proc_from,it
      integer msgtype,source
      integer msglen
      integer request(2),reqcnt
      integer taskid,np
      integer mpierr

*     **** external functions ****
      integer  Parallel2d_convert_taskid_i,Parallel_comm_world
      external Parallel2d_convert_taskid_i,Parallel_comm_world

      integer offset,chunklen
      integer tid,nthr
      integer Parallel_threadid,Parallel_nthreads
      external Parallel_threadid,Parallel_nthreads

     
      call Parallel2d_taskid_i(taskid)
      call Parallel2d_np_i(np)

*     **** allocate memory ****
      value = BA_push_get(mt_int,np,
     >                     'request',request(2),request(1))
      if (.not. value) 
     > call errquit(' D3dB_c_transpose_ijk:out of stack',0,MA_ERR)


*     **** pack A(i) array ****
      if ((op.eq.1).or.(op.eq.5)) nnfft3d = (nx(nb)/2+1)*nq1(nb)
      if ((op.eq.2).or.(op.eq.4)) nnfft3d = (ny(nb))    *nq2(nb)
      if ((op.eq.3).or.(op.eq.6)) nnfft3d = (nz(nb))    *nq3(nb)
#ifndef CRAY
!DIR$ ivdep
#endif

      tid = Parallel_threadid()
c      nthr = Parallel_nthreads()

 
!$OMP DO
      do i=1,nnfft3d
        tmp1(int_mb(h_iq_to_i1(1,op,nb)+i-1)) = A(i)
      end do
!$OMP END DO

*     **** it = 0, transpose data on same thread ****  
      msglen = int_mb(h_i2_start(1,op,nb)+2-1) 
     >       - int_mb(h_i2_start(1,op,nb)+1-1)
c      call dcopy(2*msglen,tmp1(int_mb(h_i1_start(1,op,nb))),1,
c     >                    tmp2(int_mb(h_i2_start(1,op,nb))),1)
      call Parallel_shared_vector_copy(.false.,2*msglen,
     >                    tmp1(int_mb(h_i1_start(1,op,nb))),
     >                    tmp2(int_mb(h_i2_start(1,op,nb))))
         

*     **** receive packed array data ****
      reqcnt = 0
!$OMP MASTER
      do it = 1,np-1

*        **** synchronous receive of tmp ****
         proc_from = mod((taskid)-it+np, np) 

         msgtype = 29
         source=proc_from
         msglen = (int_mb(h_i2_start(1,op,nb)+it+2-1) 
     >          -  int_mb(h_i2_start(1,op,nb)+it+1-1))
         

         if (msglen.gt.0) then
            reqcnt = reqcnt + 1
#ifdef MPI4
            stupid_msglen = msglen
            stupid_type   = msgtype
            stupid_taskid = Parallel2d_convert_taskid_i(source)
            call MPI_IRECV(tmp2(int_mb(h_i2_start(1,op,nb)+it+1-1)),
     >                    stupid_msglen,stupid_complex,
     >                    stupid_taskid,
     >                    stupid_type,stupid_world,
     >                    stupid_request,stupid_ierr) 
            int_mb(request(1)+reqcnt-1) = 0
            int_mb(request(1)+reqcnt-1) = stupid_request
#else
            call MPI_IRECV(tmp2(int_mb(h_i2_start(1,op,nb)+it+1-1)),
     >                    msglen,MPI_DOUBLE_COMPLEX,
     >                    Parallel2d_convert_taskid_i(source),
     >                    msgtype,Parallel_comm_world(),
     >                    int_mb(request(1)+reqcnt-1),mpierr) 
#endif
         end if

      end do

*     **** send packed array to other processors ****
      do it = 1,np-1

*        **** synchronous send of tmp ****
         proc_to   = mod((taskid)+it, np) 
         msglen    = (int_mb(h_i1_start(1,op,nb)+it+2-1)
     >              - int_mb(h_i1_start(1,op,nb)+it+1-1))
         msgtype   = 29
      
         if (msglen.gt.0) then
#ifdef MPI4
            stupid_msglen = msglen
            stupid_type   = msgtype
            stupid_taskid = Parallel2d_convert_taskid_i(proc_to)
            call MPI_SEND(tmp1(int_mb(h_i1_start(1,op,nb)+it+1-1)),
     >                     stupid_msglen,stupid_complex,
     >                     stupid_taskid,
     >                     stupid_type,stupid_world,stupid_ierr)
#else
            call MPI_SEND(tmp1(int_mb(h_i1_start(1,op,nb)+it+1-1)),
     >                     msglen,MPI_DOUBLE_COMPLEX,
     >                     Parallel2d_convert_taskid_i(proc_to),
     >                     msgtype,Parallel_comm_world(),mpierr)
#endif
         end if

      end do

*     **** wait for completion of mp_send, also do a sync ****
      if (np.gt.1) call Parallel_mpiWaitAll0(reqcnt,int_mb(request(1)))

!$OMP END MASTER
!$OMP BARRIER


*     **** unpack A(i) array ****
      if ((op.eq.4).or.(op.eq.6)) nnfft3d = (nx(nb)/2+1)*nq1(nb)
      if ((op.eq.1).or.(op.eq.3)) nnfft3d = (ny(nb))    *nq2(nb)
      if ((op.eq.2).or.(op.eq.5)) nnfft3d = (nz(nb))    *nq3(nb)
#ifndef CRAY
!DIR$ ivdep
#endif

!$OMP DO
      do i=1,nnfft3d
         A(i) = tmp2(int_mb(h_iq_to_i2(1,op,nb)+i-1))
      end do
!$OMP END DO


      !*** deallocate memory ***
      value =           BA_pop_stack(request(2))
      if (.not. value) 
     > call errquit(' D3dB_c_transpose_ijk:popping stack',0,MA_ERR)
      
      return
      end




*     ***********************************
*     *					*
*     *	   D3dB_t_transpose_ijk		*
*     *					*
*     ***********************************

      subroutine D3dB_t_transpose_ijk(nb,op,A,tmp1,tmp2)

*****************************************************
*                                                   *
*      This routine performs the operation          *
*               A(i,k,j) <- A(i,j,k)                * 
*                                                   *
*      np = the number of worker nodes              *
*      proc#=0...(np-1)
*                                                   *
*       this transpose uses more buffer space       *
*       then transpose2                             *
*****************************************************
      implicit none
      integer nb,op
      real*8  A(*)
      real*8  tmp1(*),tmp2(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "mpif.h"
#include "D3dB.fh"

#ifdef MPI4
#include "stupid_mpi4.fh"
#endif


*     **** indexing variables ****
      integer h_iq_to_i1(2,6,NBLOCKS)
      integer h_iq_to_i2(2,6,NBLOCKS)
      integer h_i1_start(2,6,NBLOCKS)
      integer h_i2_start(2,6,NBLOCKS)
      common / trans_blk_ijk / h_iq_to_i1,
     >                         h_iq_to_i2,
     >                         h_i1_start,
     >                         h_i2_start



*     **** local variables ***
      logical value
      integer i,nnfft3d
      integer proc_to,proc_from,it
      integer msgtype,source
c     integer msglen,type,nbytes
c     integer msgid
      integer msglen
c     integer status(MPI_STATUS_SIZE,512),request(512),reqcnt
      integer request(2),reqcnt
      integer taskid,np
      integer mpierr

*     **** external functions ****
      integer  Parallel2d_convert_taskid_i,Parallel_comm_world
      external Parallel2d_convert_taskid_i,Parallel_comm_world
      
      call Parallel2d_taskid_i(taskid)
      call Parallel2d_np_i(np)


*     **** allocate memory ****
      value = BA_push_get(mt_int,np,
     >                     'request',request(2),request(1))
      if (.not. value)
     > call errquit(' D3dB_t_transpose_ijk:out of stack',0,MA_ERR)
      

*     **** pack A(i) array ****
      if ((op.eq.1).or.(op.eq.5)) nnfft3d = (nx(nb)/2+1)*nq1(nb)
      if ((op.eq.2).or.(op.eq.4)) nnfft3d = (ny(nb))    *nq2(nb)
      if ((op.eq.3).or.(op.eq.6)) nnfft3d = (nz(nb))    *nq3(nb)
#ifndef CRAY
!DIR$ ivdep
#endif
      do i=1,nnfft3d
         tmp1(int_mb(h_iq_to_i1(1,op,nb)+i-1)) = A(i)
      end do

*     **** it = 0, transpose data on same thread ****  
      msglen = int_mb(h_i2_start(1,op,nb)+2-1) 
     >       - int_mb(h_i2_start(1,op,nb)+1-1)
      call dcopy(msglen,tmp1(int_mb(h_i1_start(1,op,nb))),1,
     >                  tmp2(int_mb(h_i2_start(1,op,nb))),1)
         

*     **** receive packed array data ****
      reqcnt = 0
      do it = 1,np-1

*        **** synchronous receive of tmp ****
         proc_from = mod((taskid)-it+np, np) 

         msgtype = 39
         source=proc_from
         msglen = (int_mb(h_i2_start(1,op,nb)+it+2-1) 
     >          -  int_mb(h_i2_start(1,op,nb)+it+1-1))
         

         if (msglen.gt.0) then
            reqcnt = reqcnt + 1
#ifdef MPI4
            stupid_msglen = msglen
            stupid_type   = msgtype
            stupid_taskid = Parallel2d_convert_taskid_i(source)
            call MPI_IRECV(tmp2(int_mb(h_i2_start(1,op,nb)+it+1-1)),
     >                    stupid_msglen,stupid_double,
     >                    stupid_taskid,
     >                    stupid_type,stupid_world,
     >                    stupid_request,stupid_ierr) 
            int_mb(request(1)+reqcnt-1) = 0
            int_mb(request(1)+reqcnt-1) = stupid_request
#else
            call MPI_IRECV(tmp2(int_mb(h_i2_start(1,op,nb)+it+1-1)),
     >                    msglen,MPI_DOUBLE_PRECISION,
     >                    Parallel2d_convert_taskid_i(source),
     >                    msgtype,Parallel_comm_world(),
     >                    int_mb(request(1)+reqcnt-1),mpierr) 
#endif
         end if

      end do

*     **** send packed array to other processors ****
      do it = 1,np-1

*        **** synchronous send of tmp ****
         proc_to   = mod((taskid)+it, np) 
         msglen    = (int_mb(h_i1_start(1,op,nb)+it+2-1)
     >              - int_mb(h_i1_start(1,op,nb)+it+1-1))
         msgtype   = 39
      
         if (msglen.gt.0) then
#ifdef MPI4
            stupid_msglen = msglen
            stupid_type   = msgtype
            stupid_taskid = Parallel2d_convert_taskid_i(proc_to)
            call MPI_SEND(tmp1(int_mb(h_i1_start(1,op,nb)+it+1-1)),
     >                     stupid_msglen,stupid_double,
     >                     stupid_taskid,
     >                     stupid_type,stupid_world,stupid_ierr)
#else
            call MPI_SEND(tmp1(int_mb(h_i1_start(1,op,nb)+it+1-1)),
     >                     msglen,MPI_DOUBLE_PRECISION,
     >                     Parallel2d_convert_taskid_i(proc_to),
     >                     msgtype,Parallel_comm_world(),mpierr)
#endif
         end if

      end do

   
*     **** wait for completion of mp_send, also do a sync ****
      if (np.gt.1) call Parallel_mpiWaitAll(reqcnt,int_mb(request(1)))


*     **** unpack A(i) array ****
      if ((op.eq.4).or.(op.eq.6)) nnfft3d = (nx(nb)/2+1)*nq1(nb)
      if ((op.eq.1).or.(op.eq.3)) nnfft3d = (ny(nb))    *nq2(nb)
      if ((op.eq.2).or.(op.eq.5)) nnfft3d = (nz(nb))    *nq3(nb)
#ifndef CRAY
!DIR$ ivdep
#endif
      do i=1,nnfft3d
         A(i) = tmp2(int_mb(h_iq_to_i2(1,op,nb)+i-1))
      end do
      

      !*** deallocate memory ***
      value =           BA_pop_stack(request(2))
      if (.not. value)
     > call errquit(' D3dB_t_transpose_ijk:popping stack',0,MA_ERR)

      return
      end



*     ***********************************
*     *                                 *
*     *         D3dB_SumAll             *
*     *                                 *
*     ***********************************

      subroutine D3dB_SumAll(sum)
c     implicit none
      real*8  sum


#ifdef MPI4
#include "stupid_mpi4.fh"
#else
#include "mpif.h"
#endif

      integer msglen,mpierr,np
      real*8 sumall

*     **** external functions ****
      integer  Parallel2d_comm_i
      external Parallel2d_comm_i

      call Parallel2d_np_i(np)
      if (np.gt.1) then
!$OMP MASTER
         msglen = 1
#ifdef MPI4
         stupid_msglen = msglen
         call MPI_Allreduce(sum,sumall,stupid_msglen,stupid_double,
     >                       stupid_sum,stupid_comm_i,stupid_ierr)
#else
         call MPI_Allreduce(sum,sumall,msglen,MPI_DOUBLE_PRECISION,
     >                       MPI_SUM,Parallel2d_comm_i(),mpierr)
#endif
         sum = sumall
!$OMP END MASTER
!$OMP BARRIER
      end if

      return
      end


*     ***********************************
*     *                                 *
*     *         D3dB_MaxAll             *
*     *                                 *
*     ***********************************
      subroutine D3dB_MaxAll(sum)
c     implicit none
      real*8  sum

#ifdef MPI4
#include "stupid_mpi4.fh"
#else
#include "mpif.h"
#endif

      integer msglen,mpierr,np
      real*8 sumall

*     **** external functions ****
      integer  Parallel2d_comm_i
      external Parallel2d_comm_i

      call Parallel2d_np_i(np)
      if (np.gt.1) then
!$OMP MASTER
         msglen = 1
#ifdef MPI4
         stupid_msglen = msglen
         call MPI_Allreduce(sum,sumall,stupid_msglen,stupid_double,
     >                       stupid_max,stupid_comm_i,stupid_ierr)
#else
         call MPI_Allreduce(sum,sumall,msglen,MPI_DOUBLE_PRECISION,
     >                       MPI_MAX,Parallel2d_comm_i(),mpierr)
#endif
         sum = sumall
!$OMP END MASTER
!$OMP BARRIER
      end if

      return
      end



*     ***********************************
*     *                                 *
*     *         D3dB_ISumAll            *
*     *                                 *
*     ***********************************

      subroutine D3dB_ISumAll(sum)
c     implicit none
      integer sum

#ifdef MPI4
#include "stupid_mpi4.fh"
#else
#include "mpif.h"
#endif

      integer msglen,mpierr,np
      integer sumall
      common /D3dB_ISumAll_tmp/ sumall

*     **** external functions ****
      integer  Parallel2d_comm_i
      external Parallel2d_comm_i

      call Parallel2d_np_i(np)
      if (np.gt.1) then
#ifdef MPI4
!$OMP MASTER
         stupid_msglen = 1
         stupid_tsum = sum
         call MPI_Allreduce(stupid_tsum,stupid_tsum2,stupid_msglen,
     >                      stupid_integer,
     >                      stupid_sum,stupid_comm_i,stupid_ierr)
!$OMP END MASTER
!$OMP BARRIER
         sum = stupid_tsum2
#else
!$OMP MASTER
         msglen = 1
         call MPI_Allreduce(sum,sumall,msglen,MPI_INTEGER,
     >                       MPI_SUM,Parallel2d_comm_i(),mpierr)
!$OMP END MASTER
!$OMP BARRIER
         sum = sumall
#endif
      end if

      return
      end





*     ***********************************
*     *                                 *
*     *         D3dB_Vector_SumAll      *
*     *                                 *
*     ***********************************

      subroutine D3dB_Vector_SumAll(n,sum)
c     implicit none
      integer n
      real*8  sum(*)

#include "bafdecls.fh"
#include "errquit.fh"

#ifdef MPI4
#include "stupid_mpi4.fh"
#else
#include "mpif.h"
#endif

      logical value
      integer msglen,mpierr

*     **** temporary workspace ****
      integer sumall(2),np

*     **** external functions ****
      integer  Parallel2d_comm_i
      external Parallel2d_comm_i


      call Parallel2d_np_i(np)
      call nwpw_timing_start(2)
      if (np.gt.1) then

*     ***** allocate temporary space ****
      value = BA_push_get(mt_dbl,n,'sumall',sumall(2),sumall(1))
      if (.not. value) 
     > call errquit('D3dB_Vector_SumAll:out of stack',0, MA_ERR)


!$OMP MASTER
#ifdef MPI4
      stupid_msglen = n
      call MPI_Allreduce(sum,dbl_mb(sumall(1)),
     >                stupid_msglen,
     >                stupid_double,
     >                stupid_sum,stupid_comm_i,stupid_ierr)
#else
      call MPI_Allreduce(sum,dbl_mb(sumall(1)),n,
     >                MPI_DOUBLE_PRECISION,
     >                MPI_SUM,Parallel2d_comm_i(),mpierr)
#endif

      call dcopy(n,dbl_mb(sumall(1)),1,sum,1)
!$OMP END MASTER
!$OMP BARRIER

      value = BA_pop_stack(sumall(2))

      end if
      call nwpw_timing_end(2)
      return
      end




*     ***********************************
*     *                                 *
*     *         D3dB_Vector_ISumAll     *
*     *                                 *
*     ***********************************

      subroutine D3dB_Vector_ISumAll(n,sum)
c     implicit none
      integer n
      integer  sum(*)

#include "bafdecls.fh"
#include "errquit.fh"

#ifdef MPI4
#include "stupid_mpi4.fh"
#else
#include "mpif.h"
#endif

      logical value
      integer msglen,mpierr

*     **** temporary workspace ****
      integer sumall(2),sumall1(2),i,taskid

*     **** external functions ****
      integer  Parallel2d_comm_i
      external Parallel2d_comm_i

      call Parallel_taskid(taskid)
      call nwpw_timing_start(2)

      call Parallel2d_np_i(np)
      if (np.gt.1) then

*     ***** allocate temporary space ****
      value = BA_push_get(mt_int,n,'sumall',sumall(2),sumall(1))
      if (.not. value) 
     > call errquit('D3dB_Vector_ISumall:out of stack',0,MA_ERR)

#ifdef MPI4
      if (.not.BA_push_get(mt_int,n,'sumall1',sumall1(2),sumall1(1)))
     > call errquit('D3dB_Vector_ISumAll:out of stack',1, MA_ERR)
      stupid_msglen = n
      call stupid_icopy8to4(n,sum,int_mb(sumall1(1)))
      call MPI_Allreduce(int_mb(sumall1(1)),
     >                   int_mb(sumall(1)),stupid_msglen,
     >                stupid_integer,
     >                stupid_sum,stupid_comm_i,stupid_ierr)
      call stupid_icopy4to8(n,int_mb(sumall(1)),sum)
      if (.not.BA_pop_stack(sumall1(2)))
     >  call errquit('D3dB_Vector_ISumAll:popping stack',0, MA_ERR)
#else
      msglen = n
      call MPI_Allreduce(sum,int_mb(sumall(1)),msglen,
     >                MPI_INTEGER,
     >                MPI_SUM,Parallel2d_comm_i(),mpierr)
      call icopy(n,int_mb(sumall(1)),1,sum,1)
#endif

      value = BA_pop_stack(sumall(2))
      if (.not. value) 
     > call errquit('D3dB_Vector_ISumAll:popping stack',0, MA_ERR)

      end if

      call nwpw_timing_end(2)
      return
      end




      
*     ***********************************
*     *					*
*     *	   D3dB_c_ptranspose1_jk	*
*     *					*
*     ***********************************

      subroutine D3dB_c_ptranspose1_jk(nbb,A,tmp1,tmp2)

*****************************************************
*                                                   *
*      This routine performs the operation          *
*               A(i,k,j) <- A(i,j,k)                * 
*                                                   *
*      np = the number of worker nodes              *
*      proc#=0...(np-1)
*                                                   *
*       this transpose uses more buffer space       *
*       then transpose2                             *
*****************************************************
      implicit none
      integer nbb
      complex*16  A(*)
      complex*16  tmp1(*),tmp2(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "mpif.h"
#include "D3dB.fh"

#ifdef MPI4
#include "stupid_mpi4.fh"
#endif

*     **** indexing variables ****
      integer iq_to_i1(2,0:3)
      integer iq_to_i2(2,0:3)
      integer iz_to_i2(2,0:3)
      integer i1_start(2,0:3)
      integer i2_start(2,0:3)
      common / ptrans_blk1 / iq_to_i1,iq_to_i2,iz_to_i2,
     >                       i1_start,i2_start


*     **** local variables ***
      logical value
      integer proc_to,proc_from,it
      integer msgtype,source
c     integer msglen,type,nbytes
c     integer msgid
      integer msglen
c     integer status(MPI_STATUS_SIZE,512),request(512),reqcnt
      integer request(2),reqcnt
      integer taskid,np
      integer mpierr,n1,n2,id

*     **** external functions ****
      integer  Parallel2d_convert_taskid_i,Parallel_comm_world
      external Parallel2d_convert_taskid_i,Parallel_comm_world
      
      !id = 2*(nbb/2) + 1
      if (nbb.lt.2) then
         id = 1
      else
         id = 3
      end if
      call Parallel2d_taskid_i(taskid)
      call Parallel2d_np_i(np)

*     **** allocate memory ****
      value = BA_push_get(mt_int,np,
     >                     'request',request(2),request(1))
      if (.not. value)
     > call errquit(' D3dB_c_ptranspose1_jk:out of stack',0,MA_ERR)


      n1 = int_mb(i1_start(1,nbb)+np) - 1
      n2 = int_mb(i2_start(1,nbb)+np) - 1

*     **** pack A(i) array ****
       call D3dB_pfft_index1_copy(n1,int_mb(iq_to_i1(1,nbb)),A,tmp1)

*     **** it = 0, transpose data on same thread ****  
      msglen = int_mb(i2_start(1,nbb)+2-1) - int_mb(i2_start(1,nbb)+1-1)
      call dcopy(2*msglen,
     >           tmp1(int_mb(i1_start(1,nbb)+1-1)),1,
     >           tmp2(int_mb(i2_start(1,nbb)+1-1)),1)
         

*     **** receive packed array data ****
      reqcnt = 0
      do it = 1,np-1

*        **** synchronous receive of tmp ****
         proc_from = mod((taskid)-it+np, np) 

         msgtype = 49
         source=proc_from
         msglen = (int_mb(i2_start(1,nbb)+it+2-1) 
     >          -  int_mb(i2_start(1,nbb)+it+1-1))
         

         if (msglen.gt.0) then
            reqcnt = reqcnt + 1
#ifdef MPI4
            stupid_msglen = msglen
            stupid_type   = msgtype
            stupid_taskid = Parallel2d_convert_taskid_i(source)
            call MPI_IRECV(tmp2(int_mb(i2_start(1,nbb)+it+1-1)),
     >                    stupid_msglen,stupid_complex,
     >                    stupid_taskid,
     >                    stupid_type,stupid_world,
     >                    stupid_request,stupid_ierr) 
            int_mb(request(1)+reqcnt-1) = 0
            int_mb(request(1)+reqcnt-1) = stupid_request
#else
            call MPI_IRECV(tmp2(int_mb(i2_start(1,nbb)+it+1-1)),
     >                    msglen,MPI_DOUBLE_COMPLEX,
     >                    Parallel2d_convert_taskid_i(source),
     >                    msgtype,Parallel_comm_world(),
     >                    int_mb(request(1)+reqcnt-1),mpierr) 
#endif
         end if

      end do

*     **** send packed array to other processors ****
      do it = 1,np-1

*        **** synchronous send of tmp ****
         proc_to   = mod((taskid)+it, np) 
         msglen    = (int_mb(i1_start(1,nbb)+it+2-1)
     >              - int_mb(i1_start(1,nbb)+it+1-1))
         msgtype   = 49
      
         if (msglen.gt.0) then
#ifdef MPI4
            stupid_msglen = msglen
            stupid_type   = msgtype
            stupid_taskid = Parallel2d_convert_taskid_i(proc_to)
            call MPI_SEND(tmp1(int_mb(i1_start(1,nbb)+it+1-1)),
     >                     stupid_msglen,stupid_complex,
     >                     stupid_taskid,
     >                     stupid_type,stupid_world,stupid_ierr)
#else
            call MPI_SEND(tmp1(int_mb(i1_start(1,nbb)+it+1-1)),
     >                     msglen,MPI_DOUBLE_COMPLEX,
     >                     Parallel2d_convert_taskid_i(proc_to),
     >                     msgtype,Parallel_comm_world(),mpierr)
#endif
         end if

      end do
   
*     **** wait for completion of mp_send, also do a sync ****
      if (np.gt.1) call Parallel_mpiWaitAll(reqcnt,int_mb(request(1)))


*     **** unpack A(i) array ****
      call D3dB_pfft_index2_copy(n2,int_mb(iq_to_i2(1,nbb)),tmp2,A)
      call D3dB_pfft_index2_zero(nfft3d(id)-n2,
     >                           int_mb(iz_to_i2(1,nbb)),A)
      

      !*** deallocate memory ***
      value =           BA_pop_stack(request(2))
      if (.not. value)
     > call errquit(' D3dB_c_ptranspose1_jk:popping stack',0,MA_ERR)

      return
      end




      
*     ***********************************
*     *					*
*     *	   D3dB_c_ptranspose2_jk	*
*     *					*
*     ***********************************

      subroutine D3dB_c_ptranspose2_jk(nbb,A,tmp1,tmp2)

*****************************************************
*                                                   *
*      This routine performs the operation          *
*               A(i,k,j) <- A(i,j,k)                * 
*                                                   *
*      np = the number of worker nodes              *
*      proc#=0...(np-1)
*                                                   *
*       this transpose uses more buffer space       *
*       then transpose2                             *
*****************************************************
      implicit none
      integer nbb
      complex*16  A(*)
      complex*16  tmp1(*),tmp2(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "mpif.h"
#include "D3dB.fh"

#ifdef MPI4
#include "stupid_mpi4.fh"
#endif


*     **** indexing variables ****
      integer iq_to_i1(2,0:3)
      integer iq_to_i2(2,0:3)
      integer iz_to_i2(2,0:3)
      integer i1_start(2,0:3)
      integer i2_start(2,0:3)
      common / ptrans_blk2 / iq_to_i1,iq_to_i2,iz_to_i2,
     >                       i1_start,i2_start


*     **** local variables ***
      logical value
      integer proc_to,proc_from,it
      integer msgtype,source
      integer msglen
      integer request(2),reqcnt
      integer taskid,np
      integer mpierr,n1,n2,id

c     **** external functions ****
      integer  Parallel2d_convert_taskid_i,Parallel_comm_world
      external Parallel2d_convert_taskid_i,Parallel_comm_world
      
      !id = 2*(nbb/2)+1
      if (nbb.lt.2) then
         id = 1
      else
         id = 3
      end if
      call Parallel2d_taskid_i(taskid)
      call Parallel2d_np_i(np)

*     **** allocate memory ****
      value = BA_push_get(mt_int,np,
     >                     'request',request(2),request(1))
      if (.not. value)
     > call errquit(' D3dB_c_ptranspose2_jk:out of stack',0,MA_ERR)

      n1 = int_mb(i1_start(1,nbb)+np) - 1
      n2 = int_mb(i2_start(1,nbb)+np) - 1

*     **** pack A(i) array ****
       call D3dB_pfft_index1_copy(n1,int_mb(iq_to_i1(1,nbb)),A,tmp1)

*     **** it = 0, transpose data on same thread ****  
      msglen = int_mb(i2_start(1,nbb)+2-1) - int_mb(i2_start(1,nbb)+1-1)
      call dcopy(2*msglen,
     >           tmp1(int_mb(i1_start(1,nbb)+1-1)),1,
     >           tmp2(int_mb(i2_start(1,nbb)+1-1)),1)
         

*     **** receive packed array data ****
      reqcnt = 0
      do it = 1,np-1

*        **** synchronous receive of tmp ****
         proc_from = mod((taskid)-it+np, np) 

         msgtype = 59
         source=proc_from
         msglen = (int_mb(i2_start(1,nbb)+it+2-1) 
     >          -  int_mb(i2_start(1,nbb)+it+1-1))
         

         if (msglen.gt.0) then
            reqcnt = reqcnt + 1
#ifdef MPI4
            stupid_msglen = msglen
            stupid_type   = msgtype
            stupid_taskid = Parallel2d_convert_taskid_i(source)
            call MPI_IRECV(tmp2(int_mb(i2_start(1,nbb)+it+1-1)),
     >                    stupid_msglen,stupid_complex,
     >                    stupid_taskid,
     >                    stupid_type,stupid_world,
     >                    stupid_request,stupid_ierr) 
            int_mb(request(1)+reqcnt-1) = 0
            int_mb(request(1)+reqcnt-1) = stupid_request
#else
            call MPI_IRECV(tmp2(int_mb(i2_start(1,nbb)+it+1-1)),
     >                    msglen,MPI_DOUBLE_COMPLEX,
     >                    Parallel2d_convert_taskid_i(source),
     >                    msgtype,Parallel_comm_world(),
     >                    int_mb(request(1)+reqcnt-1),mpierr) 
#endif
         end if

      end do

*     **** send packed array to other processors ****
      do it = 1,np-1

*        **** synchronous send of tmp ****
         proc_to   = mod((taskid)+it, np) 
         msglen    = (int_mb(i1_start(1,nbb)+it+2-1)
     >              - int_mb(i1_start(1,nbb)+it+1-1))
         msgtype   = 59
      
         if (msglen.gt.0) then
#ifdef MPI4
            stupid_msglen = msglen
            stupid_type   = msgtype
            stupid_taskid = Parallel2d_convert_taskid_i(proc_to)
            call MPI_SEND(tmp1(int_mb(i1_start(1,nbb)+it+1-1)),
     >                     stupid_msglen,stupid_complex,
     >                     stupid_taskid,
     >                     stupid_type,stupid_world,stupid_ierr)
#else
            call MPI_SEND(tmp1(int_mb(i1_start(1,nbb)+it+1-1)),
     >                     msglen,MPI_DOUBLE_COMPLEX,
     >                     Parallel2d_convert_taskid_i(proc_to),
     >                     msgtype,Parallel_comm_world(),mpierr)
#endif
         end if

      end do
   
*     **** wait for completion of mp_send, also do a sync ****
      if (np.gt.1) call Parallel_mpiWaitAll(reqcnt,int_mb(request(1)))


*     **** unpack A(i) array ****
      call D3dB_pfft_index2_copy(n2,int_mb(iq_to_i2(1,nbb)),tmp2,A)
      call D3dB_pfft_index2_zero(nfft3d(id)-n2,
     >                           int_mb(iz_to_i2(1,nbb)),A)
      

      !*** deallocate memory ***
      value =           BA_pop_stack(request(2))
      if (.not. value)
     > call errquit(' D3dB_c_ptranspose2_jk:popping stack',0,MA_ERR)

      return
      end





*     ***********************************
*     *					*
*     *	   D3dB_c_ptranspose_ijk	*
*     *					*
*     ***********************************

      subroutine D3dB_c_ptranspose_ijk(nbb,op,A,tmp1,tmp2)

*****************************************************
*                                                   *
*      This routine performs the operation          *
*               A(i,k,j) <- A(i,j,k)                * 
*                                                   *
*      np = the number of worker nodes              *
*      proc#=0...(np-1)
*                                                   *
*       this transpose uses more buffer space       *
*       then transpose2                             *
*****************************************************
      implicit none
      integer nbb,op
      complex*16  A(*)
      complex*16  tmp1(*),tmp2(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "mpif.h"
#include "D3dB.fh"

#ifdef MPI4
#include "stupid_mpi4.fh"
#endif


*     **** indexing variables ****
      integer h_iq_to_i1(2,6,0:3)
      integer h_iq_to_i2(2,6,0:3)
      integer h_iz_to_i2(2,6,0:3)
      integer h_iz_to_i2_count(6,0:3)
      integer h_i1_start(2,6,0:3)
      integer h_i2_start(2,6,0:3)
      common / ptrans_blk_ijk / h_iq_to_i1,
     >                         h_iq_to_i2,
     >                         h_iz_to_i2,
     >                         h_iz_to_i2_count,
     >                         h_i1_start,
     >                         h_i2_start



*     **** local variables ***
      logical value
      integer n1,n2,n3
      integer proc_to,proc_from,it
      integer msgtype,source
      integer msglen
      integer request(2),reqcnt
      integer taskid,np
      integer mpierr

c     **** external functions ****
      integer  Parallel2d_convert_taskid_i,Parallel_comm_world
      external Parallel2d_convert_taskid_i,Parallel_comm_world
      
      call Parallel2d_taskid_i(taskid)
      call Parallel2d_np_i(np)

*     **** allocate memory ****
      value = BA_push_get(mt_int,np,
     >                     'request',request(2),request(1))
      if (.not. value) 
     > call errquit(' D3dB_c_ptranspose_ijk:out of stack',0,MA_ERR)


      n1 = int_mb(h_i1_start(1,op,nbb)+np) - 1
      n2 = int_mb(h_i2_start(1,op,nbb)+np) - 1
      n3 = h_iz_to_i2_count(op,nbb)


*     **** pack A(i) array ****
      call D3dB_pfft_index1_copy(n1,int_mb(h_iq_to_i1(1,op,nbb)),A,tmp1)


*     **** it = 0, transpose data on same thread ****  
      msglen = int_mb(h_i2_start(1,op,nbb)+2-1) 
     >       - int_mb(h_i2_start(1,op,nbb)+1-1)
c      call dcopy(2*msglen,
c     >           tmp1(int_mb(h_i1_start(1,op,nbb))),1,
c     >           tmp2(int_mb(h_i2_start(1,op,nbb))),1)
      call Parallel_shared_vector_copy(.true.,2*msglen,
     >           tmp1(int_mb(h_i1_start(1,op,nbb))),
     >           tmp2(int_mb(h_i2_start(1,op,nbb))))


!$OMP MASTER
*     **** receive packed array data ****
      reqcnt = 0
      do it = 1,np-1

*        **** synchronous receive of tmp ****
         proc_from = mod((taskid)-it+np, np) 

         msgtype = 69
         source=proc_from
         msglen = (int_mb(h_i2_start(1,op,nbb)+it+2-1) 
     >          -  int_mb(h_i2_start(1,op,nbb)+it+1-1))
         

         if (msglen.gt.0) then
            reqcnt = reqcnt + 1
#ifdef MPI4
            stupid_msglen = msglen
            stupid_type   = msgtype
            stupid_taskid = Parallel2d_convert_taskid_i(source)
            call MPI_IRECV(tmp2(int_mb(h_i2_start(1,op,nbb)+it+1-1)),
     >                    stupid_msglen,stupid_complex,
     >                    stupid_taskid,
     >                    stupid_type,stupid_world,
     >                    stupid_request,stupid_ierr) 
            int_mb(request(1)+reqcnt-1) = 0
            int_mb(request(1)+reqcnt-1) = stupid_request
#else
            call MPI_IRECV(tmp2(int_mb(h_i2_start(1,op,nbb)+it+1-1)),
     >                    msglen,MPI_DOUBLE_COMPLEX,
     >                    Parallel2d_convert_taskid_i(source),
     >                    msgtype,Parallel_comm_world(),
     >                    int_mb(request(1)+reqcnt-1),mpierr) 
#endif
         end if

      end do

*     **** send packed array to other processors ****
      do it = 1,np-1

*        **** synchronous send of tmp ****
         proc_to   = mod((taskid)+it, np) 
         msglen    = (int_mb(h_i1_start(1,op,nbb)+it+2-1)
     >              - int_mb(h_i1_start(1,op,nbb)+it+1-1))
         msgtype   = 69
      
         if (msglen.gt.0) then
#ifdef MPI4
            stupid_msglen = msglen
            stupid_type   = msgtype
            stupid_taskid = Parallel2d_convert_taskid_i(proc_to)
            call MPI_SEND(tmp1(int_mb(h_i1_start(1,op,nbb)+it+1-1)),
     >                     stupid_msglen,stupid_complex,
     >                     stupid_taskid,
     >                     stupid_type,stupid_world,stupid_ierr)
#else
            call MPI_SEND(tmp1(int_mb(h_i1_start(1,op,nbb)+it+1-1)),
     >                     msglen,MPI_DOUBLE_COMPLEX,
     >                     Parallel2d_convert_taskid_i(proc_to),
     >                     msgtype,Parallel_comm_world(),mpierr)
#endif
         end if

      end do

   
*     **** wait for completion of mp_send, also do a sync ****
      if (np.gt.1) call Parallel_mpiWaitAll0(reqcnt,int_mb(request(1)))

!$OMP END MASTER
!$OMP BARRIER

*     **** unpack A(i) array ****
      call D3dB_pfft_index2_copy(n2,int_mb(h_iq_to_i2(1,op,nbb)),tmp2,A)
      call D3dB_pfft_index2_zero(n3,int_mb(h_iz_to_i2(1,op,nbb)),A)


      !*** deallocate memory ***
      value =           BA_pop_stack(request(2))
      if (.not. value) 
     > call errquit(' D3dB_c_transpose_ijk:popping stack',0,MA_ERR)
      
      return
      end




*     ************************************
*     *                                  *
*     *         Balance_c_balance_start  *
*     *                                  *
*     ************************************

      subroutine Balance_c_balance_start(nb,A,request,reqcnt,msgtype)
      implicit none
      integer nb
      complex*16 A(*)
      integer    request(*),reqcnt
      integer msgtype

#include "bafdecls.fh"
#include "mpif.h"
#include "balance_common.fh"

#ifdef MPI4
#include "stupid_mpi4.fh"
#endif


*     **** local variables ****
      integer j,np
      integer pto,pfrom,msglen,indx,mpierr

c     **** external functions ****
      integer  Parallel2d_convert_taskid_i,Parallel_comm_world
      external Parallel2d_convert_taskid_i,Parallel_comm_world


!$OMP MASTER
      reqcnt = 0
      if (log_mb(sender_list(1)+nb)) then
         do j=1,int_mb(npacket_list(1)+nb)
            pto    = int_mb(int_mb(proc_to_list(1)    +2*nb)+j-1)    
            msglen = int_mb(int_mb(packet_size_list(1)+2*nb)+j-1)
            indx   = int_mb(int_mb(indx_start_list(1) +2*nb)+j-1)
c            msgtype=2
c            send data....
            if (msglen.gt.0) then
            reqcnt = reqcnt + 1
#ifdef MPI4
            stupid_msglen = msglen
            stupid_type   = msgtype
            stupid_taskid = Parallel2d_convert_taskid_i(pto)
            call MPI_ISEND(A(indx),
     >                     stupid_msglen,stupid_complex,
     >                     stupid_taskid,
     >                     stupid_type,stupid_world,
     >                     stupid_request,stupid_ierr)
            request(reqcnt) = 0
            request(reqcnt) = stupid_request
#else
            call MPI_ISEND(A(indx),
     >                     msglen,MPI_DOUBLE_COMPLEX,
     >                     Parallel2d_convert_taskid_i(pto),
     >                     msgtype,Parallel_comm_world(),
     >                     request(reqcnt),mpierr)
#endif
            end if
  

         end do
      end if

      if (log_mb(receiver_list(1)+nb)) then
         do j=1,int_mb(npacket_list(1)+nb)
            pfrom  = int_mb(int_mb(proc_from_list(1)  +2*nb)+j-1)    
            msglen = int_mb(int_mb(packet_size_list(1)+2*nb)+j-1)
            indx   = int_mb(int_mb(indx_start_list(1) +2*nb)+j-1)
c            msgtype=2
c            recieve data....
            if (msglen.gt.0) then
            reqcnt = reqcnt + 1
#ifdef MPI4
            stupid_msglen = msglen
            stupid_type   = msgtype
            stupid_taskid = Parallel2d_convert_taskid_i(pfrom)
            call MPI_IRECV(A(indx),
     >                    stupid_msglen,stupid_complex,
     >                    stupid_taskid,
     >                    stupid_type,stupid_world,
     >                    stupid_request,stupid_ierr)
            request(reqcnt) = 0
            request(reqcnt) = stupid_request
#else
            call MPI_IRECV(A(indx),
     >                    msglen,MPI_DOUBLE_COMPLEX,
     >                    Parallel2d_convert_taskid_i(pfrom),
     >                    msgtype,Parallel_comm_world(),
     >                    request(reqcnt),mpierr)
#endif
            end if

         end do
      end if
!$OMP END MASTER

      return
      end




*     ************************************
*     *                                  *
*     *         Balance_c_balance_end    *
*     *                                  *
*     ************************************

      subroutine Balance_c_balance_end(nb,A,request,reqcnt)
      implicit none
      integer nb
      complex*16 A(*)
      integer    request(*),reqcnt

*     **** local variables ***
      logical value
      integer np

      call Parallel2d_np_i(np)
      if (np.gt.1) call Parallel_mpiWaitAll(reqcnt,request)
      return
      end

*     ************************************
*     *                                  *
*     *    Balance_c_unbalance_start     *
*     *                                  *
*     ************************************

      subroutine Balance_c_unbalance_start(nb,A,request,reqcnt,msgtype)
      implicit none
      integer nb
      complex*16 A(*)
      integer    request(*),reqcnt
      integer    msgtype

#include "bafdecls.fh"
#include "mpif.h"
#include "balance_common.fh"

#ifdef MPI4
#include "stupid_mpi4.fh"
#endif

*     **** local variables ****
      integer j,np
      integer pto,pfrom,msglen,indx,mpierr

c     **** external functions ****
      integer  Parallel2d_convert_taskid_i,Parallel_comm_world
      external Parallel2d_convert_taskid_i,Parallel_comm_world


!$OMP MASTER
      reqcnt=0
      if (log_mb(sender_list(1)+nb)) then
         do j=1,int_mb(npacket_list(1)+nb)
            pfrom  = int_mb(int_mb(proc_to_list(1)    +2*nb)+j-1)
            msglen = int_mb(int_mb(packet_size_list(1)+2*nb)+j-1)
            indx   = int_mb(int_mb(indx_start_list(1) +2*nb)+j-1)
c            msgtype=2
c            recieve data....
            if (msglen.gt.0) then
            reqcnt = reqcnt + 1
#ifdef MPI4
            stupid_msglen = msglen
            stupid_type   = msgtype
            stupid_taskid = Parallel2d_convert_taskid_i(pfrom)
            call MPI_IRECV(A(indx),
     >                    stupid_msglen,stupid_complex,
     >                    stupid_taskid,
     >                    stupid_type,stupid_world,
     >                    stupid_request,stupid_ierr)
            request(reqcnt) = 0
            request(reqcnt) = stupid_request
#else
            call MPI_IRECV(A(indx),
     >                    msglen,MPI_DOUBLE_COMPLEX,
     >                    Parallel2d_convert_taskid_i(pfrom),
     >                    msgtype,Parallel_comm_world(),
     >                    request(reqcnt),mpierr)
#endif
            end if
  
         end do
      end if

      if (log_mb(receiver_list(1)+nb)) then
         do j=1,int_mb(npacket_list(1)+nb)
            pto    = int_mb(int_mb(proc_from_list(1)  +2*nb)+j-1)    
            msglen = int_mb(int_mb(packet_size_list(1)+2*nb)+j-1)
            indx   = int_mb(int_mb(indx_start_list(1) +2*nb)+j-1)
c            msgtype=2
c            send data....
            if (msglen.gt.0) then
            reqcnt = reqcnt + 1
#ifdef MPI4
            stupid_msglen = msglen
            stupid_type   = msgtype
            stupid_taskid = Parallel2d_convert_taskid_i(pto)
            call MPI_ISEND(A(indx),
     >                     stupid_msglen,stupid_complex,
     >                     stupid_taskid,
     >                     stupid_type,stupid_world,
     >                     stupid_request,stupid_ierr)
            request(reqcnt) = 0
            request(reqcnt) = stupid_request
#else
            call MPI_ISEND(A(indx),
     >                     msglen,MPI_DOUBLE_COMPLEX,
     >                     Parallel2d_convert_taskid_i(pto),
     >                     msgtype,Parallel_comm_world(),
     >                     request(reqcnt),mpierr)
#endif

            end if

         end do
      end if
!$OMP END MASTER


      return
      end


*     ************************************
*     *                                  *
*     *       Balance_c_unbalance_end    *
*     *                                  *
*     ************************************

      subroutine Balance_c_unbalance_end(nb,A,request,reqcnt)
      implicit none
      integer nb
      complex*16 A(*)
      integer    request(*),reqcnt

*     **** local variables ***
      integer np

      call Parallel2d_np_i(np)
      if (np.gt.1) call Parallel_mpiWaitAll(reqcnt,request)
      return
      end




      
*     ***********************************
*     *					*
*     *	   D3dB_c_ptranspose1_jk_start	*
*     *					*
*     ***********************************

*                                           
*      This routine performs the operation   
*               A(i,k,j) <- A(i,j,k)          
*                                              
*      np = the number of worker nodes          
*      proc#=0...(np-1)
*                                                
*       this transpose uses more buffer space     
*       then transpose2                            
*

      subroutine D3dB_c_ptranspose1_jk_start(nbb,A,tmp1,tmp2,
     >                                       request,reqcnt,msgtype)
      implicit none
      integer nbb
      complex*16  A(*)
      complex*16  tmp1(*),tmp2(*)
      integer request(*),reqcnt
      integer msgtype

#include "bafdecls.fh"
#include "errquit.fh"
#include "mpif.h"
#include "D3dB.fh"

#ifdef MPI4
#include "stupid_mpi4.fh"
#endif


*     **** indexing variables ****
      integer iq_to_i1(2,0:3)
      integer iq_to_i2(2,0:3)
      integer iz_to_i2(2,0:3)
      integer i1_start(2,0:3)
      integer i2_start(2,0:3)
      common / ptrans_blk1 / iq_to_i1,iq_to_i2,iz_to_i2,
     >                       i1_start,i2_start


*     **** local variables ***
      integer proc_to,proc_from,it
      integer source
      integer msglen
      integer taskid,np
      integer mpierr,n1

c     **** external functions ****
      integer  Parallel2d_convert_taskid_i,Parallel_comm_world
      external Parallel2d_convert_taskid_i,Parallel_comm_world
      
      call Parallel2d_taskid_i(taskid)
      call Parallel2d_np_i(np)

      n1 = int_mb(i1_start(1,nbb)+np) - 1

*     **** pack A(i) array ****
       call D3dB_pfft_index1_copy(n1,int_mb(iq_to_i1(1,nbb)),A,tmp1)

*     **** it = 0, transpose data on same thread ****  
      msglen = int_mb(i2_start(1,nbb)+2-1) - int_mb(i2_start(1,nbb)+1-1)
      call dcopy(2*msglen,
     >           tmp1(int_mb(i1_start(1,nbb)+1-1)),1,
     >           tmp2(int_mb(i2_start(1,nbb)+1-1)),1)
         

*     **** receive packed array data ****
      reqcnt = 0
      do it = 1,np-1

*        **** synchronous receive of tmp ****
         proc_from = mod((taskid)-it+np, np) 

c         msgtype = 1
c         type   = 1
         source=proc_from
         msglen = (int_mb(i2_start(1,nbb)+it+2-1) 
     >          -  int_mb(i2_start(1,nbb)+it+1-1))
         

         if (msglen.gt.0) then
            reqcnt = reqcnt + 1
#ifdef MPI4
            stupid_msglen = msglen
            stupid_type   = msgtype
            stupid_taskid = Parallel2d_convert_taskid_i(source)
            call MPI_IRECV(tmp2(int_mb(i2_start(1,nbb)+it+1-1)),
     >                    stupid_msglen,stupid_complex,
     >                    stupid_taskid,
     >                    stupid_type,stupid_world,
     >                    stupid_request,stupid_ierr) 
            request(reqcnt) = 0
            request(reqcnt) = stupid_request
#else
            call MPI_IRECV(tmp2(int_mb(i2_start(1,nbb)+it+1-1)),
     >                    msglen,MPI_DOUBLE_COMPLEX,
     >                    Parallel2d_convert_taskid_i(source),
     >                    msgtype,Parallel_comm_world(),
     >                    request(reqcnt),mpierr) 
#endif
         end if

      end do

*     **** send packed array to other processors ****
      do it = 1,np-1

*        **** synchronous send of tmp ****
         proc_to   = mod((taskid)+it, np) 
         msglen    = (int_mb(i1_start(1,nbb)+it+2-1)
     >              - int_mb(i1_start(1,nbb)+it+1-1))
c         msgtype   = 1
      
         if (msglen.gt.0) then
            reqcnt = reqcnt + 1
#ifdef MPI4
            stupid_msglen = msglen
            stupid_type   = msgtype
            stupid_taskid = Parallel2d_convert_taskid_i(proc_to)
            call MPI_ISEND(tmp1(int_mb(i1_start(1,nbb)+it+1-1)),
     >                     stupid_msglen,stupid_complex,
     >                     stupid_taskid,
     >                     stupid_type,stupid_world,
     >                     stupid_request,stupid_ierr)
            request(reqcnt) = 0
            request(reqcnt) = stupid_request
#else
            call MPI_ISEND(tmp1(int_mb(i1_start(1,nbb)+it+1-1)),
     >                     msglen,MPI_DOUBLE_COMPLEX,
     >                     Parallel2d_convert_taskid_i(proc_to),
     >                     msgtype,Parallel_comm_world(),
     >                     request(reqcnt),mpierr)
#endif
         end if

      end do
   

      return
      end


*     ***********************************
*     *					*
*     *	   D3dB_c_ptranspose1_jk_end	*
*     *					*
*     ***********************************

*                                           
*      This routine performs the operation   
*               A(i,k,j) <- A(i,j,k)          
*                                              
*      np = the number of worker nodes          
*      proc#=0...(np-1)
*                                                
*       this transpose uses more buffer space     
*       then transpose2                            
*

      subroutine D3dB_c_ptranspose1_jk_end(nbb,A,tmp2,request,reqcnt)
      implicit none
      integer nbb
      complex*16  A(*)
      complex*16  tmp2(*)
      integer     request(*),reqcnt

#include "bafdecls.fh"
#include "D3dB.fh"


*     **** indexing variables ****
      integer iq_to_i1(2,0:3)
      integer iq_to_i2(2,0:3)
      integer iz_to_i2(2,0:3)
      integer i1_start(2,0:3)
      integer i2_start(2,0:3)
      common / ptrans_blk1 / iq_to_i1,iq_to_i2,iz_to_i2,
     >                       i1_start,i2_start

*     **** local variables ***
      integer np,n2,id

      !id=2*(nbb/2)+1
      if (nbb.lt.2) then
         id = 1
      else
         id = 3
      end if
      
      call Parallel2d_np_i(np)
      if (np.gt.1) call Parallel_mpiWaitAll(reqcnt,request)

*     **** unpack A(i) array ****
      n2 = int_mb(i2_start(1,nbb)+np) - 1
      call D3dB_pfft_index2_copy(n2,int_mb(iq_to_i2(1,nbb)),tmp2,A)
      call D3dB_pfft_index2_zero(nfft3d(id)-n2,
     >                           int_mb(iz_to_i2(1,nbb)),A)
      return
      end





      
*     ***********************************
*     *					*
*     *	   D3dB_c_ptranspose2_jk_start	*
*     *					*
*     ***********************************

*                                           
*      This routine performs the operation   
*               A(i,k,j) <- A(i,j,k)          
*                                              
*      np = the number of worker nodes          
*      proc#=0...(np-1)
*                                                
*       this transpose uses more buffer space     
*       then transpose2                            
*

      subroutine D3dB_c_ptranspose2_jk_start(nbb,A,tmp1,tmp2,
     >                                       request,reqcnt,msgtype)
      implicit none
      integer nbb
      complex*16  A(*)
      complex*16  tmp1(*),tmp2(*)
      integer request(*),reqcnt
      integer msgtype

#include "bafdecls.fh"
#include "errquit.fh"
#include "mpif.h"
#include "D3dB.fh"

#ifdef MPI4
#include "stupid_mpi4.fh"
#endif

*     **** indexing variables ****
      integer iq_to_i1(2,0:3)
      integer iq_to_i2(2,0:3)
      integer iz_to_i2(2,0:3)
      integer i1_start(2,0:3)
      integer i2_start(2,0:3)
      common / ptrans_blk2 / iq_to_i1,iq_to_i2,iz_to_i2,
     >                       i1_start,i2_start


*     **** local variables ***
      integer proc_to,proc_from,it
      integer source
      integer msglen
      integer taskid,np
      integer mpierr,n1

c     **** external functions ****
      integer  Parallel2d_convert_taskid_i,Parallel_comm_world
      external Parallel2d_convert_taskid_i,Parallel_comm_world
      
      call Parallel2d_taskid_i(taskid)
      call Parallel2d_np_i(np)


      n1 = int_mb(i1_start(1,nbb)+np) - 1

*     **** pack A(i) array ****
       call D3dB_pfft_index1_copy(n1,int_mb(iq_to_i1(1,nbb)),A,tmp1)

*     **** it = 0, transpose data on same thread ****  
      msglen = int_mb(i2_start(1,nbb)+2-1) - int_mb(i2_start(1,nbb)+1-1)
      call dcopy(2*msglen,
     >           tmp1(int_mb(i1_start(1,nbb)+1-1)),1,
     >           tmp2(int_mb(i2_start(1,nbb)+1-1)),1)
         

*     **** receive packed array data ****
      reqcnt = 0
      do it = 1,np-1

*        **** synchronous receive of tmp ****
         proc_from = mod((taskid)-it+np, np) 

c         msgtype = 1
c         type   = 1
         source=proc_from
         msglen = (int_mb(i2_start(1,nbb)+it+2-1) 
     >          -  int_mb(i2_start(1,nbb)+it+1-1))
         

         if (msglen.gt.0) then
            reqcnt = reqcnt + 1
#ifdef MPI4
            stupid_msglen = msglen
            stupid_type   = msgtype
            stupid_taskid = Parallel2d_convert_taskid_i(source)
            call MPI_IRECV(tmp2(int_mb(i2_start(1,nbb)+it+1-1)),
     >                    stupid_msglen,stupid_complex,
     >                    stupid_taskid,
     >                    stupid_type,stupid_world,
     >                    stupid_request,stupid_ierr)
            request(reqcnt) = 0
            request(reqcnt) = stupid_request
#else
            call MPI_IRECV(tmp2(int_mb(i2_start(1,nbb)+it+1-1)),
     >                    msglen,MPI_DOUBLE_COMPLEX,
     >                    Parallel2d_convert_taskid_i(source),
     >                    msgtype,Parallel_comm_world(),
     >                    request(reqcnt),mpierr) 
#endif
         end if

      end do

*     **** send packed array to other processors ****
      do it = 1,np-1

*        **** synchronous send of tmp ****
         proc_to   = mod((taskid)+it, np) 
         msglen    = (int_mb(i1_start(1,nbb)+it+2-1)
     >              - int_mb(i1_start(1,nbb)+it+1-1))
c         msgtype   = 1
      
         if (msglen.gt.0) then
            reqcnt = reqcnt + 1
#ifdef MPI4
            stupid_msglen = msglen
            stupid_type   = msgtype
            stupid_taskid = Parallel2d_convert_taskid_i(proc_to)
            call MPI_ISEND(tmp1(int_mb(i1_start(1,nbb)+it+1-1)),
     >                     stupid_msglen,stupid_complex,
     >                     stupid_taskid,
     >                     stupid_type,stupid_world,
     >                     stupid_request,stupid_ierr)
            request(reqcnt) = 0 
            request(reqcnt) = stupid_request
#else
            call MPI_ISEND(tmp1(int_mb(i1_start(1,nbb)+it+1-1)),
     >                     msglen,MPI_DOUBLE_COMPLEX,
     >                     Parallel2d_convert_taskid_i(proc_to),
     >                     msgtype,Parallel_comm_world(),
     >                     request(reqcnt),mpierr)
#endif
         end if

      end do
   

      return
      end


*     ***********************************
*     *					*
*     *	   D3dB_c_ptranspose2_jk_end	*
*     *					*
*     ***********************************

*                                           
*      This routine performs the operation   
*               A(i,k,j) <- A(i,j,k)          
*                                              
*      np = the number of worker nodes          
*      proc#=0...(np-1)
*                                                
*       this transpose uses more buffer space     
*       then transpose2                            
*

      subroutine D3dB_c_ptranspose2_jk_end(nbb,A,tmp2,request,reqcnt)

      implicit none
      integer nbb
      complex*16  A(*)
      complex*16  tmp2(*)
      integer     request(*),reqcnt

#include "bafdecls.fh"
#include "D3dB.fh"

*     **** indexing variables ****
      integer iq_to_i1(2,0:3)
      integer iq_to_i2(2,0:3)
      integer iz_to_i2(2,0:3)
      integer i1_start(2,0:3)
      integer i2_start(2,0:3)
      common / ptrans_blk2 / iq_to_i1,iq_to_i2,iz_to_i2,
     >                       i1_start,i2_start


*     **** local variables ***
      logical value
      integer status(2),request2(2)
      integer np,n2,id

      !id=2*(nbb/2)+1      
      if (nbb.lt.2) then
         id = 1
      else
         id = 3
      end if
      call Parallel2d_np_i(np)
      if (np.gt.1) call Parallel_mpiWaitAll(reqcnt,request)

*     **** unpack A(i) array ****
      n2 = int_mb(i2_start(1,nbb)+np) - 1
      call D3dB_pfft_index2_copy(n2,int_mb(iq_to_i2(1,nbb)),tmp2,A)
      call D3dB_pfft_index2_zero(nfft3d(id)-n2,
     >                           int_mb(iz_to_i2(1,nbb)),A)
      return
      end



*     ***********************************
*     *					*
*     *	   D3dB_c_ptranspose_ijk_start	*
*     *					*
*     ***********************************
*                                                
*      This routine performs the operation      
*               A(i,k,j) <- A(i,j,k)           
*                                             
*      np = the number of worker nodes       
*      proc#=0...(np-1)
*                                           
*       this transpose uses more buffer space 
*       then transpose2                      
*


      subroutine D3dB_c_ptranspose_ijk_start(nbb,op,A,tmp1,tmp2,
     >                                       request,reqcnt,msgtype)

      implicit none
      integer nbb,op
      complex*16  A(*)
      complex*16  tmp1(*),tmp2(*)
      integer     request(*),reqcnt
      integer     msgtype

#include "bafdecls.fh"
#include "errquit.fh"
#include "mpif.h"
#include "D3dB.fh"

#ifdef MPI4
#include "stupid_mpi4.fh"
#endif


*     **** indexing variables ****
      integer h_iq_to_i1(2,6,0:3)
      integer h_iq_to_i2(2,6,0:3)
      integer h_iz_to_i2(2,6,0:3)
      integer h_iz_to_i2_count(6,0:3)
      integer h_i1_start(2,6,0:3)
      integer h_i2_start(2,6,0:3)
      common / ptrans_blk_ijk / h_iq_to_i1,
     >                         h_iq_to_i2,
     >                         h_iz_to_i2,
     >                         h_iz_to_i2_count,
     >                         h_i1_start,
     >                         h_i2_start


*     **** local variables ***
      integer n1
      integer proc_to,proc_from,it
      integer source
      integer msglen
      integer taskid,np
      integer mpierr

c     **** external functions ****
      integer  Parallel2d_convert_taskid_i,Parallel_comm_world
      external Parallel2d_convert_taskid_i,Parallel_comm_world
      
      call Parallel2d_taskid_i(taskid)
      call Parallel2d_np_i(np)


*     **** pack A(i) array ****
      n1 = int_mb(h_i1_start(1,op,nbb)+np) - 1
      call D3dB_pfft_index1_copy(n1,int_mb(h_iq_to_i1(1,op,nbb)),A,tmp1)


*     **** it = 0, transpose data on same thread ****  
      msglen = int_mb(h_i2_start(1,op,nbb)+2-1) 
     >       - int_mb(h_i2_start(1,op,nbb)+1-1)
c      call dcopy(2*msglen,
c     >           tmp1(int_mb(h_i1_start(1,op,nbb))),1,
c     >           tmp2(int_mb(h_i2_start(1,op,nbb))),1)
      call Parallel_shared_vector_copy(.false.,2*msglen,
     >           tmp1(int_mb(h_i1_start(1,op,nbb))),
     >           tmp2(int_mb(h_i2_start(1,op,nbb))))


!$OMP MASTER
*     **** receive packed array data ****
      reqcnt = 0
      do it = 1,np-1

*        **** synchronous receive of tmp ****
         proc_from = mod((taskid)-it+np, np) 

c         msgtype = 1
c         type   = 1
         source=proc_from
         msglen = (int_mb(h_i2_start(1,op,nbb)+it+2-1) 
     >          -  int_mb(h_i2_start(1,op,nbb)+it+1-1))
         

         if (msglen.gt.0) then
            reqcnt = reqcnt + 1
#ifdef MPI4
            stupid_msglen = msglen
            stupid_type   = msgtype
            stupid_taskid = Parallel2d_convert_taskid_i(source)
            call MPI_IRECV(tmp2(int_mb(h_i2_start(1,op,nbb)+it+1-1)),
     >                    stupid_msglen,stupid_complex,
     >                    stupid_taskid,
     >                    stupid_type,stupid_world,
     >                    stupid_request,stupid_ierr) 
            request(reqcnt) = 0
            request(reqcnt) = stupid_request
#else
            call MPI_IRECV(tmp2(int_mb(h_i2_start(1,op,nbb)+it+1-1)),
     >                    msglen,MPI_DOUBLE_COMPLEX,
     >                    Parallel2d_convert_taskid_i(source),
     >                    msgtype,Parallel_comm_world(),
     >                    request(reqcnt),mpierr) 
#endif
         end if

      end do

*     **** send packed array to other processors ****
      do it = 1,np-1

*        **** synchronous send of tmp ****
         proc_to   = mod((taskid)+it, np) 
         msglen    = (int_mb(h_i1_start(1,op,nbb)+it+2-1)
     >              - int_mb(h_i1_start(1,op,nbb)+it+1-1))
c         msgtype   = 1
      
         if (msglen.gt.0) then
            reqcnt = reqcnt + 1
#ifdef MPI4
            stupid_msglen = msglen
            stupid_type   = msgtype
            stupid_taskid = Parallel2d_convert_taskid_i(proc_to)
            call MPI_ISEND(tmp1(int_mb(h_i1_start(1,op,nbb)+it+1-1)),
     >                     stupid_msglen,stupid_complex,
     >                     stupid_taskid,
     >                     stupid_type,stupid_world,
     >                     stupid_request,stupid_ierr)
            request(reqcnt) = 0
            request(reqcnt) = stupid_request
#else
            call MPI_ISEND(tmp1(int_mb(h_i1_start(1,op,nbb)+it+1-1)),
     >                     msglen,MPI_DOUBLE_COMPLEX,
     >                     Parallel2d_convert_taskid_i(proc_to),
     >                     msgtype,Parallel_comm_world(),
     >                     request(reqcnt),mpierr)
#endif
         end if

      end do
!$OMP END MASTER

      return
      end



*     ***********************************
*     *					*
*     *	   D3dB_c_ptranspose_ijk_end	*
*     *					*
*     ***********************************
*                                                  
*      This routine performs the operation        
*               A(i,k,j) <- A(i,j,k)             
*                                               
*      np = the number of worker nodes         
*      proc#=0...(np-1)                       
*                                            
*       this transpose uses more buffer space 
*       then transpose2                      
*

      subroutine D3dB_c_ptranspose_ijk_end(nbb,op,A,tmp2,
     >                                 request,reqcnt)
      implicit none
      integer nbb,op
      complex*16  A(*)
      complex*16  tmp2(*)
      integer request(*),reqcnt

#include "bafdecls.fh"
#include "errquit.fh"
#include "D3dB.fh"

*     **** indexing variables ****
      integer h_iq_to_i1(2,6,0:3)
      integer h_iq_to_i2(2,6,0:3)
      integer h_iz_to_i2(2,6,0:3)
      integer h_iz_to_i2_count(6,0:3)
      integer h_i1_start(2,6,0:3)
      integer h_i2_start(2,6,0:3)
      common / ptrans_blk_ijk / h_iq_to_i1,
     >                         h_iq_to_i2,
     >                         h_iz_to_i2,
     >                         h_iz_to_i2_count,
     >                         h_i1_start,
     >                         h_i2_start

*     **** local variables ***
      logical value
      integer n2,n3,np
      
      call Parallel2d_np_i(np)

      if (np.gt.1) call Parallel_mpiWaitAll(reqcnt,request)

*     **** unpack A(i) array ****
      n2 = int_mb(h_i2_start(1,op,nbb)+np) - 1
      n3 = h_iz_to_i2_count(op,nbb)
      call D3dB_pfft_index2_copy(n2,int_mb(h_iq_to_i2(1,op,nbb)),tmp2,A)
      call D3dB_pfft_index2_zero(n3,int_mb(h_iz_to_i2(1,op,nbb)),A)

      return
      end





*     ***********************************
*     *                                 *
*     *    D3dB_c_timereverse_start     *
*     *                                 *
*     ***********************************
*                                           
*      This routine performs the operation   
*            A(i,j,k) <- conjugate(A(i,-j,-k))
*                                              
*      np = the number of worker nodes          
*      proc#=0...(np-1)                          
*                                                 
      subroutine D3dB_c_timereverse_start(nb,A,tmp1,tmp2,
     >                                    request,reqcnt,msgtype)
      implicit none
      integer     nb
      complex*16  A(*)
      complex*16  tmp1(*)
      complex*16  tmp2(*)
      integer     request(*),reqcnt
      integer     msgtype

#include "bafdecls.fh"
#include "errquit.fh"
#include "mpif.h"
#include "D3dB.fh"

#ifdef MPI4
#include "stupid_mpi4.fh"
#endif

*     **** indexing variables ****
c     integer iq_to_i1(2**NFFT2*NSLABS)
c     integer iq_to_i2(2**NFFT2*NSLABS)
c     integer i1_start(NFFT3+1)
c     integer i2_start(NFFT3+1)
      integer iq_to_i1(2,NBLOCKS)
      integer iq_to_i2(2,NBLOCKS)
      integer i1_start(2,NBLOCKS)
      integer i2_start(2,NBLOCKS)
      common / timereverse_blk / iq_to_i1,iq_to_i2,i1_start,i2_start


*     **** local variables ***
      integer i,n1
      integer it
      integer source
      integer msglen
      integer proc_from,proc_to
      integer taskid,np
      integer index1
      integer mpierr

c     **** external functions ****
      integer  Parallel2d_convert_taskid_i,Parallel_comm_world
      external Parallel2d_convert_taskid_i,Parallel_comm_world

      call Parallel2d_taskid_i(taskid)
      call Parallel2d_np_i(np)

      n1 = int_mb(i1_start(1,nb)+np) - 1

*     **** pack A(i) array ****
c      do index1=int_mb(i1_start(1,nb)+1-1),
c     >         (int_mb(i1_start(1,nb)+np+1-1)-1)
c         tmp1(index1) = A(int_mb(iq_to_i1(1,nb)+index1-1))
c      end do
      call D3dB_pfft_index1_copy(n1,int_mb(iq_to_i1(1,nb)),A,tmp1)

*     **** it = 0, transpose data on same thread ****
      msglen = int_mb(i2_start(1,nb)+2-1) - int_mb(i2_start(1,nb)+1-1)
c      do i=1,msglen
c         tmp2(int_mb(i2_start(1,nb)+1-1)+i-1)
c     > = tmp1(int_mb(i1_start(1,nb)+1-1)+i-1)
c      end do
      call Parallel_shared_vector_copy(.false.,2*msglen,
     >                                 tmp1(int_mb(i1_start(1,nb))),
     >                                 tmp2(int_mb(i2_start(1,nb))))


!$OMP MASTER
*     **** receive packed array data ****
      reqcnt = 0
      do it = 1,np-1

*        **** synchronous receive of tmp ****
         proc_from = mod((taskid)-it+np, np)

c         msgtype = 1
c         type   = 1
         source=proc_from
         msglen = (int_mb(i2_start(1,nb)+it+2-1)
     >          -  int_mb(i2_start(1,nb)+it+1-1))


         if (msglen.gt.0) then
            reqcnt = reqcnt + 1
#ifdef MPI4
            stupid_msglen = msglen
            stupid_type   = msgtype
            stupid_taskid = Parallel2d_convert_taskid_i(source)
            call MPI_IRECV(tmp2(int_mb(i2_start(1,nb)+it+1-1)),
     >                    stupid_msglen,stupid_complex,
     >                    stupid_taskid,
     >                    stupid_type,stupid_world,
     >                    stupid_request,stupid_ierr)
            request(reqcnt) = 0
            request(reqcnt) = stupid_request
#else
            call MPI_IRECV(tmp2(int_mb(i2_start(1,nb)+it+1-1)),
     >                    msglen,MPI_DOUBLE_COMPLEX,
     >                    Parallel2d_convert_taskid_i(source),
     >                    msgtype,Parallel_comm_world(),
     >                    request(reqcnt),mpierr)
#endif
         end if

      end do


*     **** send packed array to other processors ****
      do it = 1,np-1

*        **** synchronous send of tmp ****
         proc_to   = mod((taskid)+it, np)
         msglen    = (int_mb(i1_start(1,nb)+it+2-1)
     >              - int_mb(i1_start(1,nb)+it+1-1))
c         msgtype   = 1

         if (msglen.gt.0) then
            reqcnt = reqcnt + 1
#ifdef MPI4
            stupid_msglen = msglen
            stupid_type   = msgtype
            stupid_taskid = Parallel2d_convert_taskid_i(proc_to)
            call MPI_ISEND(tmp1(int_mb(i1_start(1,nb)+it+1-1)),
     >                     stupid_msglen,stupid_complex,
     >                     stupid_taskid,
     >                     stupid_type,stupid_world,
     >                     stupid_request,stupid_ierr)
            request(reqcnt) = 0
            request(reqcnt) = stupid_request
#else
            call MPI_ISEND(tmp1(int_mb(i1_start(1,nb)+it+1-1)),
     >                     msglen,MPI_DOUBLE_COMPLEX,
     >                     Parallel2d_convert_taskid_i(proc_to),
     >                     msgtype,Parallel_comm_world(),
     >                     request(reqcnt),mpierr)
#endif
         end if

      end do
!$OMP END MASTER


      return
      end




*     ***********************************
*     *                                 *
*     *    D3dB_c_timereverse_end       *
*     *                                 *
*     ***********************************
*                                           
*      This routine performs the operation   
*            A(i,j,k) <- conjugate(A(i,-j,-k))
*                                              
*      np = the number of worker nodes          
*      proc#=0...(np-1)                          
*                                                 
      subroutine D3dB_c_timereverse_end(nb,A,tmp1,tmp2,request,reqcnt)
      implicit none
      integer     nb
      complex*16  A(*)
      complex*16  tmp1(*)
      complex*16  tmp2(*)
      integer     request(*),reqcnt

#include "bafdecls.fh"
#include "D3dB.fh"

*     **** indexing variables ****
      integer iq_to_i1(2,NBLOCKS)
      integer iq_to_i2(2,NBLOCKS)
      integer i1_start(2,NBLOCKS)
      integer i2_start(2,NBLOCKS)
      common / timereverse_blk / iq_to_i1,iq_to_i2,i1_start,i2_start


*     **** local variables ***
      logical value
      integer np,n2
      integer index2

      call Parallel2d_np_i(np)
      n2 = (int_mb(i2_start(1,nb)+np)-1)

      if (np.gt.1) call Parallel_mpiWaitAll(reqcnt,request)

c*     **** unpack A(i) array ****
c#ifndef CRAY
c!DIR$ ivdep
c#endif
c      do index2=int_mb(i2_start(1,nb)+1-1),
c     >         (int_mb(i2_start(1,nb)+np+1-1)-1)
c        A(int_mb(iq_to_i2(1,nb)+index2-1))=dconjg(tmp2(index2))
c      end do

!$OMP BARRIER

      call D3dB_pfft_index2_copy_conjg(n2,int_mb(iq_to_i2(1,nb)),tmp2,A)
      return
      end



*     ***********************************
*     *					*
*     *	       D3dB_(c,r,t)_read 	*	
*     *					*
*     ***********************************
      subroutine D3dB_c_read_pio(nb,iunit,A,tmp,jcol)
      implicit none 
      integer nb
      integer iunit
      complex*16 A(*)
      complex*16 tmp(*)
      integer    jcol

#include "bafdecls.fh"
#include "errquit.fh"
#include "mpif.h"
#include "D3dB.fh"

#ifdef MPI4
#include "stupid_mpi4.fh"
#endif

*     *** local variables ***
      logical value,fillcolumn,readcolumn
      integer MASTER,taskid,taskid_i,taskid_j
      parameter(MASTER=0)
      integer p_to, p_here,q
      integer index,k,j
      integer source,msglen,msgtype
      integer tmp1(2),tmp2(2)

c      integer taskid_j,np_j
c      integer ii,jj,jstart,jend
      integer ii
      integer status(2),mpierr

*     **** external functions ****
      integer  Parallel2d_convert_taskid_ij,Parallel2d_comm_i
      external Parallel2d_convert_taskid_ij,Parallel2d_comm_i

      call Parallel_taskid(taskid)
      call Parallel2d_taskid_i(taskid_i)
      call Parallel2d_taskid_j(taskid_j)

c      call Parallel2d_np_j(np_j)
      if (jcol.lt.0) then
c         jstart = 0
c         jend   = np_j-1
         fillcolumn = .true.
         readcolumn = .true.
      else
c         jstart = jcol
c         jend   = jcol
         fillcolumn = (taskid_j.eq.jcol)
         readcolumn = (taskid_j.eq.jcol)
      end if

      if (readcolumn) then

*     **** allocate status memory ****
      if (.not.BA_push_get(mt_int,MPI_STATUS_SIZE*4,
     >                     'status',status(2),status(1)))
     > call errquit('D3dB_c_read:out of stack',0,MA_ERR)


      !**********************
      !**** slab mapping ****
      !**********************
      if (mapping.eq.1) then
*     **** master node reads from file and distributes ****
      if (taskid_i.eq.MASTER) then
         do k=1,nz(nb)

            call dread(iunit,tmp,(nx(nb)+2)*ny(nb))

            call D3dB_ktoqp(nb,k,q,ii)
            !do jj=jstart,jend
               p_to = ii
               !p_to = Parallel2d_convert_taskid_ij(ii,jj)
               if (p_to.eq.MASTER) then
                  index = (q-1)*(nx(nb)/2+1)*ny(nb) + 1
                  call dcopy((nx(nb)+2)*ny(nb),tmp,1,A(index),1)
               else
                  msglen = (nx(nb)/2+1)*ny(nb)
                  msgtype = 79
c                  call SND(9+MSGDBL,tmp,mdtob(2*msglen),p_to,1)
#ifdef MPI4
                  stupid_msglen = msglen
                  stupid_type   = msgtype
                  stupid_taskid = p_to
                  call MPI_SEND(tmp,
     >                     stupid_msglen,stupid_complex,
     >                     stupid_taskid,
     >                     stupid_type,stupid_comm_i,stupid_ierr)
#else
                  call MPI_SEND(tmp,
     >                     msglen,MPI_DOUBLE_COMPLEX,
     >                     p_to,
     >                     msgtype,Parallel2d_comm_i(),mpierr)
#endif
               end if
            !end do
         end do

*     **** not master node ****
      else if (fillcolumn) then
         do k=1,nz(nb)
            call D3dB_ktoqp(nb,k,q,ii)
            p_here = ii
            !p_here = Parallel2d_convert_taskid_ij(ii,taskid_j)
            if (p_here.eq.taskid_i) then

               msglen  = (nx(nb)/2+1)*ny(nb)
               msgtype = 79
               source  = MASTER
c               call RCV(9+MSGDBL,tmp,mdtob(2*msglen),rcv_len,
c     >                  source,rcv_proc,1)
#ifdef MPI4
               stupid_msglen = msglen
               stupid_type   = msgtype
               stupid_taskid = source
               call MPI_RECV(tmp,
     >                    stupid_msglen,stupid_complex,
     >                    stupid_taskid,
     >                    stupid_type,stupid_comm_i,
     >                    int_mb(status(1)),stupid_ierr)
#else
               call MPI_RECV(tmp,
     >                    msglen,MPI_DOUBLE_COMPLEX,
     >                    source,
     >                    msgtype,Parallel2d_comm_i(),
     >                    int_mb(status(1)),mpierr)
#endif
               index = (q-1)*(nx(nb)/2+1)*ny(nb) + 1
               call dcopy((nx(nb)+2)*ny(nb),tmp,1,A(index),1)
               
            end if
         end do
      end if

      !*************************
      !**** hilbert mapping ****
      !*************************
      else
      
*     **** master node reads from file and distributes ****
      if (taskid_i.eq.MASTER) then
         do k=1,nz(nb)
         do j=1,ny(nb)

            call dread(iunit,tmp,(nx(nb)+2))

            q    = int_mb(q_map1(1,nb)+(j-1)+(k-1)*ny(nb))
            ii   = int_mb(p_map1(1,nb)+(j-1)+(k-1)*ny(nb))
            !do jj=jstart,jend
               p_to = ii
               !p_to = Parallel2d_convert_taskid_ij(ii,jj)

               if (p_to.eq.MASTER) then
                  index = (q-1)*(nx(nb)/2+1) + 1
                  call dcopy((nx(nb)+2),tmp,1,A(index),1)
               else
                  msglen  = (nx(nb)/2+1)
                  msgtype = 79
cc                  call SND(9+MSGDBL,tmp,mdtob(2*msglen),p_to,1)
#ifdef MPI4
                  stupid_msglen = msglen
                  stupid_type   = msgtype
                  stupid_taskid = p_to
                  call MPI_SEND(tmp,
     >                     stupid_msglen,stupid_complex,
     >                     stupid_taskid,
     >                     stupid_type,stupid_comm_i,stupid_ierr)
#else
                  call MPI_SEND(tmp,
     >                     msglen,MPI_DOUBLE_COMPLEX,
     >                     p_to,
     >                     msgtype,Parallel2d_comm_i(),mpierr)
#endif
               end if
            !end do
         end do
         end do

*     **** not master node ****
      else if (fillcolumn) then
         do k=1,nz(nb)
         do j=1,ny(nb)
            
            q      = int_mb(q_map1(1,nb)+(j-1)+(k-1)*ny(nb))
            ii     = int_mb(p_map1(1,nb)+(j-1)+(k-1)*ny(nb))
            p_here = ii
            !p_here = Parallel2d_convert_taskid_ij(ii,taskid_j)
            if (p_here.eq.taskid_i) then
               msglen = (nx(nb)/2+1)
               msgtype = 79
               source  = MASTER
c               call RCV(9+MSGDBL,tmp,mdtob(2*msglen),rcv_len,
c     >                  source,rcv_proc,1)
#ifdef MPI4
               stupid_msglen = msglen
               stupid_type   = msgtype
               stupid_taskid = source
               call MPI_RECV(tmp,
     >                    stupid_msglen,stupid_complex,
     >                    stupid_taskid,
     >                    stupid_type,stupid_comm_i,
     >                    int_mb(status(1)),stupid_ierr)
#else
               call MPI_RECV(tmp,
     >                    msglen,MPI_DOUBLE_COMPLEX,
     >                    source,
     >                    msgtype,Parallel2d_comm_i(),
     >                    int_mb(status(1)),mpierr)
#endif
               index = (q-1)*(nx(nb)/2+1) + 1
               call dcopy((nx(nb)+2),tmp,1,A(index),1)
            end if
         end do
         end do
      end if
      

*     **** allocate temporary space  ****
      value = BA_push_get(mt_dcpl,(nfft3d(nb)),'tmp1',
     >                    tmp1(2),tmp1(1))
      value = value.and.
     >      BA_push_get(mt_dcpl,(nfft3d(nb)),'tmp2',tmp2(2),tmp2(1))
      if (.not. value) 
     > call errquit('D3dB_c_read:out of stack',0,MA_ERR)
*
      call D3dB_c_transpose_ijk(nb,5,A,
     >                          dcpl_mb(tmp1(1)),
     >                          dcpl_mb(tmp2(1)))  !*** map1to3 operation ***

*     **** deallocate temporary space  ****
      value = BA_pop_stack(tmp2(2))
      value = value.and.BA_pop_stack(tmp1(2))
      if (.not. value) 
     > call errquit('D3dB_c_read:popping stack',0,MA_ERR)

      end if

*     *** deallocate status memory ***
      if (.not.BA_pop_stack(status(2)))
     > call errquit('D3dB_c_read:popping stack',0,MA_ERR)



      !*** shift filepointer by (nx(nb)+2)*ny(nb)*nz(nb) doubles ****
      else
         if (taskid_i.eq.MASTER) 
     >      call dshift_fileptr(iunit,(nx(nb)+2)*ny(nb)*nz(nb))
      end if !*** readcolumn ***


      return
      end

      subroutine D3dB_c_read(nb,iunit,A,tmp,jcol)
      implicit none 
      integer nb
      integer iunit
      complex*16 A(*)
      complex*16 tmp(*)
      integer    jcol

#include "bafdecls.fh"
#include "errquit.fh"
#include "mpif.h"
#include "D3dB.fh"

#ifdef MPI4
#include "stupid_mpi4.fh"
#endif

*     *** local variables ***
      logical value,fillcolumn
      integer MASTER,taskid
      parameter(MASTER=0)
      integer p_to, p_here,q
      integer index,k,j
      integer source,msglen,msgtype
      integer tmp1(2),tmp2(2)

      integer taskid_j,np_j
      integer ii,jj,jstart,jend
      integer status(2),mpierr

*     **** external functions ****
      integer  Parallel2d_convert_taskid_ij,Parallel_comm_world
      external Parallel2d_convert_taskid_ij,Parallel_comm_world

      call Parallel_taskid(taskid)
      call Parallel2d_taskid_j(taskid_j)

      call Parallel2d_np_j(np_j)
      if (jcol.lt.0) then
         jstart = 0
         jend   = np_j-1
         fillcolumn = .true.
      else
         jstart = jcol
         jend   = jcol
         fillcolumn = (taskid_j.eq.jcol)
      end if


*     **** allocate status memory ****
      if (.not.BA_push_get(mt_int,MPI_STATUS_SIZE*4,
     >                     'status',status(2),status(1)))
     > call errquit('D3dB_c_read:out of stack',0,MA_ERR)


      !**********************
      !**** slab mapping ****
      !**********************
      if (mapping.eq.1) then
*     **** master node reads from file and distributes ****
      if (taskid.eq.MASTER) then
         do k=1,nz(nb)

            call dread(iunit,tmp,(nx(nb)+2)*ny(nb))

            call D3dB_ktoqp(nb,k,q,ii)
            do jj=jstart,jend
               p_to = Parallel2d_convert_taskid_ij(ii,jj)
               if (p_to.eq.MASTER) then
                  index = (q-1)*(nx(nb)/2+1)*ny(nb) + 1
                  call dcopy((nx(nb)+2)*ny(nb),tmp,1,A(index),1)
               else
                  msglen = (nx(nb)/2+1)*ny(nb)
                  msgtype = 89
c                  call SND(9+MSGDBL,tmp,mdtob(2*msglen),p_to,1)
#ifdef MPI4
                  stupid_msglen = msglen
                  stupid_type   = msgtype
                  stupid_taskid = p_to
                  call MPI_SEND(tmp,
     >                     stupid_msglen,stupid_complex,
     >                     stupid_taskid,
     >                     stupid_type,stupid_world,stupid_ierr)
#else
                  call MPI_SEND(tmp,
     >                     msglen,MPI_DOUBLE_COMPLEX,
     >                     p_to,
     >                     msgtype,Parallel_comm_world(),mpierr)
#endif
               end if
            end do
         end do

*     **** not master node ****
      else if (fillcolumn) then
         do k=1,nz(nb)
            call D3dB_ktoqp(nb,k,q,ii)
            p_here = Parallel2d_convert_taskid_ij(ii,taskid_j)
            if (p_here.eq.taskid) then

               msglen  = (nx(nb)/2+1)*ny(nb)
               msgtype = 89
               source  = MASTER
c               call RCV(9+MSGDBL,tmp,mdtob(2*msglen),rcv_len,
c     >                  source,rcv_proc,1)
#ifdef MPI4
               stupid_msglen = msglen
               stupid_type   = msgtype
               stupid_taskid = source
               call MPI_RECV(tmp,
     >                    stupid_msglen,stupid_complex,
     >                    stupid_taskid,
     >                    stupid_type,stupid_world,
     >                    int_mb(status(1)),stupid_ierr)
#else
               call MPI_RECV(tmp,
     >                    msglen,MPI_DOUBLE_COMPLEX,
     >                    source,
     >                    msgtype,Parallel_comm_world(),
     >                    int_mb(status(1)),mpierr)
#endif
               index = (q-1)*(nx(nb)/2+1)*ny(nb) + 1
               call dcopy((nx(nb)+2)*ny(nb),tmp,1,A(index),1)
               
            end if
         end do
      end if

      !*************************
      !**** hilbert mapping ****
      !*************************
      else
      
*     **** master node reads from file and distributes ****
      if (taskid.eq.MASTER) then
         do k=1,nz(nb)
         do j=1,ny(nb)

            call dread(iunit,tmp,(nx(nb)+2))

            q    = int_mb(q_map1(1,nb)+(j-1)+(k-1)*ny(nb))
            ii   = int_mb(p_map1(1,nb)+(j-1)+(k-1)*ny(nb))
            do jj=jstart,jend
               p_to = Parallel2d_convert_taskid_ij(ii,jj)

               if (p_to.eq.MASTER) then
                  index = (q-1)*(nx(nb)/2+1) + 1
                  call dcopy((nx(nb)+2),tmp,1,A(index),1)
               else
                  msglen  = (nx(nb)/2+1)
                  msgtype = 89
cc                  call SND(9+MSGDBL,tmp,mdtob(2*msglen),p_to,1)
#ifdef MPI4
                  stupid_msglen = msglen
                  stupid_type   = msgtype
                  stupid_taskid = p_to
                  call MPI_SEND(tmp,
     >                     stupid_msglen,stupid_complex,
     >                     stupid_taskid,
     >                     stupid_type,stupid_world,stupid_ierr)
#else
                  call MPI_SEND(tmp,
     >                     msglen,MPI_DOUBLE_COMPLEX,
     >                     p_to,
     >                     msgtype,Parallel_comm_world(),mpierr)
#endif
               end if
            end do
         end do
         end do

*     **** not master node ****
      else if (fillcolumn) then
         do k=1,nz(nb)
         do j=1,ny(nb)
            
            q      = int_mb(q_map1(1,nb)+(j-1)+(k-1)*ny(nb))
            ii     = int_mb(p_map1(1,nb)+(j-1)+(k-1)*ny(nb))
            p_here = Parallel2d_convert_taskid_ij(ii,taskid_j)
            if (p_here.eq.taskid) then
               msglen = (nx(nb)/2+1)
               msgtype = 89
               source  = MASTER
c               call RCV(9+MSGDBL,tmp,mdtob(2*msglen),rcv_len,
c     >                  source,rcv_proc,1)
#ifdef MPI4
               stupid_msglen = msglen
               stupid_type   = msgtype
               stupid_taskid = source
               call MPI_RECV(tmp,
     >                    stupid_msglen,stupid_complex,
     >                    stupid_taskid,
     >                    stupid_type,stupid_world,
     >                    int_mb(status(1)),stupid_ierr)
#else
               call MPI_RECV(tmp,
     >                    msglen,MPI_DOUBLE_COMPLEX,
     >                    source,
     >                    msgtype,Parallel_comm_world(),
     >                    int_mb(status(1)),mpierr)
#endif
               index = (q-1)*(nx(nb)/2+1) + 1
               call dcopy((nx(nb)+2),tmp,1,A(index),1)
            end if
         end do
         end do
      end if
      

*     **** allocate temporary space  ****
      value = BA_push_get(mt_dcpl,(nfft3d(nb)),'tmp1',
     >                    tmp1(2),tmp1(1))
      value = value.and.
     >      BA_push_get(mt_dcpl,(nfft3d(nb)),'tmp2',tmp2(2),tmp2(1))
      if (.not. value) 
     > call errquit('D3dB_c_read:out of stack',0,MA_ERR)
*
      call D3dB_c_transpose_ijk(nb,5,A,
     >                          dcpl_mb(tmp1(1)),
     >                          dcpl_mb(tmp2(1)))  !*** map1to3 operation ***

*     **** deallocate temporary space  ****
      value = BA_pop_stack(tmp2(2))
      value = value.and.BA_pop_stack(tmp1(2))
      if (.not. value) 
     > call errquit('D3dB_c_read:popping stack',0,MA_ERR)

      end if

*     *** deallocate status memory ***
      if (.not.BA_pop_stack(status(2)))
     > call errquit('D3dB_c_read:popping stack',0,MA_ERR)


      return
      end


      subroutine D3dB_r_read_pio(nb,iunit,A,tmp,jcol)
      implicit none 
      integer nb
      integer iunit
      real*8  A(*)
      real*8  tmp(*)
      integer jcol

#include "bafdecls.fh"
#include "errquit.fh"
#include "mpif.h"
#include "D3dB.fh"

#ifdef MPI4
#include "stupid_mpi4.fh"
#endif


*     *** local variables ***
      logical fillcolumn,readcolumn
      integer MASTER,taskid,taskid_i,taskid_j
      parameter(MASTER=0)
      integer p_to, p_here,q
      integer j,k,index,index2,mpierr
      integer source,msglen,msgtype,status(2)

c      integer taskid_j,np_j
c      integer ii,jj,jstart,jend
      integer ii

*     **** external functions ****
      integer  Parallel2d_convert_taskid_ij,Parallel2d_comm_i
      external Parallel2d_convert_taskid_ij,Parallel2d_comm_i

      call Parallel_taskid(taskid)
      call Parallel2d_taskid_i(taskid_i)
      call Parallel2d_taskid_j(taskid_j)

c      call Parallel2d_np_j(np_j)
      if (jcol.lt.0) then
c         jstart = 0
c         jend   = np_j-1
         fillcolumn = .true.
         readcolumn = .true.
      else
c         jstart = jcol
c         jend   = jcol
         fillcolumn = (taskid_j.eq.jcol)
         readcolumn = (taskid_j.eq.jcol)
      end if

      if (readcolumn) then

c     **** allocate status memory ****
      if (.not.BA_push_get(mt_int,MPI_STATUS_SIZE*4,
     >                     'status',status(2),status(1)))
     > call errquit('D3dB_r_read:out of stack',0,MA_ERR)


      !**********************
      !**** slab mapping ****
      !**********************
      if (mapping.eq.1) then
*     **** master node reads from file and distributes ****
      if (taskid_i.eq.MASTER) then
         do k=1,nz(nb)

            call dread(iunit,tmp,(nx(nb))*ny(nb))

            call D3dB_ktoqp(nb,k,q,ii)
            !do jj=jstart,jend
               p_to = ii
               !p_to = Parallel2d_convert_taskid_ij(ii,jj)
               if (p_to.eq.MASTER) then
                  do j=1,ny(nb)
                     index = (q-1)*(nx(nb)+2)*ny(nb) 
     >                     + (j-1)*(nx(nb)+2) + 1
                     index2 = (j-1)*nx(nb) + 1
                     call dcopy(nx(nb),tmp(index2),1,A(index),1)
                     A(index+nx(nb)) = 0.0d0
                     A(index+nx(nb)+1) = 0.0d0
                  end do
               else
                  msglen = (nx(nb))*ny(nb)
                  msgtype = 99
c                  call SND(9+MSGDBL,tmp,mdtob(msglen),p_to,1)
#ifdef MPI4
                  stupid_msglen = msglen
                  stupid_type   = msgtype
                  stupid_taskid = p_to
                  call MPI_SEND(tmp,
     >                     stupid_msglen,stupid_double,
     >                     stupid_taskid,
     >                     stupid_type,stupid_comm_i,stupid_ierr)
#else
                  call MPI_SEND(tmp,
     >                     msglen,MPI_DOUBLE_PRECISION,
     >                     p_to,
     >                     msgtype,Parallel2d_comm_i,mpierr)
#endif
               end if
            !end do
         end do

*     **** not master node ****
      else if (fillcolumn) then
         do k=1,nz(nb)

            call D3dB_ktoqp(nb,k,q,ii)
            p_here = ii
            !p_here = Parallel2d_convert_taskid_ij(ii,taskid_j)
            if (p_here.eq.taskid_i) then

               msglen  = (nx(nb))*ny(nb)
               msgtype = 99
               source  = MASTER
c               call RCV(9+MSGDBL,tmp,mdtob(msglen),rcv_len,
c     >                       source,rcv_proc,1)
#ifdef MPI4
               stupid_msglen = msglen
               stupid_type   = msgtype
               stupid_taskid = source
               call MPI_RECV(tmp,
     >                    stupid_msglen,stupid_double,
     >                    stupid_taskid,
     >                    stupid_type,stupid_comm_i,
     >                    int_mb(status(1)),stupid_ierr)
#else
               call MPI_RECV(tmp,
     >                    msglen,MPI_DOUBLE_PRECISION,
     >                    source,
     >                    msgtype,Parallel2d_comm_i(),
     >                    int_mb(status(1)),mpierr)
#endif
               do j=1,ny(nb)
                  index = (q-1)*(nx(nb)+2)*ny(nb) 
     >                  + (j-1)*(nx(nb)+2) + 1
                  index2 = (j-1)*nx(nb) + 1
                  call dcopy(nx(nb),tmp(index2),1,A(index),1)
                  A(index+nx(nb)) = 0.0d0
                  A(index+nx(nb)+1) = 0.0d0
               end do
            end if
         end do
      end if



      !*************************
      !**** hilbert mapping ****
      !*************************
      else
*     **** master node reads from file and distributes ****
      if (taskid_i.eq.MASTER) then
         do k=1,nz(nb)
         do j=1,ny(nb)

            call dread(iunit,tmp,(nx(nb)))
            tmp(nx(nb)+1) = 0.0d0
            tmp(nx(nb)+2) = 0.0d0

            q    = int_mb(q_map1(1,nb)+(j-1)+(k-1)*ny(nb))
            ii   = int_mb(p_map1(1,nb)+(j-1)+(k-1)*ny(nb))
            !do jj=jstart,jend
               p_to = ii
               !p_to = Parallel2d_convert_taskid_ij(ii,jj)
               if (p_to.eq.MASTER) then
                  index = (q-1)*(nx(nb)+2) + 1
                  call dcopy((nx(nb)+2),tmp,1,A(index),1)
               else
                  msglen  = (nx(nb)+2)
                  msgtype = 99
c                  call SND(9+MSGDBL,tmp,mdtob(msglen),p_to,1)
#ifdef MPI4
                  stupid_msglen = msglen
                  stupid_type   = msgtype
                  stupid_taskid = p_to
                  call MPI_SEND(tmp,
     >                     stupid_msglen,stupid_double,
     >                     stupid_taskid,
     >                     stupid_type,stupid_comm_i,stupid_ierr)
#else
                  call MPI_SEND(tmp,
     >                     msglen,MPI_DOUBLE_PRECISION,
     >                     p_to,
     >                     msgtype,Parallel2d_comm_i(),mpierr)
#endif

               end if
            !end do
         end do
         end do

*     **** not master node ****
      else if (fillcolumn) then
         do k=1,nz(nb)
         do j=1,ny(nb)

            q      = int_mb(q_map1(1,nb)+(j-1)+(k-1)*ny(nb))
            ii     = int_mb(p_map1(1,nb)+(j-1)+(k-1)*ny(nb))
            p_here = ii
            !p_here = Parallel2d_convert_taskid_ij(ii,taskid_j)
            if (p_here.eq.taskid_i) then
               msglen  = (nx(nb)+2)
               msgtype = 99
               source  = MASTER
c               call RCV(9+MSGDBL,tmp,mdtob(msglen),rcv_len,
c     >                  source,rcv_proc,1)
#ifdef MPI4
               stupid_msglen = msglen
               stupid_type   = msgtype
               stupid_taskid = source
               call MPI_RECV(tmp,
     >                    stupid_msglen,stupid_double,
     >                    stupid_taskid,
     >                    stupid_type,stupid_comm_i,
     >                    int_mb(status(1)),stupid_ierr)
#else
               call MPI_RECV(tmp,
     >                    msglen,MPI_DOUBLE_PRECISION,
     >                    source,
     >                    msgtype,Parallel2d_comm_i(),
     >                    int_mb(status(1)),mpierr)
#endif
               index = (q-1)*(nx(nb)+2) + 1
               call dcopy((nx(nb)+2),tmp,1,A(index),1)

            end if
         end do
         end do
       end if

      end if

*     *** deallocate status memory ***
      if (.not.BA_pop_stack(status(2)))
     > call errquit('D3dB_r_read:popping stack',0,MA_ERR)


      !*** shift filepointer by nx(nb)*ny(nb)*nz(nb) doubles ****
      else
         if (taskid_i.eq.MASTER) 
     >      call dshift_fileptr(iunit,nx(nb)*ny(nb)*nz(nb))
      end if !*** readcolumn ***

      return
      end


      subroutine D3dB_r_read(nb,iunit,A,tmp,jcol)
      implicit none 
      integer nb
      integer iunit
      real*8  A(*)
      real*8  tmp(*)
      integer jcol

#include "bafdecls.fh"
#include "errquit.fh"
#include "mpif.h"
#include "D3dB.fh"

#ifdef MPI4
#include "stupid_mpi4.fh"
#endif


*     *** local variables ***
      logical fillcolumn
      integer MASTER,taskid
      parameter(MASTER=0)
      integer p_to, p_here,q
      integer j,k,index,index2,mpierr
      integer source,msglen,msgtype,status(2)

      integer taskid_j,np_j
      integer ii,jj,jstart,jend

*     **** external functions ****
      integer  Parallel2d_convert_taskid_ij,Parallel_comm_world
      external Parallel2d_convert_taskid_ij,Parallel_comm_world

      call Parallel_taskid(taskid)
      call Parallel2d_taskid_j(taskid_j)
      call Parallel2d_np_j(np_j)
      if (jcol.lt.0) then
         jstart = 0
         jend   = np_j-1
         fillcolumn = .true.
      else
         jstart = jcol
         jend   = jcol
         fillcolumn = (taskid_j.eq.jcol)
      end if


c     **** allocate status memory ****
      if (.not.BA_push_get(mt_int,MPI_STATUS_SIZE*4,
     >                     'status',status(2),status(1)))
     > call errquit('D3dB_r_read:out of stack',0,MA_ERR)


      !**********************
      !**** slab mapping ****
      !**********************
      if (mapping.eq.1) then
*     **** master node reads from file and distributes ****
      if (taskid.eq.MASTER) then
         do k=1,nz(nb)

            call dread(iunit,tmp,(nx(nb))*ny(nb))

            call D3dB_ktoqp(nb,k,q,ii)
            do jj=jstart,jend
               p_to = Parallel2d_convert_taskid_ij(ii,jj)
               if (p_to.eq.MASTER) then
                  do j=1,ny(nb)
                     index = (q-1)*(nx(nb)+2)*ny(nb) 
     >                     + (j-1)*(nx(nb)+2) + 1
                     index2 = (j-1)*nx(nb) + 1
                     call dcopy(nx(nb),tmp(index2),1,A(index),1)
                     A(index+nx(nb)) = 0.0d0
                     A(index+nx(nb)+1) = 0.0d0
                  end do
               else
                  msglen = (nx(nb))*ny(nb)
                  msgtype = 109
c                  call SND(9+MSGDBL,tmp,mdtob(msglen),p_to,1)
#ifdef MPI4
                  stupid_msglen = msglen
                  stupid_type   = msgtype
                  stupid_taskid = p_to
                  call MPI_SEND(tmp,
     >                     stupid_msglen,stupid_double,
     >                     stupid_taskid,
     >                     stupid_type,stupid_world,stupid_ierr)
#else
                  call MPI_SEND(tmp,
     >                     msglen,MPI_DOUBLE_PRECISION,
     >                     p_to,
     >                     msgtype,Parallel_comm_world(),mpierr)
#endif
               end if
            end do
         end do

*     **** not master node ****
      else if (fillcolumn) then
         do k=1,nz(nb)

            call D3dB_ktoqp(nb,k,q,ii)
            p_here = Parallel2d_convert_taskid_ij(ii,taskid_j)
            if (p_here.eq.taskid) then

               msglen  = (nx(nb))*ny(nb)
               msgtype = 109
               source  = MASTER
c               call RCV(9+MSGDBL,tmp,mdtob(msglen),rcv_len,
c     >                       source,rcv_proc,1)
#ifdef MPI4
               stupid_msglen = msglen
               stupid_type   = msgtype
               stupid_taskid = source
               call MPI_RECV(tmp,
     >                    stupid_msglen,stupid_double,
     >                    stupid_taskid,
     >                    stupid_type,stupid_world,
     >                    int_mb(status(1)),stupid_ierr)
#else
               call MPI_RECV(tmp,
     >                    msglen,MPI_DOUBLE_PRECISION,
     >                    source,
     >                    msgtype,Parallel_comm_world(),
     >                    int_mb(status(1)),mpierr)
#endif
               do j=1,ny(nb)
                  index = (q-1)*(nx(nb)+2)*ny(nb) 
     >                  + (j-1)*(nx(nb)+2) + 1
                  index2 = (j-1)*nx(nb) + 1
                  call dcopy(nx(nb),tmp(index2),1,A(index),1)
                  A(index+nx(nb)) = 0.0d0
                  A(index+nx(nb)+1) = 0.0d0
               end do
            end if
         end do
      end if



      !*************************
      !**** hilbert mapping ****
      !*************************
      else
*     **** master node reads from file and distributes ****
      if (taskid.eq.MASTER) then
         do k=1,nz(nb)
         do j=1,ny(nb)

            call dread(iunit,tmp,(nx(nb)))
            tmp(nx(nb)+1) = 0.0d0
            tmp(nx(nb)+2) = 0.0d0

            q    = int_mb(q_map1(1,nb)+(j-1)+(k-1)*ny(nb))
            ii   = int_mb(p_map1(1,nb)+(j-1)+(k-1)*ny(nb))
            do jj=jstart,jend
               p_to = Parallel2d_convert_taskid_ij(ii,jj)
               if (p_to.eq.MASTER) then
                  index = (q-1)*(nx(nb)+2) + 1
                  call dcopy((nx(nb)+2),tmp,1,A(index),1)
               else
                  msglen  = (nx(nb)+2)
                  msgtype = 109
c                  call SND(9+MSGDBL,tmp,mdtob(msglen),p_to,1)
#ifdef MPI4
                  stupid_msglen = msglen
                  stupid_type   = msgtype
                  stupid_taskid = p_to
                  call MPI_SEND(tmp,
     >                     stupid_msglen,stupid_double,
     >                     stupid_taskid,
     >                     stupid_type,stupid_world,stupid_ierr)
#else
                  call MPI_SEND(tmp,
     >                     msglen,MPI_DOUBLE_PRECISION,
     >                     p_to,
     >                     msgtype,Parallel_comm_world(),mpierr)
#endif

               end if
            end do
         end do
         end do

*     **** not master node ****
      else if (fillcolumn) then
         do k=1,nz(nb)
         do j=1,ny(nb)

            q      = int_mb(q_map1(1,nb)+(j-1)+(k-1)*ny(nb))
            ii     = int_mb(p_map1(1,nb)+(j-1)+(k-1)*ny(nb))
            p_here = Parallel2d_convert_taskid_ij(ii,taskid_j)
            if (p_here.eq.taskid) then
               msglen  = (nx(nb)+2)
               msgtype = 109
               source  = MASTER
c               call RCV(9+MSGDBL,tmp,mdtob(msglen),rcv_len,
c     >                  source,rcv_proc,1)
#ifdef MPI4
               stupid_msglen = msglen
               stupid_type   = msgtype
               stupid_taskid = source
               call MPI_RECV(tmp,
     >                    stupid_msglen,stupid_double,
     >                    stupid_taskid,
     >                    stupid_type,stupid_world,
     >                    int_mb(status(1)),stupid_ierr)
#else
               call MPI_RECV(tmp,
     >                    msglen,MPI_DOUBLE_PRECISION,
     >                    source,
     >                    msgtype,Parallel_comm_world(),
     >                    int_mb(status(1)),mpierr)
#endif
               index = (q-1)*(nx(nb)+2) + 1
               call dcopy((nx(nb)+2),tmp,1,A(index),1)

            end if
         end do
         end do
       end if

      end if

*     *** deallocate status memory ***
      if (.not.BA_pop_stack(status(2)))
     > call errquit('D3dB_r_read:popping stack',0,MA_ERR)


      return
      end




      subroutine D3dB_t_read_pio(nb,iunit,A,tmp,jcol)
      implicit none 
      integer nb
      integer iunit
      real*8  A(*)
      real*8  tmp(*)
      integer jcol

#include "bafdecls.fh"
#include "errquit.fh"
#include "mpif.h"
#include "D3dB.fh"

#ifdef MPI4
#include "stupid_mpi4.fh"
#endif

*     *** local variables ***
      integer MASTER,taskid,taskid_i,taskid_j
      parameter(MASTER=0)

      logical value,fillcolumn,readcolumn
      integer p_to, p_here,q
      integer j,k,index
      integer source,msglen,msgtype,status(2),mpierr
      integer tmp1(2),tmp2(2)

c      integer taskid_j,np_j
c      integer ii,jj,jstart,jend
      integer ii

*     **** external functions ****
      integer  Parallel2d_convert_taskid_ij,Parallel2d_comm_i
      external Parallel2d_convert_taskid_ij,Parallel2d_comm_i

      call Parallel_taskid(taskid)
      call Parallel2d_taskid_i(taskid_i)
      call Parallel2d_taskid_j(taskid_j)

c      call Parallel2d_np_j(np_j)
      if (jcol.lt.0) then
c         jstart = 0
c         jend   = np_j-1
         fillcolumn = .true.
         readcolumn = .true.
      else
c         jstart = jcol
c         jend   = jcol
         fillcolumn = (taskid_j.eq.jcol)
         readcolumn = (taskid_j.eq.jcol)
      end if


      if (readcolumn) then

*     **** allocate status memory ****
      if (.not.BA_push_get(mt_int,MPI_STATUS_SIZE*4,
     >                     'status',status(2),status(1)))
     > call errquit('D3dB_t_read:out of stack',0,MA_ERR)


      if (mapping.eq.1) then
*        **** master node reads from file and distributes ****
         if (taskid_i.eq.MASTER) then
         do k=1,nz(nb)

            call dread(iunit,tmp,(nx(nb)/2+1)*ny(nb))

            call D3dB_ktoqp(nb,k,q,ii)
            
            !do jj=jstart,jend

              p_to = ii
              !p_to = Parallel2d_convert_taskid_ij(ii,jj)
              if (p_to.eq.MASTER) then
                 index = (q-1)*(nx(nb)/2+1)*ny(nb) + 1
                 call dcopy((nx(nb)/2+1)*ny(nb),tmp,1,A(index),1)
              else
                 msglen = (nx(nb)/2+1)*ny(nb)
                 msgtype = 119
c                 call SND(9+MSGDBL,tmp,mdtob(msglen),p_to,1)
#ifdef MPI4
                 stupid_msglen = msglen
                 stupid_type   = msgtype
                 stupid_taskid = p_to
                 call MPI_SEND(tmp,
     >                     stupid_msglen,stupid_double,
     >                     stupid_taskid,
     >                     stupid_type,stupid_comm_i,stupid_ierr)
#else
                 call MPI_SEND(tmp,
     >                     msglen,MPI_DOUBLE_PRECISION,
     >                     p_to,
     >                     msgtype,Parallel2d_comm_i(),mpierr)
#endif
              end if
            !end do
         end do

*        **** not master node ****
         else if (fillcolumn) then

         do k=1,nz(nb)
            call D3dB_ktoqp(nb,k,q,ii)
            p_here = ii
            !p_here = Parallel2d_convert_taskid_ij(ii,taskid_j)
            if (p_here.eq.taskid_i) then
               msglen  = (nx(nb)/2+1)*ny(nb)
               msgtype = 119
               source  = MASTER
c               call RCV(9+MSGDBL,tmp,mdtob(msglen),rcv_len,
c     >                       source,rcv_proc,1)
#ifdef MPI4
               stupid_msglen = msglen
               stupid_type   = msgtype
               stupid_taskid = source
               call MPI_RECV(tmp,
     >                    stupid_msglen,stupid_double,
     >                    stupid_taskid,
     >                    stupid_type,stupid_comm_i,
     >                    int_mb(status(1)),stupid_ierr)
#else
               call MPI_RECV(tmp,
     >                    msglen,MPI_DOUBLE_PRECISION,
     >                    source,
     >                    msgtype,Parallel2d_comm_i(),
     >                    int_mb(status(1)),mpierr)
#endif
               index = (q-1)*(nx(nb)/2+1)*ny(nb) + 1
               call dcopy((nx(nb)/2+1)*ny(nb),tmp,1,A(index),1)
               
            end if
         end do
         end if


      !*************************
      !**** hilbert mapping ****
      !*************************
      else
*     **** master node reads from file and distributes ****
      if (taskid_i.eq.MASTER) then
         do k=1,nz(nb)
         do j=1,ny(nb)

            call dread(iunit,tmp,(nx(nb)/2+1))

            q   = int_mb(q_map1(1,nb)+(j-1)+(k-1)*ny(nb))
            ii  = int_mb(p_map1(1,nb)+(j-1)+(k-1)*ny(nb))
            !do jj=jstart,jend

               p_to = ii
               !p_to = Parallel2d_convert_taskid_ij(ii,jj)
               if (p_to.eq.MASTER) then
                  index = (q-1)*(nx(nb)/2+1) + 1
                  call dcopy((nx(nb)/2+1),tmp,1,A(index),1)
               else
                  msglen  = (nx(nb)/2+1)
                  msgtype = 119
c                  call SND(9+MSGDBL,tmp,mdtob(msglen),p_to,1)
#ifdef MPI4
                  stupid_msglen = msglen
                  stupid_type   = msgtype
                  stupid_taskid = p_to
                  call MPI_SEND(tmp,
     >                     stupid_msglen,stupid_double,
     >                     stupid_taskid,
     >                     stupid_type,stupid_comm_i,stupid_ierr)
#else
                  call MPI_SEND(tmp,
     >                     msglen,MPI_DOUBLE_PRECISION,
     >                     p_to,
     >                     msgtype,Parallel2d_comm_i(),mpierr)
#endif
               end if
            !end do
         end do
         end do

*     **** not master node ****
      else if (fillcolumn) then
         do k=1,nz(nb)
         do j=1,ny(nb)

            q      = int_mb(q_map1(1,nb)+(j-1)+(k-1)*ny(nb))
            ii     = int_mb(p_map1(1,nb)+(j-1)+(k-1)*ny(nb))
            p_here = ii
            !p_here = Parallel2d_convert_taskid_ij(ii,taskid_j)
            if (p_here.eq.taskid_i) then
               msglen  = (nx(nb)/2+1)
               msgtype = 119
               source  = MASTER
c               call RCV(9+MSGDBL,tmp,mdtob(msglen),rcv_len,
c     >                  source,rcv_proc,1)
#ifdef MPI4
               stupid_msglen = msglen
               stupid_type   = msgtype
               stupid_taskid = source
               call MPI_RECV(tmp,
     >                    stupid_msglen,stupid_double,
     >                    stupid_taskid,
     >                    stupid_type,stupid_comm_i,
     >                    int_mb(status(1)),stupid_ierr)
#else
               call MPI_RECV(tmp,
     >                    msglen,MPI_DOUBLE_PRECISION,
     >                    source,
     >                    msgtype,Parallel2d_comm_i(),
     >                    int_mb(status(1)),mpierr)
#endif
               index = (q-1)*(nx(nb)/2+1) + 1
               call dcopy((nx(nb)/2+1),tmp,1,A(index),1)

            end if
         end do
         end do
      end if

*     **** allocate temporary space  ****
      value = BA_push_get(mt_dcpl,(nfft3d(nb)),'tmp1',
     >                    tmp1(2),tmp1(1))
      value = value.and.
     >      BA_push_get(mt_dcpl,(nfft3d(nb)),'tmp2',tmp2(2),tmp2(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)
*
      call D3dB_t_transpose_ijk(nb,5,A,
     >                          dcpl_mb(tmp1(1)),
     >                          dcpl_mb(tmp2(1)))  !*** map1to3 operation ***

*     **** deallocate temporary space  ****
      value =           BA_pop_stack(tmp2(2))
      value = value.and.BA_pop_stack(tmp1(2))
      if (.not. value) call errquit('error popping stack 1',0, MA_ERR)

      end if

*     *** deallocate status memory ***
      if (.not.BA_pop_stack(status(2)))
     > call errquit('D3dB_t_read:popping stack',0,MA_ERR)


      !*** shift filepointer by (nx(nb)/2+1)*ny(nb)*nz(nb) doubles ****
      else
         if (taskid_i.eq.MASTER) 
     >      call dshift_fileptr(iunit,(nx(nb)/2+1)*ny(nb)*nz(nb))
      end if !*** readcolumn ***

      return
      end


      subroutine D3dB_t_read(nb,iunit,A,tmp,jcol)
      implicit none 
      integer nb
      integer iunit
      real*8  A(*)
      real*8  tmp(*)
      integer jcol

#include "bafdecls.fh"
#include "errquit.fh"
#include "mpif.h"
#include "D3dB.fh"

#ifdef MPI4
#include "stupid_mpi4.fh"
#endif

*     *** local variables ***
      integer MASTER,taskid
      parameter(MASTER=0)

      logical value,fillcolumn
      integer p_to, p_here,q
      integer j,k,index
      integer source,msglen,msgtype,status(2),mpierr
      integer tmp1(2),tmp2(2)

      integer taskid_j,np_j
      integer ii,jj,jstart,jend

*     **** external functions ****
      integer  Parallel2d_convert_taskid_ij,Parallel_comm_world
      external Parallel2d_convert_taskid_ij,Parallel_comm_world

      call Parallel_taskid(taskid)
      call Parallel2d_taskid_j(taskid_j)
      call Parallel2d_np_j(np_j)
      if (jcol.lt.0) then
         jstart = 0
         jend   = np_j-1
         fillcolumn = .true.
      else
         jstart = jcol
         jend   = jcol
         fillcolumn = (taskid_j.eq.jcol)
      end if


*     **** allocate status memory ****
      if (.not.BA_push_get(mt_int,MPI_STATUS_SIZE*4,
     >                     'status',status(2),status(1)))
     > call errquit('D3dB_t_read:out of stack',0,MA_ERR)


      if (mapping.eq.1) then
*        **** master node reads from file and distributes ****
         if (taskid.eq.MASTER) then
         do k=1,nz(nb)

            call dread(iunit,tmp,(nx(nb)/2+1)*ny(nb))

            call D3dB_ktoqp(nb,k,q,ii)
            
            do jj=jstart,jend

              p_to = Parallel2d_convert_taskid_ij(ii,jj)
              if (p_to.eq.MASTER) then
                 index = (q-1)*(nx(nb)/2+1)*ny(nb) + 1
                 call dcopy((nx(nb)/2+1)*ny(nb),tmp,1,A(index),1)
              else
                 msglen = (nx(nb)/2+1)*ny(nb)
                 msgtype = 129
c                 call SND(9+MSGDBL,tmp,mdtob(msglen),p_to,1)
#ifdef MPI4
                 stupid_msglen = msglen
                 stupid_type   = msgtype
                 stupid_taskid = p_to
                 call MPI_SEND(tmp,
     >                     stupid_msglen,stupid_double,
     >                     stupid_taskid,
     >                     stupid_type,stupid_world,stupid_ierr)
#else
                 call MPI_SEND(tmp,
     >                     msglen,MPI_DOUBLE_PRECISION,
     >                     p_to,
     >                     msgtype,Parallel_comm_world(),mpierr)
#endif
              end if
            end do
         end do

*        **** not master node ****
         else if (fillcolumn) then

         do k=1,nz(nb)
            call D3dB_ktoqp(nb,k,q,ii)
            p_here = Parallel2d_convert_taskid_ij(ii,taskid_j)
            if (p_here.eq.taskid) then
               msglen  = (nx(nb)/2+1)*ny(nb)
               msgtype = 129
               source  = MASTER
c               call RCV(9+MSGDBL,tmp,mdtob(msglen),rcv_len,
c     >                       source,rcv_proc,1)
#ifdef MPI4
               stupid_msglen = msglen
               stupid_type   = msgtype
               stupid_taskid = source
               call MPI_RECV(tmp,
     >                    stupid_msglen,stupid_double,
     >                    stupid_taskid,
     >                    stupid_type,stupid_world,
     >                    int_mb(status(1)),stupid_ierr)
#else
               call MPI_RECV(tmp,
     >                    msglen,MPI_DOUBLE_PRECISION,
     >                    source,
     >                    msgtype,Parallel_comm_world(),
     >                    int_mb(status(1)),mpierr)
#endif
               index = (q-1)*(nx(nb)/2+1)*ny(nb) + 1
               call dcopy((nx(nb)/2+1)*ny(nb),tmp,1,A(index),1)
               
            end if
         end do
         end if


      !*************************
      !**** hilbert mapping ****
      !*************************
      else
*     **** master node reads from file and distributes ****
      if (taskid.eq.MASTER) then
         do k=1,nz(nb)
         do j=1,ny(nb)

            call dread(iunit,tmp,(nx(nb)/2+1))

            q   = int_mb(q_map1(1,nb)+(j-1)+(k-1)*ny(nb))
            ii  = int_mb(p_map1(1,nb)+(j-1)+(k-1)*ny(nb))
            do jj=jstart,jend

               p_to = Parallel2d_convert_taskid_ij(ii,jj)
               if (p_to.eq.MASTER) then
                  index = (q-1)*(nx(nb)/2+1) + 1
                  call dcopy((nx(nb)/2+1),tmp,1,A(index),1)
               else
                  msglen  = (nx(nb)/2+1)
                  msgtype = 129
c                  call SND(9+MSGDBL,tmp,mdtob(msglen),p_to,1)
#ifdef MPI4
                  stupid_msglen = msglen
                  stupid_type   = msgtype
                  stupid_taskid = p_to
                  call MPI_SEND(tmp,
     >                     stupid_msglen,stupid_double,
     >                     stupid_taskid,
     >                     stupid_type,stupid_world,stupid_ierr)
#else
                  call MPI_SEND(tmp,
     >                     msglen,MPI_DOUBLE_PRECISION,
     >                     p_to,
     >                     msgtype,Parallel_comm_world(),mpierr)
#endif
               end if
            end do
         end do
         end do

*     **** not master node ****
      else if (fillcolumn) then
         do k=1,nz(nb)
         do j=1,ny(nb)

            q      = int_mb(q_map1(1,nb)+(j-1)+(k-1)*ny(nb))
            ii     = int_mb(p_map1(1,nb)+(j-1)+(k-1)*ny(nb))
            p_here = Parallel2d_convert_taskid_ij(ii,taskid_j)
            if (p_here.eq.taskid) then
               msglen  = (nx(nb)/2+1)
               msgtype = 129
               source  = MASTER
c               call RCV(9+MSGDBL,tmp,mdtob(msglen),rcv_len,
c     >                  source,rcv_proc,1)
#ifdef MPI4
               stupid_msglen = msglen
               stupid_type   = msgtype
               stupid_taskid = source
               call MPI_RECV(tmp,
     >                    stupid_msglen,stupid_double,
     >                    stupid_taskid,
     >                    stupid_type,stupid_world,
     >                    int_mb(status(1)),stupid_ierr)
#else
               call MPI_RECV(tmp,
     >                    msglen,MPI_DOUBLE_PRECISION,
     >                    source,
     >                    msgtype,Parallel_comm_world(),
     >                    int_mb(status(1)),mpierr)
#endif
               index = (q-1)*(nx(nb)/2+1) + 1
               call dcopy((nx(nb)/2+1),tmp,1,A(index),1)

            end if
         end do
         end do
      end if

*     **** allocate temporary space  ****
      value = BA_push_get(mt_dcpl,(nfft3d(nb)),'tmp1',
     >                    tmp1(2),tmp1(1))
      value = value.and.
     >      BA_push_get(mt_dcpl,(nfft3d(nb)),'tmp2',tmp2(2),tmp2(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)
*
      call D3dB_t_transpose_ijk(nb,5,A,
     >                          dcpl_mb(tmp1(1)),
     >                          dcpl_mb(tmp2(1)))  !*** map1to3 operation ***

*     **** deallocate temporary space  ****
      value =           BA_pop_stack(tmp2(2))
      value = value.and.BA_pop_stack(tmp1(2))
      if (.not. value) call errquit('error popping stack 1',0, MA_ERR)

      end if

*     *** deallocate status memory ***
      if (.not.BA_pop_stack(status(2)))
     > call errquit('D3dB_t_read:popping stack',0,MA_ERR)


      return
      end





*     ***********************************
*     *					*
*     *	       D3dB_(c,r,t)_write	*	
*     *					*
*     ***********************************

      subroutine D3dB_c_write_pio(nb,iunit,A,tmp,jcol)
      implicit none 
      integer nb
      integer iunit
      complex*16 A(*)
      complex*16 tmp(*)
      integer    jcol

#include "bafdecls.fh"
#include "errquit.fh"
#include "mpif.h"
#include "D3dB.fh"

#ifdef MPI4
#include "stupid_mpi4.fh"
#endif


*     *** local variables ***
      logical writecolumn
      integer MASTER,taskid,taskid_i,taskid_j
      parameter(MASTER=0)
      logical value
      integer p_from, p_here,q
      integer j,k,index
      integer dest,source,msglen,msgtype,status(2)
      integer dum,dum_msglen,dum_msgtype,mpierr
      integer tmp1(2),tmp2(2)

      integer ii

*     **** external functions ****
      integer  Parallel2d_convert_taskid_ij,Parallel2d_comm_i
      external Parallel2d_convert_taskid_ij,Parallel2d_comm_i

      call Parallel_taskid(taskid)
      call Parallel2d_taskid_i(taskid_i)
      call Parallel2d_taskid_j(taskid_j)
      writecolumn = (taskid_j.eq.jcol)

      if (writecolumn) then

*     **** allocate status memory ****
      if (.not.BA_push_get(mt_int,MPI_STATUS_SIZE*4,
     >                     'status',status(2),status(1)))
     > call errquit('D3dB_c_write:out of stack',0,MA_ERR)


      !**********************
      !**** slab mapping ****
      !**********************
      if (mapping.eq.1) then
*     **** master node gathers and writes to file ****
      if (taskid_i.eq.MASTER) then
         do k=1,nz(nb)

            call D3dB_ktoqp(nb,k,q,ii)

            p_from = ii
            !p_from = Parallel2d_convert_taskid_ij(ii,jcol)
            if (p_from.eq.MASTER) then

               index = (q-1)*(nx(nb)/2+1)*ny(nb) + 1
               call dcopy((nx(nb)+2)*ny(nb),A(index),1,tmp,1)

            else
               dum        = 99
               msglen     = (nx(nb)/2+1)*ny(nb)
               msgtype    = 139
               source     = p_from
               dum_msglen = 1
               dum_msgtype= 239
c               call SND(9+MSGINT,dum,mitob(dum_msglen),source,1)
c               call RCV(9+MSGDBL,tmp,mdtob(2*msglen),rcv_len,
c     >                       source,rcv_proc,1)
#ifdef MPI4
               stupid_msglen = dum_msglen
               stupid_type   = dum_msgtype
               stupid_taskid = source
               call MPI_SEND(dum,
     >                     stupid_msglen,stupid_integer,
     >                     stupid_taskid,
     >                     stupid_type,stupid_comm_i,stupid_ierr)
               stupid_msglen = msglen
               stupid_type   = msgtype
               call MPI_RECV(tmp,
     >                    stupid_msglen,stupid_complex,
     >                    stupid_taskid,
     >                    stupid_type,stupid_comm_i,
     >                    int_mb(status(1)),stupid_ierr)
#else
               call MPI_SEND(dum,
     >                     dum_msglen,MPI_INTEGER,
     >                     source,
     >                     dum_msgtype,Parallel2d_comm_i(),mpierr)
               call MPI_RECV(tmp,
     >                    msglen,MPI_DOUBLE_COMPLEX,
     >                    source,
     >                    msgtype,Parallel2d_comm_i(),
     >                    int_mb(status(1)),mpierr)
#endif
            end if

            call dwrite(iunit,tmp,(nx(nb)+2)*ny(nb)) 
         end do

*     **** not master node ****
      else
         do k=1,nz(nb)
         
            call D3dB_ktoqp(nb,k,q,ii)
            p_here = ii
            !p_here = Parallel2d_convert_taskid_ij(ii,jcol)
            if (p_here.eq.taskid_i) then

               index = (q-1)*(nx(nb)/2+1)*ny(nb) + 1
               call dcopy((nx(nb)+2)*ny(nb),A(index),1,tmp,1)

               msglen     = (nx(nb)/2+1)*ny(nb)
               msgtype    = 139
               dest       = MASTER
               dum_msglen = 1
               dum_msgtype= 239
c               call RCV(9+MSGINT,dum,mitob(dum_msglen),rcv_len,
c     >                       dest,rcv_proc,1)
c               call SND(9+MSGDBL,tmp,mdtob(2*msglen),dest,1)
#ifdef MPI4
                  stupid_msglen = dum_msglen
                  stupid_type   = dum_msgtype
                  stupid_taskid = dest
                  call MPI_RECV(dum,
     >                    stupid_msglen,stupid_integer,
     >                    stupid_taskid,
     >                    stupid_type,stupid_comm_i,
     >                    int_mb(status(1)),stupid_ierr)
                  stupid_msglen = msglen
                  stupid_type   = msgtype
                  call MPI_SEND(tmp,
     >                     stupid_msglen,stupid_complex,
     >                     stupid_taskid,
     >                     stupid_type,stupid_comm_i,stupid_ierr)
#else
                  call MPI_RECV(dum,
     >                    dum_msglen,MPI_INTEGER,
     >                    dest,
     >                    dum_msgtype,Parallel2d_comm_i(),
     >                    int_mb(status(1)),mpierr)
                  call MPI_SEND(tmp,
     >                     msglen,MPI_DOUBLE_COMPLEX,
     >                     dest,
     >                     msgtype,Parallel2d_comm_i(),mpierr)
#endif
            end if

         end do
      end if


      !*************************
      !**** hilbert mapping ****
      !*************************
      else

      if (taskid_j.eq.jcol) then

*       **** allocate temporary space  ****
        value = BA_push_get(mt_dcpl,(nfft3d(nb)),'tmp1',
     >                      tmp1(2),tmp1(1))
        value = value.and.
     >          BA_push_get(mt_dcpl,(nfft3d(nb)),'tmp2',tmp2(2),tmp2(1))
        if (.not. value) call errquit('out of stack memory',0, MA_ERR)
*
        call D3dB_c_transpose_ijk(nb,6,A,
     >                            dcpl_mb(tmp1(1)),
     >                            dcpl_mb(tmp2(1)))  !*** map3to1 operation ***

*       **** deallocate temporary space  ****
        value = BA_pop_stack(tmp2(2))
        value = value.and.BA_pop_stack(tmp1(2))
        if (.not. value) call errquit('error popping stack 3',0, MA_ERR)

      end if

c      call ga_sync()

*     **** master node reads from file and distributes ****
      if (taskid_i.eq.MASTER) then
         do k=1,nz(nb)
         do j=1,ny(nb)

            q      = int_mb(q_map1(1,nb)+(j-1)+(k-1)*ny(nb))
            ii     = int_mb(p_map1(1,nb)+(j-1)+(k-1)*ny(nb))
            p_from = ii
            !p_from = Parallel2d_convert_taskid_ij(ii,jcol)
            if (p_from.eq.MASTER) then

               index = (q-1)*(nx(nb)/2+1) + 1
               call dcopy((nx(nb)+2),A(index),1,tmp,1)

            else
               dum         = 99
               msglen      = (nx(nb)/2+1)
               msgtype     = 139
               source      = p_from
               dum_msglen  = 1
               dum_msgtype = 239
c               call SND(9+MSGINT,dum,mitob(dum_msglen),source,1)
c               call RCV(9+MSGDBL,tmp,mdtob(2*msglen),rcv_len,
c     >                       source,rcv_proc,1)
#ifdef MPI4
               stupid_msglen = dum_msglen
               stupid_type   = dum_msgtype
               stupid_taskid = source
               call MPI_SEND(tmp,
     >                     stupid_msglen,stupid_integer,
     >                     stupid_taskid,
     >                     stupid_type,stupid_comm_i,stupid_ierr)
               stupid_msglen = msglen
               stupid_type   = msgtype
               call MPI_RECV(tmp,
     >                    stupid_msglen,stupid_complex,
     >                    stupid_taskid,
     >                    stupid_type,stupid_comm_i,
     >                    int_mb(status(1)),stupid_ierr)
#else
               call MPI_SEND(dum,
     >                     dum_msglen,MPI_INTEGER,
     >                     source,
     >                     dum_msgtype,Parallel2d_comm_i(),mpierr)
               call MPI_RECV(tmp,
     >                    msglen,MPI_DOUBLE_COMPLEX,
     >                    source,
     >                    msgtype,Parallel2d_comm_i(),
     >                    int_mb(status(1)),mpierr)
#endif
            end if
            call dwrite(iunit,tmp,(nx(nb)+2))
         end do
         end do

*     **** not master node ****
      else
         do k=1,nz(nb)
         do j=1,ny(nb)

            q      = int_mb(q_map1(1,nb)+(j-1)+(k-1)*ny(nb))
            ii     = int_mb(p_map1(1,nb)+(j-1)+(k-1)*ny(nb))
            p_here = ii
            !p_here = Parallel2d_convert_taskid_ij(ii,jcol) 
            if (p_here.eq.taskid_i) then

               index = (q-1)*(nx(nb)/2+1) + 1
               call dcopy((nx(nb)+2),A(index),1,tmp,1)

               msglen     = (nx(nb)/2+1)
               msgtype    = 139
               dest       = MASTER
               dum_msglen = 1
               dum_msgtype= 239
c               call RCV(9+MSGINT,dum,mitob(dum_msglen),rcv_len,
c     >                       dest,rcv_proc,1)
c               call SND(9+MSGDBL,tmp,mdtob(2*msglen),dest,1)
#ifdef MPI4
                  stupid_msglen = dum_msglen
                  stupid_type   = dum_msgtype
                  stupid_taskid = dest
                  call MPI_RECV(dum,
     >                    stupid_msglen,stupid_integer,
     >                    stupid_taskid,
     >                    stupid_type,stupid_comm_i,
     >                    int_mb(status(1)),stupid_ierr)
                  stupid_msglen = msglen
                  stupid_type   = msgtype
                  call MPI_SEND(tmp,
     >                     stupid_msglen,stupid_complex,
     >                     stupid_taskid,
     >                     stupid_type,stupid_comm_i,stupid_ierr)
#else
                  call MPI_RECV(dum,
     >                    dum_msglen,MPI_INTEGER,
     >                    dest,
     >                    dum_msgtype,Parallel2d_comm_i(),
     >                    int_mb(status(1)),mpierr)
                  call MPI_SEND(tmp,
     >                     msglen,MPI_DOUBLE_COMPLEX,
     >                     dest,
     >                     msgtype,Parallel2d_comm_i(),mpierr)
#endif
            end if
         end do
         end do
      end if

      end if

*     *** deallocate status memory ***
      if (.not.BA_pop_stack(status(2)))
     > call errquit('D3dB_c_write:popping stack',0,MA_ERR)


      !*** shift filepointer by (nx(nb)+2)*ny(nb)*nz(nb) doubles ****
      else
         if (taskid_i.eq.MASTER) 
     >      call dshift_fileptr(iunit,(nx(nb)+2)*ny(nb)*nz(nb))
      end if !*** writecolumn ***

      return
      end


      subroutine D3dB_c_write(nb,iunit,A,tmp,jcol)
      implicit none 
      integer nb
      integer iunit
      complex*16 A(*)
      complex*16 tmp(*)
      integer    jcol

#include "bafdecls.fh"
#include "errquit.fh"
#include "mpif.h"
#include "D3dB.fh"

#ifdef MPI4
#include "stupid_mpi4.fh"
#endif


*     *** local variables ***
      integer MASTER,taskid,taskid_j
      parameter(MASTER=0)
      logical value
      integer p_from, p_here,q
      integer j,k,index
      integer dest,source,msglen,msgtype,status(2),mpierr
      integer dum,dum_msglen,dum_msgtype
      integer tmp1(2),tmp2(2)

      integer ii

*     **** external functions ****
      integer  Parallel2d_convert_taskid_ij,Parallel_comm_world
      external Parallel2d_convert_taskid_ij,Parallel_comm_world

      call Parallel_taskid(taskid)
      call Parallel2d_taskid_j(taskid_j)


*     **** allocate status memory ****
      if (.not.BA_push_get(mt_int,MPI_STATUS_SIZE*4,
     >                     'status',status(2),status(1)))
     > call errquit('D3dB_c_write:out of stack',0,MA_ERR)


      !**********************
      !**** slab mapping ****
      !**********************
      if (mapping.eq.1) then
*     **** master node gathers and writes to file ****
      if (taskid.eq.MASTER) then
         do k=1,nz(nb)

            call D3dB_ktoqp(nb,k,q,ii)

            p_from = Parallel2d_convert_taskid_ij(ii,jcol)
            if (p_from.eq.MASTER) then

               index = (q-1)*(nx(nb)/2+1)*ny(nb) + 1
               call dcopy((nx(nb)+2)*ny(nb),A(index),1,tmp,1)

            else
               dum        = 99
               msglen     = (nx(nb)/2+1)*ny(nb)
               msgtype    = 149
               source     = p_from
               dum_msglen = 1
               dum_msgtype= 249
c               call SND(9+MSGINT,dum,mitob(dum_msglen),source,1)
c               call RCV(9+MSGDBL,tmp,mdtob(2*msglen),rcv_len,
c     >                       source,rcv_proc,1)
#ifdef MPI4
               stupid_msglen = dum_msglen
               stupid_type   = dum_msgtype
               stupid_taskid = source
               call MPI_SEND(dum,
     >                     stupid_msglen,stupid_integer,
     >                     stupid_taskid,
     >                     stupid_type,stupid_world,stupid_ierr)
               stupid_msglen = msglen
               stupid_type   = msgtype
               call MPI_RECV(tmp,
     >                    stupid_msglen,stupid_complex,
     >                    stupid_taskid,
     >                    stupid_type,stupid_world,
     >                    int_mb(status(1)),stupid_ierr)
#else
               call MPI_SEND(dum,
     >                     dum_msglen,MPI_INTEGER,
     >                     source,
     >                     dum_msgtype,Parallel_comm_world(),mpierr)
               call MPI_RECV(tmp,
     >                    msglen,MPI_DOUBLE_COMPLEX,
     >                    source,
     >                    msgtype,Parallel_comm_world(),
     >                    int_mb(status(1)),mpierr)
#endif
            end if

            call dwrite(iunit,tmp,(nx(nb)+2)*ny(nb)) 
         end do

*     **** not master node ****
      else
         do k=1,nz(nb)
         
            call D3dB_ktoqp(nb,k,q,ii)
            p_here = Parallel2d_convert_taskid_ij(ii,jcol)
            if (p_here.eq.taskid) then

               index = (q-1)*(nx(nb)/2+1)*ny(nb) + 1
               call dcopy((nx(nb)+2)*ny(nb),A(index),1,tmp,1)

               msglen     = (nx(nb)/2+1)*ny(nb)
               msgtype    = 149
               dest       = MASTER
               dum_msglen = 1
               dum_msgtype= 249
c               call RCV(9+MSGINT,dum,mitob(dum_msglen),rcv_len,
c     >                       dest,rcv_proc,1)
c               call SND(9+MSGDBL,tmp,mdtob(2*msglen),dest,1)
#ifdef MPI4
                  stupid_msglen = dum_msglen
                  stupid_type   = dum_msgtype
                  stupid_taskid = dest
                  call MPI_RECV(dum,
     >                    stupid_msglen,stupid_integer,
     >                    stupid_taskid,
     >                    stupid_type,stupid_world,
     >                    int_mb(status(1)),stupid_ierr)
                  stupid_msglen = msglen
                  stupid_type   = msgtype
                  call MPI_SEND(tmp,
     >                     stupid_msglen,stupid_complex,
     >                     stupid_taskid,
     >                     stupid_type,stupid_world,stupid_ierr)
#else
                  call MPI_RECV(dum,
     >                    dum_msglen,MPI_INTEGER,
     >                    dest,
     >                    dum_msgtype,Parallel_comm_world(),
     >                    int_mb(status(1)),mpierr)
                  call MPI_SEND(tmp,
     >                     msglen,MPI_DOUBLE_COMPLEX,
     >                     dest,
     >                     msgtype,Parallel_comm_world(),mpierr)
#endif
            end if

         end do
      end if


      !*************************
      !**** hilbert mapping ****
      !*************************
      else

      if (taskid_j.eq.jcol) then

*       **** allocate temporary space  ****
        value = BA_push_get(mt_dcpl,(nfft3d(nb)),'tmp1',
     >                      tmp1(2),tmp1(1))
        value = value.and.
     >          BA_push_get(mt_dcpl,(nfft3d(nb)),'tmp2',tmp2(2),tmp2(1))
        if (.not. value) call errquit('out of stack memory',0, MA_ERR)
*
        call D3dB_c_transpose_ijk(nb,6,A,
     >                            dcpl_mb(tmp1(1)),
     >                            dcpl_mb(tmp2(1)))  !*** map3to1 operation ***

*       **** deallocate temporary space  ****
        value = BA_pop_stack(tmp2(2))
        value = value.and.BA_pop_stack(tmp1(2))
        if (.not. value) call errquit('error popping stack 4',0, MA_ERR)

      end if

c      call ga_sync()

*     **** master node reads from file and distributes ****
      if (taskid.eq.MASTER) then
         do k=1,nz(nb)
         do j=1,ny(nb)

            q      = int_mb(q_map1(1,nb)+(j-1)+(k-1)*ny(nb))
            ii     = int_mb(p_map1(1,nb)+(j-1)+(k-1)*ny(nb))
            p_from = Parallel2d_convert_taskid_ij(ii,jcol)
            if (p_from.eq.MASTER) then

               index = (q-1)*(nx(nb)/2+1) + 1
               call dcopy((nx(nb)+2),A(index),1,tmp,1)

            else
               dum         = 99
               msglen      = (nx(nb)/2+1)
               msgtype     = 149
               source      = p_from
               dum_msglen  = 1
               dum_msgtype = 249
c               call SND(9+MSGINT,dum,mitob(dum_msglen),source,1)
c               call RCV(9+MSGDBL,tmp,mdtob(2*msglen),rcv_len,
c     >                       source,rcv_proc,1)
#ifdef MPI4
               stupid_msglen = dum_msglen
               stupid_type   = dum_msgtype
               stupid_taskid = source
               call MPI_SEND(tmp,
     >                     stupid_msglen,stupid_integer,
     >                     stupid_taskid,
     >                     stupid_type,stupid_world,stupid_ierr)
               stupid_msglen = msglen
               stupid_type   = msgtype
               call MPI_RECV(tmp,
     >                    stupid_msglen,stupid_complex,
     >                    stupid_taskid,
     >                    stupid_type,stupid_world,
     >                    int_mb(status(1)),stupid_ierr)
#else
               call MPI_SEND(dum,
     >                     dum_msglen,MPI_INTEGER,
     >                     source,
     >                     dum_msgtype,Parallel_comm_world(),mpierr)
               call MPI_RECV(tmp,
     >                    msglen,MPI_DOUBLE_COMPLEX,
     >                    source,
     >                    msgtype,Parallel_comm_world(),
     >                    int_mb(status(1)),mpierr)
#endif
            end if
            call dwrite(iunit,tmp,(nx(nb)+2))
         end do
         end do

*     **** not master node ****
      else
         do k=1,nz(nb)
         do j=1,ny(nb)

            q      = int_mb(q_map1(1,nb)+(j-1)+(k-1)*ny(nb))
            ii     = int_mb(p_map1(1,nb)+(j-1)+(k-1)*ny(nb))
            p_here = Parallel2d_convert_taskid_ij(ii,jcol) 
            if (p_here.eq.taskid) then

               index = (q-1)*(nx(nb)/2+1) + 1
               call dcopy((nx(nb)+2),A(index),1,tmp,1)

               msglen     = (nx(nb)/2+1)
               msgtype    = 149
               dest       = MASTER
               dum_msglen = 1
               dum_msgtype= 249
c               call RCV(9+MSGINT,dum,mitob(dum_msglen),rcv_len,
c     >                       dest,rcv_proc,1)
c               call SND(9+MSGDBL,tmp,mdtob(2*msglen),dest,1)
#ifdef MPI4
                  stupid_msglen = dum_msglen
                  stupid_type   = dum_msgtype
                  stupid_taskid = dest
                  call MPI_RECV(dum,
     >                    stupid_msglen,stupid_integer,
     >                    stupid_taskid,
     >                    stupid_type,stupid_world,
     >                    int_mb(status(1)),stupid_ierr)
                  stupid_msglen = msglen
                  stupid_type   = msgtype
                  call MPI_SEND(tmp,
     >                     stupid_msglen,stupid_complex,
     >                     stupid_taskid,
     >                     stupid_type,stupid_world,stupid_ierr)
#else
                  call MPI_RECV(dum,
     >                    dum_msglen,MPI_INTEGER,
     >                    dest,
     >                    dum_msgtype,Parallel_comm_world(),
     >                    int_mb(status(1)),mpierr)
                  call MPI_SEND(tmp,
     >                     msglen,MPI_DOUBLE_COMPLEX,
     >                     dest,
     >                     msgtype,Parallel_comm_world(),mpierr)
#endif
            end if
         end do
         end do
      end if

      end if

*     *** deallocate status memory ***
      if (.not.BA_pop_stack(status(2)))
     > call errquit('D3dB_c_write:popping stack',0,MA_ERR)

      return
      end





      subroutine D3dB_r_write(nb,iunit,A,tmp,jcol)
      implicit none 
      integer nb
      integer iunit
      real*8     A(*)
      real*8     tmp(*)
      integer    jcol

#include "bafdecls.fh"
#include "errquit.fh"
#include "mpif.h"
#include "D3dB.fh"

#ifdef MPI4
#include "stupid_mpi4.fh"
#endif

*     *** local variables ***
      integer MASTER,taskid
      parameter(MASTER=0)
      integer p_from, p_here,q
      integer j,k,index,index2
      integer dest,source,msglen,msgtype,status(2),mpierr

      integer taskid_j,ii

*     **** external functions ****
      integer  Parallel2d_convert_taskid_ij,Parallel_comm_world
      external Parallel2d_convert_taskid_ij,Parallel_comm_world

      call Parallel_taskid(taskid)
      call Parallel2d_taskid_j(taskid_j)

*     **** allocate status memory ****
      if (.not.BA_push_get(mt_int,MPI_STATUS_SIZE*4,
     >                     'status',status(2),status(1)))
     > call errquit('D3dB_r_write:out of stack',0,MA_ERR)


      !**********************
      !**** slab mapping ****
      !**********************
      if (mapping.eq.1) then
*     **** master node reads from file and distributes ****
      if (taskid.eq.MASTER) then
         do k=1,nz(nb)

            call D3dB_ktoqp(nb,k,q,ii)
            p_from = Parallel2d_convert_taskid_ij(ii,jcol)
            if (p_from.eq.MASTER) then

               do j=1,ny(nb)
                 index = (q-1)*(nx(nb)+2)*ny(nb) 
     >                 + (j-1)*(nx(nb)+2) + 1
                 index2 = (j-1)*nx(nb) + 1
                 call dcopy(nx(nb),A(index),1,tmp(index2),1)
               end do

            else
               msglen  = (nx(nb))*ny(nb)
               msgtype = 159
               source  = p_from
c               call RCV(9+MSGDBL,tmp,mdtob(msglen),rcv_len,
c     >                       source,rcv_proc,1)
#ifdef MPI4
               stupid_msglen = msglen
               stupid_type   = msgtype
               stupid_taskid = source
               call MPI_RECV(tmp,
     >                    stupid_msglen,stupid_double,
     >                    stupid_taskid,
     >                    stupid_type,stupid_world,
     >                    int_mb(status(1)),stupid_ierr)
#else
               call MPI_RECV(tmp,
     >                    msglen,MPI_DOUBLE_PRECISION,
     >                    source,
     >                    msgtype,Parallel_comm_world(),
     >                    int_mb(status(1)),mpierr)
#endif
            end if

            call dwrite(iunit,tmp,(nx(nb))*ny(nb)) 
         end do

*     **** not master node ****
      else
         do k=1,nz(nb)
         
            call D3dB_ktoqp(nb,k,q,ii)
            p_here = Parallel2d_convert_taskid_ij(ii,jcol)
            if (p_here.eq.taskid) then

               do j=1,ny(nb)
                  index = (q-1)*(nx(nb)+2)*ny(nb) 
     >                  + (j-1)*(nx(nb)+2) + 1
                  index2 = (j-1)*nx(nb) + 1
                  call dcopy(nx(nb),A(index),1,tmp(index2),1)
               end do

               msglen  = (nx(nb))*ny(nb)
               msgtype = 159
               dest    = MASTER
c               call SND(9+MSGDBL,tmp,mdtob(msglen),dest,1)
#ifdef MPI4
                  stupid_msglen = msglen
                  stupid_type   = msgtype
                  stupid_taskid = dest
                  call MPI_SEND(tmp,
     >                     stupid_msglen,stupid_double,
     >                     stupid_taskid,
     >                     stupid_type,stupid_world,stupid_ierr)
#else
                  call MPI_SEND(tmp,
     >                     msglen,MPI_DOUBLE_PRECISION,
     >                     dest,
     >                     msgtype,Parallel_comm_world(),mpierr)
#endif
            end if

         end do
      end if

      !*************************
      !**** hilbert mapping ****
      !*************************
      else
*     **** master node reads from file and distributes ****
      if (taskid.eq.MASTER) then
         do k=1,nz(nb)
         do j=1,ny(nb)

            q      = int_mb(q_map1(1,nb)+(j-1)+(k-1)*ny(nb))
            ii     = int_mb(p_map1(1,nb)+(j-1)+(k-1)*ny(nb))
            p_from = Parallel2d_convert_taskid_ij(ii,jcol)      

            if (p_from.eq.MASTER) then
              index = (q-1)*(nx(nb)+2) + 1
              call dcopy(nx(nb),A(index),1,tmp,1)

            else
               msglen  = (nx(nb))
               msgtype = 159
               source  = p_from
c               call RCV(9+MSGDBL,tmp,mdtob(msglen),rcv_len,
c     >                       source,rcv_proc,1)
#ifdef MPI4
               stupid_msglen = msglen
               stupid_type   = msgtype
               stupid_taskid = source
               call MPI_RECV(tmp,
     >                    stupid_msglen,stupid_double,
     >                    stupid_taskid,
     >                    stupid_type,stupid_world,
     >                    int_mb(status(1)),stupid_ierr)
#else
               call MPI_RECV(tmp,
     >                    msglen,MPI_DOUBLE_PRECISION,
     >                    source,
     >                    msgtype,Parallel_comm_world(),
     >                    int_mb(status(1)),mpierr)
#endif


            end if

            call dwrite(iunit,tmp,(nx(nb)))
         end do
         end do

*     **** not master node ****
      else
         do k=1,nz(nb)
         do j=1,ny(nb)

            q      = int_mb(q_map1(1,nb)+(j-1)+(k-1)*ny(nb))
            ii     = int_mb(p_map1(1,nb)+(j-1)+(k-1)*ny(nb))
            p_here = Parallel2d_convert_taskid_ij(ii,jcol)
            if (p_here.eq.taskid) then

               index = (q-1)*(nx(nb)+2) + 1
               call dcopy(nx(nb),A(index),1,tmp,1)

               msglen  = nx(nb)
               msgtype = 159
               dest    = MASTER

c               call SND(9+MSGDBL,tmp,mdtob(msglen),dest,1)
#ifdef MPI4
                  stupid_msglen = msglen
                  stupid_type   = msgtype
                  stupid_taskid = dest
                  call MPI_SEND(tmp,
     >                     stupid_msglen,stupid_double,
     >                     stupid_taskid,
     >                     stupid_type,stupid_world,stupid_ierr)
#else
                  call MPI_SEND(tmp,
     >                     msglen,MPI_DOUBLE_PRECISION,
     >                     dest,
     >                     msgtype,Parallel_comm_world(),mpierr)
#endif
            end if

         end do
         end do
      end if

      end if

*     *** deallocate status memory ***
      if (.not.BA_pop_stack(status(2)))
     > call errquit('D3dB_r_write:popping stack',0,MA_ERR)

      return
      end

      subroutine D3dB_t_write(nb,iunit,A,tmp,jcol)
      implicit none 
      integer nb
      integer iunit
      real*8  A(*)
      real*8  tmp(*)
      integer jcol

#include "bafdecls.fh"
#include "errquit.fh"
#include "mpif.h"
#include "D3dB.fh"

#ifdef MPI4
#include "stupid_mpi4.fh"
#endif


*     *** local variables ***
      integer MASTER,taskid
      parameter(MASTER=0)
      logical value
      integer p_from, p_here,q
      integer j,k,index
      integer dest,source,msglen,msgtype,status(2),mpierr
      integer dum,dum_msglen,dum_msgtype
      integer tmp1(2),tmp2(2)

      integer ii,taskid_j

*     **** external functions ****
      integer  Parallel2d_convert_taskid_ij,Parallel_comm_world
      external Parallel2d_convert_taskid_ij,Parallel_comm_world

      call Parallel_taskid(taskid)
      call Parallel2d_taskid_j(taskid_j)

*     **** allocate status memory ****
      if (.not.BA_push_get(mt_int,MPI_STATUS_SIZE*4,
     >                     'status',status(2),status(1)))
     > call errquit('D3dB_t_write:out of stack',0,MA_ERR)


      !**********************
      !**** slab mapping ****
      !**********************
      if (mapping.eq.1) then
*     **** master node reads from file and distributes ****
      if (taskid.eq.MASTER) then
         do k=1,nz(nb)

            call D3dB_ktoqp(nb,k,q,ii)

            p_from = Parallel2d_convert_taskid_ij(ii,jcol)
            if (p_from.eq.MASTER) then

               index = (q-1)*(nx(nb)/2+1)*ny(nb) + 1
               call dcopy((nx(nb)/2+1)*ny(nb),A(index),1,tmp,1)

            else
               dum        = 99
               msglen     = (nx(nb)/2+1)*ny(nb)
               msgtype    = 169
               source     = p_from
               dum_msglen = 1
               dum_msgtype= 269
c               call SND(9+MSGINT,dum,mitob(dum_msglen),source,1)
c               call RCV(9+MSGDBL,tmp,mdtob(msglen),rcv_len,
c     >                       source,rcv_proc,1)
#ifdef MPI4
               stupid_msglen = dum_msglen
               stupid_type   = dum_msgtype
               stupid_taskid = source
                call MPI_SEND(dum,
     >                     stupid_msglen,stupid_integer,
     >                     stupid_taskid,
     >                     stupid_type,stupid_world,stupid_ierr)
               stupid_msglen = msglen
               stupid_type   = msgtype
               call MPI_RECV(tmp,
     >                    stupid_msglen,stupid_double,
     >                    stupid_taskid,
     >                    stupid_type,stupid_world,
     >                    int_mb(status(1)),stupid_ierr)
#else
               call MPI_SEND(dum,
     >                     dum_msglen,MPI_INTEGER,
     >                     source,
     >                     dum_msgtype,Parallel_comm_world(),mpierr)
               call MPI_RECV(tmp,
     >                    msglen,MPI_DOUBLE_PRECISION,
     >                    source,
     >                    msgtype,Parallel_comm_world(),
     >                    int_mb(status(1)),mpierr)
#endif


            end if

            call dwrite(iunit,tmp,(nx(nb)/2+1)*ny(nb)) 
         end do

*     **** not master node ****
      else
         do k=1,nz(nb)
         
            call D3dB_ktoqp(nb,k,q,ii)
            p_here = Parallel2d_convert_taskid_ij(ii,jcol)
            if (p_here.eq.taskid) then

               index = (q-1)*(nx(nb)/2+1)*ny(nb) + 1
               call dcopy((nx(nb)/2+1)*ny(nb),A(index),1,tmp,1)

               msglen     = (nx(nb)/2+1)*ny(nb)
               msgtype    = 169
               dest       = MASTER
               dum_msglen = 1
               dum_msgtype= 269
c               call RCV(9+MSGINT,dum,mitob(dum_msglen),rcv_len,
c     >                       dest,rcv_proc,1)
c               call SND(9+MSGDBL,tmp,mdtob(msglen),dest,1)
#ifdef MPI4
                  stupid_msglen = dum_msglen
                  stupid_type   = dum_msgtype
                  stupid_taskid = dest
                  call MPI_RECV(dum,
     >                    stupid_msglen,stupid_integer,
     >                    stupid_taskid,
     >                    stupid_type,stupid_world,
     >                    int_mb(status(1)),stupid_ierr)
                  stupid_msglen = msglen
                  stupid_type   = msgtype
                  call MPI_SEND(tmp,
     >                     stupid_msglen,stupid_double,
     >                     stupid_taskid,
     >                     stupid_type,stupid_world,stupid_ierr)
#else
                  call MPI_RECV(dum,
     >                    dum_msglen,MPI_INTEGER,
     >                    dest,
     >                    dum_msgtype,Parallel_comm_world(),
     >                    int_mb(status(1)),mpierr)
                  call MPI_SEND(tmp,
     >                     msglen,MPI_DOUBLE_PRECISION,
     >                     dest,
     >                     msgtype,Parallel_comm_world(),mpierr)
#endif
            end if

         end do
      end if


      !*************************
      !**** hilbert mapping ****
      !*************************
      else

      if (taskid_j.eq.jcol) then

*       **** allocate temporary space  ****
        value = BA_push_get(mt_dcpl,(nfft3d(nb)),'tmp1',
     >                      tmp1(2),tmp1(1))
        value = value.and.
     >          BA_push_get(mt_dcpl,(nfft3d(nb)),'tmp2',tmp2(2),tmp2(1))
        if (.not. value) call errquit('out of stack memory',0, MA_ERR)
*
        call D3dB_t_transpose_ijk(nb,6,A,
     >                            dcpl_mb(tmp1(1)),
     >                            dcpl_mb(tmp2(1)))  !*** map3to1 operation ***

*       **** deallocate temporary space  ****
        value = BA_pop_stack(tmp2(2))
        value = value.and.BA_pop_stack(tmp1(2))
        if (.not. value) call errquit('error popping stack 5',0, MA_ERR)

      end if


*     **** master node reads from file and distributes ****
      if (taskid.eq.MASTER) then
         do k=1,nz(nb)
         do j=1,ny(nb)

            q      = int_mb(q_map1(1,nb)+(j-1)+(k-1)*ny(nb))
            ii     = int_mb(p_map1(1,nb)+(j-1)+(k-1)*ny(nb))
            p_from = Parallel2d_convert_taskid_ij(ii,jcol)
            if (p_from.eq.MASTER) then

               index = (q-1)*(nx(nb)/2+1) + 1
               call dcopy((nx(nb)/2+1),A(index),1,tmp,1)

            else

               dum         = 99
               msglen      = (nx(nb)/2+1)
               msgtype     = 169
               source      = p_from
               dum_msglen  = 1
               dum_msgtype = 269
c               call SND(9+MSGINT,dum,mitob(dum_msglen),source,1)
c               call RCV(9+MSGDBL,tmp,mdtob(msglen),rcv_len,
c     >                       source,rcv_proc,1)
#ifdef MPI4
               stupid_msglen = dum_msglen
               stupid_type   = dum_msgtype
               stupid_taskid = source
               call MPI_SEND(dum,
     >                     stupid_msglen,stupid_integer,
     >                     stupid_taskid,
     >                     stupid_type,stupid_world,stupid_ierr)
               stupid_msglen = msglen
               stupid_type   = msgtype
               call MPI_RECV(tmp,
     >                    stupid_msglen,stupid_double,
     >                    stupid_taskid,
     >                    stupid_type,stupid_world,
     >                    int_mb(status(1)),stupid_ierr)
#else
               call MPI_SEND(dum,
     >                     dum_msglen,MPI_INTEGER,
     >                     source,
     >                     dum_msgtype,Parallel_comm_world(),mpierr)
               call MPI_RECV(tmp,
     >                    msglen,MPI_DOUBLE_PRECISION,
     >                    source,
     >                    msgtype,Parallel_comm_world(),
     >                    int_mb(status(1)),mpierr)
#endif


            end if

            call dwrite(iunit,tmp,(nx(nb)/2+1))
         end do
         end do

*     **** not master node ****
      else
         do k=1,nz(nb)
         do j=1,ny(nb)

            q      = int_mb(q_map1(1,nb)+(j-1)+(k-1)*ny(nb))
            ii     = int_mb(p_map1(1,nb)+(j-1)+(k-1)*ny(nb))
            p_here = Parallel2d_convert_taskid_ij(ii,jcol) 
            if (p_here.eq.taskid) then

               index = (q-1)*(nx(nb)/2+1) + 1
               call dcopy((nx(nb)/2+1),A(index),1,tmp,1)

               msglen     = (nx(nb)/2+1)
               msgtype    = 169
               dest       = MASTER
               dum_msglen = 1
               dum_msgtype= 269
c               call RCV(9+MSGINT,dum,mitob(dum_msglen),rcv_len,
c     >                       dest,rcv_proc,1)
c               call SND(9+MSGDBL,tmp,mdtob(msglen),dest,1)
#ifdef MPI4
                  stupid_msglen = dum_msglen
                  stupid_type   = dum_msgtype
                  stupid_taskid = dest
                  call MPI_RECV(dum,
     >                    stupid_msglen,stupid_integer,
     >                    stupid_taskid,
     >                    stupid_type,stupid_world,
     >                    int_mb(status(1)),stupid_ierr)
                  stupid_msglen = msglen
                  stupid_type   = msgtype
                  call MPI_SEND(tmp,
     >                     stupid_msglen,stupid_double,
     >                     stupid_taskid,
     >                     stupid_type,stupid_world,stupid_ierr)
#else
               call MPI_RECV(dum,
     >                    dum_msglen,MPI_INTEGER,
     >                    dest,
     >                    dum_msgtype,Parallel_comm_world(),
     >                    int_mb(status(1)),mpierr)
                  call MPI_SEND(tmp,
     >                     msglen,MPI_DOUBLE_PRECISION,
     >                     dest,
     >                     msgtype,Parallel_comm_world(),mpierr)
#endif
            end if

         end do
         end do
      end if

      end if

*     *** deallocate status memory ***
      if (.not.BA_pop_stack(status(2)))
     > call errquit('D3dB_t_write:popping stack',0,MA_ERR)

      return
      end



*     ***********************************
*     *					*
*     *	       D3dB_r_FormatWrite	*	
*     *					*
*     ***********************************

      subroutine D3dB_r_FormatWrite(nb,iunit,A,tmp)
      implicit none 
      integer nb
      integer iunit
      real*8     A(*)
      real*8     tmp(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "mpif.h"
#include "D3dB.fh"

#ifdef MPI4
#include "stupid_mpi4.fh"
#endif

*     *** local variables ***
      integer MASTER,taskid
      parameter(MASTER=0)
      integer p_from, p_here,q
      integer i,j,k,index,index2
      integer dest,source,msglen,msgtype,status(2),mpierr
      integer  Parallel_comm_world
      external Parallel_comm_world

      call Parallel_taskid(taskid)

*     **** allocate status memory ****
      if (.not.BA_push_get(mt_int,MPI_STATUS_SIZE*4,
     >                     'status',status(2),status(1)))
     > call errquit('D3dB_r_FormatWrite:out of stack',0,MA_ERR)


      !**********************
      !**** slab mapping ****
      !**********************
      if (mapping.eq.1) then
*     **** master node reads from file and distributes ****
      if (taskid.eq.MASTER) then
         do k=1,nz(nb)

            call D3dB_ktoqp(nb,k,q,p_from)

            if (p_from.eq.MASTER) then

               do j=1,ny(nb)
                 index = (q-1)*(nx(nb)+2)*ny(nb) 
     >                 + (j-1)*(nx(nb)+2) + 1
                 index2 = (j-1)*nx(nb) + 1
                 call dcopy(nx(nb),A(index),1,tmp(index2),1)
               end do

            else
               msglen  = (nx(nb))*ny(nb)
               msgtype = 179
               source  = p_from

c               call RCV(9+MSGDBL,tmp,mdtob(msglen),rcv_len,
c     >                       source,rcv_proc,1)
#ifdef MPI4
               stupid_msglen = msglen
               stupid_type   = msgtype
               stupid_taskid = source
               call MPI_RECV(tmp,
     >                    stupid_msglen,stupid_double,
     >                    stupid_taskid,
     >                    stupid_type,stupid_world,
     >                    int_mb(status(1)),stupid_ierr)
#else
               call MPI_RECV(tmp,
     >                    msglen,MPI_DOUBLE_PRECISION,
     >                    source,
     >                    msgtype,Parallel_comm_world(),
     >                    int_mb(status(1)),mpierr)
#endif
            end if

c           **** call dwrite(iunit,tmp,(nx(nb))*ny(nb)) ****
            do j=1,ny(nb)
              write(iunit,'(3E26.14)') (tmp(i+(j-1)*nx(nb)), i=1,nx(nb))
            end do

         end do

*     **** not master node ****
      else
         do k=1,nz(nb)
         
            call D3dB_ktoqp(nb,k,q,p_here)
            if (p_here.eq.taskid) then

               do j=1,ny(nb)
                  index = (q-1)*(nx(nb)+2)*ny(nb) 
     >                  + (j-1)*(nx(nb)+2) + 1
                  index2 = (j-1)*nx(nb) + 1
                  call dcopy(nx(nb),A(index),1,tmp(index2),1)
               end do

               msglen  = (nx(nb))*ny(nb)
               msgtype = 179
               dest    = MASTER
c               call SND(9+MSGDBL,tmp,mdtob(msglen),dest,1)
#ifdef MPI4
                  stupid_msglen = msglen
                  stupid_type   = msgtype
                  stupid_taskid = dest
                  call MPI_SEND(tmp,
     >                     stupid_msglen,stupid_complex,
     >                     stupid_taskid,
     >                     stupid_type,stupid_world,stupid_ierr)
#else
                  call MPI_SEND(tmp,
     >                     msglen,MPI_DOUBLE_PRECISION,
     >                     dest,
     >                     msgtype,Parallel_comm_world(),mpierr)
#endif
            end if

         end do
      end if



      !*************************
      !**** hilbert mapping ****
      !*************************
      else


*     **** master node reads from file and distributes ****
      if (taskid.eq.MASTER) then
         do k=1,nz(nb)
         do j=1,ny(nb)

            q      = int_mb(q_map1(1,nb)+(j-1)+(k-1)*ny(nb))
            p_from = int_mb(p_map1(1,nb)+(j-1)+(k-1)*ny(nb))


            if (p_from.eq.MASTER) then
              index = (q-1)*(nx(nb)+2) + 1
              call dcopy(nx(nb),A(index),1,tmp,1)

            else
               msglen  = (nx(nb))
               msgtype = 179
               source  = p_from
c               call RCV(9+MSGDBL,tmp,mdtob(msglen),rcv_len,
c     >                       source,rcv_proc,1)
#ifdef MPI4
               stupid_msglen = msglen
               stupid_type   = msgtype
               stupid_taskid = source
               call MPI_RECV(tmp,
     >                    stupid_msglen,stupid_double,
     >                    stupid_taskid,
     >                    stupid_type,stupid_world,
     >                    int_mb(status(1)),stupid_ierr)
#else
               call MPI_RECV(tmp,
     >                    msglen,MPI_DOUBLE_PRECISION,
     >                    source,
     >                    msgtype,Parallel_comm_world(),
     >                    int_mb(status(1)),mpierr)
#endif


            end if

            !call dwrite(iunit,tmp,(nx(nb)))
            write(iunit,'(3E26.14)') (tmp(i), i=1,nx(nb))
         end do
         end do

*     **** not master node ****
      else
         do k=1,nz(nb)
         do j=1,ny(nb)

            q      = int_mb(q_map1(1,nb)+(j-1)+(k-1)*ny(nb))
            p_here = int_mb(p_map1(1,nb)+(j-1)+(k-1)*ny(nb))

            if (p_here.eq.taskid) then

               index = (q-1)*(nx(nb)+2) + 1
               call dcopy(nx(nb),A(index),1,tmp,1)

               msglen  = nx(nb)
               msgtype = 179
               dest    = MASTER

c               call SND(9+MSGDBL,tmp,mdtob(msglen),dest,1)
#ifdef MPI4
                  stupid_msglen = msglen
                  stupid_type   = msgtype
                  stupid_taskid = dest
                  call MPI_SEND(tmp,
     >                     stupid_msglen,stupid_double,
     >                     stupid_taskid,
     >                     stupid_type,stupid_world,stupid_ierr)
#else
                  call MPI_SEND(tmp,
     >                     msglen,MPI_DOUBLE_PRECISION,
     >                     dest,
     >                     msgtype,Parallel_comm_world(),mpierr)
#endif
            end if

         end do
         end do
      end if

      end if

*     *** deallocate status memory ***
      if (.not.BA_pop_stack(status(2)))
     > call errquit('D3dB_r_FormatWrite:popstack',0,MA_ERR)

*     **** wait ****
      return
      end


*     *******************************************
*     *						*
*     *	       D3dB_r_FormatWrite_reverse	*	
*     *						*
*     *******************************************

      subroutine D3dB_r_FormatWrite_reverse(nb,iunit,A,tmp)
      implicit none 
      integer nb
      integer iunit
      real*8     A(*)
      real*8     tmp(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "mpif.h"
#include "D3dB.fh"

#ifdef MPI4
#include "stupid_mpi4.fh"
#endif

*     *** local variables ***
      integer MASTER,taskid
      parameter(MASTER=0)
      integer p_from, p_here,q
      integer i,j,k,index
      integer dest,source,msglen,idum,msgtype,status(2),mpierr
      integer  Parallel_comm_world
      external Parallel_comm_world

      call Parallel_taskid(taskid)

*     **** allocate status memory ****
      if (.not.BA_push_get(mt_int,MPI_STATUS_SIZE*4,
     >                     'status',status(2),status(1)))
     > call errquit('D3dB_r_FormatWrite_reverse:out stack',0,MA_ERR)


      !**********************
      !**** slab mapping ****
      !**********************
      if (mapping.eq.1) then
*     **** master node reads from file and distributes ****
      if (taskid.eq.MASTER) then
         do i=1,nx(nb)
         do j=1,ny(nb)

            do k=1,nz(nb)
              call D3dB_ktoqp(nb,k,q,p_from)
              if (p_from.eq.MASTER) then
                 index = (q-1)*(nx(nb)+2)*ny(nb) 
     >                 + (j-1)*(nx(nb)+2) + i
                 tmp(k) = A(index)
              else
                 msglen  = 1
                 msgtype = 189
                 source  = p_from

c                 call SND(9+MSGINT,idum,mitob(msglen),source,1)
c                 call RCV(9+MSGDBL,tmp(k),mdtob(msglen),rcv_len,
c     >                         source,rcv_proc,1)
#ifdef MPI4
               stupid_msglen = msglen
               stupid_type   = msgtype
               stupid_taskid = source
               call MPI_SEND(idum,
     >                     stupid_msglen,stupid_integer,
     >                     stupid_taskid,
     >                     stupid_type,stupid_world,stupid_ierr)
               call MPI_RECV(tmp(k),
     >                    stupid_msglen,stupid_double,
     >                    stupid_taskid,
     >                    stupid_type,stupid_world,
     >                    int_mb(status(1)),stupid_ierr)
#else
               call MPI_SEND(idum,
     >                     msglen,MPI_INTEGER,
     >                     source,
     >                     msgtype,Parallel_comm_world(),mpierr)
               call MPI_RECV(tmp(k),
     >                    msglen,MPI_DOUBLE_PRECISION,
     >                    source,
     >                    msgtype,Parallel_comm_world(),
     >                    int_mb(status(1)),mpierr)
#endif


              end if
            end do
            write(iunit,'(6E13.5)') (tmp(k), k=1,nz(nb))
       
         end do
         end do

*     **** not master node ****
      else
         do i=1,nx(nb)
         do j=1,ny(nb)

            do k=1,nz(nb)
              call D3dB_ktoqp(nb,k,q,p_here)
              if (p_here.eq.taskid) then

                 index = (q-1)*(nx(nb)+2)*ny(nb) 
     >                 + (j-1)*(nx(nb)+2) + i
                 tmp(1) = A(index)

                 msglen  = 1
                 msgtype = 189
                 dest    = MASTER
c                 call RCV(9+MSGINT,idum,mitob(msglen),rcv_len,
c     >                         dest,rcv_proc,1)
c                 call SND(9+MSGDBL,tmp,mdtob(msglen),dest,1)
#ifdef MPI4
                 stupid_msglen = msglen
                 stupid_type   = msgtype
                 stupid_taskid = dest
                 call MPI_RECV(idum,
     >                    stupid_msglen,stupid_integer,
     >                    stupid_taskid,
     >                    stupid_type,stupid_world,
     >                    int_mb(status(1)),stupid_ierr)
                 call MPI_SEND(tmp,
     >                     stupid_msglen,stupid_double,
     >                     stupid_taskid,
     >                     stupid_type,stupid_world,stupid_ierr)
#else
                 call MPI_RECV(idum,
     >                    msglen,MPI_INTEGER,
     >                    dest,
     >                    msgtype,Parallel_comm_world(),
     >                    int_mb(status(1)),mpierr)
                 call MPI_SEND(tmp,
     >                     msglen,MPI_DOUBLE_PRECISION,
     >                     dest,
     >                     msgtype,Parallel_comm_world(),mpierr)
#endif
              end if
            end do

         end do
         end do
      end if


      !*************************
      !**** hilbert mapping ****
      !*************************
      else

*     **** master node reads from file and distributes ****
      if (taskid.eq.MASTER) then
         do i=1,nx(nb)
         do j=1,ny(nb)

            do k=1,nz(nb)
              call D3dB_ijktoindex2p(nb,i,j,k,index,p_from)
              if (p_from.eq.MASTER) then
                 tmp(k) = A(index)
              else
                 msglen  = 1
                 msgtype = 189
                 source  = p_from
c                 call SND(9+MSGINT,idum,mitob(msglen),source,1)
c                 call RCV(9+MSGDBL,tmp(k),mdtob(msglen),rcv_len,
c     >                         source,rcv_proc,1)
#ifdef MPI4
                 stupid_msglen = msglen
                 stupid_type   = msgtype
                 stupid_taskid = source
                 call MPI_SEND(idum,
     >                     stupid_msglen,stupid_integer,
     >                     stupid_taskid,
     >                     stupid_type,stupid_world,stupid_ierr)
                 call MPI_RECV(tmp(k),
     >                    stupid_msglen,stupid_double,
     >                    stupid_taskid,
     >                    stupid_type,stupid_world,
     >                    int_mb(status(1)),stupid_ierr)
#else
                 call MPI_SEND(idum,
     >                     msglen,MPI_INTEGER,
     >                     source,
     >                     msgtype,Parallel_comm_world(),mpierr)
                 call MPI_RECV(tmp(k),
     >                    msglen,MPI_DOUBLE_PRECISION,
     >                    source,
     >                    msgtype,Parallel_comm_world(),
     >                    int_mb(status(1)),mpierr)
#endif
              end if
            end do
            write(iunit,'(6E13.5)') (tmp(k), k=1,nz(nb))

         end do
         end do

*     **** not master node ****
      else
         do i=1,nx(nb)
         do j=1,ny(nb)

            do k=1,nz(nb)
              call D3dB_ijktoindex2p(nb,i,j,k,index,p_here)
              if (p_here.eq.taskid) then

                 tmp(1) = A(index)

                 msglen  = 1
                 msgtype = 189
                 dest    = MASTER
c                 call RCV(9+MSGINT,idum,mitob(msglen),rcv_len,
c     >                         dest,rcv_proc,1)
c                 call SND(9+MSGDBL,tmp,mdtob(msglen),dest,1)
#ifdef MPI4
                 stupid_msglen = msglen
                 stupid_type   = msgtype
                 stupid_taskid = dest
                 call MPI_RECV(idum,
     >                    stupid_msglen,stupid_integer,
     >                    stupid_taskid,
     >                    stupid_type,stupid_world,
     >                    int_mb(status(1)),stupid_ierr)
                 call MPI_SEND(tmp,
     >                     stupid_msglen,stupid_double,
     >                     stupid_taskid,
     >                     stupid_type,stupid_world,stupid_ierr)
#else
                 call MPI_RECV(idum,
     >                    msglen,MPI_INTEGER,
     >                    dest,
     >                    msgtype,Parallel_comm_world(),
     >                    int_mb(status(1)),mpierr)
                 call MPI_SEND(tmp,
     >                     msglen,MPI_DOUBLE_PRECISION,
     >                     dest,
     >                     msgtype,Parallel_comm_world(),mpierr)
#endif
              end if
            end do

         end do
         end do
      end if

      end if

*     *** deallocate status memory ***
      if (.not.BA_pop_stack(status(2)))
     > call errquit('D3dB_r_FormatWrite_reverse:popstack',0,MA_ERR)

*     **** wait ****
      return
      end

*     ***********************************
*     *					*
*     *	       D3dB_r_FormatWrite_paw	*	
*     *					*
*     ***********************************

      subroutine D3dB_r_FormatWrite_paw(nb,iunit,A,tmp)
      implicit none 
      integer nb
      integer iunit
      real*8     A(*)
      real*8     tmp(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "mpif.h"
#include "D3dB.fh"

#ifdef MPI4
#include "stupid_mpi4.fh"
#endif

      
*     *** local variables ***
      integer MASTER,taskid
      parameter(MASTER=0)
      integer p_from, p_here,q
      integer i,j,k,index
      integer dest,source,msglen,msgtype,status(2),mpierr
      integer  Parallel_comm_world
      external Parallel_comm_world

      call Parallel_taskid(taskid)

*     **** allocate status memory ****
      if (.not.BA_push_get(mt_int,MPI_STATUS_SIZE*4,
     >                     'status',status(2),status(1)))
     > call errquit('D3dB_r_FormatWrite_paw:out of stack',0,MA_ERR)


*     **** master node reads from file and distributes ****
      if (taskid.eq.MASTER) then
         do k=1,nz(nb)
         do j=1,ny(nb)

            call D3dB_ktoqp(nb,j,q,p_from)

            if (p_from.eq.MASTER) then

                 index = (q-1)*(nx(nb)+2)*ny(nb) 
     >                 + (k-1)*(nx(nb)+2) + 1
                 call dcopy(nx(nb),A(index),1,tmp,1)

            else

               msglen  = (nx(nb))
               msgtype = 199
               status  = msglen
               source  = p_from

c               call RCV(9+MSGDBL,tmp,mdtob(msglen),rcv_len,
c     >                       source,rcv_proc,1)
#ifdef MPI4
               stupid_msglen = msglen
               stupid_type   = msgtype
               stupid_taskid = source
               call MPI_RECV(tmp,
     >                    stupid_msglen,stupid_double,
     >                    stupid_taskid,
     >                    stupid_type,stupid_world,
     >                    int_mb(status(1)),stupid_ierr)
#else
               call MPI_RECV(tmp,
     >                    msglen,MPI_DOUBLE_PRECISION,
     >                    source,
     >                    msgtype,Parallel_comm_world(),
     >                    int_mb(status(1)),mpierr)
#endif


            end if

c           **** call dwrite(iunit,tmp,(nx(nb))) ****
            write(iunit,'(3E26.14)') (tmp(i), i=1,nx(nb))

         end do
         end do

*     **** not master node ****
      else
         do k=1,nz(nb)
         do j=1,ny(nb)
         
            call D3dB_ktoqp(nb,j,q,p_here)
            if (p_here.eq.taskid) then

                  index = (q-1)*(nx(nb)+2)*ny(nb) 
     >                  + (k-1)*(nx(nb)+2) + 1
                  call dcopy(nx(nb),A(index),1,tmp,1)


               msglen  = (nx(nb))
               msgtype = 199
               dest    = MASTER
c               call SND(9+MSGDBL,tmp,mdtob(msglen),dest,1)
#ifdef MPI4
                  stupid_msglen = msglen
                  stupid_type   = msgtype
                  stupid_taskid = dest
                  call MPI_SEND(tmp,
     >                     stupid_msglen,stupid_double,
     >                     stupid_taskid,
     >                     stupid_type,stupid_world,stupid_ierr)
#else
                  call MPI_SEND(tmp,
     >                     msglen,MPI_DOUBLE_PRECISION,
     >                     dest,
     >                     msgtype,Parallel_comm_world(),mpierr)
#endif
            end if

         end do
         end do
      end if



*     *** deallocate status memory ***
      if (.not.BA_pop_stack(status(2)))
     > call errquit('D3dB_r_FormatWrite_paw:popstack',0,MA_ERR)

*     **** wait ****
      return
      end




*     ***********************************
*     *                                 *
*     *       D3dB_r_1to3_Contract      *
*     *                                 *
*     ***********************************
      subroutine D3dB_r_1to3_Contract(center,A,tmp1,tmp2)
      implicit none
      integer center(3)
      !real*8 A(*),B(*)
      real*8 A(*)
      real*8 tmp1(*),tmp2(*)

#include "bafdecls.fh"
#include "mpif.h"
#include "D3dB.fh"
#include "errquit.fh"

#ifdef MPI4
#include "stupid_mpi4.fh"
#endif

*     **** local variables ****
      logical value
      integer taskid,np
      integer proc_to,proc_from,source
      integer indxbig,pbig,indxsmall,psmall,itmp
      integer index1,index2,msgtype,msglen
      integer k1,k2,k3,l1,l2,l3,i,j,k,il,jl,kl,it
      integer nxh,nyh,nzh,nxh_small,nyh_small,nzh_small
      integer pindex1(2),pindex2(2),ic1_start(2),ic2_start(2)
      integer c_iq_to_i1(2),c_iq_to_i2(2)
      integer reqcnt,request(2)
      integer mpierr

*     **** external functions ****
      integer  Parallel2d_convert_taskid_i,Parallel_comm_world
      external Parallel2d_convert_taskid_i,Parallel_comm_world


      call Parallel2d_taskid_i(taskid)
      call Parallel2d_np_i(np)

      !**** allocate temporary space from stack ****
c     pindex1(np),pindex2(np),ic1_start(np+1),ic2_start(np+1)
c     c_iq_to_i1(n2ft3d(1)),c_iq_to_i2(n2ft3d(1))
      value = BA_push_get(mt_int,np,'request_1t3',request(2),request(1))
      value = value.and.
     >        BA_push_get(mt_int,np,'pindex1_1t3',pindex1(2),pindex1(1))
      value = value.and.
     >        BA_push_get(mt_int,np,'pindex2_1t3',pindex2(2),pindex2(1))
      value = value.and.
     >        BA_push_get(mt_int,np+1,'ic1_start_1t3',
     >                    ic1_start(2),ic1_start(1))
      value = value.and.
     >        BA_push_get(mt_int,np+1,'ic2_start_1t3',
     >                    ic2_start(2),ic2_start(1))
      value = value.and.
     >        BA_push_get(mt_int,n2ft3d(1),'c_iq_to_i1_1t3',
     >                    c_iq_to_i1(2),c_iq_to_i1(1))
      value = value.and.
     >        BA_push_get(mt_int,n2ft3d(1),'c_iq_to_i2_1t3',
     >                    c_iq_to_i2(2),c_iq_to_i2(1))
      if (.not.value)
     >   call errquit('D3dB_r_1to3_Contract:out of stack',0,MA_ERR)


      nxh = nx(1)/2
      nyh = ny(1)/2
      nzh = nz(1)/2

      nxh_small = nx(3)/2
      nyh_small = ny(3)/2
      nzh_small = nz(3)/2

!$OMP MASTER
      do it=0,np-1
         int_mb(pindex1(1)+it) = 0
         int_mb(pindex2(1)+it) = 0
      end do
      do k3=-nzh_small,nzh_small-1
         l3 = k3 + center(3)
         if (l3.lt.(-nzh))  l3 = l3 + nz(1)
         if (l3.gt.(nzh-1)) l3 = l3 - nz(1)
         k  = k3 + nzh_small
         kl = l3 + nzh

         do k2=-nyh_small,nyh_small-1
            l2 = k2 + center(2)
            if (l2.lt.(-nyh))  l2 = l2 + ny(1)
            if (l2.gt.(nyh-1)) l2 = l2 - ny(1)
            j  = k2 + nyh_small
            jl = l2 + nyh

            do k1=-nxh_small,nxh_small-1
               l1 = k1 + center(1)
               if (l1.lt.(-nxh))  l1 = l1 + nx(1)
               if (l1.gt.(nxh-1)) l1 = l1 - nx(1)
               i  = k1 + nxh_small
               il = l1 + nxh

               call D3dB_ijktoindex2p(1,il+1,jl+1,kl+1,indxbig,pbig)
               call D3dB_ijktoindex2p(3,i+1,j+1,k+1,indxsmall,psmall)

*              **** packing scheme **** 
               if (pbig.eq.taskid) then
                  int_mb(pindex1(1)+psmall)=int_mb(pindex1(1)+psmall)+1
               end if

*              **** unpacking scheme **** 
               if (psmall.eq.taskid) then
                  int_mb(pindex2(1)+pbig) = int_mb(pindex2(1)+pbig)+1
               end if

            end do
         end do
      end do
      index1 = 1
      index2 = 1
      do it=0,np-1
         proc_to   = mod(taskid+it,np)
         proc_from = mod(taskid-it+np,np)
         int_mb(ic1_start(1)+it) = index1
         int_mb(ic2_start(1)+it) = index2
         index1 = index1 + int_mb(pindex1(1)+proc_to)
         index2 = index2 + int_mb(pindex2(1)+proc_from)
      end do
      int_mb(ic1_start(1)+np) = index1
      int_mb(ic2_start(1)+np) = index2

      do it=0,np-1
         int_mb(pindex1(1)+it) = 0
         int_mb(pindex2(1)+it) = 0
      end do
      do k3=-nzh_small,nzh_small-1
         l3 = k3 + center(3)
         if (l3.lt.(-nzh))  l3 = l3 + nz(1)
         if (l3.gt.(nzh-1)) l3 = l3 - nz(1)
         k  = k3 + nzh_small
         kl = l3 + nzh

         do k2=-nyh_small,nyh_small-1
            l2 = k2 + center(2)
            if (l2.lt.(-nyh))  l2 = l2 + ny(1)
            if (l2.gt.(nyh-1)) l2 = l2 - ny(1)
            j  = k2 + nyh_small
            jl = l2 + nyh

            do k1=-nxh_small,nxh_small-1
               l1 = k1 + center(1)
               if (l1.lt.(-nxh))  l1 = l1 + nx(1)
               if (l1.gt.(nxh-1)) l1 = l1 - nx(1)
               i  = k1 + nxh_small
               il = l1 + nxh

               call D3dB_ijktoindex2p(1,il+1,jl+1,kl+1,indxbig,pbig)
               call D3dB_ijktoindex2p(3,i+1,j+1,k+1,indxsmall,psmall)

*              **** packing scheme **** 
               if (pbig.eq.taskid) then
                  itmp = int_mb(ic1_start(1)+psmall) 
     >                 + int_mb(pindex1(1)+psmall)
                  int_mb(c_iq_to_i1(1)+itmp-1) = indxbig
                  int_mb(pindex1(1)+psmall)=int_mb(pindex1(1)+psmall)+1
               end if

*              **** unpacking scheme **** 
               if (psmall.eq.taskid) then
                  itmp = int_mb(ic2_start(1)+pbig) 
     >                 + int_mb(pindex2(1)+pbig)
                  int_mb(c_iq_to_i2(1)+itmp-1) = indxsmall
                  int_mb(pindex2(1)+pbig) = int_mb(pindex2(1)+pbig)+1
               end if

            end do
         end do
      end do
!$OMP END MASTER
!$OMP BARRIER


*     **** pack A(i) array ****
!$OMP DO
      do i=1,index1-1
         tmp1(i) = A(int_mb(c_iq_to_i1(1)+i-1)) 
      end do
!$OMP END DO

*     **** it = 0, transpose data on same thread ****  
      msglen = int_mb(ic2_start(1)+2-1) - int_mb(ic2_start(1)+1-1)

      call Parallel_shared_vector_copy(.false.,msglen,
     >               tmp1(int_mb(ic1_start(1)+1-1)),
     >               tmp2(int_mb(ic2_start(1)+1-1)))

!$OMP MASTER

*     **** receive packed array data ****
      reqcnt = 0
      do it = 1,np-1

*        **** synchronous receive of tmp ****
         proc_from = mod((taskid)-it+np, np)

         msgtype = 93
         source=proc_from
         msglen = (int_mb(ic2_start(1)+it+2-1)
     >          -  int_mb(ic2_start(1)+it+1-1))



         if (msglen.gt.0) then
            reqcnt = reqcnt + 1
#ifdef MPI4
            stupid_msglen = msglen
            stupid_type   = msgtype
            stupid_taskid = Parallel2d_convert_taskid_i(source)
            call MPI_IRECV(tmp2(int_mb(ic2_start(1)+it+1-1)),
     >                    stupid_msglen,stupid_double,
     >                    stupid_taskid,
     >                    stupid_type,stupid_world,
     >                    stupid_request,stupid_ierr)
            int_mb(request(1)+reqcnt-1) = 0
            int_mb(request(1)+reqcnt-1) = stupid_request
#else
            call MPI_IRECV(tmp2(int_mb(ic2_start(1)+it+1-1)),
     >                    msglen,MPI_DOUBLE_PRECISION,
     >                    Parallel2d_convert_taskid_i(source),
     >                    msgtype,Parallel_comm_world(),
     >                    int_mb(request(1)+reqcnt-1),mpierr)
#endif
         end if

      end do

*     **** send packed array to other processors ****
      do it = 1,np-1

*        **** synchronous send of tmp ****
         proc_to   = mod((taskid)+it, np)
         msglen    = (int_mb(ic1_start(1)+it+2-1)
     >              - int_mb(ic1_start(1)+it+1-1))
         msgtype   = 93

         if (msglen.gt.0) then
#ifdef MPI4
            stupid_msglen = msglen
            stupid_type   = msgtype
            stupid_taskid = Parallel2d_convert_taskid_i(proc_to)
            call MPI_SEND(tmp1(int_mb(ic1_start(1)+it+1-1)),
     >                     stupid_msglen,stupid_double,
     >                     stupid_taskid,
     >                     stupid_type,stupid_world,stupid_ierr)
#else
            call MPI_SEND(tmp1(int_mb(ic1_start(1)+it+1-1)),
     >                     msglen,MPI_DOUBLE_PRECISION,
     >                     Parallel2d_convert_taskid_i(proc_to),
     >                     msgtype,Parallel_comm_world(),mpierr)
#endif
         end if

      end do


*     **** wait for completion of mp_send, also do a sync ****
      if (np.gt.1) call Parallel_mpiWaitAll0(reqcnt,int_mb(request(1)))

!$OMP END MASTER
!$OMP BARRIER

*     **** unpack B(i) array ****
!$OMP DO
      do i=1,index2-1
         A(int_mb(c_iq_to_i2(1)+i-1)) = tmp2(i)
      end do
!$OMP END DO


      !**** deallocate temporary space from stack ****
      value =           BA_pop_stack(c_iq_to_i2(2))
      value = value.and.BA_pop_stack(c_iq_to_i1(2))
      value = value.and.BA_pop_stack(ic2_start(2))
      value = value.and.BA_pop_stack(ic1_start(2))
      value = value.and.BA_pop_stack(pindex2(2))
      value = value.and.BA_pop_stack(pindex1(2))
      value = value.and.BA_pop_stack(request(2))
      if (.not.value)
     >   call errquit('D3dB_r_1to3_Contract:popping stack',0,MA_ERR)

      return
      end






*     ***********************************
*     *                                 *
*     *       D3dB_r_3to1_Expand        *
*     *                                 *
*     ***********************************
      subroutine D3dB_r_3to1_Expand(center,A,tmp1,tmp2)
      implicit none
      integer center(3)
      !real*8 A(*),B(*)
      real*8 A(*)
      real*8 tmp1(*),tmp2(*)

#include "bafdecls.fh"
#include "mpif.h"
#include "D3dB.fh"
#include "errquit.fh"

#ifdef MPI4
#include "stupid_mpi4.fh"
#endif


*     **** local variables ****
      logical value
      integer taskid,np
      integer proc_to,proc_from,source
      integer indxbig,pbig,indxsmall,psmall,itmp
      integer index1,index2,msgtype,msglen
      integer k1,k2,k3,l1,l2,l3,i,j,k,il,jl,kl,it
      integer nxh,nyh,nzh,nxh_small,nyh_small,nzh_small
      integer pindex1(2),pindex2(2),ic1_start(2),ic2_start(2)
      integer c_iq_to_i1(2),c_iq_to_i2(2)
      integer reqcnt,request(2)
      integer mpierr

*     **** external functions ****
      integer  Parallel2d_convert_taskid_i,Parallel_comm_world
      external Parallel2d_convert_taskid_i,Parallel_comm_world


      call Parallel2d_taskid_i(taskid)
      call Parallel2d_np_i(np)

      !**** allocate temporary space from stack ****
      value = BA_push_get(mt_int,np,'request',request(2),request(1))
      value = value.and.
     >        BA_push_get(mt_int,np,'pindex1',pindex1(2),pindex1(1))
      value = value.and.
     >        BA_push_get(mt_int,np,'pindex2',pindex2(2),pindex2(1))
      value = value.and.
     >        BA_push_get(mt_int,np+1,'ic1_start',
     >                    ic1_start(2),ic1_start(1))
      value = value.and.
     >        BA_push_get(mt_int,np+1,'ic2_start',
     >                    ic2_start(2),ic2_start(1))
      value = value.and.
     >        BA_push_get(mt_int,n2ft3d(1),'c_iq_to_i1',
     >                    c_iq_to_i1(2),c_iq_to_i1(1))
      value = value.and.
     >        BA_push_get(mt_int,n2ft3d(1),'c_iq_to_i2',
     >                    c_iq_to_i2(2),c_iq_to_i2(1))
      if (.not.value)
     >   call errquit('D3dB_r_3to1_Expand:out of stack',0,MA_ERR)


      nxh = nx(1)/2
      nyh = ny(1)/2
      nzh = nz(1)/2

      nxh_small = nx(3)/2
      nyh_small = ny(3)/2
      nzh_small = nz(3)/2

!$OMP MASTER
      do it=0,np-1
         int_mb(pindex1(1)+it) = 0
         int_mb(pindex2(1)+it) = 0
      end do
      do k3=-nzh_small,nzh_small-1
         l3 = k3 + center(3)
         if (l3.lt.(-nzh))  l3 = l3 + nz(1)
         if (l3.gt.(nzh-1)) l3 = l3 - nz(1)
         k  = k3 + nzh_small
         kl = l3 + nzh

         do k2=-nyh_small,nyh_small-1
            l2 = k2 + center(2)
            if (l2.lt.(-nyh))  l2 = l2 + ny(1)
            if (l2.gt.(nyh-1)) l2 = l2 - ny(1)
            j  = k2 + nyh_small
            jl = l2 + nyh

            do k1=-nxh_small,nxh_small-1
               l1 = k1 + center(1)
               if (l1.lt.(-nxh))  l1 = l1 + nx(1)
               if (l1.gt.(nxh-1)) l1 = l1 - nx(1)
               i  = k1 + nxh_small
               il = l1 + nxh

               call D3dB_ijktoindex2p(1,il+1,jl+1,kl+1,indxbig,pbig)
               call D3dB_ijktoindex2p(3,i+1,j+1,k+1,indxsmall,psmall)

*              **** packing scheme **** 
               if (psmall.eq.taskid) then
                  int_mb(pindex1(1)+pbig) = int_mb(pindex1(1)+pbig)+1
               end if

*              **** unpacking scheme **** 
               if (pbig.eq.taskid) then
                  int_mb(pindex2(1)+psmall)=int_mb(pindex2(1)+psmall)+1
               end if

            end do
         end do
      end do
      index1 = 1
      index2 = 1
      do it=0,np-1
         proc_to   = mod(taskid+it,np)
         proc_from = mod(taskid-it+np,np)
         int_mb(ic1_start(1)+it) = index1
         int_mb(ic2_start(1)+it) = index2
         index1 = index1 + int_mb(pindex1(1)+proc_to)
         index2 = index2 + int_mb(pindex2(1)+proc_from)
      end do
      int_mb(ic1_start(1)+np) = index1
      int_mb(ic2_start(1)+np) = index2

      do it=0,np-1
         int_mb(pindex1(1)+it) = 0
         int_mb(pindex2(1)+it) = 0
      end do
      do k3=-nzh_small,nzh_small-1
         l3 = k3 + center(3)
         if (l3.lt.(-nzh))  l3 = l3 + nz(1)
         if (l3.gt.(nzh-1)) l3 = l3 - nz(1)
         k  = k3 + nzh_small
         kl = l3 + nzh

         do k2=-nyh_small,nyh_small-1
            l2 = k2 + center(2)
            if (l2.lt.(-nyh))  l2 = l2 + ny(1)
            if (l2.gt.(nyh-1)) l2 = l2 - ny(1)
            j  = k2 + nyh_small
            jl = l2 + nyh

            do k1=-nxh_small,nxh_small-1
               l1 = k1 + center(1)
               if (l1.lt.(-nxh))  l1 = l1 + nx(1)
               if (l1.gt.(nxh-1)) l1 = l1 - nx(1)
               i  = k1 + nxh_small
               il = l1 + nxh

               call D3dB_ijktoindex2p(1,il+1,jl+1,kl+1,indxbig,pbig)
               call D3dB_ijktoindex2p(3,i+1,j+1,k+1,indxsmall,psmall)

*              **** packing scheme **** 
               if (psmall.eq.taskid) then
                  itmp = int_mb(ic1_start(1)+pbig) 
     >                 + int_mb(pindex1(1)+pbig)
                  int_mb(c_iq_to_i1(1)+itmp-1) = indxsmall
                  int_mb(pindex1(1)+pbig) = int_mb(pindex1(1)+pbig)+1
               end if

*              **** unpacking scheme **** 
               if (pbig.eq.taskid) then
                  itmp = int_mb(ic2_start(1)+psmall) 
     >                 + int_mb(pindex2(1)+psmall)
                  int_mb(c_iq_to_i2(1)+itmp-1) = indxbig
                  int_mb(pindex2(1)+psmall)=int_mb(pindex2(1)+psmall)+1
               end if


            end do
         end do
      end do
!$OMP END MASTER
!$OMP BARRIER


*     **** pack A(i) array ****
!$OMP DO
      do i=1,index1-1
         tmp1(i) = A(int_mb(c_iq_to_i1(1)+i-1)) 
      end do
!$OMP END DO

*     **** it = 0, transpose data on same thread ****  
      msglen = int_mb(ic2_start(1)+2-1) - int_mb(ic2_start(1)+1-1)

      call Parallel_shared_vector_copy(.false.,msglen,
     >               tmp1(int_mb(ic1_start(1)+1-1)),
     >               tmp2(int_mb(ic2_start(1)+1-1)))

      call Parallel_shared_vector_zero(.false.,n2ft3d(1),A)

!$OMP MASTER
*     **** receive packed array data ****
      reqcnt = 0
      do it = 1,np-1

*        **** synchronous receive of tmp ****
         proc_from = mod((taskid)-it+np, np)

         msgtype = 94
         source=proc_from
         msglen = (int_mb(ic2_start(1)+it+2-1)
     >          -  int_mb(ic2_start(1)+it+1-1))


         if (msglen.gt.0) then
            reqcnt = reqcnt + 1
#ifdef MPI4
            stupid_msglen = msglen
            stupid_type   = msgtype
            stupid_taskid = Parallel2d_convert_taskid_i(source)
            call MPI_IRECV(tmp2(int_mb(ic2_start(1)+it+1-1)),
     >                    stupid_msglen,stupid_double,
     >                    stupid_taskid,
     >                    stupid_type,stupid_world,
     >                    stupid_request,stupid_ierr)
            int_mb(request(1)+reqcnt-1) = 0
            int_mb(request(1)+reqcnt-1) = stupid_request
#else
            call MPI_IRECV(tmp2(int_mb(ic2_start(1)+it+1-1)),
     >                    msglen,MPI_DOUBLE_PRECISION,
     >                    Parallel2d_convert_taskid_i(source),
     >                    msgtype,Parallel_comm_world(),
     >                    int_mb(request(1)+reqcnt-1),mpierr)
#endif
         end if

      end do

*     **** send packed array to other processors ****
      do it = 1,np-1

*        **** synchronous send of tmp ****
         proc_to   = mod((taskid)+it, np)
         msglen    = (int_mb(ic1_start(1)+it+2-1)
     >              - int_mb(ic1_start(1)+it+1-1))
         msgtype   = 94

         if (msglen.gt.0) then
#ifdef MPI4
            stupid_msglen = msglen
            stupid_type   = msgtype
            stupid_taskid = Parallel2d_convert_taskid_i(proc_to)
            call MPI_SEND(tmp1(int_mb(ic1_start(1)+it+1-1)),
     >                     stupid_msglen,stupid_double,
     >                     stupid_taskid,
     >                     stupid_type,stupid_world,stupid_ierr)
#else
            call MPI_SEND(tmp1(int_mb(ic1_start(1)+it+1-1)),
     >                     msglen,MPI_DOUBLE_PRECISION,
     >                     Parallel2d_convert_taskid_i(proc_to),
     >                     msgtype,Parallel_comm_world(),mpierr)
#endif
         end if

      end do


*     **** wait for completion of mp_send, also do a sync ****
      if (np.gt.1) call Parallel_mpiWaitAll0(reqcnt,int_mb(request(1)))

!$OMP END MASTER
!$OMP BARRIER

*     **** unpack B(i) array ****
!$OMP DO
      do i=1,index2-1
         A(int_mb(c_iq_to_i2(1)+i-1)) = tmp2(i)
      end do
!$OMP END DO


      !**** deallocate temporary space from stack ****
      value =           BA_pop_stack(c_iq_to_i2(2))
      value = value.and.BA_pop_stack(c_iq_to_i1(2))
      value = value.and.BA_pop_stack(ic2_start(2))
      value = value.and.BA_pop_stack(ic1_start(2))
      value = value.and.BA_pop_stack(pindex2(2))
      value = value.and.BA_pop_stack(pindex1(2))
      value = value.and.BA_pop_stack(request(2))
      if (.not.value)
     >   call errquit('D3dB_r_3to1_Expand:popping stack',0,MA_ERR)

      return
      end



*     ***********************************
*     *                                 *
*     *         D3dB_r_1to3_init        *
*     *                                 *
*     ***********************************
      subroutine D3dB_r_1to3_init()
      implicit none

#include "bafdecls.fh"
#include "D3dB.fh"
#include "errquit.fh"

*     **** r_1to3_block common block ****
      integer maxq_1t3
      parameter (maxq_1t3=6)
      integer qhead0_1t3,qhead1_1t3,qlast_1t3
      integer i1_start_1t3(2),i2_start_1t3(2)
      integer c_iq_to_i1_1t3(2),c_iq_to_i2_1t3(2)
      integer index2_1t3(2)
      common /r_1to3_block/ i1_start_1t3,i2_start_1t3,
     >                      c_iq_to_i1_1t3,c_iq_to_i2_1t3,
     >                      index2_1t3,
     >                      qhead0_1t3,qhead1_1t3,qlast_1t3

*     **** igrid_r_block common block ****
      integer igrid_big(2),igrid_small(2),small(2),pindex1(2),pindex2(2)
      common /igrid_r_block/ igrid_big,igrid_small,small,pindex1,pindex2

*     **** local variables ****
      logical value
      integer np


      call Parallel2d_np_i(np)

      !**** allocate indexing space from heap ****
      value = BA_alloc_get(mt_int,maxq_1t3*(np+1),'i1_start_1t3',
     >                    i1_start_1t3(2),i1_start_1t3(1))
      value = value.and.
     >        BA_alloc_get(mt_int,maxq_1t3*(np+1),'i2_start_1t3',
     >                    i2_start_1t3(2),i2_start_1t3(1))
      value = value.and.
     >        BA_alloc_get(mt_int,maxq_1t3*n2ft3d(1),'c_iq_to_i1_1t3',
     >                    c_iq_to_i1_1t3(2),c_iq_to_i1_1t3(1))
      value = value.and.
     >        BA_alloc_get(mt_int,maxq_1t3*n2ft3d(1),'c_iq_to_i2_1t3',
     >                    c_iq_to_i2_1t3(2),c_iq_to_i2_1t3(1))
      value = value.and.
     >        BA_alloc_get(mt_int,maxq_1t3,'index_1t3',
     >                     index2_1t3(2),index2_1t3(1))

      value = value.and.
     >        BA_alloc_get(mt_int,4*n2ft3d(1),'igrid_big',
     >                    igrid_big(2),igrid_big(1))
      value = value.and.
     >        BA_alloc_get(mt_int,4*n2ft3d(3),'igrid_small',
     >                    igrid_small(2),igrid_small(1))
      value = value.and.
     >        BA_alloc_get(mt_int,4*n2ft3d(1),'small',
     >                    small(2),small(1))
      value = value.and.
     >        BA_alloc_get(mt_int,np,'pindex1',pindex1(2),pindex1(1))
      value = value.and.
     >        BA_alloc_get(mt_int,np,'pindex2',pindex2(2),pindex2(1))
      if (.not.value)
     >   call errquit('D3dB_r_1to3_init:out of heap',0,MA_ERR)

      call lattice_i_grid(1,int_mb(igrid_big(1)))
      call lattice_i_grid(3,int_mb(igrid_small(1)))

!$OMP MASTER
      qhead0_1t3 = 0
      qhead1_1t3 = 0
      qlast_1t3  = 0
!$OMP END MASTER
!$OMP BARRIER

      return
      end 


*     ***********************************
*     *                                 *
*     *         D3dB_r_1to3_end         *
*     *                                 *
*     ***********************************
      subroutine D3dB_r_1to3_end()
      implicit none

#include "bafdecls.fh"
#include "D3dB.fh"
#include "errquit.fh"

*     **** r_1to3_block common block ****
      integer maxq_1t3
      parameter (maxq_1t3=6)
      integer qhead0_1t3,qhead1_1t3,qlast_1t3
      integer i1_start_1t3(2),i2_start_1t3(2)
      integer c_iq_to_i1_1t3(2),c_iq_to_i2_1t3(2)
      integer index2_1t3(2)
      common /r_1to3_block/ i1_start_1t3,i2_start_1t3,
     >                      c_iq_to_i1_1t3,c_iq_to_i2_1t3,
     >                      index2_1t3,
     >                      qhead0_1t3,qhead1_1t3,qlast_1t3

*     **** igrid_r_block common block ****
      integer igrid_big(2),igrid_small(2),small(2),pindex1(2),pindex2(2)
      common /igrid_r_block/ igrid_big,igrid_small,small,pindex1,pindex2

*     **** local variables ****
      logical value

      !**** allocate indexing space from heap ****
      value =           BA_free_heap(i1_start_1t3(2))
      value = value.and.BA_free_heap(i2_start_1t3(2))
      value = value.and.BA_free_heap(c_iq_to_i1_1t3(2))
      value = value.and.BA_free_heap(c_iq_to_i2_1t3(2))
      value = value.and.BA_free_heap(index2_1t3(2))
      value = value.and.BA_free_heap(igrid_big(2))
      value = value.and.BA_free_heap(igrid_small(2))
      value = value.and.BA_free_heap(small(2))
      value = value.and.BA_free_heap(pindex1(2))
      value = value.and.BA_free_heap(pindex2(2))
      if (.not.value)
     >   call errquit('D3dB_r_1to3_end:freeing heap',0,MA_ERR)

      return
      end 



*     ***********************************
*     *                                 *
*     *    D3dB_r_1to3_Contract_start3  *
*     *                                 *
*     ***********************************
      subroutine D3dB_r_1to3_Contract_start3(center,A,tmp1,
     >                                      request,reqcnt)
      implicit none
      integer center(3)
      real*8 A(*)
      real*8 tmp1(*)
      integer request(*),reqcnt

#include "bafdecls.fh"
#include "mpif.h"
#include "D3dB.fh"
#include "errquit.fh"

#ifdef MPI4
#include "stupid_mpi4.fh"
#endif

*     **** r_1to3_block common block ****
      integer maxq_1t3
      parameter (maxq_1t3=6)
      integer qhead0_1t3,qhead1_1t3,qlast_1t3
      integer i1_start_1t3(2),i2_start_1t3(2)
      integer c_iq_to_i1_1t3(2),c_iq_to_i2_1t3(2)
      integer index2_1t3(2)
      common /r_1to3_block/ i1_start_1t3,i2_start_1t3,
     >                      c_iq_to_i1_1t3,c_iq_to_i2_1t3,
     >                      index2_1t3,
     >                      qhead0_1t3,qhead1_1t3,qlast_1t3

*     **** igrid_r_block common block ****
      integer igrid_big(2),igrid_small(2),small(2),pindex1(2),pindex2(2)
      common /igrid_r_block/ igrid_big,igrid_small,small,pindex1,pindex2

*     **** local variables ****
      logical value
      integer taskid,np
      integer proc_to,proc_from,source
      integer indxbig,pbig,indxsmall,psmall,itmp
      integer index1,index2,msgtype,msglen
      integer ismall,jsmall,ksmall,lsmall,ibig,jbig,kbig,lbig
      integer it,i,it2
      integer nxh,nyh,nzh,nxh_small,nyh_small,nzh_small
      integer small_count
      integer i1_start,i2_start,c_iq_to_i1,c_iq_to_i2
      integer shift1,shift2
c      integer nthr,thr

      call nwpw_timing_start(55)
      call Parallel2d_taskid_i(taskid)
      call Parallel2d_np_i(np)
c      call Parallel_nthreads(nthr)
c      call Parallel_threadid(thr)


      shift1 = qlast_1t3*(np+1)
      shift2 = qlast_1t3*n2ft3d(1)
      i1_start = i1_start_1t3(1) + shift1
      i2_start = i2_start_1t3(1) + shift1
      c_iq_to_i1 = c_iq_to_i1_1t3(1) + shift2
      c_iq_to_i2 = c_iq_to_i2_1t3(1) + shift2


      nxh = nx(1)/2
      nyh = ny(1)/2
      nzh = nz(1)/2

      nxh_small = nx(3)/2
      nyh_small = ny(3)/2
      nzh_small = nz(3)/2

!$OMP MASTER
      do it=0,np-1
         int_mb(pindex1(1)+it) = 0
         int_mb(pindex2(1)+it) = 0
      end do

      !*** count number of sends ****
      small_count = 0
      do i=1,n2ft3d(1)
         ibig = int_mb(igrid_big(1)+4*(i-1))
         jbig = int_mb(igrid_big(1)+4*(i-1)+1)
         kbig = int_mb(igrid_big(1)+4*(i-1)+2)
         lbig = int_mb(igrid_big(1)+4*(i-1)+3)
         if (lbig.gt.0) then

            ismall = ibig-center(1)
            if (ismall.gt.(nxh-1)) ismall = ismall - nx(1)
            if (ismall.lt.(-nxh))  ismall = ismall + nx(1)

            jsmall = jbig-center(2)
            if (jsmall.gt.(nyh-1)) jsmall = jsmall - ny(1)
            if (jsmall.lt.(-nyh))  jsmall = jsmall + ny(1)

            ksmall = kbig-center(3)
            if (ksmall.gt.(nzh-1)) ksmall = ksmall - nz(1)
            if (ksmall.lt.(-nzh))  ksmall = ksmall + nz(1)

            if ((ismall.gt.(-nxh_small-1)).and.
     >          (ismall.lt.nxh_small).and.
     >          (jsmall.gt.(-nzh_small-1)).and.
     >          (jsmall.lt.nyh_small).and.
     >          (ksmall.gt.(-nzh_small-1)).and.
     >          (ksmall.lt.nzh_small)) then
               ismall = ismall + nxh_small
               jsmall = jsmall + nyh_small
               ksmall = ksmall + nzh_small
               call D3dB_ijktoindex2p(3,ismall+1,jsmall+1,ksmall+1,
     >                             indxsmall,psmall)

               int_mb(small(1)+4*small_count)   = indxsmall
               int_mb(small(1)+4*small_count+1) = psmall
               int_mb(small(1)+4*small_count+2) = i

               int_mb(small(1)+4*small_count+3) =
     >         int_mb(pindex1(1)+psmall)

               small_count = small_count + 1
               int_mb(pindex1(1)+psmall) = int_mb(pindex1(1)+psmall) + 1
            end if
         end if
      end do

      !*** count number of recvs ****
      do i=1,n2ft3d(3)
         ismall = int_mb(igrid_small(1)+4*(i-1))
         jsmall = int_mb(igrid_small(1)+4*(i-1)+1)
         ksmall = int_mb(igrid_small(1)+4*(i-1)+2)
         lsmall = int_mb(igrid_small(1)+4*(i-1)+3)
         if (lsmall.gt.0) then

            ibig = ismall + center(1)
            if (ibig.lt.(-nxh))  ibig = ibig + nx(1)
            if (ibig.gt.(nxh-1)) ibig = ibig - nx(1)
            ibig = ibig + nxh

            jbig = jsmall + center(2)
            if (jbig.lt.(-nyh))  jbig = jbig + ny(1)
            if (jbig.gt.(nyh-1)) jbig = jbig - ny(1)
            jbig = jbig + nyh

            kbig = ksmall + center(3)
            if (kbig.lt.(-nzh))  kbig = kbig + nz(1)
            if (kbig.gt.(nzh-1)) kbig = kbig - nz(1)
            kbig = kbig + nzh

            call D3dB_ijktoindex2p(1,ibig+1,jbig+1,kbig+1,
     >                             indxbig,pbig)
            int_mb(pindex2(1)+pbig) = int_mb(pindex2(1)+pbig) + 1
         end if

      end do

      index1 = 1
      index2 = 1
      do it=0,np-1
         proc_to   = mod(taskid+it,np)
         proc_from = mod(taskid-it+np,np)
         int_mb(i1_start+it) = index1
         int_mb(i2_start+it) = index2
         index1 = index1 + int_mb(pindex1(1)+proc_to)
         index2 = index2 + int_mb(pindex2(1)+proc_from)
      end do
      int_mb(i1_start+np) = index1
      int_mb(i2_start+np) = index2

c      do it=0,np-1
c         int_mb(pindex1(1)+it) = 0
c         int_mb(pindex2(1)+it) = 0
c      end do
c
c      do i=1,n2ft3d(1)
c         ibig = int_mb(igrid_big(1)+4*(i-1))
c         jbig = int_mb(igrid_big(1)+4*(i-1)+1)
c         kbig = int_mb(igrid_big(1)+4*(i-1)+2)
c         lbig = int_mb(igrid_big(1)+4*(i-1)+3)
c         if (lbig.gt.0) then
c
c            ismall = ibig-center(1)
c            if (ismall.gt.(nxh-1)) ismall = ismall - nx(1)
c            if (ismall.lt.(-nxh))  ismall = ismall + nx(1)
c
c            jsmall = jbig-center(2)
c            if (jsmall.gt.(nyh-1)) jsmall = jsmall - ny(1)
c            if (jsmall.lt.(-nyh))  jsmall = jsmall + ny(1)
c
c            ksmall = kbig-center(3)
c            if (ksmall.gt.(nzh-1)) ksmall = ksmall - nz(1)
c            if (ksmall.lt.(-nzh))  ksmall = ksmall + nz(1)
c
c            if ((ismall.gt.(-nxh_small-1)).and.
c     >          (ismall.lt.nxh_small).and.
c     >          (jsmall.gt.(-nyh_small-1)).and.
c     >          (jsmall.lt.nyh_small).and.
c     >          (ksmall.gt.(-nzh_small-1)).and.
c     >          (ksmall.lt.nzh_small)) then
cc               ibig = ibig + nxh
cc               jbig = jbig + nyh
cc               kbig = kbig + nzh
cc               call D3dB_ijktoindex2p(1,ibig+1,jbig+1,kbig+1,
cc     >                             indxbig,pbig)
c               ismall = ismall + nxh_small
c               jsmall = jsmall + nyh_small
c               ksmall = ksmall + nzh_small
c               call D3dB_ijktoindex2p(3,ismall+1,jsmall+1,ksmall+1,
c     >                             indxsmall,psmall)
c
c               it = mod(psmall-taskid+np,np)
c               !it = mod(taskid-psmall+np,np)
c               itmp = int_mb(i1_start+it) 
c     >              + int_mb(pindex1(1)+psmall)
cc               int_mb(c_iq_to_i1+itmp-1) = indxbig
c               int_mb(c_iq_to_i1+itmp-1) = i
c               int_mb(c_iq_to_i2+itmp-1) = indxsmall
c               int_mb(pindex1(1)+psmall) = int_mb(pindex1(1)+psmall)+1
c
c            end if
c         end if
c      end do

      do i=1,small_count
         indxsmall = int_mb(small(1)+4*(i-1))
         psmall    = int_mb(small(1)+4*(i-1)+1)
         indxbig   = int_mb(small(1)+4*(i-1)+2)
         itmp      = int_mb(small(1)+4*(i-1)+3)
         it = mod(psmall-taskid+np,np)
         itmp = itmp + int_mb(i1_start+it)
cc         itmp = int_mb(i1_start+it) 
cc     >              + int_mb(pindex1(1)+psmall)
         int_mb(c_iq_to_i1+itmp-1) = indxbig
         int_mb(c_iq_to_i2+itmp-1) = indxsmall
c         int_mb(pindex1(1)+psmall) = int_mb(pindex1(1)+psmall)+1
      end do

*     **** set index2 and increment qlast_1t3 ****
      int_mb(index2_1t3(1)+qlast_1t3) = index2
      qlast_1t3 = mod(qlast_1t3+1,maxq_1t3)

!$OMP END MASTER


!$OMP DO
      do i=1,index1-1
         tmp1(i) = A(int_mb(c_iq_to_i1+i-1)) 
      end do
!$OMP END DO
      call nwpw_timing_end(55)

      return
      end



*     ***********************************
*     *                                 *
*     *    D3dB_r_1to3_Contract_mid3    *
*     *                                 *
*     ***********************************
      subroutine D3dB_r_1to3_Contract_mid3(tmp1,tmp2,request,reqcnt)
      implicit none
      real*8 tmp1(*)
      real*8 tmp2(*)
      integer request(*),reqcnt


#include "bafdecls.fh"
#include "mpif.h"
#include "D3dB.fh"
#include "errquit.fh"

#ifdef MPI4
#include "stupid_mpi4.fh"
#endif

*     **** r_1to3_block common block ****
      integer maxq_1t3
      parameter (maxq_1t3=6)
      integer qhead0_1t3,qhead1_1t3,qlast_1t3
      integer i1_start_1t3(2),i2_start_1t3(2)
      integer c_iq_to_i1_1t3(2),c_iq_to_i2_1t3(2)
      integer index2_1t3(2)
      common /r_1to3_block/ i1_start_1t3,i2_start_1t3,
     >                      c_iq_to_i1_1t3,c_iq_to_i2_1t3,
     >                      index2_1t3,
     >                      qhead0_1t3,qhead1_1t3,qlast_1t3

*     **** local variables ****
      logical value
      integer taskid,np
      integer proc_to,proc_from,source
      integer it,itmp,msgtype,msglen
      integer i1_start,i2_start,c_iq_to_i1,c_iq_to_i2
      integer shift1,shift2
      integer mpierr

c     **** external functions ****
      integer  Parallel2d_convert_taskid_i,Parallel_comm_world
      external Parallel2d_convert_taskid_i,Parallel_comm_world


      call nwpw_timing_start(55)
      call Parallel2d_taskid_i(taskid)
      call Parallel2d_np_i(np)

      shift1 = qhead0_1t3*(np+1)
      shift2 = qhead0_1t3*n2ft3d(1)
      i1_start = i1_start_1t3(1) + shift1
      i2_start = i2_start_1t3(1) + shift1
      c_iq_to_i1 = c_iq_to_i1_1t3(1) + shift2
      c_iq_to_i2 = c_iq_to_i2_1t3(1) + shift2


*     **** it = 0, transpose data on same thread ****  
      msglen = int_mb(i2_start+2-1) - int_mb(i2_start+1-1)

      call Parallel_shared_vector_copy(.false.,msglen,
     >               tmp1(int_mb(i1_start+1-1)),
     >               tmp2(int_mb(i2_start+1-1)))

      call Parallel_shared_vector_icopy(.false.,msglen,
     >               int_mb(c_iq_to_i2+int_mb(i1_start+1-1)-1),
     >               int_mb(c_iq_to_i1+int_mb(i2_start+1-1)-1))

!$OMP MASTER

*     **** receive packed array data ****
      reqcnt = 0
      do it = 1,np-1

*        **** synchronous receive of tmp ****
         proc_from = mod((taskid)-it+np, np)

         source=proc_from
         msglen = (int_mb(i2_start+it+2-1)
     >          -  int_mb(i2_start+it+1-1))

         if (msglen.gt.0) then
            msgtype = 73
            reqcnt = reqcnt + 1
#ifdef MPI4
            stupid_msglen = msglen
            stupid_type   = msgtype
            stupid_taskid = Parallel2d_convert_taskid_i(source)
            call MPI_IRECV(int_mb(c_iq_to_i1+int_mb(i2_start+it+1-1)-1),
     >                    2*stupid_msglen,stupid_integer,
     >                    stupid_taskid,
     >                    stupid_type,stupid_world,
     >                    stupid_request,stupid_ierr)
            request(reqcnt) = 0
            request(reqcnt) = stupid_request
#else
            call MPI_IRECV(int_mb(c_iq_to_i1+int_mb(i2_start+it+1-1)-1),
     >                    msglen,MPI_INTEGER,
     >                    Parallel2d_convert_taskid_i(source),
     >                    msgtype,Parallel_comm_world(),
     >                    request(reqcnt),mpierr)
#endif

            msgtype = 74
            reqcnt = reqcnt + 1
            !write(*,*) taskid," recv2=",source," msglen=",msglen
#ifdef MPI4
            stupid_msglen = msglen
            stupid_type   = msgtype
            stupid_taskid = Parallel2d_convert_taskid_i(source)
            call MPI_IRECV(tmp2(int_mb(i2_start+it+1-1)),
     >                    stupid_msglen,stupid_double,
     >                    stupid_taskid,
     >                    stupid_type,stupid_world,
     >                    stupid_request,stupid_ierr)
            request(reqcnt) = 0
            request(reqcnt) = stupid_request
#else
            call MPI_IRECV(tmp2(int_mb(i2_start+it+1-1)),
     >                    msglen,MPI_DOUBLE_PRECISION,
     >                    Parallel2d_convert_taskid_i(source),
     >                    msgtype,Parallel_comm_world(),
     >                    request(reqcnt),mpierr)
#endif
         end if

      end do
*     **** send packed array to other processors ****
      do it = 1,np-1

*        **** synchronous send of tmp ****
         proc_to   = mod((taskid)+it, np)
         msglen    = (int_mb(i1_start+it+2-1)
     >              - int_mb(i1_start+it+1-1))

         if (msglen.gt.0) then
            msgtype   = 73
            reqcnt = reqcnt + 1
#ifdef MPI4
            stupid_msglen = msglen
            stupid_type   = msgtype
            stupid_taskid = Parallel2d_convert_taskid_i(proc_to)
            call MPI_ISEND(int_mb(c_iq_to_i2+int_mb(i1_start+it+1-1)-1),
     >                     2*stupid_msglen,stupid_integer,
     >                     stupid_taskid,
     >                     stupid_type,stupid_world,
     >                     stupid_request,stupid_ierr)
            request(reqcnt) = 0
            request(reqcnt) = stupid_request
#else
            call MPI_ISEND(int_mb(c_iq_to_i2+int_mb(i1_start+it+1-1)-1),
     >                     2*msglen,MPI_INTEGER,
     >                     Parallel2d_convert_taskid_i(proc_to),
     >                     msgtype,Parallel_comm_world(),
     >                     request(reqcnt),mpierr)
#endif
            msgtype   = 74
            reqcnt = reqcnt + 1
            !write(*,*) taskid," send2=",proc_to," msglen=",msglen
#ifdef MPI4
            stupid_msglen = msglen
            stupid_type   = msgtype
            stupid_taskid = Parallel2d_convert_taskid_i(proc_to)
            call MPI_ISEND(tmp1(int_mb(i1_start+it+1-1)),
     >                     stupid_msglen,stupid_double,
     >                     stupid_taskid,
     >                     stupid_type,stupid_world,
     >                     stupid_request,stupid_ierr)
            request(reqcnt) = 0
            request(reqcnt) = stupid_request
#else
            call MPI_ISEND(tmp1(int_mb(i1_start+it+1-1)),
     >                     msglen,MPI_DOUBLE_PRECISION,
     >                     Parallel2d_convert_taskid_i(proc_to),
     >                     msgtype,Parallel_comm_world(),
     >                     request(reqcnt),mpierr)
#endif
         end if

      end do

*     *** increment qhead0_1t3
      qhead0_1t3 = mod(qhead0_1t3+1,maxq_1t3)


!$OMP END MASTER
!$OMP BARRIER
      call nwpw_timing_end(55)


      return
      end


*     ***********************************
*     *                                 *
*     *    D3dB_r_1to3_Contract_final3  *
*     *                                 *
*     ***********************************
      subroutine D3dB_r_1to3_Contract_final3(tmp2,A,request,reqcnt)
      implicit none
      real*8 tmp2(*),A(*)
      integer request(*),reqcnt

#include "bafdecls.fh"
#include "D3dB.fh"
#include "errquit.fh"

*     **** r_1to3_block common block ****
      integer maxq_1t3
      parameter (maxq_1t3=6)
      integer qhead0_1t3,qhead1_1t3,qlast_1t3
      integer i1_start_1t3(2),i2_start_1t3(2)
      integer c_iq_to_i1_1t3(2),c_iq_to_i2_1t3(2)
      integer index2_1t3(2)
      common /r_1to3_block/ i1_start_1t3,i2_start_1t3,
     >                      c_iq_to_i1_1t3,c_iq_to_i2_1t3,
     >                      index2_1t3,
     >                      qhead0_1t3,qhead1_1t3,qlast_1t3

*     **** local variables ****
      logical value
      integer taskid,np
      integer proc_to,proc_from,source
      integer itmp,msgtype,msglen,it,i
      integer i1_start,i2_start,c_iq_to_i1,c_iq_to_i2
      integer shift1,shift2,index2
      integer indxsmall,psmall
      integer ismall,jsmall,ksmall,lsmall

      call nwpw_timing_start(55)
      call Parallel2d_taskid_i(taskid)
      call Parallel2d_np_i(np)

      shift1 = qhead1_1t3*(np+1)
      shift2 = qhead1_1t3*n2ft3d(1)
      i1_start = i1_start_1t3(1) + shift1
      i2_start = i2_start_1t3(1) + shift1
      c_iq_to_i1 = c_iq_to_i1_1t3(1) + shift2
      c_iq_to_i2 = c_iq_to_i2_1t3(1) + shift2
      index2 = int_mb(index2_1t3(1)+qhead1_1t3)

!$OMP MASTER
*     **** wait for completion of mp_send, also do a sync ****
      if (np.gt.1) call Parallel_mpiWaitAll0(reqcnt,request)

*     *** increment qhead1_1t3
      qhead1_1t3 = mod(qhead1_1t3+1,maxq_1t3)

!$OMP END MASTER
!$OMP BARRIER


*     **** unpack tmp2(i) array ****
      call Parallel_shared_vector_zero(.true.,n2ft3d(3),A)
!$OMP DO
      do i=1,index2-1
         A(int_mb(c_iq_to_i1+i-1)) = tmp2(i)
      end do
!$OMP END DO
      call nwpw_timing_end(55)

      return
      end








*     ***********************************
*     *                                 *
*     *    D3dB_r_1to3_Contract_start   *
*     *                                 *
*     ***********************************
      subroutine D3dB_r_1to3_Contract_start(center,A,tmp1,
     >                                      request,reqcnt)
      implicit none
      integer center(3)
      real*8 A(*)
      real*8 tmp1(*)
      integer request(*),reqcnt

#include "bafdecls.fh"
#include "mpif.h"
#include "D3dB.fh"
#include "errquit.fh"

#ifdef MPI4
#include "stupid_mpi4.fh"
#endif

*     **** r_1to3_block common block ****
      integer maxq_1t3
      parameter (maxq_1t3=6)
      integer qhead0_1t3,qhead1_1t3,qlast_1t3
      integer i1_start_1t3(2),i2_start_1t3(2)
      integer c_iq_to_i1_1t3(2),c_iq_to_i2_1t3(2)
      integer index2_1t3(2)
      common /r_1to3_block/ i1_start_1t3,i2_start_1t3,
     >                      c_iq_to_i1_1t3,c_iq_to_i2_1t3,
     >                      index2_1t3,
     >                      qhead0_1t3,qhead1_1t3,qlast_1t3

*     **** local variables ****
      logical value
      integer taskid,np
      integer proc_to,proc_from,source
      integer indxbig,pbig,indxsmall,psmall,itmp
      integer index1,index2,msgtype,msglen
      integer k1,k2,k3,l1,l2,l3,i,j,k,il,jl,kl,it
      integer nxh,nyh,nzh,nxh_small,nyh_small,nzh_small
      integer pindex1(2),pindex2(2)
      integer i1_start,i2_start,c_iq_to_i1,c_iq_to_i2
      integer shift1,shift2
      

      call nwpw_timing_start(55)
      call Parallel2d_taskid_i(taskid)
      call Parallel2d_np_i(np)

      !**** allocate temporary space from stack ****
      value = BA_push_get(mt_int,np,'pindex1_1t3',pindex1(2),pindex1(1))
      value = value.and.
     >        BA_push_get(mt_int,np,'pindex2_1t3',pindex2(2),pindex2(1))
      if (.not.value)
     >  call errquit('D3dB_r_1to3_Contract_start:out of stack',0,MA_ERR)


      shift1 = qlast_1t3*(np+1)
      shift2 = qlast_1t3*n2ft3d(1)
      i1_start = i1_start_1t3(1) + shift1
      i2_start = i2_start_1t3(1) + shift1
      c_iq_to_i1 = c_iq_to_i1_1t3(1) + shift2
      c_iq_to_i2 = c_iq_to_i2_1t3(1) + shift2

      nxh = nx(1)/2
      nyh = ny(1)/2
      nzh = nz(1)/2

      nxh_small = nx(3)/2
      nyh_small = ny(3)/2
      nzh_small = nz(3)/2

!$OMP MASTER
      do it=0,np-1
         int_mb(pindex1(1)+it) = 0
         int_mb(pindex2(1)+it) = 0
      end do
      do k3=-nzh_small,nzh_small-1
         l3 = k3 + center(3)
         if (l3.lt.(-nzh))  l3 = l3 + nz(1)
         if (l3.gt.(nzh-1)) l3 = l3 - nz(1)
         k  = k3 + nzh_small
         kl = l3 + nzh

         do k2=-nyh_small,nyh_small-1
            l2 = k2 + center(2)
            if (l2.lt.(-nyh))  l2 = l2 + ny(1)
            if (l2.gt.(nyh-1)) l2 = l2 - ny(1)
            j  = k2 + nyh_small
            jl = l2 + nyh

            do k1=-nxh_small,nxh_small-1
               l1 = k1 + center(1)
               if (l1.lt.(-nxh))  l1 = l1 + nx(1)
               if (l1.gt.(nxh-1)) l1 = l1 - nx(1)
               i  = k1 + nxh_small
               il = l1 + nxh

               call D3dB_ijktoindex2p(1,il+1,jl+1,kl+1,indxbig,pbig)
               call D3dB_ijktoindex2p(3,i+1,j+1,k+1,indxsmall,psmall)


*              **** packing scheme **** 
               if (pbig.eq.taskid) then
                  int_mb(pindex1(1)+psmall)=int_mb(pindex1(1)+psmall)+1
               end if

*              **** unpacking scheme **** 
               if (psmall.eq.taskid) then
                  int_mb(pindex2(1)+pbig) = int_mb(pindex2(1)+pbig)+1
               end if

            end do
         end do
      end do
      index1 = 1
      index2 = 1
      do it=0,np-1
         proc_to   = mod(taskid+it,np)
         proc_from = mod(taskid-it+np,np)
         int_mb(i1_start+it) = index1
         int_mb(i2_start+it) = index2
         index1 = index1 + int_mb(pindex1(1)+proc_to)
         index2 = index2 + int_mb(pindex2(1)+proc_from)
      end do
      int_mb(i1_start+np) = index1
      int_mb(i2_start+np) = index2


      do it=0,np-1
         int_mb(pindex1(1)+it) = 0
         int_mb(pindex2(1)+it) = 0
      end do
      do k3=-nzh_small,nzh_small-1
         l3 = k3 + center(3)
         if (l3.lt.(-nzh))  l3 = l3 + nz(1)
         if (l3.gt.(nzh-1)) l3 = l3 - nz(1)
         k  = k3 + nzh_small
         kl = l3 + nzh

         do k2=-nyh_small,nyh_small-1
            l2 = k2 + center(2)
            if (l2.lt.(-nyh))  l2 = l2 + ny(1)
            if (l2.gt.(nyh-1)) l2 = l2 - ny(1)
            j  = k2 + nyh_small
            jl = l2 + nyh

            do k1=-nxh_small,nxh_small-1
               l1 = k1 + center(1)
               if (l1.lt.(-nxh))  l1 = l1 + nx(1)
               if (l1.gt.(nxh-1)) l1 = l1 - nx(1)
               i  = k1 + nxh_small
               il = l1 + nxh

               call D3dB_ijktoindex2p(1,il+1,jl+1,kl+1,indxbig,pbig)
               call D3dB_ijktoindex2p(3,i+1,j+1,k+1,indxsmall,psmall)


*              **** packing scheme **** 
               if (pbig.eq.taskid) then
                  it = mod(psmall-taskid+np,np)
                  itmp = int_mb(i1_start+it) 
     >                 + int_mb(pindex1(1)+psmall)
                  int_mb(c_iq_to_i1+itmp-1) = indxbig
                  int_mb(pindex1(1)+psmall) =int_mb(pindex1(1)+psmall)+1
               end if

*              **** unpacking scheme **** 
               if (psmall.eq.taskid) then
                  it = mod(taskid-pbig+np,np)
                  itmp = int_mb(i2_start+it) 
     >                 + int_mb(pindex2(1)+pbig)
                  int_mb(c_iq_to_i2+itmp-1) = indxsmall
                  int_mb(pindex2(1)+pbig)   = int_mb(pindex2(1)+pbig)+1
               end if

            end do
         end do
      end do

*     **** pop stack ****
      value =           BA_pop_stack(pindex2(2))
      value = value.and.BA_pop_stack(pindex1(2))
      if (.not.value)
     >   call errquit('D3dB_r_1to3_Contract_start:popping stack',
     >                0,MA_ERR)


*     **** set index2 and increment qlast_1t3 ****
      int_mb(index2_1t3(1)+qlast_1t3) = index2
      qlast_1t3 = mod(qlast_1t3+1,maxq_1t3)

!$OMP END MASTER
!$OMP BARRIER


*     **** pack A(i) array ****
!$OMP DO
      do i=1,index1-1
         tmp1(i) = A(int_mb(c_iq_to_i1+i-1)) 
      end do
!$OMP END DO
      call nwpw_timing_end(55)

      return
      end


*     ***********************************
*     *                                 *
*     *    D3dB_r_1to3_Contract_mid     *
*     *                                 *
*     ***********************************
      subroutine D3dB_r_1to3_Contract_mid(tmp1,tmp2,request,reqcnt)
      implicit none
      real*8 tmp1(*)
      real*8 tmp2(*)
      integer request(*),reqcnt
      

#include "bafdecls.fh"
#include "mpif.h"
#include "D3dB.fh"
#include "errquit.fh"

#ifdef MPI4
#include "stupid_mpi4.fh"
#endif

*     **** r_1to3_block common block ****
      integer maxq_1t3
      parameter (maxq_1t3=6)
      integer qhead0_1t3,qhead1_1t3,qlast_1t3
      integer i1_start_1t3(2),i2_start_1t3(2)
      integer c_iq_to_i1_1t3(2),c_iq_to_i2_1t3(2)
      integer index2_1t3(2)
      common /r_1to3_block/ i1_start_1t3,i2_start_1t3,
     >                      c_iq_to_i1_1t3,c_iq_to_i2_1t3,
     >                      index2_1t3,
     >                      qhead0_1t3,qhead1_1t3,qlast_1t3

*     **** local variables ****
      logical value
      integer taskid,np
      integer proc_to,proc_from,source
      integer it,itmp,msgtype,msglen
      integer i1_start,i2_start,c_iq_to_i1,c_iq_to_i2
      integer shift1,shift2
      integer mpierr

c     **** external functions ****
      integer  Parallel2d_convert_taskid_i,Parallel_comm_world
      external Parallel2d_convert_taskid_i,Parallel_comm_world


      call nwpw_timing_start(55)
      call Parallel2d_taskid_i(taskid)
      call Parallel2d_np_i(np)

      shift1 = qhead0_1t3*(np+1)
      shift2 = qhead0_1t3*n2ft3d(1)
      i1_start = i1_start_1t3(1) + shift1
      i2_start = i2_start_1t3(1) + shift1
      c_iq_to_i1 = c_iq_to_i1_1t3(1) + shift2
      c_iq_to_i2 = c_iq_to_i2_1t3(1) + shift2


*     **** it = 0, transpose data on same thread ****  
      msglen = int_mb(i2_start+2-1) - int_mb(i2_start+1-1)

      call Parallel_shared_vector_copy(.false.,msglen,
     >               tmp1(int_mb(i1_start+1-1)),
     >               tmp2(int_mb(i2_start+1-1)))

!$OMP MASTER

*     **** receive packed array data ****
      reqcnt = 0
      do it = 1,np-1

*        **** synchronous receive of tmp ****
         proc_from = mod((taskid)-it+np, np)

         msgtype = 73
         source=proc_from
         msglen = (int_mb(i2_start+it+2-1)
     >          -  int_mb(i2_start+it+1-1))

         if (msglen.gt.0) then
            reqcnt = reqcnt + 1
#ifdef MPI4
            stupid_msglen = msglen
            stupid_type   = msgtype
            stupid_taskid = Parallel2d_convert_taskid_i(source)
            call MPI_IRECV(tmp2(int_mb(i2_start+it+1-1)),
     >                    stupid_msglen,stupid_double,
     >                    stupid_taskid,
     >                    stupid_type,stupid_world,
     >                    stupid_request,stupid_ierr)
            request(reqcnt) = 0
            request(reqcnt) = stupid_request
#else
            call MPI_IRECV(tmp2(int_mb(i2_start+it+1-1)),
     >                    msglen,MPI_DOUBLE_PRECISION,
     >                    Parallel2d_convert_taskid_i(source),
     >                    msgtype,Parallel_comm_world(),
     >                    request(reqcnt),mpierr)
#endif
         end if

      end do

*     **** send packed array to other processors ****
      do it = 1,np-1

*        **** synchronous send of tmp ****
         proc_to   = mod((taskid)+it, np)
         msglen    = (int_mb(i1_start+it+2-1)
     >              - int_mb(i1_start+it+1-1))
         msgtype   = 73

         if (msglen.gt.0) then
            reqcnt = reqcnt + 1
#ifdef MPI4
            stupid_msglen = msglen
            stupid_type   = msgtype
            stupid_taskid = Parallel2d_convert_taskid_i(proc_to)
            call MPI_ISEND(tmp1(int_mb(i1_start+it+1-1)),
     >                     stupid_msglen,stupid_double,
     >                     stupid_taskid,
     >                     stupid_type,stupid_world,
     >                     stupid_request,stupid_ierr)
            request(reqcnt) = 0
            request(reqcnt) = stupid_request
#else
            call MPI_ISEND(tmp1(int_mb(i1_start+it+1-1)),
     >                     msglen,MPI_DOUBLE_PRECISION,
     >                     Parallel2d_convert_taskid_i(proc_to),
     >                     msgtype,Parallel_comm_world(),
     >                     request(reqcnt),mpierr)
#endif
         end if

      end do

*     *** increment qhead0_1t3
      qhead0_1t3 = mod(qhead0_1t3+1,maxq_1t3)


!$OMP END MASTER
!$OMP BARRIER
      call nwpw_timing_end(55)


      return
      end


*     ***********************************
*     *                                 *
*     *    D3dB_r_1to3_Contract_final   *
*     *                                 *
*     ***********************************
      subroutine D3dB_r_1to3_Contract_final(tmp2,A,request,reqcnt)
      implicit none
      real*8 tmp2(*),A(*)
      integer request(*),reqcnt

#include "bafdecls.fh"
#include "D3dB.fh"
#include "errquit.fh"

*     **** r_1to3_block common block ****
      integer maxq_1t3
      parameter (maxq_1t3=6)
      integer qhead0_1t3,qhead1_1t3,qlast_1t3
      integer i1_start_1t3(2),i2_start_1t3(2)
      integer c_iq_to_i1_1t3(2),c_iq_to_i2_1t3(2)
      integer index2_1t3(2)
      common /r_1to3_block/ i1_start_1t3,i2_start_1t3,
     >                      c_iq_to_i1_1t3,c_iq_to_i2_1t3,
     >                      index2_1t3,
     >                      qhead0_1t3,qhead1_1t3,qlast_1t3

*     **** local variables ****
      logical value
      integer taskid,np
      integer proc_to,proc_from,source
      integer itmp,msgtype,msglen,it,i
      integer i1_start,i2_start,c_iq_to_i1,c_iq_to_i2
      integer shift1,shift2,index2
      integer indxsmall,psmall


      call nwpw_timing_start(55)
      call Parallel2d_taskid_i(taskid)
      call Parallel2d_np_i(np)

      
      shift1 = qhead1_1t3*(np+1)
      shift2 = qhead1_1t3*n2ft3d(1)
      i1_start = i1_start_1t3(1) + shift1
      i2_start = i2_start_1t3(1) + shift1
      c_iq_to_i1 = c_iq_to_i1_1t3(1) + shift2
      c_iq_to_i2 = c_iq_to_i2_1t3(1) + shift2
      index2 = int_mb(index2_1t3(1)+qhead1_1t3)

!$OMP MASTER
*     **** wait for completion of mp_send, also do a sync ****
      if (np.gt.1) call Parallel_mpiWaitAll0(reqcnt,request)

*     *** increment qhead1_1t3
      qhead1_1t3 = mod(qhead1_1t3+1,maxq_1t3)

!$OMP END MASTER
!$OMP BARRIER


*     **** unpack tmp2(i) array ****
      call Parallel_shared_vector_zero(.true.,n2ft3d(3),A)
!$OMP DO
      do i=1,index2-1
         A(int_mb(c_iq_to_i2+i-1)) = tmp2(i)
      end do
!$OMP END DO
      call nwpw_timing_end(55)

      return
      end



*     ***********************************
*     *                                 *
*     *         D3dB_r_3to1_init        *
*     *                                 *
*     ***********************************
      subroutine D3dB_r_3to1_init()
      implicit none

#include "bafdecls.fh"
#include "D3dB.fh"
#include "errquit.fh"

*     **** r_3to1_block common block ****
      integer maxq_3t1
      parameter (maxq_3t1=6)
      integer qhead0_3t1,qhead1_3t1,qlast_3t1
      integer i1_start_3t1(2),i2_start_3t1(2)
      integer c_iq_to_i1_3t1(2),c_iq_to_i2_3t1(2)
      integer index2_3t1(2)
      common /r_3to1_block/ i1_start_3t1,i2_start_3t1,
     >                      c_iq_to_i1_3t1,c_iq_to_i2_3t1,
     >                      index2_3t1,
     >                      qhead0_3t1,qhead1_3t1,qlast_3t1

*     **** local variables ****
      logical value
      integer np


      call Parallel2d_np_i(np)

      !**** allocate indexing space from heap ****
      value = BA_alloc_get(mt_int,maxq_3t1*(np+1),'i1_start_3t1',
     >                    i1_start_3t1(2),i1_start_3t1(1))
      value = value.and.
     >        BA_alloc_get(mt_int,maxq_3t1*(np+1),'i2_start_3t1',
     >                    i2_start_3t1(2),i2_start_3t1(1))
      value = value.and.
     >        BA_alloc_get(mt_int,maxq_3t1*n2ft3d(1),'c_iq_to_i1_3t1',
     >                    c_iq_to_i1_3t1(2),c_iq_to_i1_3t1(1))
      value = value.and.
     >        BA_alloc_get(mt_int,maxq_3t1*n2ft3d(1),'c_iq_to_i2_3t1',
     >                    c_iq_to_i2_3t1(2),c_iq_to_i2_3t1(1))
      value = value.and.
     >        BA_alloc_get(mt_int,maxq_3t1,'index_3t1',
     >                     index2_3t1(2),index2_3t1(1))
      if (.not.value)
     >   call errquit('D3dB_r_3to1_init:out of heap',0,MA_ERR)

!$OMP MASTER
      qhead0_3t1 = 0
      qhead1_3t1 = 0
      qlast_3t1  = 0
!$OMP END MASTER
!$OMP BARRIER

      return
      end 



*     ***********************************
*     *                                 *
*     *         D3dB_r_3to1_end         *
*     *                                 *
*     ***********************************
      subroutine D3dB_r_3to1_end()
      implicit none

#include "bafdecls.fh"
#include "D3dB.fh"
#include "errquit.fh"

*     **** r_3to1_block common block ****
      integer maxq_3t1
      parameter (maxq_3t1=6)
      integer qhead0_3t1,qhead1_3t1,qlast_3t1
      integer i1_start_3t1(2),i2_start_3t1(2)
      integer c_iq_to_i1_3t1(2),c_iq_to_i2_3t1(2)
      integer index2_3t1(2)
      common /r_3to1_block/ i1_start_3t1,i2_start_3t1,
     >                      c_iq_to_i1_3t1,c_iq_to_i2_3t1,
     >                      index2_3t1,
     >                      qhead0_3t1,qhead1_3t1,qlast_3t1

*     **** local variables ****
      logical value

      !**** freeing indexing space from heap ****
      value =           BA_free_heap(i1_start_3t1(2))
      value = value.and.BA_free_heap(i2_start_3t1(2))
      value = value.and.BA_free_heap(c_iq_to_i1_3t1(2))
      value = value.and.BA_free_heap(c_iq_to_i2_3t1(2))
      value = value.and.BA_free_heap(index2_3t1(2))
      if (.not.value)
     >   call errquit('D3dB_r_3to1_end:freeing heap',0,MA_ERR)

      return
      end





*     ***********************************
*     *                                 *
*     *    D3dB_r_3to1_Expand_start     *
*     *                                 *
*     ***********************************
      subroutine D3dB_r_3to1_Expand_start(center,A,tmp1,
     >                                      request,reqcnt)
      implicit none
      integer center(3)
      real*8 A(*)
      real*8 tmp1(*)
      integer request(*),reqcnt

#include "bafdecls.fh"
#include "mpif.h"
#include "D3dB.fh"
#include "errquit.fh"

#ifdef MPI4
#include "stupid_mpi4.fh"
#endif

*     **** r_3to1_block common block ****
      integer maxq_3t1
      parameter (maxq_3t1=6)
      integer qhead0_3t1,qhead1_3t1,qlast_3t1
      integer i1_start_3t1(2),i2_start_3t1(2)
      integer c_iq_to_i1_3t1(2),c_iq_to_i2_3t1(2)
      integer index2_3t1(2)
      common /r_3to1_block/ i1_start_3t1,i2_start_3t1,
     >                      c_iq_to_i1_3t1,c_iq_to_i2_3t1,
     >                      index2_3t1,
     >                      qhead0_3t1,qhead1_3t1,qlast_3t1

*     **** local variables ****
      logical value
      integer taskid,np
      integer proc_to,proc_from,source
      integer indxbig,pbig,indxsmall,psmall,itmp
      integer index1,index2,msgtype,msglen
      integer k1,k2,k3,l1,l2,l3,i,j,k,il,jl,kl,it
      integer nxh,nyh,nzh,nxh_small,nyh_small,nzh_small
      integer pindex1(2),pindex2(2)
      integer i1_start,i2_start,c_iq_to_i1,c_iq_to_i2
      integer shift1,shift2
      

      call nwpw_timing_start(55)
      call Parallel2d_taskid_i(taskid)
      call Parallel2d_np_i(np)

      !**** allocate temporary space from stack ****
      value = BA_push_get(mt_int,np,'pindex1_3t1',pindex1(2),pindex1(1))
      value = value.and.
     >        BA_push_get(mt_int,np,'pindex2_3t1',pindex2(2),pindex2(1))
      if (.not.value)
     >  call errquit('D3dB_r_3to1_Expand_start:out of stack',0,MA_ERR)

      shift1 = qlast_3t1*(np+1)
      shift2 = qlast_3t1*n2ft3d(1)
      i1_start = i1_start_3t1(1) + shift1
      i2_start = i2_start_3t1(1) + shift1
      c_iq_to_i1 = c_iq_to_i1_3t1(1) + shift2
      c_iq_to_i2 = c_iq_to_i2_3t1(1) + shift2

      nxh = nx(1)/2
      nyh = ny(1)/2
      nzh = nz(1)/2

      nxh_small = nx(3)/2
      nyh_small = ny(3)/2
      nzh_small = nz(3)/2

!$OMP MASTER
      do it=0,np-1
         int_mb(pindex1(1)+it) = 0
         int_mb(pindex2(1)+it) = 0
      end do
      do k3=-nzh_small,nzh_small-1
         l3 = k3 + center(3)
         if (l3.lt.(-nzh))  l3 = l3 + nz(1)
         if (l3.gt.(nzh-1)) l3 = l3 - nz(1)
         k  = k3 + nzh_small
         kl = l3 + nzh

         do k2=-nyh_small,nyh_small-1
            l2 = k2 + center(2)
            if (l2.lt.(-nyh))  l2 = l2 + ny(1)
            if (l2.gt.(nyh-1)) l2 = l2 - ny(1)
            j  = k2 + nyh_small
            jl = l2 + nyh

            do k1=-nxh_small,nxh_small-1
               l1 = k1 + center(1)
               if (l1.lt.(-nxh))  l1 = l1 + nx(1)
               if (l1.gt.(nxh-1)) l1 = l1 - nx(1)
               i  = k1 + nxh_small
               il = l1 + nxh

               call D3dB_ijktoindex2p(1,il+1,jl+1,kl+1,indxbig,pbig)
               call D3dB_ijktoindex2p(3,i+1,j+1,k+1,indxsmall,psmall)

*              **** packing scheme **** 
               if (psmall.eq.taskid) then
                  int_mb(pindex1(1)+pbig) = int_mb(pindex1(1)+pbig)+1
               end if

*              **** unpacking scheme **** 
               if (pbig.eq.taskid) then
                  int_mb(pindex2(1)+psmall)=int_mb(pindex2(1)+psmall)+1
               end if


            end do
         end do
      end do
      index1 = 1
      index2 = 1
      do it=0,np-1
         proc_to   = mod(taskid+it,np)
         proc_from = mod(taskid-it+np,np)
         int_mb(i1_start+it) = index1
         int_mb(i2_start+it) = index2
         index1 = index1 + int_mb(pindex1(1)+proc_to)
         index2 = index2 + int_mb(pindex2(1)+proc_from)
      end do
      int_mb(i1_start+np) = index1
      int_mb(i2_start+np) = index2

      do it=0,np-1
         int_mb(pindex1(1)+it) = 0
         int_mb(pindex2(1)+it) = 0
      end do
      do k3=-nzh_small,nzh_small-1
         l3 = k3 + center(3)
         if (l3.lt.(-nzh))  l3 = l3 + nz(1)
         if (l3.gt.(nzh-1)) l3 = l3 - nz(1)
         k  = k3 + nzh_small
         kl = l3 + nzh

         do k2=-nyh_small,nyh_small-1
            l2 = k2 + center(2)
            if (l2.lt.(-nyh))  l2 = l2 + ny(1)
            if (l2.gt.(nyh-1)) l2 = l2 - ny(1)
            j  = k2 + nyh_small
            jl = l2 + nyh

            do k1=-nxh_small,nxh_small-1
               l1 = k1 + center(1)
               if (l1.lt.(-nxh))  l1 = l1 + nx(1)
               if (l1.gt.(nxh-1)) l1 = l1 - nx(1)
               i  = k1 + nxh_small
               il = l1 + nxh

               call D3dB_ijktoindex2p(1,il+1,jl+1,kl+1,indxbig,pbig)
               call D3dB_ijktoindex2p(3,i+1,j+1,k+1,indxsmall,psmall)

*              **** packing scheme **** 
               if (psmall.eq.taskid) then
                  it = mod(pbig-taskid+np,np)
                  itmp = int_mb(i1_start+it)
     >                 + int_mb(pindex1(1)+pbig)
                  int_mb(c_iq_to_i1+itmp-1) = indxsmall
                  int_mb(pindex1(1)+pbig) = int_mb(pindex1(1)+pbig)+1
               end if

*              **** unpacking scheme **** 
               if (pbig.eq.taskid) then
                  it = mod(taskid-psmall+np,np)
                  itmp = int_mb(i2_start+it)
     >                 + int_mb(pindex2(1)+psmall)
                  int_mb(c_iq_to_i2+itmp-1) = indxbig
                  int_mb(pindex2(1)+psmall)=int_mb(pindex2(1)+psmall)+1
               end if


            end do
         end do
      end do

*     **** pop stack ****
      value =           BA_pop_stack(pindex2(2))
      value = value.and.BA_pop_stack(pindex1(2))
      if (.not.value)
     >   call errquit('D3dB_r_1to3_Contract_start:popping stack',
     >                0,MA_ERR)


*     **** set index2 and increment qlast_3t1 ****
      int_mb(index2_3t1(1)+qlast_3t1) = index2
      qlast_3t1 = mod(qlast_3t1+1,maxq_3t1)

!$OMP END MASTER
!$OMP BARRIER


*     **** pack A(i) array ****
!$OMP DO
      do i=1,index1-1
         tmp1(i) = A(int_mb(c_iq_to_i1+i-1)) 
      end do
!$OMP END DO

      call nwpw_timing_end(55)

      return
      end


*     ***********************************
*     *                                 *
*     *    D3dB_r_3to1_Expand_mid       *
*     *                                 *
*     ***********************************
      subroutine D3dB_r_3to1_Expand_mid(tmp1,tmp2,request,reqcnt)
      implicit none
      real*8 tmp1(*)
      real*8 tmp2(*)
      integer request(*),reqcnt
      

#include "bafdecls.fh"
#include "mpif.h"
#include "D3dB.fh"
#include "errquit.fh"

#ifdef MPI4
#include "stupid_mpi4.fh"
#endif

*     **** r_3to1_block common block ****
      integer maxq_3t1
      parameter (maxq_3t1=6)
      integer qhead0_3t1,qhead1_3t1,qlast_3t1
      integer i1_start_3t1(2),i2_start_3t1(2)
      integer c_iq_to_i1_3t1(2),c_iq_to_i2_3t1(2)
      integer index2_3t1(2)
      common /r_3to1_block/ i1_start_3t1,i2_start_3t1,
     >                      c_iq_to_i1_3t1,c_iq_to_i2_3t1,
     >                      index2_3t1,
     >                      qhead0_3t1,qhead1_3t1,qlast_3t1

*     **** local variables ****
      logical value
      integer taskid,np
      integer proc_to,proc_from,source
      integer it,itmp,msgtype,msglen
      integer i1_start,i2_start,c_iq_to_i1,c_iq_to_i2
      integer shift1,shift2
      integer mpierr

c     **** external functions ****
      integer  Parallel2d_convert_taskid_i,Parallel_comm_world
      external Parallel2d_convert_taskid_i,Parallel_comm_world


      call nwpw_timing_start(55)
      call Parallel2d_taskid_i(taskid)
      call Parallel2d_np_i(np)

      shift1 = qhead0_3t1*(np+1)
      shift2 = qhead0_3t1*n2ft3d(1)
      i1_start = i1_start_3t1(1) + shift1
      i2_start = i2_start_3t1(1) + shift1
      c_iq_to_i1 = c_iq_to_i1_3t1(1) + shift2
      c_iq_to_i2 = c_iq_to_i2_3t1(1) + shift2


*     **** it = 0, transpose data on same thread ****  
      msglen = int_mb(i2_start+2-1) - int_mb(i2_start+1-1)

      call Parallel_shared_vector_copy(.false.,msglen,
     >               tmp1(int_mb(i1_start+1-1)),
     >               tmp2(int_mb(i2_start+1-1)))



!$OMP MASTER

*     **** receive packed array data ****
      reqcnt = 0
      do it = 1,np-1

*        **** synchronous receive of tmp ****
         proc_from = mod((taskid)-it+np, np)

         msgtype = 54
         source=proc_from
         msglen = (int_mb(i2_start+it+2-1)
     >          -  int_mb(i2_start+it+1-1))

         if (msglen.gt.0) then
            reqcnt = reqcnt + 1
#ifdef MPI4
            stupid_msglen = msglen
            stupid_type   = msgtype
            stupid_taskid = Parallel2d_convert_taskid_i(source)
            call MPI_IRECV(tmp2(int_mb(i2_start+it+1-1)),
     >                    stupid_msglen,stupid_double,
     >                    stupid_taskid,
     >                    stupid_type,stupid_world,
     >                    stupid_request,stupid_ierr)
            request(reqcnt) = 0
            request(reqcnt) = stupid_request
#else
            call MPI_IRECV(tmp2(int_mb(i2_start+it+1-1)),
     >                    msglen,MPI_DOUBLE_PRECISION,
     >                    Parallel2d_convert_taskid_i(source),
     >                    msgtype,Parallel_comm_world(),
     >                    request(reqcnt),mpierr)
#endif
         end if

      end do

*     **** send packed array to other processors ****
      do it = 1,np-1

*        **** synchronous send of tmp ****
         proc_to   = mod((taskid)+it, np)
         msglen    = (int_mb(i1_start+it+2-1)
     >              - int_mb(i1_start+it+1-1))
         msgtype   = 54

         if (msglen.gt.0) then
            reqcnt = reqcnt + 1
#ifdef MPI4
            stupid_msglen = msglen
            stupid_type   = msgtype
            stupid_taskid = Parallel2d_convert_taskid_i(proc_to)
            call MPI_ISEND(tmp1(int_mb(i1_start+it+1-1)),
     >                     stupid_msglen,stupid_double,
     >                     stupid_taskid,
     >                     stupid_type,stupid_world,
     >                     stupid_request,stupid_ierr)
            request(reqcnt) = 0
            request(reqcnt) = stupid_request
#else
            call MPI_ISEND(tmp1(int_mb(i1_start+it+1-1)),
     >                     msglen,MPI_DOUBLE_PRECISION,
     >                     Parallel2d_convert_taskid_i(proc_to),
     >                     msgtype,Parallel_comm_world(),
     >                     request(reqcnt),mpierr)
#endif
         end if

      end do
      

*     *** increment qhead0_3t1
      qhead0_3t1 = mod(qhead0_3t1+1,maxq_3t1)

c      if (np.gt.1) call Parallel_mpiWaitAll0(reqcnt,request)

!$OMP END MASTER
!$OMP BARRIER

      call nwpw_timing_end(55)

      return
      end


*     ***********************************
*     *                                 *
*     *    D3dB_r_3to1_Expand_final     *
*     *                                 *
*     ***********************************
      subroutine D3dB_r_3to1_Expand_final(tmp2,A,request,reqcnt)
      implicit none
      real*8 tmp2(*),A(*)
      integer request(*),reqcnt

#include "bafdecls.fh"
#include "D3dB.fh"
#include "errquit.fh"

*     **** r_3to1_block common block ****
      integer maxq_3t1
      parameter (maxq_3t1=6)
      integer qhead0_3t1,qhead1_3t1,qlast_3t1
      integer i1_start_3t1(2),i2_start_3t1(2)
      integer c_iq_to_i1_3t1(2),c_iq_to_i2_3t1(2)
      integer index2_3t1(2)
      common /r_3to1_block/ i1_start_3t1,i2_start_3t1,
     >                      c_iq_to_i1_3t1,c_iq_to_i2_3t1,
     >                      index2_3t1,
     >                      qhead0_3t1,qhead1_3t1,qlast_3t1

*     **** local variables ****
      logical value
      integer taskid,np
      integer proc_to,proc_from,source
      integer itmp,msgtype,msglen,it,i
      integer i1_start,i2_start,c_iq_to_i1,c_iq_to_i2
      integer shift1,shift2,index2


      call nwpw_timing_start(55)
      call Parallel2d_taskid_i(taskid)
      call Parallel2d_np_i(np)

      
      shift1 = qhead1_3t1*(np+1)
      shift2 = qhead1_3t1*n2ft3d(1)
      i1_start = i1_start_3t1(1) + shift1
      i2_start = i2_start_3t1(1) + shift1
      c_iq_to_i1 = c_iq_to_i1_3t1(1) + shift2
      c_iq_to_i2 = c_iq_to_i2_3t1(1) + shift2
      index2 = int_mb(index2_3t1(1)+qhead1_3t1)


!$OMP MASTER
*     **** wait for completion of mp_send, also do a sync ****
      if (np.gt.1) call Parallel_mpiWaitAll0(reqcnt,request)

*     *** increment qhead1_3t1
      qhead1_3t1 = mod(qhead1_3t1+1,maxq_3t1)

!$OMP END MASTER
!$OMP BARRIER

*     **** unpack tmp2(i) array ****
      call Parallel_shared_vector_zero(.true.,n2ft3d(1),A)
!$OMP DO
      do i=1,index2-1
         A(int_mb(c_iq_to_i2+i-1)) = tmp2(i)
      end do
!$OMP END DO

      call nwpw_timing_end(55)

      return
      end


*     ***********************************
*     *                                 *
*     *    D3dB_r_3to1_Expand_start3    *
*     *                                 *
*     ***********************************
      subroutine D3dB_r_3to1_Expand_start3(center,A,tmp1,
     >                                      request,reqcnt)
      implicit none
      integer center(3)
      real*8 A(*)
      real*8 tmp1(*)
      integer request(*),reqcnt

#include "bafdecls.fh"
#include "mpif.h"
#include "D3dB.fh"
#include "errquit.fh"

#ifdef MPI4
#include "stupid_mpi4.fh"
#endif

*     **** r_3to1_block common block ****
      integer maxq_3t1
      parameter (maxq_3t1=6)
      integer qhead0_3t1,qhead1_3t1,qlast_3t1
      integer i1_start_3t1(2),i2_start_3t1(2)
      integer c_iq_to_i1_3t1(2),c_iq_to_i2_3t1(2)
      integer index2_3t1(2)
      common /r_3to1_block/ i1_start_3t1,i2_start_3t1,
     >                      c_iq_to_i1_3t1,c_iq_to_i2_3t1,
     >                      index2_3t1,
     >                      qhead0_3t1,qhead1_3t1,qlast_3t1

*     **** igrid_r_block common block ****
      integer igrid_big(2),igrid_small(2),small(2),pindex1(2),pindex2(2)
      common /igrid_r_block/ igrid_big,igrid_small,small,pindex1,pindex2

*     **** local variables ****
      logical value
      integer taskid,np
      integer proc_to,proc_from,source
      integer indxbig,pbig,indxsmall,psmall,itmp
      integer index1,index2,msgtype,msglen
      integer ismall,jsmall,ksmall,lsmall,ibig,jbig,kbig,lbig
      integer i,it
      integer nxh,nyh,nzh,nxh_small,nyh_small,nzh_small
      integer small_count
      integer i1_start,i2_start,c_iq_to_i1,c_iq_to_i2
      integer shift1,shift2


      call nwpw_timing_start(55)
      call Parallel2d_taskid_i(taskid)
      call Parallel2d_np_i(np)

      shift1 = qlast_3t1*(np+1)
      shift2 = qlast_3t1*n2ft3d(1)
      i1_start = i1_start_3t1(1) + shift1
      i2_start = i2_start_3t1(1) + shift1
      c_iq_to_i1 = c_iq_to_i1_3t1(1) + shift2
      c_iq_to_i2 = c_iq_to_i2_3t1(1) + shift2

      nxh = nx(1)/2
      nyh = ny(1)/2
      nzh = nz(1)/2

      nxh_small = nx(3)/2
      nyh_small = ny(3)/2
      nzh_small = nz(3)/2

!$OMP MASTER
      do it=0,np-1
         int_mb(pindex1(1)+it) = 0
         int_mb(pindex2(1)+it) = 0
      end do
      !*** count number of sends ****
      small_count = 0
      do i=1,n2ft3d(3)
         ismall = int_mb(igrid_small(1)+4*(i-1))
         jsmall = int_mb(igrid_small(1)+4*(i-1)+1)
         ksmall = int_mb(igrid_small(1)+4*(i-1)+2)
         lsmall = int_mb(igrid_small(1)+4*(i-1)+3)
         if (lsmall.gt.0) then

            ibig = ismall+center(1)
            if (ibig.gt.(nxh-1)) ibig = ibig - nx(1)
            if (ibig.lt.(-nxh))  ibig = ibig + nx(1)

            jbig = jsmall+center(2)
            if (jbig.gt.(nyh-1)) jbig = jbig - ny(1)
            if (jbig.lt.(-nyh))  jbig = jbig + ny(1)

            kbig = ksmall+center(3)
            if (kbig.gt.(nzh-1)) kbig = kbig - nz(1)
            if (kbig.lt.(-nzh))  kbig = kbig + nz(1)

            ibig = ibig + nxh
            jbig = jbig + nyh
            kbig = kbig + nzh
            call D3dB_ijktoindex2p(1,ibig+1,jbig+1,kbig+1,
     >                             indxbig,pbig)

            int_mb(small(1)+4*small_count)   = indxbig
            int_mb(small(1)+4*small_count+1) = pbig
            int_mb(small(1)+4*small_count+2) = i
            int_mb(small(1)+4*small_count+3) = int_mb(pindex1(1)+pbig)
            small_count = small_count + 1

            int_mb(pindex1(1)+pbig) = int_mb(pindex1(1)+pbig) + 1
         end if
      end do

      !*** count number of recvs ****
      do i=1,n2ft3d(1)
         ibig = int_mb(igrid_big(1)+4*(i-1))
         jbig = int_mb(igrid_big(1)+4*(i-1)+1)
         kbig = int_mb(igrid_big(1)+4*(i-1)+2)
         lbig = int_mb(igrid_big(1)+4*(i-1)+3)
         if (lbig.gt.0) then

            ismall = ibig - center(1)
            if (ismall.lt.(-nxh))  ismall = ismall + nx(1)
            if (ismall.gt.(nxh-1)) ismall = ismall - nx(1)

            jsmall = jbig - center(2)
            if (jsmall.lt.(-nyh))  jsmall = jsmall + ny(1)
            if (jsmall.gt.(nyh-1)) jsmall = jsmall - ny(1)

            ksmall = kbig - center(3)
            if (ksmall.lt.(-nzh))  ksmall = ksmall + nz(1)
            if (ksmall.gt.(nzh-1)) ksmall = ksmall - nz(1)

            if ((ismall.gt.(-nxh_small-1)).and.
     >          (ismall.lt.nxh_small).and.
     >          (jsmall.gt.(-nzh_small-1)).and.
     >          (jsmall.lt.nyh_small).and.
     >          (ksmall.gt.(-nzh_small-1)).and.
     >          (ksmall.lt.nzh_small)) then
               ismall = ismall + nxh_small
               jsmall = jsmall + nyh_small
               ksmall = ksmall + nzh_small
               call D3dB_ijktoindex2p(3,ismall+1,jsmall+1,ksmall+1,
     >                             indxsmall,psmall)
               int_mb(pindex2(1)+psmall) = int_mb(pindex2(1)+psmall) + 1
            end if
         end if
      end do


      index1 = 1
      index2 = 1
      do it=0,np-1
         proc_to   = mod(taskid+it,np)
         proc_from = mod(taskid-it+np,np)
         int_mb(i1_start+it) = index1
         int_mb(i2_start+it) = index2
         index1 = index1 + int_mb(pindex1(1)+proc_to)
         index2 = index2 + int_mb(pindex2(1)+proc_from)
      end do
      int_mb(i1_start+np) = index1
      int_mb(i2_start+np) = index2

c      do it=0,np-1
c         int_mb(pindex1(1)+it) = 0
c         int_mb(pindex2(1)+it) = 0
c      end do
c
c      do i=1,n2ft3d(3)
c         ismall = int_mb(igrid_small(1)+4*(i-1))
c         jsmall = int_mb(igrid_small(1)+4*(i-1)+1)
c         ksmall = int_mb(igrid_small(1)+4*(i-1)+2)
c         lsmall = int_mb(igrid_small(1)+4*(i-1)+3)
c         if (lsmall.gt.0) then
c
c            ibig = ismall + center(1)
c            if (ibig.gt.(nxh-1)) ibig = ibig - nx(1)
c            if (ibig.lt.(-nxh))  ibig = ibig + nx(1)
c
c            jbig = jsmall + center(2)
c            if (jbig.gt.(nyh-1)) jbig = jbig - ny(1)
c            if (jbig.lt.(-nyh))  jbig = jbig + ny(1)
c
c            kbig = ksmall + center(3)
c            if (kbig.gt.(nzh-1)) kbig = kbig - nz(1)
c            if (kbig.lt.(-nzh))  kbig = kbig + nz(1)
c
c            ibig = ibig + nxh
c            jbig = jbig + nyh
c            kbig = kbig + nzh
cc            ismall = ismall + nxh_small
cc            jsmall = jsmall + nyh_small
cc            ksmall = ksmall + nzh_small
c            call D3dB_ijktoindex2p(1,ibig+1,jbig+1,kbig+1,
c     >                          indxbig,pbig)
cc            call D3dB_ijktoindex2p(3,ismall+1,jsmall+1,ksmall+1,
cc     >                          indxsmall,psmall)
c
c            it = mod(pbig-taskid+np,np)
c            itmp = int_mb(i1_start+it)
c     >           + int_mb(pindex1(1)+pbig)
cc            int_mb(c_iq_to_i1+itmp-1) = indxsmall
c            int_mb(c_iq_to_i1+itmp-1) = i
c            int_mb(c_iq_to_i2+itmp-1) = indxbig
c            int_mb(pindex1(1)+pbig) = int_mb(pindex1(1)+pbig)+1
c
c         end if
c      end do

      do i=1,small_count
         indxbig   = int_mb(small(1)+4*(i-1))
         pbig      = int_mb(small(1)+4*(i-1)+1)
         indxsmall = int_mb(small(1)+4*(i-1)+2)
         itmp      = int_mb(small(1)+4*(i-1)+3)
         it = mod(pbig-taskid+np,np)
         itmp = itmp + int_mb(i1_start+it)
cc         itmp = int_mb(i1_start+it) 
cc     >              + int_mb(pindex1(1)+psmall)
         int_mb(c_iq_to_i1+itmp-1) = indxsmall
         int_mb(c_iq_to_i2+itmp-1) = indxbig
c         int_mb(pindex1(1)+pbig) = int_mb(pindex1(1)+pbig)+1
      end do





*     **** set index2 and increment qlast_1t3 ****
      int_mb(index2_3t1(1)+qlast_3t1) = index2
      qlast_3t1 = mod(qlast_3t1+1,maxq_3t1)

!$OMP END MASTER


!$OMP DO
      do i=1,index1-1
         tmp1(i) = A(int_mb(c_iq_to_i1+i-1))
      end do
!$OMP END DO
      call nwpw_timing_end(55)

      return
      end 



*     ***********************************
*     *                                 *
*     *    D3dB_r_3to1_Expand_mid3      *
*     *                                 *
*     ***********************************
      subroutine D3dB_r_3to1_Expand_mid3(tmp1,tmp2,request,reqcnt)
      implicit none
      real*8 tmp1(*)
      real*8 tmp2(*)
      integer request(*),reqcnt


#include "bafdecls.fh"
#include "mpif.h"
#include "D3dB.fh"
#include "errquit.fh"

#ifdef MPI4
#include "stupid_mpi4.fh"
#endif

*     **** r_3to1_block common block ****
      integer maxq_3t1
      parameter (maxq_3t1=6)
      integer qhead0_3t1,qhead1_3t1,qlast_3t1
      integer i1_start_3t1(2),i2_start_3t1(2)
      integer c_iq_to_i1_3t1(2),c_iq_to_i2_3t1(2)
      integer index2_3t1(2)
      common /r_3to1_block/ i1_start_3t1,i2_start_3t1,
     >                      c_iq_to_i1_3t1,c_iq_to_i2_3t1,
     >                      index2_3t1,
     >                      qhead0_3t1,qhead1_3t1,qlast_3t1

*     **** local variables ****
      logical value
      integer taskid,np
      integer proc_to,proc_from,source
      integer it,itmp,msgtype,msglen
      integer i1_start,i2_start,c_iq_to_i1,c_iq_to_i2
      integer shift1,shift2
      integer mpierr

c     **** external functions ****
      integer  Parallel2d_convert_taskid_i,Parallel_comm_world
      external Parallel2d_convert_taskid_i,Parallel_comm_world

      call nwpw_timing_start(55)
      call Parallel2d_taskid_i(taskid)
      call Parallel2d_np_i(np)

      shift1 = qhead0_3t1*(np+1)
      shift2 = qhead0_3t1*n2ft3d(1)
      i1_start = i1_start_3t1(1) + shift1
      i2_start = i2_start_3t1(1) + shift1
      c_iq_to_i1 = c_iq_to_i1_3t1(1) + shift2
      c_iq_to_i2 = c_iq_to_i2_3t1(1) + shift2


*     **** it = 0, transpose data on same thread ****  
      msglen = int_mb(i2_start+2-1) - int_mb(i2_start+1-1)

      call Parallel_shared_vector_copy(.false.,msglen,
     >               tmp1(int_mb(i1_start+1-1)),
     >               tmp2(int_mb(i2_start+1-1)))

      call Parallel_shared_vector_icopy(.false.,msglen,
     >               int_mb(c_iq_to_i2+int_mb(i1_start+1-1)-1),
     >               int_mb(c_iq_to_i1+int_mb(i2_start+1-1)-1))

!$OMP MASTER

*     **** receive packed array data ****
      reqcnt = 0
      do it = 1,np-1

*        **** synchronous receive of tmp ****
         proc_from = mod((taskid)-it+np, np)

         source=proc_from
         msglen = (int_mb(i2_start+it+2-1)
     >          -  int_mb(i2_start+it+1-1))

         if (msglen.gt.0) then
            msgtype = 53
            reqcnt = reqcnt + 1
#ifdef MPI4
            stupid_msglen = msglen
            stupid_type   = msgtype
            stupid_taskid = Parallel2d_convert_taskid_i(source)
            call MPI_IRECV(int_mb(c_iq_to_i1+int_mb(i2_start+it+1-1)-1),
     >                    2*stupid_msglen,stupid_integer,
     >                    stupid_taskid,
     >                    stupid_type,stupid_world,
     >                    stupid_request,stupid_ierr)
            request(reqcnt) = 0
            request(reqcnt) = stupid_request
#else
            call MPI_IRECV(int_mb(c_iq_to_i1+int_mb(i2_start+it+1-1)-1),
     >                    2*msglen,MPI_INTEGER,
     >                    Parallel2d_convert_taskid_i(source),
     >                    msgtype,Parallel_comm_world(),
     >                    request(reqcnt),mpierr)
#endif
         end if

         if (msglen.gt.0) then
            msgtype = 54
            reqcnt = reqcnt + 1
#ifdef MPI4
            stupid_msglen = msglen
            stupid_type   = msgtype
            stupid_taskid = Parallel2d_convert_taskid_i(source)
            call MPI_IRECV(tmp2(int_mb(i2_start+it+1-1)),
     >                    stupid_msglen,stupid_double,
     >                    stupid_taskid,
     >                    stupid_type,stupid_world,
     >                    stupid_request,stupid_ierr)
            request(reqcnt) = 0
            request(reqcnt) = stupid_request
#else
            call MPI_IRECV(tmp2(int_mb(i2_start+it+1-1)),
     >                    msglen,MPI_DOUBLE_PRECISION,
     >                    Parallel2d_convert_taskid_i(source),
     >                    msgtype,Parallel_comm_world(),
     >                    request(reqcnt),mpierr)
#endif
         end if

      end do


*     **** send packed array to other processors ****
      do it = 1,np-1

*        **** synchronous send of tmp ****
         proc_to   = mod((taskid)+it, np)
         msglen    = (int_mb(i1_start+it+2-1)
     >              - int_mb(i1_start+it+1-1))

         if (msglen.gt.0) then
            msgtype   = 53
            reqcnt = reqcnt + 1
#ifdef MPI4
            stupid_msglen = msglen
            stupid_type   = msgtype
            stupid_taskid = Parallel2d_convert_taskid_i(proc_to)
            call MPI_ISEND(int_mb(c_iq_to_i2+int_mb(i1_start+it+1-1)-1),
     >                     2*stupid_msglen,stupid_integer,
     >                     stupid_taskid,
     >                     stupid_type,stupid_world,
     >                     stupid_request,stupid_ierr)
            request(reqcnt) = 0
            request(reqcnt) = stupid_request
#else
            call MPI_ISEND(int_mb(c_iq_to_i2+int_mb(i1_start+it+1-1)-1),
     >                     2*msglen,MPI_INTEGER,
     >                     Parallel2d_convert_taskid_i(proc_to),
     >                     msgtype,Parallel_comm_world(),
     >                     request(reqcnt),mpierr)
#endif
         end if

         if (msglen.gt.0) then
            msgtype   = 54
            reqcnt = reqcnt + 1
#ifdef MPI4
            stupid_msglen = msglen
            stupid_type   = msgtype
            stupid_taskid = Parallel2d_convert_taskid_i(proc_to)
            call MPI_ISEND(tmp1(int_mb(i1_start+it+1-1)),
     >                     stupid_msglen,stupid_double,
     >                     stupid_taskid,
     >                     stupid_type,stupid_world,
     >                     stupid_request,stupid_ierr)
            request(reqcnt) = 0
            request(reqcnt) = stupid_request
#else
            call MPI_ISEND(tmp1(int_mb(i1_start+it+1-1)),
     >                     msglen,MPI_DOUBLE_PRECISION,
     >                     Parallel2d_convert_taskid_i(proc_to),
     >                     msgtype,Parallel_comm_world(),
     >                     request(reqcnt),mpierr)
#endif
         end if
      end do


*     *** increment qhead0_3t1
      qhead0_3t1 = mod(qhead0_3t1+1,maxq_3t1)

c      if (np.gt.1) call Parallel_mpiWaitAll0(reqcnt,request)

!$OMP END MASTER
!$OMP BARRIER

      call nwpw_timing_end(55)

      return
      end

*     ***********************************
*     *                                 *
*     *    D3dB_r_3to1_Expand_final3    *
*     *                                 *
*     ***********************************
      subroutine D3dB_r_3to1_Expand_final3(tmp2,A,request,reqcnt)
      implicit none
      real*8 tmp2(*),A(*)
      integer request(*),reqcnt

#include "bafdecls.fh"
#include "D3dB.fh"
#include "errquit.fh"

*     **** r_3to1_block common block ****
      integer maxq_3t1
      parameter (maxq_3t1=6)
      integer qhead0_3t1,qhead1_3t1,qlast_3t1
      integer i1_start_3t1(2),i2_start_3t1(2)
      integer c_iq_to_i1_3t1(2),c_iq_to_i2_3t1(2)
      integer index2_3t1(2)
      common /r_3to1_block/ i1_start_3t1,i2_start_3t1,
     >                      c_iq_to_i1_3t1,c_iq_to_i2_3t1,
     >                      index2_3t1,
     >                      qhead0_3t1,qhead1_3t1,qlast_3t1

*     **** local variables ****
      logical value
      integer taskid,np
      integer proc_to,proc_from,source
      integer itmp,msgtype,msglen,it,i
      integer i1_start,i2_start,c_iq_to_i1,c_iq_to_i2
      integer shift1,shift2,index2

      call nwpw_timing_start(55)
      call Parallel2d_taskid_i(taskid)
      call Parallel2d_np_i(np)

      shift1 = qhead1_3t1*(np+1)
      shift2 = qhead1_3t1*n2ft3d(1)
      i1_start = i1_start_3t1(1) + shift1
      i2_start = i2_start_3t1(1) + shift1
      c_iq_to_i1 = c_iq_to_i1_3t1(1) + shift2
      c_iq_to_i2 = c_iq_to_i2_3t1(1) + shift2
      index2 = int_mb(index2_3t1(1)+qhead1_3t1)

!$OMP MASTER
*     **** wait for completion of mp_send, also do a sync ****
      if (np.gt.1) call Parallel_mpiWaitAll0(reqcnt,request)

*     *** increment qhead1_3t1
      qhead1_3t1 = mod(qhead1_3t1+1,maxq_3t1)

!$OMP END MASTER
!$OMP BARRIER

*     **** unpack tmp2(i) array ****
      call Parallel_shared_vector_zero(.true.,n2ft3d(1),A)
!$OMP DO
      do i=1,index2-1
         A(int_mb(c_iq_to_i1+i-1)) = tmp2(i)
      end do
!$OMP END DO

      call nwpw_timing_end(55)

      return
      end












*     ***********************************
*     *                                 *
*     *   D3dB_r_1to3_Contract_start2   *
*     *                                 *
*     ***********************************
      subroutine D3dB_r_1to3_Contract_start2(center,A,tmp1,
     >                                      request,reqcnt)
      implicit none
      integer center(3)
      real*8 A(*)
      real*8 tmp1(*)
      integer request(*),reqcnt

#include "bafdecls.fh"
#include "mpif.h"
#include "D3dB.fh"
#include "errquit.fh"

#ifdef MPI4
#include "stupid_mpi4.fh"
#endif

*     **** r_1to3_block common block ****
      integer maxq_1t3
      parameter (maxq_1t3=6)
      integer qhead0_1t3,qhead1_1t3,qlast_1t3
      integer i1_start_1t3(2),i2_start_1t3(2)
      integer c_iq_to_i1_1t3(2),c_iq_to_i2_1t3(2)
      integer index2_1t3(2)
      common /r_1to3_block/ i1_start_1t3,i2_start_1t3,
     >                      c_iq_to_i1_1t3,c_iq_to_i2_1t3,
     >                      index2_1t3,
     >                      qhead0_1t3,qhead1_1t3,qlast_1t3

*     **** local variables ****
      logical value
      integer taskid,np
      integer proc_to,proc_from,source
      integer indxbig,pbig,indxsmall,psmall,itmp
      integer index1,index2,msgtype,msglen
      integer k1,k2,k3,l1,l2,l3,i,j,k,il,jl,kl,it
      integer nxh,nyh,nzh,nxh_small,nyh_small,nzh_small
      integer pindex1(2),pindex2(2)
      integer i1_start,i2_start,c_iq_to_i1,c_iq_to_i2
      integer shift1,shift2


      call Parallel2d_taskid_i(taskid)
      call Parallel2d_np_i(np)
      call Parallel_shared_vector_zero(.false.,n2ft3d(1),tmp1)

      shift1 = qlast_1t3*(np+1)
      shift2 = qlast_1t3*n2ft3d(1)
      i1_start = i1_start_1t3(1) + shift1
      i2_start = i2_start_1t3(1) + shift1
      c_iq_to_i1 = c_iq_to_i1_1t3(1) + shift2
      c_iq_to_i2 = c_iq_to_i2_1t3(1) + shift2

      nxh = nx(1)/2
      nyh = ny(1)/2
      nzh = nz(1)/2

      nxh_small = nx(3)/2
      nyh_small = ny(3)/2
      nzh_small = nz(3)/2

!$OMP MASTER

      index1 = 1
      index2 = 1
      do it=0,np-1
         proc_to   = mod(taskid+it,np)
         proc_from = mod(taskid-it+np,np)

         int_mb(i1_start+it) = index1
         int_mb(i2_start+it) = index2

         do k3=-nzh_small,nzh_small-1
            l3 = k3 + center(3)
            if (l3.lt.(-nzh))  l3 = l3 + nz(1)
            if (l3.gt.(nzh-1)) l3 = l3 - nz(1)
            k  = k3 + nzh_small
            kl = l3 + nzh

            do k2=-nyh_small,nyh_small-1
               l2 = k2 + center(2)
               if (l2.lt.(-nyh))  l2 = l2 + ny(1)
               if (l2.gt.(nyh-1)) l2 = l2 - ny(1)
               j  = k2 + nyh_small
               jl = l2 + nyh

               do k1=-nxh_small,nxh_small-1
                  l1 = k1 + center(1)
                  if (l1.lt.(-nxh))  l1 = l1 + nx(1)
                  if (l1.gt.(nxh-1)) l1 = l1 - nx(1)
                  i  = k1 + nxh_small
                  il = l1 + nxh

                  call D3dB_ijktoindex2p(1,il+1,jl+1,kl+1,indxbig,pbig)
                  call D3dB_ijktoindex2p(3,i+1,j+1,k+1,indxsmall,psmall)

*                 **** packing scheme **** 
                  if ((pbig.eq.taskid).and.(psmall.eq.proc_to)) then
                     int_mb(c_iq_to_i1+index1-1) = indxbig
                     index1 = index1 + 1
                  end if

*                 **** unpacking scheme **** 
                  if ((psmall.eq.taskid).and.(pbig.eq.proc_from)) then
                     int_mb(c_iq_to_i2+index2-1) = indxsmall
                     index2 = index2 + 1
                  end if
               end do
            end do
         end do
      end do
      int_mb(i1_start+np) = index1
      int_mb(i2_start+np) = index2

*     **** set index2 and increment qlast_1t3 ****
      int_mb(index2_1t3(1)+qlast_1t3) = index2
      qlast_1t3 = mod(qlast_1t3+1,maxq_1t3)

!$OMP END MASTER
!$OMP BARRIER


*     **** pack A(i) array ****
!$OMP DO
      do i=1,index1-1
         tmp1(i) = A(int_mb(c_iq_to_i1+i-1))
      end do
!$OMP END DO



      return
      end





*     ***********************************
*     *                                 *
*     *    D3dB_r_3to1_Expand_start2    *
*     *                                 *
*     ***********************************
      subroutine D3dB_r_3to1_Expand_start2(center,A,tmp1,
     >                                     request,reqcnt)
      implicit none
      integer center(3)
      real*8 A(*)
      real*8 tmp1(*)
      integer request(*),reqcnt

#include "bafdecls.fh"
#include "mpif.h"
#include "D3dB.fh"
#include "errquit.fh"

#ifdef MPI4
#include "stupid_mpi4.fh"
#endif

*     **** r_3to1_block common block ****
      integer maxq_3t1
      parameter (maxq_3t1=6)
      integer qhead0_3t1,qhead1_3t1,qlast_3t1
      integer i1_start_3t1(2),i2_start_3t1(2)
      integer c_iq_to_i1_3t1(2),c_iq_to_i2_3t1(2)
      integer index2_3t1(2)
      common /r_3to1_block/ i1_start_3t1,i2_start_3t1,
     >                      c_iq_to_i1_3t1,c_iq_to_i2_3t1,
     >                      index2_3t1,
     >                      qhead0_3t1,qhead1_3t1,qlast_3t1

*     **** local variables ****
      logical value
      integer taskid,np
      integer proc_to,proc_from,source
      integer indxbig,pbig,indxsmall,psmall,itmp
      integer index1,index2,msgtype,msglen
      integer k1,k2,k3,l1,l2,l3,i,j,k,il,jl,kl,it
      integer nxh,nyh,nzh,nxh_small,nyh_small,nzh_small
      integer pindex1(2),pindex2(2)
      integer i1_start,i2_start,c_iq_to_i1,c_iq_to_i2
      integer shift1,shift2


      call Parallel2d_taskid_i(taskid)
      call Parallel2d_np_i(np)
      call Parallel_shared_vector_zero(.false.,n2ft3d(1),tmp1)

      shift1 = qlast_3t1*(np+1)
      shift2 = qlast_3t1*n2ft3d(1)
      i1_start = i1_start_3t1(1) + shift1
      i2_start = i2_start_3t1(1) + shift1
      c_iq_to_i1 = c_iq_to_i1_3t1(1) + shift2
      c_iq_to_i2 = c_iq_to_i2_3t1(1) + shift2

      nxh = nx(1)/2
      nyh = ny(1)/2
      nzh = nz(1)/2

      nxh_small = nx(3)/2
      nyh_small = ny(3)/2
      nzh_small = nz(3)/2

!$OMP MASTER

      index1 = 1
      index2 = 1
      do it=0,np-1
         proc_to   = mod(taskid+it,np)
         proc_from = mod(taskid-it+np,np)

         int_mb(i1_start+it) = index1
         int_mb(i2_start+it) = index2

         do k3=-nzh_small,nzh_small-1
            l3 = k3 + center(3)
            if (l3.lt.(-nzh))  l3 = l3 + nz(1)
            if (l3.gt.(nzh-1)) l3 = l3 - nz(1)
            k  = k3 + nzh_small
            kl = l3 + nzh

            do k2=-nyh_small,nyh_small-1
               l2 = k2 + center(2)
               if (l2.lt.(-nyh))  l2 = l2 + ny(1)
               if (l2.gt.(nyh-1)) l2 = l2 - ny(1)
               j  = k2 + nyh_small
               jl = l2 + nyh

               do k1=-nxh_small,nxh_small-1
                  l1 = k1 + center(1)
                  if (l1.lt.(-nxh))  l1 = l1 + nx(1)
                  if (l1.gt.(nxh-1)) l1 = l1 - nx(1)
                  i  = k1 + nxh_small
                  il = l1 + nxh

                  call D3dB_ijktoindex2p(3,i+1,j+1,k+1,indxsmall,psmall)
                  call D3dB_ijktoindex2p(1,il+1,jl+1,kl+1,indxbig,pbig)

*                 **** packing scheme **** 
                  if ((psmall.eq.taskid).and.(pbig.eq.proc_to)) then
                     int_mb(c_iq_to_i1+index1-1) = indxsmall
                     index1 = index1 + 1
                  end if

*                 **** unpacking scheme **** 
                  if ((pbig.eq.taskid).and.(psmall.eq.proc_from)) then
                     int_mb(c_iq_to_i2+index2-1) = indxbig
                     index2 = index2 + 1
                  end if
               end do
            end do
         end do
      end do
      int_mb(i1_start+np) = index1
      int_mb(i2_start+np) = index2

*     **** set index2 and increment qlast_3t1 ****
      int_mb(index2_3t1(1)+qlast_3t1) = index2
      qlast_3t1 = mod(qlast_3t1+1,maxq_3t1)

!$OMP END MASTER
!$OMP BARRIER


*     **** pack A(i) array ****
!$OMP DO
      do i=1,index1-1
         tmp1(i) = A(int_mb(c_iq_to_i1+i-1))
      end do
!$OMP END DO


      return
      end








