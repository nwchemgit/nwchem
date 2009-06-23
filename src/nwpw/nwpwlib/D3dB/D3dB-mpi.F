#define NBLOCKS 2

*
* $Id: D3dB-mpi.F,v 1.21 2007-04-12 17:10:13 bylaska Exp $
*

#ifdef MPI
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

#include "mafdecls.fh"
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
      integer msglen,type
      integer request(2),reqcnt
      integer taskid,np
      integer mpierr

*     **** external functions ****
      integer  Parallel2d_convert_taskid_i
      external Parallel2d_convert_taskid_i
      
      call Parallel2d_taskid_i(taskid)
      call Parallel2d_np_i(np)

*     **** allocate memory ****
      value = MA_push_get(mt_int,np,
     >                     'request',request(2),request(1))
      if (.not. value)
     > call errquit(' D3dB_c_transpose_jk:out of stack',0,MA_ERR)


*     **** pack A(i) array ****
#ifndef CRAY
!DIR$ ivdep
#endif
      do i=1,nfft3d(nb) 
         tmp1(int_mb(iq_to_i1(1,nb)+i-1)) = A(i)
      end do

*     **** it = 0, transpose data on same thread ****  
      msglen = int_mb(i2_start(1,nb)+2-1) - int_mb(i2_start(1,nb)+1-1)
#ifndef CRAY
!DIR$ ivdep
#endif
      do i=1,msglen
         tmp2(int_mb(i2_start(1,nb)+1-1)+i-1) 
     > = tmp1(int_mb(i1_start(1,nb)+1-1)+i-1)
      end do
         

*     **** receive packed array data ****
      reqcnt = 0
      do it = 1,np-1

*        **** synchronous receive of tmp ****
         proc_from = mod((taskid)-it+np, np) 

         msgtype = 1
         source=proc_from
         type   = 1
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
     >                    msgtype,MPI_COMM_WORLD,
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
         msgtype   = 1
      
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
     >                     msgtype,MPI_COMM_WORLD,mpierr)
#endif
         end if

      end do

   
*     **** wait for completion of mp_send, also do a sync ****
      if (np.gt.1) call Parallel_mpiWaitAll(reqcnt,int_mb(request(1)))


*     **** unpack A(i) array ****
#ifndef CRAY
!DIR$ ivdep
#endif
      do i=1,nfft3d(nb) !***(nx(nb)/2+1)*ny(nb)*nq(nb)
         A(i) = tmp2(int_mb(iq_to_i2(1,nb)+i-1))
      end do
      

      !*** deallocate memory ***
      value =           MA_pop_stack(request(2))
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

#include "mafdecls.fh"
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
      integer msglen,msgtype,type
      integer proc_from,proc_to
      integer request(2),reqcnt
      integer taskid,np
      integer index1,index2
      integer mpierr

*     **** external functions ****
      integer  Parallel2d_convert_taskid_i
      external Parallel2d_convert_taskid_i

      call Parallel2d_taskid_i(taskid)
      call Parallel2d_np_i(np)

*     **** allocate memory ****
      value = MA_push_get(mt_int,np,
     >                     'request',request(2),request(1))
      if (.not. value)
     > call errquit(' D3dB_c_timereverse:out of stack',0,MA_ERR)


*     **** pack A(i) array ****
#ifndef CRAY
!DIR$ ivdep
#endif
      do index1=int_mb(i1_start(1,nb)+1-1),
     >         (int_mb(i1_start(1,nb)+np+1-1)-1)
         tmp1(index1) = A(int_mb(iq_to_i1(1,nb)+index1-1))
      end do

*     **** it = 0, transpose data on same thread ****
      msglen = int_mb(i2_start(1,nb)+2-1) - int_mb(i2_start(1,nb)+1-1)
      call dcopy(2*msglen,tmp1(int_mb(i1_start(1,nb))),1,
     >                    tmp2(int_mb(i2_start(1,nb))),1)


*     **** receive packed array data ****
      reqcnt = 0
      do it = 1,np-1

*        **** synchronous receive of tmp ****
         proc_from = mod((taskid)-it+np, np)

         msgtype = 1
         source=proc_from
         type   = 1
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
     >                    msgtype,MPI_COMM_WORLD,
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
         msgtype   = 1

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
     >                     msgtype,MPI_COMM_WORLD,mpierr)
#endif
         end if

      end do


*     **** wait for completion of mp_send, also do a sync ****
      if (np.gt.1) call Parallel_mpiWaitAll(reqcnt,int_mb(request(1)))


*     **** unpack A(i) array ****
#ifndef CRAY
!DIR$ ivdep
#endif
      do index2=int_mb(i2_start(1,nb)+1-1),
     >         (int_mb(i2_start(1,nb)+np+1-1)-1)
        A(int_mb(iq_to_i2(1,nb)+index2-1))=dconjg(tmp2(index2))
      end do


      !*** deallocate memory ***
      value =           MA_pop_stack(request(2))
      if (.not. value)
     > call errquit(' D3dB_c_timereverse:popping stack',0,MA_ERR)


      return
      end



#endif





#ifdef MPI
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

#include "mafdecls.fh"
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
      integer msglen,type
      integer request(2),reqcnt
      integer taskid,np
      integer mpierr

*     **** external functions ****
      integer  Parallel2d_convert_taskid_i
      external Parallel2d_convert_taskid_i

      
      call Parallel2d_taskid_i(taskid)
      call Parallel2d_np_i(np)

*     **** allocate memory ****
      value = MA_push_get(mt_int,np,
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
      do i=1,nnfft3d
         tmp1(int_mb(h_iq_to_i1(1,op,nb)+i-1)) = A(i)
      end do

*     **** it = 0, transpose data on same thread ****  
      msglen = int_mb(h_i2_start(1,op,nb)+2-1) 
     >       - int_mb(h_i2_start(1,op,nb)+1-1)
      call dcopy(2*msglen,tmp1(int_mb(h_i1_start(1,op,nb))),1,
     >                    tmp2(int_mb(h_i2_start(1,op,nb))),1)
         

*     **** receive packed array data ****
      reqcnt = 0
      do it = 1,np-1

*        **** synchronous receive of tmp ****
         proc_from = mod((taskid)-it+np, np) 

         msgtype = 1
         source=proc_from
         type   = 1
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
     >                    msgtype,MPI_COMM_WORLD,
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
         msgtype   = 1
      
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
     >                     msgtype,MPI_COMM_WORLD,mpierr)
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
      value =           MA_pop_stack(request(2))
      if (.not. value) 
     > call errquit(' D3dB_c_transpose_ijk:popping stack',0,MA_ERR)
      
      return
      end
#endif




#ifdef MPI
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

#include "mafdecls.fh"
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
      integer msglen,type
c     integer status(MPI_STATUS_SIZE,512),request(512),reqcnt
      integer request(2),reqcnt
      integer taskid,np
      integer mpierr

*     **** external functions ****
      integer  Parallel2d_convert_taskid_i
      external Parallel2d_convert_taskid_i
      
      call Parallel2d_taskid_i(taskid)
      call Parallel2d_np_i(np)


*     **** allocate memory ****
      value = MA_push_get(mt_int,np,
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

         msgtype = 1
         source=proc_from
         type   = 1
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
     >                    msgtype,MPI_COMM_WORLD,
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
         msgtype   = 1
      
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
     >                     msgtype,MPI_COMM_WORLD,mpierr)
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
      value =           MA_pop_stack(request(2))
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

*     **** external functions ****
      integer  Parallel2d_comm_i
      external Parallel2d_comm_i

      call Parallel2d_np_i(np)
      if (np.gt.1) then

#ifdef MPI4
         stupid_msglen = 1
         stupid_tsum = sum
         call MPI_Allreduce(stupid_tsum,stupid_tsum2,stupid_msglen,
     >                      stupid_integer,
     >                      stupid_sum,stupid_comm_i,stupid_ierr)
         sum = stupid_tsum2
#else
         msglen = 1
         call MPI_Allreduce(sum,sumall,msglen,MPI_INTEGER,
     >                       MPI_SUM,Parallel2d_comm_i(),mpierr)
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

#include "mafdecls.fh"
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
      value = MA_push_get(mt_dbl,n,'sumall',sumall(2),sumall(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)

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
      value = MA_pop_stack(sumall(2))

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

#include "mafdecls.fh"
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
     >                stupid_sum,stupid_comm_i,stupid_ierr)
      call stupid_icopy4to8(n,int_mb(sumall(1)),sum)
      if (.not.MA_pop_stack(sumall1(2)))
     >  call errquit('error popping stack',0, MA_ERR)
#else
      msglen = n
      call MPI_Allreduce(sum,int_mb(sumall(1)),msglen,
     >                MPI_INTEGER,
     >                MPI_SUM,Parallel2d_comm_i(),mpierr)
      call icopy(n,int_mb(sumall(1)),1,sum,1)
#endif

      value = MA_pop_stack(sumall(2))
      if (.not. value) call errquit('error popping stack',0, MA_ERR)

      end if

      call nwpw_timing_end(2)
      return
      end

#endif


#ifdef MPI

      
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

#include "mafdecls.fh"
#include "errquit.fh"
#include "mpif.h"
#include "D3dB.fh"

#ifdef MPI4
#include "stupid_mpi4.fh"
#endif

*     **** indexing variables ****
      integer iq_to_i1(2,0:1)
      integer iq_to_i2(2,0:1)
      integer iz_to_i2(2,0:1)
      integer i1_start(2,0:1)
      integer i2_start(2,0:1)
      common / ptrans_blk1 / iq_to_i1,iq_to_i2,iz_to_i2,
     >                       i1_start,i2_start


*     **** local variables ***
      logical value
      integer proc_to,proc_from,it
      integer msgtype,source
c     integer msglen,type,nbytes
c     integer msgid
      integer msglen,type
c     integer status(MPI_STATUS_SIZE,512),request(512),reqcnt
      integer request(2),reqcnt
      integer taskid,np
      integer mpierr,n1,n2

*     **** external functions ****
      integer  Parallel2d_convert_taskid_i
      external Parallel2d_convert_taskid_i
      
      call Parallel2d_taskid_i(taskid)
      call Parallel2d_np_i(np)

*     **** allocate memory ****
      value = MA_push_get(mt_int,np,
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

         msgtype = 1
         source=proc_from
         type   = 1
         msglen = (int_mb(i2_start(1,nbb)+it+2-1) 
     >          -  int_mb(i2_start(1,nbb)+it+1-1))
         

         if (msglen.gt.0) then
            reqcnt = reqcnt + 1
#ifdef MPI4
            stupid_msglen = msglen
            stupid_type   = type
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
     >                    type,MPI_COMM_WORLD,
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
         msgtype   = 1
      
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
     >                     msgtype,MPI_COMM_WORLD,mpierr)
#endif
         end if

      end do
   
*     **** wait for completion of mp_send, also do a sync ****
      if (np.gt.1) call Parallel_mpiWaitAll(reqcnt,int_mb(request(1)))


*     **** unpack A(i) array ****
      call D3dB_pfft_index2_copy(n2,int_mb(iq_to_i2(1,nbb)),tmp2,A)
      call D3dB_pfft_index2_zero(nfft3d(1)-n2,int_mb(iz_to_i2(1,nbb)),A)
      

      !*** deallocate memory ***
      value =           MA_pop_stack(request(2))
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

#include "mafdecls.fh"
#include "errquit.fh"
#include "mpif.h"
#include "D3dB.fh"

#ifdef MPI4
#include "stupid_mpi4.fh"
#endif


*     **** indexing variables ****
      integer iq_to_i1(2,0:1)
      integer iq_to_i2(2,0:1)
      integer iz_to_i2(2,0:1)
      integer i1_start(2,0:1)
      integer i2_start(2,0:1)
      common / ptrans_blk2 / iq_to_i1,iq_to_i2,iz_to_i2,
     >                       i1_start,i2_start


*     **** local variables ***
      logical value
      integer proc_to,proc_from,it
      integer msgtype,source
      integer msglen,type
      integer request(2),reqcnt
      integer taskid,np
      integer mpierr,n1,n2

c     **** external functions ****
      integer  Parallel2d_convert_taskid_i
      external Parallel2d_convert_taskid_i
      
      call Parallel2d_taskid_i(taskid)
      call Parallel2d_np_i(np)

*     **** allocate memory ****
      value = MA_push_get(mt_int,np,
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

         msgtype = 1
         source=proc_from
         type   = 1
         msglen = (int_mb(i2_start(1,nbb)+it+2-1) 
     >          -  int_mb(i2_start(1,nbb)+it+1-1))
         

         if (msglen.gt.0) then
            reqcnt = reqcnt + 1
#ifdef MPI4
            stupid_msglen = msglen
            stupid_type   = type
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
     >                    type,MPI_COMM_WORLD,
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
         msgtype   = 1
      
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
     >                     msgtype,MPI_COMM_WORLD,mpierr)
#endif
         end if

      end do
   
*     **** wait for completion of mp_send, also do a sync ****
      if (np.gt.1) call Parallel_mpiWaitAll(reqcnt,int_mb(request(1)))


*     **** unpack A(i) array ****
      call D3dB_pfft_index2_copy(n2,int_mb(iq_to_i2(1,nbb)),tmp2,A)
      call D3dB_pfft_index2_zero(nfft3d(1)-n2,int_mb(iz_to_i2(1,nbb)),A)
      

      !*** deallocate memory ***
      value =           MA_pop_stack(request(2))
      if (.not. value)
     > call errquit(' D3dB_c_ptranspose2_jk:popping stack',0,MA_ERR)

      return
      end


#endif



#ifdef MPI
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

#include "mafdecls.fh"
#include "errquit.fh"
#include "mpif.h"
#include "D3dB.fh"

#ifdef MPI4
#include "stupid_mpi4.fh"
#endif


*     **** indexing variables ****
      integer h_iq_to_i1(2,6,0:1)
      integer h_iq_to_i2(2,6,0:1)
      integer h_iz_to_i2(2,6,0:1)
      integer h_iz_to_i2_count(6,0:1)
      integer h_i1_start(2,6,0:1)
      integer h_i2_start(2,6,0:1)
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
      integer msglen,type
      integer request(2),reqcnt
      integer taskid,np
      integer mpierr

c     **** external functions ****
      integer  Parallel2d_convert_taskid_i
      external Parallel2d_convert_taskid_i
      
      call Parallel2d_taskid_i(taskid)
      call Parallel2d_np_i(np)

*     **** allocate memory ****
      value = MA_push_get(mt_int,np,
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
      call dcopy(2*msglen,
     >           tmp1(int_mb(h_i1_start(1,op,nbb))),1,
     >           tmp2(int_mb(h_i2_start(1,op,nbb))),1)


*     **** receive packed array data ****
      reqcnt = 0
      do it = 1,np-1

*        **** synchronous receive of tmp ****
         proc_from = mod((taskid)-it+np, np) 

         msgtype = 1
         source=proc_from
         type   = 1
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
     >                    msgtype,MPI_COMM_WORLD,
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
         msgtype   = 1
      
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
     >                     msgtype,MPI_COMM_WORLD,mpierr)
#endif
         end if

      end do

   
*     **** wait for completion of mp_send, also do a sync ****
      if (np.gt.1) call Parallel_mpiWaitAll(reqcnt,int_mb(request(1)))


*     **** unpack A(i) array ****
      call D3dB_pfft_index2_copy(n2,int_mb(h_iq_to_i2(1,op,nbb)),tmp2,A)
      call D3dB_pfft_index2_zero(n3,int_mb(h_iz_to_i2(1,op,nbb)),A)


      !*** deallocate memory ***
      value =           MA_pop_stack(request(2))
      if (.not. value) 
     > call errquit(' D3dB_c_transpose_ijk:popping stack',0,MA_ERR)
      
      return
      end
#endif




#ifdef MPI
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

#include "mafdecls.fh"
#include "mpif.h"
#include "balance_common.fh"

#ifdef MPI4
#include "stupid_mpi4.fh"
#endif


*     **** local variables ****
      integer j,np
      integer pto,pfrom,msglen,indx,mpierr

c     **** external functions ****
      integer  Parallel2d_convert_taskid_i
      external Parallel2d_convert_taskid_i


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
     >                     msgtype,MPI_COMM_WORLD,
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
     >                    msgtype,MPI_COMM_WORLD,
     >                    request(reqcnt),mpierr)
#endif
            end if

         end do
      end if

      return
      end
#endif




#ifdef MPI
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
#endif

#ifdef MPI
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

#include "mafdecls.fh"
#include "mpif.h"
#include "balance_common.fh"

#ifdef MPI4
#include "stupid_mpi4.fh"
#endif

*     **** local variables ****
      integer j,np
      integer pto,pfrom,msglen,indx,mpierr

c     **** external functions ****
      integer  Parallel2d_convert_taskid_i
      external Parallel2d_convert_taskid_i


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
     >                    msgtype,MPI_COMM_WORLD,
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
     >                     msgtype,MPI_COMM_WORLD,
     >                     request(reqcnt),mpierr)
#endif

            end if

         end do
      end if

      return
      end
#endif


#ifdef MPI
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

#endif



#ifdef MPI
      
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

#include "mafdecls.fh"
#include "errquit.fh"
#include "mpif.h"
#include "D3dB.fh"

#ifdef MPI4
#include "stupid_mpi4.fh"
#endif


*     **** indexing variables ****
      integer iq_to_i1(2,0:1)
      integer iq_to_i2(2,0:1)
      integer iz_to_i2(2,0:1)
      integer i1_start(2,0:1)
      integer i2_start(2,0:1)
      common / ptrans_blk1 / iq_to_i1,iq_to_i2,iz_to_i2,
     >                       i1_start,i2_start


*     **** local variables ***
      integer proc_to,proc_from,it
      integer source
      integer msglen
      integer taskid,np
      integer mpierr,n1

c     **** external functions ****
      integer  Parallel2d_convert_taskid_i
      external Parallel2d_convert_taskid_i
      
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
     >                    msgtype,MPI_COMM_WORLD,
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
     >                     msgtype,MPI_COMM_WORLD,
     >                     request(reqcnt),mpierr)
#endif
         end if

      end do
   

      return
      end

#endif

#ifdef MPI
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

#include "mafdecls.fh"
#include "D3dB.fh"


*     **** indexing variables ****
      integer iq_to_i1(2,0:1)
      integer iq_to_i2(2,0:1)
      integer iz_to_i2(2,0:1)
      integer i1_start(2,0:1)
      integer i2_start(2,0:1)
      common / ptrans_blk1 / iq_to_i1,iq_to_i2,iz_to_i2,
     >                       i1_start,i2_start

*     **** local variables ***
      integer np,n2
      
      call Parallel2d_np_i(np)
      if (np.gt.1) call Parallel_mpiWaitAll(reqcnt,request)

*     **** unpack A(i) array ****
      n2 = int_mb(i2_start(1,nbb)+np) - 1
      call D3dB_pfft_index2_copy(n2,int_mb(iq_to_i2(1,nbb)),tmp2,A)
      call D3dB_pfft_index2_zero(nfft3d(1)-n2,int_mb(iz_to_i2(1,nbb)),A)
      return
      end

#endif



#ifdef MPI

      
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

#include "mafdecls.fh"
#include "errquit.fh"
#include "mpif.h"
#include "D3dB.fh"

#ifdef MPI4
#include "stupid_mpi4.fh"
#endif

*     **** indexing variables ****
      integer iq_to_i1(2,0:1)
      integer iq_to_i2(2,0:1)
      integer iz_to_i2(2,0:1)
      integer i1_start(2,0:1)
      integer i2_start(2,0:1)
      common / ptrans_blk2 / iq_to_i1,iq_to_i2,iz_to_i2,
     >                       i1_start,i2_start


*     **** local variables ***
      integer proc_to,proc_from,it
      integer source
      integer msglen
      integer taskid,np
      integer mpierr,n1

c     **** external functions ****
      integer  Parallel2d_convert_taskid_i
      external Parallel2d_convert_taskid_i
      
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
     >                    msgtype,MPI_COMM_WORLD,
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
     >                     msgtype,MPI_COMM_WORLD,
     >                     request(reqcnt),mpierr)
#endif
         end if

      end do
   

      return
      end

#endif

#ifdef MPI
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

#include "mafdecls.fh"
#include "D3dB.fh"

*     **** indexing variables ****
      integer iq_to_i1(2,0:1)
      integer iq_to_i2(2,0:1)
      integer iz_to_i2(2,0:1)
      integer i1_start(2,0:1)
      integer i2_start(2,0:1)
      common / ptrans_blk2 / iq_to_i1,iq_to_i2,iz_to_i2,
     >                       i1_start,i2_start


*     **** local variables ***
      logical value
      integer status(2),request2(2)
      integer np,n2
      
      call Parallel2d_np_i(np)
      if (np.gt.1) call Parallel_mpiWaitAll(reqcnt,request)

*     **** unpack A(i) array ****
      n2 = int_mb(i2_start(1,nbb)+np) - 1
      call D3dB_pfft_index2_copy(n2,int_mb(iq_to_i2(1,nbb)),tmp2,A)
      call D3dB_pfft_index2_zero(nfft3d(1)-n2,int_mb(iz_to_i2(1,nbb)),A)
      return
      end

#endif


#ifdef MPI
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

#include "mafdecls.fh"
#include "errquit.fh"
#include "mpif.h"
#include "D3dB.fh"

#ifdef MPI4
#include "stupid_mpi4.fh"
#endif


*     **** indexing variables ****
      integer h_iq_to_i1(2,6,0:1)
      integer h_iq_to_i2(2,6,0:1)
      integer h_iz_to_i2(2,6,0:1)
      integer h_iz_to_i2_count(6,0:1)
      integer h_i1_start(2,6,0:1)
      integer h_i2_start(2,6,0:1)
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
      integer  Parallel2d_convert_taskid_i
      external Parallel2d_convert_taskid_i
      
      call Parallel2d_taskid_i(taskid)
      call Parallel2d_np_i(np)


*     **** pack A(i) array ****
      n1 = int_mb(h_i1_start(1,op,nbb)+np) - 1
      call D3dB_pfft_index1_copy(n1,int_mb(h_iq_to_i1(1,op,nbb)),A,tmp1)


*     **** it = 0, transpose data on same thread ****  
      msglen = int_mb(h_i2_start(1,op,nbb)+2-1) 
     >       - int_mb(h_i2_start(1,op,nbb)+1-1)
      call dcopy(2*msglen,
     >           tmp1(int_mb(h_i1_start(1,op,nbb))),1,
     >           tmp2(int_mb(h_i2_start(1,op,nbb))),1)


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
     >                    msgtype,MPI_COMM_WORLD,
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
     >                     msgtype,MPI_COMM_WORLD,
     >                     request(reqcnt),mpierr)
#endif
         end if

      end do

      return
      end

#endif


#ifdef MPI
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

#include "mafdecls.fh"
#include "errquit.fh"
#include "D3dB.fh"

*     **** indexing variables ****
      integer h_iq_to_i1(2,6,0:1)
      integer h_iq_to_i2(2,6,0:1)
      integer h_iz_to_i2(2,6,0:1)
      integer h_iz_to_i2_count(6,0:1)
      integer h_i1_start(2,6,0:1)
      integer h_i2_start(2,6,0:1)
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

#endif




#ifdef MPI
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

#include "mafdecls.fh"
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
      integer i
      integer it
      integer source
      integer msglen
      integer proc_from,proc_to
      integer taskid,np
      integer index1
      integer mpierr

c     **** external functions ****
      integer  Parallel2d_convert_taskid_i
      external Parallel2d_convert_taskid_i

      call Parallel2d_taskid_i(taskid)
      call Parallel2d_np_i(np)


*     **** pack A(i) array ****
#ifndef CRAY
!DIR$ ivdep
#endif
      do index1=int_mb(i1_start(1,nb)+1-1),
     >         (int_mb(i1_start(1,nb)+np+1-1)-1)
         tmp1(index1) = A(int_mb(iq_to_i1(1,nb)+index1-1))
      end do

*     **** it = 0, transpose data on same thread ****
      msglen = int_mb(i2_start(1,nb)+2-1) - int_mb(i2_start(1,nb)+1-1)
#ifndef CRAY
!DIR$ ivdep
#endif
      do i=1,msglen
         tmp2(int_mb(i2_start(1,nb)+1-1)+i-1)
     > = tmp1(int_mb(i1_start(1,nb)+1-1)+i-1)
      end do


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
     >                    msgtype,MPI_COMM_WORLD,
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
     >                     msgtype,MPI_COMM_WORLD,
     >                     request(reqcnt),mpierr)
#endif
         end if

      end do

      return
      end

#endif




#ifdef MPI
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

#include "mafdecls.fh"
#include "D3dB.fh"

*     **** indexing variables ****
      integer iq_to_i1(2,NBLOCKS)
      integer iq_to_i2(2,NBLOCKS)
      integer i1_start(2,NBLOCKS)
      integer i2_start(2,NBLOCKS)
      common / timereverse_blk / iq_to_i1,iq_to_i2,i1_start,i2_start


*     **** local variables ***
      logical value
      integer np
      integer index2

      call Parallel2d_np_i(np)

      if (np.gt.1) call Parallel_mpiWaitAll(reqcnt,request)

*     **** unpack A(i) array ****
#ifndef CRAY
!DIR$ ivdep
#endif
      do index2=int_mb(i2_start(1,nb)+1-1),
     >         (int_mb(i2_start(1,nb)+np+1-1)-1)
        A(int_mb(iq_to_i2(1,nb)+index2-1))=dconjg(tmp2(index2))
      end do

      return
      end

#endif

