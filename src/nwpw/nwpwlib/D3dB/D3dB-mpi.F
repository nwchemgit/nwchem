#define NBLOCKS 2

*
* $Id: D3dB-mpi.F,v 1.9 2005-02-02 09:10:39 bylaska Exp $
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
      integer i,q
      integer proc_to,proc_from,it
      integer size,msgtype,source
c     integer msglen,type,nbytes
c     integer msgid
      integer msglen,type
c     integer status(MPI_STATUS_SIZE,512),request(512),reqcnt
      integer status(2),request(2),reqcnt
      integer taskid,np
      integer mpierr
      
      call Parallel_taskid(taskid)
      call Parallel_np(np)

*     **** allocate memory ****
      value = MA_push_get(mt_int,MPI_STATUS_SIZE*np,
     >                     'status',status(2),status(1))
      value = value.and.
     >        MA_push_get(mt_int,np,
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
            call MPI_IRECV(tmp2(int_mb(i2_start(1,nb)+it+1-1)),
     >                    msglen,MPI_DOUBLE_COMPLEX,
     >                    source,type,MPI_COMM_WORLD,
     >                    int_mb(request(1)+reqcnt-1),mpierr) 
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
            call MPI_SEND(tmp1(int_mb(i1_start(1,nb)+it+1-1)),
     >                     msglen,MPI_DOUBLE_COMPLEX,
     >                     proc_to,msgtype,MPI_COMM_WORLD,mpierr)
         end if

      end do

   
*     **** wait for completion of mp_send, also do a sync ****
      call MPI_WAITALL(reqcnt,int_mb(request(1)),
     >                        int_mb(status(1)),mpierr)

*     **** unpack A(i) array ****
#ifndef CRAY
!DIR$ ivdep
#endif
      do i=1,nfft3d(nb) !***(nx(nb)/2+1)*ny(nb)*nq(nb)
         A(i) = tmp2(int_mb(iq_to_i2(1,nb)+i-1))
      end do
      

      !*** deallocate memory ***
      value =           MA_pop_stack(request(2))
      value = value.and.MA_pop_stack(status(2))
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
      integer  rcv_len,rcv_proc
      integer i,c
      integer it
      integer source
      integer msglen,msgtype,type
      integer proc_from,proc_to
c     integer status(MPI_STATUS_SIZE,512),request(512),reqcnt
      integer status(2),request(2),reqcnt
      integer taskid,np
      integer index1,index2
      integer mpierr

      call Parallel_taskid(taskid)
      call Parallel_np(np)

*     **** allocate memory ****
      value = MA_push_get(mt_int,MPI_STATUS_SIZE*np,
     >                     'status',status(2),status(1))
      value = value.and.
     >        MA_push_get(mt_int,np,
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
            call MPI_IRECV(tmp2(int_mb(i2_start(1,nb)+it+1-1)),
     >                    msglen,MPI_DOUBLE_COMPLEX,
     >                    source,type,MPI_COMM_WORLD,
     >                    int_mb(request(1)+reqcnt-1),mpierr)
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
            call MPI_SEND(tmp1(int_mb(i1_start(1,nb)+it+1-1)),
     >                     msglen,MPI_DOUBLE_COMPLEX,
     >                     proc_to,msgtype,MPI_COMM_WORLD,mpierr)
         end if

      end do


*     **** wait for completion of mp_send, also do a sync ****
      call MPI_WAITALL(reqcnt,int_mb(request(1)),
     >                        int_mb(status(1)),mpierr)


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
      value = value.and.MA_pop_stack(status(2))
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
#include "mpif.h"
#include "D3dB.fh"


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
      integer i,q,nnfft3d
      integer proc_to,proc_from,it
      integer size,msgtype,source
c     integer msglen,type,nbytes
c     integer msgid
      integer msglen,type
c      integer status(MPI_STATUS_SIZE,512),request(512),reqcnt
      integer status(2),request(2),reqcnt
      integer taskid,np
      integer mpierr
      
      call Parallel_taskid(taskid)
      call Parallel_np(np)

*     **** allocate memory ****
      value = MA_push_get(mt_int,MPI_STATUS_SIZE*np,
     >                     'status',status(2),status(1))
      value = value.and.
     >        MA_push_get(mt_int,np,
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
#ifndef CRAY
!DIR$ ivdep
#endif
      do i=1,msglen
         tmp2(int_mb(h_i2_start(1,op,nb)+1-1)+i-1) 
     > = tmp1(int_mb(h_i1_start(1,op,nb)+1-1)+i-1)
      end do
         

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
            call MPI_IRECV(tmp2(int_mb(h_i2_start(1,op,nb)+it+1-1)),
     >                    msglen,MPI_DOUBLE_COMPLEX,
     >                    source,type,MPI_COMM_WORLD,
     >                    int_mb(request(1)+reqcnt-1),mpierr) 
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
            call MPI_SEND(tmp1(int_mb(h_i1_start(1,op,nb)+it+1-1)),
     >                     msglen,MPI_DOUBLE_COMPLEX,
     >                     proc_to,msgtype,MPI_COMM_WORLD,mpierr)
         end if

      end do

   
*     **** wait for completion of mp_send, also do a sync ****
      call MPI_WAITALL(reqcnt,int_mb(request(1)),
     >                        int_mb(status(1)),mpierr)

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
      value = value.and.MA_pop_stack(status(2))
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
      integer i,q,nnfft3d
      integer proc_to,proc_from,it
      integer size,msgtype,source
c     integer msglen,type,nbytes
c     integer msgid
      integer msglen,type
c     integer status(MPI_STATUS_SIZE,512),request(512),reqcnt
      integer status(2),request(2),reqcnt
      integer taskid,np
      integer mpierr
      
      call Parallel_taskid(taskid)
      call Parallel_np(np)


*     **** allocate memory ****
      value = MA_push_get(mt_int,MPI_STATUS_SIZE*np,
     >                     'status',status(2),status(1))
      value = value.and.
     >        MA_push_get(mt_int,np,
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
#ifndef CRAY
!DIR$ ivdep
#endif
      do i=1,msglen
         tmp2(int_mb(h_i2_start(1,op,nb)+1-1)+i-1) 
     > = tmp1(int_mb(h_i1_start(1,op,nb)+1-1)+i-1)
      end do
         

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
            call MPI_IRECV(tmp2(int_mb(h_i2_start(1,op,nb)+it+1-1)),
     >                    msglen,MPI_DOUBLE_PRECISION,
     >                    source,type,MPI_COMM_WORLD,
     >                    int_mb(request(1)+reqcnt-1),mpierr) 
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
            call MPI_SEND(tmp1(int_mb(h_i1_start(1,op,nb)+it+1-1)),
     >                     msglen,MPI_DOUBLE_PRECISION,
     >                     proc_to,msgtype,MPI_COMM_WORLD,mpierr)
         end if

      end do

   
*     **** wait for completion of mp_send, also do a sync ****
      call MPI_WAITALL(reqcnt,int_mb(request(1)),
     >                        int_mb(status(1)),mpierr)

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
      value = value.and.MA_pop_stack(status(2))
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


#include "mpif.h"

      integer MASTER
      parameter (MASTER=0)
      integer msglen,mpierr
      real*8 sumall,sumt

c     msglen = 8
      msglen = 1

c      sumt = sum

      call MPI_Allreduce(sum,sumall,msglen,MPI_DOUBLE_PRECISION,
     >                       MPI_SUM,MPI_COMM_WORLD,mpierr)
c      call MPI_BCAST(sumall,msglen,MPI_DOUBLE_PRECISION,
c     >               MASTER,MPI_COMM_WORLD,mpierr)

      sum = sumall
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
#include "mpif.h"
#include "errquit.fh"



      logical value
      integer MASTER
      parameter (MASTER=0)
      integer msglen

*     **** temporary workspace ****
      integer sumall(2),np


      call Parallel_np(np)
      call nwpw_timing_start(2)
      if (np.gt.1) then

*     ***** allocate temporary space ****
      value = MA_push_get(mt_dbl,n,'sumall',sumall(2),sumall(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)

      msglen = n


      call MPI_Allreduce(sum,dbl_mb(sumall(1)),msglen,
     >                MPI_DOUBLE_PRECISION,
     >                MPI_SUM,MPI_COMM_WORLD,mpierr)
c      call MPI_BCAST(dbl_mb(sumall(1)),msglen,
c     >               MPI_DOUBLE_PRECISION,
c     >               MASTER,MPI_COMM_WORLD,mpierr)



      call dcopy(n,dbl_mb(sumall(1)),1,sum,1)
      value = MA_pop_stack(sumall(2))

      end if
      call nwpw_timing_end(2)
      return
      end



#endif


