#define NBLOCKS 2

*
* $Id: D3dB-mpi.F,v 1.3 2004-04-21 15:09:47 bylaska Exp $
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


#include "mpif.h"
      integer mpierr

*     **** local variables ***
      integer i,q
      integer proc_to,proc_from,it
      integer size,msgtype,source
c     integer msglen,type,nbytes
c     integer msgid
      integer msglen,type
      integer status(MPI_STATUS_SIZE,512),request(512),reqcnt
      integer taskid,np
      
      call Parallel_taskid(taskid)
      call Parallel_np(np)

*     **** pack A(i) array ****
!DIR$ ivdep
      do i=1,nfft3d(nb) 
         tmp1(int_mb(iq_to_i1(1,nb)+i-1)) = A(i)
      end do

*     **** it = 0, transpose data on same thread ****  
      msglen = int_mb(i2_start(1,nb)+2-1) - int_mb(i2_start(1,nb)+1-1)
!DIR$ ivdep
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
     >                    request(reqcnt),mpierr) 
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
      call MPI_WAITALL(reqcnt,request,status,mpierr)

*     **** unpack A(i) array ****
!DIR$ ivdep
      do i=1,nfft3d(nb) !***(nx(nb)/2+1)*ny(nb)*nq(nb)
         A(i) = tmp2(int_mb(iq_to_i2(1,nb)+i-1))
      end do
      
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


#include "mpif.h"
         integer mpierr

*     **** local variables ***
      integer  rcv_len,rcv_proc
      integer i,c
      integer it
      integer source
      integer msglen,msgtype,type
      integer proc_from,proc_to
      integer status(MPI_STATUS_SIZE,512),request(512),reqcnt
      integer taskid,np
      integer index1,index2

      call Parallel_taskid(taskid)
      call Parallel_np(np)

*     **** pack A(i) array ****
!DIR$ ivdep
      do index1=int_mb(i1_start(1,nb)+1-1),
     >         (int_mb(i1_start(1,nb)+np+1-1)-1)
         tmp1(index1) = A(int_mb(iq_to_i1(1,nb)+index1-1))
      end do

*     **** it = 0, transpose data on same thread ****
      msglen = int_mb(i2_start(1,nb)+2-1) - int_mb(i2_start(1,nb)+1-1)
!DIR$ ivdep
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
     >                    request(reqcnt),mpierr)
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
      call MPI_WAITALL(reqcnt,request,status,mpierr)


*     **** unpack A(i) array ****
!DIR$ ivdep
      do index2=int_mb(i2_start(1,nb)+1-1),
     >         (int_mb(i2_start(1,nb)+np+1-1)-1)
        A(int_mb(iq_to_i2(1,nb)+index2-1))=dconjg(tmp2(index2))
      end do

      return
      end



#endif


