#define NBLOCKS 2


#ifdef TCGMSG
*     ***********************************
*     *					*
*     *	   C3dB_c_transpose_jk		*
*     *					*
*     ***********************************

      subroutine C3dB_c_transpose_jk(nb,A,tmp1,tmp2)

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
      integer     nb
      complex*16  A(*)
      complex*16  tmp1(*),tmp2(*)

#include "mafdecls.fh"
#include "C3dB.fh"


*     **** indexing variables ****
c     integer iq_to_i1((NFFT1/2+1)*NFFT2*NSLABS)
c     integer iq_to_i2((NFFT1/2+1)*NFFT2*NSLABS)
c     integer i1_start(NFFT3+1)
c     integer i2_start(NFFT3+1)
      integer iq_to_i1(2,NBLOCKS)
      integer iq_to_i2(2,NBLOCKS)
      integer i1_start(2,NBLOCKS)
      integer i2_start(2,NBLOCKS)
      common / c_trans_blk / iq_to_i1,iq_to_i2,i1_start,i2_start

*     **** Used to avoid asynchronous communications ****
      integer Nchannels(NBLOCKS)
      integer channel_proc(2,NBLOCKS)
      integer channel_type(2,NBLOCKS)
      common / c_channel_blk / channel_proc,channel_type,Nchannels

#include "tcgmsg.fh"
#include "msgtypesf.h"
      integer  rcv_len,rcv_proc

*     **** local variables ***
      integer i,c
      integer it
      integer source
      integer msglen
      integer pfrom,pto
      integer taskid,np
      
      call Parallel_taskid(taskid)
      call Parallel_np(np)

*     **** pack A(i) array ****
      do i=1,nx(nb)*ny(nb)*nq(nb)
         tmp1(int_mb(iq_to_i1(1,nb)+i-1)) = A(i)
      end do

*     **** it = 0, transpose data on same thread ****  
      msglen = int_mb(i2_start(1,nb)+2-1) - int_mb(i2_start(1,nb)+1-1)
      do i=1,msglen
         tmp2(int_mb(i2_start(1,nb)+1-1)+i-1) 
     > = tmp1(int_mb(i1_start(1,nb)+1-1)+i-1)
      end do
         

      do c=1,Nchannels(nb)
*        **** receive packed array data ****
         if (int_mb(channel_type(1,nb)+c-1) .eq. 1) then
            pfrom=int_mb(channel_proc(1,nb)+c-1)
            it = mod((taskid+np-pfrom),np)

            source=pfrom
            msglen = (int_mb(i2_start(1,nb)+it+2-1) 
     >             -  int_mb(i2_start(1,nb)+it+1-1))
         
            if (msglen.gt.0) then
               call RCV(9+MSGDBL,
     >                  tmp2(int_mb(i2_start(1,nb)+it+1-1)),
     >                  mdtob(2*msglen),rcv_len,
     >                  source,rcv_proc,1)
            end if
         end if

*        **** send packed array to other processors ****
         if (int_mb(channel_type(1,nb)+c-1) .eq. 0) then
            pto=int_mb(channel_proc(1,nb)+c-1)
            it = mod((pto-taskid+np),np)

            msglen    = (int_mb(i1_start(1,nb)+it+2-1)
     >                - int_mb(i1_start(1,nb)+it+1-1))
      
            if (msglen.gt.0) then
               call SND(9+MSGDBL,
     >                  tmp1(int_mb(i1_start(1,nb)+it+1-1)),
     >                  mdtob(2*msglen),pto,1)
            end if
         end if

      end do

   
*     **** unpack A(i) array ****
      do i=1,nx(nb)*ny(nb)*nq(nb)
         A(i) = tmp2(int_mb(iq_to_i2(1,nb)+i-1))
      end do
      
      return
      end


             
*     ***********************************
*     *					*
*     *	   C3dB_r_transpose_jk		*
*     *					*
*     ***********************************

      subroutine C3dB_r_transpose_jk(nb,A,tmp1,tmp2)

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
      integer     nb
      real*8      A(*)
      real*8      tmp1(*),tmp2(*)

#include "mafdecls.fh"
#include "C3dB.fh"


*     **** indexing variables ****
c     integer iq_to_i1((NFFT1/2+1)*NFFT2*NSLABS)
c     integer iq_to_i2((NFFT1/2+1)*NFFT2*NSLABS)
c     integer i1_start(NFFT3+1)
c     integer i2_start(NFFT3+1)
      integer iq_to_i1(2,NBLOCKS)
      integer iq_to_i2(2,NBLOCKS)
      integer i1_start(2,NBLOCKS)
      integer i2_start(2,NBLOCKS)
      common / c_trans_blk / iq_to_i1,iq_to_i2,i1_start,i2_start

*     **** Used to avoid asynchronous communications ****
      integer Nchannels(NBLOCKS)
      integer channel_proc(2,NBLOCKS)
      integer channel_type(2,NBLOCKS)
      common / c_channel_blk / channel_proc,channel_type,Nchannels

#include "tcgmsg.fh"
#include "msgtypesf.h"
      integer  rcv_len,rcv_proc

*     **** local variables ***
      integer i,c
      integer it
      integer source
      integer msglen
      integer pfrom,pto
      integer taskid,np
      
      call Parallel_taskid(taskid)
      call Parallel_np(np)

*     **** pack A(i) array ****
      do i=1,nx(nb)*ny(nb)*nq(nb)
         tmp1(int_mb(iq_to_i1(1,nb)+i-1)) = A(i)
      end do

*     **** it = 0, transpose data on same thread ****  
      msglen = int_mb(i2_start(1,nb)+2-1) - int_mb(i2_start(1,nb)+1-1)
      do i=1,msglen
         tmp2(int_mb(i2_start(1,nb)+1-1)+i-1) 
     > = tmp1(int_mb(i1_start(1,nb)+1-1)+i-1)
      end do
         

      do c=1,Nchannels(nb)
*        **** receive packed array data ****
         if (int_mb(channel_type(1,nb)+c-1) .eq. 1) then
            pfrom=int_mb(channel_proc(1,nb)+c-1)
            it = mod((taskid+np-pfrom),np)

            source=pfrom
            msglen = (int_mb(i2_start(1,nb)+it+2-1) 
     >             -  int_mb(i2_start(1,nb)+it+1-1))
         
            if (msglen.gt.0) then
               call RCV(9+MSGDBL,
     >                  tmp2(int_mb(i2_start(1,nb)+it+1-1)),
     >                  mdtob(msglen),rcv_len,
     >                  source,rcv_proc,1)
            end if
         end if

*        **** send packed array to other processors ****
         if (int_mb(channel_type(1,nb)+c-1) .eq. 0) then
            pto=int_mb(channel_proc(1,nb)+c-1)
            it = mod((pto-taskid+np),np)

            msglen    = (int_mb(i1_start(1,nb)+it+2-1)
     >                - int_mb(i1_start(1,nb)+it+1-1))
      
            if (msglen.gt.0) then
               call SND(9+MSGDBL,
     >                  tmp1(int_mb(i1_start(1,nb)+it+1-1)),
     >                  mdtob(msglen),pto,1)
            end if
         end if

      end do

   
*     **** unpack A(i) array ****
      do i=1,nx(nb)*ny(nb)*nq(nb)
         A(i) = tmp2(int_mb(iq_to_i2(1,nb)+i-1))
      end do
      
      return
      end
#endif





#ifdef TCGMSG
*     ***********************************
*     *					*
*     *	   C3dB_c_transpose_ijk		*
*     *					*
*     ***********************************

      subroutine C3dB_c_transpose_ijk(nb,op,A,tmp1,tmp2)

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
      integer     nb,op
      complex*16  A(*)
      complex*16  tmp1(*),tmp2(*)

#include "mafdecls.fh"

#include "C3dB.fh"


*     **** indexing variables ****
      integer h_iq_to_i1(2,6,NBLOCKS)
      integer h_iq_to_i2(2,6,NBLOCKS)
      integer h_i1_start(2,6,NBLOCKS)
      integer h_i2_start(2,6,NBLOCKS)
      common / c_trans_blk_ijk / h_iq_to_i1,
     >                         h_iq_to_i2,
     >                         h_i1_start,
     >                         h_i2_start

*     **** Used to avoid asynchronous communications ****
      integer Nchannels(NBLOCKS)
      integer channel_proc(2,NBLOCKS)
      integer channel_type(2,NBLOCKS)
      common / c_channel_blk / channel_proc,channel_type,Nchannels

#include "tcgmsg.fh"
#include "msgtypesf.h"
      integer  rcv_len,rcv_proc

*     **** local variables ***
      integer i,c,nnfft3d
      integer it
      integer source
      integer msglen
      integer pfrom,pto
      integer taskid,np
      
      call Parallel_taskid(taskid)
      call Parallel_np(np)


*     **** pack A(i) array ****
      if ((op.eq.1).or.(op.eq.5)) nnfft3d = (nx(nb))*nq1(nb)
      if ((op.eq.2).or.(op.eq.4)) nnfft3d = (ny(nb))*nq2(nb)
      if ((op.eq.3).or.(op.eq.6)) nnfft3d = (nz(nb))*nq3(nb)
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
         

      do c=1,Nchannels(nb)
*        **** receive packed array data ****
         if (int_mb(channel_type(1,nb)+c-1) .eq. 1) then
            pfrom=int_mb(channel_proc(1,nb)+c-1)
            it = mod((taskid+np-pfrom),np)

            source=pfrom
            msglen = (int_mb(h_i2_start(1,op,nb)+it+2-1) 
     >             -  int_mb(h_i2_start(1,op,nb)+it+1-1))
         
            if (msglen.gt.0) then
               call RCV(9+MSGDBL,
     >                  tmp2(int_mb(h_i2_start(1,op,nb)+it+1-1)),
     >                  mdtob(2*msglen),rcv_len,
     >                  source,rcv_proc,1)
            end if
         end if

*        **** send packed array to other processors ****
         if (int_mb(channel_type(1,nb)+c-1) .eq. 0) then
            pto=int_mb(channel_proc(1,nb)+c-1)
            it = mod((pto-taskid+np),np)

            msglen    = (int_mb(h_i1_start(1,op,nb)+it+2-1)
     >                -  int_mb(h_i1_start(1,op,nb)+it+1-1))
      
            if (msglen.gt.0) then
               call SND(9+MSGDBL,
     >                  tmp1(int_mb(h_i1_start(1,op,nb)+it+1-1)),
     >                  mdtob(2*msglen),pto,1)
            end if
         end if

      end do


*     **** unpack A(i) array ****
      if ((op.eq.4).or.(op.eq.6)) nnfft3d = (nx(nb))*nq1(nb)
      if ((op.eq.1).or.(op.eq.3)) nnfft3d = (ny(nb))*nq2(nb)
      if ((op.eq.2).or.(op.eq.5)) nnfft3d = (nz(nb))*nq3(nb)
#ifndef CRAY
!DIR$ ivdep
#endif
      do i=1,nnfft3d
         A(i) = tmp2(int_mb(h_iq_to_i2(1,op,nb)+i-1))
      end do

      
      return
      end



*     ***********************************
*     *					*
*     *	   C3dB_r_transpose_ijk		*
*     *					*
*     ***********************************

      subroutine C3dB_r_transpose_ijk(nb,op,A,tmp1,tmp2)

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

#include "C3dB.fh"


*     **** indexing variables ****
      integer h_iq_to_i1(2,6,NBLOCKS)
      integer h_iq_to_i2(2,6,NBLOCKS)
      integer h_i1_start(2,6,NBLOCKS)
      integer h_i2_start(2,6,NBLOCKS)
      common / c_trans_blk_ijk / h_iq_to_i1,
     >                         h_iq_to_i2,
     >                         h_i1_start,
     >                         h_i2_start

*     **** Used to avoid asynchronous communications ****
      integer Nchannels(NBLOCKS)
      integer channel_proc(2,NBLOCKS)
      integer channel_type(2,NBLOCKS)
      common / c_channel_blk / channel_proc,channel_type,Nchannels

#include "tcgmsg.fh"
#include "msgtypesf.h"
      integer  rcv_len,rcv_proc

*     **** local variables ***
      integer i,c,nnfft3d
      integer it
      integer source
      integer msglen
      integer pfrom,pto
      integer taskid,np
      
      call Parallel_taskid(taskid)
      call Parallel_np(np)


*     **** pack A(i) array ****
      if ((op.eq.1).or.(op.eq.5)) nnfft3d = (nx(nb))*nq1(nb)
      if ((op.eq.2).or.(op.eq.4)) nnfft3d = (ny(nb))*nq2(nb)
      if ((op.eq.3).or.(op.eq.6)) nnfft3d = (nz(nb))*nq3(nb)
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
         

      do c=1,Nchannels(nb)
*        **** receive packed array data ****
         if (int_mb(channel_type(1,nb)+c-1) .eq. 1) then
            pfrom=int_mb(channel_proc(1,nb)+c-1)
            it = mod((taskid+np-pfrom),np)

            source=pfrom
            msglen = (int_mb(h_i2_start(1,op,nb)+it+2-1) 
     >             -  int_mb(h_i2_start(1,op,nb)+it+1-1))
         
            if (msglen.gt.0) then
               call RCV(9+MSGDBL,
     >                  tmp2(int_mb(h_i2_start(1,op,nb)+it+1-1)),
     >                  mdtob(msglen),rcv_len,
     >                  source,rcv_proc,1)
            end if
         end if

*        **** send packed array to other processors ****
         if (int_mb(channel_type(1,nb)+c-1) .eq. 0) then
            pto=int_mb(channel_proc(1,nb)+c-1)
            it = mod((pto-taskid+np),np)

            msglen    = (int_mb(h_i1_start(1,op,nb)+it+2-1)
     >                -  int_mb(h_i1_start(1,op,nb)+it+1-1))
      
            if (msglen.gt.0) then
               call SND(9+MSGDBL,
     >                  tmp1(int_mb(h_i1_start(1,op,nb)+it+1-1)),
     >                  mdtob(msglen),pto,1)
            end if
         end if

      end do


*     **** unpack A(i) array ****
      if ((op.eq.4).or.(op.eq.6)) nnfft3d = (nx(nb))*nq1(nb)
      if ((op.eq.1).or.(op.eq.3)) nnfft3d = (ny(nb))*nq2(nb)
      if ((op.eq.2).or.(op.eq.5)) nnfft3d = (nz(nb))*nq3(nb)
#ifndef CRAY
!DIR$ ivdep
#endif
      do i=1,nnfft3d
         A(i) = tmp2(int_mb(h_iq_to_i2(1,op,nb)+i-1))
      end do

      
      return
      end


#endif

