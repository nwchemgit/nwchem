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
      integer     nb
      complex*16  A(*)
      complex*16  tmp1(*),tmp2(*)

#include "bafdecls.fh"

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

*     **** Used to avoid asynchronous communications ****
      integer Nchannels(NBLOCKS)
      integer channel_proc(2,NBLOCKS)
      integer channel_type(2,NBLOCKS)
      common / channel_blk / channel_proc,channel_type,Nchannels

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

*     **** external functions ****
      integer  Parallel2d_convert_taskid_i
      external Parallel2d_convert_taskid_i
      
      call Parallel2d_taskid_i(taskid)
      call Parallel2d_np_i(np)

*     **** pack A(i) array ****
#ifndef CRAY
!DIR$ ivdep
#endif
      do i=1,nfft3d(nb)  !(nx(nb)/2+1)*ny(nb)*nq(nb)
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
     >                  Parallel2d_convert_taskid_i(source),
     >                  rcv_proc,1)
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
     >                  mdtob(2*msglen),
     >                  Parallel2d_convert_taskid_i(pto),
     >                  1)
            end if
         end if

      end do

   
*     **** unpack A(i) array ****
#ifndef CRAY
!DIR$ ivdep
#endif
      do i=1,nfft3d(nb)  !(nx(nb)/2+1)*ny(nb)*nq(nb)
         A(i) = tmp2(int_mb(iq_to_i2(1,nb)+i-1))
      end do

      
      return
      end

c*     ***********************************
c*     *					*
c*     *	   D3dB_nc_transpose_jk		*
c*     *					*
c*     ***********************************
c
c      subroutine D3dB_nc_transpose_jk(nb,ne,A,tmp1,tmp2)
c
c*****************************************************
c*                                                   *
c*      This routine performs the operation          *
c*               A(i,k,j,n) <- A(i,j,k,n)                * 
c*                                                   *
c*      np = the number of worker nodes              *
c*      proc#=0...(np-1)
c*                                                   *
c*       this transpose uses more buffer space       *
c*       then transpose2                             *
c*****************************************************
c      implicit none
c      integer     nb,ne
c      complex*16  A(*)
c      complex*16  tmp1(*),tmp2(*)
c
c#include "bafdecls.fh"
c
c#include "D3dB.fh"
c
c
c*     **** indexing variables ****
cc     integer iq_to_i1((NFFT1/2+1)*NFFT2*NSLABS)
cc     integer iq_to_i2((NFFT1/2+1)*NFFT2*NSLABS)
cc     integer i1_start(NFFT3+1)
cc     integer i2_start(NFFT3+1)
c      integer iq_to_i1(2,NBLOCKS)
c      integer iq_to_i2(2,NBLOCKS)
c      integer i1_start(2,NBLOCKS)
c      integer i2_start(2,NBLOCKS)
c      common / trans_blk / iq_to_i1,iq_to_i2,i1_start,i2_start
c
c*     **** Used to avoid asynchronous communications ****
c      integer Nchannels(NBLOCKS)
c      integer channel_proc(2,NBLOCKS)
c      integer channel_type(2,NBLOCKS)
c      common / channel_blk / channel_proc,channel_type,Nchannels
c
c#include "tcgmsg.fh"
c#include "msgtypesf.h"
c      integer  rcv_len,rcv_proc
c
c*     **** local variables ***
c      integer i,c,n,nnfft3d
c      integer it
c      integer source
c      integer msglen
c      integer pfrom,pto
c      integer taskid,np
c      
c      call Parallel_taskid(taskid)
c      call Parallel_np(np)
c
c      !nnfft3d = (nx(nb)/2+1)*ny(nb)*nq(nb)
c
c*     **** pack A(i) array ****
c      do i=1,nfft3d(nb)
c#ifndef CRAY
c!DIR$ ivdep
c#endif
c      do n=1,ne
c         tmp1(n+(int_mb(iq_to_i1(1,nb)+i-1)-1)*ne) 
c     >   = A(i+(n-1)*nfft3d(nb))
c      end do
c      end do
c
c*     **** it = 0, transpose data on same thread ****  
c      msglen = int_mb(i2_start(1,nb)+2-1) - int_mb(i2_start(1,nb)+1-1)
c      do i=1,msglen
c#ifndef CRAY
c!DIR$ ivdep
c#endif
c      do n=1,ne
c         tmp2(n+(int_mb(i2_start(1,nb))+i-2)*ne)
c     > = tmp1(n+(int_mb(i1_start(1,nb))+i-2)*ne)
c      end do
c      end do
c         
c
c      do c=1,Nchannels(nb)
c*        **** receive packed array data ****
c         if (int_mb(channel_type(1,nb)+c-1) .eq. 1) then
c            pfrom=int_mb(channel_proc(1,nb)+c-1)
c            it = mod((taskid+np-pfrom),np)
c
c            source=pfrom
c            msglen = (int_mb(i2_start(1,nb)+it+2-1) 
c     >             -  int_mb(i2_start(1,nb)+it+1-1))*ne
c         
c            if (msglen.gt.0) then
c               call RCV(9+MSGDBL,
c     >                  tmp2(1+(int_mb(i2_start(1,nb)+it+1-1)-1)*ne),
c     >                  mdtob(2*msglen),rcv_len,
c     >                  source,rcv_proc,1)
c            end if
c         end if
c
c*        **** send packed array to other processors ****
c         if (int_mb(channel_type(1,nb)+c-1) .eq. 0) then
c            pto=int_mb(channel_proc(1,nb)+c-1)
c            it = mod((pto-taskid+np),np)
c
c            msglen    = (int_mb(i1_start(1,nb)+it+2-1)
c     >                - int_mb(i1_start(1,nb)+it+1-1))*ne
c      
c            if (msglen.gt.0) then
c               call SND(9+MSGDBL,
c     >                  tmp1(1+(int_mb(i1_start(1,nb)+it+1-1)-1)*ne),
c     >                  mdtob(2*msglen),pto,1)
c            end if
c         end if
c
c      end do
c
c   
c*     **** unpack A(i) array ****
c      do i=1,nfft3d(nb)
c#ifndef CRAY
c!DIR$ ivdep
c#endif
c      do n=1,ne
c         A(i+(n-1)*nfft3d(nb)) 
c     >   = tmp2(n+(int_mb(iq_to_i2(1,nb)+i-1)-1)*ne)
c      end do
c      end do
c
c      
c      return
c      end
c




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

*     **** Used to avoid asynchronous communications ****
      integer Nchannels(NBLOCKS)
      integer channel_proc(2,NBLOCKS)
      integer channel_type(2,NBLOCKS)
      common / channel_blk / channel_proc,channel_type,Nchannels

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
      integer index1,index2

*     **** external functions ****
      integer  Parallel2d_convert_taskid_i
      external Parallel2d_convert_taskid_i

      call Parallel2d_taskid_i(taskid)
      call Parallel2d_np_i(np)

*     **** pack A(i) array ****
      do index1=int_mb(i1_start(1,nb)+1-1),
     >         (int_mb(i1_start(1,nb)+np+1-1)-1)
         tmp1(index1) = A(int_mb(iq_to_i1(1,nb)+index1-1))
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
     >                  Parallel2d_convert_taskid_i(source),
     >                  rcv_proc,1)
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
     >                  mdtob(2*msglen),
     >                  Parallel2d_convert_taskid_i(pto),
     >                  1)
            end if
         end if

      end do

*     **** unpack A(i) array ****
      do index2=int_mb(i2_start(1,nb)+1-1),
     >         (int_mb(i2_start(1,nb)+np+1-1)-1)
        A(int_mb(iq_to_i2(1,nb)+index2-1))=dconjg(tmp2(index2))
      end do

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
      integer     nb,op
      complex*16  A(*)
      complex*16  tmp1(*),tmp2(*)

#include "bafdecls.fh"

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

*     **** Used to avoid asynchronous communications ****
      integer Nchannels(NBLOCKS)
      integer channel_proc(2,NBLOCKS)
      integer channel_type(2,NBLOCKS)
      common / channel_blk / channel_proc,channel_type,Nchannels

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

*     **** external functions ****
      integer  Parallel2d_convert_taskid_i
      external Parallel2d_convert_taskid_i
      
      call Parallel2d_taskid_i(taskid)
      call Parallel2d_np_i(np)


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
     >                  Parallel2d_convert_taskid_i(source),
     >                  rcv_proc,1)
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
     >                  mdtob(2*msglen),
     >                  Parallel2d_convert_taskid_i(pto),1)
            end if
         end if

      end do


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

*     **** Used to avoid asynchronous communications ****
      integer Nchannels(NBLOCKS)
      integer channel_proc(2,NBLOCKS)
      integer channel_type(2,NBLOCKS)
      common / channel_blk / channel_proc,channel_type,Nchannels

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

*     **** external functions ****
      integer  Parallel2d_convert_taskid_i
      external Parallel2d_convert_taskid_i
      
      call Parallel2d_taskid_i(taskid)
      call Parallel2d_np_i(np)


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
     >                  Parallel2d_convert_taskid_i(source),
     >                  rcv_proc,1)
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
     >                  mdtob(msglen),
     >                  Parallel2d_convert_taskid_i(pto),
     >                  1)
            end if
         end if

      end do


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

#include "tcgmsg.fh"
#include "msgtypesf.h"

*     **** local variables ****
      integer np_i

*     **** external functions ****
      integer  Parallel2d_comm_i
      external Parallel2d_comm_i

      call Parallel2d_np_i(np_i)
      if (np_i.gt.1) then
         call GA_PGROUP_DGOP(Parallel2d_comm_i(),
     >                       9+MSGDBL,sum,1,'+')
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
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "errquit.fh"


*     **** temporary workspace ****
      integer np_i

*     **** external functions ****
      integer  Parallel2d_comm_i
      external Parallel2d_comm_i

      call nwpw_timing_start(2)

      call Parallel2d_np_i(np_i)
      if (np_i.gt.1) then
         call GA_PGROUP_DGOP(Parallel2d_comm_i(),
     >                       9+MSGDBL,sum,n,'+')
      end if

      call nwpw_timing_end(2)

      return
      end




*     ***********************************
*     *                                 *
*     *         D3dB_ISumAll            *
*     *                                 *
*     ***********************************

      subroutine D3dB_ISumAll(sum)
c     implicit none
      integer  sum

#include "tcgmsg.fh"
#include "msgtypesf.h"

*     **** local variables ****
      integer np_i

*     **** external functions ****
      integer  Parallel2d_comm_i
      external Parallel2d_comm_i


      call Parallel2d_np_i(np_i)
      if (np_i.gt.1) then
         call GA_PGROUP_IGOP(Parallel2d_comm_i(),
     >                       9+MSGINT,sum,1,'+')
      end if

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
#include "tcgmsg.fh"
#include "msgtypesf.h"

*     **** local variables ****
      integer np_i

*     **** external functions ****
      integer  Parallel2d_comm_i
      external Parallel2d_comm_i

      call nwpw_timing_start(2)

      call Parallel2d_np_i(np_i)
      if (np_i.gt.1) then
         call GA_PGROUP_IGOP(Parallel2d_comm_i(),
     >                       9+MSGINT,sum,n,'+')
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
      integer     nbb
      complex*16  A(*)
      complex*16  tmp1(*),tmp2(*)

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



*     **** Used to avoid asynchronous communications ****
      integer Nchannels(NBLOCKS)
      integer channel_proc(2,NBLOCKS)
      integer channel_type(2,NBLOCKS)
      common / channel_blk / channel_proc,channel_type,Nchannels

#include "tcgmsg.fh"
#include "msgtypesf.h"
      integer  rcv_len,rcv_proc

*     **** local variables ***
      integer c,it
      integer source
      integer msglen
      integer pfrom,pto
      integer taskid,np
      integer n1,n2,id

*     **** external functions ****
      integer  Parallel2d_convert_taskid_i
      external Parallel2d_convert_taskid_i
      

      id = 2*(nbb/2) + 1

      call Parallel2d_taskid_i(taskid)
      call Parallel2d_np_i(np)


      n1 = int_mb(i1_start(1,nbb)+np) - 1
      n2 = int_mb(i2_start(1,nbb)+np) - 1

*     **** pack A(i) array ****
      call D3dB_pfft_index1_copy(n1,int_mb(iq_to_i1(1,nbb)),A,tmp1)


*     **** it = 0, transpose data on same thread ****  
      msglen = int_mb(i2_start(1,nbb)+2-1) - int_mb(i2_start(1,nbb)+1-1)
      call dcopy(2*msglen,
     >           tmp1(int_mb(i1_start(1,nbb)+1-1)),1,
     >           tmp2(int_mb(i2_start(1,nbb)+1-1)),1)


      do c=1,Nchannels(1)
*        **** receive packed array data ****
         if (int_mb(channel_type(1,1)+c-1) .eq. 1) then
            pfrom=int_mb(channel_proc(1,1)+c-1)
            it = mod((taskid+np-pfrom),np)

            source=pfrom
            msglen = (int_mb(i2_start(1,nbb)+it+2-1) 
     >             -  int_mb(i2_start(1,nbb)+it+1-1))
         
            if (msglen.gt.0) then
               call RCV(9+MSGDBL,
     >                  tmp2(int_mb(i2_start(1,nbb)+it+1-1)),
     >                  mdtob(2*msglen),rcv_len,
     >                  Parallel2d_convert_taskid_i(source),
     >                  rcv_proc,1)
            end if
         end if

*        **** send packed array to other processors ****
         if (int_mb(channel_type(1,1)+c-1) .eq. 0) then
            pto=int_mb(channel_proc(1,1)+c-1)
            it = mod((pto-taskid+np),np)

            msglen    = (int_mb(i1_start(1,nbb)+it+2-1)
     >                - int_mb(i1_start(1,nbb)+it+1-1))
      
            if (msglen.gt.0) then
               call SND(9+MSGDBL,
     >                  tmp1(int_mb(i1_start(1,nbb)+it+1-1)),
     >                  mdtob(2*msglen),
     >                  Parallel2d_convert_taskid_i(pto),
     >                  1)
            end if
         end if

      end do

   
*     **** unpack A(i) array ****
      call D3dB_pfft_index2_copy(n2,int_mb(iq_to_i2(1,nbb)),tmp2,A)
      call D3dB_pfft_index2_zero(nfft3d(id)-n2,
     >                           int_mb(iz_to_i2(1,nbb)),A)

      
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
      integer     nbb
      complex*16  A(*)
      complex*16  tmp1(*),tmp2(*)

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



*     **** Used to avoid asynchronous communications ****
      integer Nchannels(NBLOCKS)
      integer channel_proc(2,NBLOCKS)
      integer channel_type(2,NBLOCKS)
      common / channel_blk / channel_proc,channel_type,Nchannels

#include "tcgmsg.fh"
#include "msgtypesf.h"
      integer  rcv_len,rcv_proc

*     **** local variables ***
      integer c,it
      integer source
      integer msglen
      integer pfrom,pto
      integer taskid,np
      integer n1,n2,id

*     **** external functions ****
      integer  Parallel2d_convert_taskid_i
      external Parallel2d_convert_taskid_i
      
      id = 2*(nbb/2) + 1

      call Parallel2d_taskid_i(taskid)
      call Parallel2d_np_i(np)


      n1 = int_mb(i1_start(1,nbb)+np) - 1
      n2 = int_mb(i2_start(1,nbb)+np) - 1

*     **** pack A(i) array ****
      call D3dB_pfft_index1_copy(n1,int_mb(iq_to_i1(1,nbb)),A,tmp1)


*     **** it = 0, transpose data on same thread ****  
      msglen = int_mb(i2_start(1,nbb)+2-1) - int_mb(i2_start(1,nbb)+1-1)
      call dcopy(2*msglen,
     >           tmp1(int_mb(i1_start(1,nbb)+1-1)),1,
     >           tmp2(int_mb(i2_start(1,nbb)+1-1)),1)


      do c=1,Nchannels(1)
*        **** receive packed array data ****
         if (int_mb(channel_type(1,1)+c-1) .eq. 1) then
            pfrom=int_mb(channel_proc(1,1)+c-1)
            it = mod((taskid+np-pfrom),np)

            source=pfrom
            msglen = (int_mb(i2_start(1,nbb)+it+2-1) 
     >             -  int_mb(i2_start(1,nbb)+it+1-1))
         
            if (msglen.gt.0) then
               call RCV(9+MSGDBL,
     >                  tmp2(int_mb(i2_start(1,nbb)+it+1-1)),
     >                  mdtob(2*msglen),rcv_len,
     >                  Parallel2d_convert_taskid_i(source),
     >                  rcv_proc,1)
            end if
         end if

*        **** send packed array to other processors ****
         if (int_mb(channel_type(1,1)+c-1) .eq. 0) then
            pto=int_mb(channel_proc(1,1)+c-1)
            it = mod((pto-taskid+np),np)

            msglen    = (int_mb(i1_start(1,nbb)+it+2-1)
     >                - int_mb(i1_start(1,nbb)+it+1-1))
      
            if (msglen.gt.0) then
               call SND(9+MSGDBL,
     >                  tmp1(int_mb(i1_start(1,nbb)+it+1-1)),
     >                  mdtob(2*msglen),
     >                  Parallel2d_convert_taskid_i(pto),
     >                  1)
            end if
         end if

      end do

*     **** unpack A(i) array ****
      call D3dB_pfft_index2_copy(n2,int_mb(iq_to_i2(1,nbb)),tmp2,A)
      call D3dB_pfft_index2_zero(nfft3d(id)-n2,
     >                           int_mb(iz_to_i2(1,nbb)),A)
      
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
      subroutine D3dB_c_timereverse_start(nb,A,tmp1,tmp2,request,reqcnt)
      implicit none
      integer     nb
      complex*16  A(*)
      complex*16  tmp1(*)
      complex*16  tmp2(*)
      integer     request(*),reqcnt

#include "bafdecls.fh"
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

*     **** Used to avoid asynchronous communications ****
      integer Nchannels(NBLOCKS)
      integer channel_proc(2,NBLOCKS)
      integer channel_type(2,NBLOCKS)
      common / channel_blk / channel_proc,channel_type,Nchannels

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
      integer index1

*     **** external functions ****
      integer  Parallel2d_convert_taskid_i
      external Parallel2d_convert_taskid_i


      call Parallel2d_taskid_i(taskid)
      call Parallel2d_np_i(np)

*     **** pack A(i) array ****
      do index1=int_mb(i1_start(1,nb)+1-1),
     >         (int_mb(i1_start(1,nb)+np+1-1)-1)
         tmp1(index1) = A(int_mb(iq_to_i1(1,nb)+index1-1))
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
     >                  Parallel2d_convert_taskid_i(source),
     >                  rcv_proc,1)
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
     >                  mdtob(2*msglen),
     >                  Parallel2d_convert_taskid_i(pto),1)
            end if
         end if

      end do


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
#include "errquit.fh"


*     **** indexing variables ****
      integer iq_to_i1(2,NBLOCKS)
      integer iq_to_i2(2,NBLOCKS)
      integer i1_start(2,NBLOCKS)
      integer i2_start(2,NBLOCKS)
      common / timereverse_blk / iq_to_i1,iq_to_i2,i1_start,i2_start


*     **** local variables ***
      integer np
      integer index2


      call Parallel2d_np_i(np)

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


*     ************************************
*     *                                  *
*     *         Balance_c_balance_start  *
*     *                                  *
*     ************************************

      subroutine Balance_c_balance_start(nb,A,request,reqcnt,msgtype)
      implicit none
      integer nb
      complex*16 A(*)
      integer    request(*),reqcnt,msgtype

      call Balance_c_balance(nb,A)
      return
      end

*     ************************************
*     *                                  *
*     *         Balance_c_balance_end    *
*     *                                  *
*     ************************************
* 
      subroutine Balance_c_balance_end(nb,A,request,reqcnt)
      implicit none
      integer nb
      complex*16 A(*)
      integer    request(*),reqcnt

*     *** dummy routine ***
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
      integer    request(*),reqcnt,msgtype

      call Balance_c_unbalance(nb,A)
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

*     *** dummy routine ***
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
      integer request(*),reqcnt,msgtype

#include "bafdecls.fh"
#include "errquit.fh"
#include "D3dB.fh"

*     **** indexing variables ****
      integer iq_to_i1(2,0:3)
      integer iq_to_i2(2,0:3)
      integer iz_to_i2(2,0:3)
      integer i1_start(2,0:3)
      integer i2_start(2,0:3)
      common / ptrans_blk1 / iq_to_i1,iq_to_i2,iz_to_i2,
     >                       i1_start,i2_start

*     **** Used to avoid asynchronous communications ****
      integer Nchannels(NBLOCKS)
      integer channel_proc(2,NBLOCKS)
      integer channel_type(2,NBLOCKS)
      common / channel_blk / channel_proc,channel_type,Nchannels

#include "tcgmsg.fh"
#include "msgtypesf.h"
      integer  rcv_len,rcv_proc


*     **** local variables ***
      integer c,it
      integer source
      integer msglen
      integer pfrom,pto
      integer taskid,np
      integer n1

*     **** external functions ****
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
         

      do c=1,Nchannels(1)
*        **** receive packed array data ****
         if (int_mb(channel_type(1,1)+c-1) .eq. 1) then
            pfrom=int_mb(channel_proc(1,1)+c-1)
            it = mod((taskid+np-pfrom),np)

            source=pfrom
            msglen = (int_mb(i2_start(1,nbb)+it+2-1)
     >             -  int_mb(i2_start(1,nbb)+it+1-1))

            if (msglen.gt.0) then
               call RCV(9+MSGDBL,
     >                  tmp2(int_mb(i2_start(1,nbb)+it+1-1)),
     >                  mdtob(2*msglen),rcv_len,
     >                  Parallel2d_convert_taskid_i(source),
     >                  rcv_proc,1)
            end if
         end if

*        **** send packed array to other processors ****
         if (int_mb(channel_type(1,1)+c-1) .eq. 0) then
            pto=int_mb(channel_proc(1,1)+c-1)
            it = mod((pto-taskid+np),np)

            msglen    = (int_mb(i1_start(1,nbb)+it+2-1)
     >                - int_mb(i1_start(1,nbb)+it+1-1))

            if (msglen.gt.0) then
               call SND(9+MSGDBL,
     >                  tmp1(int_mb(i1_start(1,nbb)+it+1-1)),
     >                  mdtob(2*msglen),
     >                  Parallel2d_convert_taskid_i(pto),1)
            end if
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
#include "errquit.fh"
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
      
      id = 2*(nbb/2) + 1

      call Parallel2d_np_i(np)

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
      integer request(*),reqcnt,msgtype

#include "bafdecls.fh"
#include "errquit.fh"
#include "D3dB.fh"


*     **** indexing variables ****
      integer iq_to_i1(2,0:3)
      integer iq_to_i2(2,0:3)
      integer iz_to_i2(2,0:3)
      integer i1_start(2,0:3)
      integer i2_start(2,0:3)
      common / ptrans_blk2 / iq_to_i1,iq_to_i2,iz_to_i2,
     >                       i1_start,i2_start



*     **** Used to avoid asynchronous communications ****
      integer Nchannels(NBLOCKS)
      integer channel_proc(2,NBLOCKS)
      integer channel_type(2,NBLOCKS)
      common / channel_blk / channel_proc,channel_type,Nchannels

#include "tcgmsg.fh"
#include "msgtypesf.h"
      integer  rcv_len,rcv_proc


*     **** local variables ***
      integer c,it
      integer source
      integer msglen
      integer pfrom,pto
      integer taskid,np
      integer n1

*     **** external functions ****
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
         

      do c=1,Nchannels(1)
*        **** receive packed array data ****
         if (int_mb(channel_type(1,1)+c-1) .eq. 1) then
            pfrom=int_mb(channel_proc(1,1)+c-1)
            it = mod((taskid+np-pfrom),np)

            source=pfrom
            msglen = (int_mb(i2_start(1,nbb)+it+2-1)
     >             -  int_mb(i2_start(1,nbb)+it+1-1))

            if (msglen.gt.0) then
               call RCV(9+MSGDBL,
     >                  tmp2(int_mb(i2_start(1,nbb)+it+1-1)),
     >                  mdtob(2*msglen),rcv_len,
     >                  Parallel2d_convert_taskid_i(source),
     >                  rcv_proc,1)
            end if
         end if

*        **** send packed array to other processors ****
         if (int_mb(channel_type(1,1)+c-1) .eq. 0) then
            pto=int_mb(channel_proc(1,1)+c-1)
            it = mod((pto-taskid+np),np)

            msglen    = (int_mb(i1_start(1,nbb)+it+2-1)
     >                - int_mb(i1_start(1,nbb)+it+1-1))

            if (msglen.gt.0) then
               call SND(9+MSGDBL,
     >                  tmp1(int_mb(i1_start(1,nbb)+it+1-1)),
     >                  mdtob(2*msglen),
     >                  Parallel2d_convert_taskid_i(pto),
     >                  1)
            end if
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
#include "errquit.fh"
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
      integer np,n2,id
      
      id = 2*(nbb/2) + 1

      call Parallel2d_np_i(np)

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
      integer     request(*),reqcnt,msgtype

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

*     **** Used to avoid asynchronous communications ****
      integer Nchannels(NBLOCKS)
      integer channel_proc(2,NBLOCKS)
      integer channel_type(2,NBLOCKS)
      common / channel_blk / channel_proc,channel_type,Nchannels

#include "tcgmsg.fh"
#include "msgtypesf.h"
      integer  rcv_len,rcv_proc


*     **** local variables ***
      integer c,n1
      integer it
      integer source
      integer msglen
      integer pfrom,pto
      integer taskid,np

*     **** external functions ****
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



      do c=1,Nchannels(1)
*        **** receive packed array data ****
         if (int_mb(channel_type(1,1)+c-1) .eq. 1) then
            pfrom=int_mb(channel_proc(1,1)+c-1)
            it = mod((taskid+np-pfrom),np)

            source=pfrom
            msglen = (int_mb(h_i2_start(1,op,nbb)+it+2-1)
     >             -  int_mb(h_i2_start(1,op,nbb)+it+1-1))

            if (msglen.gt.0) then
               call RCV(9+MSGDBL,
     >                  tmp2(int_mb(h_i2_start(1,op,nbb)+it+1-1)),
     >                  mdtob(2*msglen),rcv_len,
     >                  Parallel2d_convert_taskid_i(source),
     >                  rcv_proc,1)
            end if
         end if

*        **** send packed array to other processors ****
         if (int_mb(channel_type(1,1)+c-1) .eq. 0) then
            pto=int_mb(channel_proc(1,1)+c-1)
            it = mod((pto-taskid+np),np)

            msglen    = (int_mb(h_i1_start(1,op,nbb)+it+2-1)
     >                -  int_mb(h_i1_start(1,op,nbb)+it+1-1))

            if (msglen.gt.0) then
               call SND(9+MSGDBL,
     >                  tmp1(int_mb(h_i1_start(1,op,nbb)+it+1-1)),
     >                  mdtob(2*msglen),
     >                  Parallel2d_convert_taskid_i(pto),
     >                  1)
            end if
         end if

      end do

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
      integer n2,n3,np

      
      call Parallel2d_np_i(np)

*     **** unpack A(i) array ****
      n2 = int_mb(h_i2_start(1,op,nbb)+np) - 1
      n3 = h_iz_to_i2_count(op,nbb)
      call D3dB_pfft_index2_copy(n2,int_mb(h_iq_to_i2(1,op,nbb)),tmp2,A)
      call D3dB_pfft_index2_zero(n3,int_mb(h_iz_to_i2(1,op,nbb)),A)


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

*     **** Used to avoid asynchronous communications ****
      integer Nchannels(NBLOCKS)
      integer channel_proc(2,NBLOCKS)
      integer channel_type(2,NBLOCKS)
      common / channel_blk / channel_proc,channel_type,Nchannels

#include "tcgmsg.fh"
#include "msgtypesf.h"
      integer  rcv_len,rcv_proc


*     **** local variables ***
      integer c,n1,n2,n3
      integer it
      integer source
      integer msglen
      integer pfrom,pto
      integer taskid,np

*     **** external functions ****
      integer  Parallel2d_convert_taskid_i
      external Parallel2d_convert_taskid_i


      call Parallel2d_taskid_i(taskid)
      call Parallel2d_np_i(np)

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


      do c=1,Nchannels(1)
*        **** receive packed array data ****
         if (int_mb(channel_type(1,1)+c-1) .eq. 1) then
            pfrom=int_mb(channel_proc(1,1)+c-1)
            it = mod((taskid+np-pfrom),np)

            source=pfrom
            msglen = (int_mb(h_i2_start(1,op,nbb)+it+2-1)
     >             -  int_mb(h_i2_start(1,op,nbb)+it+1-1))

            if (msglen.gt.0) then
               call RCV(9+MSGDBL,
     >                  tmp2(int_mb(h_i2_start(1,op,nbb)+it+1-1)),
     >                  mdtob(2*msglen),rcv_len,
     >                  Parallel2d_convert_taskid_i(source),
     >                  rcv_proc,1)
            end if
         end if

*        **** send packed array to other processors ****
         if (int_mb(channel_type(1,1)+c-1) .eq. 0) then
            pto=int_mb(channel_proc(1,1)+c-1)
            it = mod((pto-taskid+np),np)

            msglen    = (int_mb(h_i1_start(1,op,nbb)+it+2-1)
     >                -  int_mb(h_i1_start(1,op,nbb)+it+1-1))

            if (msglen.gt.0) then
               call SND(9+MSGDBL,
     >                  tmp1(int_mb(h_i1_start(1,op,nbb)+it+1-1)),
     >                  mdtob(2*msglen),
     >                  Parallel2d_convert_taskid_i(pto),
     >                  1)
            end if
         end if

      end do

*     **** unpack A(i) array ****
      call D3dB_pfft_index2_copy(n2,int_mb(h_iq_to_i2(1,op,nbb)),tmp2,A)
      call D3dB_pfft_index2_zero(n3,int_mb(h_iz_to_i2(1,op,nbb)),A)

      
      return
      end


*     ***********************************
*     *                                 *
*     *      D3dB_channel_init          *
*     *                                 *
*     ***********************************

      subroutine D3dB_channel_init(nb)
      implicit none
      integer nb

#include "bafdecls.fh"
#include "errquit.fh"
#include "D3dB.fh"


*     **** Used to avoid asynchronous communications ****
      integer Nchannels(NBLOCKS)
      integer channel_proc(2,NBLOCKS)
      integer channel_type(2,NBLOCKS)
      common / channel_blk / channel_proc,channel_type,Nchannels
      integer pair1(2),pair2(2)
      integer pair_step(2),pair_tmp(2)
      integer step,Nstep,icount,jcount,jcount_max

*     **** local variables ****
      integer np,taskid
      integer i,j,k
      logical value


      call Parallel2d_taskid_i(taskid)
      call Parallel2d_np_i(np)
*
*     **** Define Channels - which are used to avoid ****
*     **** asynchronous communications               ****

      value = BA_alloc_get(mt_int,(2*np),
     >        'channel_proc',channel_proc(2,nb),channel_proc(1,nb))
      value = value.and.
     >        BA_alloc_get(mt_int,(2*np),
     >        'channel_type',channel_type(2,nb),channel_type(1,nb))
      if (.not. value) 
     > call errquit('D3dB_channel_init:out of heap memory',0, MA_ERR)



      value = BA_push_get(mt_int,(np*(np-1)/2),
     >                    'pair1',pair1(2),pair1(1))
      value = value.and.
     >        BA_push_get(mt_int,(np*(np-1)/2),
     >                    'pair2',pair2(2),pair2(1))
      value = value.and.
     >        BA_push_get(mt_int,(np*(np-1)/2),
     >                'pair_step',pair_step(2),pair_step(1))
      value = value.and.
     >        BA_push_get(mt_int,(np*(np-1)/2),
     >                'pair_tmp',pair_tmp(2),pair_tmp(1))
      if (.not. value) 
     >  call errquit('D3dB_channel_init:out of stack memory',0, MA_ERR)

*     *** define pair1,pair2 ****
      icount = 0
      do i=0,     (np-1)
      do j=(i+1), (np-1)
         icount = icount + 1
         int_mb(pair1(1)+icount-1) = i
         int_mb(pair2(1)+icount-1) = j
      end do
      end do

*     **** define pair_step ****
      do i=1,(np*(np-1)/2)
         int_mb(pair_step(1)+i-1) = (-1)
      end do

      step=0
      jcount = 0
      jcount_max = np*(np-1)/2
c      do while(.not. full_ps(int_mb(pair_step(1)),np))
      do while(jcount.lt.jcount_max)
         step=step+1

         icount = 0
         do i=1, (np*(np-1)/2)
            if (int_mb(pair_step(1)+i-1).eq.(-1)) then
               value=.true.
               do k=1,icount
                 j = int_mb(pair_tmp(1)+k-1)

                 if
     >        ((int_mb(pair1(1)+i-1).eq.int_mb(pair1(1)+j-1)).or.
     >         (int_mb(pair1(1)+i-1).eq.int_mb(pair2(1)+j-1)).or.
     >         (int_mb(pair2(1)+i-1).eq.int_mb(pair1(1)+j-1)).or.
     >         (int_mb(pair2(1)+i-1).eq.int_mb(pair2(1)+j-1)))
     >           value=.false.

               end do
               if (value) then
                 int_mb(pair_step(1)+i-1) = step
                 icount = icount + 1
                 jcount = jcount + 1
                 int_mb(pair_tmp(1)+icount-1) = i
               end if
            end if
         end do
      end do
      Nstep=step


*     **** define channels ***
      Nchannels(nb)=0
      do step=1,Nstep
         do i=1,(np*(np-1)/2)
            if (int_mb(pair_step(1)+i-1).eq.step) then
*              **** send then recv ****
               if (int_mb(pair1(1)+i-1).eq.taskid) then
                  Nchannels(nb)=Nchannels(nb)+1
                  int_mb(channel_proc(1,nb)+Nchannels(nb)-1)
     >             = int_mb(pair2(1)+i-1)
                  int_mb(channel_type(1,nb)+Nchannels(nb)-1) = 0
                  Nchannels(nb)=Nchannels(nb)+1
                  int_mb(channel_proc(1,nb)+Nchannels(nb)-1)
     >            = int_mb(pair2(1)+i-1)
                  int_mb(channel_type(1,nb)+Nchannels(nb)-1) = 1
               end if

*              **** recv then send ****
               if (int_mb(pair2(1)+i-1).eq.taskid) then
                  Nchannels(nb)=Nchannels(nb)+1
                  int_mb(channel_proc(1,nb)+Nchannels(nb)-1)
     >            = int_mb(pair1(1)+i-1)
                  int_mb(channel_type(1,nb)+Nchannels(nb)-1) = 1
                  Nchannels(nb)=Nchannels(nb)+1
                  int_mb(channel_proc(1,nb)+Nchannels(nb)-1)
     >            = int_mb(pair1(1)+i-1)
                  int_mb(channel_type(1,nb)+Nchannels(nb)-1) = 0
               end if
            end if
         end do
      end do


      value=          BA_pop_stack(pair_tmp(2))
      value=value.and.BA_pop_stack(pair_step(2))
      value=value.and.BA_pop_stack(pair2(2))
      value=value.and.BA_pop_stack(pair1(2))
      if (.not. value) 
     >  call errquit('D3dB_channel_init:popping stack',0, MA_ERR)

      return
      end

c      logical function full_ps(ps,np)
c      implicit none
c      integer ps(*)
c      integer np
c
c      integer i
c      logical value
c
c      value=.true.
c      do i=1,(np*(np-1)/2)
c         if (ps(i).eq.(-1)) value=.false.
c      end do
c
c      full_ps=value
c      return
c      end

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

#include "D3dB.fh"

#include "tcgmsg.fh"
#include "msgtypesf.h"
      integer rcv_len,rcv_proc

      
*     *** local variables ***
      logical value,fillcolumn,readcolumn
      integer MASTER,taskid,taskid_i
      parameter(MASTER=0)
      integer p_to, p_here,q
      integer index,k,j
      integer source,msglen
      integer tmp1(2),tmp2(2)

      integer taskid_j,np_j
      integer ii,jj,jstart,jend

*     **** external functions ****
      integer  Parallel2d_convert_taskid_ij
      external Parallel2d_convert_taskid_ij

      call Parallel_taskid(taskid)
      call Parallel2d_taskid_i(taskid_i)
      call Parallel2d_taskid_j(taskid_j)

c      call Parallel2d_np_j(np_j)
c      call Parallel2d_taskid_j(taskid_j)
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
               p_to = Parallel2d_convert_taskid_ij(ii,taskid_j)
               !p_to = Parallel2d_convert_taskid_ij(ii,jj)
               if (ii.eq.MASTER) then
                  index = (q-1)*(nx(nb)/2+1)*ny(nb) + 1
                  call dcopy((nx(nb)+2)*ny(nb),tmp,1,A(index),1)
               else
                  msglen = (nx(nb)/2+1)*ny(nb)
                  call SND(9+MSGDBL,tmp,mdtob(2*msglen),p_to,1)
               end if
            !end do
         end do

*     **** not master node ****
      else if (fillcolumn) then
         do k=1,nz(nb)
            call D3dB_ktoqp(nb,k,q,ii)
            p_here = Parallel2d_convert_taskid_ij(ii,taskid_j)
            if (p_here.eq.taskid) then

               msglen = (nx(nb)/2+1)*ny(nb)
               !source = MASTER
               source = Parallel2d_convert_taskid_ij(MASTER,taskid_j)

               call RCV(9+MSGDBL,tmp,mdtob(2*msglen),rcv_len,
     >                  source,rcv_proc,1)


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
               p_to = Parallel2d_convert_taskid_ij(ii,taskid_j)
               !p_to = Parallel2d_convert_taskid_ij(ii,jj)

               if (ii.eq.MASTER) then
                  index = (q-1)*(nx(nb)/2+1) + 1
                  call dcopy((nx(nb)+2),tmp,1,A(index),1)
               else
                  msglen = (nx(nb)/2+1)
                  call SND(9+MSGDBL,tmp,mdtob(2*msglen),p_to,1)
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
            p_here = Parallel2d_convert_taskid_ij(ii,taskid_j)
            if (p_here.eq.taskid) then
               msglen = (nx(nb)/2+1)
               !source = MASTER
               source = Parallel2d_convert_taskid_ij(MASTER,taskid_j)
               call RCV(9+MSGDBL,tmp,mdtob(2*msglen),rcv_len,
     >                  source,rcv_proc,1)

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
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)
*
      call D3dB_c_transpose_ijk(nb,5,A,
     >                          dcpl_mb(tmp1(1)),
     >                          dcpl_mb(tmp2(1)))  !*** map1to3 operation ***

*     **** deallocate temporary space  ****
      value = BA_pop_stack(tmp2(2))
      value = value.and.BA_pop_stack(tmp1(2))
      if (.not. value) call errquit('error popping stack',0, MA_ERR)

      end if


      !*** shift filepointer by (nx(nb)+2)*ny(nb)*nz(nb) doubles ****
      else
         if (taskid_i.eq.MASTER)
     >      call dshift_fileptr(iunit,(nx(nb)+2)*ny(nb)*nz(nb))
      end if !*** readcolumn ***



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
#include "D3dB.fh"

#include "tcgmsg.fh"
#include "msgtypesf.h"
      integer rcv_len,rcv_proc

      
*     *** local variables ***
      logical fillcolumn,readcolumn
      integer MASTER,taskid,taskid_i
      parameter(MASTER=0)
      integer p_to, p_here,q
      integer j,k,index,index2
      integer source,msglen

      integer taskid_j,np_j
      integer ii,jj,jstart,jend

*     **** external functions ****
      integer  Parallel2d_convert_taskid_ij
      external Parallel2d_convert_taskid_ij

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
               p_to = Parallel2d_convert_taskid_ij(ii,taskid_j)
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
                  call SND(9+MSGDBL,tmp,mdtob(msglen),p_to,1)
               end if
            !end do
         end do

*     **** not master node ****
      else if (fillcolumn) then
         do k=1,nz(nb)

            call D3dB_ktoqp(nb,k,q,ii)
            p_here = Parallel2d_convert_taskid_ij(ii,taskid_j)
            if (p_here.eq.taskid) then

               msglen = (nx(nb))*ny(nb)
               !source  = MASTER
               source  = Parallel2d_convert_taskid_ij(MASTER,taskid_j)

               call RCV(9+MSGDBL,tmp,mdtob(msglen),rcv_len,
     >                       source,rcv_proc,1)
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

            q    = int_mb(q_map1(1,nb)+(j-1)+(k-1)*ny(nb))
            ii   = int_mb(p_map1(1,nb)+(j-1)+(k-1)*ny(nb))
            !do jj=jstart,jend
               p_to = Parallel2d_convert_taskid_ij(ii,taskid_j)
               !p_to = Parallel2d_convert_taskid_ij(ii,jj)
               if (p_to.eq.MASTER) then
                  index = (q-1)*(nx(nb)+2) + 1
                  call dcopy((nx(nb)+2),tmp,1,A(index),1)
               else
                  msglen = (nx(nb)+2)
                  call SND(9+MSGDBL,tmp,mdtob(msglen),p_to,1)
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
            p_here = Parallel2d_convert_taskid_ij(ii,taskid_j)
            if (p_here.eq.taskid) then

               msglen = (nx(nb)+2)
               !source  = MASTER
               source  = Parallel2d_convert_taskid_ij(MASTER,taskid_j)

               call RCV(9+MSGDBL,tmp,mdtob(msglen),rcv_len,
     >                  source,rcv_proc,1)

               index = (q-1)*(nx(nb)+2) + 1
               call dcopy((nx(nb)+2),tmp,1,A(index),1)

            end if
         end do
         end do
       end if

      end if

      !*** shift filepointer by nx(nb)*ny(nb)*nz(nb) doubles ****
      else
         if (taskid_i.eq.MASTER)
     >      call dshift_fileptr(iunit,nx(nb)*ny(nb)*nz(nb))
      end if !*** readcolumn ***

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
#include "D3dB.fh"



#include "tcgmsg.fh"
#include "msgtypesf.h"
      integer rcv_len,rcv_proc

      
*     *** local variables ***
      integer MASTER,taskid,taskid_i
      parameter(MASTER=0)

      logical value,fillcolumn,readcolumn
      integer p_to, p_here,q
      integer j,k,index
      integer source,msglen
      integer tmp1(2),tmp2(2)

      integer taskid_j,np_j
      integer ii,jj,jstart,jend

*     **** external functions ****
      integer  Parallel2d_convert_taskid_ij
      external Parallel2d_convert_taskid_ij

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


      if (mapping.eq.1) then
*        **** master node reads from file and distributes ****
         if (taskid_i.eq.MASTER) then
         do k=1,nz(nb)

            call dread(iunit,tmp,(nx(nb)/2+1)*ny(nb))

            call D3dB_ktoqp(nb,k,q,ii)
            !do jj=jstart,jend

              p_to = Parallel2d_convert_taskid_ij(ii,taskid_j)
              !p_to = Parallel2d_convert_taskid_ij(ii,jj)
              if (p_to.eq.MASTER) then
                 index = (q-1)*(nx(nb)/2+1)*ny(nb) + 1
                 call dcopy((nx(nb)/2+1)*ny(nb),tmp,1,A(index),1)
              else
                 msglen = (nx(nb)/2+1)*ny(nb)
                 call SND(9+MSGDBL,tmp,mdtob(msglen),p_to,1)
              end if
            !end do
         end do

*        **** not master node ****
         else if (fillcolumn) then

         do k=1,nz(nb)
            call D3dB_ktoqp(nb,k,q,ii)
            p_here = Parallel2d_convert_taskid_ij(ii,taskid_j)
            if (p_here.eq.taskid) then

               msglen = (nx(nb)/2+1)*ny(nb)
               !source  = MASTER
               source = Parallel2d_convert_taskid_ij(MASTER,taskid_j)

               call RCV(9+MSGDBL,tmp,mdtob(msglen),rcv_len,
     >                       source,rcv_proc,1)


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

               p_to = Parallel2d_convert_taskid_ij(ii,taskid_j)
               !p_to = Parallel2d_convert_taskid_ij(ii,jj)
               if (ii.eq.MASTER) then
                  index = (q-1)*(nx(nb)/2+1) + 1
                  call dcopy((nx(nb)/2+1),tmp,1,A(index),1)
               else
                  msglen = (nx(nb)/2+1)
                  call SND(9+MSGDBL,tmp,mdtob(msglen),p_to,1)
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
            p_here = Parallel2d_convert_taskid_ij(ii,taskid_j)
            if (p_here.eq.taskid) then

               msglen = (nx(nb)/2+1)
               !source  = MASTER
               source = Parallel2d_convert_taskid_ij(MASTER,taskid_j)

               call RCV(9+MSGDBL,tmp,mdtob(msglen),rcv_len,
     >                  source,rcv_proc,1)


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
      if (.not. value) call errquit('error popping stack',0, MA_ERR)



      end if


      !*** shift filepointer by (nx(nb)/2+1)*ny(nb)*nz(nb) doubles ****
      else
         if (taskid_i.eq.MASTER)
     >      call dshift_fileptr(iunit,(nx(nb)/2+1)*ny(nb)*nz(nb))
      end if !*** readcolumn ***

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
#include "D3dB.fh"


#include "tcgmsg.fh"
#include "msgtypesf.h"
      integer rcv_len,rcv_proc

      
*     *** local variables ***
      logical writecolumn
      integer MASTER,taskid,taskid_i
      parameter(MASTER=0)
      logical value
      integer p_from, p_here,q
      integer j,k,index
      integer dest,source,status,msglen
      integer dum,dum_msglen
      integer tmp1(2),tmp2(2)

      integer ii,taskid_j

*     **** external functions ****
      integer  Parallel2d_convert_taskid_ij
      external Parallel2d_convert_taskid_ij

      call Parallel_taskid(taskid)
      call Parallel2d_taskid_i(taskid_i)
      call Parallel2d_taskid_j(taskid_j)
      writecolumn = (taskid_j.eq.jcol)

      if (writecolumn) then

      !**********************
      !**** slab mapping ****
      !**********************
      if (mapping.eq.1) then
*     **** master node reads from file and distributes ****
      if (taskid_i.eq.MASTER) then
         do k=1,nz(nb)

            call D3dB_ktoqp(nb,k,q,ii)

            p_from = Parallel2d_convert_taskid_ij(ii,taskid_j)
            !p_from = Parallel2d_convert_taskid_ij(ii,jcol)
            if (ii.eq.MASTER) then

               index = (q-1)*(nx(nb)/2+1)*ny(nb) + 1
               call dcopy((nx(nb)+2)*ny(nb),A(index),1,tmp,1)

            else

               msglen  = (nx(nb)/2+1)*ny(nb)
               status  = msglen
               source  = p_from

               dum = 99
               dum_msglen = 1
               call SND(9+MSGINT,dum,mitob(dum_msglen),source,1)
               call RCV(9+MSGDBL,tmp,mdtob(2*msglen),rcv_len,
     >                       source,rcv_proc,1)

            end if

            call dwrite(iunit,tmp,(nx(nb)+2)*ny(nb)) 
         end do

*     **** not master node ****
      else
         do k=1,nz(nb)
         
            call D3dB_ktoqp(nb,k,q,ii)
            p_here = Parallel2d_convert_taskid_ij(ii,taskid_j)
            if (p_here.eq.taskid) then

               index = (q-1)*(nx(nb)/2+1)*ny(nb) + 1
               call dcopy((nx(nb)+2)*ny(nb),A(index),1,tmp,1)


               msglen = (nx(nb)/2+1)*ny(nb)
               !dest   = MASTER
               dest   = Parallel2d_convert_taskid_ij(MASTER,taskid_j)

               dum_msglen = 1
               call RCV(9+MSGINT,dum,mitob(dum_msglen),rcv_len,
     >                       dest,rcv_proc,1)
               call SND(9+MSGDBL,tmp,mdtob(2*msglen),dest,1)

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
        if (.not. value) call errquit('error popping stack',0, MA_ERR)

      end if


*     **** master node reads from file and distributes ****
      if (taskid_i.eq.MASTER) then
         do k=1,nz(nb)
         do j=1,ny(nb)

            q      = int_mb(q_map1(1,nb)+(j-1)+(k-1)*ny(nb))
            ii     = int_mb(p_map1(1,nb)+(j-1)+(k-1)*ny(nb))
            p_from = Parallel2d_convert_taskid_ij(ii,taskid_j)
            !p_from = Parallel2d_convert_taskid_ij(ii,jcol)
            if (ii.eq.MASTER) then

               index = (q-1)*(nx(nb)/2+1) + 1
               call dcopy((nx(nb)+2),A(index),1,tmp,1)

            else

               msglen  = (nx(nb)/2+1)
               status  = msglen
               source  = p_from

               dum = 99
               dum_msglen = 1
               call SND(9+MSGINT,dum,mitob(dum_msglen),source,1)
               call RCV(9+MSGDBL,tmp,mdtob(2*msglen),rcv_len,
     >                       source,rcv_proc,1)

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
            p_here = Parallel2d_convert_taskid_ij(ii,taskid_j) 
            if (p_here.eq.taskid) then

               index = (q-1)*(nx(nb)/2+1) + 1
               call dcopy((nx(nb)+2),A(index),1,tmp,1)


               msglen = (nx(nb)/2+1)
               !dest   = MASTER
               dest   = Parallel2d_convert_taskid_ij(MASTER,taskid_j) 

               dum_msglen = 1
               call RCV(9+MSGINT,dum,mitob(dum_msglen),rcv_len,
     >                       dest,rcv_proc,1)
               call SND(9+MSGDBL,tmp,mdtob(2*msglen),dest,1)

            end if

         end do
         end do
      end if

      end if


      !*** shift filepointer by (nx(nb)+2)*ny(nb)*nz(nb) doubles ****
      else
         if (taskid_i.eq.MASTER)
     >      call dshift_fileptr(iunit,(nx(nb)+2)*ny(nb)*nz(nb))
      end if !*** writecolumn ***


      return
      end

      subroutine D3dB_r_write_pio(nb,iunit,A,tmp,jcol)
      implicit none 
      integer nb
      integer iunit
      real*8     A(*)
      real*8     tmp(*)
      integer    jcol

#include "bafdecls.fh"
#include "D3dB.fh"



#include "tcgmsg.fh"
#include "msgtypesf.h"
      integer rcv_len,rcv_proc

      
*     *** local variables ***
      logical writecolumn
      integer MASTER,taskid,taskid_i
      parameter(MASTER=0)
      integer p_from, p_here,q
      integer j,k,index,index2
      integer dest,source,status,msglen

      integer taskid_j,ii

*     **** external functions ****
      integer  Parallel2d_convert_taskid_ij
      external Parallel2d_convert_taskid_ij

      call Parallel_taskid(taskid)
      call Parallel2d_taskid_i(taskid_i)
      call Parallel2d_taskid_j(taskid_j)
      writecolumn = (taskid_j.eq.jcol)

      if (writecolumn) then

      !**********************
      !**** slab mapping ****
      !**********************
      if (mapping.eq.1) then
*     **** master node reads from file and distributes ****
      if (taskid_i.eq.MASTER) then
         do k=1,nz(nb)

            call D3dB_ktoqp(nb,k,q,ii)
            p_from = Parallel2d_convert_taskid_ij(ii,taskid_j)
            if (ii.eq.MASTER) then

               do j=1,ny(nb)
                 index = (q-1)*(nx(nb)+2)*ny(nb) 
     >                 + (j-1)*(nx(nb)+2) + 1
                 index2 = (j-1)*nx(nb) + 1
                 call dcopy(nx(nb),A(index),1,tmp(index2),1)
               end do

            else

               msglen  = (nx(nb))*ny(nb)
               status  = msglen
               source  = p_from

               call RCV(9+MSGDBL,tmp,mdtob(msglen),rcv_len,
     >                       source,rcv_proc,1)

            end if

            call dwrite(iunit,tmp,(nx(nb))*ny(nb)) 
         end do

*     **** not master node ****
      else
         do k=1,nz(nb)
         
            call D3dB_ktoqp(nb,k,q,ii)
            p_here = Parallel2d_convert_taskid_ij(ii,taskid_j)
            !p_here = Parallel2d_convert_taskid_ij(ii,jcol)
            if (p_here.eq.taskid) then

               do j=1,ny(nb)
                  index = (q-1)*(nx(nb)+2)*ny(nb) 
     >                  + (j-1)*(nx(nb)+2) + 1
                  index2 = (j-1)*nx(nb) + 1
                  call dcopy(nx(nb),A(index),1,tmp(index2),1)
               end do


               msglen = (nx(nb))*ny(nb)
               !dest    = MASTER
               dest   = Parallel2d_convert_taskid_ij(MASTER,taskid_j)

               call SND(9+MSGDBL,tmp,mdtob(msglen),dest,1)

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

            q      = int_mb(q_map1(1,nb)+(j-1)+(k-1)*ny(nb))
            ii     = int_mb(p_map1(1,nb)+(j-1)+(k-1)*ny(nb))
            p_from = Parallel2d_convert_taskid_ij(ii,taskid_j)      

            if (ii.eq.MASTER) then
              index = (q-1)*(nx(nb)+2) + 1
              call dcopy(nx(nb),A(index),1,tmp,1)

            else

               msglen  = (nx(nb))
               status  = msglen
               source  = p_from

               call RCV(9+MSGDBL,tmp,mdtob(msglen),rcv_len,
     >                       source,rcv_proc,1)

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
            p_here = Parallel2d_convert_taskid_ij(ii,taskid_j)
            if (p_here.eq.taskid) then

               index = (q-1)*(nx(nb)+2) + 1
               call dcopy(nx(nb),A(index),1,tmp,1)

               msglen  = nx(nb)
               dest    = MASTER

               call SND(9+MSGDBL,tmp,mdtob(msglen),dest,1)

            end if

         end do
         end do
      end if

      end if

      !*** shift filepointer by nx(nb)*ny(nb)*nz(nb) doubles ****
      else
         if (taskid_i.eq.MASTER)
     >      call dshift_fileptr(iunit,nx(nb)*ny(nb)*nz(nb))
      end if !*** writecolumn ***

      return
      end

      subroutine D3dB_t_write_pio(nb,iunit,A,tmp,jcol)
      implicit none 
      integer nb
      integer iunit
      real*8  A(*)
      real*8  tmp(*)
      integer jcol

#include "bafdecls.fh"
#include "errquit.fh"
#include "D3dB.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"

      integer rcv_len,rcv_proc

      
*     *** local variables ***
      logical writecolumn
      integer MASTER,taskid,taskid_i
      parameter(MASTER=0)
      logical value
      integer p_from, p_here,q
      integer j,k,index
      integer dest,source,status,msglen
      integer dum,dum_msglen
      integer tmp1(2),tmp2(2)

      integer ii,taskid_j

*     **** external functions ****
      integer  Parallel2d_convert_taskid_ij
      external Parallel2d_convert_taskid_ij

      call Parallel_taskid(taskid)
      call Parallel2d_taskid_i(taskid_i)
      call Parallel2d_taskid_j(taskid_j)
      writecolumn = (taskid_j.eq.jcol)

      if (writecolumn) then

      !**********************
      !**** slab mapping ****
      !**********************
      if (mapping.eq.1) then
*     **** master node reads from file and distributes ****
      if (taskid_i.eq.MASTER) then
         do k=1,nz(nb)

            call D3dB_ktoqp(nb,k,q,ii)

            p_from = Parallel2d_convert_taskid_ij(ii,taskid_j)
            if (ii.eq.MASTER) then

               index = (q-1)*(nx(nb)/2+1)*ny(nb) + 1
               call dcopy((nx(nb)/2+1)*ny(nb),A(index),1,tmp,1)

            else

               msglen  = (nx(nb)/2+1)*ny(nb)
               status  = msglen
               source  = p_from

               dum = 99
               dum_msglen = 1
               call SND(9+MSGINT,dum,mitob(dum_msglen),source,1)
               call RCV(9+MSGDBL,tmp,mdtob(msglen),rcv_len,
     >                       source,rcv_proc,1)

            end if

            call dwrite(iunit,tmp,(nx(nb)/2+1)*ny(nb)) 
         end do

*     **** not master node ****
      else
         do k=1,nz(nb)
         
            call D3dB_ktoqp(nb,k,q,ii)
            p_here = Parallel2d_convert_taskid_ij(ii,taskid_j)
            if (p_here.eq.taskid) then

               index = (q-1)*(nx(nb)/2+1)*ny(nb) + 1
               call dcopy((nx(nb)/2+1)*ny(nb),A(index),1,tmp,1)


               msglen = (nx(nb)/2+1)*ny(nb)
               !dest   = MASTER
               dest   = Parallel2d_convert_taskid_ij(MASTER,taskid_j)

               dum_msglen = 1
               call RCV(9+MSGINT,dum,mitob(dum_msglen),rcv_len,
     >                       dest,rcv_proc,1)
               call SND(9+MSGDBL,tmp,mdtob(msglen),dest,1)

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
        if (.not. value) call errquit('error popping stack',0, MA_ERR)

      end if


*     **** master node reads from file and distributes ****
      if (taskid_i.eq.MASTER) then
         do k=1,nz(nb)
         do j=1,ny(nb)

            q      = int_mb(q_map1(1,nb)+(j-1)+(k-1)*ny(nb))
            ii     = int_mb(p_map1(1,nb)+(j-1)+(k-1)*ny(nb))
            p_from = Parallel2d_convert_taskid_ij(ii,taskid_j)
            if (ii.eq.MASTER) then

               index = (q-1)*(nx(nb)/2+1) + 1
               call dcopy((nx(nb)/2+1),A(index),1,tmp,1)

            else

               msglen  = (nx(nb)/2+1)
               status  = msglen
               source  = p_from

               dum = 99
               dum_msglen = 1
               call SND(9+MSGINT,dum,mitob(dum_msglen),source,1)
               call RCV(9+MSGDBL,tmp,mdtob(msglen),rcv_len,
     >                       source,rcv_proc,1)

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
            p_here = Parallel2d_convert_taskid_ij(ii,taskid_j) 
            if (p_here.eq.taskid) then

               index = (q-1)*(nx(nb)/2+1) + 1
               call dcopy((nx(nb)/2+1),A(index),1,tmp,1)


               msglen = (nx(nb)/2+1)
               !dest   = MASTER
               dest   = Parallel2d_convert_taskid_ij(MASTER,taskid_j) 

               dum_msglen = 1
               call RCV(9+MSGINT,dum,mitob(dum_msglen),rcv_len,
     >                       dest,rcv_proc,1)
               call SND(9+MSGDBL,tmp,mdtob(msglen),dest,1)

            end if

         end do
         end do
      end if
      end if

      !*** shift filepointer by (nx(nb)/2+1)*ny(nb)*nz(nb) doubles ****
      else
         if (taskid_i.eq.MASTER)
     >      call dshift_fileptr(iunit,(nx(nb)/2+1)*ny(nb)*nz(nb))
      end if !*** writecolumn ***


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

#include "D3dB.fh"

#include "tcgmsg.fh"
#include "msgtypesf.h"
      integer rcv_len,rcv_proc

      
*     *** local variables ***
      logical value,fillcolumn
      integer MASTER,taskid
      parameter(MASTER=0)
      integer p_to, p_here,q
      integer index,k,j
      integer source,msglen
      integer tmp1(2),tmp2(2)

      integer taskid_j,np_j
      integer ii,jj,jstart,jend

*     **** external functions ****
      integer  Parallel2d_convert_taskid_ij
      external Parallel2d_convert_taskid_ij

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
                  call SND(9+MSGDBL,tmp,mdtob(2*msglen),p_to,1)
               end if
            end do
         end do

*     **** not master node ****
      else if (fillcolumn) then
         do k=1,nz(nb)
            call D3dB_ktoqp(nb,k,q,ii)
            p_here = Parallel2d_convert_taskid_ij(ii,taskid_j)
            if (p_here.eq.taskid) then

               msglen = (nx(nb)/2+1)*ny(nb)
               source = MASTER

               call RCV(9+MSGDBL,tmp,mdtob(2*msglen),rcv_len,
     >                  source,rcv_proc,1)
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
                  msglen = (nx(nb)/2+1)
                  call SND(9+MSGDBL,tmp,mdtob(2*msglen),p_to,1)
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
               source = MASTER
               call RCV(9+MSGDBL,tmp,mdtob(2*msglen),rcv_len,
     >                  source,rcv_proc,1)

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
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)
*
      call D3dB_c_transpose_ijk(nb,5,A,
     >                          dcpl_mb(tmp1(1)),
     >                          dcpl_mb(tmp2(1)))  !*** map1to3 operation ***

*     **** deallocate temporary space  ****
      value = BA_pop_stack(tmp2(2))
      value = value.and.BA_pop_stack(tmp1(2))
      if (.not. value) call errquit('error popping stack',0, MA_ERR)

      end if


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
#include "D3dB.fh"

#include "tcgmsg.fh"
#include "msgtypesf.h"
      integer rcv_len,rcv_proc

      
*     *** local variables ***
      logical fillcolumn
      integer MASTER,taskid
      parameter(MASTER=0)
      integer p_to, p_here,q
      integer j,k,index,index2
      integer source,msglen

      integer taskid_j,np_j
      integer ii,jj,jstart,jend

*     **** external functions ****
      integer  Parallel2d_convert_taskid_ij
      external Parallel2d_convert_taskid_ij

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
                  call SND(9+MSGDBL,tmp,mdtob(msglen),p_to,1)
               end if
            end do
         end do

*     **** not master node ****
      else if (fillcolumn) then
         do k=1,nz(nb)

            call D3dB_ktoqp(nb,k,q,ii)
            p_here = Parallel2d_convert_taskid_ij(ii,taskid_j)
            if (p_here.eq.taskid) then

               msglen = (nx(nb))*ny(nb)
               source  = MASTER

               call RCV(9+MSGDBL,tmp,mdtob(msglen),rcv_len,
     >                       source,rcv_proc,1)
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

            q    = int_mb(q_map1(1,nb)+(j-1)+(k-1)*ny(nb))
            ii   = int_mb(p_map1(1,nb)+(j-1)+(k-1)*ny(nb))
            do jj=jstart,jend
               p_to = Parallel2d_convert_taskid_ij(ii,jj)
               if (p_to.eq.MASTER) then
                  index = (q-1)*(nx(nb)+2) + 1
                  call dcopy((nx(nb)+2),tmp,1,A(index),1)
               else
                  msglen = (nx(nb)+2)
                  call SND(9+MSGDBL,tmp,mdtob(msglen),p_to,1)
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

               msglen = (nx(nb)+2)
               source  = MASTER

               call RCV(9+MSGDBL,tmp,mdtob(msglen),rcv_len,
     >                  source,rcv_proc,1)

               index = (q-1)*(nx(nb)+2) + 1
               call dcopy((nx(nb)+2),tmp,1,A(index),1)

            end if
         end do
         end do
       end if

      end if

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
#include "D3dB.fh"



#include "tcgmsg.fh"
#include "msgtypesf.h"
      integer rcv_len,rcv_proc

      
*     *** local variables ***
      integer MASTER,taskid
      parameter(MASTER=0)

      logical value,fillcolumn
      integer p_to, p_here,q
      integer j,k,index
      integer source,msglen
      integer tmp1(2),tmp2(2)

      integer taskid_j,np_j
      integer ii,jj,jstart,jend

*     **** external functions ****
      integer  Parallel2d_convert_taskid_ij
      external Parallel2d_convert_taskid_ij

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



      if (mapping.eq.1) then
*        **** master node reads from file and distributes ****
         if (taskid.eq.MASTER) then
         do k=1,nz(nb)

            call dread(iunit,tmp,(nx(nb)/2+1)*ny(nb))

            call D3dB_ktoqp(nb,k,q,ii)
            !do jj=jstart,jend

              p_to = Parallel2d_convert_taskid_ij(ii,jj)
              if (p_to.eq.MASTER) then
                 index = (q-1)*(nx(nb)/2+1)*ny(nb) + 1
                 call dcopy((nx(nb)/2+1)*ny(nb),tmp,1,A(index),1)
              else
                 msglen = (nx(nb)/2+1)*ny(nb)
                 call SND(9+MSGDBL,tmp,mdtob(msglen),p_to,1)
              end if
            !end do
         end do

*        **** not master node ****
         else if (fillcolumn) then

         do k=1,nz(nb)
            call D3dB_ktoqp(nb,k,q,ii)
            p_here = Parallel2d_convert_taskid_ij(ii,taskid_j)
            if (p_here.eq.taskid) then

               msglen = (nx(nb)/2+1)*ny(nb)
               source  = MASTER
               call RCV(9+MSGDBL,tmp,mdtob(msglen),rcv_len,
     >                       source,rcv_proc,1)


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
                  msglen = (nx(nb)/2+1)
                  call SND(9+MSGDBL,tmp,mdtob(msglen),p_to,1)
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
               source  = MASTER

               call RCV(9+MSGDBL,tmp,mdtob(msglen),rcv_len,
     >                  source,rcv_proc,1)

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
      if (.not. value) call errquit('error popping stack',0, MA_ERR)



      end if


      return
      end




*     ***********************************
*     *					*
*     *	       D3dB_(c,r,t)_write	*	
*     *					*
*     ***********************************

      subroutine D3dB_c_write(nb,iunit,A,tmp,jcol)
      implicit none 
      integer nb
      integer iunit
      complex*16 A(*)
      complex*16 tmp(*)
      integer    jcol


#include "bafdecls.fh"
#include "errquit.fh"
#include "D3dB.fh"


#include "tcgmsg.fh"
#include "msgtypesf.h"
      integer rcv_len,rcv_proc

      
*     *** local variables ***
      integer MASTER,taskid
      parameter(MASTER=0)
      logical value
      integer p_from, p_here,q
      integer j,k,index
      integer dest,source,status,msglen
      integer dum,dum_msglen
      integer tmp1(2),tmp2(2)

      integer ii,taskid_j

*     **** external functions ****
      integer  Parallel2d_convert_taskid_ij
      external Parallel2d_convert_taskid_ij

      call Parallel_taskid(taskid)
      call Parallel2d_taskid_j(taskid_j)


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
               call dcopy((nx(nb)+2)*ny(nb),A(index),1,tmp,1)

            else

               msglen  = (nx(nb)/2+1)*ny(nb)
               status  = msglen
               source  = p_from

               dum = 99
               dum_msglen = 1
               call SND(9+MSGINT,dum,mitob(dum_msglen),source,1)
               call RCV(9+MSGDBL,tmp,mdtob(2*msglen),rcv_len,
     >                       source,rcv_proc,1)

            end if

            call dwrite(iunit,tmp,(nx(nb)+2)*ny(nb)) 
         end do

*     **** not master node ****
      else
         do k=1,nz(nb)
         
            call D3dB_ktoqp(nb,k,q,ii)
            p_here = Parallel2d_convert_taskid_ij(ii,taskid_j)
            if (p_here.eq.taskid) then

               index = (q-1)*(nx(nb)/2+1)*ny(nb) + 1
               call dcopy((nx(nb)+2)*ny(nb),A(index),1,tmp,1)


               msglen = (nx(nb)/2+1)*ny(nb)
               dest   = MASTER

               dum_msglen = 1
               call RCV(9+MSGINT,dum,mitob(dum_msglen),rcv_len,
     >                       dest,rcv_proc,1)
               call SND(9+MSGDBL,tmp,mdtob(2*msglen),dest,1)

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
        if (.not. value) call errquit('error popping stack',0, MA_ERR)

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
               call dcopy((nx(nb)+2),A(index),1,tmp,1)

            else

               msglen  = (nx(nb)/2+1)
               status  = msglen
               source  = p_from

               dum = 99
               dum_msglen = 1
               call SND(9+MSGINT,dum,mitob(dum_msglen),source,1)
               call RCV(9+MSGDBL,tmp,mdtob(2*msglen),rcv_len,
     >                       source,rcv_proc,1)

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
            p_here = Parallel2d_convert_taskid_ij(ii,taskid_j) 
            if (p_here.eq.taskid) then

               index = (q-1)*(nx(nb)/2+1) + 1
               call dcopy((nx(nb)+2),A(index),1,tmp,1)

               msglen = (nx(nb)/2+1)
               dest   = MASTER

               dum_msglen = 1
               call RCV(9+MSGINT,dum,mitob(dum_msglen),rcv_len,
     >                       dest,rcv_proc,1)
               call SND(9+MSGDBL,tmp,mdtob(2*msglen),dest,1)

            end if

         end do
         end do
      end if

      end if



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
#include "D3dB.fh"

#include "tcgmsg.fh"
#include "msgtypesf.h"
      integer rcv_len,rcv_proc

      
*     *** local variables ***
      integer MASTER,taskid
      parameter(MASTER=0)
      integer p_from, p_here,q
      integer j,k,index,index2
      integer dest,source,status,msglen

      integer taskid_j,ii

*     **** external functions ****
      integer  Parallel2d_convert_taskid_ij
      external Parallel2d_convert_taskid_ij

      call Parallel_taskid(taskid)
      call Parallel2d_taskid_j(taskid_j)


      !**********************
      !**** slab mapping ****
      !**********************
      if (mapping.eq.1) then
*     **** master node reads from file and distributes ****
      if (taskid.eq.MASTER) then
         do k=1,nz(nb)

            call D3dB_ktoqp(nb,k,q,ii)
            p_from = Parallel2d_convert_taskid_ij(ii,taskid_j)
            if (p_from.eq.MASTER) then

               do j=1,ny(nb)
                 index = (q-1)*(nx(nb)+2)*ny(nb) 
     >                 + (j-1)*(nx(nb)+2) + 1
                 index2 = (j-1)*nx(nb) + 1
                 call dcopy(nx(nb),A(index),1,tmp(index2),1)
               end do

            else

               msglen  = (nx(nb))*ny(nb)
               status  = msglen
               source  = p_from

               call RCV(9+MSGDBL,tmp,mdtob(msglen),rcv_len,
     >                       source,rcv_proc,1)

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


               msglen = (nx(nb))*ny(nb)
               dest    = MASTER

               call SND(9+MSGDBL,tmp,mdtob(msglen),dest,1)

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
            p_from = Parallel2d_convert_taskid_ij(ii,taskid_j)      

            if (p_from.eq.MASTER) then
              index = (q-1)*(nx(nb)+2) + 1
              call dcopy(nx(nb),A(index),1,tmp,1)

            else

               msglen  = (nx(nb))
               status  = msglen
               source  = p_from

               call RCV(9+MSGDBL,tmp,mdtob(msglen),rcv_len,
     >                       source,rcv_proc,1)

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
            p_here = Parallel2d_convert_taskid_ij(ii,taskid_j)
            if (p_here.eq.taskid) then

               index = (q-1)*(nx(nb)+2) + 1
               call dcopy(nx(nb),A(index),1,tmp,1)

               msglen  = nx(nb)
               dest    = MASTER

               call SND(9+MSGDBL,tmp,mdtob(msglen),dest,1)

            end if

         end do
         end do
      end if

      end if

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
#include "D3dB.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"

      integer rcv_len,rcv_proc

      
*     *** local variables ***
      integer MASTER,taskid
      parameter(MASTER=0)
      logical value
      integer p_from, p_here,q
      integer j,k,index
      integer dest,source,status,msglen
      integer dum,dum_msglen
      integer tmp1(2),tmp2(2)

      integer ii,taskid_j

*     **** external functions ****
      integer  Parallel2d_convert_taskid_ij
      external Parallel2d_convert_taskid_ij

      call Parallel_taskid(taskid)
      call Parallel2d_taskid_j(taskid_j)


      !**********************
      !**** slab mapping ****
      !**********************
      if (mapping.eq.1) then
*     **** master node reads from file and distributes ****
      if (taskid.eq.MASTER) then
         do k=1,nz(nb)

            call D3dB_ktoqp(nb,k,q,ii)

            p_from = Parallel2d_convert_taskid_ij(ii,taskid_j)
            if (p_from.eq.MASTER) then

               index = (q-1)*(nx(nb)/2+1)*ny(nb) + 1
               call dcopy((nx(nb)/2+1)*ny(nb),A(index),1,tmp,1)

            else

               msglen  = (nx(nb)/2+1)*ny(nb)
               status  = msglen
               source  = p_from

               dum = 99
               dum_msglen = 1
               call SND(9+MSGINT,dum,mitob(dum_msglen),source,1)
               call RCV(9+MSGDBL,tmp,mdtob(msglen),rcv_len,
     >                       source,rcv_proc,1)

            end if

            call dwrite(iunit,tmp,(nx(nb)/2+1)*ny(nb)) 
         end do

*     **** not master node ****
      else
         do k=1,nz(nb)
         
            call D3dB_ktoqp(nb,k,q,ii)
            p_here = Parallel2d_convert_taskid_ij(ii,taskid_j)
            if (p_here.eq.taskid) then

               index = (q-1)*(nx(nb)/2+1)*ny(nb) + 1
               call dcopy((nx(nb)/2+1)*ny(nb),A(index),1,tmp,1)


               msglen = (nx(nb)/2+1)*ny(nb)
               dest   = MASTER

               dum_msglen = 1
               call RCV(9+MSGINT,dum,mitob(dum_msglen),rcv_len,
     >                       dest,rcv_proc,1)
               call SND(9+MSGDBL,tmp,mdtob(msglen),dest,1)

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
        if (.not. value) call errquit('error popping stack',0, MA_ERR)

      end if


*     **** master node reads from file and distributes ****
      if (taskid.eq.MASTER) then
         do k=1,nz(nb)
         do j=1,ny(nb)

            q      = int_mb(q_map1(1,nb)+(j-1)+(k-1)*ny(nb))
            ii     = int_mb(p_map1(1,nb)+(j-1)+(k-1)*ny(nb))
            p_from = Parallel2d_convert_taskid_ij(ii,taskid_j)
            if (p_from.eq.MASTER) then

               index = (q-1)*(nx(nb)/2+1) + 1
               call dcopy((nx(nb)/2+1),A(index),1,tmp,1)

            else

               msglen  = (nx(nb)/2+1)
               status  = msglen
               source  = p_from

               dum = 99
               dum_msglen = 1
               call SND(9+MSGINT,dum,mitob(dum_msglen),source,1)
               call RCV(9+MSGDBL,tmp,mdtob(msglen),rcv_len,
     >                       source,rcv_proc,1)

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
            p_here = Parallel2d_convert_taskid_ij(ii,taskid_j) 
            if (p_here.eq.taskid) then

               index = (q-1)*(nx(nb)/2+1) + 1
               call dcopy((nx(nb)/2+1),A(index),1,tmp,1)


               msglen = (nx(nb)/2+1)
               dest   = MASTER

               dum_msglen = 1
               call RCV(9+MSGINT,dum,mitob(dum_msglen),rcv_len,
     >                       dest,rcv_proc,1)
               call SND(9+MSGDBL,tmp,mdtob(msglen),dest,1)

            end if

         end do
         end do
      end if
      end if


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
#include "D3dB.fh"


#include "tcgmsg.fh"
#include "msgtypesf.h"
      integer rcv_len,rcv_proc

      
*     *** local variables ***
      integer MASTER,taskid
      parameter(MASTER=0)
      integer p_from, p_here,q
      integer i,j,k,index,index2
      integer dest,source,status,msglen

      call Parallel_taskid(taskid)

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
               status  = msglen
               source  = p_from

               call RCV(9+MSGDBL,tmp,mdtob(msglen),rcv_len,
     >                       source,rcv_proc,1)

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
               dest    = MASTER

               call SND(9+MSGDBL,tmp,mdtob(msglen),dest,1)

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
               status  = msglen
               source  = p_from

               call RCV(9+MSGDBL,tmp,mdtob(msglen),rcv_len,
     >                       source,rcv_proc,1)

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
               dest    = MASTER

               call SND(9+MSGDBL,tmp,mdtob(msglen),dest,1)

            end if

         end do
         end do
      end if



      end if

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

#include "D3dB.fh"

#include "tcgmsg.fh"
#include "msgtypesf.h"
      integer rcv_len,rcv_proc

      
*     *** local variables ***
      integer MASTER,taskid
      parameter(MASTER=0)
      integer p_from, p_here,q
      integer i,j,k,index
      integer dest,source,status,msglen,idum

      call Parallel_taskid(taskid)

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
                 status  = msglen
                 source  = p_from

                 call SND(9+MSGINT,idum,mitob(msglen),source,1)
                 call RCV(9+MSGDBL,tmp(k),mdtob(msglen),rcv_len,
     >                         source,rcv_proc,1)

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
                 dest    = MASTER

                 call RCV(9+MSGINT,idum,mitob(msglen),rcv_len,
     >                         dest,rcv_proc,1)
                 call SND(9+MSGDBL,tmp,mdtob(msglen),dest,1)

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
                 status  = msglen
                 source  = p_from

                 call SND(9+MSGINT,idum,mitob(msglen),source,1)
                 call RCV(9+MSGDBL,tmp(k),mdtob(msglen),rcv_len,
     >                         source,rcv_proc,1)

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
                 dest    = MASTER

                 call RCV(9+MSGINT,idum,mitob(msglen),rcv_len,
     >                         dest,rcv_proc,1)
                 call SND(9+MSGDBL,tmp,mdtob(msglen),dest,1)

              end if
            end do

         end do
         end do
      end if

      end if

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

#include "D3dB.fh"


#include "tcgmsg.fh"
#include "msgtypesf.h"
      integer rcv_len,rcv_proc

      
*     *** local variables ***
      integer MASTER,taskid
      parameter(MASTER=0)
      integer p_from, p_here,q
      integer i,j,k,index
      integer dest,source,status,msglen

      call Parallel_taskid(taskid)

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
               status  = msglen
               source  = p_from

               call RCV(9+MSGDBL,tmp,mdtob(msglen),rcv_len,
     >                       source,rcv_proc,1)

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
               dest    = MASTER

               call SND(9+MSGDBL,tmp,mdtob(msglen),dest,1)

            end if

         end do
         end do
      end if

*     **** wait ****
      return
      end



