#define NBLOCKS 2

*
* $Id: D3dB-tcgmsg.F,v 1.6 2004-05-05 19:42:38 edo Exp $
*

#ifdef TCGMSG
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
      
      call Parallel_taskid(taskid)
      call Parallel_np(np)

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
c#include "mafdecls.fh"
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

      call Parallel_taskid(taskid)
      call Parallel_np(np)

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

#include "mafdecls.fh"

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
      
      call Parallel_taskid(taskid)
      call Parallel_np(np)


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

#include "mafdecls.fh"

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
      
      call Parallel_taskid(taskid)
      call Parallel_np(np)


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




#endif

