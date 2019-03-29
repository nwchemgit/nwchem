#define NBLOCKS 2



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

#include "bafdecls.fh"
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

*     **** external functions ****
      integer  Parallel3d_convert_taskid_i
      external Parallel3d_convert_taskid_i

      
      call Parallel3d_taskid_i(taskid)
      call Parallel3d_np_i(np)

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
     >                  Parallel3d_convert_taskid_i(source),
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
     >                  Parallel3d_convert_taskid_i(pto),1)
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

#include "bafdecls.fh"
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

*     **** external functions ****
      integer  Parallel3d_convert_taskid_i
      external Parallel3d_convert_taskid_i
      
      call Parallel3d_taskid_i(taskid)
      call Parallel3d_np_i(np)

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
     >                  Parallel3d_convert_taskid_i(source),
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
     >                  mdtob(msglen),
     >                  Parallel3d_convert_taskid_i(pto),1)
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

#include "bafdecls.fh"

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

*     **** external functions ****
      integer  Parallel3d_convert_taskid_i
      external Parallel3d_convert_taskid_i
      
      call Parallel3d_taskid_i(taskid)
      call Parallel3d_np_i(np)


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
     >                  Parallel3d_convert_taskid_i(source),
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
     >                  Parallel3d_convert_taskid_i(pto),1)
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

#include "bafdecls.fh"

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

*     **** external functions ****
      integer  Parallel3d_convert_taskid_i
      external Parallel3d_convert_taskid_i
      
      call Parallel3d_taskid_i(taskid)
      call Parallel3d_np_i(np)


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
     >                  Parallel3d_convert_taskid_i(source),
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
     >                  Parallel3d_convert_taskid_i(pto),1)
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
*     *	   C3dB_c_ptranspose1_jk	*
*     *					*
*     ***********************************

      subroutine C3dB_c_ptranspose1_jk(fb,A,tmp1,tmp2)

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
      integer     fb
      complex*16  A(*)
      complex*16  tmp1(*),tmp2(*)

#include "bafdecls.fh"
#include "C3dB.fh"

*     **** indexing variables ****
      integer iq_to_i1(2,0:1)
      integer iq_to_i2(2,0:1)
      integer iz_to_i2(2,0:1)
      integer i1_start(2,0:1)
      integer i2_start(2,0:1)
      common / c_ptrans_blk1 / iq_to_i1,iq_to_i2,iz_to_i2,
     >                         i1_start,i2_start

*     **** Used to avoid asynchronous communications ****
      integer Nchannels(NBLOCKS)
      integer channel_proc(2,NBLOCKS)
      integer channel_type(2,NBLOCKS)
      common / c_channel_blk / channel_proc,channel_type,Nchannels

#include "tcgmsg.fh"
#include "msgtypesf.h"
      integer  rcv_len,rcv_proc

*     **** local variables ***
      integer c,it
      integer source
      integer msglen
      integer pfrom,pto
      integer taskid,np
      integer n1,n2

*     **** external functions ****
      integer  Parallel3d_convert_taskid_i
      external Parallel3d_convert_taskid_i

      call Parallel3d_taskid_i(taskid)
      call Parallel3d_np_i(np)

      n1 = int_mb(i1_start(1,fb)+np) - 1
      n2 = int_mb(i2_start(1,fb)+np) - 1

*     **** pack A(i) array ****
      call C3dB_pfft_index1_copy(n1,int_mb(iq_to_i1(1,fb)),A,tmp1)


*     **** it = 0, transpose data on same thread ****  
      msglen = int_mb(i2_start(1,fb)+2-1) - int_mb(i2_start(1,fb)+1-1)
      call dcopy(2*msglen,
     >           tmp1(int_mb(i1_start(1,fb)+1-1)),1,
     >           tmp2(int_mb(i2_start(1,fb)+1-1)),1)


      do c=1,Nchannels(1)
*        **** receive packed array data ****
         if (int_mb(channel_type(1,1)+c-1) .eq. 1) then
            pfrom=int_mb(channel_proc(1,1)+c-1)
            it = mod((taskid+np-pfrom),np)

            source=pfrom
            msglen = (int_mb(i2_start(1,fb)+it+2-1) 
     >             -  int_mb(i2_start(1,fb)+it+1-1))
         
            if (msglen.gt.0) then
               call RCV(9+MSGDBL,
     >                  tmp2(int_mb(i2_start(1,fb)+it+1-1)),
     >                  mdtob(2*msglen),rcv_len,
     >                  Parallel3d_convert_taskid_i(source),
     >                  rcv_proc,1)
            end if
         end if

*        **** send packed array to other processors ****
         if (int_mb(channel_type(1,1)+c-1) .eq. 0) then
            pto=int_mb(channel_proc(1,1)+c-1)
            it = mod((pto-taskid+np),np)

            msglen    = (int_mb(i1_start(1,fb)+it+2-1)
     >                 - int_mb(i1_start(1,fb)+it+1-1))
      
            if (msglen.gt.0) then
               call SND(9+MSGDBL,
     >                  tmp1(int_mb(i1_start(1,fb)+it+1-1)),
     >                  mdtob(2*msglen),
     >                  Parallel3d_convert_taskid_i(pto),
     >                  1)
            end if
         end if

      end do

*     **** unpack A(i) array ****
      call C3dB_pfft_index2_copy(n2,int_mb(iq_to_i2(1,fb)),tmp2,A)
      call C3dB_pfft_index2_zero(nfft3d(1)-n2,int_mb(iz_to_i2(1,fb)),A)

      
      return
      end



*     ***********************************
*     *					*
*     *	   C3dB_c_ptranspose2_jk	*
*     *					*
*     ***********************************

      subroutine C3dB_c_ptranspose2_jk(fb,A,tmp1,tmp2)

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
      integer     fb
      complex*16  A(*)
      complex*16  tmp1(*),tmp2(*)

#include "bafdecls.fh"
#include "C3dB.fh"

*     **** indexing variables ****
      integer iq_to_i1(2,0:1)
      integer iq_to_i2(2,0:1)
      integer iz_to_i2(2,0:1)
      integer i1_start(2,0:1)
      integer i2_start(2,0:1)
      common / c_ptrans_blk2 / iq_to_i1,iq_to_i2,iz_to_i2,
     >                       i1_start,i2_start



*     **** Used to avoid asynchronous communications ****
      integer Nchannels(NBLOCKS)
      integer channel_proc(2,NBLOCKS)
      integer channel_type(2,NBLOCKS)
      common / c_channel_blk / channel_proc,channel_type,Nchannels

#include "tcgmsg.fh"
#include "msgtypesf.h"
      integer  rcv_len,rcv_proc

*     **** local variables ***
      integer c,it
      integer source
      integer msglen
      integer pfrom,pto
      integer taskid,np
      integer n1,n2

*     **** external functions ****
      integer  Parallel3d_convert_taskid_i
      external Parallel3d_convert_taskid_i

      call Parallel3d_taskid_i(taskid)
      call Parallel3d_np_i(np)

      n1 = int_mb(i1_start(1,fb)+np) - 1
      n2 = int_mb(i2_start(1,fb)+np) - 1

*     **** pack A(i) array ****
      call C3dB_pfft_index1_copy(n1,int_mb(iq_to_i1(1,fb)),A,tmp1)


*     **** it = 0, transpose data on same thread ****  
      msglen = int_mb(i2_start(1,fb)+2-1) - int_mb(i2_start(1,fb)+1-1)
      call dcopy(2*msglen,
     >           tmp1(int_mb(i1_start(1,fb)+1-1)),1,
     >           tmp2(int_mb(i2_start(1,fb)+1-1)),1)


      do c=1,Nchannels(1)
*        **** receive packed array data ****
         if (int_mb(channel_type(1,1)+c-1) .eq. 1) then
            pfrom=int_mb(channel_proc(1,1)+c-1)
            it = mod((taskid+np-pfrom),np)

            source=pfrom
            msglen = (int_mb(i2_start(1,fb)+it+2-1) 
     >             -  int_mb(i2_start(1,fb)+it+1-1))
         
            if (msglen.gt.0) then
               call RCV(9+MSGDBL,
     >                  tmp2(int_mb(i2_start(1,fb)+it+1-1)),
     >                  mdtob(2*msglen),rcv_len,
     >                  Parallel3d_convert_taskid_i(source),
     >                  rcv_proc,1)
            end if
         end if

*        **** send packed array to other processors ****
         if (int_mb(channel_type(1,1)+c-1) .eq. 0) then
            pto=int_mb(channel_proc(1,1)+c-1)
            it = mod((pto-taskid+np),np)

            msglen    = (int_mb(i1_start(1,fb)+it+2-1)
     >                 - int_mb(i1_start(1,fb)+it+1-1))
      
            if (msglen.gt.0) then
               call SND(9+MSGDBL,
     >                  tmp1(int_mb(i1_start(1,fb)+it+1-1)),
     >                  mdtob(2*msglen),
     >                  Parallel3d_convert_taskid_i(pto),
     >                  1)
            end if
         end if

      end do

*     **** unpack A(i) array ****
      call C3dB_pfft_index2_copy(n2,int_mb(iq_to_i2(1,fb)),tmp2,A)
      call C3dB_pfft_index2_zero(nfft3d(1)-n2,int_mb(iz_to_i2(1,fb)),A)
      
      return
      end




*     ***********************************
*     *					*
*     *	   C3dB_c_ptranspose1_jk_start	*
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

      subroutine C3dB_c_ptranspose1_jk_start(fb,A,tmp1,tmp2,
     >                                       request,reqcnt,msgtype)
      implicit none
      integer fb
      complex*16  A(*)
      complex*16  tmp1(*),tmp2(*)
      integer request(*),reqcnt,msgtype

#include "bafdecls.fh"
#include "errquit.fh"
#include "C3dB.fh"

*     **** indexing variables ****
      integer iq_to_i1(2,0:1)
      integer iq_to_i2(2,0:1)
      integer iz_to_i2(2,0:1)
      integer i1_start(2,0:1)
      integer i2_start(2,0:1)
      common / c_ptrans_blk1 / iq_to_i1,iq_to_i2,iz_to_i2,
     >                       i1_start,i2_start

*     **** Used to avoid asynchronous communications ****
      integer Nchannels(NBLOCKS)
      integer channel_proc(2,NBLOCKS)
      integer channel_type(2,NBLOCKS)
      common / c_channel_blk / channel_proc,channel_type,Nchannels

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
      integer  Parallel3d_convert_taskid_i
      external Parallel3d_convert_taskid_i

      
      call Parallel3d_taskid_i(taskid)
      call Parallel3d_np_i(np)

      n1 = int_mb(i1_start(1,fb)+np) - 1

*     **** pack A(i) array ****
       call C3dB_pfft_index1_copy(n1,int_mb(iq_to_i1(1,fb)),A,tmp1)

*     **** it = 0, transpose data on same thread ****  
      msglen = int_mb(i2_start(1,fb)+2-1) - int_mb(i2_start(1,fb)+1-1)
      call dcopy(2*msglen,
     >           tmp1(int_mb(i1_start(1,fb)+1-1)),1,
     >           tmp2(int_mb(i2_start(1,fb)+1-1)),1)
         

      do c=1,Nchannels(1)
*        **** receive packed array data ****
         if (int_mb(channel_type(1,1)+c-1) .eq. 1) then
            pfrom=int_mb(channel_proc(1,1)+c-1)
            it = mod((taskid+np-pfrom),np)

            source=pfrom
            msglen = (int_mb(i2_start(1,fb)+it+2-1)
     >             -  int_mb(i2_start(1,fb)+it+1-1))

            if (msglen.gt.0) then
               call RCV(9+MSGDBL,
     >                  tmp2(int_mb(i2_start(1,fb)+it+1-1)),
     >                  mdtob(2*msglen),rcv_len,
     >                  Parallel3d_convert_taskid_i(source),
     >                  rcv_proc,1)
            end if
         end if

*        **** send packed array to other processors ****
         if (int_mb(channel_type(1,1)+c-1) .eq. 0) then
            pto=int_mb(channel_proc(1,1)+c-1)
            it = mod((pto-taskid+np),np)

            msglen    = (int_mb(i1_start(1,fb)+it+2-1)
     >                 - int_mb(i1_start(1,fb)+it+1-1))

            if (msglen.gt.0) then
               call SND(9+MSGDBL,
     >                  tmp1(int_mb(i1_start(1,fb)+it+1-1)),
     >                  mdtob(2*msglen),
     >                  Parallel3d_convert_taskid_i(pto),1)
            end if
         end if

      end do

      return
      end


*     ***********************************
*     *					*
*     *	   C3dB_c_ptranspose1_jk_end	*
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

      subroutine C3dB_c_ptranspose1_jk_end(fb,A,tmp2,request,reqcnt)

      implicit none
      integer fb
      complex*16  A(*)
      complex*16  tmp2(*)
      integer     request(*),reqcnt

#include "bafdecls.fh"
#include "errquit.fh"
#include "C3dB.fh"


*     **** indexing variables ****
      integer iq_to_i1(2,0:1)
      integer iq_to_i2(2,0:1)
      integer iz_to_i2(2,0:1)
      integer i1_start(2,0:1)
      integer i2_start(2,0:1)
      common / c_ptrans_blk1 / iq_to_i1,iq_to_i2,iz_to_i2,
     >                       i1_start,i2_start

*     **** local variables ***
      integer np,n2
      
      call Parallel3d_np_i(np)

*     **** unpack A(i) array ****
      n2 = int_mb(i2_start(1,fb)+np) - 1
      call C3dB_pfft_index2_copy(n2,int_mb(iq_to_i2(1,fb)),tmp2,A)
      call C3dB_pfft_index2_zero(nfft3d(1)-n2,int_mb(iz_to_i2(1,fb)),A)

      return
      end


      
*     ***********************************
*     *					*
*     *	   C3dB_c_ptranspose2_jk_start	*
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

      subroutine C3dB_c_ptranspose2_jk_start(fb,A,tmp1,tmp2,
     >                                       request,reqcnt,msgtype)
      implicit none
      integer fb
      complex*16  A(*)
      complex*16  tmp1(*),tmp2(*)
      integer request(*),reqcnt,msgtype

#include "bafdecls.fh"
#include "errquit.fh"
#include "C3dB.fh"


*     **** indexing variables ****
      integer iq_to_i1(2,0:1)
      integer iq_to_i2(2,0:1)
      integer iz_to_i2(2,0:1)
      integer i1_start(2,0:1)
      integer i2_start(2,0:1)
      common / c_ptrans_blk2 / iq_to_i1,iq_to_i2,iz_to_i2,
     >                       i1_start,i2_start



*     **** Used to avoid asynchronous communications ****
      integer Nchannels(NBLOCKS)
      integer channel_proc(2,NBLOCKS)
      integer channel_type(2,NBLOCKS)
      common / c_channel_blk / channel_proc,channel_type,Nchannels

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
      integer  Parallel3d_convert_taskid_i
      external Parallel3d_convert_taskid_i

      
      call Parallel3d_taskid_i(taskid)
      call Parallel3d_np_i(np)


      n1 = int_mb(i1_start(1,fb)+np) - 1

*     **** pack A(i) array ****
       call C3dB_pfft_index1_copy(n1,int_mb(iq_to_i1(1,fb)),A,tmp1)

*     **** it = 0, transpose data on same thread ****  
      msglen = int_mb(i2_start(1,fb)+2-1) - int_mb(i2_start(1,fb)+1-1)
      call dcopy(2*msglen,
     >           tmp1(int_mb(i1_start(1,fb)+1-1)),1,
     >           tmp2(int_mb(i2_start(1,fb)+1-1)),1)
         

      do c=1,Nchannels(1)
*        **** receive packed array data ****
         if (int_mb(channel_type(1,1)+c-1) .eq. 1) then
            pfrom=int_mb(channel_proc(1,1)+c-1)
            it = mod((taskid+np-pfrom),np)

            source=pfrom
            msglen = (int_mb(i2_start(1,fb)+it+2-1)
     >             -  int_mb(i2_start(1,fb)+it+1-1))

            if (msglen.gt.0) then
               call RCV(9+MSGDBL,
     >                  tmp2(int_mb(i2_start(1,fb)+it+1-1)),
     >                  mdtob(2*msglen),rcv_len,
     >                  Parallel3d_convert_taskid_i(source),
     >                  rcv_proc,1)
            end if
         end if

*        **** send packed array to other processors ****
         if (int_mb(channel_type(1,1)+c-1) .eq. 0) then
            pto=int_mb(channel_proc(1,1)+c-1)
            it = mod((pto-taskid+np),np)

            msglen    = (int_mb(i1_start(1,fb)+it+2-1)
     >                 - int_mb(i1_start(1,fb)+it+1-1))

            if (msglen.gt.0) then
               call SND(9+MSGDBL,
     >                  tmp1(int_mb(i1_start(1,fb)+it+1-1)),
     >                  mdtob(2*msglen),
     >                  Parallel3d_convert_taskid_i(pto),
     >                  1)
            end if
         end if

      end do


      return
      end


*     ***********************************
*     *					*
*     *	   C3dB_c_ptranspose2_jk_end	*
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

      subroutine C3dB_c_ptranspose2_jk_end(fb,A,tmp2,request,reqcnt)

      implicit none
      integer fb
      complex*16  A(*)
      complex*16  tmp2(*)
      integer     request(*),reqcnt

#include "bafdecls.fh"
#include "errquit.fh"
#include "C3dB.fh"


*     **** indexing variables ****
      integer iq_to_i1(2,0:1)
      integer iq_to_i2(2,0:1)
      integer iz_to_i2(2,0:1)
      integer i1_start(2,0:1)
      integer i2_start(2,0:1)
      common / c_ptrans_blk2 / iq_to_i1,iq_to_i2,iz_to_i2,
     >                       i1_start,i2_start


*     **** local variables ***
      integer np,n2
      
      call Parallel3d_np_i(np)

*     **** unpack A(i) array ****
      n2 = int_mb(i2_start(1,fb)+np) - 1
      call C3dB_pfft_index2_copy(n2,int_mb(iq_to_i2(1,fb)),tmp2,A)
      call C3dB_pfft_index2_zero(nfft3d(1)-n2,int_mb(iz_to_i2(1,fb)),A)

      return
      end



*     ***********************************
*     *					*
*     *	   C3dB_c_ptranspose_ijk_start	*
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


      subroutine C3dB_c_ptranspose_ijk_start(fb,op,A,tmp1,tmp2,
     >                                       request,reqcnt,msgtype)

      implicit none
      integer fb,op
      complex*16  A(*)
      complex*16  tmp1(*),tmp2(*)
      integer     request(*),reqcnt,msgtype

#include "bafdecls.fh"
#include "errquit.fh"
#include "C3dB.fh"


*     **** indexing variables ****
      integer h_iq_to_i1(2,6,0:1)
      integer h_iq_to_i2(2,6,0:1)
      integer h_iz_to_i2(2,6,0:1)
      integer h_iz_to_i2_count(6,0:1)
      integer h_i1_start(2,6,0:1)
      integer h_i2_start(2,6,0:1)
      common / c_ptrans_blk_ijk / h_iq_to_i1,
     >                         h_iq_to_i2,
     >                         h_iz_to_i2,
     >                         h_iz_to_i2_count,
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
      integer c,n1
      integer it
      integer source
      integer msglen
      integer pfrom,pto
      integer taskid,np

*     **** external functions ****
      integer  Parallel3d_convert_taskid_i
      external Parallel3d_convert_taskid_i

      
      call Parallel3d_taskid_i(taskid)
      call Parallel3d_np_i(np)

*     **** pack A(i) array ****
      n1 = int_mb(h_i1_start(1,op,fb)+np) - 1
      call C3dB_pfft_index1_copy(n1,int_mb(h_iq_to_i1(1,op,fb)),A,tmp1)


*     **** it = 0, transpose data on same thread ****  
      msglen = int_mb(h_i2_start(1,op,fb)+2-1) 
     >       - int_mb(h_i2_start(1,op,fb)+1-1)
      call dcopy(2*msglen,
     >           tmp1(int_mb(h_i1_start(1,op,fb))),1,
     >           tmp2(int_mb(h_i2_start(1,op,fb))),1)



      do c=1,Nchannels(1)
*        **** receive packed array data ****
         if (int_mb(channel_type(1,1)+c-1) .eq. 1) then
            pfrom=int_mb(channel_proc(1,1)+c-1)
            it = mod((taskid+np-pfrom),np)

            source=pfrom
            msglen = (int_mb(h_i2_start(1,op,fb)+it+2-1)
     >             -  int_mb(h_i2_start(1,op,fb)+it+1-1))

            if (msglen.gt.0) then
               call RCV(9+MSGDBL,
     >                  tmp2(int_mb(h_i2_start(1,op,fb)+it+1-1)),
     >                  mdtob(2*msglen),rcv_len,
     >                  Parallel3d_convert_taskid_i(source),
     >                  rcv_proc,1)
            end if
         end if

*        **** send packed array to other processors ****
         if (int_mb(channel_type(1,1)+c-1) .eq. 0) then
            pto=int_mb(channel_proc(1,1)+c-1)
            it = mod((pto-taskid+np),np)

            msglen    = (int_mb(h_i1_start(1,op,fb)+it+2-1)
     >                -  int_mb(h_i1_start(1,op,fb)+it+1-1))

            if (msglen.gt.0) then
               call SND(9+MSGDBL,
     >                  tmp1(int_mb(h_i1_start(1,op,fb)+it+1-1)),
     >                  mdtob(2*msglen),
     >                  Parallel3d_convert_taskid_i(pto),
     >                  1)
            end if
         end if

      end do

      return
      end



*     ***********************************
*     *					*
*     *	   C3dB_c_ptranspose_ijk_end	*
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

      subroutine C3dB_c_ptranspose_ijk_end(fb,op,A,tmp2,
     >                                 request,reqcnt)

      implicit none
      integer fb,op
      complex*16  A(*)
      complex*16  tmp2(*)
      integer request(*),reqcnt

#include "bafdecls.fh"
#include "errquit.fh"
#include "C3dB.fh"


*     **** indexing variables ****
      integer h_iq_to_i1(2,6,0:1)
      integer h_iq_to_i2(2,6,0:1)
      integer h_iz_to_i2(2,6,0:1)
      integer h_iz_to_i2_count(6,0:1)
      integer h_i1_start(2,6,0:1)
      integer h_i2_start(2,6,0:1)
      common / c_ptrans_blk_ijk / h_iq_to_i1,
     >                         h_iq_to_i2,
     >                         h_iz_to_i2,
     >                         h_iz_to_i2_count,
     >                         h_i1_start,
     >                         h_i2_start

*     **** local variables ***
      integer n2,n3,np

      
      call Parallel3d_np_i(np)

*     **** unpack A(i) array ****
      n2 = int_mb(h_i2_start(1,op,fb)+np) - 1
      n3 = h_iz_to_i2_count(op,fb)
      call C3dB_pfft_index2_copy(n2,int_mb(h_iq_to_i2(1,op,fb)),tmp2,A)
      call C3dB_pfft_index2_zero(n3,int_mb(h_iz_to_i2(1,op,fb)),A)

      return
      end



*     ***********************************
*     *					*
*     *	   C3dB_c_ptranspose_ijk	*
*     *					*
*     ***********************************

      subroutine C3dB_c_ptranspose_ijk(fb,op,A,tmp1,tmp2)

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
      integer fb,op
      complex*16  A(*)
      complex*16  tmp1(*),tmp2(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "C3dB.fh"


*     **** indexing variables ****
      integer h_iq_to_i1(2,6,0:1)
      integer h_iq_to_i2(2,6,0:1)
      integer h_iz_to_i2(2,6,0:1)
      integer h_iz_to_i2_count(6,0:1)
      integer h_i1_start(2,6,0:1)
      integer h_i2_start(2,6,0:1)
      common / c_ptrans_blk_ijk / h_iq_to_i1,
     >                         h_iq_to_i2,
     >                         h_iz_to_i2,
     >                         h_iz_to_i2_count,
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
      integer c,n1,n2,n3
      integer it
      integer source
      integer msglen
      integer pfrom,pto
      integer taskid,np

*     **** external functions ****
      integer  Parallel3d_convert_taskid_i
      external Parallel3d_convert_taskid_i

      call Parallel3d_taskid_i(taskid)
      call Parallel3d_np_i(np)

      n1 = int_mb(h_i1_start(1,op,fb)+np) - 1
      n2 = int_mb(h_i2_start(1,op,fb)+np) - 1
      n3 = h_iz_to_i2_count(op,fb)


*     **** pack A(i) array ****
      call C3dB_pfft_index1_copy(n1,int_mb(h_iq_to_i1(1,op,fb)),A,tmp1)

*     **** it = 0, transpose data on same thread ****  
      msglen = int_mb(h_i2_start(1,op,fb)+2-1) 
     >       - int_mb(h_i2_start(1,op,fb)+1-1)
      call dcopy(2*msglen,
     >           tmp1(int_mb(h_i1_start(1,op,fb))),1,
     >           tmp2(int_mb(h_i2_start(1,op,fb))),1)


      do c=1,Nchannels(1)
*        **** receive packed array data ****
         if (int_mb(channel_type(1,1)+c-1) .eq. 1) then
            pfrom=int_mb(channel_proc(1,1)+c-1)
            it = mod((taskid+np-pfrom),np)

            source=pfrom
            msglen = (int_mb(h_i2_start(1,op,fb)+it+2-1)
     >             -  int_mb(h_i2_start(1,op,fb)+it+1-1))

            if (msglen.gt.0) then
               call RCV(9+MSGDBL,
     >                  tmp2(int_mb(h_i2_start(1,op,fb)+it+1-1)),
     >                  mdtob(2*msglen),rcv_len,
     >                  Parallel3d_convert_taskid_i(source),
     >                  rcv_proc,1)
            end if
         end if

*        **** send packed array to other processors ****
         if (int_mb(channel_type(1,1)+c-1) .eq. 0) then
            pto=int_mb(channel_proc(1,1)+c-1)
            it = mod((pto-taskid+np),np)

            msglen    = (int_mb(h_i1_start(1,op,fb)+it+2-1)
     >                -  int_mb(h_i1_start(1,op,fb)+it+1-1))

            if (msglen.gt.0) then
               call SND(9+MSGDBL,
     >                  tmp1(int_mb(h_i1_start(1,op,fb)+it+1-1)),
     >                  mdtob(2*msglen),
     >                  Parallel3d_convert_taskid_i(pto),
     >                  1)
            end if
         end if

      end do

*     **** unpack A(i) array ****
      call C3dB_pfft_index2_copy(n2,int_mb(h_iq_to_i2(1,op,fb)),tmp2,A)
      call C3dB_pfft_index2_zero(n3,int_mb(h_iz_to_i2(1,op,fb)),A)

      return
      end




*     ***********************************
*     *					*
*     *	   C3dB_channel_init	        *
*     *					*
*     ***********************************

      subroutine C3dB_channel_init(nb)
      implicit none
      integer nb

#include "bafdecls.fh"
#include "errquit.fh"
#include "C3dB.fh"

*     **** Used to avoid asynchronous communications ****
      integer Nchannels(NBLOCKS)
      integer channel_proc(2,NBLOCKS)
      integer channel_type(2,NBLOCKS)
      common / c_channel_blk / channel_proc,channel_type,Nchannels
      integer pair1(2),pair2(2)
      integer pair_step(2)
      integer step,Nstep,icount


*     **** local variables ****
      integer proc_to,proc_from
      integer pto,qto,np,taskid
      integer pfrom,qfrom
      integer phere,qhere
      integer index1,index2,itmp
      integer i,j,k,it
      logical value

*     **** external functions ****
      logical  full_ps_C3dB
      external full_ps_C3dB


      call Parallel3d_taskid_i(taskid)
      call Parallel3d_np_i(np)

*     **** Define Channels - which are used to avoid ****
*     **** asynchronous communications               ****
      value = BA_alloc_get(mt_int,(2*np),
     >        'channel_proc',channel_proc(2,nb),channel_proc(1,nb))
      value = value.and.
     >        BA_alloc_get(mt_int,(2*np),
     >        'channel_type',channel_type(2,nb),channel_type(1,nb))
      if (.not. value) call errquit('out of heap memory',0, MA_ERR)



      value = BA_push_get(mt_int,(np*(np-1)/2),
     >                    'pair1',pair1(2),pair1(1))
      value = value.and.
     >        BA_push_get(mt_int,(np*(np-1)/2),
     >                    'pair2',pair2(2),pair2(1))
      value = value.and.
     >        BA_push_get(mt_int,(np*(np-1)/2),
     >                'pair_step',pair_step(2),pair_step(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)

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
      do while(.not. full_ps_C3dB(int_mb(pair_step(1)),np))
         step=step+1
         do i=1, (np*(np-1)/2)
            if (int_mb(pair_step(1)+i-1).eq.(-1)) then
               value=.true.
               do j=1,(i-1)
                  if (int_mb(pair_step(1)+j-1).eq.step) then
                     if (int_mb(pair1(1)+i-1)
     >               .eq.int_mb(pair1(1)+j-1)) value=.false.
                     if (int_mb(pair1(1)+i-1)
     >               .eq.int_mb(pair2(1)+j-1)) value=.false.
                     if (int_mb(pair2(1)+i-1)
     >               .eq.int_mb(pair1(1)+j-1)) value=.false.
                     if (int_mb(pair2(1)+i-1)
     >               .eq.int_mb(pair2(1)+j-1)) value=.false.
                  end if
               end do
               if (value) int_mb(pair_step(1)+i-1) = step
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



      value=BA_pop_stack(pair_step(2))
      value=BA_pop_stack(pair2(2))
      value=BA_pop_stack(pair1(2))
      return
      end


      
      logical function full_ps_C3dB(ps,np)
      implicit none
      integer ps(*)
      integer np

      integer i
      logical value

      value=.true.
      do i=1,(np*(np-1)/2)
         if (ps(i).eq.(-1)) value=.false.
      end do

      full_ps_C3dB=value
      return
      end
             

*     ************************************
*     *                                  *
*     *      c_Balance_c_balance_start   *
*     *                                  *
*     ************************************

      subroutine c_Balance_c_balance_start(nb,A,request,reqcnt,msgtype)
      implicit none
      integer nb
      complex*16 A(*)
      integer    request(*),reqcnt,msgtype

      call c_Balance_c_balance(nb,A)
      return
      end

*     ************************************
*     *                                  *
*     *       c_Balance_c_balance_end    *
*     *                                  *
*     ************************************
*
      subroutine c_Balance_c_balance_end(nb,A,request,reqcnt)
      implicit none
      integer nb
      complex*16 A(*)
      integer    request(*),reqcnt

*     *** dummy routine ***
      return
      end

*     ************************************
*     *                                  *
*     *  c_Balance_c_unbalance_start     *
*     *                                  *
*     ************************************

      subroutine c_Balance_c_unbalance_start(nb,A,
     >                                       request,reqcnt,msgtype)
      implicit none
      integer nb
      complex*16 A(*)
      integer    request(*),reqcnt,msgtype

      call c_Balance_c_unbalance(nb,A)
      return
      end

*     ************************************
*     *                                  *
*     *     c_Balance_c_unbalance_end    *
*     *                                  *
*     ************************************

      subroutine c_Balance_c_unbalance_end(nb,A,request,reqcnt)
      implicit none
      integer nb
      complex*16 A(*)
      integer    request(*),reqcnt

*     *** dummy routine ***
      return
      end



*     ***********************************
*     *					*
*     *	        C3dB_SumAll		*	
*     *					*
*     ***********************************

      subroutine C3dB_SumAll(sum)
c     implicit none
      real*8  sum

#include "tcgmsg.fh"
#include "msgtypesf.h"

*     **** local variables ****
      integer np_i

*     **** external functions ****
      integer  Parallel3d_comm_i
      external Parallel3d_comm_i

      call Parallel3d_np_i(np_i)
      if (np_i.gt.1) then
         call GA_PGROUP_DGOP(Parallel3d_comm_i(),
     >                       9+MSGDBL,sum,1,'+')
      end if
      return
      end

*     ***********************************
*     *					*
*     *	        C3dB_ISumAll		*	
*     *					*
*     ***********************************

      subroutine C3dB_ISumAll(sum)
c     implicit none
      integer  sum

#include "tcgmsg.fh"
#include "msgtypesf.h"


*     **** local variables ****
      integer np_i

*     **** external functions ****
      integer  Parallel3d_comm_i
      external Parallel3d_comm_i

      call Parallel3d_np_i(np_i)
      if (np_i.gt.1) then
         call GA_PGROUP_IGOP(Parallel3d_comm_i(),
     >                       9+MSGINT,sum,1,'+')
      end if

      return
      end



*     ***********************************
*     *					*
*     *	        C3dB_Vector_SumAll	*	
*     *					*
*     ***********************************

      subroutine C3dB_Vector_SumAll(n,sum1)
c     implicit none
      integer n
      real*8  sum1(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"

*     **** temporary workspace ****
      integer np_i

*     **** external functions ****
      integer  Parallel3d_comm_i
      external Parallel3d_comm_i

      call nwpw_timing_start(2)

      call Parallel3d_np_i(np_i)
      if (np_i.gt.1) then
         call GA_PGROUP_DGOP(Parallel3d_comm_i(),
     >                       9+MSGDBL,sum1,n,'+')
      end if

      call nwpw_timing_end(2)

      return
      end


*     ***********************************
*     *					*
*     *	        C3dB_Vector_ISumAll	*	
*     *					*
*     ***********************************

      subroutine C3dB_Vector_ISumAll(n,sum)
c     implicit none
      integer n
      integer  sum(*)




#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "bafdecls.fh"
#include "errquit.fh"

*     **** local variables ****
      integer np_i

*     **** external functions ****
      integer  Parallel3d_comm_i
      external Parallel3d_comm_i

      call nwpw_timing_start(2)

      call Parallel3d_np_i(np_i)
      if (np_i.gt.1) then
         call GA_PGROUP_IGOP(Parallel3d_comm_i(),
     >                       9+MSGINT,sum,n,'+')
      end if

      call nwpw_timing_end(2)
      return
      end




*     ***********************************
*     *					*
*     *	       C3dB_(c,r)_read  	*	
*     *					*
*     ***********************************

      subroutine C3dB_c_read(nb,iunit,A,tmp,jcol,kcol)
      implicit none 
      integer nb
      integer iunit
      complex*16 A(*)
      complex*16 tmp(*)
      integer jcol,kcol

#include "bafdecls.fh"
#include "errquit.fh"

#include "C3dB.fh"

#include "tcgmsg.fh"
#include "msgtypesf.h"
      integer rcv_len,rcv_proc

      
*     *** local variables ***
      logical value,fillcolumn,fillzone
      integer MASTER,taskid
      parameter(MASTER=0)
      integer p_to, p_here,q
      integer index,j,k
      integer source,msglen
      integer tmp1(2),tmp2(2)

      integer taskid_j,np_j
      integer taskid_k,np_k
      integer ii,jj,kk,jstart,jend,kstart,kend

*     **** external functions ****
      integer  Parallel3d_convert_taskid_ijk
      external Parallel3d_convert_taskid_ijk

      call Parallel_taskid(taskid)

      call Parallel3d_np_j(np_j)
      call Parallel3d_taskid_j(taskid_j)
      if (jcol.lt.0) then
         jstart = 0
         jend = np_j-1
         fillcolumn = .true.
      else
         jstart = jcol
         jend   = jcol
         fillcolumn = (taskid_j.eq.jcol)
      endif
      call Parallel3d_np_k(np_k)
      call Parallel3d_taskid_k(taskid_k)
      if (kcol.lt.0) then
         kstart = 0
         kend = np_k-1
         fillzone = .true.
      else
         kstart = kcol
         kend   = kcol
         fillzone = (taskid_k.eq.kcol)
      endif


      !**********************
      !**** slab mapping ****
      !**********************
      if (mapping.eq.1) then
*     **** master node reads from file and distributes ****
      if (taskid.eq.MASTER) then
         do k=1,nz(nb)

            call dread(iunit,tmp,2*nx(nb)*ny(nb))

            call C3dB_ktoqp(nb,k,q,ii)
            do kk=kstart,kend
            do jj=jstart,jend
               p_to = Parallel3d_convert_taskid_ijk(ii,jj,kk)
               if (p_to.eq.MASTER) then
                  index = 1 + (q-1)*nx(nb)*ny(nb)
                  call zcopy(nx(nb)*ny(nb),tmp,1,A(index),1)
               else
                  msglen = nx(nb)*ny(nb)
                  call SND(9+MSGDBL,tmp,mdtob(2*msglen),p_to,1)
               end if
            end do
            end do
         end do

*     **** not master node ****
      else if (fillcolumn.and.fillzone) then
         do k=1,nz(nb)
            call C3dB_ktoqp(nb,k,q,ii)
            p_here = Parallel3d_convert_taskid_ijk(ii,taskid_j,taskid_k)
            if (p_here.eq.taskid) then
               msglen = nx(nb)*ny(nb)
               source  = MASTER
               call RCV(9+MSGDBL,tmp,mdtob(2*msglen),rcv_len,
     >                  source,rcv_proc,1)
               index = 1 + (q-1)*nx(nb)*ny(nb) 
               call zcopy(nx(nb)*ny(nb),tmp,1,A(index),1)
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

            call dread(iunit,tmp,2*nx(nb))

            q    = int_mb(q_map1(1,nb)+(j-1)+(k-1)*ny(nb))
            ii   = int_mb(p_map1(1,nb)+(j-1)+(k-1)*ny(nb))
            do kk=kstart,kend
            do jj=jstart,jend
               p_to = Parallel3d_convert_taskid_ijk(ii,jj,kk)
               if (p_to.eq.MASTER) then
                  index = 1 + (q-1)*nx(nb)
                  call zcopy(nx(nb),tmp,1,A(index),1)
               else
                  msglen = nx(nb)
                  call SND(9+MSGDBL,tmp,mdtob(2*msglen),p_to,1)
               end if
            end do
            end do
         end do
         end do

*     **** not master node ****
      else if (fillcolumn.and.fillzone) then
         do k=1,nz(nb)
         do j=1,ny(nb)

            q      = int_mb(q_map1(1,nb)+(j-1)+(k-1)*ny(nb))
            ii     = int_mb(p_map1(1,nb)+(j-1)+(k-1)*ny(nb))
            p_here = Parallel3d_convert_taskid_ijk(ii,taskid_j,taskid_k)
            if (p_here.eq.taskid) then
               msglen = nx(nb)
               source  = MASTER
               call RCV(9+MSGDBL,tmp,mdtob(2*msglen),rcv_len,
     >                  source,rcv_proc,1)
               index = 1 + (q-1)*nx(nb)
               call zcopy(nx(nb),tmp,1,A(index),1)
            end if
         end do
         end do
      end if

      if (fillcolumn.and.fillzone) then
*        **** allocate temporary space  ****
         value =BA_push_get(mt_dcpl,(nfft3d(nb)),'tmp1',tmp1(2),tmp1(1))
         value =value.and.
     >          BA_push_get(mt_dcpl,(nfft3d(nb)),'tmp2',tmp2(2),tmp2(1))
         if (.not. value) 
     >      call errquit('C3dB_c_read:out of stack memory',0,MA_ERR)
*
         call C3dB_c_transpose_ijk(nb,5,A,
     >                          dcpl_mb(tmp1(1)),
     >                          dcpl_mb(tmp2(1)))  !*** map1to3 operation ***

*        **** deallocate temporary space  ****
         value =           BA_pop_stack(tmp2(2))
         value = value.and.BA_pop_stack(tmp1(2))
         if (.not. value) 
     >      call errquit('C3dB_c_read:error popping stack',0,MA_ERR)
      endif

      end if

*     **** WAIT ****
c      call ga_sync()

      return
      end

      subroutine C3dB_r_read(nb,iunit,A,tmp,jcol,kcol,dotrans)
      implicit none 
      integer nb
      integer iunit
      real*8  A(*)
      real*8  tmp(*)
      integer jcol,kcol
      logical dotrans

#include "bafdecls.fh"
#include "errquit.fh"
#include "C3dB.fh"

c#include "mpif.h"

#include "tcgmsg.fh"
#include "msgtypesf.h"
      integer rcv_len,rcv_proc

      
*     *** local variables ***
      integer MASTER,taskid
      parameter(MASTER=0)

      logical value,fillcolumn,fillzone
      integer p_to, p_here,q
      integer j,k,index,index2
      integer source,msglen
      integer tmp1(2),tmp2(2)

      integer taskid_j,np_j
      integer taskid_k,np_k
      integer ii,jj,kk,jstart,jend,kstart,kend
      integer msgtype
      real*8 dum
c      integer mpierr,mstatus(MPI_STATUS_SIZE)

*     **** external functions ****
      integer  Parallel3d_convert_taskid_ijk
      external Parallel3d_convert_taskid_ijk

      call Parallel_taskid(taskid)

      call Parallel3d_np_j(np_j)
      call Parallel3d_taskid_j(taskid_j)
      if (jcol.lt.0) then
         jstart = 0
         jend = np_j-1
         fillcolumn = .true.
      else
         jstart = jcol
         jend   = jcol
         fillcolumn = (taskid_j.eq.jcol)
      endif
      call Parallel3d_np_k(np_k)
      call Parallel3d_taskid_k(taskid_k)
      if (kcol.lt.0) then
         kstart = 0
         kend = np_k-1
         fillzone = .true.
      else
         kstart = kcol
         kend   = kcol
         fillzone = (taskid_k.eq.kcol)
      endif

      !**********************
      !**** slab mapping ****
      !**********************
      if (mapping.eq.1) then
*     **** master node reads from file and distributes ****
      if (taskid.eq.MASTER) then
         do k=1,nz(nb)

            call dread(iunit,tmp,(nx(nb))*ny(nb))

            call C3dB_ktoqp(nb,k,q,ii)
            do kk=kstart,kend
            do jj=jstart,jend
               p_to = Parallel3d_convert_taskid_ijk(ii,jj,kk)
               if (p_to.eq.MASTER) then
                  index = 1 + (q-1)*nx(nb)*ny(nb)
                  call dcopy(nx(nb)*ny(nb),tmp,1,A(index),1)
               else
                  msglen = nx(nb)*ny(nb)
                  call SND(9+MSGDBL,tmp,mdtob(msglen),p_to,1)
               end if
            end do
            end do

         end do

*     **** not master node ****
      else if (fillcolumn.and.fillzone) then
         do k=1,nz(nb)
            call C3dB_ktoqp(nb,k,q,ii)
            p_here = Parallel3d_convert_taskid_ijk(ii,taskid_j,taskid_k)
            if (p_here.eq.taskid) then
               msglen  = nx(nb)*ny(nb)
               source   = MASTER
               call RCV(9+MSGDBL,tmp,mdtob(msglen),rcv_len,
     >                  source,rcv_proc,1)
               index = 1 + (q-1)*nx(nb)*ny(nb)
               call dcopy(nx(nb)*ny(nb),tmp,1,A(index),1)
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
            do kk=kstart,kend
            do jj=jstart,jend
               p_to = Parallel3d_convert_taskid_ijk(ii,jj,kk)
               if (p_to.eq.MASTER) then
                  index = (q-1)*nx(nb) + 1
                  call dcopy(nx(nb),tmp,1,A(index),1)
               else
                  msglen = nx(nb)
                  call SND(9+MSGDBL,tmp,mdtob(msglen),p_to,1)
               end if
            end do
            end do
         end do
         end do

*     **** not master node ****
      else if (fillcolumn.and.fillzone) then
         do k=1,nz(nb)
         do j=1,ny(nb)

            q      = int_mb(q_map1(1,nb)+(j-1)+(k-1)*ny(nb))
            ii     = int_mb(p_map1(1,nb)+(j-1)+(k-1)*ny(nb))
            p_here = Parallel3d_convert_taskid_ijk(ii,taskid_j,taskid_k)
            if (p_here.eq.taskid) then
               msglen = nx(nb)
               source  = MASTER
               call RCV(9+MSGDBL,tmp,mdtob(msglen),rcv_len,
     >                  source,rcv_proc,1)
               index = (q-1)*nx(nb) + 1
               call dcopy(nx(nb),tmp,1,A(index),1)
            end if
         end do
         end do
       end if

      if (fillcolumn.and.fillzone.and.dotrans) then
*        **** allocate temporary space  ****
         value =BA_push_get(mt_dcpl,(nfft3d(nb)),'tmp1',tmp1(2),tmp1(1))
         value =value.and.
     >          BA_push_get(mt_dcpl,(nfft3d(nb)),'tmp2',tmp2(2),tmp2(1))
         if (.not. value) 
     >      call errquit('C3dB_r_read:out of stack memory',0,MA_ERR)
         call C3dB_r_transpose_ijk(nb,5,A,
     >                          dcpl_mb(tmp1(1)),
     >                          dcpl_mb(tmp2(1)))  !*** map1to3 operation ***
*        **** deallocate temporary space  ****
         value =           BA_pop_stack(tmp2(2))
         value = value.and.BA_pop_stack(tmp1(2))
         if (.not. value) 
     >      call errquit('C3dB_r_read:error popping stack',0,MA_ERR)
      end if

      end if

*     **** WAIT ****
c      call ga_sync()

      return
      end


*     ***********************************
*     *					*
*     *	       C3dB_(c,r)_write	        *	
*     *					*
*     ***********************************

      subroutine C3dB_c_write(nb,iunit,A,tmp,jcol,kcol)
      implicit none 
      integer nb
      integer iunit
      complex*16 A(*)
      complex*16 tmp(*)
      integer jcol,kcol

#include "bafdecls.fh"
#include "errquit.fh"
#include "C3dB.fh"

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

      integer ii,taskid_j,taskid_k

*     **** external functions ****
      integer  Parallel3d_convert_taskid_ijk
      external Parallel3d_convert_taskid_ijk

      call Parallel_taskid(taskid)
      call Parallel3d_taskid_j(taskid_j)
      call Parallel3d_taskid_k(taskid_k)


      !**********************
      !**** slab mapping ****
      !**********************
      if (mapping.eq.1) then
*     **** master node reads from file and distributes ****
      if (taskid.eq.MASTER) then
         do k=1,nz(nb)

            call C3dB_ktoqp(nb,k,q,ii)
            p_from = Parallel3d_convert_taskid_ijk(ii,jcol,kcol)

            if (p_from.eq.MASTER) then
               index = 1 + (q-1)*nx(nb)*ny(nb) 
               call zcopy(nx(nb)*ny(nb),A(index),1,tmp,1)
            else
               msglen  = nx(nb)*ny(nb)
               status  = msglen
               source  = p_from
               dum = 99
               dum_msglen = 1
               call SND(9+MSGINT,dum,mitob(dum_msglen),source,1)
               call RCV(9+MSGDBL,tmp,mdtob(2*msglen),rcv_len,
     >                       source,rcv_proc,1)
            end if
            call dwrite(iunit,tmp,2*nx(nb)*ny(nb)) 
         end do

*     **** not master node ****
      else
         do k=1,nz(nb)
         
            call C3dB_ktoqp(nb,k,q,ii)
            p_here = Parallel3d_convert_taskid_ijk(ii,jcol,kcol)
            if (p_here.eq.taskid) then

               index = 1 + (q-1)*nx(nb)*ny(nb)
               call zcopy(nx(nb)*ny(nb),A(index),1,tmp,1)

               msglen  = nx(nb)*ny(nb)
               dest    = MASTER
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

      if ((taskid_j.eq.jcol).and.(taskid_k.eq.kcol)) then
*        **** allocate temporary space  ****
         value = BA_push_get(mt_dcpl,(nfft3d(nb)),'tmp1',
     >                       tmp1(2),tmp1(1))
         value = value.and.
     >         BA_push_get(mt_dcpl,(nfft3d(nb)),'tmp2',tmp2(2),tmp2(1))
         if (.not. value) 
     >      call errquit('C3dB_c_write:out of stack memory',0,MA_ERR)
*
         call C3dB_c_transpose_ijk(nb,6,A,
     >                          dcpl_mb(tmp1(1)),
     >                          dcpl_mb(tmp2(1)))  !*** map3to1 operation ***

*        **** deallocate temporary space  ****
         value = BA_pop_stack(tmp2(2))
         value = value.and.BA_pop_stack(tmp1(2))
         if (.not. value) 
     >      call errquit('C3dB_c_write:error popping stack',0, MA_ERR)
      end if

c      call ga_sync()

*     **** master node reads from file and distributes ****
      if (taskid.eq.MASTER) then
         do k=1,nz(nb)
         do j=1,ny(nb)

            q      = int_mb(q_map1(1,nb)+(j-1)+(k-1)*ny(nb))
            ii     = int_mb(p_map1(1,nb)+(j-1)+(k-1)*ny(nb))
            p_from = Parallel3d_convert_taskid_ijk(ii,jcol,kcol)
            if (p_from.eq.MASTER) then

               index = 1 + (q-1)*nx(nb)
               call zcopy(nx(nb),A(index),1,tmp,1)
            else
               msglen  = nx(nb)
               status  = msglen
               source  = p_from

               dum = 99
               dum_msglen = 1
               call SND(9+MSGINT,dum,mitob(dum_msglen),source,1)
               call RCV(9+MSGDBL,tmp,mdtob(2*msglen),rcv_len,
     >                       source,rcv_proc,1)
            end if

            call dwrite(iunit,tmp,2*nx(nb))
         end do
         end do

*     **** not master node ****
      else
         do k=1,nz(nb)
         do j=1,ny(nb)

            q      = int_mb(q_map1(1,nb)+(j-1)+(k-1)*ny(nb))
            ii     = int_mb(p_map1(1,nb)+(j-1)+(k-1)*ny(nb))
            p_here = Parallel3d_convert_taskid_ijk(ii,jcol,kcol)
            if (p_here.eq.taskid) then

               index = 1 + (q-1)*nx(nb)
               call zcopy(nx(nb),A(index),1,tmp,1)

               msglen  = nx(nb)
               dest    = MASTER
               dum_msglen = 1
               call RCV(9+MSGINT,dum,mitob(dum_msglen),rcv_len,
     >                       dest,rcv_proc,1)
               call SND(9+MSGDBL,tmp,mdtob(2*msglen),dest,1)
            end if

         end do
         end do
      end if

      end if

*     **** WAIT ****
c      call ga_sync()

      return
      end

      subroutine C3dB_r_write(nb,iunit,A,tmp,jcol,kcol)
      implicit none 
      integer nb
      integer iunit
      real*8     A(*)
      real*8     tmp(*)
      integer jcol,kcol

#include "bafdecls.fh"
#include "C3dB.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"

      integer rcv_len,rcv_proc
      
*     *** local variables ***
      integer MASTER,taskid
      parameter(MASTER=0)
      integer p_from, p_here,q
      integer j,k,index,index2
      integer dest,source,status,msglen

      integer ii,taskid_j,taskid_k

*     **** external functions ****
      integer  Parallel3d_convert_taskid_ijk
      external Parallel3d_convert_taskid_ijk

      call Parallel_taskid(taskid)
      call Parallel3d_taskid_j(taskid_j)
      call Parallel3d_taskid_k(taskid_k)


      !**********************
      !**** slab mapping ****
      !**********************
      if (mapping.eq.1) then
*     **** master node reads from file and distributes ****
      if (taskid.eq.MASTER) then
         do k=1,nz(nb)

            call C3dB_ktoqp(nb,k,q,ii)
            p_from = Parallel3d_convert_taskid_ijk(ii,jcol,kcol)
            if (p_from.eq.MASTER) then
               index = 1 + (q-1)*nx(nb)*ny(nb)
               call dcopy(nx(nb)*ny(nb),A(index),1,tmp,1)
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
         
            call C3dB_ktoqp(nb,k,q,ii)
            p_here = Parallel3d_convert_taskid_ijk(ii,jcol,kcol)
            if (p_here.eq.taskid) then

               index = 1 + (q-1)*nx(nb)*ny(nb)
               call dcopy(nx(nb)*ny(nb),A(index),1,tmp,1)
 
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
            ii     = int_mb(p_map1(1,nb)+(j-1)+(k-1)*ny(nb))
            p_from = Parallel3d_convert_taskid_ijk(ii,jcol,kcol)

            if (p_from.eq.MASTER) then
              index = 1 + (q-1)*nx(nb)
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
            p_here = Parallel3d_convert_taskid_ijk(ii,jcol,kcol)

            if (p_here.eq.taskid) then

               index = 1 + (q-1)*nx(nb)
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
c      call ga_sync()
      return
      end




c $Id$
