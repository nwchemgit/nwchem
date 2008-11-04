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

#include "mafdecls.fh"
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

#include "mafdecls.fh"
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

#include "mafdecls.fh"
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

#include "mafdecls.fh"
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

#include "mafdecls.fh"
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

#include "mafdecls.fh"
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

#include "mafdecls.fh"
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

#include "mafdecls.fh"
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

#include "mafdecls.fh"
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

#include "mafdecls.fh"
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
      value = MA_alloc_get(mt_int,(2*np),
     >        'channel_proc',channel_proc(2,nb),channel_proc(1,nb))
      value = value.and.
     >        MA_alloc_get(mt_int,(2*np),
     >        'channel_type',channel_type(2,nb),channel_type(1,nb))
      if (.not. value) call errquit('out of heap memory',0, MA_ERR)



      value = MA_push_get(mt_int,(np*(np-1)/2),
     >                    'pair1',pair1(2),pair1(1))
      value = value.and.
     >        MA_push_get(mt_int,(np*(np-1)/2),
     >                    'pair2',pair2(2),pair2(1))
      value = value.and.
     >        MA_push_get(mt_int,(np*(np-1)/2),
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



      value=MA_pop_stack(pair_step(2))
      value=MA_pop_stack(pair2(2))
      value=MA_pop_stack(pair1(2))
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

#include "mafdecls.fh"
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
#include "mafdecls.fh"
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


