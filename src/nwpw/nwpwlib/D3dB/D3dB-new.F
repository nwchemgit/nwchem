#define TCGMSG
#define NBLOCKS 2

*
* $Id: D3dB-new.F,v 1.1 2001-08-30 17:41:40 bylaska Exp $
*

*     ***********************************************************
*     *								*
*     *   		   D3dB library				*
*     *			(MPI implemenation)			*
*     *								*
*     *   Author - Eric Bylaska					*
*     *   date   - 3/23/96					*
*     *								*
*     ***********************************************************

*	The D3dB (distributed three-dimensional block) library is to 
* be used for handling three kinds of data structures.  The first 
* data structure, denoted by "r", is a double precision array of
* length (nx(nb)+2)*ny(nb)*nz.  The second data structure, denoted by "c", is
* a double complex array of length of (nx(nb)/2+1)*ny(nb)*nz.  The third data
* structure, denoted by "t", is a dobule precision array of length
* (nx(nb)/2+1)*ny(nb)*nz.
*
*	The three data structures are distributed across threads, p, in
* the k (i.e. nz(nb)) dimension using a cyclic decomposition.  So that
* a "r" array A is defined as double precision A(nx(nb)+2,ny(nb),nq(nb)) on
* each thread.  
*
*	Where 
*		np = number of threads
*		nq(nb) = ceil(nz(nb)/np).
*		0 <= p < np
*		1 <= q <= nq(nb)
*		1 <= k <= nz(nb)
* 
* 	The mapping of k -> q is defined as:
*	
*		k = ((q-1)*np + p) + 1
*		q = ((k-1) - p)/np + 1
*		p = (k-1) mod np
*
*  Libraries used: mpi, blas, fftpack, and compressed_io
*
*  common blocks used in this library:
*
*       integer nq,nx(NBLOCKS),ny,nz(nb)
*	common	/ D3dB / nq,nx,ny,nz
*
*	integer q_map(NFFT3),p_map(NFFT3),k_map(NFFT3)
* 	common /D3dB_mapping / q_map,p_map,k_map
*
*     integer iq_to_i1((NFFT1/2+1)*NFFT2*NSLABS)
*     integer iq_to_i2((NFFT1/2+1)*NFFT2*NSLABS)
*     integer i1_start(NPROCS+1)
*     integer i2_start(NPROCS+1)
*     common / trans_blk / iq_to_i1,iq_to_i2,i1_start,i2_start

*     **** local variables ****

*     ***********************************
*     *					*
*     *	       Mapping_Init		*	
*     *					*
*     ***********************************

      subroutine Mapping_Init(nb)
      implicit none
      integer nb
     
#include "mafdecls.fh"

      integer nq(NBLOCKS),nx(NBLOCKS),ny(NBLOCKS),nz(NBLOCKS)
      common	/ D3dB / nq,nx,ny,nz

c     integer q_map(NFFT3),p_map(NFFT3),k_map(NFFT3)
      integer q_map(2,NBLOCKS),p_map(2,NBLOCKS),k_map(2,NBLOCKS)
      common /D3dB_mapping / q_map,p_map,k_map

      integer k,q,p
*     integer kn
      integer taskid,np
      logical value
      

*     **** allocate q_map,p_map,k_map
      value = MA_alloc_get(mt_int,nz(nb),'q_map',q_map(2,nb),
     >                                       q_map(1,nb))
      value = MA_alloc_get(mt_int,nz(nb),'p_map',p_map(2,nb),
     >                                       p_map(1,nb))
      value = MA_alloc_get(mt_int,nz(nb),'k_map',k_map(2,nb),
     >                                       k_map(1,nb))

      call Parallel_np(np)
      call Parallel_taskid(taskid)


*     ****************************
*     ****** cyclic mapping ******
*     ****************************
      p = 0
      q = 1
      do k=1,nz(nb)
         int_mb(q_map(1,nb)+k-1) = q
         int_mb(p_map(1,nb)+k-1) = p
         if (p .eq. taskid) nq(nb) = q

         p        = p+1
         if (p .ge. np) then
            p = 0
            q = q + 1
         end if
      end do

*     **************************
*     ***** block mapping ******
*     **************************
*
*     **** make sure nz(nb) is a multiple of np ****
*     kn = mod(nz(nb),np)
*     if (kn.ne.0) then
*        kn=(nz(nb)/np)+1
*     else
*        kn=(nz(nb)/np)
*     end if
*
*     p=0
*     q=1
*     do k=1,nz(nb)
*        int_mb(q_map(1,nb)+k-1) = q
*        int_mb(p_map(1,nb)+k-1) = p
*        if (p .eq. taskid) nq(nb) = q
*
*        q=q+1
*        if (q .gt. (kn)) then
*           q = 1
*           p = p + 1
*        end if
*     end do


      do k=1,nz(nb)
         if (int_mb(p_map(1,nb)+k-1) .eq. taskid) then
c           k_map(q_map(k)) = k
            int_mb(k_map(1,nb)+int_mb(q_map(1,nb)+k-1)-1) = k
         end if
      end do
 
      return
      end

*     ***********************************
*     *					*
*     *	          D3dB_end   		*	
*     *					*
*     ***********************************
      subroutine D3dB_end(nb)
      implicit none
      integer nb

#include "mafdecls.fh"

c     integer q_map(NFFT3),p_map(NFFT3),k_map(NFFT3)
      integer q_map(2,NBLOCKS),p_map(2,NBLOCKS),k_map(2,NBLOCKS)
      common /D3dB_mapping / q_map,p_map,k_map

c     integer iq_to_i1((NFFT1/2+1)*NFFT2*NSLABS)
c     integer iq_to_i2((NFFT1/2+1)*NFFT2*NSLABS)
c     integer i1_start(NFFT3+1)
c     integer i2_start(NFFT3+1)
      integer iq_to_i1(2,NBLOCKS)
      integer iq_to_i2(2,NBLOCKS)
      integer i1_start(2,NBLOCKS)
      integer i2_start(2,NBLOCKS)
      common / trans_blk / iq_to_i1,iq_to_i2,i1_start,i2_start

      integer Nchannels(NBLOCKS)
      integer channel_proc(2,NBLOCKS)
      integer channel_type(2,NBLOCKS)
      common / channel_blk / channel_proc,channel_type,Nchannels

      logical value

      call D3dB_timereverse_end(nb)

      value = MA_free_heap(q_map(2,nb))
      value = MA_free_heap(p_map(2,nb))
      value = MA_free_heap(k_map(2,nb))
      value = MA_free_heap(i1_start(2,nb))
      value = MA_free_heap(i2_start(2,nb))
      value = MA_free_heap(iq_to_i1(2,nb))
      value = MA_free_heap(iq_to_i2(2,nb))
      value = MA_free_heap(channel_proc(2,nb))
      value = MA_free_heap(channel_type(2,nb))

      return
      end
         
*     ***********************************
*     *					*
*     *	          D3dB_qtok   		*	
*     *					*
*     ***********************************

      subroutine D3dB_qtok(nb,q,k)      
      implicit none
      integer nb
      integer q,k

#include "mafdecls.fh"
      
c     integer q_map(NFFT3),p_map(NFFT3),k_map(NFFT3)
      integer q_map(2,NBLOCKS),p_map(2,NBLOCKS),k_map(2,NBLOCKS)
      common /D3dB_mapping / q_map,p_map,k_map
      
c     k = k_map(q)
      k = int_mb(k_map(1,nb)+q-1)

      return
      end

*     ***********************************
*     *					*
*     *	          D3dB_ktoqp  		*	
*     *					*
*     ***********************************

      subroutine D3dB_ktoqp(nb,k,q,p)      
      implicit none
      integer nb
      integer k,q,p

#include "mafdecls.fh"
      
c     integer q_map(NFFT3),p_map(NFFT3),k_map(NFFT3)
      integer q_map(2,NBLOCKS),p_map(2,NBLOCKS),k_map(2,NBLOCKS)
      common /D3dB_mapping / q_map,p_map,k_map
      
c     q = q_map(k)
c     p = p_map(k)

      q = int_mb(q_map(1,nb)+k-1)
      p = int_mb(p_map(1,nb)+k-1)
      return
      end

*     ***********************************
*     *					*
*     *	        D3dB_nfft3d		*
*     *					*
*     ***********************************

      subroutine D3dB_nfft3d(nb,nfft3d)
      implicit none
      integer nb
      integer nfft3d

      integer nq(NBLOCKS),nx(NBLOCKS),ny(NBLOCKS),nz(NBLOCKS)
      common	/ D3dB / nq,nx,ny,nz

      nfft3d = (nx(nb)/2+1)*ny(nb)*nq(nb)
      return
      end

*     ***********************************
*     *					*
*     *	        D3dB_n2ft3d		*
*     *					*
*     ***********************************

      subroutine D3dB_n2ft3d(nb,n2ft3d)
      implicit none
      integer nb
      integer n2ft3d

      integer nq(NBLOCKS),nx(NBLOCKS),ny(NBLOCKS),nz(NBLOCKS)
      common	/ D3dB / nq,nx,ny,nz

      n2ft3d = (nx(nb)+2)*ny(nb)*nq(nb)
      return
      end


*     ***********************************
*     *					*
*     *	        D3dB_nq			*	
*     *					*
*     ***********************************

      subroutine D3dB_nq(nb,nqtmp)
      implicit none
      integer nb
      integer nqtmp

      integer nq(NBLOCKS),nx(NBLOCKS),ny(NBLOCKS),nz(NBLOCKS)
      common	/ D3dB / nq,nx,ny,nz

      nqtmp = nq(nb)

      return 
      end

*     ***********************************
*     *					*
*     *	        D3dB_nx			*	
*     *					*
*     ***********************************
     
      subroutine D3dB_nx(nb,nxtmp)
      implicit none
      integer nb
      integer nxtmp
      
      integer nq(NBLOCKS),nx(NBLOCKS),ny(NBLOCKS),nz(NBLOCKS)
      common	/ D3dB / nq,nx,ny,nz

      nxtmp = nx(nb)
      return
      end

*     ***********************************
*     *					*
*     *	        D3dB_ny			*	
*     *					*
*     ***********************************

      subroutine D3dB_ny(nb,nytmp)
      implicit none
      integer nb
      integer nytmp
      
      integer nq(NBLOCKS),nx(NBLOCKS),ny(NBLOCKS),nz(NBLOCKS)
      common	/ D3dB / nq,nx,ny,nz

      nytmp = ny(nb)
      return
      end

*     ***********************************
*     *					*
*     *	        D3dB_nz			*	
*     *					*
*     ***********************************

      subroutine D3dB_nz(nb,nztmp)
      implicit none
      integer nb
      integer nztmp
      
      integer nq(NBLOCKS),nx(NBLOCKS),ny(NBLOCKS),nz(NBLOCKS)
      common	/ D3dB / nq,nx,ny,nz

      nztmp = nz(nb)
      return
      end


*     ***********************************
*     *					*
*     *	        D3dB_Init		*	
*     *					*
*     ***********************************

      subroutine D3dB_Init(nb,nx_in,ny_in,nz_in)
      implicit none
      integer nb
      integer nx_in,ny_in,nz_in

      integer nq(NBLOCKS),nx(NBLOCKS),ny(NBLOCKS),nz(NBLOCKS)
      common	/ D3dB / nq,nx,ny,nz

     

*     ***** initialize D3dB common block *****
      nx(nb)     = nx_in
      ny(nb)     = ny_in
      nz(nb)     = nz_in


*     **** do other initializations ****
      call Mapping_Init(nb)
      call D3dB_c_transpose_jk_init(nb)
      call D3dB_c_timereverse_init(nb)

      return
      end


*     ***********************************
*     *					*
*     *	        D3dB_SumAll		*	
*     *					*
*     ***********************************

      subroutine D3dB_SumAll(sum)
c     implicit none
      real*8  sum

#ifdef MPI
      include 'mpif.h'
      integer mpierr
#endif
#ifdef TCGMSG
#include "tcgmsg.fh"
#include "msgtypesf.h"
#endif

      integer MASTER
      parameter (MASTER=0)
      integer msglen
      real*8 sumall,sumt
    
c     msglen = 8
      msglen = 1

      sumt = sum
#ifdef MPI
c     call mp_combine(sumt,sumall,msglen,d_vadd,allgrp)
      call MPI_REDUCE(sumt,sumall,msglen,MPI_DOUBLE_PRECISION,
     >                       MPI_SUM,MASTER,MPI_COMM_WORLD,mpierr)
      call MPI_BCAST(sumall,msglen,MPI_DOUBLE_PRECISION,
     >               MASTER,MPI_COMM_WORLD,mpierr)
#endif
#ifdef TCGMSG
      call GA_DGOP(9+MSGDBL,sumt,1,'+')
      sumall=sumt
#endif

      sum = sumall
      return
      end


*     ***********************************
*     *					*
*     *	        D3dB_ISumAll		*	
*     *					*
*     ***********************************

      subroutine D3dB_ISumAll(sum)
c     implicit none
      integer  sum

#ifdef MPI
      include 'mpif.h'
      integer mpierr
#endif
#ifdef TCGMSG
#include "tcgmsg.fh"
#include "msgtypesf.h"
#endif

      integer MASTER
      parameter (MASTER=0)
      integer msglen
      integer sumall,sumt
    
c     msglen = 8
      msglen = 1

      sumt = sum
#ifdef MPI
c     call mp_combine(sumt,sumall,msglen,d_vadd,allgrp)
      call MPI_REDUCE(sumt,sumall,msglen,MPI_INTEGER,
     >                       MPI_SUM,MASTER,MPI_COMM_WORLD,mpierr)
      call MPI_BCAST(sumall,msglen,MPI_INTEGER,
     >               MASTER,MPI_COMM_WORLD,mpierr)
#endif
#ifdef TCGMSG
      call GA_IGOP(9+MSGINT,sumt,1,'+')
      sumall=sumt
#endif

      sum = sumall
      return
      end



*     ***********************************
*     *					*
*     *	        D3dB_(c,r,t)_Zero	*	
*     *					*
*     ***********************************

      subroutine D3dB_c_Zero(nb,A)     
      implicit none 
      integer nb
      complex*16 A(*)

      integer nq(NBLOCKS),nx(NBLOCKS),ny(NBLOCKS),nz(NBLOCKS)
      common	/ D3dB / nq,nx,ny,nz

      call dcopy((nx(nb)+2)*ny(nb)*nq(nb),0.0d0,0,A,1)
      return
      end


      subroutine D3dB_r_Zero(nb,A)     
      implicit none 
      integer nb
      real*8  A(*)

      integer nq(NBLOCKS),nx(NBLOCKS),ny(NBLOCKS),nz(NBLOCKS)
      common	/ D3dB / nq,nx,ny,nz

      call dcopy((nx(nb)+2)*ny(nb)*nq(nb),0.0d0,0,A,1)
      return
      end



*     ***********************************
*     *					*
*     *	        D3dB_(c,r,t)_Copy	*	
*     *					*
*     ***********************************

      subroutine D3dB_c_Copy(nb,A,B)     
      implicit none 
      integer nb
      complex*16 A(*)
      complex*16 B(*)

      integer nq(NBLOCKS),nx(NBLOCKS),ny(NBLOCKS),nz(NBLOCKS)
      common	/ D3dB / nq,nx,ny,nz

      call dcopy((nx(nb)+2)*ny(nb)*nq(nb),A,1,B,1)
      return
      end

      subroutine D3dB_r_Copy(nb,A,B)     
      implicit none 
      integer nb
      real*8 A(*)
      real*8 B(*)

      integer nq(NBLOCKS),nx(NBLOCKS),ny(NBLOCKS),nz(NBLOCKS)
      common	/ D3dB / nq,nx,ny,nz

      call dcopy((nx(nb)+2)*ny(nb)*nq(nb),A,1,B,1)
      return
      end

      subroutine D3dB_t_Copy(nb,A,B)     
      implicit none 
      integer nb
      real*8 A(*)
      real*8 B(*)

      integer nq(NBLOCKS),nx(NBLOCKS),ny(NBLOCKS),nz(NBLOCKS)
      common	/ D3dB / nq,nx,ny,nz

      call dcopy((nx(nb)/2+1)*ny(nb)*nq(nb),A,1,B,1)
      return
      end

*     ***********************************
*     *					*
*     *	        D3dB_cr_fft3b		*
*     *					*
*     ***********************************

      subroutine D3dB_cr_fft3b(nb,A)

*****************************************************
*                                                   *
*      This routine performs the operation of       *
*      a three dimensional complex to complex       *
*      inverse fft                                  *
*           A(nx,ny(nb),nz(nb)) <- FFT3^(-1)[A(kx,ky,kz)]   * 
*                                                   *
*      Entry - 					    *
*              A: a column distribuded 3d block     *
*              tmp: tempory work space must be at   *
*                    least the size of (complex)    *
*                    (nfft*nfft + 1) + 10*nfft      * 
*                                                   *
*       Exit - A is transformed and the imaginary   *
*              part of A is set to zero             *
*       uses - D3dB_c_transpose_jk, dcopy           *
*                                                   *
*****************************************************

      implicit none
      integer nb
      complex*16  A(*)

#include "mafdecls.fh"


      integer nq(NBLOCKS),nx(NBLOCKS),ny(NBLOCKS),nz(NBLOCKS)
      common	/ D3dB / nq,nx,ny,nz

*     *** local variables ***
      integer i,j,k,q,indx

c     complex*16  tmp1(*)
c     complex*16  tmp2(*)
c     real*8      tmp3(*)
      integer nfft3d
      integer tmp1(2),tmp2(2),tmp3(2)
      logical value



*     ***** timing variables *****
*     1 - fft_time
*     2 - ABC_time
*     3 - L_time
*     4 - vlocal_time
*     5 - vnl_time
*     6 - vh_time
      real*8 tim1,tim2
      real*8 times(50)
      common / timing / times


      call current_second(tim1)

*     ***** allocate temporary space ****
      call D3dB_nfft3d(nb,nfft3d)
      value = MA_push_get(mt_dcpl,(nfft3d), 'ffttmp1',tmp1(2),tmp1(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,(nfft3d), 'ffttmp2',tmp2(2),tmp2(1))
      value = value.and.
     >        MA_push_get(mt_dbl,(2*nfft3d),'ffttmp3',tmp3(2),tmp3(1))
      if (.not. value) call errquit('out of stack memory',0)

*     ********************************************
*     ***         Do a transpose of A          ***
*     ***      A(kx,kz,ky) <- A(kx,ky,kz)      ***
*     ********************************************
c     call D3dB_c_transpose_jk(nb,A,dcpl_mb(tmp2(1)),dbl_mb(tmp3(1)))

*     *************************************************
*     ***     do fft along kz dimension             ***
*     ***   A(kx,nz(nb),ky) <- fft1d^(-1)[A(kx,kz,ky)]  ***
*     *************************************************
      call dcffti(nz(nb),dcpl_mb(tmp1(1)))
      do q=1,nq(nb)
      do i=1,(nx(nb)/2+1)
         do k=1,nz(nb)
            indx = i + (k-1)*(nx(nb)/2+1) + (q-1)*(nx(nb)/2+1)*nz(nb)
            dcpl_mb(tmp2(1)+k-1) = A(indx)
         end do
         call dcfftb(nz(nb),dcpl_mb(tmp2(1)),dcpl_mb(tmp1(1)))
         do k=1,nz(nb)
            indx = i + (k-1)*(nx(nb)/2+1) + (q-1)*(nx(nb)/2+1)*nz(nb)
            A(indx) = dcpl_mb(tmp2(1)+k-1)
         end do
      end do
      end do

*     ********************************************
*     ***         Do a transpose of A          ***
*     ***      A(kx,ky,nz(nb)) <- A(kx,nz(nb),ky)      ***
*     ********************************************
      call D3dB_c_transpose_jk(nb,A,dcpl_mb(tmp2(1)),dbl_mb(tmp3(1)))

*     *************************************************
*     ***     do fft along ky dimension             ***
*     ***   A(kx,ny(nb),nz(nb)) <- fft1d^(-1)[A(kx,ky,nz(nb))]  ***
*     *************************************************
      call dcffti(ny(nb),dcpl_mb(tmp1(1)))
      do q=1,nq(nb)
      do i=1,(nx(nb)/2+1)
         do j=1,ny(nb)
            indx = i + (j-1)*(nx(nb)/2+1) + (q-1)*(nx(nb)/2+1)*ny(nb)
            dcpl_mb(tmp2(1)+j-1) = A(indx)
         end do
         call dcfftb(ny(nb),dcpl_mb(tmp2(1)),dcpl_mb(tmp1(1)))
         do j=1,ny(nb)
            indx = i + (j-1)*(nx(nb)/2+1) + (q-1)*(nx(nb)/2+1)*ny(nb)
            A(indx) = dcpl_mb(tmp2(1)+j-1)
         end do
      end do
      end do

*     *************************************************
*     ***     do fft along kx dimension             ***
*     ***   A(nx(nb),ny(nb),nz(nb)) <- fft1d^(-1)[A(kx,ny(nb),nz(nb))]  ***
*     *************************************************
      call drffti(nx(nb),dcpl_mb(tmp1(1)))
      do q=1,nq(nb)
      do j=1,ny(nb)
         indx = 1 + (j-1)*(nx(nb)/2+1) + (q-1)*(nx(nb)/2+1)*ny(nb)
         call dcopy((nx(nb)+2),A(indx),1,dbl_mb(tmp3(1)),1)
         do i=2,nx(nb)
            dbl_mb(tmp3(1)+i-1) = dbl_mb(tmp3(1)+i)
         end do
         call drfftb(nx(nb),dbl_mb(tmp3(1)),dcpl_mb(tmp1(1)))
         dbl_mb(tmp3(1)+nx(nb)) = 0.0d0
         dbl_mb(tmp3(1)+nx(nb)+1) = 0.0d0
         call dcopy((nx(nb)+2),dbl_mb(tmp3(1)),1,A(indx),1)
      end do
      end do

    
*     **** deallocate temporary space  ****
      value = MA_pop_stack(tmp3(2))
      value = MA_pop_stack(tmp2(2))
      value = MA_pop_stack(tmp1(2))

      call current_second(tim2)
      times(1) = times(1) + (tim2-tim1)
      return
      end





*     ***********************************
*     *					*
*     *	        D3dB_rc_fft3f		*
*     *					*
*     ***********************************

      subroutine D3dB_rc_fft3f(nb,A)

*****************************************************
*                                                   *
*      This routine performs the operation of       *
*      a three dimensional complex to complex fft   *
*           A(kx,ky,kz) <- FFT3[A(nx(nb),ny(nb),nz(nb))]        * 
*                                                   *
*      Entry - 					    *
*              A: a column distribuded 3d block     *
*              tmp: tempory work space must be at   *
*                    least the size of (complex)    *
*                    (nfft*nfft + 1) + 10*nfft      * 
*                                                   *
*       Exit - A is transformed                     *
*                                                   *
*       uses - transpose1 subroutine                *
*                                                   *
*****************************************************

      implicit none
      integer nb
      complex*16  A(*)

#include "mafdecls.fh"

      integer nq(NBLOCKS),nx(NBLOCKS),ny(NBLOCKS),nz(NBLOCKS)
      common	/ D3dB / nq,nx,ny,nz

*     *** local variables ***
      integer i,j,k,q,indx

      integer nfft3d
      integer tmp1(2),tmp2(2),tmp3(2)
      logical value


*     ***** timing variables *****
*     1 - fft_time
*     2 - ABC_time
*     3 - L_time
*     4 - vlocal_time
*     5 - vnl_time
*     6 - vh_time
      real*8 tim1,tim2
      real*8 times(50)
      common / timing / times

      call current_second(tim1)

*     ***** allocate temporary space ****
      call D3dB_nfft3d(nb,nfft3d)
      value = MA_push_get(mt_dcpl,(nfft3d),'tmp1',tmp1(2),tmp1(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,(nfft3d),'tmp2',tmp2(2),tmp2(1))
      value = value.and.
     >        MA_push_get(mt_dbl,(2*nfft3d),'tmp3',tmp3(2),tmp3(1))
      if (.not. value) call errquit('out of stack memory',0)


*     ********************************************
*     ***     do fft along nx(nb) dimension        ***
*     ***   A(kx,ny(nb),nz(nb)) <- fft1d[A(nx(nb),ny(nb),nz(nb))]  ***
*     ********************************************
      call drffti(nx(nb),dcpl_mb(tmp1(1)))
      do q=1,nq(nb)
      do j=1,ny(nb)
         indx = 1 + (j-1)*(nx(nb)/2+1) + (q-1)*(nx(nb)/2+1)*ny(nb)
         call dcopy((nx(nb)+2),A(indx),1,dbl_mb(tmp3(1)),1)
         call drfftf(nx(nb),dbl_mb(tmp3(1)),dcpl_mb(tmp1(1)))
         do i=nx(nb),2,-1
            dbl_mb(tmp3(1)+i) = dbl_mb(tmp3(1)+i-1)
         end do
         dbl_mb(tmp3(1)+1)    = 0.0d0
         dbl_mb(tmp3(1)+nx(nb)+1) = 0.0d0
         call dcopy((nx(nb)+2),dbl_mb(tmp3(1)),1,A(indx),1)
         
      end do
      end do

*     ********************************************
*     ***     do fft along ny(nb) dimension        ***
*     ***   A(kx,ky,nz(nb)) <- fft1d[A(kx,ny(nb),nz(nb))]  ***
*     ********************************************
      call dcffti(ny(nb),dcpl_mb(tmp1(1)))
      do q=1,nq(nb)
      do i=1,(nx(nb)/2+1)
         do j=1,ny(nb)
            indx = i + (j-1)*(nx(nb)/2+1) + (q-1)*(nx(nb)/2+1)*ny(nb)
            dcpl_mb(tmp2(1)+j-1) = A(indx)
         end do
         call dcfftf(ny(nb),dcpl_mb(tmp2(1)),dcpl_mb(tmp1(1)))
         do j=1,ny(nb)
            indx = i + (j-1)*(nx(nb)/2+1) + (q-1)*(nx(nb)/2+1)*ny(nb)
            A(indx) = dcpl_mb(tmp2(1)+j-1)
         end do
      end do
      end do


*     ********************************************
*     ***         Do a transpose of A          ***
*     ***      A(ky,nz(nb),ky) <- A(kx,ky,nz(nb))      ***
*     ********************************************
      call D3dB_c_transpose_jk(nb,A,dcpl_mb(tmp2(1)),dbl_mb(tmp3(1)))


*     ********************************************
*     ***     do fft along nz(nb) dimension        ***
*     ***   A(kx,kz,ky) <- fft1d[A(kx,nz(nb),ky)]  ***
*     ********************************************
      call dcffti(nz(nb),dcpl_mb(tmp1(1)))
      do q=1,nq(nb)
      do i=1,(nx(nb)/2+1)
         do k=1,nz(nb)
            indx = i + (k-1)*(nx(nb)/2+1) + (q-1)*(nx(nb)/2+1)*ny(nb)
            dcpl_mb(tmp2(1)+k-1) = A(indx)
         end do
         call dcfftf(nz(nb),dcpl_mb(tmp2(1)),dcpl_mb(tmp1(1)))
         do k=1,nz(nb)
            indx = i + (k-1)*(nx(nb)/2+1) + (q-1)*(nx(nb)/2+1)*ny(nb)
            A(indx) = dcpl_mb(tmp2(1)+k-1)
         end do
      end do
      end do

*     ********************************************
*     ***         Do a transpose of A          ***
*     ***      A(kx,ky,kz) <- A(kx,kz,ky)      ***
*     ********************************************
c     call D3dB_c_transpose_jk(nb,A,dcpl_mb(tmp2(1)),dbl_mb(tmp3(1)))


*     **** deallocate temporary space  ****
      value = MA_pop_stack(tmp3(2))
      value = MA_pop_stack(tmp2(2))
      value = MA_pop_stack(tmp1(2))

      call current_second(tim2)
      times(1) = times(1) + (tim2-tim1)
      return
      end




*     ***********************************
*     *					*
*     *	       D3dB_(c,r,t)_Read 	*	
*     *					*
*     ***********************************

      subroutine D3dB_c_Read(nb,iunit,A,tmp)
      implicit none 
      integer nb
      integer iunit
      complex*16 A(*)
      complex*16 tmp(*)

      integer nq(NBLOCKS),nx(NBLOCKS),ny(NBLOCKS),nz(NBLOCKS)
      common	/ D3dB / nq,nx,ny,nz

#ifdef MPI
      include 'mpif.h'
      integer mpierr
      integer msgtype,status
#endif
#ifdef TCGMSG
#include "tcgmsg.fh"
#include "msgtypesf.h"
      integer rcv_len,rcv_proc
#endif
      
*     *** local variables ***
      integer MASTER,taskid
      parameter(MASTER=0)
      integer p_to, p_here,q
      integer index,k
      integer source,msglen

      call Parallel_taskid(taskid)

*     **** master node reads from file and distributes ****
      if (taskid.eq.MASTER) then
         do k=1,nz(nb)

            call dread(iunit,tmp,(nx(nb)+2)*ny(nb))

            call D3dB_ktoqp(nb,k,q,p_to)
            if (p_to.eq.MASTER) then
               index = (q-1)*(nx(nb)/2+1)*ny(nb) + 1
               call dcopy((nx(nb)+2)*ny(nb),tmp,1,A(index),1)
            else
               msglen = (nx(nb)/2+1)*ny(nb)
#ifdef MPI
               msgtype= 33
               call MPI_SEND(tmp,msglen,MPI_DOUBLE_COMPLEX,
     >                       p_to,msgtype,MPI_COMM_WORLD,mpierr)
#endif
#ifdef TCGMSG
               call SND(9+MSGDBL,tmp,mdtob(2*msglen),p_to,1)
#endif
            end if
         end do

*     **** not master node ****
      else
         do k=1,nz(nb)
            call D3dB_ktoqp(nb,k,q,p_here)
            if (p_here.eq.taskid) then

               msglen = (nx(nb)/2+1)*ny(nb)
               source  = MASTER
#ifdef MPI
               msgtype = 33
               call MPI_RECV(tmp,msglen,MPI_DOUBLE_COMPLEX,
     >                       source,msgtype,MPI_COMM_WORLD,
     >                       status,mpierr)
#endif
#ifdef TCGMSG
               call RCV(9+MSGDBL,tmp,mdtob(2*msglen),rcv_len,
     >                  source,rcv_proc,1)
#endif

               index = (q-1)*(nx(nb)/2+1)*ny(nb) + 1
               call dcopy((nx(nb)+2)*ny(nb),tmp,1,A(index),1)
               
            end if
         end do
      end if

*     **** wait ****
#ifdef MPI
      call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
#endif

      return
      end

      subroutine D3dB_r_Read(nb,iunit,A,tmp)
      implicit none 
      integer nb
      integer iunit
      real*8  A(*)
      real*8  tmp(*)

      integer nq(NBLOCKS),nx(NBLOCKS),ny(NBLOCKS),nz(NBLOCKS)
      common	/ D3dB / nq,nx,ny,nz

#ifdef MPI
      include 'mpif.h'
      integer mpierr
      integer msgtype,status
#endif
#ifdef TCGMSG
#include "tcgmsg.fh"
#include "msgtypesf.h"
      integer rcv_len,rcv_proc
#endif
      
*     *** local variables ***
      integer MASTER,taskid
      parameter(MASTER=0)
      integer p_to, p_here,q
      integer j,k,index,index2
      integer source,msglen

      call Parallel_taskid(taskid)

*     **** master node reads from file and distributes ****
      if (taskid.eq.MASTER) then
         do k=1,nz(nb)

            call dread(iunit,tmp,(nx(nb))*ny(nb))

            call D3dB_ktoqp(nb,k,q,p_to)
            if (p_to.eq.MASTER) then

               do j=1,ny(nb)
                  index = (q-1)*(nx(nb)+2)*ny(nb) 
     >                  + (j-1)*(nx(nb)+2) + 1
                  index2 = (j-1)*nx(nb) + 1
                  call dcopy(nx(nb),tmp(index2),1,A(index),1)
                  A(index+nx(nb)) = 0.0d0
                  A(index+nx(nb)+1) = 0.0d0
               end do
              
            else
               msglen = (nx(nb))*ny(nb)
#ifdef MPI
               msgtype= 33
               call MPI_SEND(tmp,msglen,MPI_DOUBLE_PRECISION,
     >                       p_to,msgtype,MPI_COMM_WORLD,mpierr)
#endif
#ifdef TCGMSG
               call SND(9+MSGDBL,tmp,mdtob(msglen),p_to,1)
#endif
            end if
         end do

*     **** not master node ****
      else
         do k=1,nz(nb)
            call D3dB_ktoqp(nb,k,q,p_here)
            if (p_here.eq.taskid) then

               msglen = (nx(nb))*ny(nb)
               source  = MASTER
#ifdef MPI
               msgtype = 33
               call MPI_RECV(tmp,msglen,MPI_DOUBLE_PRECISION,
     >                       source,msgtype,MPI_COMM_WORLD,
     >                       status,mpierr)
#endif
#ifdef TCGMSG
               call RCV(9+MSGDBL,tmp,mdtob(msglen),rcv_len,
     >                       source,rcv_proc,1)
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

*     **** wait ****
#ifdef MPI
      call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
#endif
      return
      end

      subroutine D3dB_t_Read(nb,iunit,A,tmp)
      implicit none 
      integer nb
      integer iunit
      real*8  A(*)
      real*8  tmp(*)

      integer nq(NBLOCKS),nx(NBLOCKS),ny(NBLOCKS),nz(NBLOCKS)
      common	/ D3dB / nq,nx,ny,nz

#ifdef MPI
      include 'mpif.h'
      integer mpierr
      integer msgtype,status
#endif
#ifdef TCGMSG
#include "tcgmsg.fh"
#include "msgtypesf.h"
      integer rcv_len,rcv_proc
#endif

      
*     *** local variables ***
      integer MASTER,taskid
      parameter(MASTER=0)
      integer p_to, p_here,q
      integer k,index
      integer source,msglen

      call Parallel_taskid(taskid)

*     **** master node reads from file and distributes ****
      if (taskid.eq.MASTER) then
         do k=1,nz(nb)

            call dread(iunit,tmp,(nx(nb)/2+1)*ny(nb))

            call D3dB_ktoqp(nb,k,q,p_to)
            if (p_to.eq.MASTER) then
               index = (q-1)*(nx(nb)/2+1)*ny(nb) + 1
               call dcopy((nx(nb)/2+1)*ny(nb),tmp,1,A(index),1)
            else
               msglen = (nx(nb)/2+1)*ny(nb)
#ifdef MPI
               msgtype= 33
               call MPI_SEND(tmp,msglen,MPI_DOUBLE_PRECISION,
     >                       p_to,msgtype,MPI_COMM_WORLD,mpierr)
#endif
#ifdef TCGMSG
               call SND(9+MSGDBL,tmp,mdtob(msglen),p_to,1)
#endif
            end if
         end do

*     **** not master node ****
      else
         do k=1,nz(nb)
            call D3dB_ktoqp(nb,k,q,p_here)
            if (p_here.eq.taskid) then

               msglen = (nx(nb)/2+1)*ny(nb)
               source  = MASTER
#ifdef MPI
               msgtype = 33
               call MPI_RECV(tmp,msglen,MPI_DOUBLE_PRECISION,
     >                       source,msgtype,MPI_COMM_WORLD,
     >                       status,mpierr)
#endif
#ifdef TCGMSG
               call RCV(9+MSGDBL,tmp,mdtob(msglen),rcv_len,
     >                       source,rcv_proc,1)
#endif

               index = (q-1)*(nx(nb)/2+1)*ny(nb) + 1
               call dcopy((nx(nb)/2+1)*ny(nb),tmp,1,A(index),1)
               
            end if
         end do
      end if

*     **** wait ****
#ifdef MPI
      call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
#endif
      return
      end



*     ***********************************
*     *					*
*     *	       D3dB_(c,r,t)_SMul 	*	
*     *					*
*     ***********************************

*  This routine performs the operation	C = scale * A
* where scale is a real*8 number.

      subroutine D3dB_c_SMul(nb,scale,A,C)     
      implicit none 
      integer    nb
      real*8     scale
      complex*16 A(*)
      complex*16 C(*)

      integer nq(NBLOCKS),nx(NBLOCKS),ny(NBLOCKS),nz(NBLOCKS)
      common	/ D3dB / nq,nx,ny,nz

      integer i

      do i=1,(nx(nb)/2+1)*ny(nb)*nq(nb)
         C(i) = scale*A(i)
      end do
      return
      end


      subroutine D3dB_r_SMul(nb,scale,A,C)     
      implicit none 
      integer nb
      real*8     scale
      real*8 A(*)
      real*8 C(*)

      integer nq(NBLOCKS),nx(NBLOCKS),ny(NBLOCKS),nz(NBLOCKS)
      common	/ D3dB / nq,nx,ny,nz

      integer i

      do i=1,(nx(nb)+2)*ny(nb)*nq(nb)
         C(i) = scale*A(i)
      end do
      return
      end 


      subroutine D3dB_t_SMul(nb,scale,A,C)     
      implicit none 
      integer nb
      real*8 scale
      real*8 A(*)
      real*8 C(*)

      integer nq(NBLOCKS),nx(NBLOCKS),ny(NBLOCKS),nz(NBLOCKS)
      common	/ D3dB / nq,nx,ny,nz

      integer i

      do i=1,(nx(nb)/2+1)*ny(nb)*nq(nb)
         C(i) = scale*A(i)
      end do
      return
      end 

      subroutine D3dB_c_ZMul(nb,scale,A,C)     
      implicit none 
      integer    nb
      complex*16 scale
      complex*16 A(*)
      complex*16 C(*)

      integer nq(NBLOCKS),nx(NBLOCKS),ny(NBLOCKS),nz(NBLOCKS)
      common	/ D3dB / nq,nx,ny,nz

      integer i

      do i=1,(nx(nb)/2+1)*ny(nb)*nq(nb)
         C(i) = scale*A(i)
      end do
      return
      end



*     ***********************************
*     *					*
*     *	       D3dB_ct_Sqr	 	*	
*     *					*
*     ***********************************

*  This routine performs the operation	C = A * A

      subroutine D3dB_ct_Sqr(nb,A,C)     
      implicit none 
      integer    nb
      complex*16 A(*)
      real*8     C(*)

      integer nq(NBLOCKS),nx(NBLOCKS),ny(NBLOCKS),nz(NBLOCKS)
      common	/ D3dB / nq,nx,ny,nz

      integer i

      do i=1,(nx(nb)/2+1)*ny(nb)*nq(nb)
         C(i) = dble(A(i))**2 + dimag(A(i))**2
      end do
      return
      end

*     ***********************************
*     *					*
*     *	       D3dB_rr_Sqr	 	*	
*     *					*
*     ***********************************

      subroutine D3dB_rr_Sqr(nb,A,C)     
      implicit none 
      integer nb
      real*8 A(*)
      real*8 C(*)

      integer nq(NBLOCKS),nx(NBLOCKS),ny(NBLOCKS),nz(NBLOCKS)
      common	/ D3dB / nq,nx,ny,nz

      integer i

      do i=1,(nx(nb)+2)*ny(nb)*nq(nb)
         C(i) = A(i)**2
      end do
      return
      end

*     ***********************************
*     *					*
*     *	       D3dB_rr_Sqrt	 	*	
*     *					*
*     ***********************************

      subroutine D3dB_rr_Sqrt(nb,A,C)     
      implicit none 
      integer nb
      real*8 A(*)
      real*8 C(*)

      integer nq(NBLOCKS),nx(NBLOCKS),ny(NBLOCKS),nz(NBLOCKS)
      common	/ D3dB / nq,nx,ny,nz

      integer i

      do i=1,(nx(nb)+2)*ny(nb)*nq(nb)
         C(i) = dsqrt(A(i))
      end do
      return
      end


*     ***********************************
*     *					*
*     *	       D3dB_tt_Sqr	 	*	
*     *					*
*     ***********************************

      subroutine D3dB_tt_Sqr(nb,A,C)     
      implicit none 
      integer nb
      real*8 A(*)
      real*8 C(*)

      integer nq(NBLOCKS),nx(NBLOCKS),ny(NBLOCKS),nz(NBLOCKS)
      common	/ D3dB / nq,nx,ny,nz

      integer i

      do i=1,(nx(nb)/2+1)*ny(nb)*nq(nb)
         C(i) = A(i)**2
      end do
      return
      end

#ifdef MPI
*     ***********************************
*     *					*
*     *	   D3dB_c_transpose_jk_init	*
*     *					*
*     ***********************************

      subroutine D3dB_c_transpose_jk_init(nb)
      implicit none
      integer nb

#include "mafdecls.fh"

      integer nq(NBLOCKS),nx(NBLOCKS),ny(NBLOCKS),nz(NBLOCKS)
      common	/ D3dB / nq,nx,ny,nz

c     integer iq_to_i1((NFFT1/2+1)*NFFT2*NSLABS)
c     integer iq_to_i2((NFFT1/2+1)*NFFT2*NSLABS)
c     integer i1_start(NFFT3+1)
c     integer i2_start(NFFT3+1)
      integer iq_to_i1(2,NBLOCKS)
      integer iq_to_i2(2,NBLOCKS)
      integer i1_start(2,NBLOCKS)
      integer i2_start(2,NBLOCKS)
      common / trans_blk / iq_to_i1,iq_to_i2,i1_start,i2_start

*     **** local variables ****
      integer proc_to,proc_from
      integer pto,qto,np,taskid
      integer pfrom,qfrom
      integer phere,qhere
      integer index1,index2,itmp
      integer i,j,k,it
      logical value

*     **** allocate trans_blk common block ****
      value = MA_alloc_get(mt_int,((nx(nb)/2+1)*ny(nb)*nq(nb)),
     >                     'iq_to_i1',iq_to_i1(2,nb),iq_to_i1(1,nb))
      value = MA_alloc_get(mt_int,((nx(nb)/2+1)*ny(nb)*nq(nb)),
     >                     'iq_to_i2',iq_to_i2(2,nb),iq_to_i2(1,nb))

      value = MA_alloc_get(mt_int,(nz(nb)+1),
     >                     'i1_start',i1_start(2,nb),i1_start(1,nb))
      value = MA_alloc_get(mt_int,(nz(nb)+1),
     >                     'i2_start',i2_start(2,nb),i2_start(1,nb))

      call Parallel_taskid(taskid)
      call Parallel_np(np)

      index1 = 1 
      index2 = 1
      do it=0,np-1
         proc_to   = mod(taskid+it,np)
         proc_from = mod(taskid-it+np,np)
c        i1_start(it+1) = index1
c        i2_start(it+1) = index2
         int_mb(i1_start(1,nb)+it) = index1
         int_mb(i2_start(1,nb)+it) = index2

         do k=1,nz(nb)
         do j=1,ny(nb)

*           **** packing scheme **** 
            call D3dB_ktoqp(nb,k,qhere,phere)
            call D3dB_ktoqp(nb,j,qto,pto)
            if ((phere.eq.taskid).and.(pto.eq.proc_to)) then
               do i=1,(nx(nb)/2+1)
                  itmp = i + (j-1)*(nx(nb)/2+1) 
     >                     + (qhere-1)*(nx(nb)/2+1)*ny(nb)
c                 iq_to_i1(itmp) = index1
                  int_mb(iq_to_i1(1,nb)+itmp-1) = index1
                  index1 = index1 + 1
               end do
            end if
             
*           **** unpacking scheme ****
            call D3dB_ktoqp(nb,j,qhere,phere)
            call D3dB_ktoqp(nb,k,qfrom,pfrom)
            if ((phere.eq.taskid).and.(pfrom.eq.proc_from)) then
               do i=1,(nx(nb)/2+1)
                  itmp = i + (k-1)*(nx(nb)/2+1) 
     >                     + (qhere-1)*(nx(nb)/2+1)*ny(nb)
c                 iq_to_i2(itmp) = index2
                  int_mb(iq_to_i2(1,nb)+itmp-1) = index2
                  index2 = index2 + 1
               end do
            end if
         end do
         end do
      end do
c     i1_start(np+1) = index1
c     i2_start(np+1) = index2
      int_mb(i1_start(1,nb)+np) = index1
      int_mb(i2_start(1,nb)+np) = index2

      return
      end
             
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

      integer nq(NBLOCKS),nx(NBLOCKS),ny(NBLOCKS),nz(NBLOCKS)
      common	/ D3dB / nq,nx,ny,nz


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


      include 'mpif.h'
      integer mpierr

*     **** local variables ***
      integer i,q
      integer proc_to,proc_from,it
      integer size,msgtype,source
c     integer msglen,type,nbytes
c     integer msgid
      integer msglen,type
      integer status(512),request(512),reqcnt
      integer taskid,np
      
      call Parallel_taskid(taskid)
      call Parallel_np(np)

*     **** pack A(i) array ****
      do i=1,(nx(nb)/2+1)*ny(nb)*nq(nb)
         tmp1(int_mb(iq_to_i1(1,nb)+i-1)) = A(i)
      end do

*     **** it = 0, transpose data on same thread ****  
      msglen = int_mb(i2_start(1,nb)+2-1) - int_mb(i2_start(1,nb)+1-1)
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
*           reqcnt = reqcnt + 1
            call MPI_SEND(tmp1(int_mb(i1_start(1,nb)+it+1-1)),
     >                     msglen,MPI_DOUBLE_COMPLEX,
     >                     proc_to,msgtype,MPI_COMM_WORLD,mpierr)
         end if

      end do

   
*     **** wait for completion of mp_send, also do a sync ****
      call MPI_WAITALL(reqcnt,request,status,mpierr)

*     **** unpack A(i) array ****
      do i=1,(nx(nb)/2+1)*ny(nb)*nq(nb)
         A(i) = tmp2(int_mb(iq_to_i2(1,nb)+i-1))
      end do
      
      return
      end
#endif


#ifdef TCGMSG
*     ***********************************
*     *					*
*     *	   D3dB_c_transpose_jk_init	*
*     *					*
*     ***********************************

      subroutine D3dB_c_transpose_jk_init(nb)
      implicit none
      integer nb

#include "mafdecls.fh"

      integer nq(NBLOCKS),nx(NBLOCKS),ny(NBLOCKS),nz(NBLOCKS)
      common	/ D3dB / nq,nx,ny,nz

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
      logical  full_ps
      external full_ps


*     **** allocate trans_blk common block ****
      value = MA_alloc_get(mt_int,((nx(nb)/2+1)*ny(nb)*nq(nb)),
     >                     'iq_to_i1',iq_to_i1(2,nb),iq_to_i1(1,nb))
      value = MA_alloc_get(mt_int,((nx(nb)/2+1)*ny(nb)*nq(nb)),
     >                     'iq_to_i2',iq_to_i2(2,nb),iq_to_i2(1,nb))

      value = MA_alloc_get(mt_int,(nz(nb)+1),
     >                     'i1_start',i1_start(2,nb),i1_start(1,nb))
      value = MA_alloc_get(mt_int,(nz(nb)+1),
     >                     'i2_start',i2_start(2,nb),i2_start(1,nb))

      call Parallel_taskid(taskid)
      call Parallel_np(np)

      index1 = 1 
      index2 = 1
      do it=0,np-1
         proc_to   = mod(taskid+it,np)
         proc_from = mod(taskid-it+np,np)
c        i1_start(it+1) = index1
c        i2_start(it+1) = index2
         int_mb(i1_start(1,nb)+it) = index1
         int_mb(i2_start(1,nb)+it) = index2

         do k=1,nz(nb)
         do j=1,ny(nb)

*           **** packing scheme **** 
            call D3dB_ktoqp(nb,k,qhere,phere)
            call D3dB_ktoqp(nb,j,qto,pto)
            if ((phere.eq.taskid).and.(pto.eq.proc_to)) then
               do i=1,(nx(nb)/2+1)
                  itmp = i + (j-1)*(nx(nb)/2+1) 
     >                     + (qhere-1)*(nx(nb)/2+1)*ny(nb)
c                 iq_to_i1(itmp) = index1
                  int_mb(iq_to_i1(1,nb)+itmp-1) = index1
                  index1 = index1 + 1
               end do
            end if
             
*           **** unpacking scheme ****
            call D3dB_ktoqp(nb,j,qhere,phere)
            call D3dB_ktoqp(nb,k,qfrom,pfrom)
            if ((phere.eq.taskid).and.(pfrom.eq.proc_from)) then
               do i=1,(nx(nb)/2+1)
                  itmp = i + (k-1)*(nx(nb)/2+1) 
     >                     + (qhere-1)*(nx(nb)/2+1)*ny(nb)
c                 iq_to_i2(itmp) = index2
                  int_mb(iq_to_i2(1,nb)+itmp-1) = index2
                  index2 = index2 + 1
               end do
            end if
         end do
         end do
      end do
c     i1_start(np+1) = index1
c     i2_start(np+1) = index2
      int_mb(i1_start(1,nb)+np) = index1
      int_mb(i2_start(1,nb)+np) = index2


*
*     **** Define Channels - which are used to avoid ****
*     **** asynchronous communications               ****
      
      value = MA_alloc_get(mt_int,(2*np),
     >        'channel_proc',channel_proc(2,nb),channel_proc(1,nb))
      value = value.and.
     >        MA_alloc_get(mt_int,(2*np),
     >        'channel_type',channel_type(2,nb),channel_type(1,nb))
      if (.not. value) call errquit('out of heap memory',0)



      value = MA_push_get(mt_int,(np*(np-1)/2),
     >                    'pair1',pair1(2),pair1(1))
      value = value.and.
     >        MA_push_get(mt_int,(np*(np-1)/2),
     >                    'pair2',pair2(2),pair2(1))
      value = value.and.
     >        MA_push_get(mt_int,(np*(np-1)/2),
     >                'pair_step',pair_step(2),pair_step(1))
      if (.not. value) call errquit('out of stack memory',0)

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
      do while(.not. full_ps(int_mb(pair_step(1)),np))
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

      
      logical function full_ps(ps,np)
      implicit none
      integer ps(*)
      integer np

      integer i
      logical value

      value=.true.
      do i=1,(np*(np-1)/2)
         if (ps(i).eq.(-1)) value=.false.
      end do

      full_ps=value
      return
      end
             
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

      integer nq(NBLOCKS),nx(NBLOCKS),ny(NBLOCKS),nz(NBLOCKS)
      common	/ D3dB / nq,nx,ny,nz


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
      do i=1,(nx(nb)/2+1)*ny(nb)*nq(nb)
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
      do i=1,(nx(nb)/2+1)*ny(nb)*nq(nb)
         A(i) = tmp2(int_mb(iq_to_i2(1,nb)+i-1))
      end do
      
      return
      end
#endif





*     ***********************************
*     *					*
*     *	       D3dB_(c,r,t)_Write	*	
*     *					*
*     ***********************************

      subroutine D3dB_c_Write(nb,iunit,A,tmp)
      implicit none 
      integer nb
      integer iunit
      complex*16 A(*)
      complex*16 tmp(*)

      integer nq(NBLOCKS),nx(NBLOCKS),ny(NBLOCKS),nz(NBLOCKS)
      common	/ D3dB / nq,nx,ny,nz

#ifdef MPI
      include 'mpif.h'
      integer mpierr
      integer msgtype
#endif
#ifdef TCGMSG
#include "tcgmsg.fh"
#include "msgtypesf.h"
      integer rcv_len,rcv_proc
#endif
      
*     *** local variables ***
      integer MASTER,taskid
      parameter(MASTER=0)
      integer p_from, p_here,q
      integer k,index
      integer dest,source,status,msglen

      call Parallel_taskid(taskid)

*     **** master node reads from file and distributes ****
      if (taskid.eq.MASTER) then
         do k=1,nz(nb)

            call D3dB_ktoqp(nb,k,q,p_from)

            if (p_from.eq.MASTER) then

               index = (q-1)*(nx(nb)/2+1)*ny(nb) + 1
               call dcopy((nx(nb)+2)*ny(nb),A(index),1,tmp,1)

            else

               msglen  = (nx(nb)/2+1)*ny(nb)
               status  = msglen
               source  = p_from
#ifdef MPI
               msgtype = 33
               call MPI_RECV(tmp,msglen,MPI_DOUBLE_COMPLEX,
     >                       source,msgtype,MPI_COMM_WORLD,
     >                       status,mpierr)
#endif
#ifdef TCGMSG
               call RCV(9+MSGDBL,tmp,mdtob(2*msglen),rcv_len,
     >                       source,rcv_proc,1)
#endif
            end if

            call dwrite(iunit,tmp,(nx(nb)+2)*ny(nb)) 
         end do

*     **** not master node ****
      else
         do k=1,nz(nb)
         
            call D3dB_ktoqp(nb,k,q,p_here)
            if (p_here.eq.taskid) then

               index = (q-1)*(nx(nb)/2+1)*ny(nb) + 1
               call dcopy((nx(nb)+2)*ny(nb),A(index),1,tmp,1)


               msglen  = (nx(nb)/2+1)*ny(nb)
               dest    = MASTER
#ifdef MPI
               msgtype = 33
               call MPI_SEND(tmp,msglen,MPI_DOUBLE_COMPLEX,
     >                       dest,msgtype,MPI_COMM_WORLD,mpierr)
#endif
#ifdef TCGMSG
               call SND(9+MSGDBL,tmp,mdtob(2*msglen),dest,1)
#endif
            end if

         end do
      end if

*     **** wait ****
#ifdef MPI
      call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
#endif

      return
      end

      subroutine D3dB_r_Write(nb,iunit,A,tmp)
      implicit none 
      integer nb
      integer iunit
      real*8     A(*)
      real*8     tmp(*)

      integer nq(NBLOCKS),nx(NBLOCKS),ny(NBLOCKS),nz(NBLOCKS)
      common	/ D3dB / nq,nx,ny,nz

#ifdef MPI
      include 'mpif.h'
      integer mpierr
      integer msgtype
#endif
#ifdef TCGMSG
#include "tcgmsg.fh"
#include "msgtypesf.h"
      integer rcv_len,rcv_proc
#endif
      
*     *** local variables ***
      integer MASTER,taskid
      parameter(MASTER=0)
      integer p_from, p_here,q
      integer j,k,index,index2
      integer dest,source,status,msglen

      call Parallel_taskid(taskid)

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
#ifdef MPI
               msgtype = 33
               call MPI_RECV(tmp,msglen,MPI_DOUBLE_PRECISION,
     >                       source,msgtype,MPI_COMM_WORLD,
     >                       status,mpierr)
#endif
#ifdef TCGMSG
               call RCV(9+MSGDBL,tmp,mdtob(msglen),rcv_len,
     >                       source,rcv_proc,1)
#endif
            end if

            call dwrite(iunit,tmp,(nx(nb))*ny(nb)) 
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
#ifdef MPI
               msgtype = 33
               call MPI_SEND(tmp,msglen,MPI_DOUBLE_PRECISION,
     >                       dest,msgtype,MPI_COMM_WORLD,mpierr)
#endif
#ifdef TCGMSG
               call SND(9+MSGDBL,tmp,mdtob(msglen),dest,1)
#endif
            end if

         end do
      end if

*     **** wait ****
#ifdef MPI
      call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
#endif
      return
      end

*     ***********************************
*     *									*
*     *	       D3dB_r_FormatWrite		*	
*     *									*
*     ***********************************

      subroutine D3dB_r_FormatWrite(nb,iunit,A,tmp)
      implicit none 
      integer nb
      integer iunit
      real*8     A(*)
      real*8     tmp(*)

      integer nq(NBLOCKS),nx(NBLOCKS),ny(NBLOCKS),nz(NBLOCKS)
      common	/ D3dB / nq,nx,ny,nz

#ifdef MPI
      include 'mpif.h'
      integer mpierr
      integer msgtype
#endif
#ifdef TCGMSG
#include "tcgmsg.fh"
#include "msgtypesf.h"
      integer rcv_len,rcv_proc
#endif
      
*     *** local variables ***
      integer MASTER,taskid
      parameter(MASTER=0)
      integer p_from, p_here,q
      integer i,j,k,index,index2
      integer dest,source,status,msglen

      call Parallel_taskid(taskid)

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
#ifdef MPI
               msgtype = 33
               call MPI_RECV(tmp,msglen,MPI_DOUBLE_PRECISION,
     >                       source,msgtype,MPI_COMM_WORLD,
     >                       status,mpierr)
#endif
#ifdef TCGMSG
               call RCV(9+MSGDBL,tmp,mdtob(msglen),rcv_len,
     >                       source,rcv_proc,1)
#endif
            end if

c           **** call dwrite(iunit,tmp,(nx(nb))*ny(nb)) ****
            do j=1,ny(nb)
              write(iunit,'(6E13.5)') (tmp(i+(j-1)*nx(nb)), i=1,nx(nb))
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
#ifdef MPI
               msgtype = 33
               call MPI_SEND(tmp,msglen,MPI_DOUBLE_PRECISION,
     >                       dest,msgtype,MPI_COMM_WORLD,mpierr)
#endif
#ifdef TCGMSG
               call SND(9+MSGDBL,tmp,mdtob(msglen),dest,1)
#endif
            end if

         end do
      end if

*     **** wait ****
#ifdef MPI
      call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
#endif
      return
      end

*     *******************************************
*     *											*
*     *	       D3dB_r_FormatWrite_reverse		*	
*     *											*
*     *******************************************

      subroutine D3dB_r_FormatWrite_reverse(nb,iunit,A,tmp)
      implicit none 
      integer nb
      integer iunit
      real*8     A(*)
      real*8     tmp(*)

      integer nq(NBLOCKS),nx(NBLOCKS),ny(NBLOCKS),nz(NBLOCKS)
      common	/ D3dB / nq,nx,ny,nz

#ifdef MPI
      include 'mpif.h'
      integer mpierr
      integer msgtype
#endif
#ifdef TCGMSG
#include "tcgmsg.fh"
#include "msgtypesf.h"
      integer rcv_len,rcv_proc
#endif
      
*     *** local variables ***
      integer MASTER,taskid
      parameter(MASTER=0)
      integer p_from, p_here,q
      integer i,j,k,index,index2
      integer dest,source,status,msglen

      call Parallel_taskid(taskid)

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
#ifdef MPI
                 msgtype = 33
                 call MPI_RECV(tmp(k),msglen,MPI_DOUBLE_PRECISION,
     >                         source,msgtype,MPI_COMM_WORLD,
     >                         status,mpierr)
#endif
#ifdef TCGMSG
                 call RCV(9+MSGDBL,tmp(k),mdtob(msglen),rcv_len,
     >                         source,rcv_proc,1)
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
                 dest    = MASTER
#ifdef MPI
                 msgtype = 33
                 call MPI_SEND(tmp,msglen,MPI_DOUBLE_PRECISION,
     >                       dest,msgtype,MPI_COMM_WORLD,mpierr)
#endif
#ifdef TCGMSG
                 call SND(9+MSGDBL,tmp,mdtob(msglen),dest,1)
#endif
              end if
            end do

         end do
         end do
      end if

*     **** wait ****
#ifdef MPI
      call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
#endif
      return
      end



*     ***********************************
*     *					*
*     *	         D3dB_cc_dot  	 	*	
*     *					*
*     ***********************************

      subroutine D3dB_cc_dot(nb,A,B,sumall)     
      implicit none
      integer nb 
      complex*16 A(*)
      complex*16 B(*)
      real*8     sumall

      integer nq(NBLOCKS),nx(NBLOCKS),ny(NBLOCKS),nz(NBLOCKS)
      common	/ D3dB / nq,nx,ny,nz

      integer i,j,q,index,np
      real*8  sum

*     ***** timing variables *****
*     1 - fft_time
*     2 - ABC_time
*     3 - L_time
*     4 - vlocal_time
*     5 - vnl_time
*     6 - vh_time
      real*8 tim1,tim2
      real*8 times(50)
      common / timing / times

      call current_second(tim1)

      call Parallel_np(np)

*     **** sum up dot product on this node ****
      sum = 0.0d0

*     ***** kx!=0 plane, so double count *****
      do q=1,nq(nb)
         do j=1,ny(nb)
         do i=2,(nx(nb)/2+1)
            index = (q-1)*(nx(nb)/2+1)*ny(nb) 
     >            + (j-1)*(nx(nb)/2+1) + i
            sum = sum + dble(A(index))  * dble(B(index))
     >                + dimag(A(index)) * dimag(B(index))
         end do
         end do
      end do
      sum = sum*2.0d0

*     ***** kx==0 plane, so single count *****
      do q=1,nq(nb)
         do j=1,ny(nb)
            i=1
            index = (q-1)*(nx(nb)/2+1)*ny(nb) + (j-1)*(nx(nb)/2+1) + 1
            sum = sum + dble(A(index))  * dble(B(index))
     >                + dimag(A(index)) * dimag(B(index))
         end do
      end do
      

*     **** add up sums from other nodes ****
      if (np.gt.1) then
         call D3dB_SumAll(sum)
      end if

      call current_second(tim2)
      times(2) = times(2) + (tim2-tim1)

      sumall = sum
      return
      end

*     ***********************************
*     *					*
*     *	         D3dB_cc_idot  	 	*	
*     *					*
*     ***********************************

      subroutine D3dB_cc_idot(nb,A,B,sumall)     
      implicit none
      integer nb 
      complex*16 A(*)
      complex*16 B(*)
      real*8     sumall

      integer nq(NBLOCKS),nx(NBLOCKS),ny(NBLOCKS),nz(NBLOCKS)
      common	/ D3dB / nq,nx,ny,nz

      integer i,j,q,index,np
      real*8  sum

*     ***** timing variables *****
*     1 - fft_time
*     2 - ABC_time
*     3 - L_time
*     4 - vlocal_time
*     5 - vnl_time
*     6 - vh_time
      real*8 tim1,tim2
      real*8 times(50)
      common / timing / times

      call current_second(tim1)

      call Parallel_np(np)

*     **** sum up dot product on this node ****
      sum = 0.0d0

*     ***** kx!=0 plane, so double count *****
      do q=1,nq(nb)
         do j=1,ny(nb)
         do i=2,(nx(nb)/2+1)
            index = (q-1)*(nx(nb)/2+1)*ny(nb) 
     >            + (j-1)*(nx(nb)/2+1) + i
            sum = sum + dble(A(index))  * dble(B(index))
     >                + dimag(A(index)) * dimag(B(index))
         end do
         end do
      end do
      sum = sum*2.0d0

*     ***** kx==0 plane, so single count *****
      do q=1,nq(nb)
         do j=1,ny(nb)
            i=1
            index = (q-1)*(nx(nb)/2+1)*ny(nb) + (j-1)*(nx(nb)/2+1) + 1
            sum = sum + dble(A(index))  * dble(B(index))
     >                + dimag(A(index)) * dimag(B(index))
         end do
      end do
      

*     **** do not add up sums from other nodes ****
       
      call current_second(tim2)
      times(2) = times(2) + (tim2-tim1)

      sumall = sum
      return
      end

*     ***********************************
*     *					*
*     *	         D3dB_tt_dot  	 	*	
*     *					*
*     ***********************************

      subroutine D3dB_tt_dot(nb,A,B,sumall)     
      implicit none 
      integer nb
      real*8 A(*)
      real*8 B(*)
      real*8 sumall

      integer nq(NBLOCKS),nx(NBLOCKS),ny(NBLOCKS),nz(NBLOCKS)
      common	/ D3dB / nq,nx,ny,nz

      integer i,j,q,index,np,nxh
      real*8  sum

      nxh=nx(nb)/2
      call Parallel_np(np)

*     **** sum up dot product on this node ****
      sum = 0.0d0

*     ***** k!=0 plane, so double count *****
      do q=1,nq(nb)
         do j=1,ny(nb)
         do i=2,(nxh+1)
            index = (q-1)*(nxh+1)*ny(nb) + (j-1)*(nxh+1) + i
            sum = sum + A(index)*B(index)
         end do
         end do
      end do
      sum = sum*2.0d0

*     **** kx==0 plane, so single count *****
      do q=1,nq(nb)
         do j=1,ny(nb)
            i=1
            index = (q-1)*(nxh+1)*ny(nb) + (j-1)*(nxh+1) + 1
            sum = sum + A(index)*B(index)
         end do
      end do
      

*     **** add up sums from other nodes ****
      if (np.gt.1) then
         call D3dB_SumAll(sum)
      end if

      sumall = sum
      return
      end


*     ***********************************
*     *					*
*     *	         D3dB_tt_idot  	 	*	
*     *					*
*     ***********************************

      subroutine D3dB_tt_idot(nb,A,B,sumall)     
      implicit none 
      integer nb
      real*8 A(*)
      real*8 B(*)
      real*8 sumall

      integer nq(NBLOCKS),nx(NBLOCKS),ny(NBLOCKS),nz(NBLOCKS)
      common	/ D3dB / nq,nx,ny,nz

      integer i,j,q,index,np,nxh
      real*8  sum

*     ***** timing variables *****
*     1 - fft_time
*     2 - dot_time
*     3 - ABC_time
*     4 - L_time
*     5 - vlocal_time
*     6 - vnl_time
*     7 - vh_time
*     8 - strfac_time
*     9 - maskC_time
      real*8 tim1,tim2
      real*8 times(50)
      common / timing / times


      call current_second(tim1)

      nxh=nx(nb)/2
      call Parallel_np(np)

*     **** sum up dot product on this node ****
      sum = 0.0d0

*     ***** k!=0 plane, so double count *****
      do q=1,nq(nb)
         do j=1,ny(nb)
         do i=2,(nxh+1)
            index = (q-1)*(nxh+1)*ny(nb) + (j-1)*(nxh+1) + i
            sum = sum + A(index)*B(index)
         end do
         end do
      end do
      sum = sum*2.0d0

*     **** kx==0 plane, so single count *****
      do q=1,nq(nb)
         do j=1,ny(nb)
            i=1
            index = (q-1)*(nxh+1)*ny(nb) + (j-1)*(nxh+1) + 1
            sum = sum + A(index)*B(index)
         end do
      end do
      

*     **** !!!! do not add up sums from other nodes ****

      call current_second(tim2)
      times(2) = times(2) + (tim2-tim1)

      sumall = sum
      return
      end






*     ***********************************
*     *					*
*     *	         D3dB_rr_dot  	 	*	
*     *					*
*     ***********************************

      subroutine D3dB_rr_dot(nb,A,B,sumall)     
      implicit none 
      integer nb
      real*8  A(*)
      real*8  B(*)
      real*8  sumall

      integer nq(NBLOCKS),nx(NBLOCKS),ny(NBLOCKS),nz(NBLOCKS)
      common	/ D3dB / nq,nx,ny,nz

      integer i,np
      real*8  sum

      call Parallel_np(np)

*     **** sum up dot product on this node ****
      sum = 0.0d0
      do i=1,(nx(nb)+2)*ny(nb)*nq(nb)
         sum = sum + A(i)*B(i)
      end do

*     **** add up sums from other nodes ****
      if (np.gt.1) then
         call D3dB_SumAll(sum)
      end if

      sumall = sum
      return
      end

*     ***********************************
*     *					*
*     *	         D3dB_rr_idot  	 	*	
*     *					*
*     ***********************************

      subroutine D3dB_rr_idot(nb,A,B,sumall)     
      implicit none 
      integer nb
      real*8  A(*)
      real*8  B(*)
      real*8  sumall

      integer nq(NBLOCKS),nx(NBLOCKS),ny(NBLOCKS),nz(NBLOCKS)
      common	/ D3dB / nq,nx,ny,nz

      integer i,np
      real*8  sum

      call Parallel_np(np)

*     **** sum up dot product on this node ****
      sum = 0.0d0
      do i=1,(nx(nb)+2)*ny(nb)*nq(nb)
         sum = sum + A(i)*B(i)
      end do

*     **** add up sums from other nodes ****
*     if (np.gt.1) then
*        call D3dB_SumAll(sum)
*     end if

      sumall = sum
      return
      end



*     ***********************************
*     *					*
*     *	         D3dB_cc_Mul  	 	*	
*     *					*
*     ***********************************

      subroutine D3dB_cc_Mul(nb,A,B,C)     
      implicit none 
      integer    nb
      complex*16 A(*)
      complex*16 B(*)
      complex*16 C(*)

      integer nq(NBLOCKS),nx(NBLOCKS),ny(NBLOCKS),nz(NBLOCKS)
      common	/ D3dB / nq,nx,ny,nz

      integer i

      do i=1,(nx(nb)/2+1)*ny(nb)*nq(nb)
            C(i) = dconjg(A(i)) * B(i)
         end do

      return
      end

*     ***********************************
*     *					*
*     *	         D3dB_lc_Mask  	 	*	
*     *					*
*     ***********************************

      subroutine D3dB_lc_Mask(nb,masker,A)
      implicit none 
      integer    nb
      logical    masker(*)
      complex*16 A(*)

      integer nq(NBLOCKS),nx(NBLOCKS),ny(NBLOCKS),nz(NBLOCKS)
      common	/ D3dB / nq,nx,ny,nz

      integer i

      do i=1,(nx(nb)/2+1)*ny(nb)*nq(nb)
         if (masker(i)) A(i) = dcmplx(0.0d0,0.0d0)
      end do
      return
      end

*     ***********************************
*     *					*
*     *	         D3dB_lr_Mask  	 	*	
*     *					*
*     ***********************************

      subroutine D3dB_lr_Mask(nb,masker,A)
      implicit none 
      integer   nb
      logical   masker(*)
      real*8    A(*)

      integer nq(NBLOCKS),nx(NBLOCKS),ny(NBLOCKS),nz(NBLOCKS)
      common	/ D3dB / nq,nx,ny,nz

      integer i

      do i=1,(nx(nb)/2+1)*ny(nb)*nq(nb)
         if (masker(i)) A(i) = 0.0d0
      end do
      return
      end


*     ***********************************
*     *					*
*     *	         D3dB_tc_Mul  	 	*	
*     *					*
*     ***********************************

      subroutine D3dB_tc_Mul(nb,A,B,C)     
      implicit none 
      integer    nb
      real*8     A(*)
      complex*16 B(*)
      complex*16 C(*)

      integer nq(NBLOCKS),nx(NBLOCKS),ny(NBLOCKS),nz(NBLOCKS)
      common	/ D3dB / nq,nx,ny,nz

      integer i

      do i=1,(nx(nb)/2+1)*ny(nb)*nq(nb)
            C(i) = A(i) * B(i)
      end do

      return
      end

*     ***********************************
*     *					*
*     *	         D3dB_rr_Mul  	 	*	
*     *					*
*     ***********************************

      subroutine D3dB_rr_Mul(nb,A,B,C)     
      implicit none 
      integer nb
      real*8 A(*)
      real*8 B(*)
      real*8 C(*)

      integer nq(NBLOCKS),nx(NBLOCKS),ny(NBLOCKS),nz(NBLOCKS)
      common	/ D3dB / nq,nx,ny,nz

      integer i

      do i=1,(nx(nb)+2)*ny(nb)*nq(nb)
            C(i) = A(i) * B(i)
         end do

      return
      end



*     ***********************************
*     *					*
*     *	         D3dB_cc_Sum  	 	*	
*     *					*
*     ***********************************

      subroutine D3dB_cc_Sum(nb,A,B,C)     
      implicit none 
      integer    nb
      complex*16 A(*)
      complex*16 B(*)
      complex*16 C(*)

      integer nq(NBLOCKS),nx(NBLOCKS),ny(NBLOCKS),nz(NBLOCKS)
      common	/ D3dB / nq,nx,ny,nz

      integer i

      do i=1,(nx(nb)/2+1)*ny(nb)*nq(nb)
         C(i) = A(i) + B(i)
      end do

      return
      end


*     ***********************************
*     *					*
*     *	         D3dB_rr_Sum  	 	*	
*     *					*
*     ***********************************

      subroutine D3dB_rr_Sum(nb,A,B,C)     
      implicit none 
      integer nb
      real*8  A(*)
      real*8  B(*)
      real*8  C(*)

      integer nq(NBLOCKS),nx(NBLOCKS),ny(NBLOCKS),nz(NBLOCKS)
      common	/ D3dB / nq,nx,ny,nz

      integer i

      do i=1,(nx(nb)+2)*ny(nb)*nq(nb)
         C(i) = A(i) + B(i)
      end do

      return
      end

*     ***********************************
*     *					*
*     *	         D3dB_cc_Sub  	 	*	
*     *					*
*     ***********************************

      subroutine D3dB_cc_Sub(nb,A,B,C)     
      implicit none 
      integer    nb
      complex*16 A(*)
      complex*16 B(*)
      complex*16 C(*)

      integer nq(NBLOCKS),nx(NBLOCKS),ny(NBLOCKS),nz(NBLOCKS)
      common	/ D3dB / nq,nx,ny,nz

      integer i

      do i=1,(nx(nb)/2+1)*ny(nb)*nq(nb)
         C(i) = A(i) - B(i)
      end do

      return
      end


*     ***********************************
*     *									*
*     *	         D3dB_rr_Sub  	 		*	
*     *									*
*     ***********************************

      subroutine D3dB_rr_Sub(nb,A,B,C)     
      implicit none 
      integer nb
      real*8  A(*)
      real*8  B(*)
      real*8  C(*)

      integer nq(NBLOCKS),nx(NBLOCKS),ny(NBLOCKS),nz(NBLOCKS)
      common	/ D3dB / nq,nx,ny,nz

      integer i

      do i=1,(nx(nb)+2)*ny(nb)*nq(nb)
         C(i) = A(i) - B(i)
      end do

      return
      end




*     ***********************************
*     *					*
*     *	         D3dB_cc_zaxpy 	 	*	
*     *					*
*     ***********************************

      subroutine D3dB_cc_zaxpy(nb,alpha,A,B)     
      implicit none 
      integer    nb
      complex*16 alpha
      complex*16 A(*)
      complex*16 B(*)

      integer nq(NBLOCKS),nx(NBLOCKS),ny(NBLOCKS),nz(NBLOCKS)
      common	/ D3dB / nq,nx,ny,nz

      integer i

      do i=1,(nx(nb)/2+1)*ny(nb)*nq(nb)
         B(i) = B(i) + alpha*A(i)
      end do

      return
      end



*     ***********************************
*     *					*
*     *	         D3dB_cc_daxpy 	 	*	
*     *					*
*     ***********************************

      subroutine D3dB_cc_daxpy(nb,alpha,A,B)     
      implicit none 
      integer    nb
      real*8     alpha
      complex*16 A(*)
      complex*16 B(*)

      integer nq(NBLOCKS),nx(NBLOCKS),ny(NBLOCKS),nz(NBLOCKS)
      common	/ D3dB / nq,nx,ny,nz

      integer i

      do i=1,(nx(nb)/2+1)*ny(nb)*nq(nb)
         B(i) = B(i) + alpha*A(i)
      end do

      return
      end

*     ***********************************
*     *					*
*     *	         D3dB_rr_daxpy 	 	*	
*     *					*
*     ***********************************

      subroutine D3dB_rr_daxpy(nb,alpha,A,B)     
      implicit none 
      integer nb
      real*8  alpha
      real*8  A(*)
      real*8  B(*)

      integer nq(NBLOCKS),nx(NBLOCKS),ny(NBLOCKS),nz(NBLOCKS)
      common	/ D3dB / nq,nx,ny,nz

      integer i

      do i=1,(nx(nb)+2)*ny(nb)*nq(nb)
         B(i) = B(i) + alpha* B(i)
      end do

      return
      end

*     ***********************************
*     *                                 *
*     *          D3dB_rr_Divide         *
*     *                                 *
*     ***********************************

      subroutine D3dB_rr_Divide(nb,A,B,C)
      implicit none
      integer nb
      real*8 A(*)
      real*8 B(*)
      real*8 C(*)

      integer nq(NBLOCKS),nx(NBLOCKS),ny(NBLOCKS),nz(NBLOCKS)
      common	/ D3dB / nq,nx,ny,nz

      real*8 eta
      parameter (eta=1.0d-9)

      integer i,j,q,index

      do q=1,nq(nb)
      do j=1,ny(nb)
      do i=1,nx(nb)
         index = i + (j-1)*(nx(nb)+2) + (q-1)*(nx(nb)+2)*ny(nb)
         if (dabs(B(index)) .le. eta) then
           C(index) = 0.0d0
         else
           C(index) = A(index) / B(index)
         end if
      end do
      end do
      end do

      return
      end

*     ***********************************
*     *                                 *
*     *          D3dB_rr_Minus          *
*     *                                 *
*     ***********************************
      subroutine D3dB_rr_Minus(nb,A,B,C)
      implicit none
      integer nb
      real*8  A(*)
      real*8  B(*)
      real*8  C(*)

      integer nq(NBLOCKS),nx(NBLOCKS),ny(NBLOCKS),nz(NBLOCKS)
      common	/ D3dB / nq,nx,ny,nz

      integer i

      do i=1,(nx(nb)+2)*ny(nb)*nq(nb)
         C(i) = A(i) - B(i)
      end do

      return
      end




*     ***********************************
*     *					*
*     *          D3dB_r_Zero_Ends	*
*     *					*
*     ***********************************

      subroutine D3dB_r_Zero_Ends(nb,A)
      integer nb
      real*8 A(*) 

      integer nq(NBLOCKS),nx(NBLOCKS),ny(NBLOCKS),nz(NBLOCKS)
      common	/ D3dB / nq,nx,ny,nz

      integer j,q,index
      do q=1,nq(nb)
         do j=1,ny(nb)
            index = (nx(nb)+1) + (j-1)*(nx(nb)+2) 
     >                         + (q-1)*(nx(nb)+2)*(ny(nb))
            A(index)   = 0.0d0
            A(index+1) = 0.0d0
         end do
      end do

      return
      end


*     ***********************************
*     *					*
*     *	         D3dB_r_dsum  	 	*	
*     *					*
*     ***********************************

      subroutine D3dB_r_dsum(nb,A,sumall)     
      implicit none
      integer nb 
      real*8  A(*)
      real*8  sumall

      integer nq(NBLOCKS),nx(NBLOCKS),ny(NBLOCKS),nz(NBLOCKS)
      common	/ D3dB / nq,nx,ny,nz

      integer i,np
      real*8 sum

      call Parallel_np(np)

*     **** sum up dot product on this node ****
      sum = 0.0d0
      do i=1,(nx(nb)+2)*ny(nb)*nq(nb)
         sum = sum + A(i)
      end do

*     **** add up sums from other nodes ****
      if (np.gt.1) then
        call D3dB_SumAll(sum)
      end if

      sumall = sum

      return
      end

*     ***********************************
*     *					*
*     *	         D3dB_t_dsum  	 	*	
*     *					*
*     ***********************************

      subroutine D3dB_t_dsum(nb,A,sumall)     
      implicit none 
      integer nb
      real*8  A(*)
      real*8  sumall

      integer nq(NBLOCKS),nx(NBLOCKS),ny(NBLOCKS),nz(NBLOCKS)
      common	/ D3dB / nq,nx,ny,nz

      integer i,j,q,np,nxh,index
      real*8 sum

      nxh = nx(nb)/2
      call Parallel_np(np)

*     **** sum up dot product on this node ****
      sum = 0.0d0

*     ***** k!=0 plane so double count *****
      do q=1,nq(nb)
      do j=1,ny(nb)
         do i=2,(nxh+1)
            index = (q-1)*(nxh+1)*ny(nb) + (j-1)*(nxh+1) + i
             sum = sum + A(index)
         end do
      end do
      end do
      sum = sum*2.0d0

*     ***** k==0 plane, so single count *****
      do q=1,nq(nb)
      do j=1,ny(nb)
            index = (q-1)*(nxh+1)*ny(nb) + (j-1)*(nxh+1) + 1
             sum = sum + A(index)
      end do
      end do

*     **** add up sums from other nodes ****
      if (np.gt.1) then
        call D3dB_SumAll(sum)
      end if

      sumall = sum

      return
      end


*     ***********************************
*     *					*
*     *	     D3dB_cc_Vector_dot 	*	
*     *					*
*     ***********************************

      subroutine D3dB_cc_Vector_dot(nb,nfft3d,nn,ne,A,B,sumall)     
      implicit none 
      integer    nb
      integer    nfft3d,nn,ne
      complex*16 A(*)
      complex*16 B(*)
      real*8     sumall(nn,nn)

      integer nq(NBLOCKS),nx(NBLOCKS),ny(NBLOCKS),nz(NBLOCKS)
      common	/ D3dB / nq,nx,ny,nz

      integer i,j,q,index,np
      integer index1,index2
      integer n,m,shift1,shift2
      real*8  sum

*     ***** timing variables *****
*     1 - fft_time
*     2 - dot_time
      real*8 tim1,tim2
      real*8 times(50)
      common / timing / times

      call current_second(tim1)

      call Parallel_np(np)

*     **** sum up dot product on this node ****
      do n=1,ne
      do m=n,ne

        shift1 = (n-1)*nfft3d
        shift2 = (m-1)*nfft3d
        sum    = 0.0d0
*       ***** kx!=0 plane, so double count *****
        do q=1,nq(nb)
           do j=1,ny(nb)
           do i=2,(nx(nb)/2+1)
              index = (q-1)*(nx(nb)/2+1)*ny(nb) 
     >              + (j-1)*(nx(nb)/2+1) + i
              index1 = index+shift1
              index2 = index+shift2
              sum = sum + dble(A(index1))  * dble(B(index2))
     >                  + dimag(A(index1)) * dimag(B(index2))
           end do
           end do
        end do
        sum = sum*2.0d0

*       ***** kx==0 plane, so single count *****
        do q=1,nq(nb)
           do j=1,ny(nb)
              i=1
              index = (q-1)*(nx(nb)/2+1)*ny(nb) 
     >              + (j-1)*(nx(nb)/2+1) + 1
              index1 = index+shift1
              index2 = index+shift2
              sum = sum + dble(A(index1))  * dble(B(index2))
     >                  + dimag(A(index1)) * dimag(B(index2))
           end do
        end do

         sumall(n,m) = sum
         sumall(m,n) = sum
      end do
      end do

      call current_second(tim2)
      times(2) = times(2) + (tim2-tim1)

*     **** add up sums from other nodes ****
      if (np.gt.1) then
         call D3dB_Vector_SumAll(nn*ne,sumall)
      end if


      return
      end



*     ***********************************
*     *					*
*     *	     D3dB_cc_Vector_ndot 	*	
*     *					*
*     ***********************************

      subroutine D3dB_cc_Vector_ndot(nb,nfft3d,ne,A,B,sumall)     
      implicit none 
      integer    nb
      integer    nfft3d,ne
      complex*16 A(*)
      complex*16 B(*)
      real*8     sumall(ne)

      integer nq(NBLOCKS),nx(NBLOCKS),ny(NBLOCKS),nz(NBLOCKS)
      common	/ D3dB / nq,nx,ny,nz

      integer i,j,q,index,np
      integer index1,index2
      integer n,shift1,shift2
      real*8  sum

*     ***** timing variables *****
*     1 - fft_time
*     2 - ABC_time
*     3 - L_time
*     4 - vlocal_time
*     5 - vnl_time
*     6 - vh_time
      real*8 tim1,tim2
      real*8 times(50)
      common / timing / times

      call current_second(tim1)

      call Parallel_np(np)

*     **** sum up dot product on this node ****
      do n=1,ne

        shift1 = (n-1)*nfft3d
        shift2 = 0
        sum    = 0.0d0
*       ***** kx!=0 plane, so double count *****
        do q=1,nq(nb)
           do j=1,ny(nb)
           do i=2,(nx(nb)/2+1)
              index = (q-1)*(nx(nb)/2+1)*ny(nb) 
     >              + (j-1)*(nx(nb)/2+1) + i
              index1 = index+shift1
              index2 = index+shift2
              sum = sum + dble(A(index1))  * dble(B(index2))
     >                  + dimag(A(index1)) * dimag(B(index2))
           end do
           end do
        end do
        sum = sum*2.0d0

*       ***** kx==0 plane, so single count *****
        do q=1,nq(nb)
           do j=1,ny(nb)
              i=1
              index = (q-1)*(nx(nb)/2+1)*ny(nb) 
     >              + (j-1)*(nx(nb)/2+1) + 1
              index1 = index+shift1
              index2 = index+shift2
              sum = sum + dble(A(index1))  * dble(B(index2))
     >                  + dimag(A(index1)) * dimag(B(index2))
           end do
        end do

         sumall(n) = sum
      end do

      call current_second(tim2)
      times(2) = times(2) + (tim2-tim1)

*     **** add up sums from other nodes ****
      if (np.gt.1) then
         call D3dB_Vector_SumAll(ne,sumall)
      end if


      return
      end





*     ***********************************
*     *					*
*     *	        D3dB_Vector_SumAll	*	
*     *					*
*     ***********************************

      subroutine D3dB_Vector_SumAll(n,sum)
c     implicit none
      integer n
      real*8  sum(*)

#include "mafdecls.fh"


#ifdef MPI
      include 'mpif.h'
      integer mpierr
#endif
#ifdef TCGMSG
#include "tcgmsg.fh"
#include "msgtypesf.h"
#endif

      integer MASTER
      parameter (MASTER=0)
      integer msglen
      logical value

*     **** temporary workspace ****
      integer sumall(2)

*     ***** timing variables *****
*     1 - fft_time
*     2 - dot_time
      real*8 tim1,tim2
      real*8 times(50)
      common / timing / times

      call current_second(tim1)

*     ***** allocate temporary space ****
      value = MA_push_get(mt_dbl,n,'sumall',sumall(2),sumall(1))
      if (.not. value) call errquit('out of stack memory',0)

      msglen = n

#ifdef MPI
      call dcopy(n,0.0d0,0,dbl_mb(sumall(1)),1)
      call MPI_REDUCE(sum,dbl_mb(sumall(1)),msglen,
     >                     MPI_DOUBLE_PRECISION,
     >                     MPI_SUM,MASTER,MPI_COMM_WORLD,mpierr)
      call MPI_BCAST(dbl_mb(sumall(1)),msglen,
     >                     MPI_DOUBLE_PRECISION,
     >                     MASTER,MPI_COMM_WORLD,mpierr)
#endif
#ifdef TCGMSG
      call dcopy(n,sum,1,dbl_mb(sumall(1)),1)
      call GA_DGOP(9+MSGDBL,dbl_mb(sumall(1)),n,'+')
#endif

      call dcopy(n,dbl_mb(sumall(1)),1,sum,1)
      value = MA_pop_stack(sumall(2))

      call current_second(tim2)
      times(2) = times(2) + (tim2-tim1)
      return
      end


*     ***********************************
*     *					*
*     *	        D3dB_Vector_ISumAll	*	
*     *					*
*     ***********************************

      subroutine D3dB_Vector_ISumAll(n,sum)
c     implicit none
      integer n
      integer  sum(*)

#include "mafdecls.fh"


#ifdef MPI
      include 'mpif.h'
      integer mpierr
#endif
#ifdef TCGMSG
#include "tcgmsg.fh"
#include "msgtypesf.h"
#endif

      integer MASTER
      parameter (MASTER=0)
      integer msglen
      logical value

*     **** temporary workspace ****
      integer sumall(2)

*     ***** timing variables *****
*     1 - fft_time
*     2 - dot_time
      real*8 tim1,tim2
      real*8 times(50)
      common / timing / times

      call current_second(tim1)

*     ***** allocate temporary space ****
      value = MA_push_get(mt_int,n,'sumall',sumall(2),sumall(1))
      if (.not. value) call errquit('out of stack memory',0)

      msglen = n

#ifdef MPI
      call icopy(n,0.0d0,0,int_mb(sumall(1)),1)
      call MPI_REDUCE(sum,int_mb(sumall(1)),msglen,
     >                     MPI_INTEGER,
     >                     MPI_SUM,MASTER,MPI_COMM_WORLD,mpierr)
      call MPI_BCAST(int_mb(sumall(1)),msglen,
     >                     MPI_INTEGER,
     >                     MASTER,MPI_COMM_WORLD,mpierr)
#endif
#ifdef TCGMSG
      call icopy(n,sum,1,int_mb(sumall(1)),1)
      call GA_IGOP(9+MSGINT,int_mb(sumall(1)),n,'+')
#endif

      call icopy(n,int_mb(sumall(1)),1,sum,1)
      value = MA_pop_stack(sumall(2))

      call current_second(tim2)
      times(2) = times(2) + (tim2-tim1)
      return
      end

c *** icopy define in src/util directory!!!
c      subroutine icopy(n,a,inca,b,incb)
c      integer n
c      integer a(*),inca
c      integer b(*),incb
c
c      integer i,shifta,shiftb
c
c      shifta = 1
c      shiftb = 1
c      do i=1,n
c        b(shiftb)=a(shifta)
c        shifta=shifta+inca
c        shiftb=shiftb+incb
c      end do
c
c      return
c      end


*     ***********************************
*     *                                 *
*     *          D3dB_ic_Mul            *
*     *                                 *
*     ***********************************

      subroutine D3dB_ic_Mul(nb,A,B,C)
      implicit none
      integer    nb
      real*8     A(*)
      complex*16 B(*)
      complex*16 C(*)

      integer nq(NBLOCKS),nx(NBLOCKS),ny(NBLOCKS),nz(NBLOCKS)
      common	/ D3dB / nq,nx,ny,nz

      integer i

      do i=1,(nx(nb)/2+1)*ny(nb)*nq(nb)
            C(i) = dcmplx(0.0d0,A(i)) * B(i)
      end do

      return
      end

*     ***********************************
*     *                                 *
*     *          D3dB_r_Expand          *
*     *                                 *
*     ***********************************

      subroutine D3dB_r_Expand(nb1,A,nb2,B)
      implicit none
      integer nb1
      real*8  A(*)
      integer nb2
      real*8  B(*)

      integer nq(NBLOCKS),nx(NBLOCKS),ny(NBLOCKS),nz(NBLOCKS)
      common	/ D3dB / nq,nx,ny,nz

      integer i,j,q,index1,index2

      call dcopy(nq(nb2)*ny(nb2)*(nx(nb2)+2),0.0d0,0,B,1)
      do q=1,nq(nb1)
      do j=1,ny(nb1)
      do i=1,nx(nb1)
         index1 = i + (j-1)*(nx(nb1)+2) + (q-1)*(nx(nb1)+2)*ny(nb1)
         index2 = i + (j-1)*(nx(nb2)+2) + (q-1)*(nx(nb2)+2)*ny(nb2)
         B(index2) = A(index1)
      end do
      end do
      end do

      return
      end

*     ***********************************
*     *                                 *
*     *          D3dB_r_Contract        *
*     *                                 *
*     ***********************************

      subroutine D3dB_r_Contract(nb2,B,nb1,A)
      implicit none
      integer nb2
      real*8  B(*)
      integer nb1
      real*8  A(*)

      integer nq(NBLOCKS),nx(NBLOCKS),ny(NBLOCKS),nz(NBLOCKS)
      common	/ D3dB / nq,nx,ny,nz

      integer i,j,q,index1,index2

      call dcopy(nq(nb1)*ny(nb1)*(nx(nb1)+2),0.0d0,0,A,1)
      do q=1,nq(nb1)
      do j=1,ny(nb1)
      do i=1,nx(nb1)
         index1 = i + (j-1)*(nx(nb1)+2) + (q-1)*(nx(nb1)+2)*ny(nb1)
         index2 = i + (j-1)*(nx(nb2)+2) + (q-1)*(nx(nb2)+2)*ny(nb2)
         A(index1) = B(index2)
      end do
      end do
      end do

      return
      end


#ifdef TCGMSG
*     ***********************************
*     *					*
*     *	   D3dB_c_timereverse_init	*
*     *					*
*     ***********************************

      subroutine D3dB_c_timereverse_init(nb)
      implicit none
      integer nb

#include "mafdecls.fh"

      integer nq(NBLOCKS),nx(NBLOCKS),ny(NBLOCKS),nz(NBLOCKS)
      common	/ D3dB / nq,nx,ny,nz

c     integer iq_to_i1(2*NFFT2*NSLABS)
c     integer iq_to_i2(2*NFFT2*NSLABS)
c     integer i1_start(NFFT3+1)
c     integer i2_start(NFFT3+1)
      integer iq_to_i1(2,NBLOCKS)
      integer iq_to_i2(2,NBLOCKS)
      integer i1_start(2,NBLOCKS)
      integer i2_start(2,NBLOCKS)
      common / timereverse_blk / iq_to_i1,iq_to_i2,i1_start,i2_start

*     **** local variables ****
      integer proc_to,proc_from
      integer pto,qto,np,taskid
      integer pfrom,qfrom
      integer phere,qhere
      integer index1,index2,itmp
      integer it
      integer i2,j2,k2
      integer i3,j3,k3
      integer nyh,nzh
      logical value


*     **** allocate timereverse_blk common block ****
      value = MA_alloc_get(mt_int,(2*ny(nb)*nq(nb)),
     >                     'iq_to_i1',iq_to_i1(2,nb),iq_to_i1(1,nb))
      value = value.and.
     >        MA_alloc_get(mt_int,(2*ny(nb)*nq(nb)),
     >                     'iq_to_i2',iq_to_i2(2,nb),iq_to_i2(1,nb))

      value = value.and.
     >        MA_alloc_get(mt_int,(nz(nb)+1),
     >                     'i1_start',i1_start(2,nb),i1_start(1,nb))
      value = value.and.
     >        MA_alloc_get(mt_int,(nz(nb)+1),
     >                     'i2_start',i2_start(2,nb),i2_start(1,nb))
      if (.not.value) call errquit('out of heap memory',0)

      call Parallel_taskid(taskid)
      call Parallel_np(np)

  
      nyh = ny(nb)/2
      nzh = nz(nb)/2

      index1 = 1 
      index2 = 1
      do it=0,np-1
         proc_to   = mod(taskid+it,np)
         proc_from = mod(taskid-it+np,np)
c        i1_start(it+1) = index1
c        i2_start(it+1) = index2
         int_mb(i1_start(1,nb)+it) = index1
         int_mb(i2_start(1,nb)+it) = index2

*        *********************
*        **** K=(0,0,k3)  ****
*        *********************
         do k3=1,(nzh-1)
            i3 =  k3
            j3 = -k3
            if (i3.lt.0) i3 = i3 + nz(nb)
            if (j3.lt.0) j3 = j3 + nz(nb)
            i2 = 1
            i3 = i3 + 1
            j2 = 1
            j3 = j3 + 1

*           **** packing scheme ****
            call D3dB_ktoqp(nb,i3,qhere,phere)
            call D3dB_ktoqp(nb,j3,qto,pto)
            if ((phere.eq.taskid).and.(pto.eq.proc_to)) then
               itmp = 1 + (i2-1)*(nx(nb)/2+1) 
     >                  + (qhere-1)*(nx(nb)/2+1)*ny(nb)
c               iq_to_i1(index1) = itmp
               int_mb(iq_to_i1(1,nb)+index1-1)=itmp
               index1 = index1 + 1
            end if

*           **** unpacking scheme ****
            call D3dB_ktoqp(nb,j3,qhere,phere)
            call D3dB_ktoqp(nb,i3,qfrom,pfrom)
            if ((phere.eq.taskid).and.(pfrom.eq.proc_from)) then
               itmp = 1 + (j2-1)*(nx(nb)/2+1) 
     >                  + (qhere-1)*(nx(nb)/2+1)*ny(nb)
c               iq_to_i2(index2) = itmp
               int_mb(iq_to_i2(1,nb)+index2-1) = itmp
               index2 = index2 + 1
            end if

         end do

*        *********************
*        **** k=(0,k2,k3) ****
*        *********************
         do k3=(-nzh+1),(nzh-1)
         do k2=1,(nyh-1)
            i2 =  k2
            i3 =  k3
            j2 = -k2
            j3 = -k3
            if (i2.lt.0) i2 = i2 + ny(nb)
            if (i3.lt.0) i3 = i3 + nz(nb)
            if (j2.lt.0) j2 = j2 + ny(nb)
            if (j3.lt.0) j3 = j3 + nz(nb)
            i2 = i2 + 1
            i3 = i3 + 1
            j2 = j2 + 1
            j3 = j3 + 1

*           **** packing scheme ****
            call D3dB_ktoqp(nb,i3,qhere,phere)
            call D3dB_ktoqp(nb,j3,qto,pto)
            if ((phere.eq.taskid).and.(pto.eq.proc_to)) then
               itmp = 1 + (i2-1)*(nx(nb)/2+1) 
     >                  + (qhere-1)*(nx(nb)/2+1)*ny(nb)
c               iq_to_i1(index1) = itmp
               int_mb(iq_to_i1(1,nb)+index1-1) = itmp
               index1 = index1 + 1
            end if

*           **** unpacking scheme ****
            call D3dB_ktoqp(nb,j3,qhere,phere)
            call D3dB_ktoqp(nb,i3,qfrom,pfrom)
            if ((phere.eq.taskid).and.(pfrom.eq.proc_from)) then
               itmp = 1 + (j2-1)*(nx(nb)/2+1) 
     >                  + (qhere-1)*(nx(nb)/2+1)*ny(nb)
c               iq_to_i2(index2) = itmp
               int_mb(iq_to_i2(1,nb)+index2-1)=itmp
               index2 = index2 + 1
            end if
         end do
         end do

      end do
c     i1_start(np+1) = index1
c     i2_start(np+1) = index2
      int_mb(i1_start(1,nb)+np) = index1
      int_mb(i2_start(1,nb)+np) = index2

      return
      end


             
*     ***********************************
*     *					*
*     *	   D3dB_c_timereverse 		*
*     *					*
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

      integer nq(NBLOCKS),nx(NBLOCKS),ny(NBLOCKS),nz(NBLOCKS)
      common	/ D3dB / nq,nx,ny,nz


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
*     *	     D3dB_timereverse_end 	*	
*     *					*
*     ***********************************
      subroutine D3dB_timereverse_end(nb)
      implicit none
      integer nb

#include "mafdecls.fh"

c     integer iq_to_i1((NFFT1/2+1)*NFFT2*NSLABS)
c     integer iq_to_i2((NFFT1/2+1)*NFFT2*NSLABS)
c     integer i1_start(NFFT3+1)
c     integer i2_start(NFFT3+1)
      integer iq_to_i1(2,NBLOCKS)
      integer iq_to_i2(2,NBLOCKS)
      integer i1_start(2,NBLOCKS)
      integer i2_start(2,NBLOCKS)
      common / timereverse_blk / iq_to_i1,iq_to_i2,i1_start,i2_start

      logical value

      value = MA_free_heap(i1_start(2,nb))
      value = value.and.
     >        MA_free_heap(i2_start(2,nb))
      value = value.and.
     >        MA_free_heap(iq_to_i1(2,nb))
      value = value.and.
     >        MA_free_heap(iq_to_i2(2,nb))
      if (.not.value) call errquit('error freeing heap',0)
      return
      end



#endif



















