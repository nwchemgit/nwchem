ccccccccccc  #define TCGMSG
#define NBLOCKS 3

*
* $Id$
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
     
#include "bafdecls.fh"
#include "errquit.fh"
#include "D3dB.fh"


      integer k,q,p
*     integer kn
      integer taskid,np
      logical value
      integer  tid,Parallel_threadid
      external Parallel_threadid
      
      call Parallel2d_np_i(np)
      call Parallel2d_taskid_i(taskid)



*     **************************
*     ****** Slab mapping ******
*     **************************
      if (mapping.eq.1) then

*     **** allocate q_map,p_map,k_map
      value = BA_alloc_get(mt_int,nz(nb),'q_map',q_map(2,nb),
     >                                       q_map(1,nb))
      value = value.and.BA_alloc_get(mt_int,nz(nb),'p_map',p_map(2,nb),
     >                                       p_map(1,nb))
      value = value.and.BA_alloc_get(mt_int,nz(nb),'k_map',k_map(2,nb),
     >                                       k_map(1,nb))
      if (.not. value)
     > call errquit('Mapping_init:out of heap memory',0, MA_ERR)



*     ****** cyclic ******
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


c      if (nb.eq.1) then
c       p = 0
c       q = 1
c       do k=1,nz(nb)
c         int_mb(q_map(1,nb)+k-1) = q
c         int_mb(p_map(1,nb)+k-1) = p
c         if (p .eq. taskid) nq(nb) = q
c
c         p        = p+1
c         if (p .ge. np) then
c            p = 0
c            q = q + 1
c         end if
c       end do
c      else if (nb.eq.2) then
c       p = 0
c       q = 1
c       do k=1,nz(1)
c         int_mb(q_map(1,nb)+k-1)       = int_mb(q_map(1,1)+k-1)
c         int_mb(q_map(1,nb)+k+nz(1)-1) = int_mb(q_map(1,1)+k-1)+nq(1)
c         int_mb(p_map(1,nb)+k-1)       = int_mb(p_map(1,1)+k-1)
c         int_mb(p_map(1,nb)+k+nz(1)-1) = int_mb(p_map(1,1)+k-1)
c       end do
c
c      end if

*     ***** block  ******
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

      !*** not used anymore!! ****
      do k=1,nz(nb)
         if (int_mb(p_map(1,nb)+k-1) .eq. taskid) then
c           k_map(q_map(k)) = k
            int_mb(k_map(1,nb)+int_mb(q_map(1,nb)+k-1)-1) = k
         end if
      end do

      nfft3d(nb)     = (nx(nb)/2+1)*ny(nb)*nq(nb)
      n2ft3d(nb)     = 2*nfft3d(nb)
      nfft3d_map(nb) = nfft3d(nb)
      n2ft3d_map(nb) = n2ft3d(nb)


*     ******************************
*     ****** Hilbert mappings ******
*     ******************************
      else





*     **** allocate q_map1,p_map1,q_map2,p_map2,q_map3,p_map3 ****
      value =           BA_alloc_get(mt_int,ny(nb)*nz(nb),
     >                              'q_map1',
     >                               q_map1(2,nb),
     >                               q_map1(1,nb))


      value = value.and.BA_alloc_get(mt_int,ny(nb)*nz(nb),
     >                              'p_map1',
     >                               p_map1(2,nb),
     >                               p_map1(1,nb))

      value = value.and.BA_alloc_get(mt_int,nz(nb)*(nx(nb)/2+1),
     >                              'q_map2',
     >                               q_map2(2,nb),
     >                               q_map2(1,nb))
      value = value.and.BA_alloc_get(mt_int,nz(nb)*(nx(nb)/2+1),
     >                              'p_map2',
     >                               p_map2(2,nb),
     >                               p_map2(1,nb))

      value = value.and.BA_alloc_get(mt_int,ny(nb)*(nx(nb)/2+1),
     >                              'q_map3',
     >                               q_map3(2,nb),
     >                               q_map3(1,nb))
      value = value.and.BA_alloc_get(mt_int,ny(nb)*(nx(nb)/2+1),
     >                              'p_map3',
     >                               p_map3(2,nb),
     >                               p_map3(1,nb))
      if (.not. value)
     > call errquit('Mapping_init:out of heap memory',1, MA_ERR)


      !**** double grid map1 defined wrt to single grid         ****
      !**** makes expand and contract routines trivial parallel ****

!MATHIAS
      if (mapping2d.eq.1) then
         if ((nb.eq.1).or.(nb.gt.2)) then
           call hilbert2d_map(ny(nb),nz(nb),int_mb(p_map1(1,nb)))
         end if
         call hilbert2d_map(nz(nb),(nx(nb)/2+1),int_mb(p_map2(1,nb)))
         call hilbert2d_map((nx(nb)/2+1),ny(nb),int_mb(p_map3(1,nb)))
      else
         if ((nb.eq.1).or.(nb.gt.2)) then
           call hcurve_map(ny(nb),nz(nb),int_mb(p_map1(1,nb)))
         end if
        
         call hcurve_map(nz(nb),(nx(nb)/2+1),int_mb(p_map2(1,nb)))
         call hcurve_map((nx(nb)/2+1),ny(nb),int_mb(p_map3(1,nb)))
      end if



c!$OMP critical
c      write(*,*) "checking p_map1,q_map1 ",Parallel_threadid()
c      do k=1,ny(nb)*nz(nb)
c         write(*,*) Parallel_threadid(),k,
c     >              int_mb(p_map1(1,nb)+k-1)
c      end do
c!$OMP end critical
c!$OMP critical
c      write(*,*) "checking p_map2,q_map2 ",Parallel_threadid()
c      do k=1,(nx(nb)/2+1)*nz(nb)
c         write(*,*) Parallel_threadid(),k,
c     >              int_mb(p_map2(1,nb)+k-1)
c      end do
c!$OMP end critical
c!$OMP critical
c      write(*,*) "checking p_map2,q_map2 ",Parallel_threadid()
c      do k=1,(nx(nb)/2+1)*ny(nb)
c         write(*,*) Parallel_threadid(),k,
c     >              int_mb(p_map3(1,nb)+k-1)
c      end do
c!$OMP end critical


      !**** double grid map1 defined wrt to single grid         ****
      !**** makes expand and contract routines trivial parallel ****
      if ((nb.eq.1).or.(nb.gt.2)) then
      call generate_map_indexes(taskid,np,
     >                          ny(nb),nz(nb),
     >                          int_mb(p_map1(1,nb)),
     >                          int_mb(q_map1(1,nb)),nq1(nb))
      else
        nq1(2) = 4*nq1(1)
        call expand_hilbert2d(np,ny(1),nz(1),
     >                        int_mb(p_map1(1,1)),int_mb(q_map1(1,1)),
     >                        int_mb(p_map1(1,2)),int_mb(q_map1(1,2)))
      end if
      call generate_map_indexes(taskid,np,
     >                          nz(nb),nx(nb)/2+1,
     >                          int_mb(p_map2(1,nb)),
     >                          int_mb(q_map2(1,nb)),nq2(nb))
      call generate_map_indexes(taskid,np,
     >                          nx(nb)/2+1,ny(nb),
     >                          int_mb(p_map3(1,nb)),
     >                          int_mb(q_map3(1,nb)),nq3(nb))

c      if (taskid.eq.0) then
c      write(*,*) taskid,"nq2=",nq2(nb), ny(nb)*nq2(nb)
c      write(*,*) taskid,"nq1=",nq1(nb),(nx(nb)/2+1)*nq1(nb)
c      write(*,*) taskid,"nq3=",nq3(nb), nz(nb)*nq3(nb)
c
c      write(*,*) 'hilbert map nb=',nb
c      do j=0,ny(nb)-1
c        write(*,'(A,80I4)') 'hilbert map:',
c     >   j,(int_mb(p_map3(1,nb)+j*(nx(nb)/2+1)))
c      end do
c      write(*,*)
c      
c      write(*,*) 'hilbert map nb=',nb
c      do j=0,ny(nb)-1
c        write(*,'(A,80I4)') 'hilbert map:',
c     >   (int_mb(p_map3(1,nb)+k+j*(nx(nb)/2+1)), k=0,nx(nb)/2)
c      end do
c      write(*,*)
c      end if

      nfft3d(nb) = (nx(nb)/2+1)*nq1(nb)
      if ((ny(nb)*nq2(nb)).gt.nfft3d(nb)) nfft3d(nb) = ny(nb)*nq2(nb)
      if ((nz(nb)*nq3(nb)).gt.nfft3d(nb)) nfft3d(nb) = nz(nb)*nq3(nb)
      n2ft3d(nb) = 2*nfft3d(nb)

      nfft3d_map(nb) = nz(nb)*nq3(nb)
      n2ft3d_map(nb) = (nx(nb)+2)*nq1(nb)


      
      end if

 
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

#include "bafdecls.fh"
#include "errquit.fh"
#include "D3dB.fh"


*     *** hilbert tranpose data structure ****
      integer h_iq_to_i1(2,6,NBLOCKS)
      integer h_iq_to_i2(2,6,NBLOCKS)
      integer h_i1_start(2,6,NBLOCKS)
      integer h_i2_start(2,6,NBLOCKS)
      common / trans_blk_ijk / h_iq_to_i1,
     >                         h_iq_to_i2,
     >                         h_i1_start,
     >                         h_i2_start


c     integer iq_to_i1((NFFT1/2+1)*NFFT2*NSLABS)
c     integer iq_to_i2((NFFT1/2+1)*NFFT2*NSLABS)
c     integer i1_start(NFFT3+1)
c     integer i2_start(NFFT3+1)
      integer iq_to_i1(2,NBLOCKS)
      integer iq_to_i2(2,NBLOCKS)
      integer i1_start(2,NBLOCKS)
      integer i2_start(2,NBLOCKS)
      common / trans_blk / iq_to_i1,iq_to_i2,i1_start,i2_start

#ifndef MPI
      integer Nchannels(NBLOCKS)
      integer channel_proc(2,NBLOCKS)
      integer channel_type(2,NBLOCKS)
      common / channel_blk / channel_proc,channel_type,Nchannels
#endif

      logical value
      integer i
      logical  control_single_precision_on
      external control_single_precision_on

      call D3dB_timereverse_end(nb)
      call D3dB_fft_end(nb)
      if (control_single_precision_on()) call D3dBs_fft_end(nb)

      value=.true.

      !**** slab mapping ****
      if (mapping.eq.1) then
      value = value.and.BA_free_heap(q_map(2,nb))
      value = value.and.BA_free_heap(p_map(2,nb))
      value = value.and.BA_free_heap(k_map(2,nb))
      end if

      !**** hilbert mappings ****
      if (mapping.eq.2) then
      value = value.and.BA_free_heap(q_map1(2,nb))
      value = value.and.BA_free_heap(p_map1(2,nb))
      value = value.and.BA_free_heap(q_map2(2,nb))
      value = value.and.BA_free_heap(p_map2(2,nb))
      value = value.and.BA_free_heap(q_map3(2,nb))
      value = value.and.BA_free_heap(p_map3(2,nb))
      end if

      !**** slab transpose mappings ****
      if (mapping.eq.1) then
      value = value.and.BA_free_heap(i1_start(2,nb))
      value = value.and.BA_free_heap(i2_start(2,nb))
      value = value.and.BA_free_heap(iq_to_i1(2,nb))
      value = value.and.BA_free_heap(iq_to_i2(2,nb))
      end if

      !**** hilbert transpose mappings ****
      if (mapping.eq.2) then
      do i=1,6
      value = value.and.BA_free_heap(h_i1_start(2,i,nb))
      value = value.and.BA_free_heap(h_i2_start(2,i,nb))
      value = value.and.BA_free_heap(h_iq_to_i1(2,i,nb))
      value = value.and.BA_free_heap(h_iq_to_i2(2,i,nb))
      end do
      end if

#ifndef MPI
      value = value.and.BA_free_heap(channel_proc(2,nb))
      value = value.and.BA_free_heap(channel_type(2,nb))
#endif

      if (.not. value)
     > call errquit('D3dB_end:freeing heap memory',0, MA_ERR)




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

#include "bafdecls.fh"
#include "D3dB.fh"

      
      
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

#include "bafdecls.fh"
#include "D3dB.fh"

      
      
c     q = q_map(k)
c     p = p_map(k)

      q = int_mb(q_map(1,nb)+k-1)
      p = int_mb(p_map(1,nb)+k-1)
      return
      end


*     ***********************************
*     *					*
*     *	          D3dB_ijktoindexp	*	
*     *					*
*     ***********************************

      subroutine D3dB_ijktoindexp(nb,i,j,k,indx,p)      
      implicit none
      integer nb
      integer i,j,k
      integer indx,p

#include "bafdecls.fh"
#include "D3dB.fh"

      integer q

      !**** slab mapping ***
      if (mapping.eq.1) then
      q = int_mb(q_map(1,nb)+k-1)
      p = int_mb(p_map(1,nb)+k-1)

      indx = i + (j-1)*(nx(nb)/2+1) + (q-1)*(nx(nb)/2+1)*ny(nb)

      !**** hilbert mapping ****
      else
      q = int_mb(q_map3(1,nb)+(i-1)+(j-1)*(nx(nb)/2+1))
      p = int_mb(p_map3(1,nb)+(i-1)+(j-1)*(nx(nb)/2+1))
       
      indx = k + (q-1)*nz(nb)
      end if

      return
      end



*     ***********************************
*     *                                 *
*     *           D3dB_ijktoindex1p     *
*     *                                 *
*     ***********************************

      subroutine D3dB_ijktoindex1p(nb,i,j,k,indx,p)
      implicit none
      integer nb
      integer i,j,k
      integer indx,p

#include "bafdecls.fh"
#include "D3dB.fh"

      integer q

      !**** slab mapping ***
      if (mapping.eq.1) then
      q = int_mb(q_map(1,nb)+j-1)
      p = int_mb(p_map(1,nb)+j-1)

      indx = i + (k-1)*(nx(nb)/2+1) + (q-1)*(nx(nb)/2+1)*nz(nb)

      !**** hilbert mapping ****
      else
      q = int_mb(q_map2(1,nb)+(k-1)+(i-1)*(nz(nb)))
      p = int_mb(p_map2(1,nb)+(k-1)+(i-1)*(nz(nb)))

      indx = j + (q-1)*ny(nb)
      end if

      return
      end




*     ***********************************
*     *                                 *
*     *           D3dB_ijktoindex2p     *
*     *                                 *
*     ***********************************

      subroutine D3dB_ijktoindex2p(nb,i,j,k,indx,p)
      implicit none
      integer nb
      integer i,j,k
      integer indx,p

#include "bafdecls.fh"
#include "D3dB.fh"


      integer q

      !**** slab mapping ****
      if (mapping.eq.1) then
      q = int_mb(q_map(1,nb)+j-1)
      p = int_mb(p_map(1,nb)+j-1)

      indx = i + (k-1)*(nx(nb)+2) + (q-1)*(nx(nb)+2)*ny(nb)

      !**** hilbert mapping ****
      else
      q = int_mb(q_map1(1,nb)+(j-1)+(k-1)*ny(nb))
      p = int_mb(p_map1(1,nb)+(j-1)+(k-1)*ny(nb))
      
      indx = i + (q-1)*(nx(nb)+2)
      end if

      return
      end





*     ***********************************
*     *					*
*     *	        D3dB_nfft3d		*
*     *					*
*     ***********************************

      subroutine D3dB_nfft3d(nb,nfft3d_out)
      implicit none
      integer nb
      integer nfft3d_out

#include "D3dB.fh"
      nfft3d_out = nfft3d(nb)
      return
      end


*     ***********************************
*     *                                 *
*     *         D3dB_nfft3d_map         *
*     *                                 *
*     ***********************************

      subroutine D3dB_nfft3d_map(nb,nfft3d_out)
      implicit none
      integer nb
      integer nfft3d_out

#include "D3dB.fh"

      nfft3d_out = nfft3d_map(nb)
      return
      end


*     ***********************************
*     *					*
*     *	        D3dB_n2ft3d		*
*     *					*
*     ***********************************

      subroutine D3dB_n2ft3d(nb,n2ft3d_out)
      implicit none
      integer nb
      integer n2ft3d_out

#include "D3dB.fh"

      n2ft3d_out = n2ft3d(nb)
      return
      end


*     ***********************************
*     *                                 *
*     *         D3dB_n2ft3d_map         *
*     *                                 *
*     ***********************************

      subroutine D3dB_n2ft3d_map(nb,n2ft3d_out)
      implicit none
      integer nb
      integer n2ft3d_out

#include "D3dB.fh"

      n2ft3d_out = n2ft3d_map(nb)
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

#include "D3dB.fh"


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

#include "D3dB.fh"

      
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
      
#include "D3dB.fh"


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

#include "D3dB.fh"

      
      nztmp = nz(nb)
      return
      end

*     ***********************************
*     *                                 *
*     *         D3dB_zplane_size        *
*     *                                 *
*     ***********************************

      subroutine D3dB_zplane_size(nb,zplane_sizetmp)
      implicit none
      integer nb
      integer zplane_sizetmp

#include "D3dB.fh"

      zplane_sizetmp = zplane_size(nb)
      return
      end

*     ***********************************
*     *					*
*     *	        D3dB_Init		*	
*     *					*
*     ***********************************

      subroutine D3dB_Init(nb,nx_in,ny_in,nz_in,map_in)
      implicit none
      integer nb
      integer nx_in,ny_in,nz_in
      integer map_in
      logical value, MA_verify_allocator_stuff
      external MA_verify_allocator_stuff

#include "D3dB.fh"


      !**** local variables ****
      integer MASTER
      parameter (MASTER=0)
      integer taskid,np
      integer  Parallel_threadid
      external Parallel_threadid
      logical  control_single_precision_on
      external control_single_precision_on





      call Parallel2d_np_i(np)
      call Parallel_taskid(taskid)

      !**** Make sure ngrid is consistent with mapping ***
      if (map_in.eq.1) then
        if ((np.gt.nz_in).or.(ny_in.ne.nz_in)) then
          if (taskid.eq.MASTER) then
            write(6,*) 'Error: for slab decomposition the',
     >                 ' number of processors must ',
     >                 ' be in the range ( 1 ...ngrid(3)=',
     >                   nz_in,')'
           write(6,*) ' and ngrid(2) == ngrid(3), ',
     >                ' ngrid(2)=',ny_in,
     >                ' ngrid(3)=',nz_in
          end if
          call errquit('D3dB_Init: mapping error',0,0)
        end if
        if (mod(nx_in,2).ne.0) then
          if (taskid.eq.MASTER) then
           write(6,*)
     >      'Error: ngrid(1) must be even (ngrid(1) mod 2 == 0)'
           write(6,*) 'Error: ngrid(1)=',nx_in
          end if
          call errquit('D3dB_Init: slab mapping error',0,0)
        end if
      end if

      if (map_in.ge.2) then
        if (np.gt.(ny_in*nz_in)) then
          if (taskid.eq.MASTER) then
           write(6,*) 'Error: np > MIN(ngrid(2)*ngrid(3),',
     >                ' (ngrid(1)/2+1)*ngrid(2),',
     >                ' (ngrid(1)/2+1)*ngrid(3))'
           write(6,*) 'Error: np > ngrid(2)*ngrid(3)'
           write(6,*) 'Error: for the Hilbert decomposition the',
     >                 ' the number of processors must ',
     >                 ' be in the range ( 1 ...',
     >                   ny_in*nz_in,')'
          end if
          call errquit('D3dB_Init: Hilbert mapping error',0,0)
        end if
        if (np.gt.((nx_in/2+1)*ny_in)) then
          if (taskid.eq.MASTER) then
           write(6,*) 'Error: np > MIN(ngrid(2)*ngrid(3),',
     >                ' (ngrid(1)/2+1)*ngrid(2),',
     >                ' (ngrid(1)/2+1)*ngrid(3))'
           write(6,*) 'Error: np > (ngrid(1)/2+1)*ngrid(2)'
           write(6,*) 'Error: for the Hilbert decomposition the',
     >                 ' the number of processors must ',
     >                 ' be in the range ( 1 ...',
     >                   (nx_in/2+1)*ny_in,')'
          end if
          call errquit('D3dB_Init: Hilbert mapping error',0,0)
        end if
        if (np.gt.((nx_in/2+1)*nz_in)) then
          if (taskid.eq.MASTER) then
           write(6,*) 'Error: np > MIN(ngrid(2)*ngrid(3),',
     >                ' (ngrid(1)/2+1)*ngrid(2),',
     >                ' (ngrid(1)/2+1)*ngrid(3))'
           write(6,*) 'Error: np > (ngrid(1)/2+1)*ngrid(3)'
           write(6,*) 'Error: for the Hilbert decomposition the', 
     >                 ' the number of processors must ',
     >                 ' be in the range ( 1 ...',
     >                   (nx_in/2+1)*nz_in,')'
          end if
          call errquit('D3dB_Init: Hilbert mapping error',0,0)
        end if
        if (mod(nx_in,2).ne.0) then
          if (taskid.eq.MASTER) then
           write(6,*)
     >      'Error: ngrid(1) must be even (ngrid(1) mod 2 == 0)'
           write(6,*) 'Error: ngrid(1)=',nx_in
          end if
          call errquit('D3dB_Init: Hilbert mapping error',0,0)
        end if
      end if


*     ***** initialize D3dB common block *****
      nx(nb)     = nx_in
      ny(nb)     = ny_in
      nz(nb)     = nz_in
      mapping    = map_in
      mapping2d  = 1
      if (mapping.eq.3) then
         mapping   = 2
         mapping2d = 2
      end if
 
      

*     **** do other initializations ****
      call Mapping_Init(nb)

      if (mapping.eq.1) call D3dB_c_transpose_jk_init(nb)
      if (mapping.eq.2) call D3dB_c_transpose_ijk_init(nb)

#ifndef MPI
      call D3dB_channel_init(nb)
#endif

      call D3dB_c_timereverse_init(nb)
      call D3dB_fft_init(nb)
      if (control_single_precision_on()) call D3dBs_fft_init(nb)
      
      return
      end


c*     ***********************************
c*     *					*
c*     *	        D3dB_SumAll		*	
c*     *					*
c*     ***********************************
c
c      subroutine D3dB_SumAll(sum)
cc     implicit none
c      real*8  sum
c
c
c#include "tcgmsg.fh"
c#include "msgtypesf.h"
c
c
c      integer MASTER
c      parameter (MASTER=0)
c      integer msglen
c      real*8 sumall,sumt
c
c*     **** external functions ****
c      integer  Parallel2d_comm_i
c      external Parallel2d_comm_i
c    
cc     msglen = 8
c      msglen = 1
c
c      sumt = sum
c
c      call GA_PGROUP_DGOP(Parallel2d_comm_i(),
c     >                    9+MSGDBL,sumt,1,'+')
c
cc      call GA_DGOP(9+MSGDBL,sumt,1,'+')
cc     call DGOP(9+MSGDBL,sumt,1,'+')
c      sumall=sumt
c
c
c      sum = sumall
c      return
c      end


c*     ***********************************
c*     *					*
c*     *	        D3dB_ISumAll		*	
c*     *					*
c*     ***********************************
c
c      subroutine D3dB_ISumAll(sum)
cc     implicit none
c      integer  sum
c
c
c#include "tcgmsg.fh"
c#include "msgtypesf.h"
c
c
c      integer MASTER
c      parameter (MASTER=0)
c      integer msglen
c      integer sumall,sumt
c
c*     **** external functions ****
c      integer  Parallel2d_comm_i
c      external Parallel2d_comm_i
c
c    
cc     msglen = 8
c      msglen = 1
c
c      sumt = sum
c
c      call GA_PGROUP_IGOP(Parallel2d_comm_i(),
c     >                    9+MSGINT,sumt,1,'+')
cc     call GA_IGOP(9+MSGINT,sumt,1,'+')
c      sumall=sumt
c
c
c      sum = sumall
c      return
c      end



*     ***********************************
*     *					*
*     *	        D3dB_(c,r,t)_Zero	*	
*     *					*
*     ***********************************

      subroutine D3dB_c_Zero(nb,A)     
      implicit none 
      integer nb
      complex*16 A(*)

#include "D3dB.fh"


      !call dcopy(n2ft3d_map(nb),0.0d0,0,A,1)
      call Parallel_shared_vector_zero(.true.,n2ft3d_map(nb),A)
      return
      end

      subroutine D3dB_c_nZero(nb,n,A)
      implicit none
      integer nb,n
      complex*16 A(*)

#include "D3dB.fh"


      !call dcopy(n*n2ft3d_map(nb),0.0d0,0,A,1)
      call Parallel_shared_vector_zero(.true.,n*n2ft3d_map(nb),A)
      return
      end


      subroutine D3dB_r_Zero(nb,A)     
      implicit none 
      integer nb
      real*8  A(*)

#include "D3dB.fh"

c      call dcopy(n2ft3d_map(nb),0.0d0,0,A,1)
      call Parallel_shared_vector_zero(.true.,n2ft3d_map(nb),A)
      return
      end


      subroutine D3dB_r_nZero(nb,n,A)
      implicit none
      integer nb,n
      real*8  A(*)

#include "D3dB.fh"

      !call dcopy(n*n2ft3d_map(nb),0.0d0,0,A,1)
      call Parallel_shared_vector_zero(.true.,n*n2ft3d_map(nb),A)
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

#include "D3dB.fh"

      !call dcopy(2*nfft3d_map(nb),A,1,B,1)
      call Parallel_shared_vector_copy(.true.,2*nfft3d_map(nb),A,B)
      return
      end

      subroutine D3dB_c_nCopy(nb,n,A,B)
      implicit none
      integer nb,n
      complex*16 A(*)
      complex*16 B(*)

#include "D3dB.fh"

      !call dcopy(n*2*nfft3d_map(nb),A,1,B,1)
      call Parallel_shared_vector_copy(.true.,n*2*nfft3d_map(nb),A,B)
      return
      end


      subroutine D3dB_r_Copy(nb,A,B)     
      implicit none 
      integer nb
      real*8 A(*)
      real*8 B(*)

#include "D3dB.fh"
      integer i
      !call dcopy(n2ft3d_map(nb),A,1,B,1)
      call Parallel_shared_vector_copy(.true.,n2ft3d_map(nb),A,B)

      return
      end


      subroutine D3dB_r_nCopy(nb,n,A,B)
      implicit none
      integer nb,n
      real*8 A(*)
      real*8 B(*)

#include "D3dB.fh"

      !call dcopy(n*n2ft3d_map(nb),A,1,B,1)
      call Parallel_shared_vector_copy(.true.,n*n2ft3d_map(nb),A,B)
      return
      end


      subroutine D3dB_t_Copy(nb,A,B)     
      implicit none 
      integer nb
      real*8 A(*)
      real*8 B(*)

#include "D3dB.fh"

      !call dcopy(nfft3d_map(nb),A,1,B,1)
      call Parallel_shared_vector_copy(.true.,nfft3d_map(nb),A,B)
      return
      end


      subroutine D3dB_t_nCopy(nb,n,A,B)
      implicit none
      integer nb,n
      real*8 A(*)
      real*8 B(*)

#include "D3dB.fh"

      !call dcopy(n*nfft3d_map(nb),A,1,B,1)
      call Parallel_shared_vector_copy(.true.,n*nfft3d_map(nb),A,B)
      return
      end


      subroutine D3dB_tc_Copy(nb,A,B)
      implicit none
      integer nb
      real*8     A(*)
      complex*16 B(*)

#include "D3dB.fh"

      integer i

!$OMP DO
      do i=1,nfft3d_map(nb)
         B(i) = dcmplx(A(i),0.0d0)
      end do
!$OMP END DO

      return
      end

      subroutine D3dB_ct_Copy(nb,A,B)
      implicit none
      integer nb
      complex*16 A(*)
      real*8     B(*)

#include "D3dB.fh"

      integer i
!$OMP DO
      do i=1,nfft3d_map(nb)
         B(i) = dble(A(i))
      end do
!$OMP END DO

      return
      end




*     ***********************************
*     *					*
*     *	        D3dB_fft_init		*
*     *					*
*     ***********************************

      subroutine D3dB_fft_init(nb)
      implicit none
      integer nb

#include "bafdecls.fh"
#include "errquit.fh"

#include "D3dB.fh"

      integer tmpx(2,NBLOCKS),tmpy(2,NBLOCKS),tmpz(2,NBLOCKS)
      common    / D3dB_fft / tmpx,tmpy,tmpz

#ifdef FFTW3
#include "fftw3.fh"
       integer Atest(2),nxh,nxhz,nxhy
#endif

      logical value
      integer i,mthr, Parallel_maxthreads
      external Parallel_maxthreads

      mthr = Parallel_maxthreads()
      

      !call D3dB_nfft3d(nb,nfft3d)

c      value = BA_alloc_get(mt_dcpl,(nfft3d(nb)),
c     >        'fttmpx',tmpx(2,nb),tmpx(1,nb))
c      value = value.and.
c     >        BA_alloc_get(mt_dcpl,(nfft3d(nb)),
c     >        'fttmpy',tmpy(2,nb),tmpy(1,nb))
c      value = value.and.
c     >        BA_alloc_get(mt_dcpl,(nfft3d(nb)),
c     >        'fttmpz',tmpz(2,nb),tmpz(1,nb))
      value = BA_alloc_get(mt_dcpl,mthr*(2*nx(nb)+15),
     >        'fttmpx',tmpx(2,nb),tmpx(1,nb))
      value = value.and.
     >        BA_alloc_get(mt_dcpl,mthr*(2*ny(nb)+15),
     >        'fttmpy',tmpy(2,nb),tmpy(1,nb))
      value = value.and.
     >        BA_alloc_get(mt_dcpl,mthr*(2*nz(nb)+15),
     >        'fttmpz',tmpz(2,nb),tmpz(1,nb))
      if (.not. value) 
     >   call errquit('D3dB_fft_init:out of heap memory',0, MA_ERR)


#ifdef MLIB
      call drc1ft(dcpl_mb(tmpx(1,nb)),nx(nb),
     >            dcpl_mb(tmpx(1,nb)),-3,ierr)
      call z1dfft(dcpl_mb(tmpx(1,nb)),ny(nb),
     >            dcpl_mb(tmpy(1,nb)),-3,ierr)
      call z1dfft(dcpl_mb(tmpx(1,nb)),nz(nb),
     >            dcpl_mb(tmpz(1,nb)),-3,ierr)

#else
      do i=1,mthr
        !write(*,*) "DEBUG init fft arrays of thread ", i-1
        call drffti(nx(nb),dcpl_mb(tmpx(1,nb)+(i-1)*(2*nx(nb)+15)))
        call dcffti(ny(nb),dcpl_mb(tmpy(1,nb)+(i-1)*(2*ny(nb)+15)))
        call dcffti(nz(nb),dcpl_mb(tmpz(1,nb)+(i-1)*(2*nz(nb)+15)))
      end do
#endif

#ifdef FFTW3

c       call dfftw_init_threads()
c       call dfftw_plan_with_nthreads(2)

       iforward  = FFTW_FORWARD
       ibackward = FFTW_BACKWARD
c       iestimate = FFTW_PATIENT
c       iestimate = FFTW_MEASURE
c       iestimate = FFTW_ESTIMATE
c       iestimate = FFTW_EXHAUSTIVE
       call icopy(nplans*NBLOCKS,0,0,plans,1)
      if (mapping.eq.1) then
         nxh = (nx(nb)/2+1)
         nxhz = nxh*nz(nb)
         nxhy = nxh*ny(nb)
         if (.not.BA_alloc_get(mt_dcpl,nx(nb)*ny(nb)*nq(nb),
     >                       'Atest',Atest(2),Atest(1)))
     >     call errquit('D3dB_fft_init:out of heap memory',0,MA_ERR)

         call dfftw_plan_many_dft(plans(1,nb),1,nz(nb),nxh,
     >                        dcpl_mb(Atest(1)),nxhz,nxh,1,
     >                        dcpl_mb(Atest(1)),nxhz,nxh,1,
     >                        ibackward,FFTW_EXHAUSTIVE)
         call dfftw_plan_many_dft(plans(2,nb),1,ny(nb),nxh,
     >                        dcpl_mb(Atest(1)),nxhy,nxh,1,
     >                        dcpl_mb(Atest(1)),nxhy,nxh,1, 
     >                        ibackward,FFTW_EXHAUSTIVE)
         call dfftw_plan_many_dft_c2r(plans(3,nb),1,nx(nb),
     >        ny(nb)*nq(nb),
     >        dcpl_mb(Atest(1)),nxh       *ny(nb)*nq(nb),1,nxh,
     >        dcpl_mb(Atest(1)),(nx(nb)+2)*ny(nb)*nq(nb),1,nx(nb)+2, 
     >        FFTW_EXHAUSTIVE)

         call dfftw_plan_many_dft_r2c(plans(4,nb),1,nx(nb),
     >        ny(nb)*nq(nb),
     >        dcpl_mb(Atest(1)),(nx(nb)+2)*ny(nb)*nq(nb),1,nx(nb)+2,
     >        dcpl_mb(Atest(1)),nxh       *ny(nb)*nq(nb),1,nxh,
     >        FFTW_EXHAUSTIVE)
         call dfftw_plan_many_dft(plans(5,nb),1,ny(nb),nxh,
     >                        dcpl_mb(Atest(1)),nxhy,nxh,1,
     >                        dcpl_mb(Atest(1)),nxhy,nxh,1, 
     >                        iforward,FFTW_EXHAUSTIVE)
         call dfftw_plan_many_dft(plans(6,nb),1,nz(nb),nxh,
     >                        dcpl_mb(Atest(1)),nxhz,nxh,1,
     >                        dcpl_mb(Atest(1)),nxhz,nxh,1,
     >                        iforward,FFTW_EXHAUSTIVE)

         call dfftw_plan_many_dft(plans(7,nb),1,nz(nb),1,
     >                        dcpl_mb(Atest(1)),nxhz,nxh,1,
     >                        dcpl_mb(Atest(1)),nxhz,nxh,1,
     >                        ibackward,FFTW_EXHAUSTIVE)
         call dfftw_plan_many_dft(plans(8,nb),1,ny(nb),1,
     >                        dcpl_mb(Atest(1)),nxhy,nxh,1,
     >                        dcpl_mb(Atest(1)),nxhy,nxh,1,
     >                        ibackward,FFTW_EXHAUSTIVE)

          call dfftw_plan_many_dft(plans(9,nb),1,ny(nb),1,
     >                        dcpl_mb(Atest(1)),nxhy,nxh,1,
     >                        dcpl_mb(Atest(1)),nxhy,nxh,1,
     >                        iforward,FFTW_EXHAUSTIVE)
          call dfftw_plan_many_dft(plans(10,nb),1,nz(nb),1,
     >                        dcpl_mb(Atest(1)),nxhz,nxh,1,
     >                        dcpl_mb(Atest(1)),nxhz,nxh,1,
     >                        iforward,FFTW_EXHAUSTIVE)


         if (.not.BA_free_heap(Atest(2)))
     >   call errquit('D3dB_fft_init:freeing heap',0,MA_ERR)

      else

         nxh  = (nx(nb)/2+1)
         if (.not.BA_alloc_get(mt_dcpl,nfft3d(nb),'Atest',
     >                         Atest(2),Atest(1)))
     >     call errquit('D3dB_fft_init:out of heap memory',0,MA_ERR)

         call dfftw_plan_many_dft(plans(11,nb),1,nz(nb),nq3(nb),
     >                        dcpl_mb(Atest(1)),nz(nb)*nq3(nb),1,nz(nb),
     >                        dcpl_mb(Atest(1)),nz(nb)*nq3(nb),1,nz(nb),
     >                        ibackward,FFTW_EXHAUSTIVE)
         call dfftw_plan_many_dft(plans(12,nb),1,ny(nb),nq2(nb),
     >                        dcpl_mb(Atest(1)),ny(nb)*nq2(nb),1,ny(nb),
     >                        dcpl_mb(Atest(1)),ny(nb)*nq2(nb),1,ny(nb),
     >                        ibackward,FFTW_EXHAUSTIVE)
         call dfftw_plan_many_dft_c2r(plans(13,nb),1,nx(nb),
     >        nq1(nb),
     >        dcpl_mb(Atest(1)),nfft3d(nb),1,nxh,
     >        dcpl_mb(Atest(1)),n2ft3d(nb),1,(nx(nb)+2), 
     >        FFTW_EXHAUSTIVE)
         call dfftw_plan_many_dft_r2c(plans(14,nb),1,nx(nb),
     >        nq1(nb),
     >        dcpl_mb(Atest(1)),n2ft3d(nb),1,nx(nb)+2,
     >        dcpl_mb(Atest(1)),nfft3d(nb),1,nxh,
     >        FFTW_EXHAUSTIVE)

         call dfftw_plan_many_dft(plans(15,nb),1,ny(nb),nq2(nb),
     >                        dcpl_mb(Atest(1)),ny(nb)*nq2(nb),1,ny(nb),
     >                        dcpl_mb(Atest(1)),ny(nb)*nq2(nb),1,ny(nb),
     >                        iforward,FFTW_EXHAUSTIVE)
         call dfftw_plan_many_dft(plans(16,nb),1,nz(nb),nq3(nb),
     >                        dcpl_mb(Atest(1)),nz(nb)*nq3(nb),1,nz(nb),
     >                        dcpl_mb(Atest(1)),nz(nb)*nq3(nb),1,nz(nb),
     >                        iforward,FFTW_EXHAUSTIVE)

         call dfftw_plan_many_dft(plans(17,nb),1,nz(nb),1,
     >                        dcpl_mb(Atest(1)),nz(nb),1,1,
     >                        dcpl_mb(Atest(1)),nz(nb),1,1,
     >                        ibackward,FFTW_EXHAUSTIVE)
         call dfftw_plan_many_dft(plans(18,nb),1,ny(nb),1,
     >                        dcpl_mb(Atest(1)),ny(nb),1,1,
     >                        dcpl_mb(Atest(1)),ny(nb),1,1,
     >                        ibackward,FFTW_EXHAUSTIVE)

          call dfftw_plan_many_dft(plans(19,nb),1,ny(nb),1,
     >                        dcpl_mb(Atest(1)),ny(nb),1,1,
     >                        dcpl_mb(Atest(1)),ny(nb),1,1,
     >                        iforward,FFTW_EXHAUSTIVE)
          call dfftw_plan_many_dft(plans(20,nb),1,nz(nb),1,
     >                        dcpl_mb(Atest(1)),nz(nb),1,1,
     >                        dcpl_mb(Atest(1)),nz(nb),1,1,
     >                        iforward,FFTW_EXHAUSTIVE)

         if (.not.BA_free_heap(Atest(2)))
     >   call errquit('D3dB_fft_init:freeing heap',0,MA_ERR)

      end if
#endif

      return
      end

*     ***********************************
*     *                                 *
*     *         D3dB_fft_end            *
*     *                                 *
*     ***********************************

      subroutine D3dB_fft_end(nb)
      implicit none
      integer nb

#include "bafdecls.fh"
#include "errquit.fh"

#include "D3dB.fh"


      integer tmpx(2,NBLOCKS),tmpy(2,NBLOCKS),tmpz(2,NBLOCKS)
      common    / D3dB_fft / tmpx,tmpy,tmpz
#ifdef USE_OPENMP
      logical is_par
      common    / Debug_openmp / is_par
#endif

      logical value
      integer i

#ifdef FFTW3
      do i=1,nplans
         if (plans(i,nb).ne.0) call dfftw_destroy_plan(plans(i,nb))
      end do
c      call dfftw_cleanup_threads()
#endif

      value =           BA_free_heap(tmpx(2,nb))
      value = value.and.BA_free_heap(tmpy(2,nb))
      value = value.and.BA_free_heap(tmpz(2,nb))
      if (.not.value)
     >   call errquit(
     >   'D3dB_fft_end:error deallocatingof heap memory',0, MA_ERR)

      return
      end







c*     ***********************************
c*     *                                 *
c*     *         D3dB_n_fft_init         *
c*     *                                 *
c*     ***********************************
c
c      subroutine D3dB_n_fft_init(nb,ne)
c      implicit none
c      integer nb,ne
c
c#include "bafdecls.fh"
c#include "errquit.fh"
c
c      integer tmp2(2,NBLOCKS),tmp3(2,NBLOCKS)
c      common    / D3dB_n_fft / tmp2,tmp3
c
c      logical value
c      integer nfft3d
c
c      call D3dB_nfft3d(nb,nfft3d)
c      value = BA_alloc_get(mt_dcpl,(ne*nfft3d),
c     >        'fttmp2_h',tmp2(2,nb),tmp2(1,nb))
c      value = value.and.
c     >        BA_alloc_get(mt_dbl,(2*ne*nfft3d),
c     >        'fttmp3_h',tmp3(2,nb),tmp3(1,nb))
c      if (.not. value)
c     > call errquit('D3dB_n_fft_init:out of heap memory',0, MA_ERR)
c
c      return
c      end
c
c*     ***********************************
c*     *                                 *
c*     *         D3dB_n_fft_end          *
c*     *                                 *
c*     ***********************************
c
c      subroutine D3dB_n_fft_end(nb)
c      implicit none
c      integer nb
c
c#include "bafdecls.fh"
c#include "errquit.fh"
c
c      integer tmp2(2,NBLOCKS),tmp3(2,NBLOCKS)
c      common    / D3dB_n_fft / tmp2,tmp3
c
c      logical value
c
c      value =           BA_free_heap(tmp2(2,nb))
c      value = value.and.BA_free_heap(tmp3(2,nb))
c      if (.not.value)
c     > call errquit(
c     > 'D3dB_n_fft_end:error deallocatingof heap memory',0,MA_ERR)
c
c      return
c      end



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


#include "bafdecls.fh"
#include "errquit.fh"

#include "D3dB.fh"


      integer tmpx(2,NBLOCKS),tmpy(2,NBLOCKS),tmpz(2,NBLOCKS)
      common    / D3dB_fft / tmpx,tmpy,tmpz

*     *** local variables ***
      integer i,j,k,q,indx
      integer nxh,nxhy,nxhz,indx0,indx1

      
      !integer tmp1(2),tmp2(2),tmp3(2)
      integer tmp2(2),tmp3(2)
      logical value

      integer  tid,nthr,offset
      integer  Parallel_threadid,Parallel_nthreads
      external Parallel_threadid,Parallel_nthreads


      tid  = Parallel_threadid()
      nthr = Parallel_nthreads()

      call nwpw_timing_start(1)

*     ***** allocate temporary space ****
      !call D3dB_nfft3d(nb,nfft3d)
      value = BA_push_get(mt_dcpl,(nfft3d(nb)),'ffttmp2',
     >                    tmp2(2),tmp2(1))
      value = value.and.
     >      BA_push_get(mt_dbl,(n2ft3d(nb)),'ffttmp3',tmp3(2),tmp3(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)

       nxh = (nx(nb)/2+1)
       nxhz = nxh*nz(nb)
       nxhy = nxh*ny(nb)
      
      !**********************
      !**** slab mapping ****
      !**********************
      if (mapping.eq.1) then
*     ********************************************
*     ***         Do a transpose of A          ***
*     ***      A(kx,kz,ky) <- A(kx,ky,kz)      ***
*     ********************************************
c     call D3dB_c_transpose_jk(nb,A,dcpl_mb(tmp2(1)),dbl_mb(tmp3(1)))

*     *************************************************
*     ***     do fft along kz dimension             ***
*     ***   A(kx,nz(nb),ky) <- fft1d^(-1)[A(kx,kz,ky)]  ***
*     *************************************************
#ifdef MLIB
      !call z1dfft(dbl_mb(tmp3(1)),nz(nb),dcpl_mb(tmpz(1)),-3,ierr)
      do q=1,nq(nb)
      do i=1,(nx(nb)/2+1)
         do k=1,nz(nb)
            indx = i + (k-1)*(nx(nb)/2+1) + (q-1)*(nx(nb)/2+1)*nz(nb)
            dcpl_mb(tmp2(1)+k-1) = A(indx)
         end do
         call z1dfft(dcpl_mb(tmp2(1)),nz(nb),
     >               dcpl_mb(tmpz(1,nb)),-2,ierr)
         do k=1,nz(nb)
            indx = i + (k-1)*(nx(nb)/2+1) + (q-1)*(nx(nb)/2+1)*nz(nb)
            A(indx) = dcpl_mb(tmp2(1)+k-1)
         end do
      end do
      end do
      !call dscal((nx(nb)+2)*ny(nb)*nq(nb),dble(nz(nb)),A,1)

#else

#ifdef FFTW3
      do q=1,nq(nb)
        indx = 1+(q-1)*nxhz
        call dfftw_execute_dft(plans(1,nb),A(indx),A(indx))
      end do
#else
      !call dcffti(nz(nb),dcpl_mb(tmp1(1)))
      indx0 = 0
      do q=1,nq(nb)
      do i=1,nxh

         indx  = i + indx0
         indx1 = indx
         do k=1,nz(nb)
            dcpl_mb(tmp2(1)+k-1) = A(indx)
            indx = indx + nxh
         end do
         call dcfftb(nz(nb),dcpl_mb(tmp2(1)),dcpl_mb(tmpz(1,nb)))
         do k=1,nz(nb)
            A(indx1) = dcpl_mb(tmp2(1)+k-1)
            indx1 = indx1 + nxh
         end do

      end do
      indx0 = indx0 + nxhz
      end do
#endif
#endif

*     ********************************************
*     ***         Do a transpose of A          ***
*     ***      A(kx,ky,nz(nb)) <- A(kx,nz(nb),ky)      ***
*     ********************************************
      call D3dB_c_transpose_jk(nb,A,dcpl_mb(tmp2(1)),dbl_mb(tmp3(1)))

*     *************************************************
*     ***     do fft along ky dimension             ***
*     ***   A(kx,ny(nb),nz(nb)) <- fft1d^(-1)[A(kx,ky,nz(nb))]  ***
*     *************************************************
#ifdef MLIB
      !call z1dfft(dbl_mb(tmp3(1)),ny(nb),dcpl_mb(tmp1(1)),-3,ierr)
      do q=1,nq(nb)
      do i=1,(nx(nb)/2+1)
         do j=1,ny(nb)
            indx = i + (j-1)*(nx(nb)/2+1) + (q-1)*(nx(nb)/2+1)*ny(nb)
            dcpl_mb(tmp2(1)+j-1) = A(indx)
         end do
         call z1dfft(dcpl_mb(tmp2(1)),ny(nb),
     >               dcpl_mb(tmpy(1,nb)),-2,ierr)
         do j=1,ny(nb)
            indx = i + (j-1)*(nx(nb)/2+1) + (q-1)*(nx(nb)/2+1)*ny(nb)
            A(indx) = dcpl_mb(tmp2(1)+j-1)
         end do
      end do
      end do
      !call dscal((nx(nb)+2)*ny(nb)*nq(nb),dble(ny(nb)),A,1)
#else

#ifdef FFTW3
      do q=1,nq(nb)
         indx = 1+(q-1)*nxhy
         call dfftw_execute_dft(plans(2,nb),A(indx),A(indx))
      end do
#else
      !call dcffti(ny(nb),dcpl_mb(tmp1(1)))
      indx0 = 0
      do q=1,nq(nb)
      do i=1,(nx(nb)/2+1)
      
         indx  = i + indx0
         indx1 = indx
         do j=1,ny(nb)
            dcpl_mb(tmp2(1)+j-1) = A(indx)
            indx = indx + nxh
         end do
         call dcfftb(ny(nb),dcpl_mb(tmp2(1)),dcpl_mb(tmpy(1,nb)))
         do j=1,ny(nb)
            A(indx1) = dcpl_mb(tmp2(1)+j-1)
            indx1 = indx1 + nxh
         end do

      end do
      indx0 = indx0 + nxhy
      end do
#endif
#endif

*     *************************************************
*     ***     do fft along kx dimension             ***
*     ***   A(nx(nb),ny(nb),nz(nb)) <- fft1d^(-1)[A(kx,ny(nb),nz(nb))]  ***
*     *************************************************
#ifdef MLIB
      !call drc1ft (dbl_mb(tmp3(1)),nx(nb),dcpl_mb(tmp1(1)),-3,ierr)
      do q=1,nq(nb)
      do j=1,ny(nb)
         indx = 1 + (j-1)*(nx(nb)/2+1) + (q-1)*(nx(nb)/2+1)*ny(nb)
         call drc1ft(A(indx),nx(nb),dcpl_mb(tmpx(1,nb)),-2,ierr)
      end do
      end do
c     call drcfts(A,nx(nb),1,ny(nb)*nq(nb),
c    >                  nx(nb)+2,-2,ierr)
c     call dscal((nx(nb)+2)*ny(nb)*nq(nb),dble(nx(nb)),A,1)

#else

#ifdef FFTW3
      call dfftw_execute_dft_c2r(plans(3,nb),A,A)

#else
      !call drffti(nx(nb),dcpl_mb(tmp1(1)))

c      do q=1,nq(nb)
c      do j=1,ny(nb)
c         indx = 1 + (j-1)*(nx(nb)/2+1) + (q-1)*(nx(nb)/2+1)*ny(nb)
c         call dcopy((nx(nb)+2),A(indx),1,dbl_mb(tmp3(1)),1)
c         do i=2,nx(nb)
c            dbl_mb(tmp3(1)+i-1) = dbl_mb(tmp3(1)+i)
c         end do
c         call drfftb(nx(nb),dbl_mb(tmp3(1)),dcpl_mb(tmp1(1)))
c         dbl_mb(tmp3(1)+nx(nb)) = 0.0d0
c         dbl_mb(tmp3(1)+nx(nb)+1) = 0.0d0
c         call dcopy((nx(nb)+2),dbl_mb(tmp3(1)),1,A(indx),1)
c      end do
c      end do

      call cshift1_fftb(nx(nb),ny(nb),nq(nb),1,A)
      indx = 1
      do q=1,nq(nb)
      do j=1,ny(nb)
         !indx = 1 + (j-1)*(nx(nb)/2+1) + (q-1)*(nx(nb)/2+1)*ny(nb)
         call drfftb(nx(nb),A(indx),dcpl_mb(tmpx(1,nb)))
         indx = indx + nxh
      end do
      end do
      call zeroend_fftb(nx(nb),ny(nb),nq(nb),1,A)
#endif

#endif


      !*************************
      !**** hilbert mapping ****
      !*************************
      else


*     *************************************************
*     ***     do fft along kz dimension             ***
*     ***   A(nz(nb),kx,ky) <- fft1d^(-1)[A(kz,kx,ky)]  ***
*     *************************************************
#ifdef MLIB
      indx = 1
      do q=1,nq3(nb)
         !indx = 1 + (q-1)*nz(nb)
         call z1dfft(A(indx),nz(nb),dcpl_mb(tmpz(1,nb)),-2,ierr)
         indx = indx + nz(nb)
      end do
#else

#ifdef FFTW3
      call dfftw_execute_dft(plans(11,nb),A,A)

#else

      offset=tid*(2*nz(nb)+15)
      do i=tid+1,nq3(nb),nthr
        call dcfftb(nz(nb),A(1+(i-1)*nz(nb)),dcpl_mb(tmpz(1,nb)+offset))
      end do
!$OMP BARRIER

#endif
#endif

      call D3dB_c_transpose_ijk(nb,3,A,dcpl_mb(tmp2(1)),dbl_mb(tmp3(1)))

*     *************************************************
*     ***     do fft along ky dimension             ***
*     ***   A(ny(nb),nz(nb),kx) <- fft1d^(-1)[A(ky,nz(nb),kx)]  ***
*     *************************************************
#ifdef MLIB
      indx = 1
      do q=1,nq2(nb)
         !indx = 1 + (q-1)*ny(nb)
         call z1dfft(A(indx),ny(nb),dcpl_mb(tmpy(1,nb)),-2,ierr)
         indx = indx + ny(nb)
      end do
#else

#ifdef FFTW3
      call dfftw_execute_dft(plans(12,nb),A,A)

#else
      offset=tid*(2*ny(nb)+15)
      do i=tid+1,nq2(nb),nthr
        call dcfftb(ny(nb),A(1+(i-1)*ny(nb)),dcpl_mb(tmpy(1,nb)+offset))
      end do
!$OMP BARRIER
#endif
#endif

      call D3dB_c_transpose_ijk(nb,4,A,dcpl_mb(tmp2(1)),dbl_mb(tmp3(1)))

*     *************************************************
*     ***     do fft along kx dimension             ***
*     ***   A(nx(nb),ny(nb),nz(nb)) <- fft1d^(-1)[A(kx,ny(nb),nz(nb))]  ***
*     *************************************************
#ifdef MLIB
      indx = 1
      do q=1,nq1(nb)
         !indx = 1 + (q-1)*(nx(nb)/2+1)
         call drc1ft(A(indx),nx(nb),dcpl_mb(tmpx(1,nb)),-2,ierr)
         indx = indx + nxh
      end do
#else

#ifdef FFTW3
      call dfftw_execute_dft_c2r(plans(13,nb),A,A)

#else
      offset=tid*(2*nx(nb)+15)
      call cshift1_fftb(nx(nb),nq1(nb),1,1,A)
      do i=tid+1,nq1(nb),nthr
        call drfftb(nx(nb),A(1+(i-1)*nxh),dcpl_mb(tmpx(1,nb)+offset))
      end do
!$OMP BARRIER

      call zeroend_fftb(nx(nb),nq1(nb),1,1,A)

#endif
#endif

      end if

    
*     **** deallocate temporary space  ****
      value = BA_pop_stack(tmp3(2))
      value = value.and.BA_pop_stack(tmp2(2))
      !value = BA_pop_stack(tmp1(2))
      if (.not. value) call errquit('popping stack memory',0,MA_ERR)

      call nwpw_timing_end(1)

      return
      end

      subroutine D3dB_fftbx_sub(n,nx,nxh,tmpx,A)
      implicit none
      integer n,nx,nxh
      real*8     tmpx(2*nx+15)
      complex*16 A(nxh,n)
      integer i



      do i=1,n
         call drfftb(nx,A(1,i),tmpx)
      end do


      return
      end

      subroutine D3dB_fftby_sub2(n,ny,tmpy,A)
      implicit none
      integer n,ny
      real*8     tmpy(4*ny+15)
      complex*16 A(ny,n)
      integer i



      do i=1,n
         call dcfftb(ny,A(1,i),tmpy)
      end do


      return
      end

      subroutine D3dB_fftbz_sub2(n,nz,tmpz,A)
      implicit none
      integer n,nz
      real*8     tmpz(4*nz+15)
      complex*16 A(nz,n)
      integer i



      do i=1,n
         call dcfftb(nz,A(1,i),tmpz)
      end do


      return
      end

      subroutine cshift1_fftb(nx,ny,nq,ne,A)
      implicit none
      integer nx,ny,nq,ne
      real*8 A(*)

      integer i,j,indx

!$OMP DO
      do j=1,(ny*nq*ne)
        indx = 1+(j-1)*(nx+2)
c        indx = 1 + (j-1)*(nx+2) + (q-1)*(nx+2)*ny
c    >            + (n-1)*(nx+2)*ny*nq
         do i=2,nx
            A(indx+i-1) = A(indx+i)
         end do
      end do
!$OMP END DO

c     end do
c     end do
      return
      end 


      subroutine cshift1_fftb1(nx,A)
      implicit none
      integer nx
      real*8 A(*)
      integer i
      do i=2,nx
         A(i) = A(i+1)
      end do 
      return
      end


      subroutine zeroend_fftb1(nx,A)
      implicit none
      integer nx
      real*8 A(*)
      integer i
      A(nx+1) = 0.0d0
      A(nx+2) = 0.0d0
      return
      end


      subroutine zeroend_fftb(nx,ny,nq,ne,A)
      implicit none
      integer nx,ny,nq,ne
      real*8 A(*)

      integer i,indx

!$OMP DO
      do i=1,(ny*nq*ne)
         indx = nx+1+(i-1)*(nx+2)
         A(indx)   = 0.0d0
         A(indx+1) = 0.0d0
      end do
!$OMP END DO

      return
      end 

c*     ***********************************
c*     *					*
c*     *	        D3dB_ncr_fft3b		*
c*     *					*
c*     ***********************************
c
c      subroutine D3dB_ncr_fft3b(nb,ne,A)
c
c*****************************************************
c*                                                   *
c*      This routine performs the operation of       *
c*      a three dimensional complex to complex       *
c*      inverse fft                                  *
c* A(nx,ny(nb),nz(nb),n) <- FFT3^(-1)[A(kx,ky,kz,n)] * 
c*                                                   *
c*      Entry - 					    *
c*              A: a column distribuded 3d block     *
c*              tmp: tempory work space must be at   *
c*                    least the size of (complex)    *
c*                    (nfft*nfft + 1) + 10*nfft      * 
c*                                                   *
*       Exit - A is transformed and the imaginary   *
cc*              part of A is set to zero             *
c*       uses - D3dB_c_transpose_jk, dcopy           *
c*                                                   *
c*****************************************************
c
c      implicit none
c      integer nb,ne
c      complex*16  A(*)
c
c#include "bafdecls.fh"
c#include "errquit.fh"
c
c#include "D3dB.fh"
c
c
c      integer tmpx(2,NBLOCKS),tmpy(2,NBLOCKS),tmpz(2,NBLOCKS)
c      common    / D3dB_fft / tmpx,tmpy,tmpz
c
c      integer tmp2(2,NBLOCKS),tmp3(2,NBLOCKS)
c      common    / D3dB_n_fft / tmp2,tmp3
c
c*     *** local variables ***
c      integer i,j,k,q,n,indx,ierr
c
cc     complex*16  tmp1(*)
cc     complex*16  tmp2(*)
cc     real*8      tmp3(*)
c      !integer tmp1(2),tmp2(2),tmp3(2)
c      logical value
c
c
c      call nwpw_timing_start(1)
c
c*     ***** allocate temporary space ****
c      !call D3dB_nfft3d(nb,nfft3d)
c
c
c*     ********************************************
c*     ***         Do a transpose of A          ***
c*     ***      A(kx,kz,ky) <- A(kx,ky,kz)      ***
c*     ********************************************
cc     call D3dB_c_transpose_jk(nb,A,dcpl_mb(tmp2(1)),dbl_mb(tmp3(1,nb)))
c
c*     *************************************************
c*     ***     do fft along kz dimension             ***
c*     ***   A(kx,nz(nb),ky) <- fft1d^(-1)[A(kx,kz,ky)]  ***
c*     *************************************************
c#ifdef MLIB
c      do n=1,ne
c      do q=1,nq(nb)
c      do i=1,(nx(nb)/2+1)
c         do k=1,nz(nb)
c            indx = i + (k-1)*(nx(nb)/2+1) + (q-1)*(nx(nb)/2+1)*nz(nb)
c     >               + (n-1)*nfft3d(nb)
c            dcpl_mb(tmp2(1,nb)+k-1) = A(indx)
c         end do
c         call z1dfft(dcpl_mb(tmp2(1,nb)),nz(nb),
c     >               dcpl_mb(tmpz(1,nb)),-2,ierr)
c         do k=1,nz(nb)
c            indx = i + (k-1)*(nx(nb)/2+1) + (q-1)*(nx(nb)/2+1)*nz(nb)
c     >               + (n-1)*nfft3d(nb)
c            A(indx) = dcpl_mb(tmp2(1,nb)+k-1)
c         end do
c      end do
c      end do
c      end do
c
c#else
c      do n=1,ne
c      do q=1,nq(nb)
c      do i=1,(nx(nb)/2+1)
c         do k=1,nz(nb)
c            indx = i + (k-1)*(nx(nb)/2+1) + (q-1)*(nx(nb)/2+1)*nz(nb)
c     >               + (n-1)*nfft3d(nb)
c            dcpl_mb(tmp2(1,nb)+k-1) = A(indx)
c         end do
c         call dcfftb(nz(nb),dcpl_mb(tmp2(1,nb)),dcpl_mb(tmpz(1,nb)))
c         do k=1,nz(nb)
c            indx = i + (k-1)*(nx(nb)/2+1) + (q-1)*(nx(nb)/2+1)*nz(nb)
c     >               + (n-1)*nfft3d(nb)
c            A(indx) = dcpl_mb(tmp2(1,nb)+k-1)
c         end do
c      end do
c      end do
c      end do
c#endif
c
c*     ********************************************
c*     ***         Do a transpose of A          ***
c*     ***      A(kx,ky,nz(nb)) <- A(kx,nz(nb),ky)      ***
c*     ********************************************
c      do n=1,ne
c        indx = 1 + (n-1)*nfft3d(nb)
c        call D3dB_c_transpose_jk(nb,A(indx),
c     >                           dcpl_mb(tmp2(1,nb)),
c     >                           dbl_mb(tmp3(1,nb)))
c      end do
cc     call D3dB_nc_transpose_jk(nb,ne,A,
cc    >                          dcpl_mb(tmp2(1,nb)),
cc    >                          dbl_mb(tmp3(1,nb)))
c
c*     *************************************************
c*     ***     do fft along ky dimension             ***
c*     ***   A(kx,ny(nb),nz(nb)) <- fft1d^(-1)[A(kx,ky,nz(nb))]  ***
c*     *************************************************
c#ifdef MLIB
c      do n=1,ne
c      do q=1,nq(nb)
c      do i=1,(nx(nb)/2+1)
c         do j=1,ny(nb)
c            indx = i + (j-1)*(nx(nb)/2+1) + (q-1)*(nx(nb)/2+1)*ny(nb)
c     >               + (n-1)*nfft3d(nb)
c            dcpl_mb(tmp2(1,nb)+j-1) = A(indx)
c         end do
c         call z1dfft(dcpl_mb(tmp2(1,nb)),ny(nb),
c     >               dcpl_mb(tmpy(1,nb)),-2,ierr)
c         do j=1,ny(nb)
c            indx = i + (j-1)*(nx(nb)/2+1) + (q-1)*(nx(nb)/2+1)*ny(nb)
c     >               + (n-1)*nfft3d(nb)
c            A(indx) = dcpl_mb(tmp2(1,nb)+j-1)
c         end do
c      end do
c      end do
c      end do
c      !call dscal((nx(nb)+2)*ny(nb)*nq(nb),dble(ny(nb)),A,1)
c#else
c      do n=1,ne
c      do q=1,nq(nb)
c      do i=1,(nx(nb)/2+1)
c         do j=1,ny(nb)
c            indx = i + (j-1)*(nx(nb)/2+1) + (q-1)*(nx(nb)/2+1)*ny(nb)
c     >               + (n-1)*nfft3d(nb)
c            dcpl_mb(tmp2(1,nb)+j-1) = A(indx)
c         end do
c         call dcfftb(ny(nb),dcpl_mb(tmp2(1,nb)),dcpl_mb(tmpy(1,nb)))
c         do j=1,ny(nb)
c            indx = i + (j-1)*(nx(nb)/2+1) + (q-1)*(nx(nb)/2+1)*ny(nb)
c     >               + (n-1)*nfft3d(nb)
c            A(indx) = dcpl_mb(tmp2(1,nb)+j-1)
c         end do
c      end do
c      end do
c      end do
c#endif
c
c*     *************************************************
c*     ***     do fft along kx dimension             ***
c*     ***   A(nx(nb),ny(nb),nz(nb)) <- fft1d^(-1)[A(kx,ny(nb),nz(nb))]  ***
c*     *************************************************
c#ifdef MLIB
c      do n=1,ne
c      do q=1,nq(nb)
c      do j=1,ny(nb)
c         indx = 1 + (j-1)*(nx(nb)/2+1) + (q-1)*(nx(nb)/2+1)*ny(nb)
c     >            + (n-1)*nfft3d(nb)
c         call drc1ft(A(indx),nx(nb),dcpl_mb(tmpx(1,nb)),-2,ierr)
c      end do
c      end do
c      end do
c
c#else
c
c      call cshift1_fftb(nx(nb),ny(nb),nq(nb),ne,A)
c      do n=1,ne
c      do q=1,nq(nb)
c      do j=1,ny(nb)
c         indx = 1 + (j-1)*(nx(nb)/2+1) + (q-1)*(nx(nb)/2+1)*ny(nb)
c     >            + (n-1)*nfft3d(nb)
c         call drfftb(nx(nb),A(indx),dcpl_mb(tmpx(1,nb)))
c      end do
c      end do
c      end do
c      call zeroend_fftb(nx(nb),ny(nb),nq(nb),ne,A)
c#endif
c    
c      call nwpw_timing_end(1)
c      return
c      end




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

#include "bafdecls.fh"
#include "errquit.fh"

#include "D3dB.fh"

      integer tmpx(2,NBLOCKS),tmpy(2,NBLOCKS),tmpz(2,NBLOCKS)
      common    / D3dB_fft / tmpx,tmpy,tmpz


*     *** local variables ***
      integer i,j,k,q,indx,indx1
      integer nxh,nxhy,nxhz

      !integer tmp1(2),tmp2(2),tmp3(2)
      integer tmp2(2),tmp3(2)
      logical value

      integer  tid,nthr,offset,nnn
      integer  Parallel_threadid,Parallel_nthreads
      external Parallel_threadid,Parallel_nthreads

      tid  = Parallel_threadid()
      nthr = Parallel_nthreads()

      call nwpw_timing_start(1)


*     ***** allocate temporary space ****
      !call D3dB_nfft3d(nb,nfft3d)
      nnn = nfft3d(nb)
      if ((nthr*ny(nb)).gt.nnn) nnn = nthr*ny(nb)
      if ((nthr*nz(nb)).gt.nnn) nnn = nthr*nz(nb)
      value = BA_push_get(mt_dcpl,nnn,'tmp2',tmp2(2),tmp2(1))
c      value = BA_push_get(mt_dcpl,(nfft3d(nb)),'tmp2',tmp2(2),tmp2(1))
      value = value.and.
     >        BA_push_get(mt_dbl,(n2ft3d(nb)),'tmp3',tmp3(2),tmp3(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)

      nxh = (nx(nb)/2+1)
      nxhz = nxh*nz(nb)
      nxhy = nxh*ny(nb)

      !**********************
      !**** slab mapping ****
      !**********************
      if (mapping.eq.1) then
*     ********************************************
*     ***     do fft along nx(nb) dimension        ***
*     ***   A(kx,ny(nb),nz(nb)) <- fft1d[A(nx(nb),ny(nb),nz(nb))]  ***
*     ********************************************
#ifdef MLIB
c     do q=1,nq(nb)
c     do j=1,ny(nb)
c        indx = 1 + (j-1)*(nx(nb)/2+1) + (q-1)*(nx(nb)/2+1)*ny(nb)
c        call drc1ft(A(indx),nx(nb),dcpl_mb(tmpx(1,nb)),1,ierr)
c     end do
c     end do
      call drcfts(A,nx(nb),1,ny(nb)*nq(nb),
     >                  nx(nb)+2,1,ierr)

#else

#ifdef FFTW3
      call dfftw_execute_dft_r2c(plans(4,nb),A,A)

#else
      offset=tid*(2*nx(nb)+15)
      do j=tid+1,ny(nb)*nq(nb),nthr
         call drfftf(nx(nb),A((j-1)*nxh+1),dcpl_mb(tmpx(1,nb)+offset))
      end do
      call cshift_fftf(nx(nb),ny(nb),nq(nb),1,A)
#endif
#endif


*     ********************************************
*     ***     do fft along ny(nb) dimension        ***
*     ***   A(kx,ky,nz(nb)) <- fft1d[A(kx,ny(nb),nz(nb))]  ***
*     ********************************************

#ifdef MLIB
      !call z1dfft(dbl_mb(tmp3(1)),ny(nb),dcpl_mb(tmp1(1)),-3,ierr)
      do q=1,nq(nb)
      do i=1,(nx(nb)/2+1)
         do j=1,ny(nb)
            indx = i + (j-1)*(nx(nb)/2+1) + (q-1)*(nx(nb)/2+1)*ny(nb)
            dcpl_mb(tmp2(1)+j-1) = A(indx)
         end do
         !indx = i + (q-1)*(nx(nb)/2+1)*ny(nb)
         !call zcopy(ny(nb),A(indx),(nx(nb)/2+1),dcpl_mb(tmp2(1)),1)
         call z1dfft(dcpl_mb(tmp2(1)),ny(nb),dcpl_mb(tmpy(1,nb)),1,ierr)
         do j=1,ny(nb)
            indx = i + (j-1)*(nx(nb)/2+1) + (q-1)*(nx(nb)/2+1)*ny(nb)
            A(indx) = dcpl_mb(tmp2(1)+j-1)
         end do
         !call zcopy(ny(nb),dcpl_mb(tmp2(1)),1,A(indx),(nx(nb)/2+1))
      end do
      end do
ccc   *** this should be faster but isn't  ***
c      do i=1,(nx(nb)/2+1)
c        !indx = 1 + (q-1)*(nx(nb)/2+1)*ny(nb)
c        call zffts(A(i),ny(nb),(nx(nb)/2+1),nq(nb),
c     >           (nx(nb)/2+1)*nq(nb),1,ierr)
c      end do
#else

#ifdef FFTW3
      do q=1,nq(nb)
         indx = 1+(q-1)*nxhy
         call dfftw_execute_dft(plans(5,nb),A(indx),A(indx))
      end do

#else
      offset=tid*(2*ny(nb)+15)
      do i=tid+1,nxh,nthr
      indx = i
      indx1= i
      do q=1,nq(nb)
         do j=1,ny(nb)
            !indx = i + (j-1)*(nx(nb)/2+1) + (q-1)*(nx(nb)/2+1)*ny(nb)
            dcpl_mb(tmp2(1)+j-1+tid*ny(nb)) = A(indx)
            indx = indx + nxh
         end do

         call dcfftf(ny(nb),dcpl_mb(tmp2(1)+tid*ny(nb)),
     >                      dcpl_mb(tmpy(1,nb)+offset))

         do j=1,ny(nb)
            !indx = i + (j-1)*(nx(nb)/2+1) + (q-1)*(nx(nb)/2+1)*ny(nb)
            A(indx1) = dcpl_mb(tmp2(1)+j-1+tid*ny(nb))
            indx1 = indx1 + nxh
         end do
      end do
      end do

#endif
#endif

*     ********************************************
*     ***         Do a transpose of A          ***
*     ***      A(ky,nz(nb),ky) <- A(kx,ky,nz(nb))      ***
*     ********************************************
      call D3dB_c_transpose_jk(nb,A,dcpl_mb(tmp2(1)),dbl_mb(tmp3(1)))


*     ********************************************
*     ***     do fft along nz(nb) dimension        ***
*     ***   A(kx,kz,ky) <- fft1d[A(kx,nz(nb),ky)]  ***
*     ********************************************
#ifdef MLIB
      !call z1dfft(dbl_mb(tmp3(1)),nz(nb),dcpl_mb(tmp1(1)),-3,ierr)
      do q=1,nq(nb)
      do i=1,(nx(nb)/2+1)
         do k=1,nz(nb)
            indx = i + (k-1)*(nx(nb)/2+1) + (q-1)*(nx(nb)/2+1)*ny(nb)
            dcpl_mb(tmp2(1)+k-1) = A(indx)
         end do
         call z1dfft(dcpl_mb(tmp2(1)),nz(nb),dcpl_mb(tmpz(1,nb)),1,ierr)
         do k=1,nz(nb)
            indx = i + (k-1)*(nx(nb)/2+1) + (q-1)*(nx(nb)/2+1)*ny(nb)
            A(indx) = dcpl_mb(tmp2(1)+k-1)
         end do
      end do
      end do
#else

#ifdef FFTW3
      do q=1,nq(nb)
         indx = 1+(q-1)*nxhz
         call dfftw_execute_dft(plans(6,nb),A(indx),A(indx))
      end do

#else
      offset=tid*(2*nz(nb)+15)
      do i=tid+1,nxh,nthr
      indx  = i
      indx1 = i
      do q=1,nq(nb)

         do k=1,nz(nb)
            !indx = i + (k-1)*(nx(nb)/2+1) + (q-1)*(nx(nb)/2+1)*ny(nb)
            dcpl_mb(tmp2(1)+k-1+tid*nz(nb)) = A(indx)
            indx = indx + nxh
         end do
         call dcfftf(nz(nb),dcpl_mb(tmp2(1)+tid*nz(nb)),
     >                      dcpl_mb(tmpz(1,nb)+offset))
         do k=1,nz(nb)
            !indx = i + (k-1)*(nx(nb)/2+1) + (q-1)*(nx(nb)/2+1)*ny(nb)
            A(indx1) = dcpl_mb(tmp2(1)+k-1+tid*nz(nb))
            indx1 = indx1 + nxh
         end do

      end do
      end do

#endif
#endif

*     ********************************************
*     ***         Do a transpose of A          ***
*     ***      A(kx,ky,kz) <- A(kx,kz,ky)      ***
*     ********************************************
c     call D3dB_c_transpose_jk(nb,A,dcpl_mb(tmp2(1)),dbl_mb(tmp3(1)))





      !*************************
      !**** hilbert mapping ****
      !*************************
      else

*     ********************************************
*     ***     do fft along nx(nb) dimension        ***
*     ***   A(kx,ny(nb),nz(nb)) <- fft1d[A(nx(nb),ny(nb),nz(nb))]  ***
*     ********************************************
#ifdef MLIB
      call drcfts(A,nx(nb),1,nq1(nb),
     >                  nx(nb)+2,1,ierr)
#else

#ifdef FFTW3
      call dfftw_execute_dft_r2c(plans(14,nb),A,A)
#else

      offset=tid*(2*nx(nb)+15)
      do i=tid+1,nq1(nb),nthr
        call drfftf(nx(nb),A(1+(i-1)*nxh),dcpl_mb(tmpx(1,nb)+offset))
      end do
      call cshift_fftf(nx(nb),nq1(nb),1,1,A)


#endif
#endif

      call D3dB_c_transpose_ijk(nb,1,A,dcpl_mb(tmp2(1)),dbl_mb(tmp3(1)))

*     ********************************************
*     ***     do fft along ny(nb) dimension        ***
*     ***   A(ky,nz(nb),kx) <- fft1d[A(ny(nb),nz(nb),kx)]  ***
*     ********************************************
#ifdef MLIB
      indx = 1
      do q=1,nq2(nb)
         !indx = 1 + (q-1)*ny(nb)
         call z1dfft(A(indx),ny(nb),dcpl_mb(tmpy(1,nb)),1,ierr)
         indx = indx + ny(nb)
      end do
#else

#ifdef FFTW3
      call dfftw_execute_dft(plans(15,nb),A,A)

#else
      offset=tid*(2*ny(nb)+15)
      do i=tid+1,nq2(nb),nthr
        call dcfftf(ny(nb),A(1+(i-1)*ny(nb)),dcpl_mb(tmpy(1,nb)+offset))
      end do
!$OMP BARRIER
#endif
#endif

      call D3dB_c_transpose_ijk(nb,2,A,dcpl_mb(tmp2(1)),dbl_mb(tmp3(1)))

*     ********************************************
*     ***     do fft along nz(nb) dimension        ***
*     ***   A(kz,kx,ky) <- fft1d[A(nz(nb),kx,ky)]  ***
*     ********************************************
#ifdef MLIB
      indx = 1
      do q=1,nq3(nb)
         !indx = 1 + (q-1)*nz(nb)
         call z1dfft(A(indx),nz(nb),dcpl_mb(tmpz(1,nb)),1,ierr)
         indx = indx + nz(nb)
      end do
#else

#ifdef FFTW3
      call dfftw_execute_dft(plans(16,nb),A,A)
#else
      offset=tid*(2*nz(nb)+15)
      do i=tid+1,nq3(nb),nthr
        call dcfftf(nz(nb),A(1+(i-1)*nz(nb)),dcpl_mb(tmpz(1,nb)+offset))
      end do
!$OMP BARRIER

#endif
#endif

      end if

*     **** deallocate temporary space  ****
      value = BA_pop_stack(tmp3(2))
      value = BA_pop_stack(tmp2(2))
      !value = BA_pop_stack(tmp1(2))

      call nwpw_timing_end(1)

      return
      end

      subroutine D3dB_fftfx_sub(n,nx,nxh,tmpx,A)
      implicit none
      integer n,nx,nxh
      real*8     tmpx(2*nx+15)
      complex*16 A(nxh,n)
      integer i



      do i=1,n
         call drfftf(nx,A(1,i),tmpx)
      end do


      return
      end

      subroutine  D3dB_fftfy_sub2(n,ny,tmpy,A)
      implicit none
      integer n,ny
      real*8     tmpy(4*ny+15)
      complex*16 A(ny,n)
      integer i,indx



      do i=1,n
         call dcfftf(ny,A(1,i),tmpy)
      end do


      return
      end

      subroutine  D3dB_fftfz_sub2(n,nz,tmpz,A)
      implicit none
      integer n,nz
      real*8     tmpz(4*nz+15)
      complex*16 A(nz,n)
      integer i



      do i=1,n
         call dcfftf(nz,A(1,i),tmpz)
      end do


      return
      end





      subroutine cshift_fftf(nx,ny,nq,ne,A)
      implicit none
      integer nx,ny,nq,ne
      real*8 A(*)

      integer i,j,indx

!$OMP BARRIER
!$OMP DO 
      do j=1,(ny*nq*ne)
        indx = 1+(j-1)*(nx+2)
c        indx = 1 + (j-1)*(nx+2) + (q-1)*(nx+2)*ny
c     >        + (n-1)*(nx+2)*ny*nq

         do i=nx,2,-1
            A(indx+i) = A(indx+i-1)
         end do
         A(indx+1)    = 0.0d0
         A(indx+nx+1) = 0.0d0
!         indx = indx + (nx+2)
      end do
!$OMP END DO

      return
      end 


      subroutine cshift_fftf_ab(nx,ny,nq,ne,A,B)
      implicit none
      integer nx,ny,nq,ne
      real*8 A(*)
      real*8 B(*)

      integer i,j,indx

      indx = 1
      do j=1,(ny*nq*ne)
CDIR$ NOVECTOR
         do i=nx,2,-1
            B(indx+i) = A(indx+i-1)
         end do
         B(indx+1)    = 0.0d0
         B(indx+nx+1) = 0.0d0
         indx = indx + (nx+2)
      end do

      return
      end




c*     ***********************************
c*     *					*
c*     *	        D3dB_nrc_fft3f		*
c*     *					*
c*     ***********************************
c
c      subroutine D3dB_nrc_fft3f(nb,ne,A)
c
c*****************************************************
c*                                                   *
c*      This routine performs the operation of       *
c*      a three dimensional complex to complex fft   *
c*           A(kx,ky,kz) <- FFT3[A(nx(nb),ny(nb),nz(nb))]        * 
c*                                                   *
c*      Entry - 					    *
c*              A: a column distribuded 3d block     *
c*              tmp: tempory work space must be at   *
c*                    least the size of (complex)    *
c*                    (nfft*nfft + 1) + 10*nfft      * 
c*                                                   *
c*       Exit - A is transformed                     *
c*                                                   *
c*       uses - transpose1 subroutine                *
c*                                                   *
c*****************************************************
c
c      implicit none
c      integer nb,ne
c      complex*16  A(*)
c
c#include "bafdecls.fh"
c#include "errquit.fh"
c
c#include "D3dB.fh"
c
c      integer tmpx(2,NBLOCKS),tmpy(2,NBLOCKS),tmpz(2,NBLOCKS)
c      common    / D3dB_fft / tmpx,tmpy,tmpz
c
c*     *** local variables ***
c      integer i,j,k,q,n,indx,ierr
c
c      !integer tmp1(2),tmp2(2),tmp3(2)
c      integer tmp2(2),tmp3(2)
c      logical value
c
c
c      call nwpw_timing_start(1)
c
c
c*     ***** allocate temporary space ****
c      !call D3dB_nfft3d(nb,nfft3d)
c      value = BA_push_get(mt_dcpl,(nfft3d(nb)),'tmp2',tmp2(2),tmp2(1))
c      value = value.and.
c     >        BA_push_get(mt_dbl,(n2ft3d(nb)),'tmp3',tmp3(2),tmp3(1))
c      if (.not. value) call errquit('out of stack memory',0, MA_ERR)
c
c
c*     ********************************************
c*     ***     do fft along nx(nb) dimension        ***
c*     ***   A(kx,ny(nb),nz(nb)) <- fft1d[A(nx(nb),ny(nb),nz(nb))]  ***
c*     ********************************************
c#ifdef MLIB
c      call drcfts(A,nx(nb),1,ny(nb)*nq(nb)*ne,
c     >                  nx(nb)+2,1,ierr)
c
c#else
c      !call drffti(nx(nb),dcpl_mb(tmp1(1)))
c      do n=1,ne
c      do q=1,nq(nb)
c      do j=1,ny(nb)
c         indx = 1 + (j-1)*(nx(nb)/2+1) + (q-1)*(nx(nb)/2+1)*ny(nb)
c     >        + (n-1)*nfft3d(nb)
c
c
c         call drfftf(nx(nb),A(indx),dcpl_mb(tmpx(1,nb)))
c      end do
c      end do
c      end do
c      call cshift_fftf(nx(nb),ny(nb),nq(nb),ne,A)
c#endif
c
c
c*     ********************************************
c*     ***     do fft along ny(nb) dimension        ***
c*     ***   A(kx,ky,nz(nb)) <- fft1d[A(kx,ny(nb),nz(nb))]  ***
c*     ********************************************
c
c#ifdef MLIB
c      !call z1dfft(dbl_mb(tmp3(1)),ny(nb),dcpl_mb(tmp1(1)),-3,ierr)
c      do n=1,ne
c      do q=1,nq(nb)
c      do i=1,(nx(nb)/2+1)
c         do j=1,ny(nb)
c            indx = i + (j-1)*(nx(nb)/2+1) + (q-1)*(nx(nb)/2+1)*ny(nb)
c     >               + (n-1)*nfft3d(nb)
c            dcpl_mb(tmp2(1)+j-1) = A(indx)
c         end do
c         !indx = i + (q-1)*(nx(nb)/2+1)*ny(nb)
c         !call zcopy(ny(nb),A(indx),(nx(nb)/2+1),dcpl_mb(tmp2(1)),1)
c         call z1dfft(dcpl_mb(tmp2(1)),ny(nb),dcpl_mb(tmpy(1,nb)),1,ierr)
c         do j=1,ny(nb)
c            indx = i + (j-1)*(nx(nb)/2+1) + (q-1)*(nx(nb)/2+1)*ny(nb)
c     >               + (n-1)*nfft3d(nb)
c            A(indx) = dcpl_mb(tmp2(1)+j-1)
c         end do
c         !call zcopy(ny(nb),dcpl_mb(tmp2(1)),1,A(indx),(nx(nb)/2+1))
c      end do
c      end do
c      end do
cccc   *** this should be faster but isn't  ***
cc      do i=1,(nx(nb)/2+1)
cc        !indx = 1 + (q-1)*(nx(nb)/2+1)*ny(nb)
cc        call zffts(A(i),ny(nb),(nx(nb)/2+1),nq(nb),
cc     >           (nx(nb)/2+1)*nq(nb),1,ierr)
cc      end do
c#else
c      !call dcffti(ny(nb),dcpl_mb(tmp1(1)))
c      do n=1,ne
c      do q=1,nq(nb)
c      do i=1,(nx(nb)/2+1)
c         do j=1,ny(nb)
c            indx = i + (j-1)*(nx(nb)/2+1) + (q-1)*(nx(nb)/2+1)*ny(nb)
c     >               + (n-1)*nfft3d(nb)
c            dcpl_mb(tmp2(1)+j-1) = A(indx)
c         end do
c         call dcfftf(ny(nb),dcpl_mb(tmp2(1)),dcpl_mb(tmpy(1,nb)))
c         do j=1,ny(nb)
c            indx = i + (j-1)*(nx(nb)/2+1) + (q-1)*(nx(nb)/2+1)*ny(nb)
c     >               + (n-1)*nfft3d(nb)
c            A(indx) = dcpl_mb(tmp2(1)+j-1)
c         end do
c      end do
c      end do
c      end do
c#endif
c
c
c*     ********************************************
c*     ***         Do a transpose of A          ***
c*     ***      A(ky,nz(nb),ky) <- A(kx,ky,nz(nb))      ***
c*     ********************************************
c      do n=1,ne
c        indx = 1 + (n-1)*nfft3d(nb)
c        call D3dB_c_transpose_jk(nb,A(indx),
c     >                           dcpl_mb(tmp2(1)),
c     >                           dbl_mb(tmp3(1)))
c      end do
c
c
c*     ********************************************
c*     ***     do fft along nz(nb) dimension        ***
c*     ***   A(kx,kz,ky) <- fft1d[A(kx,nz(nb),ky)]  ***
c*     ********************************************
c#ifdef MLIB
c      !call z1dfft(dbl_mb(tmp3(1)),nz(nb),dcpl_mb(tmp1(1)),-3,ierr)
c      do n=1,ne
c      do q=1,nq(nb)
c      do i=1,(nx(nb)/2+1)
c         do k=1,nz(nb)
c            indx = i + (k-1)*(nx(nb)/2+1) + (q-1)*(nx(nb)/2+1)*ny(nb)
c     >               + (n-1)*nfft3d(nb)
c            dcpl_mb(tmp2(1)+k-1) = A(indx)
c         end do
c         call z1dfft(dcpl_mb(tmp2(1)),nz(nb),dcpl_mb(tmpz(1,nb)),1,ierr)
c         do k=1,nz(nb)
c            indx = i + (k-1)*(nx(nb)/2+1) + (q-1)*(nx(nb)/2+1)*ny(nb)
c     >               + (n-1)*nfft3d(nb)
c            A(indx) = dcpl_mb(tmp2(1)+k-1)
c         end do
c      end do
c      end do
c      end do
c#else
c      !call dcffti(nz(nb),dcpl_mb(tmp1(1)))
c      do n=1,ne
c      do q=1,nq(nb)
c      do i=1,(nx(nb)/2+1)
c         do k=1,nz(nb)
c            indx = i + (k-1)*(nx(nb)/2+1) + (q-1)*(nx(nb)/2+1)*ny(nb)
c     >               + (n-1)*nfft3d(nb)
c            dcpl_mb(tmp2(1)+k-1) = A(indx)
c         end do
c         call dcfftf(nz(nb),dcpl_mb(tmp2(1)),dcpl_mb(tmpz(1,nb)))
c         do k=1,nz(nb)
c            indx = i + (k-1)*(nx(nb)/2+1) + (q-1)*(nx(nb)/2+1)*ny(nb)
c     >               + (n-1)*nfft3d(nb)
c            A(indx) = dcpl_mb(tmp2(1)+k-1)
c         end do
c      end do
c      end do
c      end do
c#endif
c
c*     ********************************************
c*     ***         Do a transpose of A          ***
c*     ***      A(kx,ky,kz) <- A(kx,kz,ky)      ***
c*     ********************************************
cc     call D3dB_c_transpose_jk(nb,A,dcpl_mb(tmp2(1)),dbl_mb(tmp3(1)))
c
c
c*     **** deallocate temporary space  ****
c      value = BA_pop_stack(tmp3(2))
c      value = BA_pop_stack(tmp2(2))
c      !value = BA_pop_stack(tmp1(2))

c      call nwpw_timing_end(1)
c      return
c      end




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

#include "D3dB.fh"

      integer i

      do i=1,nfft3d_map(nb)
         C(i) = scale*A(i)
      end do
      return
      end


      subroutine D3dB_c_SMul1(nb,scale,A)
      implicit none 
      integer    nb
      real*8     scale
      complex*16 A(*)
      
#include "D3dB.fh"

      integer i
      
      do i=1,nfft3d_map(nb)
         A(i) = scale*A(i)
      end do 
      return 
      end



      subroutine D3dB_r_SMul(nb,scale,A,C)     
      implicit none 
      integer nb
      real*8     scale
      real*8 A(*)
      real*8 C(*)

#include "D3dB.fh"

      integer i

!$OMP DO
      do i=1,n2ft3d_map(nb)
         C(i) = scale*A(i)
      end do
!$OMP END DO
      return
      end 

      subroutine D3dB_r_SMul1(nb,scale,A)
      implicit none
      integer nb
      real*8     scale
      real*8 A(*)

#include "D3dB.fh"

      integer i

!$OMP DO
      do i=1,n2ft3d_map(nb)
         A(i) = scale*A(i)
      end do
!$OMP END DO
      return
      end


      subroutine D3dB_t_SMul(nb,scale,A,C)     
      implicit none 
      integer nb
      real*8 scale
      real*8 A(*)
      real*8 C(*)

#include "D3dB.fh"

      integer i

      do i=1,nfft3d_map(nb)
         C(i) = scale*A(i)
      end do
      return
      end 


      subroutine D3dB_t_SMul1(nb,scale,A)
      implicit none
      integer nb
      real*8 scale
      real*8 A(*)

#include "D3dB.fh"

      integer i

      do i=1,nfft3d_map(nb)
         A(i) = scale*A(i)
      end do
      return
      end



      subroutine D3dB_c_ZMul(nb,scale,A,C)     
      implicit none 
      integer    nb
      complex*16 scale
      complex*16 A(*)
      complex*16 C(*)

#include "D3dB.fh"


      integer i

      do i=1,nfft3d_map(nb)
         C(i) = scale*A(i)
      end do
      return
      end

      subroutine D3dB_r_Power1(nb,y,A)
      implicit none
      integer nb
      real*8  y
      real*8 A(*)

#include "D3dB.fh"

      integer ii

      do ii=1,n2ft3d_map(nb)
         A(ii) = A(ii)**y
      end do
      
      return
      end


      subroutine D3dB_rr_Power(nb,y,A,B)
      implicit none
      integer nb
      real*8  y
      real*8 A(*)
      real*8 B(*)

#include "D3dB.fh"

      integer i

      do i=1,n2ft3d_map(nb)
         B(i) = A(i)**y
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

#include "D3dB.fh"

      integer i

      do i=1,nfft3d_map(nb)
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

#include "D3dB.fh"


      integer i

!$OMP DO
      do i=1,n2ft3d_map(nb)
         C(i) = A(i)**2
      end do
!$OMP END DO
      return
      end


*     ***********************************
*     *                                 *
*     *        D3dB_rr_SqrAdd           *
*     *                                 *
*     ***********************************

      subroutine D3dB_rr_SqrAdd(nb,A,C)
      implicit none
      integer nb
      real*8 A(*)
      real*8 C(*)

#include "D3dB.fh"


      integer i

!$OMP DO
      do i=1,n2ft3d_map(nb)
         C(i) = C(i) + A(i)**2
      end do
!$OMP END DO
      return
      end




*     ***********************************
*     *                                 *
*     *        D3dB_rr_Sqr1             *
*     *                                 *
*     ***********************************

      subroutine D3dB_rr_Sqr1(nb,A)
      implicit none
      integer nb
      real*8 A(*)

#include "D3dB.fh"

      integer i

!$OMP DO
      do i=1,n2ft3d_map(nb)
         A(i) = A(i)**2
      end do
!$OMP END DO
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

#include "D3dB.fh"

      integer i

      do i=1,n2ft3d_map(nb)
         C(i) = dsqrt(A(i))
      end do
      return
      end



*     ***********************************
*     *                                 *
*     *        D3dB_rr_Sqrt1            *
*     *                                 *
*     ***********************************

      subroutine D3dB_rr_Sqrt1(nb,A)
      implicit none
      integer nb
      real*8 A(*)

#include "D3dB.fh"

      integer i

!$OMP DO
      do i=1,n2ft3d_map(nb)
         A(i) = dsqrt(A(i))
      end do
!$OMP END DO
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

#include "D3dB.fh"

      integer i

      do i=1,nfft3d_map(nb)
         C(i) = A(i)**2
      end do
      return
      end

             

*     ***********************************
*     *					*
*     *	   D3dB_c_transpose_jk_init	*
*     *					*
*     ***********************************

      subroutine D3dB_c_transpose_jk_init(nb)
      implicit none
      integer nb

#include "bafdecls.fh"
#include "errquit.fh"
#include "D3dB.fh"


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
      value = BA_alloc_get(mt_int,((nx(nb)/2+1)*ny(nb)*nq(nb)),
     >                     'iq_to_i1',iq_to_i1(2,nb),iq_to_i1(1,nb))
      value=value.and.BA_alloc_get(mt_int,((nx(nb)/2+1)*ny(nb)*nq(nb)),
     >                     'iq_to_i2',iq_to_i2(2,nb),iq_to_i2(1,nb))
      value = value.and.BA_alloc_get(mt_int,(nz(nb)+1),
     >                     'i1_start',i1_start(2,nb),i1_start(1,nb))
      value = value.and.BA_alloc_get(mt_int,(nz(nb)+1),
     >                     'i2_start',i2_start(2,nb),i2_start(1,nb))
      if (.not. value) 
     > call errquit('D3dB_transpose_jk_init:out of heap',0,MA_ERR)

      call Parallel2d_taskid_i(taskid)
      call Parallel2d_np_i(np)

!MATHIAS
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
*     *	         D3dB_cc_dot  	 	*	
*     *					*
*     ***********************************

      subroutine D3dB_cc_dot(nb,A,B,sumall)     
      implicit none
      integer nb 
      complex*16 A(*)
      complex*16 B(*)
      real*8     sumall


#include "D3dB.fh"

      integer i,j,k,q,index,np,taskid,p
      real*8  sum


      call nwpw_timing_start(2)

      call Parallel2d_np_i(np)

*     **** sum up dot product on this node ****
      sum = 0.0d0

      !**********************
      !**** slab mapping ****
      !**********************
      if (mapping.eq.1) then
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

      !*************************
      !**** hilbert mapping ****
      !*************************
      else
      call Parallel2d_taskid_i(taskid)
*     ***** kx!=0 plane, so double count *****
      do index=1,nfft3d_map(nb)
            sum = sum + dble(A(index))  * dble(B(index))
     >                + dimag(A(index)) * dimag(B(index))
      end do
      sum = sum*2.0d0

*     ***** kx==0 plane, so single count *****
      do k=1,nz(nb)
      do j=1,ny(nb)
         i=1
         call D3dB_ijktoindexp(1,i,j,k,index,p)
         if (p.eq.taskid) then
         sum = sum - dble(A(index))  * dble(B(index))
     >             - dimag(A(index)) * dimag(B(index))
         end if
      end do
      end do

      end if
      

*     **** add up sums from other nodes ****
      if (np.gt.1) then
         call D3dB_SumAll(sum)
      end if

      call nwpw_timing_end(2)

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


#include "D3dB.fh"

      integer i,j,k,q,index,np,taskid,p
      real*8  sum


      call nwpw_timing_start(2)

c      call Parallel2d_np_i(np)

*     **** sum up dot product on this node ****
      sum = 0.0d0

      !**********************
      !**** slab mapping ****
      !**********************
      if (mapping.eq.1) then
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
      

      !*************************
      !**** hilbert mapping ****
      !*************************
      else
      call Parallel2d_taskid_i(taskid)
*     ***** kx!=0 plane, so double count *****
      do index=1,nfft3d_map(nb)
            sum = sum + dble(A(index))  * dble(B(index))
     >                + dimag(A(index)) * dimag(B(index))
      end do
      sum = sum*2.0d0

*     ***** kx==0 plane, so single count *****
      do k=1,nz(nb)
      do j=1,ny(nb)
         i=1
         call D3dB_ijktoindexp(1,i,j,k,index,p)
         if (p.eq.taskid) then
         sum = sum - dble(A(index))  * dble(B(index))
     >             - dimag(A(index)) * dimag(B(index))
         end if
      end do
      end do

      end if


*     **** do not add up sums from other nodes ****
       
      call nwpw_timing_end(2)

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

#include "D3dB.fh"

      integer i,j,k,q,index,np,nxh,taskid,p
      real*8  sum

      nxh=nx(nb)/2
      call Parallel2d_np_i(np)

*     **** sum up dot product on this node ****
      sum = 0.0d0

      !**********************
      !**** slab mapping ****
      !**********************
      if (mapping.eq.1) then
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

      
      !*************************
      !**** hilbert mapping ****
      !*************************
      else
      call Parallel2d_taskid_i(taskid)
*     ***** kx!=0 plane, so double count *****
      do index=1,nfft3d_map(nb)
            sum = sum + A(index)*B(index)
      end do
      sum = sum*2.0d0

*     ***** kx==0 plane, so single count *****
      do k=1,nz(nb)
      do j=1,ny(nb)
         i=1
         call D3dB_ijktoindexp(1,i,j,k,index,p)
         if (p.eq.taskid) then
         sum = sum - A(index)*B(index)
         end if
      end do
      end do

      end if
      

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


#include "D3dB.fh"

      integer i,j,k,q,index,np,nxh,taskid,p
      real*8  sum


      call nwpw_timing_start(2)

      nxh=nx(nb)/2
c      call Parallel2d_np_i(np)

*     **** sum up dot product on this node ****
      sum = 0.0d0

      !**********************
      !**** slab mapping ****
      !**********************
      if (mapping.eq.1) then
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
      


      !*************************
      !**** hilbert mapping ****
      !*************************
      else
      call Parallel2d_taskid_i(taskid)
*     ***** kx!=0 plane, so double count *****
      do index=1,nfft3d_map(nb)
            sum = sum + A(index)*B(index)
      end do
      sum = sum*2.0d0

*     ***** kx==0 plane, so single count *****
      do k=1,nz(nb)
      do j=1,ny(nb)
         i=1
         call D3dB_ijktoindexp(1,i,j,k,index,p)
         if (p.eq.taskid) then
         sum = sum - A(index)*B(index)
         end if
      end do
      end do

      end if



*     **** !!!! do not add up sums from other nodes ****

      call nwpw_timing_end(2)

      sumall = sum
      return
      end



*     ***********************************
*     *					*
*     *	         D3dB_rr_dot  	 	*	
*     *					*
*     ***********************************
*     shared memory output
*     - sumall

      subroutine D3dB_rr_dot(nb,A,B,sumall)     
      implicit none 
      integer nb
      real*8  A(*)
      real*8  B(*)
      real*8  sumall

#include "D3dB.fh"

      integer i,np
      real*8  sum
      common /D3dB_RR_TSUM/ sum

      call Parallel2d_np_i(np)

*     **** sum up dot product on this node ****
!$OMP MASTER
      sum = 0.0d0
!$OMP END MASTER
!$OMP BARRIER
!$OMP DO REDUCTION(+:sum)
      do i=1,n2ft3d_map(nb)
         sum = sum + A(i)*B(i)
      end do
!$OMP END DO


*     **** add up sums from other nodes ****
      if (np.gt.1) then
         call D3dB_SumAll(sum)
      end if

!$OMP MASTER
      sumall = sum
!$OMP END MASTER
!$OMP BARRIER
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

#include "D3dB.fh"

      integer i,np
      real*8  sum
      common /D3dB_RR_TSUM/ sum


*     **** sum up dot product on this node ****
!$OMP MASTER
      sum = 0.0d0
!$OMP END MASTER
!$OMP BARRIER
!$OMP DO reduction(+:sum)
      do i=1,n2ft3d_map(nb)
         sum = sum + A(i)*B(i)
      end do
!$OMP END DO

*     **** add up sums from other nodes ****
*     if (np.gt.1) then
*        call D3dB_SumAll(sum)
*     end if

!$OMP MASTER
      sumall = sum
!$OMP END MASTER
!$OMP BARRIER
      return
      end

*     ***********************************
*     *                                 *
*     *          D3dB_rrm_sym_dot       *       
*     *                                 *
*     ***********************************

      subroutine D3dB_rrm_sym_dot(nb,n,A,B,matrix)
      implicit none
      integer nb,n
      real*8  A(*)
      real*8  B(*)
      real*8  matrix(n,n)

#include "D3dB.fh"

*     **** local variables ****
      integer j,k
      integer np

      call nwpw_timing_start(2)
      call Parallel2d_np_i(np)

      do k=1,n
        call DGEMM_OMP('T','N',k,1,n2ft3d(nb),
     >             1.0d0,
     >             A,n2ft3d(nb),
     >             B(1+(k-1)*n2ft3d(nb)),n2ft3d(nb),
     >             0.0d0,
     >             matrix(1,k),k)
      end do

!$OMP DO
      do k=1,n
      do j=k+1,n
        matrix(j,k) = matrix(k,j)
      end do
      end do
!$OMP END DO

      if (np.gt.1) call D3dB_Vector_SumAll(n*n,matrix)
      call nwpw_timing_end(2)
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

#include "D3dB.fh"

      integer i

      do i=1,nfft3d_map(nb)
            C(i) = dconjg(A(i)) * B(i)
         end do

      return
      end

*     ***********************************
*     *                                 *
*     *          D3dB_cc_Mul2           *       
*     *                                 *
*     ***********************************

      subroutine D3dB_cc_Mul2(nb,A,B)
      implicit none
      integer    nb
      complex*16 A(*)
      complex*16 B(*)

#include "D3dB.fh"

      integer i

      do i=1,nfft3d_map(nb)
         B(i) = A(i) * B(i)
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

#include "D3dB.fh"

      integer i

      do i=1,nfft3d_map(nb)
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

#include "D3dB.fh"

      integer i

      do i=1,nfft3d_map(nb)
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

#include "D3dB.fh"

      integer i

      do i=1,nfft3d_map(nb)
            C(i) = A(i) * B(i)
      end do

      return
      end



*     ***********************************
*     *                                 *
*     *          D3dB_tc_Mul2           *
*     *                                 *
*     ***********************************

      subroutine D3dB_tc_Mul2(nb,A,B)
      implicit none
      integer    nb
      real*8     A(*)
      complex*16 B(*)

#include "D3dB.fh"

      integer i

!$OMP DO
      do i=1,nfft3d_map(nb)
            B(i) = B(i) * A(i)
      end do
!$OMP END DO

      return
      end

*     ***********************************
*     *					*
*     *	         D3dB_rr_Mul  	 	*	
*     *					*
*     ***********************************

      subroutine D3dB_rr_Mul(nb,A,B,C)     
      implicit none 

#include "D3dB.fh"

      integer nb
      real*8 A(*)
      real*8 B(*)
      real*8 C(*)
      integer i,n

!$OMP DO
      do i=1,n2ft3d_map(nb)
         C(i) = A(i) * B(i)
      end do
!$OMP END DO

      return
      end

*     ***********************************
*     *                                 *
*     *          D3dB_rr_Mul2           *
*     *                                 *
*     ***********************************

      subroutine D3dB_rr_Mul2(nb,A,B)
      implicit none
      integer nb
      real*8 A(*)
      real*8 B(*)

#include "D3dB.fh"

      integer i

!$OMP DO
      do i=1,n2ft3d_map(nb)
         B(i) = B(i) * A(i)
      end do
!$OMP END DO

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

#include "D3dB.fh"

      integer i

!$OMP DO
      do i=1,nfft3d_map(nb)
         C(i) = A(i) + B(i)
      end do
!$OMP END DO

      return
      end


*     ***********************************
*     *                                 *
*     *          D3dB_cc_Sum2           *
*     *                                 *
*     ***********************************

      subroutine D3dB_cc_Sum2(nb,A,B)
      implicit none
      integer    nb
      complex*16 A(*)
      complex*16 B(*)

#include "D3dB.fh"

      integer i

!$OMP DO
      do i=1,nfft3d_map(nb)
         B(i) = B(i) + A(i)
      end do
!$OMP END DO

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

#include "D3dB.fh"

      integer i

!$OMP DO
      do i=1,n2ft3d_map(nb)
         C(i) = B(i)+A(i)
      end do
!$OMP END DO

      return
      end


*     ***********************************
*     *                                 *
*     *          D3dB_rr_Sum2           *
*     *                                 *
*     ***********************************
      subroutine D3dB_rr_Sum2(nb,A,B)
      implicit none
      integer nb
      real*8  A(*)
      real*8  B(*)

#include "D3dB.fh"

      integer i

!$OMP DO
      do i=1,n2ft3d_map(nb)
         B(i) = B(i) + A(i)
      end do
!$OMP END DO
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

#include "D3dB.fh"

      integer i

      do i=1,nfft3d_map(nb)
         C(i) = A(i) - B(i)
      end do

      return
      end


*     ***********************************
*     *					*
*     *	         D3dB_rr_Sub   		*	
*     *					*
*     ***********************************

      subroutine D3dB_rr_Sub(nb,A,B,C)     
      implicit none 

#include "D3dB.fh"

      integer nb
      real*8  A(n2ft3d_map(nb))
      real*8  B(n2ft3d_map(nb))
      real*8  C(n2ft3d_map(nb))
      integer i


!$OMP DO
      do i=1,n2ft3d_map(nb)
         C(i) = A(i) - B(i)
      end do
!$OMP END DO

      return
      end



*     ***********************************
*     *                                 *
*     *          D3dB_rr_Sub2           *
*     *                                 *
*     ***********************************

      subroutine D3dB_rr_Sub2(nb,A,B)
      implicit none

#include "D3dB.fh"

      integer nb
      real*8  A(n2ft3d_map(nb))
      real*8  B(n2ft3d_map(nb))

      integer i

!$OMP DO
      do i=1,n2ft3d_map(nb)
         B(i) = B(i) - A(i)
      end do
!$OMP END DO


      return
      end

*     ***********************************
*     *                                 *
*     *          D3dB_rr_Multiply2      *       
*     *                                 *
*     ***********************************
      subroutine D3dB_rr_Multiply2(nb,A,B)
      implicit none 
      integer nb
      real*8  A(*)
      real*8  B(*)

#include "D3dB.fh"

      integer i

!$OMP DO
      do i=1,n2ft3d_map(nb)
         B(i) = B(i)*A(i)
      end do
!$OMP END DO

      return
      end



*     ***********************************
*     *                                 *
*     *       D3dB_rrr_MultiplyAdd      *       
*     *                                 *
*     ***********************************
      subroutine D3dB_rrr_MultiplyAdd(nb,A,B,C)
      implicit none
      integer nb
      real*8  A(*)
      real*8  B(*)
      real*8  C(*)

#include "D3dB.fh"

      integer i

!$OMP DO
      do i=1,n2ft3d_map(nb)
         C(i) = C(i) + B(i)*A(i)
      end do
!$OMP END DO

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

#include "D3dB.fh"

      integer i

!$OMP DO
      do i=1,nfft3d_map(nb)
         B(i) = B(i) + alpha*A(i)
      end do
!$OMP END DO

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

#include "D3dB.fh"

      integer i

!$OMP DO
      do i=1,nfft3d_map(nb)
         B(i) = B(i) + alpha*A(i)
      end do
!$OMP END DO

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

#include "D3dB.fh"

      integer i

!$OMP DO
      do i=1,n2ft3d_map(nb)
         B(i) = B(i) + alpha* A(i)
      end do
!$OMP END DO

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

#include "D3dB.fh"

      real*8 eta
      parameter (eta=1.0d-9)

      integer index

      !do q=1,nq(nb)
      !do j=1,ny(nb)
      !do i=1,nx(nb)
CDIR$ NOVECTOR
!$OMP DO
      do index = 1,n2ft3d_map(nb)
         !index = i + (j-1)*(nx(nb)+2) + (q-1)*(nx(nb)+2)*ny(nb)
         if (dabs(B(index)) .le. eta) then
           C(index) = 0.0d0
         else
           C(index) = A(index) / B(index)
         end if
      end do
!$OMP END DO
      !end do
      !end do
      !end do

      return
      end



*     ***********************************
*     *                                 *
*     *          D3dB_rr_Divide2         *
*     *                                 *
*     ***********************************

      subroutine D3dB_rr_Divide2(nb,A,B)
      implicit none
      integer nb
      real*8 A(*)
      real*8 B(*)

#include "D3dB.fh"

      real*8 eta
      parameter (eta=1.0d-9)

      integer index

      !do q=1,nq(nb)
      !do j=1,ny(nb)
      !do i=1,nx(nb)
CDIR$ NOVECTOR
!$OMP DO
      do index = 1,n2ft3d_map(nb)
         !index = i + (j-1)*(nx(nb)+2) + (q-1)*(nx(nb)+2)*ny(nb)
         if (dabs(A(index)) .le. eta) then
           B(index) = 0.0d0
         else
           B(index) = B(index) / A(index)
         end if
      end do
!$OMP END DO
      !end do
      !end do
      !end do

      return
      end



*     ***********************************
*     *                                 *
*     *          D3dB_r_ABS             *
*     *                                 *
*     ***********************************

      subroutine D3dB_r_ABS(nb,A,C)
      implicit none
      integer nb
      real*8 A(*)
      real*8 C(*)

#include "D3dB.fh"


      integer index

      !do q=1,nq(nb)
      !do j=1,ny(nb)
      !do i=1,nx(nb)
!$OMP DO
      do index=1,n2ft3d_map(nb)
         !index = i + (j-1)*(nx(nb)+2) + (q-1)*(nx(nb)+2)*ny(nb)
         C(index) = dabs(A(index))
      end do
!$OMP END DO
      !end do
      !end do
      !end do

      return
      end

      subroutine D3dB_r_abs1(nb,A)
      implicit none
      integer nb
      real*8 A(*)

#include "D3dB.fh"

      integer index

!$OMP DO
      do index=1,n2ft3d_map(nb)
         !index = i + (j-1)*(nx(nb)+2) + (q-1)*(nx(nb)+2)*ny(nb)
         A(index) = dabs(A(index))
      end do
!$OMP END DO
      return
      end

*     ***********************************
*     *                                 *
*     *          D3dB_r_ABSMAX          *
*     *                                 *
*     ***********************************

      subroutine D3dB_r_ABSMAX(nb,A,amax)
      implicit none
      integer nb
      real*8 A(*)
      real*8 amax

#include "D3dB.fh"

      integer index
      real*8 aa

      amax = 0.0d0
      do index=1,n2ft3d_map(nb)
         !index = i + (j-1)*(nx(nb)+2) + (q-1)*(nx(nb)+2)*ny(nb)
         aa = dabs(A(index))
         if (aa.gt.amax) amax = aa
      end do
      call D3dB_MaxAll(amax)

      return
      end


*     ***********************************
*     *                                 *
*     *          D3dB_r_ZeroNegative    *
*     *                                 *
*     ***********************************

      subroutine D3dB_r_ZeroNegative(nb,A)
      implicit none
      integer nb
      real*8 A(*)

#include "D3dB.fh"

      integer i

      do i=1,n2ft3d_map(nb)
         if (A(i).lt.0.0d0) A(i) = 0.0d0
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

#include "D3dB.fh"

      integer i

!$OMP DO
      do i=1,n2ft3d_map(nb)
         C(i) = A(i) - B(i)
      end do
!$OMP END DO

      return
      end

*     ***********************************
*     *					*
*     *          D3dB_r_Zero_Ends0	*
*     *					*
*     ***********************************

      subroutine D3dB_r_Zero_Ends0(nb,A)
      integer nb
      real*8 A(*) 

#include "D3dB.fh"

      integer j,k,q,index,taskid,p

      !**** slab mapping ****
      if (mapping.eq.1) then
      do q=1,nq(nb)
         do j=1,ny(nb)
            index = (nx(nb)+1) + (j-1)*(nx(nb)+2) 
     >                         + (q-1)*(nx(nb)+2)*(ny(nb))

            A(index)   = 0.0d0
            A(index+1) = 0.0d0
         end do
      end do


      !**** hilbert mapping ****
      else
        call Parallel2d_taskid_i(taskid)
        do k=1,nz(nb)
        do j=1,ny(nb)

         call D3dB_ijktoindex2p(nb,(nx(nb)+1),j,k,index,p)
         if (p.eq.taskid) A(index) = 0.0d0

         call D3dB_ijktoindex2p(nb,(nx(nb)+2),j,k,index,p)
         if (p.eq.taskid) A(index) = 0.0d0

        end do
        end do
      end if

      if (n2ft3d_map(nb).lt.n2ft3d(nb)) then
         call dcopy((n2ft3d(nb)-n2ft3d_map(nb)),
     >              0.0d0,0,A(n2ft3d_map(nb)+1),1)
      end if

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

#include "D3dB.fh"

      integer j,k,q,index,taskid,p

      !**** slab mapping ****
      if (mapping.eq.1) then
!$OMP DO
      do q=1,nq(nb)
         do j=1,ny(nb)
            index = (nx(nb)+1) + (j-1)*(nx(nb)+2) 
     >                         + (q-1)*(nx(nb)+2)*(ny(nb))

            A(index)   = 0.0d0
            A(index+1) = 0.0d0
         end do
      end do
!$OMP END DO


      !**** hilbert mapping ****
      else
        call Parallel2d_taskid_i(taskid)
!$OMP DO
        do k=1,nz(nb)
        do j=1,ny(nb)

         call D3dB_ijktoindex2p(nb,(nx(nb)+1),j,k,index,p)
         if (p.eq.taskid) A(index) = 0.0d0

         call D3dB_ijktoindex2p(nb,(nx(nb)+2),j,k,index,p)
         if (p.eq.taskid) A(index) = 0.0d0

        end do
        end do
!$OMP end DO
      end if

!$OMP MASTER
      if (n2ft3d_map(nb).lt.n2ft3d(nb)) then
         call dcopy((n2ft3d(nb)-n2ft3d_map(nb)),
     >              0.0d0,0,A(n2ft3d_map(nb)+1),1)
      end if
!$OMP END MASTER
!$OMP BARRIER

      return
      end



*     ***********************************
*     *                                 *
*     *          D3dB_r_notZero_Ends    *
*     *                                 *
*     ***********************************

      subroutine D3dB_r_notZero_Ends(nb,A)
      integer nb
      real*8 A(*) 

#include "D3dB.fh"

      integer j,k,q,index,taskid,p

      !**** slab mapping ****
      if (mapping.eq.1) then
      do q=1,nq(nb)
         do j=1,ny(nb)
            index = (nx(nb)+1) + (j-1)*(nx(nb)+2)
     >                         + (q-1)*(nx(nb)+2)*(ny(nb))
            A(index)   = 1.0d0
            A(index+1) = 1.0d0
         end do
      end do


      !**** hilbert mapping ****
      else
        call Parallel2d_taskid_i(taskid)
        do k=1,nz(nb)
        do j=1,ny(nb)

         call D3dB_ijktoindex2p(nb,(nx(nb)+1),j,k,index,p)
         if (p.eq.taskid) A(index) = 1.0d0

         call D3dB_ijktoindex2p(nb,(nx(nb)+2),j,k,index,p)
         if (p.eq.taskid) A(index) = 1.0d0

        end do
        end do
      end if

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

#include "D3dB.fh"

      integer i,np
      real*8 sum

      call Parallel2d_np_i(np)

*     **** sum up dot product on this node ****
      sum = 0.0d0
      do i=1,n2ft3d_map(nb)
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

#include "D3dB.fh"

      integer i,j,k,q,np,nxh,index,taskid,p
      real*8 sum

      nxh = nx(nb)/2
      call Parallel2d_np_i(np)

*     **** sum up dot product on this node ****
      sum = 0.0d0

      !**********************
      !**** slab mapping ****
      !**********************
      if (mapping.eq.1) then
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


      !*************************
      !**** hilbert mapping ****
      !*************************
      else
      call Parallel2d_taskid_i(taskid)
*     ***** kx!=0 plane, so double count *****
      do index=1,nfft3d_map(nb)
            sum = sum + A(index)
      end do
      sum = sum*2.0d0

*     ***** kx==0 plane, so single count *****
      do k=1,nz(nb)
      do j=1,ny(nb)
         i=1
         call D3dB_ijktoindexp(nb,i,j,k,index,p)
         if (p.eq.taskid) then
         sum = sum - A(index)
         end if
      end do
      end do

      end if




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

      subroutine D3dB_cc_Vector_dot(nb,nnfft3d,nn,ne,A,B,sumall)     
      implicit none 
      integer    nb
      integer    nnfft3d,nn,ne
      complex*16 A(*)
      complex*16 B(*)
      real*8     sumall(nn,nn)


#include "D3dB.fh"

      integer i,j,k,q,index,np,taskid,p
      integer index1,index2
      integer n,m,shift1,shift2
      real*8  sum


      call nwpw_timing_start(2)

      call Parallel2d_np_i(np)

      !**********************
      !**** slab mapping ****
      !**********************
      if (mapping.eq.1) then
*     **** sum up dot product on this node ****
      do n=1,ne
      do m=n,ne

        shift1 = (n-1)*nnfft3d
        shift2 = (m-1)*nnfft3d
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


      !*************************
      !**** hilbert mapping ****
      !*************************
      else
      call Parallel2d_taskid_i(taskid)
*     **** sum up dot product on this node ****
      do n=1,ne
      do m=n,ne

        shift1 = (n-1)*nnfft3d
        shift2 = (m-1)*nnfft3d
        sum    = 0.0d0

*       ***** kx!=0 plane, so double count *****
        do index=1,nfft3d_map(nb)
            index1 = index+shift1
            index2 = index+shift2
            sum = sum + dble(A(index1))  * dble(B(index2))
     >                + dimag(A(index1)) * dimag(B(index2))
        end do
        sum = sum*2.0d0

*       ***** kx==0 plane, so single count *****
        do k=1,nz(nb)
        do j=1,ny(nb)
         i=1
         call D3dB_ijktoindexp(nb,i,j,k,index,p)
         if (p.eq.taskid) then
         index1 = index+shift1
         index2 = index+shift2
         sum = sum - dble(A(index1))  * dble(B(index2))
     >             - dimag(A(index1)) * dimag(B(index2))
         end if
        end do
        end do


        sumall(n,m) = sum
        sumall(m,n) = sum
      end do
      end do

      end if


*     **** add up sums from other nodes ****
      if (np.gt.1) then
         call D3dB_Vector_SumAll(nn*ne,sumall)
      end if

      call nwpw_timing_end(2)

      return
      end



*     ***********************************
*     *					*
*     *	     D3dB_cc_Vector_ndot 	*	
*     *					*
*     ***********************************

      subroutine D3dB_cc_Vector_ndot(nb,nnfft3d,ne,A,B,sumall)     
      implicit none 
      integer    nb
      integer    nnfft3d,ne
      complex*16 A(*)
      complex*16 B(*)
      real*8     sumall(ne)


#include "D3dB.fh"

      integer i,j,k,q,index,np,taskid,p
      integer index1,index2
      integer n,shift1,shift2
      real*8  sum


      call nwpw_timing_start(2)

      call Parallel2d_np_i(np)

      !**********************
      !**** slab mapping ****
      !**********************
      if (mapping.eq.1) then
*     **** sum up dot product on this node ****
      do n=1,ne

        shift1 = (n-1)*nnfft3d
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


      !*************************
      !**** hilbert mapping ****
      !*************************
      else
      call Parallel2d_taskid_i(taskid)
*     **** sum up dot product on this node ****
      do n=1,ne

        shift1 = (n-1)*nnfft3d
        shift2 = 0
        sum    = 0.0d0

*       ***** kx!=0 plane, so double count *****
        do index=1,nfft3d_map(nb)
            index1 = index+shift1
            index2 = index+shift2
            sum = sum + dble(A(index1))  * dble(B(index2))
     >                + dimag(A(index1)) * dimag(B(index2))
        end do
        sum = sum*2.0d0

*       ***** kx==0 plane, so single count *****
        do k=1,nz(nb)
        do j=1,ny(nb)
         i=1
         call D3dB_ijktoindexp(nb,i,j,k,index,p)
         if (p.eq.taskid) then
         index1 = index+shift1
         index2 = index+shift2
         sum = sum - dble(A(index1))  * dble(B(index2))
     >             - dimag(A(index1)) * dimag(B(index2))
         end if
        end do
        end do

        sumall(n) = sum
      end do

      end if


*     **** add up sums from other nodes ****
      if (np.gt.1) then
         call D3dB_Vector_SumAll(ne,sumall)
      end if

      call nwpw_timing_end(2)

      return
      end





c*     ***********************************
c*     *					*
c*     *	        D3dB_Vector_SumAll	*	
c*     *					*
c*     ***********************************
c
c      subroutine D3dB_Vector_SumAll(n,sum)
cc     implicit none
c      integer n
c      real*8  sum(*)
c
c#include "bafdecls.fh"
c
c#include "tcgmsg.fh"
c#include "msgtypesf.h"
c#include "errquit.fh"
c
c
c
c      logical value
c      integer MASTER
c      parameter (MASTER=0)
c      integer msglen
c
c*     **** temporary workspace ****
c      integer sumall(2),np
c
c*     **** external functions ****
c      integer  Parallel2d_comm_i
c      external Parallel2d_comm_i
c
c
c
c      call Parallel_np(np)
c      call nwpw_timing_start(2)
c      if (np.gt.1) then
c
c*     ***** allocate temporary space ****
c      value = BA_push_get(mt_dbl,n,'sumall',sumall(2),sumall(1))
c      if (.not. value) call errquit('out of stack memory',0, MA_ERR)
c
c      msglen = n
c
c
c      call dcopy(n,sum,1,dbl_mb(sumall(1)),1)
c
c      call GA_PGROUP_DGOP(Parallel2d_comm_i(),
c     >                    9+MSGDBL,dbl_mb(sumall(1)),n,'+')
cc     call GA_DGOP(9+MSGDBL,dbl_mb(sumall(1)),n,'+')
cc     call DGOP(9+MSGDBL,dbl_mb(sumall(1)),n,'+')
c
c
c      call dcopy(n,dbl_mb(sumall(1)),1,sum,1)
c      value = BA_pop_stack(sumall(2))
c
c      end if
c      call nwpw_timing_end(2)
c      return
c      end
c
c
c*     ***********************************
c*     *					*
c*     *	        D3dB_Vector_ISumAll	*	
c*     *					*
c*     ***********************************
c
c      subroutine D3dB_Vector_ISumAll(n,sum)
cc     implicit none
c      integer n
c      integer  sum(*)
c
c#include "bafdecls.fh"
c#include "errquit.fh"
c
c
c
c#include "tcgmsg.fh"
c#include "msgtypesf.h"
c
c
c
c      integer MASTER
c      parameter (MASTER=0)
c      integer msglen
c      logical value
c
c*     **** temporary workspace ****
c      integer sumall(2)
c
c*     **** external functions ****
c      integer  Parallel2d_comm_i
c      external Parallel2d_comm_i
c
c
c      call nwpw_timing_start(2)
c
c*     ***** allocate temporary space ****
c      value = BA_push_get(mt_int,n,'sumall',sumall(2),sumall(1))
c      if (.not. value) call errquit('out of stack memory',0, MA_ERR)
c
c      msglen = n
c
c
c      call icopy(n,sum,1,int_mb(sumall(1)),1)
c      call GA_PGROUP_IGOP(Parallel2d_comm_i(),
c     >                    9+MSGINT,int_mb(sumall(1)),n,'+')
cc     call GA_IGOP(9+MSGINT,int_mb(sumall(1)),n,'+')
c
c
c      call icopy(n,int_mb(sumall(1)),1,sum,1)
c      value = BA_pop_stack(sumall(2))
c
c      call nwpw_timing_end(2)
c      return
c      end


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
cpgi$r opt=1
      subroutine D3dB_ic_Mul(nb,A,B,C)
      implicit none
      integer    nb
      real*8     A(*)
      complex*16 B(*)
      complex*16 C(*)

#include "D3dB.fh"

      integer i

!$OMP DO
      do i=1,nfft3d_map(nb)
            C(i) = dcmplx(0.0d0,A(i)) * B(i)
      end do
!$OMP END DO

      return
      end

*     ***********************************
*     *                                 *
*     *          D3dB_ic_Mul2           *
*     *                                 *
*     ***********************************
cpgi$r opt=1
      subroutine D3dB_ic_Mul2(nb,A,B)
      implicit none
      integer    nb
      real*8     A(*)
      complex*16 B(*)

#include "D3dB.fh"

      integer i

!$OMP DO
      do i=1,nfft3d_map(nb)
            B(i) = dcmplx(0.0d0,A(i)) * B(i)
      end do
!$OMP END DO

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

#include "D3dB.fh"

      integer i,j,q,index1,index2

      if (mapping.eq.1) then
c      call dcopy(nq(nb2)*ny(nb2)*(nx(nb2)+2),0.0d0,0,B,1)
      call Parallel_shared_vector_zero(.true.,
     >                                 nq(nb2)*ny(nb2)*(nx(nb2)+2),B)
!$OMP DO
      do j=1,ny(nb1)
      do q=1,nq(nb1)
      do i=1,nx(nb1)
         index1 = i + (j-1)*(nx(nb1)+2) + (q-1)*(nx(nb1)+2)*ny(nb1)
         index2 = i + (j-1)*(nx(nb2)+2) + (q-1)*(nx(nb2)+2)*ny(nb2)
         B(index2) = A(index1)
      end do
      end do
      end do
!$OMP END DO

      else
c      call dcopy(n2ft3d(nb2),0.0d0,0,B,1)
      call Parallel_shared_vector_zero(.true.,n2ft3d(nb2),B)
!$OMP DO
      do q=1,nq1(nb1)
      do i=1,nx(nb1)
         index1 = i + (q-1)*(nx(nb1)+2)
         index2 = i + (q-1)*(nx(nb2)+2)
         B(index2) = A(index1)
      end do
      end do
!$OMP END DO
      end if
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

#include "D3dB.fh"

      integer i,j,q,index1,index2

      if (mapping.eq.1) then
c      call dcopy(nq(nb1)*ny(nb1)*(nx(nb1)+2),0.0d0,0,A,1)
      call Parallel_shared_vector_zero(.true.,
     >                                 nq(nb1)*ny(nb1)*(nx(nb1)+2),A)
!$OMP DO
      do j=1,ny(nb1)
      do q=1,nq(nb1)
      do i=1,nx(nb1)
         index1 = i + (j-1)*(nx(nb1)+2) + (q-1)*(nx(nb1)+2)*ny(nb1)
         index2 = i + (j-1)*(nx(nb2)+2) + (q-1)*(nx(nb2)+2)*ny(nb2)
         A(index1) = B(index2)
      end do
      end do
      end do
!$OMP END DO

      else
c      call dcopy(n2ft3d(nb1),0.0d0,0,A,1)
      call Parallel_shared_vector_zero(.true.,n2ft3d(nb1),A)
!$OMP DO
      do q=1,nq1(nb1)
      do i=1,nx(nb1)
         index1 = i + (q-1)*(nx(nb1)+2)
         index2 = i + (q-1)*(nx(nb2)+2)
         A(index1) = B(index2)
      end do
      end do
!$OMP END DO
      end if

      return
      end

*     ***********************************
*     *					*
*     *	   D3dB_timereverse_size	*
*     *					*
*     ***********************************

      integer function D3dB_timereverse_size(nb)
      implicit none
      integer nb

#include "D3dB.fh"

*     **** local variables ****
      integer proc_to,proc_from
      integer pto,np,taskid
      integer phere,itmp,itmp2
      integer index1,index2
      integer it,size
      integer i2,j2,k2
      integer i3,j3,k3
      integer nyh,nzh

      call Parallel2d_taskid_i(taskid)
      call Parallel2d_np_i(np)

      nyh = ny(nb)/2
      nzh = nz(nb)/2

      index1 = 1
      index2 = 1
      do it=0,np-1
         proc_to   = mod(taskid+it,np)
         proc_from = mod(taskid-it+np,np)

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
            call D3dB_ijktoindexp(nb,1,1,i3,itmp,phere)
            call D3dB_ijktoindexp(nb,1,1,j3,itmp2,pto)
            if ((phere.eq.taskid).and.(pto.eq.proc_to)) then
               index1 = index1 + 1
            end if

*           **** unpacking scheme ****
            if ((pto.eq.taskid).and.(phere.eq.proc_from)) then
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
            call D3dB_ijktoindexp(nb,1,i2,i3,itmp,phere)
            call D3dB_ijktoindexp(nb,1,j2,j3,itmp2,pto)
            if ((phere.eq.taskid).and.(pto.eq.proc_to)) then
               index1 = index1 + 1
            end if

*           **** unpacking scheme ****
            if ((pto.eq.taskid).and.(phere.eq.proc_from)) then
               index2 = index2 + 1
            end if
         end do
         end do

      end do
      size = index1
      if (size.lt.index2) size = index2

      D3dB_timereverse_size = size
      return
      end

*     ***********************************
*     *					*
*     *	   D3dB_c_timereverse_init	*
*     *					*
*     ***********************************

      subroutine D3dB_c_timereverse_init(nb)
      implicit none
      integer nb

#include "bafdecls.fh"
#include "D3dB.fh"
#include "errquit.fh"

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
      integer pto,np,taskid
      integer phere
      integer index1,index2,itmp,itmp2
      integer it
      integer i2,j2,k2
      integer i3,j3,k3
      integer nyh,nzh
      logical value

      !**** external functions ****
      integer  D3dB_timereverse_size
      external D3dB_timereverse_size

      call Parallel2d_taskid_i(taskid)
      call Parallel2d_np_i(np)

      nyh = ny(nb)/2
      nzh = nz(nb)/2

*     **** set zplane_size ****
      zplane_size(nb) = D3dB_timereverse_size(nb)
      
*     **** allocate timereverse_blk common block ****
      value = BA_alloc_get(mt_int,zplane_size(nb),
     >                     'iq_to_i1',iq_to_i1(2,nb),iq_to_i1(1,nb))
      value = value.and.
     >        BA_alloc_get(mt_int,zplane_size(nb),
     >                     'iq_to_i2',iq_to_i2(2,nb),iq_to_i2(1,nb))

      value = value.and.
     >        BA_alloc_get(mt_int,(np+1),
     >                     'i1_start',i1_start(2,nb),i1_start(1,nb))
      value = value.and.
     >        BA_alloc_get(mt_int,(np+1),
     >                     'i2_start',i2_start(2,nb),i2_start(1,nb))
      if (.not.value) call errquit('out of heap memory',0, MA_ERR)


  
!MATHIAS
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
            call D3dB_ijktoindexp(nb,1,1,i3,itmp,phere)
            call D3dB_ijktoindexp(nb,1,1,j3,itmp2,pto)
            !call D3dB_ktoqp(nb,i3,qhere,phere)
            !call D3dB_ktoqp(nb,j3,qto,pto)
            if ((phere.eq.taskid).and.(pto.eq.proc_to)) then
c               itmp = 1 + (i2-1)*(nx(nb)/2+1) 
c     >                  + (qhere-1)*(nx(nb)/2+1)*ny(nb)
c               iq_to_i1(index1) = itmp
               int_mb(iq_to_i1(1,nb)+index1-1)=itmp
               index1 = index1 + 1
            end if

*           **** unpacking scheme ****
            !call D3dB_ktoqp(nb,j3,qhere,phere)
            !call D3dB_ktoqp(nb,i3,qfrom,pfrom)
            if ((pto.eq.taskid).and.(phere.eq.proc_from)) then
c               itmp = 1 + (j2-1)*(nx(nb)/2+1) 
c     >                  + (qhere-1)*(nx(nb)/2+1)*ny(nb)
c               iq_to_i2(index2) = itmp
               int_mb(iq_to_i2(1,nb)+index2-1) = itmp2
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
            call D3dB_ijktoindexp(nb,1,i2,i3,itmp,phere)
            call D3dB_ijktoindexp(nb,1,j2,j3,itmp2,pto)
            !call D3dB_ktoqp(nb,i3,qhere,phere)
            !call D3dB_ktoqp(nb,j3,qto,pto)
            if ((phere.eq.taskid).and.(pto.eq.proc_to)) then
c               itmp = 1 + (i2-1)*(nx(nb)/2+1) 
c     >                  + (qhere-1)*(nx(nb)/2+1)*ny(nb)
c               iq_to_i1(index1) = itmp
               int_mb(iq_to_i1(1,nb)+index1-1) = itmp
               index1 = index1 + 1
            end if

*           **** unpacking scheme ****
            !call D3dB_ktoqp(nb,j3,qhere,phere)
            !call D3dB_ktoqp(nb,i3,qfrom,pfrom)
            if ((pto.eq.taskid).and.(phere.eq.proc_from)) then
c              itmp = 1 + (j2-1)*(nx(nb)/2+1) 
c    >                  + (qhere-1)*(nx(nb)/2+1)*ny(nb)
c               iq_to_i2(index2) = itmp
               int_mb(iq_to_i2(1,nb)+index2-1)=itmp2
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
*     *	     D3dB_timereverse_end 	*	
*     *					*
*     ***********************************
      subroutine D3dB_timereverse_end(nb)
      implicit none
      integer nb

#include "bafdecls.fh"
#include "errquit.fh"

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

      value = BA_free_heap(i1_start(2,nb))
      value = value.and.
     >        BA_free_heap(i2_start(2,nb))
      value = value.and.
     >        BA_free_heap(iq_to_i1(2,nb))
      value = value.and.
     >        BA_free_heap(iq_to_i2(2,nb))
      if (.not.value) call errquit('error freeing heap',0, MA_ERR)
      return
      end



      subroutine D3dB_pfft_index1_copy(n,index,a,b)
      implicit none
      integer n
      integer index(*)
      complex*16  a(*),b(*)
      integer i

#ifndef CRAY
!DIR$ ivdep
#endif
!$OMP DO
      do i=1,n
        b(i) = a(index(i))
      end do
!$OMP END DO

      return
      end

      subroutine D3dB_pfft_index2_copy(n,index,a,b)
      implicit none
      integer n
      integer index(*)
      complex*16  a(*),b(*)
      integer i
#ifndef CRAY
!DIR$ ivdep
#endif
!$OMP DO 
      do i=1,n
        b(index(i)) = a(i)
      end do
!$OMP END DO NOWAIT
      return
      end


      subroutine D3dB_pfft_index2_zero(n,index,a)
      implicit none
      integer n
      integer index(*)
      complex*16  a(*)
      integer i
#ifndef CRAY
!DIR$ ivdep
#endif
!$OMP DO
      do i=1,n
        a(index(i)) = 0.0d0
      end do
!$OMP END DO
      return
      end



      subroutine D3dB_pfft_index2_copy_conjg(n,index,a,b)
      implicit none
      integer n
      integer index(*)
      complex*16  a(*),b(*)
      integer i
#ifndef CRAY
!DIR$ ivdep
#endif
!$OMP DO 
      do i=1,n
        b(index(i)) = dconjg(a(i))
      end do
!$OMP END DO NOWAIT
      return
      end



