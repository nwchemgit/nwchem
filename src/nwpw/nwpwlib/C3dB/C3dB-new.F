#define NBLOCKS 2
*
* $Id$
*

*     ***********************************************************
*     *								*
*     *   		   C3dB library				*
*     *		(NWChem implemenation, version 0.1)	        *
*     *								*
*     *   Author - Eric Bylaska					*
*     *   date   - 11/16/01					*
*     *								*
*     ***********************************************************
*   The C3dB (full complex distributed three-dimensional block) library
*is to be used for handling three kinds of data structures.  The first
* data structure, denoted by "r", is a double precision array of
* length (nx)*ny*nz.  The second data structure, denoted by "c", is
* a double complex array of length of (nx)*ny*nz.
*
*   The two data structures are distributed across threads, p, in
* the k (i.e. nz) dimension using a cyclic decomposition.  So that
* a "r" array A is defined as double precision A(nx,ny,nq) on
* each thread.
*
*   Where
*       np = number of threads
*       nq = ceil(nz/np).
*       0 <= p < np
*       1 <= q <= nq
*       1 <= k <= nz
*
*   The mapping of k -> q is defined as:
*
*       k = ((q-1)*np + p) + 1
*       q = ((k-1) - p)/np + 1
*       p = (k-1) mod np
*
*  Libraries used: mpi, blas, fftpack, and compressed_io
*

*  common blocks used in this library:
*
*       integer nq,nx,ny,nz
*   common  / C3dB / nq,nx,ny,nz
*
*   integer q_map(NFFT3),p_map(NFFT3),k_map(NFFT3)
*   common /C3dB_mapping / q_map,p_map,k_map
*
*     integer iq_to_i1((NFFT1)*NFFT2*NSLABS)
*     integer iq_to_i2((NFFT1)*NFFT2*NSLABS)
*     integer i1_start(NPROCS+1)
*     integer i2_start(NPROCS+1)
*     common / c_trans_blk / iq_to_i1,iq_to_i2,i1_start,i2_start


*     ***********************************
*     *					*
*     *	       Mapping_Init_C3dB	*	
*     *					*
*     ***********************************

      subroutine Mapping_Init_C3dB(nb)
      implicit none
      integer nb
     
#include "bafdecls.fh"
#include "errquit.fh"
#include "C3dB.fh"


      integer k,q,p,j
*     integer kn
      integer taskid,np
      logical value


      call Parallel3d_np_i(np)
      call Parallel3d_taskid_i(taskid)


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
     > call errquit('Mapping_init:out of heap memory',1, MA_ERR)



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

      do k=1,nz(nb)
         if (int_mb(p_map(1,nb)+k-1) .eq. taskid) then
            int_mb(k_map(1,nb)+int_mb(q_map(1,nb)+k-1)-1) = k
         end if
      end do

      nfft3d(nb)     = nx(nb)*ny(nb)*nq(nb)
      n2ft3d(nb)     = nfft3d(nb)
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

      value = value.and.BA_alloc_get(mt_int,nz(nb)*nx(nb),
     >                              'q_map2',
     >                               q_map2(2,nb),
     >                               q_map2(1,nb))
      value = value.and.BA_alloc_get(mt_int,nz(nb)*nx(nb),
     >                              'p_map2',
     >                               p_map2(2,nb),
     >                               p_map2(1,nb))

      value = value.and.BA_alloc_get(mt_int,ny(nb)*nx(nb),
     >                              'q_map3',
     >                               q_map3(2,nb),
     >                               q_map3(1,nb))
      value = value.and.BA_alloc_get(mt_int,ny(nb)*nx(nb),
     >                              'p_map3',
     >                               p_map3(2,nb),
     >                               p_map3(1,nb))
      if (.not. value)
     > call errquit('Mapping_init:out of heap memory',1, MA_ERR)


      !**** double grid map1 defined wrt to single grid         ****
      !**** makes expand and contract routines trivial parallel ****
      if (mapping2d.eq.1) then
         if (nb.eq.1) then
           call hilbert2d_map(ny(nb),nz(nb),int_mb(p_map1(1,nb)))
         end if
         call hilbert2d_map(nz(nb),nx(nb),int_mb(p_map2(1,nb)))
         call hilbert2d_map(nx(nb),ny(nb),int_mb(p_map3(1,nb)))
      else
         if (nb.eq.1) then
           call hcurve_map(ny(nb),nz(nb),int_mb(p_map1(1,nb)))
         end if
         call hcurve_map(nz(nb),nx(nb),int_mb(p_map2(1,nb)))
         call hcurve_map(nx(nb),ny(nb),int_mb(p_map3(1,nb)))
      end if



      !**** double grid map1 defined wrt to single grid         ****
      !**** makes expand and contract routines trivial parallel ****
      if (nb.eq.1) then
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
     >                          nz(nb),nx(nb),
     >                          int_mb(p_map2(1,nb)),
     >                          int_mb(q_map2(1,nb)),nq2(nb))
      call generate_map_indexes(taskid,np,
     >                          nx(nb),ny(nb),
     >                          int_mb(p_map3(1,nb)),
     >                          int_mb(q_map3(1,nb)),nq3(nb))

c      if (taskid.eq.0) then
c      write(*,*) taskid,"nq2=",nq2(nb), ny(nb)*nq2(nb)
c      write(*,*) taskid,"nq1=",nq1(nb), nx(nb)*nq1(nb)
c      write(*,*) taskid,"nq3=",nq3(nb), nz(nb)*nq3(nb)
c      write(*,*) 'hilbert map1 nb=',nb
c      do j=0,nz(nb)-1
c        write(*,'(A,80I4)') 'hilbert map:',
c     >   (int_mb(p_map1(1,nb)+k+j*ny(nb)), k=0,ny(nb)-1)
c      end do
c      write(*,*)
c      write(*,*) 'hilbert map2 nb=',nb
c      do j=0,nx(nb)-1
c        write(*,'(A,80I4)') 'hilbert map:',
c     >   (int_mb(p_map2(1,nb)+k+j*nz(nb)), k=0,nz(nb)-1)
c      end do
c      write(*,*)
c      write(*,*) 'hilbert map3 nb=',nb
c      do j=0,ny(nb)-1
c        write(*,'(A,80I4)') 'hilbert map:',
c     >   (int_mb(p_map3(1,nb)+k+j*nx(nb)), k=0,nx(nb)-1)
c      end do
c      write(*,*)
c      end if

      nfft3d(nb) = nx(nb)*nq1(nb)
      if ((ny(nb)*nq2(nb)).gt.nfft3d(nb)) nfft3d(nb) = ny(nb)*nq2(nb)
      if ((nz(nb)*nq3(nb)).gt.nfft3d(nb)) nfft3d(nb) = nz(nb)*nq3(nb)
      n2ft3d(nb) = nfft3d(nb)

      nfft3d_map(nb) = nz(nb)*nq3(nb)
      n2ft3d_map(nb) = nx(nb)*nq1(nb)

      end if

      return
      end

*     ***********************************
*     *					*
*     *	          C3dB_end   		*	
*     *					*
*     ***********************************
      subroutine C3dB_end(nb)
      implicit none
      integer nb

#include "bafdecls.fh"
#include "errquit.fh"
#include "C3dB.fh"


*     *** hilbert tranpose data structure ****
      integer h_iq_to_i1(2,6,NBLOCKS)
      integer h_iq_to_i2(2,6,NBLOCKS)
      integer h_i1_start(2,6,NBLOCKS)
      integer h_i2_start(2,6,NBLOCKS)
      common / c_trans_blk_ijk / h_iq_to_i1,
     >                           h_iq_to_i2,
     >                           h_i1_start,
     >                           h_i2_start

      integer iq_to_i1(2,NBLOCKS)
      integer iq_to_i2(2,NBLOCKS)
      integer i1_start(2,NBLOCKS)
      integer i2_start(2,NBLOCKS)
      common / c_trans_blk / iq_to_i1,iq_to_i2,i1_start,i2_start

#ifndef MPI
      integer Nchannels(NBLOCKS)
      integer channel_proc(2,NBLOCKS)
      integer channel_type(2,NBLOCKS)
      common / c_channel_blk / channel_proc,channel_type,Nchannels
#endif

      logical value
      integer i


      call C3dB_fft_end(nb)
      value = .true.

      !**** slab mappings ****
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
     > call errquit('C3dB_end:freeing heap memory',0, MA_ERR)
      return
      end
         
*     ***********************************
*     *					*
*     *	          C3dB_qtok   		*	
*     *					*
*     ***********************************

      subroutine C3dB_qtok(nb,q,k)      
      implicit none
      integer nb
      integer q,k

#include "bafdecls.fh"
#include "C3dB.fh"
      
      k = int_mb(k_map(1,nb)+q-1)

      return
      end

*     ***********************************
*     *					*
*     *	          C3dB_ktoqp  		*	
*     *					*
*     ***********************************

      subroutine C3dB_ktoqp(nb,k,q,p)      
      implicit none
      integer nb
      integer k,q,p

#include "bafdecls.fh"
#include "C3dB.fh"
      
      q = int_mb(q_map(1,nb)+k-1)
      p = int_mb(p_map(1,nb)+k-1)
      return
      end


*     ***********************************
*     *                                 *
*     *           C3dB_ijktoindexp      *
*     *                                 *
*     ***********************************

      subroutine C3dB_ijktoindexp(nb,i,j,k,indx,p)
      implicit none
      integer nb
      integer i,j,k
      integer indx,p

#include "bafdecls.fh"
#include "C3dB.fh"

      integer q

      !**** slab mapping ****
      if (mapping.eq.1) then
      q = int_mb(q_map(1,nb)+k-1)
      p = int_mb(p_map(1,nb)+k-1)

      indx = i + (j-1)*(nx(nb)) + (q-1)*(nx(nb))*ny(nb)

      !**** hilbert mapping ****
      else
      q = int_mb(q_map3(1,nb)+(i-1)+(j-1)*nx(nb))
      p = int_mb(p_map3(1,nb)+(i-1)+(j-1)*nx(nb))

      indx = k + (q-1)*nz(nb)


      end if

      return
      end



*     ***********************************
*     *                                 *
*     *           C3dB_ijktoindex1p     *
*     *                                 *
*     ***********************************

      subroutine C3dB_ijktoindex1p(nb,i,j,k,indx,p)
      implicit none
      integer nb
      integer i,j,k
      integer indx,p

#include "bafdecls.fh"
#include "C3dB.fh"

      integer q

      !**** slab mapping ***
      if (mapping.eq.1) then
      q = int_mb(q_map(1,nb)+j-1)
      p = int_mb(p_map(1,nb)+j-1)

      indx = i + (k-1)*nx(nb) + (q-1)*nx(nb)*nz(nb)

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
*     *           C3dB_ijktoindex2p     *
*     *                                 *
*     ***********************************

      subroutine C3dB_ijktoindex2p(nb,i,j,k,indx,p)
      implicit none
      integer nb
      integer i,j,k
      integer indx,p

#include "bafdecls.fh"
#include "C3dB.fh"


      integer q

      !**** slab mapping ****
      if (mapping.eq.1) then
      q = int_mb(q_map(1,nb)+j-1)
      p = int_mb(p_map(1,nb)+j-1)

      indx = i + (k-1)*(nx(nb)) + (q-1)*(nx(nb))*ny(nb)

      !**** hilbert mapping ****
      else
      q = int_mb(q_map1(1,nb)+(j-1)+(k-1)*ny(nb))
      p = int_mb(p_map1(1,nb)+(j-1)+(k-1)*ny(nb))

      indx = i + (q-1)*nx(nb)
      end if


      return
      end





*     ***********************************
*     *                                 *
*     *         C3dB_ijk_to_srqp        *
*     *                                 *
*     ***********************************

      subroutine C3dB_ijk_to_srqp(nb,i,j,k, s,r,q,p)
      implicit none
      integer nb
      integer i,j,k
      integer s,r,q,p

#include "bafdecls.fh"
#include "C3dB.fh"


c     q = q_map(k)
c     p = p_map(k)

      s = i
      r = j
      q = int_mb(q_map(1,nb)+k-1)
      p = int_mb(p_map(1,nb)+k-1)
      return
      end




*     ***********************************
*     *					*
*     *	        C3dB_nfft3d		*
*     *					*
*     ***********************************

      subroutine C3dB_nfft3d(nb,nfft3d_out)
      implicit none
      integer nb
      integer nfft3d_out

#include "C3dB.fh"

      nfft3d_out = nfft3d(nb)
      return
      end



*     ***********************************
*     *                                 *
*     *         C3dB_nfft3d_map         *
*     *                                 *
*     ***********************************

      subroutine C3dB_nfft3d_map(nb,nfft3d_out)
      implicit none
      integer nb
      integer nfft3d_out

#include "C3dB.fh"

      nfft3d_out = nfft3d_map(nb)
      return
      end


*     ***********************************
*     *					*
*     *	        C3dB_n2ft3d		*
*     *					*
*     ***********************************

      subroutine C3dB_n2ft3d(nb,n2ft3d_out)
      implicit none
      integer nb
      integer n2ft3d_out

#include "C3dB.fh"

      n2ft3d_out = n2ft3d(nb)
      return
      end


*     ***********************************
*     *                                 *
*     *         C3dB_n2ft3d_map         *
*     *                                 *
*     ***********************************

      subroutine C3dB_n2ft3d_map(nb,n2ft3d_out)
      implicit none
      integer nb
      integer n2ft3d_out

#include "C3dB.fh"

      n2ft3d_out = n2ft3d_map(nb)
      return
      end


*     ***********************************
*     *                                 *
*     *         C3dB_nqq                *
*     *                                 *
*     ***********************************

      subroutine C3dB_nqq(nb,nqtmp)
      implicit none
      integer nb
      integer nqtmp

#include "C3dB.fh"

      !**** slab mapping ****
      if (mapping.eq.1) then
         nqtmp = ny(nb)*nq(nb)
      !**** hilbert mapping ****
      else
         nqtmp = nq(nb)
      end if

      return
      end



*     ***********************************
*     *					*
*     *	        C3dB_nq			*	
*     *					*
*     ***********************************

      subroutine C3dB_nq(nb,nqtmp)
      implicit none
      integer nb
      integer nqtmp

#include "C3dB.fh"

      nqtmp = nq(nb)

      return 
      end

*     ***********************************
*     *					*
*     *	        C3dB_nx			*	
*     *					*
*     ***********************************
     
      subroutine C3dB_nx(nb,nxtmp)
      implicit none
      integer nb
      integer nxtmp
      
#include "C3dB.fh"

      nxtmp = nx(nb)
      return
      end

*     ***********************************
*     *					*
*     *	        C3dB_ny			*	
*     *					*
*     ***********************************

      subroutine C3dB_ny(nb,nytmp)
      implicit none
      integer nb
      integer nytmp
      
#include "C3dB.fh"

      nytmp = ny(nb)
      return
      end

*     ***********************************
*     *					*
*     *	        C3dB_nz			*	
*     *					*
*     ***********************************

      subroutine C3dB_nz(nb,nztmp)
      implicit none
      integer nb
      integer nztmp
      
#include "C3dB.fh"

      nztmp = nz(nb)
      return
      end


*     ***********************************
*     *					*
*     *	        C3dB_Init		*	
*     *					*
*     ***********************************

      subroutine C3dB_Init(nb,nx_in,ny_in,nz_in,map_in)
      implicit none
      integer nb
      integer nx_in,ny_in,nz_in
      integer map_in

#include "C3dB.fh"

      !**** local variables ****
      integer MASTER
      parameter (MASTER=0)
      integer taskid,np

      call Parallel3d_np_i(np)
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
          call errquit('C3dB_Init: mapping error',0,0)
        end if
        if (mod(nx_in,2).ne.0) then
          if (taskid.eq.MASTER) then
           write(6,*)
     >      'Error: ngrid(1) must be even (ngrid(1) mod 2 == 0)'
           write(6,*) 'Error: ngrid(1)=',nx_in
          end if
          call errquit('C3dB_Init: slab mapping error',0,0)
        end if
      end if


      if (map_in.ge.2) then
        if (np.gt.(ny_in*nz_in)) then
          if (taskid.eq.MASTER) then
           write(6,*) 'Error: np > MIN(ngrid(2)*ngrid(3),',
     >                ' ngrid(1)*ngrid(2),',
     >                ' ngrid(1)*ngrid(3))'
           write(6,*) 'Error: np > ngrid(2)*ngrid(3)'
           write(6,*) 'Error: for the Hilbert decomposition the',
     >                 ' the number of processors must ',
     >                 ' be in the range ( 1 ...',
     >                   ny_in*nz_in,')'
          end if
          call errquit('C3dB_Init: Hilbert mapping error',0,0)
        end if
        if (np.gt.(nx_in*ny_in)) then
          if (taskid.eq.MASTER) then
           write(6,*) 'Error: np > MIN(ngrid(2)*ngrid(3),',
     >                ' ngrid(1)*ngrid(2),',
     >                ' ngrid(1)*ngrid(3))'
           write(6,*) 'Error: np > ngrid(1)*ngrid(2)'
           write(6,*) 'Error: for the Hilbert decomposition the',
     >                 ' the number of processors must ',
     >                 ' be in the range ( 1 ...',
     >                   nx_in*ny_in,')'
          end if
          call errquit('C3dB_Init: Hilbert mapping error',0,0)
        end if
        if (np.gt.(nx_in*nz_in)) then
          if (taskid.eq.MASTER) then
           write(6,*) 'Error: np > MIN(ngrid(2)*ngrid(3),',
     >                ' ngrid(1)*ngrid(2),',
     >                ' ngrid(1)*ngrid(3))'
           write(6,*) 'Error: np > ngrid(1)*ngrid(3)'
           write(6,*) 'Error: for the Hilbert decomposition the',
     >                 ' the number of processors must ',
     >                 ' be in the range ( 1 ...',
     >                   nx_in*nz_in,')'
          end if
          call errquit('C3dB_Init: Hilbert mapping error',0,0)
        end if
        if (mod(nx_in,2).ne.0) then
          if (taskid.eq.MASTER) then
           write(6,*)
     >      'Error: ngrid(1) must be even (ngrid(1) mod 2 == 0)'
           write(6,*) 'Error: ngrid(1)=',nx_in
          end if
          call errquit('C3dB_Init: Hilbert mapping error',0,0)
        end if
      end if


*     ***** initialize C3dB common block *****
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
      call Mapping_Init_C3dB(nb)
      if (mapping.eq.1) call C3dB_c_transpose_jk_init(nb)
      if (mapping.eq.2) call C3dB_c_transpose_ijk_init(nb)

#ifndef MPI
      call C3dB_channel_init(nb)
#endif

      call C3dB_fft_init(nb)

      return
      end


*     ***********************************
*     *					*
*     *	        C3dB_(c,r)_Zero  	*	
*     *					*
*     ***********************************

      subroutine C3dB_c_Zero(nb,A)     
      implicit none 
      integer nb
      complex*16 A(*)

#include "C3dB.fh"

      call dcopy(2*nfft3d_map(nb),0.0d0,0,A,1)
      return
      end


      subroutine C3dB_r_Zero(nb,A)     
      implicit none 
      integer nb
      real*8  A(*)

  
#include "C3dB.fh"

      call dcopy(n2ft3d_map(nb),0.0d0,0,A,1)
      return
      end



*     ***********************************
*     *					*
*     *	        C3dB_(c,r)_Copy	*	
*     *					*
*     ***********************************

      subroutine C3dB_c_Copy(nb,A,B)     
      implicit none 
      integer nb
      complex*16 A(*)
      complex*16 B(*)

#include "C3dB.fh"

      call dcopy(2*nfft3d_map(nb),A,1,B,1)
      return
      end

      subroutine C3dB_r_Copy(nb,A,B)     
      implicit none 
      integer nb
      real*8 A(*)
      real*8 B(*)

#include "C3dB.fh"

      call dcopy(n2ft3d_map(nb),A,1,B,1)
      return
      end

      subroutine C3dB_t_Copy(nb,A,B)
      implicit none
      integer nb
      real*8 A(*)
      real*8 B(*)

#include "C3dB.fh"

      call dcopy(nfft3d_map(nb),A,1,B,1)
      return
      end

*     ***********************************
*     *					*
*     *	        C3dB_ct_Copy      	*	
*     *					*
*     ***********************************

      subroutine C3dB_ct_Copy(nb,A,B)     
      implicit none 
      integer nb
      complex*16 A(*)
      real*8     B(*)
      
#include "C3dB.fh"
      
      integer i

      do i=1,nfft3d_map(nb)
         B(i) = dble(A(i))
      end do
      return
      end

*     ***********************************
*     *                                 *
*     *         C3dB_tc_Copy            *       
*     *                                 *
*     ***********************************

      subroutine C3dB_tc_Copy(nb,A,B)
      implicit none
      integer nb
      real*8     A(*)
      complex*16 B(*)

#include "C3dB.fh"

      integer i

!$OMP DO
      do i=1,nfft3d_map(nb)
         B(i) = dcmplx(A(i),0.0d0)
      end do
!$OMP END DO
      return
      end




*     ***********************************
*     *                                 *
*     *         C3dB_fft_init           *
*     *                                 *
*     ***********************************

      subroutine C3dB_fft_init(nb)
      implicit none
      integer nb

#include "bafdecls.fh"
#include "errquit.fh"

#include "C3dB.fh"


      integer tmpx(2,NBLOCKS),tmpy(2,NBLOCKS),tmpz(2,NBLOCKS)
      common    / C3dB_fft / tmpx,tmpy,tmpz

      logical value


      value = BA_alloc_get(mt_dcpl,(nfft3d(nb)),
     >        'fttmpx',tmpx(2,nb),tmpx(1,nb))
      value = value.and.
     >        BA_alloc_get(mt_dcpl,(nfft3d(nb)),
     >        'fttmpy',tmpy(2,nb),tmpy(1,nb))
      value = value.and.
     >        BA_alloc_get(mt_dcpl,(nfft3d(nb)),
     >        'fttmpz',tmpz(2,nb),tmpz(1,nb))
      if (.not. value)
     >   call errquit('C3dB_fft_init:out of heap memory',0, MA_ERR)


#ifdef MLIB
      call z1dfft(dcpl_mb(tmpx(1,nb)),nx(nb),
     >            dcpl_mb(tmpx(1,nb)),-3,ierr)
      call z1dfft(dcpl_mb(tmpx(1,nb)),ny(nb),
     >            dcpl_mb(tmpy(1,nb)),-3,ierr)
      call z1dfft(dcpl_mb(tmpx(1,nb)),nz(nb),
     >            dcpl_mb(tmpz(1,nb)),-3,ierr)

#else
      call dcffti(nx(nb),dcpl_mb(tmpx(1,nb)))
      call dcffti(ny(nb),dcpl_mb(tmpy(1,nb)))
      call dcffti(nz(nb),dcpl_mb(tmpz(1,nb)))
#endif

      return
      end


*     ***********************************
*     *                                 *
*     *         C3dB_fft_end            *
*     *                                 *
*     ***********************************

      subroutine C3dB_fft_end(nb)
      implicit none
      integer nb

#include "bafdecls.fh"
#include "errquit.fh"

#include "C3dB.fh"


      integer tmpx(2,NBLOCKS),tmpy(2,NBLOCKS),tmpz(2,NBLOCKS)
      common    / C3dB_fft / tmpx,tmpy,tmpz

      logical value

      value =           BA_free_heap(tmpx(2,nb))
      value = value.and.BA_free_heap(tmpy(2,nb))
      value = value.and.BA_free_heap(tmpz(2,nb))
      if (.not.value)
     >   call errquit(
     >   'C3dB_fft_end:error deallocatingof heap memory',0, MA_ERR)

      return
      end




     
*     ***********************************
*     *					*
*     *	        C3dB_cr_fft3b		*
*     *					*
*     ***********************************

      subroutine C3dB_cr_fft3b(nb,A)

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
*       uses - C3dB_c_transpose_jk, dcopy           *
*                                                   *
*****************************************************

      implicit none
      integer nb
      complex*16  A(*)

#include "bafdecls.fh"
#include "C3dB.fh"
#include "errquit.fh"


      integer tmpx(2,NBLOCKS),tmpy(2,NBLOCKS),tmpz(2,NBLOCKS)
      common    / C3dB_fft / tmpx,tmpy,tmpz


*     *** local variables ***
      integer i,j,q,indx

c     complex*16  tmp1(*)
c     complex*16  tmp2(*)
c     real*8      tmp3(*)
      !integer nfft3d
      integer tmp1(2),tmp2(2),ierr
      logical value



      call nwpw_timing_start(1)

*     ***** allocate temporary space ****
      !call C3dB_nfft3d(nb,nfft3d)
      value = BA_push_get(mt_dcpl,(nfft3d(nb)), 
     >                    'ffttmp1',tmp1(2),tmp1(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,(nfft3d(nb)), 
     >                   'ffttmp2',tmp2(2),tmp2(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)




      !**********************
      !**** slab mapping ****
      !**********************
      if (mapping.eq.1) then

*     ********************************************
*     ***         Do a transpose of A          ***
*     ***      A(kx,kz,ky) <- A(kx,ky,kz)      ***
*     ********************************************
c     call C3dB_c_transpose_jk(nb,A,dcpl_mb(tmp2(1)),dbl_mb(tmp3(1)))

*     *************************************************
*     ***     do fft along kz dimension             ***
*     ***   A(kx,nz(nb),ky) <- fft1d^(-1)[A(kx,kz,ky)]  ***
*     *************************************************
#ifdef MLIB
      !call dcffti(nz(nb),dcpl_mb(tmp1(1)))
      do q=1,nq(nb)
      do i=1,nx(nb)
         indx = i + (q-1)*nx(nb)*nz(nb)
         call zcopy(nz(nb),A(indx),nx(nb),dcpl_mb(tmp2(1)),1)
         call z1dfft(dcpl_mb(tmp2(1)),nz(nb),
     >               dcpl_mb(tmpz(1,nb)),-2,ierr)
         call zcopy(nz(nb),dcpl_mb(tmp2(1)),1,A(indx),nx(nb))
      end do
      end do
#else
      !call dcffti(nz(nb),dcpl_mb(tmp1(1)))
      do q=1,nq(nb)
      do i=1,nx(nb)
         indx = i + (q-1)*nx(nb)*nz(nb)
         call zcopy(nz(nb),A(indx),nx(nb),dcpl_mb(tmp2(1)),1)
         call dcfftb(nz(nb),dcpl_mb(tmp2(1)),dcpl_mb(tmpz(1,nb)))
         call zcopy(nz(nb),dcpl_mb(tmp2(1)),1,A(indx),nx(nb))
      end do
      end do
#endif

*     ********************************************
*     ***         Do a transpose of A          ***
*     ***      A(kx,ky,nz(nb)) <- A(kx,nz(nb),ky)      ***
*     ********************************************
      call C3dB_c_transpose_jk(nb,A,dcpl_mb(tmp2(1)),dcpl_mb(tmp1(1)))

*     *************************************************
*     ***     do fft along ky dimension             ***
*     ***   A(kx,ny(nb),nz(nb)) <- fft1d^(-1)[A(kx,ky,nz(nb))]  ***
*     *************************************************
#ifdef MLIB
      do q=1,nq(nb)
      do i=1,nx(nb)
         indx = i + (q-1)*nx(nb)*ny(nb)
         call zcopy(ny(nb),A(indx),nx(nb),dcpl_mb(tmp2(1)),1)
         call z1dfft(dcpl_mb(tmp2(1)),ny(nb),
     >               dcpl_mb(tmpy(1,nb)),-2,ierr)
         call zcopy(ny(nb),dcpl_mb(tmp2(1)),1,A(indx),nx(nb))
      end do
      end do
#else
      !call dcffti(ny(nb),dcpl_mb(tmp1(1)))
      do q=1,nq(nb)
      do i=1,nx(nb)
         indx = i + (q-1)*nx(nb)*ny(nb)
         call zcopy(ny(nb),A(indx),nx(nb),dcpl_mb(tmp2(1)),1)
         call dcfftb(ny(nb),dcpl_mb(tmp2(1)),dcpl_mb(tmpy(1,nb)))
         call zcopy(ny(nb),dcpl_mb(tmp2(1)),1,A(indx),nx(nb))
      end do
      end do
#endif

*     *************************************************
*     ***     do fft along kx dimension             ***
*     ***   A(nx(nb),ny(nb),nz(nb)) <- fft1d^(-1)[A(kx,ny(nb),nz(nb))]  ***
*     *************************************************
#ifdef MLIB
      indx = 1
      do q=1,nq(nb)
      do j=1,ny(nb)
         call z1dfft(A(indx),nx(nb),
     >               dcpl_mb(tmpx(1,nb)),-2,ierr)
         indx = indx + nx(nb)
      end do
      end do
#else
      !call dcffti(nx(nb),dcpl_mb(tmp1(1)))
      do q=1,nq(nb)
      do j=1,ny(nb)
         indx = 1 + (j-1)*nx(nb) + (q-1)*nx(nb)*ny(nb)
         call zcopy(nx(nb),A(indx),1,dcpl_mb(tmp2(1)),1)
         call dcfftb(nx(nb),dcpl_mb(tmp2(1)),dcpl_mb(tmpx(1,nb)))
         call zcopy(nx(nb),dcpl_mb(tmp2(1)),1,A(indx),1)
      end do
      end do
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
      indx = 1
      do q=1,nq3(nb)
         !indx = 1 + (q-1)*nz(nb)
         call dcfftb(nz(nb),A(indx),dcpl_mb(tmpz(1,nb)))
         indx = indx + nz(nb)
      end do
#endif

      call C3dB_c_transpose_ijk(nb,3,A,dcpl_mb(tmp1(1)),
     >                                 dcpl_mb(tmp2(1)))


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
      indx = 1
      do q=1,nq2(nb)
         !indx = 1 + (q-1)*ny(nb)
         call dcfftb(ny(nb),A(indx),dcpl_mb(tmpy(1,nb)))
         indx = indx + ny(nb)
      end do
#endif

      call C3dB_c_transpose_ijk(nb,4,A,dcpl_mb(tmp1(1)),
     >                                 dcpl_mb(tmp2(1)))

*     *************************************************
*     ***     do fft along kx dimension             ***
*     ***   A(nx(nb),ny(nb),nz(nb)) <- fft1d^(-1)[A(kx,ny(nb),nz(nb))]  ***
*     *************************************************
#ifdef MLIB
      indx = 1
      do q=1,nq1(nb)
         !indx = 1 + (q-1)*nx(nb)
         call z1dfft(A(indx),nx(nb),dcpl_mb(tmpx(1,nb)),-2,ierr)
         indx = indx + nx(nb)
      end do
#else
      indx = 1
      do q=1,nq1(nb)
         !indx = 1 + (q-1)*ny(nb)
         call dcfftb(nx(nb),A(indx),dcpl_mb(tmpx(1,nb)))
         indx = indx + nx(nb)
      end do
#endif



      end if
    
*     **** deallocate temporary space  ****
      value = BA_pop_stack(tmp2(2))
      value = BA_pop_stack(tmp1(2))

      call nwpw_timing_end(1)
      return
      end





*     ***********************************
*     *					*
*     *	        C3dB_rc_fft3f		*
*     *					*
*     ***********************************

      subroutine C3dB_rc_fft3f(nb,A)

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
#include "C3dB.fh"
#include "errquit.fh"


      integer tmpx(2,NBLOCKS),tmpy(2,NBLOCKS),tmpz(2,NBLOCKS)
      common    / C3dB_fft / tmpx,tmpy,tmpz


*     *** local variables ***
      integer i,j,q,indx

      !integer nfft3d
      integer tmp1(2),tmp2(2)
      logical value


      call nwpw_timing_start(1)

*     ***** allocate temporary space ****
      !call C3dB_nfft3d(nb,nfft3d)
      value = BA_push_get(mt_dcpl,(nfft3d(nb)),'tmp1',tmp1(2),tmp1(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,(nfft3d(nb)),'tmp2',tmp2(2),tmp2(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)



      !**********************
      !**** slab mapping ****
      !**********************
      if (mapping.eq.1) then

*     ********************************************
*     ***     do fft along nx(nb) dimension        ***
*     ***   A(kx,ny(nb),nz(nb)) <- fft1d[A(nx(nb),ny(nb),nz(nb))]  ***
*     ********************************************
      !call dcffti(nx(nb),dcpl_mb(tmp1(1)))
      do q=1,nq(nb)
      do j=1,ny(nb)
         indx = 1 + (j-1)*nx(nb) + (q-1)*nx(nb)*ny(nb)
         call zcopy((nx(nb)),A(indx),1,dcpl_mb(tmp2(1)),1)
         call dcfftf(nx(nb),dcpl_mb(tmp2(1)),dcpl_mb(tmpx(1,nb)))
         call zcopy(nx(nb),dcpl_mb(tmp2(1)),1,A(indx),1)
      end do
      end do

*     ********************************************
*     ***     do fft along ny(nb) dimension        ***
*     ***   A(kx,ky,nz(nb)) <- fft1d[A(kx,ny(nb),nz(nb))]  ***
*     ********************************************
      !call dcffti(ny(nb),dcpl_mb(tmp1(1)))
      do q=1,nq(nb)
      do i=1,nx(nb)
         indx = i + (q-1)*nx(nb)*ny(nb)
         call zcopy(ny(nb),A(indx),nx(nb),dcpl_mb(tmp2(1)),1)
         call dcfftf(ny(nb),dcpl_mb(tmp2(1)),dcpl_mb(tmpy(1,nb)))
         call zcopy(ny(nb),dcpl_mb(tmp2(1)),1,A(indx),nx(nb))
      end do
      end do


*     ********************************************
*     ***         Do a transpose of A          ***
*     ***      A(ky,nz(nb),ky) <- A(kx,ky,nz(nb))      ***
*     ********************************************
      call C3dB_c_transpose_jk(nb,A,dcpl_mb(tmp2(1)),dcpl_mb(tmp1(1)))


*     ********************************************
*     ***     do fft along nz(nb) dimension        ***
*     ***   A(kx,kz,ky) <- fft1d[A(kx,nz(nb),ky)]  ***
*     ********************************************
      !call dcffti(nz(nb),dcpl_mb(tmp1(1)))
      do q=1,nq(nb)
      do i=1,nx(nb)
         indx = i + (q-1)*nx(nb)*ny(nb)
         call zcopy(nz(nb),A(indx),nx(nb),dcpl_mb(tmp2(1)),1)
         call dcfftf(nz(nb),dcpl_mb(tmp2(1)),dcpl_mb(tmpz(1,nb)))
         call zcopy(nz(nb),dcpl_mb(tmp2(1)),1,A(indx),nx(nb))
      end do
      end do

*     ********************************************
*     ***         Do a transpose of A          ***
*     ***      A(kx,ky,kz) <- A(kx,kz,ky)      ***
*     ********************************************
c     call C3dB_c_transpose_jk(nb,A,dcpl_mb(tmp2(1)),dcpl_mb(tmp1(1)))



      !*************************
      !**** hilbert mapping ****
      !*************************
      else

*     ********************************************
*     ***     do fft along nx(nb) dimension        ***
*     ***   A(kx,ny(nb),nz(nb)) <- fft1d[A(nx(nb),ny(nb),nz(nb))]  ***
*     ********************************************
#ifdef MLIB
      indx = 1
      do q=1,nq1(nb)
         !indx = 1 + (q-1)*nx(nb)
         call z1dfft(A(indx),nx(nb),dcpl_mb(tmpx(1,nb)),1,ierr)
         indx = indx + nx(nb)
      end do
#else
      indx = 1
      do q=1,nq1(nb)
         !indx = 1 + (q-1)*nx(nb)
         call dcfftf(nx(nb),A(indx),dcpl_mb(tmpx(1,nb)))
         indx = indx + nx(nb)
      end do
#endif

      call C3dB_c_transpose_ijk(nb,1,A,dcpl_mb(tmp1(1)),
     >                                 dcpl_mb(tmp2(1)))

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
      indx = 1
      do q=1,nq2(nb)
         !indx = 1 + (q-1)*ny(nb)
         call dcfftf(ny(nb),A(indx),dcpl_mb(tmpy(1,nb)))
         indx = indx + ny(nb)
      end do
#endif

      call C3dB_c_transpose_ijk(nb,2,A,dcpl_mb(tmp1(1)),
     >                                 dcpl_mb(tmp2(1)))

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
      indx = 1
      do q=1,nq3(nb)
         !indx = 1 + (q-1)*nz(nb)
         call dcfftf(nz(nb),A(indx),dcpl_mb(tmpz(1,nb)))
         indx = indx + nz(nb)
      end do
#endif



      end if

*     **** deallocate temporary space  ****
      value = BA_pop_stack(tmp2(2))
      value = BA_pop_stack(tmp1(2))

      call nwpw_timing_end(1)
      return
      end





*     ***********************************
*     *					*
*     *	       C3dB_(c,r)_SMul 	        *	
*     *					*
*     ***********************************

*  This routine performs the operation	C = scale * A
* where scale is a real*8 number.

      subroutine C3dB_c_SMul(nb,scale,A,C)     
      implicit none 
      integer    nb
      real*8     scale
      complex*16 A(*)
      complex*16 C(*)

#include "C3dB.fh"

      integer i

!$OMP DO 
      do i=1,nfft3d_map(nb)
         C(i) = scale*A(i)
      end do
!$OMP END DO 
      return
      end


      subroutine C3dB_c_SMul1(nb,scale,A)
      implicit none
      integer    nb
      real*8     scale
      complex*16 A(*)

#include "C3dB.fh"

      integer i

!$OMP DO
      do i=1,nfft3d_map(nb)
         A(i) = scale*A(i)
      end do
!$OMP END DO
      return
      end


      subroutine C3dB_b_SMul1(nb,scale,A)
      implicit none
      integer    nb        
      real*8     scale
      complex*16 A(*)
#include "C3dB.fh"
      integer i     

!$OMP DO
      do i=1,n2ft3d_map(nb)
         A(i) = scale*A(i)
      end do      
!$OMP END DO
      return
      end



      subroutine C3dB_r_SMul(nb,scale,A,C)     
      implicit none 
      integer nb
      real*8     scale
      real*8 A(*)
      real*8 C(*)

#include "C3dB.fh"

      integer i

!$OMP DO 
      do i=1,n2ft3d_map(nb)
         C(i) = scale*A(i)
      end do
!$OMP END DO 
      return
      end 

      subroutine C3dB_r_SMul1(nb,scale,A)
      implicit none
      integer nb
      real*8     scale
      real*8 A(*)

#include "C3dB.fh"

      integer i

!$OMP DO 
      do i=1,n2ft3d_map(nb)
         A(i) = scale*A(i)
      end do
!$OMP END DO 
      return
      end


      subroutine C3dB_c_ZMul(nb,scale,A,C)     
      implicit none 
      integer    nb
      complex*16 scale
      complex*16 A(*)
      complex*16 C(*)

#include "C3dB.fh"

      integer i

!$OMP DO 
      do i=1,nfft3d_map(nb)
         C(i) = scale*A(i)
      end do
!$OMP END DO 
      return
      end


*     ***********************************
*     *                                 *
*     *        C3dB_rc_SMul             *
*     *                                 *
*     ***********************************

*  This routine performs the operation  C = scale * A
* where scale and A are real*8 numbers.

      subroutine C3dB_rc_SMul(nb,scale,A,C)
      implicit none
      integer    nb
      real*8     scale
      real*8 A(*)
      complex*16 C(*)

#include "C3dB.fh"

      integer i

!$OMP DO 
      do i=1,n2ft3d_map(nb)
         C(i) = dcmplx(scale*A(i),0.0d0)
      end do
!$OMP END DO 
      return
      end

*     ***********************************
*     *					*
*     *	       C3dB_cr_aSqrpy	 	*	
*     *					*
*     ***********************************

*  This routine performs the operation	C = C + w*A * A

      subroutine C3dB_cr_aSqrpy(nb,w,A,C)     
      implicit none 
      integer    nb
      real*8     w
      complex*16 A(*)
      real*8     C(*),ar,ai

#include "C3dB.fh"

      integer i

!$OMP DO 
      do i=1,n2ft3d_map(nb)
         ar=dble(A(i))
         ai=dimag(A(i))
         C(i) = C(i) + w*(ar*ar+ai*ai)
      end do
!$OMP END DO 
      return
      end




*     ***********************************
*     *					*
*     *	       C3dB_cr_Sqr	 	*	
*     *					*
*     ***********************************

*  This routine performs the operation	C = A * A

      subroutine C3dB_cr_Sqr(nb,A,C)     
      implicit none 
      integer    nb
      complex*16 A(*)
      real*8     C(*),ar,ai

#include "C3dB.fh"

      integer i

!$OMP DO 
      do i=1,n2ft3d_map(nb)
         ar=dble(A(i))
         ai=dimag(A(i))
         C(i) = ar*ar + ai*ai
      end do
!$OMP END DO 
      return
      end

      
*     ***********************************
*     *                                 *
*     *        C3dB_cc_absSqrt1         *       
*     *                                 *
*     ***********************************

*  This routine performs the operation  A = sqrt(abs(A))
      
      subroutine C3dB_cc_absSqrt1(nb,A)
      implicit none 
      integer    nb
      complex*16 A(*)
      real*8     ar,ai

#include "C3dB.fh"
      
      integer i

!$OMP DO 
      do i=1,n2ft3d_map(nb)
         ar=dble(A(i))
         !ai=dimag(A(i))
         !A(i) = dcmplx(dsqrt(dsqrt(ar*ar + ai*ai)),0.0d0)
         A(i) = dcmplx(dsqrt(dabs(ar)),0.0d0)
      end do
!$OMP END DO 
      return
      end 



*     ***********************************
*     *                                 *
*     *        C3dB_cr_real             *
*     *                                 *
*     ***********************************

*  This routine performs the operation  C = real(A)

      subroutine C3dB_cr_real(nb,A,C)
      implicit none
      integer    nb
      complex*16 A(*)
      real*8     C(*)

#include "C3dB.fh"

      integer i

!$OMP DO 
      do i=1,n2ft3d_map(nb)
         C(i) = dble(A(i))
      end do
!$OMP END DO 
      return
      end



*     ***********************************
*     *                                 *
*     *        C3dB_cr_imag             *
*     *                                 *
*     ***********************************

*  This routine performs the operation  C = imag(A)

      subroutine C3dB_cr_imag(nb,A,C)
      implicit none
      integer    nb
      complex*16 A(*)
      real*8     C(*)

#include "C3dB.fh"

      integer i

!$OMP DO 
      do i=1,n2ft3d_map(nb)
         C(i) = dimag(A(i))
      end do
!$OMP END DO 
      return
      end




*     ***********************************
*     *                                 *
*     *        C3dB_ccr_Mul             *
*     *                                 *
*     ***********************************

*  This routine performs the operation  C = dble(A * B)

      subroutine C3dB_ccr_Mul(nb,A,B,C)
      implicit none
      integer    nb
      complex*16 A(*)
      complex*16 B(*)
      real*8     C(*)

#include "C3dB.fh"

      integer i

!$OMP DO 
      do i=1,n2ft3d_map(nb)
         C(i) = dble(A(i))*dble(B(i)) + dimag(A(i))*dimag(B(i))
      end do
!$OMP END DO 
      return
      end



*     ***********************************
*     *					*
*     *	       C3dB_rr_Sqr	 	*	
*     *					*
*     ***********************************

      subroutine C3dB_rr_Sqr(nb,A,C)     
      implicit none 
      integer nb
      real*8 A(*)
      real*8 C(*)

#include "C3dB.fh"

      integer i

!$OMP DO 
      do i=1,n2ft3d_map(nb)
         C(i) = A(i)**2
      end do
!$OMP END DO 
      return
      end

*     ***********************************
*     *					*
*     *	       C3dB_rr_Sqrt	 	*	
*     *					*
*     ***********************************

      subroutine C3dB_rr_Sqrt(nb,A,C)     
      implicit none 
      integer nb
      real*8 A(*)
      real*8 C(*)

#include "C3dB.fh"

      integer i

!$OMP DO 
      do i=1,n2ft3d_map(nb)
         C(i) = dsqrt(A(i))
      end do
!$OMP END DO 
      return
      end




*     ***********************************
*     *					*
*     *	   C3dB_c_transpose_jk_init	*
*     *					*
*     ***********************************

      subroutine C3dB_c_transpose_jk_init(nb)
      implicit none
      integer nb

#include "bafdecls.fh"
#include "errquit.fh"
#include "C3dB.fh"

c     integer iq_to_i1((NFFT1/2+1)*NFFT2*NSLABS)
c     integer iq_to_i2((NFFT1/2+1)*NFFT2*NSLABS)
c     integer i1_start(NFFT3+1)
c     integer i2_start(NFFT3+1)
      integer iq_to_i1(2,NBLOCKS)
      integer iq_to_i2(2,NBLOCKS)
      integer i1_start(2,NBLOCKS)
      integer i2_start(2,NBLOCKS)
      common / c_trans_blk / iq_to_i1,iq_to_i2,i1_start,i2_start



*     **** local variables ****
      integer proc_to,proc_from
      integer pto,qto,np,taskid
      integer pfrom,qfrom
      integer phere,qhere
      integer index1,index2,itmp
      integer i,j,k,it
      logical value

*     **** external functions ****

*     **** allocate c_trans_blk common block ****
      value = BA_alloc_get(mt_int,(nx(nb)*ny(nb)*nq(nb)),
     >                     'iq_to_i1',iq_to_i1(2,nb),iq_to_i1(1,nb))
      value = BA_alloc_get(mt_int,(nx(nb)*ny(nb)*nq(nb)),
     >                     'iq_to_i2',iq_to_i2(2,nb),iq_to_i2(1,nb))

      value = BA_alloc_get(mt_int,(nz(nb)+1),
     >                     'i1_start',i1_start(2,nb),i1_start(1,nb))
      value = BA_alloc_get(mt_int,(nz(nb)+1),
     >                     'i2_start',i2_start(2,nb),i2_start(1,nb))

      call Parallel3d_taskid_i(taskid)
      call Parallel3d_np_i(np)

      index1 = 1 
      index2 = 1
      do it=0,np-1
         proc_to   = mod(taskid+it,np)
         proc_from = mod(taskid-it+np,np)
         int_mb(i1_start(1,nb)+it) = index1
         int_mb(i2_start(1,nb)+it) = index2

         do k=1,nz(nb)
         do j=1,ny(nb)

*           **** packing scheme **** 
            call C3dB_ktoqp(nb,k,qhere,phere)
            call C3dB_ktoqp(nb,j,qto,pto)
            if ((phere.eq.taskid).and.(pto.eq.proc_to)) then
               do i=1,nx(nb)
                  itmp = i + (j-1)*nx(nb) 
     >                     + (qhere-1)*nx(nb)*ny(nb)
                  int_mb(iq_to_i1(1,nb)+itmp-1) = index1
                  index1 = index1 + 1
               end do
            end if
             
*           **** unpacking scheme ****
            call C3dB_ktoqp(nb,j,qhere,phere)
            call C3dB_ktoqp(nb,k,qfrom,pfrom)
            if ((phere.eq.taskid).and.(pfrom.eq.proc_from)) then
               do i=1,nx(nb)
                  itmp = i + (k-1)*nx(nb) 
     >                     + (qhere-1)*nx(nb)*ny(nb)
                  int_mb(iq_to_i2(1,nb)+itmp-1) = index2
                  index2 = index2 + 1
               end do
            end if
         end do
         end do
      end do
      int_mb(i1_start(1,nb)+np) = index1
      int_mb(i2_start(1,nb)+np) = index2

      return
      end


*     *******************************************
*     *						*
*     *	       C3dB_r_FormatWrite_reverse	*	
*     *						*
*     *******************************************

      subroutine C3dB_r_FormatWrite_reverse(nb,iunit,A,tmp)
      implicit none 
      integer nb
      integer iunit
      real*8     A(*)
      real*8     tmp(*)

#include "C3dB.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"

      integer rcv_len,rcv_proc

      
*     *** local variables ***
      integer MASTER,taskid
      parameter(MASTER=0)
      integer p_from, p_here,q
      integer i,j,k,index
      integer dest,source,status,msglen,idum

      call Parallel3d_taskid_i(taskid)

      !**********************
      !**** slab mapping ****
      !**********************
      if (mapping.eq.1) then

*     **** master node reads from file and distributes ****
      if (taskid.eq.MASTER) then
         do i=1,nx(nb)
         do j=1,ny(nb)

            do k=1,nz(nb)
              call C3dB_ktoqp(nb,k,q,p_from)
              if (p_from.eq.MASTER) then
                 index = (q-1)*(nx(nb))*ny(nb) 
     >                 + (j-1)*(nx(nb)) + i
                 tmp(k) = A(index)
              else
                 msglen  = 1
                 status  = msglen
                 source  = p_from
                 idum = -999

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
              call C3dB_ktoqp(nb,k,q,p_here)
              if (p_here.eq.taskid) then

                 index = (q-1)*(nx(nb))*ny(nb) 
     >                 + (j-1)*(nx(nb)) + i
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
              call C3dB_ijktoindex2p(nb,i,j,k,index,p_from)
              if (p_from.eq.MASTER) then
                 tmp(k) = A(index)
              else
                 msglen  = 1
                 status  = msglen
                 source  = p_from
                 idum = -999

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
              call C3dB_ijktoindex2p(nb,i,j,k,index,p_here)

              if (p_here.eq.taskid) then

                 tmp(1) = A(index)

                 msglen  = 1
                 dest    = MASTER
                 rcv_proc = MASTER
                 rcv_len  = mitob(1)
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
*     *	         C3dB_cc_dot  	 	*	
*     *					*
*     ***********************************

      subroutine C3dB_cc_dot(nb,A,B,sumall)     
      implicit none
      integer nb 
      real*8 A(*)
      real*8 B(*)
      real*8     sumall

#include "C3dB.fh"

      integer np
      real*8  sum


*     **** external functions ****
      real*8 ddot
      external ddot

      call nwpw_timing_start(2)

      call Parallel3d_np_i(np)

*     **** sum up dot product on this node ****
      sum = ddot(nfft3d_map(nb),A(1),2,B(1),2)
     >    + ddot(nfft3d_map(nb),A(2),2,B(2),2)



*     **** add up sums from other nodes ****
      if (np.gt.1) then
         call C3dB_SumAll(sum)
      end if

      call nwpw_timing_end(2)

      sumall = sum
      return
      end

*     ***********************************
*     *					*
*     *	         C3dB_cc_idot  	 	*	
*     *					*
*     ***********************************

      subroutine C3dB_cc_idot(nb,A,B,sumall)     
      implicit none
      integer nb 
      real*8 A(*)
      real*8 B(*)
      real*8     sumall

#include "C3dB.fh"

      real*8  sum


      real*8   ddot
      external ddot

      call nwpw_timing_start(2)

*     **** sum up dot product on this node ****
      sum = ddot(nfft3d_map(nb),A(1),2,B(1),2)
     >    + ddot(nfft3d_map(nb),A(2),2,B(2),2)

*     **** do not add up sums from other nodes ****       
      call nwpw_timing_end(2)

      sumall = sum
      return
      end


      subroutine C3dB_bb_idot(nb,A,B,sumall)
      implicit none
      integer nb
      real*8 A(*)
      real*8 B(*)
      real*8     sumall

#include "C3dB.fh"   
               
      real*8  sum     

               
      real*8   ddot
      external ddot

      call nwpw_timing_start(2)

*     **** sum up dot product on this node ****
      sum = ddot(n2ft3d_map(nb),A(1),2,B(1),2)
     >    + ddot(n2ft3d_map(nb),A(2),2,B(2),2)
         
*     **** do not add up sums from other nodes ****       
      call nwpw_timing_end(2)

      sumall = sum              
      return
      end



*     ***********************************
*     *					*
*     *	         C3dB_rr_dot  	 	*	
*     *					*
*     ***********************************

      subroutine C3dB_rr_dot(nb,A,B,sumall)     
      implicit none 
      integer nb
      real*8  A(*)
      real*8  B(*)
      real*8  sumall

#include "C3dB.fh"

      integer np
      real*8  sum

      real*8   ddot
      external ddot

      call Parallel3d_np_i(np)

*     **** sum up dot product on this node ****
      sum = ddot(n2ft3d_map(nb),A,1,B,1)

*     **** add up sums from other nodes ****
      if (np.gt.1) then
         call C3dB_SumAll(sum)
      end if

      sumall = sum
      return
      end


*     ***********************************
*     *                                 *
*     *          C3dB_tt_dot            *
*     *                                 *
*     ***********************************

      subroutine C3dB_tt_dot(nb,A,B,sumall)
      implicit none
      integer nb
      real*8  A(*)
      real*8  B(*)
      real*8  sumall

#include "C3dB.fh"

      integer np
      real*8  sum

      real*8   ddot
      external ddot

      call Parallel3d_np_i(np)

*     **** sum up dot product on this node ****
      sum = ddot(nfft3d_map(nb),A,1,B,1)

*     **** add up sums from other nodes ****
      if (np.gt.1) then
         call C3dB_SumAll(sum)
      end if

      sumall = sum
      return
      end

*     ***********************************
*     *                                 *
*     *          C3dB_tt_idot           *
*     *                                 *
*     ***********************************

      subroutine C3dB_tt_idot(nb,A,B,sumall)
      implicit none
      integer nb
      real*8  A(*)
      real*8  B(*)
      real*8  sumall

#include "C3dB.fh"

      integer np
      real*8  sum

      real*8   ddot
      external ddot

      call Parallel3d_np_i(np)

*     **** sum up dot product on this node ****
      sum = ddot(nfft3d_map(nb),A,1,B,1)

c*     **** add up sums from other nodes ****
c      if (np.gt.1) then
c         call C3dB_SumAll(sum)
c      end if

      sumall = sum
      return
      end


*     ***********************************
*     *					*
*     *	         C3dB_rr_idot  	 	*	
*     *					*
*     ***********************************

      subroutine C3dB_rr_idot(nb,A,B,sumall)     
      implicit none 
      integer nb
      real*8  A(*)
      real*8  B(*)
      real*8  sumall

#include "C3dB.fh"

      integer np
      real*8  sum

      real*8   ddot
      external ddot

c      call Parallel3d_np_i(np)

*     **** sum up dot product on this node ****
      sum = ddot(n2ft3d_map(nb),A,1,B,1)

*     **** add up sums from other nodes ****
*     if (np.gt.1) then
*        call C3dB_SumAll(sum)
*     end if

      sumall = sum
      return
      end



*     ***********************************
*     *					*
*     *	         C3dB_cc_Mul  	 	*	
*     *					*
*     ***********************************

      subroutine C3dB_cc_Mul(nb,A,B,C)     
      implicit none 
      integer    nb
      complex*16 A(*)
      complex*16 B(*)
      complex*16 C(*)

#include "C3dB.fh"

      integer i

!$OMP DO 
      do i=1,nfft3d_map(nb)
         C(i) = dconjg(A(i)) * B(i)
      end do
!$OMP END DO 

      return
      end



      subroutine C3dB_bb_Mul(nb,A,B,C)
      implicit none 
      integer    nb
      complex*16 A(*)
      complex*16 B(*)
      complex*16 C(*)
#include "C3dB.fh"
      integer i
!$OMP DO 
      do i=1,n2ft3d_map(nb)
         C(i) = dconjg(A(i)) * B(i)
      end do
!$OMP END DO 
      return
      end


      subroutine C3dB_bb_ncMul(nb,A,B,C)
      implicit none
      integer    nb
      complex*16 A(*)
      complex*16 B(*)
      complex*16 C(*)
#include "C3dB.fh"
      integer i
!$OMP DO 
      do i=1,n2ft3d_map(nb)
         C(i) = A(i) * B(i)
      end do
!$OMP END DO 
      return
      end                  



*     ***********************************
*     *                                 *
*     *          C3dB_cc_Mul2           *
*     *                                 *
*     ***********************************

      subroutine C3dB_cc_Mul2(nb,A,B)
      implicit none
      integer    nb
      complex*16 A(*)
      complex*16 B(*)

#include "C3dB.fh"

      integer i

!$OMP DO 
      do i=1,nfft3d_map(nb)
         B(i) = dconjg(A(i)) * B(i)
      end do
!$OMP END DO 

      return
      end


*     ***********************************
*     *                                 *
*     *          C3dB_cc_Mul2c          *
*     *                                 *
*     ***********************************

      subroutine C3dB_cc_Mul2c(nb,A,B)
      implicit none
      integer    nb
      complex*16 A(*)
      complex*16 B(*)

#include "C3dB.fh"

      integer i

!$OMP DO 
      do i=1,nfft3d_map(nb)
         B(i) = dconjg(B(i)) * A(i)
      end do
!$OMP END DO 

      return
      end




      subroutine C3dB_bb_Mul2c(nb,A,B)
      implicit none          
      integer    nb          
      complex*16 A(*)        
      complex*16 B(*)

#include "C3dB.fh"
            
      integer i    

!$OMP DO 
      do i=1,n2ft3d_map(nb)
         B(i) = dconjg(B(i)) * A(i)
      end do
!$OMP END DO 
            
      return
      end


*     ***********************************
*     *					*
*     *	         C3dB_lc_Mask  	 	*	
*     *					*
*     ***********************************

      subroutine C3dB_lc_Mask(nb,masker,A)
      implicit none 
      integer    nb
      logical    masker(*)
      complex*16 A(*)

#include "C3dB.fh"

      integer i

!$OMP DO 
      do i=1,nfft3d_map(nb)
         if (masker(i)) A(i) = dcmplx(0.0d0,0.0d0)
      end do
!$OMP END DO 
      return
      end

*     ***********************************
*     *					*
*     *	         C3dB_lr_Mask  	 	*	
*     *					*
*     ***********************************

      subroutine C3dB_lr_Mask(nb,masker,A)
      implicit none 
      integer   nb
      logical   masker(*)
      real*8    A(*)

#include "C3dB.fh"

      integer i

!$OMP DO 
      do i=1,nfft3d_map(nb)
         if (masker(i)) A(i) = 0.0d0
      end do
!$OMP END DO 
      return
      end


*     ***********************************
*     *					*
*     *	         C3dB_rc_Mul  	 	*	
*     *					*
*     ***********************************

      subroutine C3dB_rc_Mul(nb,A,B,C)     
      implicit none 
      integer    nb
      real*8     A(*)
      complex*16 B(*)
      complex*16 C(*)

#include "C3dB.fh"

      integer i

!$OMP DO 
      do i=1,nfft3d_map(nb)
         C(i) = A(i) * B(i)
      end do
!$OMP END DO 

      return
      end

*     ***********************************
*     *					*
*     *	         C3dB_rr_Mul  	 	*	
*     *					*
*     ***********************************

      subroutine C3dB_rr_Mul(nb,A,B,C)     
      implicit none 
      integer nb
      real*8 A(*)
      real*8 B(*)
      real*8 C(*)

#include "C3dB.fh"

      integer i

!$OMP DO 
      do i=1,n2ft3d_map(nb)
         C(i) = A(i) * B(i)
      end do
!$OMP END DO 

      return
      end

*     ***********************************
*     *                                 *
*     *          C3dB_rr_Mul2           *
*     *                                 *
*     ***********************************

      subroutine C3dB_rr_Mul2(nb,A,B)
      implicit none
      integer nb
      real*8 A(*)
      real*8 B(*)

#include "C3dB.fh"

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
*     *	         C3dB_cc_Sum  	 	*	
*     *					*
*     ***********************************

      subroutine C3dB_cc_Sum(nb,A,B,C)     
      implicit none 
      integer    nb
      complex*16 A(*)
      complex*16 B(*)
      complex*16 C(*)

#include "C3dB.fh"

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
*     *          C3dB_cc_Sum2            *
*     *                                 *
*     ***********************************

      subroutine C3dB_cc_Sum2(nb,A,B)
      implicit none
      integer    nb
      complex*16 A(*)
      complex*16 B(*)

#include "C3dB.fh"

      integer i

!$OMP DO 
      do i=1,nfft3d_map(nb)
         B(i) = B(i) + A(i)
      end do
!$OMP END DO 

      return
      end

           
      subroutine C3dB_bb_Sum2(nb,A,B)
      implicit none
      integer    nb
      complex*16 A(*)
      complex*16 B(*)
            
#include "C3dB.fh"
            
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
*     *	         C3dB_rc_Sum  	 	*	
*     *					*
*     ***********************************

      subroutine C3dB_rc_Sum(nb,A,B,C)     
      implicit none 
      integer    nb
      real*8     A(*)
      complex*16 B(*)
      complex*16 C(*)

#include "C3dB.fh"

      integer i

!$OMP DO 
      do i=1,n2ft3d_map(nb)
         C(i) = A(i) + B(i)
      end do
!$OMP END DO 

      return
      end

*     ***********************************
*     *                                 *
*     *          C3dB_rc_Sum2           *
*     *                                 *
*     ***********************************

      subroutine C3dB_rc_Sum2(nb,A,B)
      implicit none
      integer    nb
      real*8     A(*)
      complex*16 B(*)

#include "C3dB.fh"

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
*     *	         C3dB_rr_Sum  	 	*	
*     *					*
*     ***********************************

      subroutine C3dB_rr_Sum(nb,A,B,C)     
      implicit none 
      integer nb
      real*8  A(*)
      real*8  B(*)
      real*8  C(*)

#include "C3dB.fh"

      integer i

!$OMP DO 
      do i=1,n2ft3d_map(nb)
         C(i) = A(i) + B(i)
      end do
!$OMP END DO 

      return
      end


*     ***********************************
*     *                                 *
*     *          C3dB_rr_Sum2           *
*     *                                 *
*     ***********************************

      subroutine C3dB_rr_Sum2(nb,A,B)
      implicit none
      integer nb
      real*8  A(*)
      real*8  B(*)

#include "C3dB.fh"

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
*     *	         C3dB_rrc_Sum  	 	*	
*     *					*
*     ***********************************

      subroutine C3dB_rrc_Sum(nb,A,B,C)     
      implicit none 
      integer    nb
      real*8     A(*)
      real*8     B(*)
      complex*16 C(*)

#include "C3dB.fh"

      integer i

!$OMP DO 
      do i=1,n2ft3d_map(nb)
         C(i) = dcmplx((A(i) + B(i)),0.0d0)
      end do
!$OMP END DO 

      return
      end


*     ***********************************
*     *                                 *
*     *          C3dB_cc_Sub2           *       
*     *                                 *
*     ***********************************

      subroutine C3dB_cc_Sub2(nb,B,C)
      implicit none 
      integer    nb
      complex*16 B(*)
      complex*16 C(*)

#include "C3dB.fh"

      integer i

      do i=1,nfft3d_map(nb)
         C(i) = C(i) - B(i)
      end do

      return
      end

      subroutine C3dB_bb_Sub2(nb,B,C)
      implicit none
      integer    nb         
      complex*16 B(*)
      complex*16 C(*)

#include "C3dB.fh"

      integer i

!$OMP DO 
      do i=1,n2ft3d_map(nb)
         C(i) = C(i) - B(i) 
      end do                           
!$OMP END DO 
      return
      end




*     ***********************************
*     *					*
*     *	         C3dB_cc_Sub  	 	*	
*     *					*
*     ***********************************

      subroutine C3dB_cc_Sub(nb,A,B,C)     
      implicit none 
      integer    nb
      complex*16 A(*)
      complex*16 B(*)
      complex*16 C(*)

#include "C3dB.fh"

      integer i

!$OMP DO 
      do i=1,nfft3d_map(nb)
         C(i) = A(i) - B(i)
      end do
!$OMP END DO 

      return
      end


*     ***********************************
*     *					*
*     *	         C3dB_rr_Sub  	 	*	
*     *					*
*     ***********************************

      subroutine C3dB_rr_Sub(nb,A,B,C)     
      implicit none 
      integer nb
      real*8  A(*)
      real*8  B(*)
      real*8  C(*)

#include "C3dB.fh"

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
*     *	         C3dB_cc_zaxpy 	 	*	
*     *					*
*     ***********************************

      subroutine C3dB_cc_zaxpy(nb,alpha,A,B)     
      implicit none 
      integer    nb
      complex*16 alpha
      complex*16 A(*)
      complex*16 B(*)

#include "C3dB.fh"

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
*     *	         C3dB_cc_daxpy 	 	*	
*     *					*
*     ***********************************

      subroutine C3dB_cc_daxpy(nb,alpha,A,B)     
      implicit none 
      integer    nb
      real*8     alpha
      complex*16 A(*)
      complex*16 B(*)

#include "C3dB.fh"

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
*     *	         C3dB_rr_daxpy 	 	*	
*     *					*
*     ***********************************

      subroutine C3dB_rr_daxpy(nb,alpha,A,B)     
      implicit none 
      integer nb
      real*8  alpha
      real*8  A(*)
      real*8  B(*)

#include "C3dB.fh"

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
*     *          C3dB_rr_Divide         *
*     *                                 *
*     ***********************************

      subroutine C3dB_rr_Divide(nb,A,B,C)
      implicit none
      integer nb
      real*8 A(*)
      real*8 B(*)
      real*8 C(*)

#include "C3dB.fh"

      real*8 eta
      parameter (eta=1.0d-9)

      integer i


!$OMP DO 
      do i=1,n2ft3d_map(nb)
         if (dabs(B(i)) .le. eta) then
           C(i) = 0.0d0
         else
           C(i) = A(i) / B(i)
         end if
      end do
!$OMP END DO 

      return
      end


*     ***********************************
*     *                                 *
*     *          C3dB_cc_Divide2        *
*     *                                 *
*     ***********************************

      subroutine C3dB_cc_Divide2(nb,A,B)
      implicit none
      integer nb
      complex*16 A(*)
      complex*16 B(*)

#include "C3dB.fh"

      real*8 eta,ar,br
      parameter (eta=1.0d-9)

      integer i

!$OMP DO 
      do i=1,nfft3d_map(nb)
         ar = dble(A(i))
         br = dble(B(i))
         if (dabs(ar) .le. eta) then
           B(i) = dcmplx(0.0d0,0.0d0)
         else
           B(i) = dcmplx(br/ar,0.0d0)
         end if
      end do
!$OMP END DO 

      return
      end


*     ***********************************
*     *                                 *
*     *          C3dB_r_abs1             *
*     *                                 *
*     ***********************************
      subroutine C3dB_r_abs1(nb,A)
      implicit none
      integer nb
      real*8 A(*)

#include "C3dB.fh"

      integer i

!$OMP DO 
      do i=1,n2ft3d_map(nb)
         A(i) = dabs(A(i))
      end do
!$OMP END DO 

      return
      end


*     ***********************************
*     *                                 *
*     *         C3dB_c_abs1             *
*     *                                 *
*     ***********************************
      subroutine C3dB_c_abs1(nb,A)
      implicit none
      integer nb
      complex*16 A(*)

#include "C3dB.fh"

      integer i

!$OMP DO 
      do i=1,nfft3d_map(nb)
         A(i) = dcmplx(dabs(dble(A(i))),dabs(dimag(A(i))))
      end do
!$OMP END DO 

      return
      end


*     ***********************************
*     *                                 *
*     *          C3dB_r_thresh          *
*     *                                 *
*     ***********************************

      subroutine C3dB_r_thresh(nb,eta,A)
      implicit none
      integer nb
      real*8 eta
      real*8 A(*)

#include "C3dB.fh"

      integer i

!$OMP DO 
      do i=1,n2ft3d_map(nb)
         if (dabs(A(i)).lt.eta) A(i) = 0.0d0
      end do
!$OMP END DO 

      return
      end



*     ***********************************
*     *                                 *
*     *          C3dB_rr_Minus          *
*     *                                 *
*     ***********************************
      subroutine C3dB_rr_Minus(nb,A,B,C)
      implicit none
      integer nb
      real*8  A(*)
      real*8  B(*)
      real*8  C(*)

#include "C3dB.fh"

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
*     *	         C3dB_r_dsum  	 	*	
*     *					*
*     ***********************************

      subroutine C3dB_r_dsum(nb,A,sumall)     
      implicit none
      integer nb 
      real*8  A(*)
      real*8  sumall

#include "C3dB.fh"

      integer i,np
      real*8 sum

      call Parallel3d_np_i(np)

*     **** sum up dot product on this node ****
      sum = 0.0d0
      do i=1,n2ft3d_map(nb)
         sum = sum + A(i)
      end do

*     **** add up sums from other nodes ****
      if (np.gt.1) then
        call C3dB_SumAll(sum)
      end if

      sumall = sum

      return
      end

*     ***********************************
*     *					*
*     *	         C3dB_c_dsum  	 	*	
*     *					*
*     ***********************************

      subroutine C3dB_c_dsum(nb,A,sumall)     
      implicit none
      integer nb 
      complex*16  A(*)
      complex*16 sumall

#include "C3dB.fh"

      integer i,np
      complex*16 sum

      call Parallel3d_np_i(np)

*     **** sum up dot product on this node ****
      sum = dcmplx(0.0d0,0.0d0)
      do i=1,nfft3d_map(nb)
         sum = sum + A(i)
      end do

*     **** add up sums from other nodes ****
      if (np.gt.1) then
        call C3dB_Vector_SumAll(2,sum)
      end if

      sumall = sum

      return
      end




*     ***********************************
*     *					*
*     *	     C3dB_cc_Vector_dot 	*	
*     *					*
*     ***********************************

      subroutine C3dB_cc_Vector_dot(nb,nnfft3d,nn,ne,A,B,sumall)     
      implicit none 
      integer    nb
      integer    nnfft3d,nn,ne
      real*8 A(*)
      real*8 B(*)
      real*8     sumall(nn,nn)

#include "C3dB.fh"

      integer np
      integer n,m,shift1,shift2
      real*8  sum

      real*8   ddot
      external ddot

      call nwpw_timing_start(2)

      call Parallel3d_np_i(np)

*     **** sum up dot product on this node ****
      do n=1,ne
      do m=n,ne

        shift1 = 1 + (n-1)*nnfft3d*2
        shift2 = 1 + (m-1)*nnfft3d*2

        sum = ddot(nfft3d_map(nb),A(shift1),2,B(shift2),2)
     >      + ddot(nfft3d_map(nb),A(shift1+1),2,B(shift2+1),2)

         sumall(n,m) = sum
         sumall(m,n) = sum
      end do
      end do


*     **** add up sums from other nodes ****
      if (np.gt.1) then
         call C3dB_Vector_SumAll(nn*ne,sumall)
      end if

      call nwpw_timing_end(2)

      return
      end



*     ***********************************
*     *					*
*     *	     C3dB_cc_Vector_ndot 	*	
*     *					*
*     ***********************************

      subroutine C3dB_cc_Vector_ndot(nb,nnfft3d,ne,A,B,sumall)     
      implicit none 
      integer    nb
      integer    nnfft3d,ne
      real*8 A(*)
      real*8 B(*)
      real*8     sumall(ne)

#include "C3dB.fh"

      integer np
      integer n,shift1
      real*8  sum


      real*8   ddot
      external ddot

      call nwpw_timing_start(2)

      call Parallel3d_np_i(np)

*     **** sum up dot product on this node ****
      do n=1,ne

        shift1 = 1 + (n-1)*nnfft3d*2
        sum = ddot(nfft3d_map(nb),A(shift1),2,B(1),2)
     >      + ddot(nfft3d_map(nb),A(shift1+1),2,B(2),2)

        sumall(n) = sum
      end do

*     **** add up sums from other nodes ****
      if (np.gt.1) then
         call C3dB_Vector_SumAll(ne,sumall)
      end if

      call nwpw_timing_end(2)
      return
      end



*     ***********************************
*     *                                 *
*     *          C3dB_ic_Mul            *
*     *                                 *
*     ***********************************

      subroutine C3dB_ic_Mul(nb,A,B,C)
      implicit none
      integer    nb
      real*8     A(*)
      complex*16 B(*)
      complex*16 C(*)

#include "C3dB.fh"

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
*     *         C3dB_D3dB_r_Copy        *
*     *                                 *
*     ***********************************

      subroutine C3dB_D3dB_r_Copy(nb,A_c3db,B_d3db)
      implicit none
      integer nb
      real*8     A_c3db(*)
      real*8     B_d3db(*)

#include "C3dB.fh"

      integer q,indx1,indx2,nqq

      !**** slab mapping ****
      if (mapping.eq.1)  then
        nqq = nq(nb)*ny(nb) 
      !**** hilbert mapping ****
      else
        nqq = nq1(nb)
      end if
        
      indx1 = 1
      indx2 = 1
      do q=1,nqq
         call dcopy(nx(nb),A_c3db(indx1),1,B_d3dB(indx2),1)
         indx1 = indx1 + nx(nb)
         indx2 = indx2 + (nx(nb)+2)
      end do
      call dcopy(nqq,0.0d0,0,B_d3db(nx(nb)+1),nx(nb)+2)
      call dcopy(nqq,0.0d0,0,B_d3db(nx(nb)+2),nx(nb)+2)

      return
      end


*     ***********************************
*     *                                 *
*     *         D3dB_C3dB_r_Copy        *
*     *                                 *
*     ***********************************

      subroutine D3dB_C3dB_r_Copy(nb,A_d3db,B_c3db)
      implicit none
      integer nb
      real*8     A_d3db(*)
      real*8     B_c3db(*)

#include "C3dB.fh"

      integer q,indx1,indx2,nqq

      !**** slab mapping ****
      if (mapping.eq.1)  then
        nqq = nq(nb)*ny(nb) 
      !**** hilbert mapping ****
      else
        nqq = nq1(nb)
      end if

      indx1 = 1
      indx2 = 1
      do q=1,nqq
         call dcopy(nx(nb),A_d3db(indx1),1,B_c3dB(indx2),1)
         indx1 = indx1 + (nx(nb)+2)
         indx2 = indx2 + nx(nb)
      end do
      return
      end


      subroutine C3dB_pfft_index1_copy(n,index,a,b)
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

      subroutine C3dB_pfft_index2_copy(n,index,a,b)
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
!$OMP END DO 
      return
      end

      subroutine C3dB_pfft_index2_zero(n,index,a)
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

