
#define NBLOCKS 4


*
* $Id$
*

*     ***********************************************************
*     *								*
*     *   		   D1dB library				*
*     *			(MPI implemenation)			*
*     *								*
*     *   Author - Eric Bylaska					*
*     *   date   - 3/23/96					*
*     *								*
*     ***********************************************************

*	The D1dB (distributed one-dimensional block) library is to 
* be used for handling three kinds of data structures.  The first 
* data structure, denoted by "r", is a double precision array of
* length (nx(nb)+2)*ny(nb)*nz.  The second data structure, denoted by "c", is
* a double complex array of length of (nx(nb)/2+1)*ny(nb)*nz.  The third data

* (nx(nb)/2+1)*ny(nb)*nz.
*
*	The data structures is distributed across threads, p, in
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
*       integer nq,nz(NBLOCKS)
*	common	/ D3dB / nq,nz

*     **** local variables ****

*     ***********************************
*     *					*
*     *	       Mapping1d_init		*	
*     *					*
*     ***********************************

      subroutine Mapping1d_init(nb)
      implicit none
      integer nb
     
#include "bafdecls.fh"
#include "errquit.fh"
#include "D1dB.fh"


      integer k,q,p
*     integer kn
      integer taskid
      logical value
      
c      call Parallel2d_np_j(np_j)
      call Parallel2d_taskid_j(taskid)


*     **** allocate q_map,p_map,k_map
      value = BA_alloc_get(mt_int,nz(nb),'q_map',q_map(2,nb),
     >                                       q_map(1,nb))
      value = value.and.BA_alloc_get(mt_int,nz(nb),'p_map',p_map(2,nb),
     >                                       p_map(1,nb))
      value = value.and.BA_alloc_get(mt_int,nz(nb),'k_map',k_map(2,nb),
     >                                       k_map(1,nb))
      value = value.and.BA_alloc_get(mt_int,np_j,'nqarray',
     >                               nqarray(2,nb),nqarray(1,nb))
      if (.not. value)
     > call errquit('Mapping_init:out of heap memory',0, MA_ERR)


*     ****************************
*     ****** Cyclic Mapping ******
*     ****************************
      if (mapping.eq.0) then
         p = 0
         q = 1
         do k=1,nz(nb)
            int_mb(q_map(1,nb)+k-1) = q
            int_mb(p_map(1,nb)+k-1) = p
            if (p .eq. taskid) nq(nb) = q

            p = p+1
            if (p .ge. np_j) then
               p = 0
               q = q + 1
            end if
         end do

*     ***********************
*     **** block mapping ****
*     ***********************
      else

         do p=0,np_j-1
           int_mb(nqarray(1,nb)+p) = 0
         end do

         p = 0
         do k=1,nz(nb)
           int_mb(nqarray(1,nb)+p) = int_mb(nqarray(1,nb)+p) + 1
           p = mod(p+1,np_j)
         end do

         k = 1
         do p=0,np_j-1
           do q=1,int_mb(nqarray(1,nb)+p)
              int_mb(q_map(1,nb)+k-1) = q
              int_mb(p_map(1,nb)+k-1) = p
              k = k + 1
           end do
         end do

         nq(nb) = int_mb(nqarray(1,nb)+taskid)

      end if
      !*** probably not needed!! ****
      do k=1,nz(nb)
         if (int_mb(p_map(1,nb)+k-1) .eq. taskid) then
            int_mb(k_map(1,nb)+int_mb(q_map(1,nb)+k-1)-1) = k
         end if
      end do

      return
      end

*     ***********************************
*     *					*
*     *	          D1dB_end   		*	
*     *					*
*     ***********************************
      subroutine D1dB_end(nb)
      implicit none
      integer nb

#include "bafdecls.fh"
#include "errquit.fh"
#include "D1dB.fh"


      logical value

      value =.true.
      value = value.and.BA_free_heap(q_map(2,nb))
      value = value.and.BA_free_heap(p_map(2,nb))
      value = value.and.BA_free_heap(k_map(2,nb))
      value = value.and.BA_free_heap(nqarray(2,nb))
      
      if (.not. value)
     > call errquit('D1dB_end:freeing heap memory',0, MA_ERR)

      return
      end
         
*     ***********************************
*     *					*
*     *	          D1dB_qtok   		*	
*     *					*
*     ***********************************

      subroutine D1dB_qtok(nb,q,k)      
      implicit none
      integer nb
      integer q,k

#include "bafdecls.fh"
#include "D1dB.fh"

      k = int_mb(k_map(1,nb)+q-1)
      return
      end

*     ***********************************
*     *					*
*     *	          D1dB_ktoqp  		*	
*     *					*
*     ***********************************

      subroutine D1dB_ktoqp(nb,k,q,p)      
      implicit none
      integer nb
      integer k,q,p

#include "bafdecls.fh"
#include "D1dB.fh"

      q = int_mb(q_map(1,nb)+k-1)
      p = int_mb(p_map(1,nb)+k-1)
      return
      end



*     ***********************************
*     *					*
*     *	        D1dB_nq			*	
*     *					*
*     ***********************************

      subroutine D1dB_nq(nb,nqtmp)
      implicit none
      integer nb
      integer nqtmp

#include "D1dB.fh"


      nqtmp = nq(nb)
      return 
      end


*     ***********************************
*     *					*
*     *	        D1dB_nz			*	
*     *					*
*     ***********************************

      subroutine D1dB_nz(nb,nztmp)
      implicit none
      integer nb
      integer nztmp

#include "D1dB.fh"
    
      nztmp = nz(nb)
      return
      end



*     ***********************************
*     *					*
*     *	        D1dB_init		*	
*     *					*
*     ***********************************

      subroutine D1dB_init(nb,nz_in,map_in)
      implicit none
      integer nb
      integer nz_in
      integer map_in

#include "D1dB.fh"

      !**** local variables ****
      integer MASTER
      parameter (MASTER=0)
      integer taskid

      call Parallel2d_np_j(np_j)
      call Parallel_taskid(taskid)


      !**** Make sure ngrid is consistent with mapping ***
      if (np_j.gt.nz_in) then
          if (taskid.eq.MASTER) then
            write(6,*) 'Error: for 1d decomposition the',
     >                 ' number of columns (np_j processors) must ',
     >                 ' be in the range ( 1 ...nz=',
     >                   nz_in,')'
          end if
          call errquit('D1dB_init: mapping error',0,0)
      end if

*     ***** initialize D1dB common block *****
      nz(nb)     = nz_in
      mapping    = map_in

*     **** do other initializations ****
      call Mapping1d_init(nb)
      
      return
      end



