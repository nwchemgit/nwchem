*
* $Id$
*

*     ***********************************************************
*     *								*
*     *   		   K1dB library				*
*     *			(MPI implemenation)			*
*     *								*
*     *   Author - Eric Bylaska					*
*     *   date   - 3/23/96					*
*     *								*
*     ***********************************************************

*	The K1dB (distributed one-dimensional block) library is to 
* be used for handling a parallel 1d data structure.  
*
*	The data structures is distributed across threads, p, in
* the k (i.e. nz) dimension using a cyclic decomposition.  
*
*	Where 
*		np = number of threads
*		nq = ceil(nz/np).
*		0 <= p < np
*		1 <= q <= nq
*		1 <= k <= nz
* 
* 	The mapping of k -> q is defined as:
*	
*		k = ((q-1)*np + p) + 1
*		q = ((k-1) - p)/np + 1
*		p = (k-1) mod np
*


*     ***********************************
*     *					*
*     *	       k_Mapping1d_init		*	
*     *					*
*     ***********************************

      subroutine k_Mapping1d_init()
      implicit none
     
#include "bafdecls.fh"
#include "errquit.fh"
#include "K1dB.fh"


      integer k,q,p,taskid
*     integer kn
      logical value

      call Parallel3d_taskid_k(taskid)

*     **** allocate q_map,p_map,k_map
      value=BA_alloc_get(mt_int,nz,'q_map',q_map(2),q_map(1))
      value=value.and.BA_alloc_get(mt_int,nz,'p_map',p_map(2),p_map(1))
      value=value.and.BA_alloc_get(mt_int,nz,'k_map',k_map(2),k_map(1))
      value=value.and.BA_alloc_get(mt_int,np_k,'nqarray',
     >                               nqarray(2),nqarray(1))
      if (.not. value)
     > call errquit('k_Mapping_init:out of heap memory',0, MA_ERR)


*     ****************************
*     ****** Cyclic Mapping ******
*     ****************************
      if (mapping.eq.0) then
         p = 0
         q = 1
         do k=1,nz
            int_mb(q_map(1)+k-1) = q
            int_mb(p_map(1)+k-1) = p
            if (p .eq. taskid) nq = q

            p = p+1
            if (p .ge. np_k) then
               p = 0
               q = q + 1
            end if
         end do

*     ***********************
*     **** block mapping ****
*     ***********************
      else

         do p=0,np_k-1
           int_mb(nqarray(1)+p) = 0
         end do

         p = 0
         do k=1,nz
           int_mb(nqarray(1)+p) = int_mb(nqarray(1)+p) + 1
           p = mod(p+1,np_k)
         end do

         k = 1
         do p=0,np_k-1
           do q=1,int_mb(nqarray(1)+p)
              int_mb(q_map(1)+k-1) = q
              int_mb(p_map(1)+k-1) = p
              k = k + 1
           end do
         end do

         nq = int_mb(nqarray(1)+taskid)

      end if
      !*** probably not needed!! ****
      do k=1,nz
         if (int_mb(p_map(1)+k-1) .eq. taskid) then
            int_mb(k_map(1)+int_mb(q_map(1)+k-1)-1) = k
         end if
      end do

      return
      end

*     ***********************************
*     *					*
*     *	          K1dB_end   		*	
*     *					*
*     ***********************************
      subroutine K1dB_end()
      implicit none

#include "bafdecls.fh"
#include "errquit.fh"
#include "K1dB.fh"


      logical value

      value =.true.
      value = value.and.BA_free_heap(q_map(2))
      value = value.and.BA_free_heap(p_map(2))
      value = value.and.BA_free_heap(k_map(2))
      value = value.and.BA_free_heap(nqarray(2))
      
      if (.not. value)
     > call errquit('K1dB_end:freeing heap memory',0, MA_ERR)

      return
      end
         
*     ***********************************
*     *					*
*     *	          K1dB_qtok   		*	
*     *					*
*     ***********************************

      subroutine K1dB_qtok(q,k)      
      implicit none
      integer q,k

#include "bafdecls.fh"
#include "K1dB.fh"

      k = int_mb(k_map(1)+q-1)
      return
      end

*     ***********************************
*     *					*
*     *	          K1dB_ktoqp  		*	
*     *					*
*     ***********************************

      subroutine K1dB_ktoqp(k,q,p)      
      implicit none
      integer k,q,p

#include "bafdecls.fh"
#include "K1dB.fh"

      q = int_mb(q_map(1)+k-1)
      p = int_mb(p_map(1)+k-1)
      return
      end



*     ***********************************
*     *					*
*     *	        K1dB_nq			*	
*     *					*
*     ***********************************

      subroutine K1dB_nq(nqtmp)
      implicit none
      integer nqtmp

#include "K1dB.fh"


      nqtmp = nq
      return 
      end


*     ***********************************
*     *					*
*     *	        K1dB_nz			*	
*     *					*
*     ***********************************

      subroutine K1dB_nz(nztmp)
      implicit none
      integer nztmp

#include "K1dB.fh"
    
      nztmp = nz
      return
      end



*     ***********************************
*     *					*
*     *	        K1dB_init		*	
*     *					*
*     ***********************************

      subroutine K1dB_init(nz_in,map_in)
      implicit none
      integer nz_in
      integer map_in

#include "K1dB.fh"

      !**** local variables ****
      integer MASTER
      parameter (MASTER=0)
      integer taskid

      call Parallel3d_np_k(np_k)
      call Parallel_taskid(taskid)

      !**** Make sure ngrid is consistent with mapping ***
      if (np_k.gt.nz_in) then
          if (taskid.eq.MASTER) then
            write(6,*) 'Possible Error?: for 1d decomposition the',
     >                 ' number of columns (np_k processors) must ',
     >                 ' be in the range ( 1 ...nz=',
     >                   nz_in,')'
            call util_flush(6)
          end if
c          call errquit('K1dB_init: mapping error',0,0)
      end if

*     ***** initialize K1dB common block *****
      nz      = nz_in
      mapping = map_in

*     **** do other initializations ****
      call k_Mapping1d_init()
      
      return
      end



