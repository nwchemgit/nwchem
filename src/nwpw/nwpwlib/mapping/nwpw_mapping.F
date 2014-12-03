*
* $Id$
*

*     ***********************************************************
*     *															*
*     *   		   nwpw_mapping library							*
*     *															*
*     *   Author - Eric Bylaska									*
*     *   date   - 10/09/01										*
*     *															*
*     ***********************************************************

*	The nwpw_mapping  library is to  be used for handling parallel
* mapping of three-dimensional grids {A(n1(nb),n2(nb),n3(nb))}.
*
*	The three-dimensional grids are distributed across threads, p, in
* the k {i.e. n3(nb)} dimension using a cyclic decomposition.  So that
* three-dimensional grid A is defined as  A(n1(nb),n2(nb),nq(nb)) on
* each thread.  
*
*	Where 
*		np = number of threads
*		nq(nb) = ceil(nz(nb)/np).
*		0 <= p < np
*		1 <= qq <= nq(nb)
*		1 <= k <= n3(nb)
*
*		1 <= i <= n1(nb)
*		1 <= j <= n2(nb)
*       1 <= q <= (n1(nb)*n2(nb)*nq(nb))
* 
* 	The mapping of i,j,k -> q,p is defined as:
*	
*       
*		k  = ((qq-1)*np + p) + 1
*		qq = ((k-1) - p)/np + 1
*		p  = (k-1) mod np
*
*       q  = i + (j-1)*n1(nb) + (qq-1)*n1(nb)*n2(nb)
*
*

*     ***********************************
*     *									*
*     *	       nwpw_mapping_start		*	
*     *									*
*     ***********************************
      subroutine nwpw_mapping_start()
      implicit none

***** nwpw_mapping_block ****
#include "nwpw_mapping.fh"

*     **** initialize lists ****
      call nwpw_ilist_init(mem_n1)
      call nwpw_ilist_init(mem_n2)
      call nwpw_ilist_init(mem_n3)
      call nwpw_ilist_init(mem_nq)

      call nwpw_ilist_init(mem_qmap1)
      call nwpw_ilist_init(mem_qmap2)

      call nwpw_ilist_init(mem_pmap1)
      call nwpw_ilist_init(mem_pmap2)

      call nwpw_ilist_init(mem_kmap1)
      call nwpw_ilist_init(mem_kmap2)

      return
      end


*     ***********************************
*     *									*
*     *	       nwpw_mappint_init		*	
*     *									*
*     ***********************************

      subroutine nwpw_mapping_init(nb,nx,ny,nz)
      implicit none
#include "errquit.fh"
      integer nb
      integer nx,ny,nz
     
#include "bafdecls.fh"

***** nwpw_mapping_block ****
#include "nwpw_mapping.fh"


      integer qmap(2),pmap(2),kmap(2)
      integer k,q,p
      integer taskid,np,nq
      logical value
      character*16 id
      

*     **** allocate q_map,p_map,k_map *****
      id='qmap'//CHAR(ICHAR('a')+nb)
      value = BA_alloc_get(mt_int,nz,id,qmap(2),qmap(1))

      id='pmap'//CHAR(ICHAR('a')+nb)
      value = value.and.
     >        BA_alloc_get(mt_int,nz,id,pmap(2),pmap(1))

      id='kmap'//CHAR(ICHAR('a')+nb)
      value = BA_alloc_get(mt_int,nz,id,kmap(2),kmap(1))
      if (.not. value) then
        call errquit('nwpw_mapping_init: out of heap memory',0, MA_ERR)
      end if


      call Parallel_np(np)
      call Parallel_taskid(taskid)

*     ****************************
*     ****** cyclic mapping ******
*     ****************************
      p = 0
      q = 1
      do k=1,nz
         int_mb(qmap(1)+k-1) = q
         int_mb(pmap(1)+k-1) = p
         if (p .eq. taskid) nq = q

         p        = p+1
         if (p .ge. np) then
            p = 0
            q = q + 1
         end if
      end do


      do k=1,nz
         if (int_mb(pmap(1)+k-1) .eq. taskid) then
            int_mb(kmap(1)+int_mb(qmap(1)+k-1)-1) = k
         end if
      end do



*     *** add to lists ****
      call nwpw_ilist_add(mem_n1, nb,nx)
      call nwpw_ilist_add(mem_n2, nb,ny)
      call nwpw_ilist_add(mem_n3, nb,nz)
      call nwpw_ilist_add(mem_nq, nb,nq)

      call nwpw_ilist_add(mem_qmap1, nb,qmap(1))
      call nwpw_ilist_add(mem_qmap2, nb,qmap(2))

      call nwpw_ilist_add(mem_pmap1, nb,pmap(1))
      call nwpw_ilist_add(mem_pmap2, nb,pmap(2))

      call nwpw_ilist_add(mem_kmap1, nb,kmap(1))
      call nwpw_ilist_add(mem_kmap2, nb,kmap(2))

      return
      end

*     ***********************************
*     *									*
*     *	        nwpw_mapping_end   		*	
*     *									*
*     ***********************************
      subroutine nwpw_mapping_end(nb)
      implicit none
#include "errquit.fh"
      integer nb

#include "bafdecls.fh"

***** nwpw_mapping_block ****
#include "nwpw_mapping.fh"


*     **** local variables ****
      logical value
      integer qmap2,pmap2,kmap2


*     **** free maps from heap ****
      call nwpw_ilist_get(mem_qmap2, nb,qmap2)
      call nwpw_ilist_get(mem_pmap2, nb,pmap2)
      call nwpw_ilist_get(mem_kmap2, nb,kmap2)
      value = BA_free_heap(qmap2)
      value = value.and.
     >        BA_free_heap(pmap2)
      value = value.and.
     >        BA_free_heap(kmap2)
      if (.not. value) then
        call errquit('nwpw_mapping_end: error freeing heap memory',0,
     &       MA_ERR)
      end if

*     **** delete id from lists ****
      call nwpw_ilist_delete(mem_n1, nb)
      call nwpw_ilist_delete(mem_n2, nb)
      call nwpw_ilist_delete(mem_n3, nb)
      call nwpw_ilist_delete(mem_nq, nb)

      call nwpw_ilist_delete(mem_qmap1, nb)
      call nwpw_ilist_delete(mem_qmap2, nb)

      call nwpw_ilist_delete(mem_pmap1, nb)
      call nwpw_ilist_delete(mem_pmap2, nb)

      call nwpw_ilist_delete(mem_kmap1, nb)
      call nwpw_ilist_delete(mem_kmap2, nb)

      return
      end

         
*     ***********************************
*     *									*
*     *	       nwpw_mapping_pqtoijk		*
*     *									*
*     ***********************************

      subroutine nwpw_mapping_qptoijk(nb,q,p,i,j,k)      
      implicit none
      integer nb
      integer q,p
      integer i,j,k

#include "bafdecls.fh"
      
***** nwpw_mapping_block ****
#include "nwpw_mapping.fh"
      
*     ***** local variables ****
      integer kmap
      integer n1,n2,qq,tq


      call nwpw_ilist_get(mem_kmap1, nb,kmap)
      call nwpw_ilist_get(mem_n1, nb,n1)
      call nwpw_ilist_get(mem_n2, nb,n2)

      tq = q

      qq = tq/(n1*n2) + 1
      tq = tq - (qq-1)*n1*n2

      j  = tq/n1 + 1
      i = tq - (j-1)*n1

      k = int_mb(kmap+qq-1)

      return
      end


*     ***********************************
*     *									*
*     *	        nwpw_mapping_ijktoqp	*	
*     *									*
*     ***********************************

      subroutine nwpw_mapping_ijktoqp(nb,i,j,k,q,p)      
      implicit none
      integer nb
      integer i,j,k,q,p

#include "bafdecls.fh"

***** nwpw_mapping_block ****
#include "nwpw_mapping.fh"

*     ***** local variables ****
      integer qmap,pmap
      integer n1,n2,qq

      call nwpw_ilist_get(mem_n1, nb,n1)
      call nwpw_ilist_get(mem_n2, nb,n2)
      call nwpw_ilist_get(mem_qmap1, nb,qmap)
      call nwpw_ilist_get(mem_pmap1, nb,pmap)
      
      p  = int_mb(pmap+k-1)
      qq = int_mb(qmap+k-1)

      q = i + (j-1)*n1 + (qq-1)*n1*n2

      return
      end

*     ***********************************
*     *									*
*     *	      nwpw_mapping_nfft3d		*
*     *									*
*     ***********************************

      subroutine nwpw_mapping_nfft3d(nb,nfft3d)
      implicit none
      integer nb
      integer nfft3d

#include "bafdecls.fh"

***** nwpw_mapping_block ****
#include "nwpw_mapping.fh"

*     ***** local variables ****
      integer n1,n2,nq

      call nwpw_ilist_get(mem_n1, nb,n1)
      call nwpw_ilist_get(mem_n2, nb,n2)
      call nwpw_ilist_get(mem_nq, nb,nq)

      nfft3d = n1*n2*nq
      return
      end








