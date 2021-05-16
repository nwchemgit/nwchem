*
*     $Id$ 
*
      subroutine nwpw_list_start(nb,size_list0,nsize0)
      implicit none
      integer nb
      integer size_list0,nsize0

#include "bafdecls.fh"
#include "nwpw_list_common.fh"

*     **** local variables ****
      logical value
      integer i

!$OMP MASTER
      size_list(nb) = size_list0
      nsize(nb) = nsize0
!$OMP END MASTER
!$OMP BARRIER

      
      value = .true.
      do i=1,size_list(nb)
!$OMP MASTER
         indx(i,nb) = i
!$OMP END MASTER
         value = value.and.
     >           BA_alloc_get(mt_dbl,nsize(nb),
     >            'Grsm_list1',grsm_mem(2,i,nb),grsm_mem(1,i,nb))
      end do
      if (.not. value) call errquit('nwpw_list_start:get heap',0,0)
         
      return
      end

      subroutine nwpw_list_end(nb)
      implicit none
      integer nb

#include "bafdecls.fh"
#include "nwpw_list_common.fh"

*     **** local variables ****
      logical value
      integer i

      
      value = .true.
      do i=1,size_list(nb)
         value = value.and.
     >           BA_free_heap(grsm_mem(2,i,nb))
      end do
      if (.not. value) call errquit('nwpw_list_end:free heap',0,0)

      return
      end


*     ************************************
*     *                                  *
*     *        nwpw_list_ptr             *
*     *                                  *
*     ************************************
*
*  Warning - should not be used if Grassmann
*          list is stored on disk.
*
      subroutine nwpw_list_ptr(nb,m,ptr)
      implicit none
      integer nb
      integer m
      integer ptr

#include "bafdecls.fh"
#include "nwpw_list_common.fh"

      ptr = grsm_mem(1,indx(m,nb),nb)
      return
      end


*     ************************************
*     *                                  *
*     *          nwpw_list_load          *
*     *                                  *
*     ************************************
      subroutine nwpw_list_load(nb,m,A)
      implicit none
      integer nb
      integer m
      real*8 A(*)

#include "bafdecls.fh"
#include "nwpw_list_common.fh"

      !call dcopy(nsize(nb),dbl_mb(grsm_mem(1,indx(m,nb),nb)),1,A,1)
      call Parallel_shared_vector_copy(.true.,nsize(nb),
     >                    dbl_mb(grsm_mem(1,indx(m,nb),nb)),A)
      
      return
      end


*     ************************************
*     *                                  *
*     *        nwpw_list_store           *
*     *                                  *
*     ************************************
      subroutine nwpw_list_store(nb,m,A)
      implicit none
      integer nb
      integer m
      real*8 A(*)

#include "bafdecls.fh"
#include "nwpw_list_common.fh"

      !call dcopy(nsize(nb),A,1,dbl_mb(grsm_mem(1,indx(m,nb),nb)),1)
      call Parallel_shared_vector_copy(.true.,nsize(nb),
     >                    A,dbl_mb(grsm_mem(1,indx(m,nb),nb)))
      return
      end


*     ************************************
*     *                                  *
*     *        nwpw_list_shift           *
*     *                                  *
*     ************************************
      subroutine nwpw_list_shift(nb)
      implicit none
      integer nb

#include "nwpw_list_common.fh"

*     **** local variables ****
      integer i,tmp
      
!$OMP MASTER
      tmp = indx(1,nb)
      do i=1,size_list(nb)-1
         indx(i,nb) = indx(i+1,nb)
      end do
      indx(size_list(nb),nb) = tmp
!$OMP END MASTER
      
      return
      end

*     ************************************
*     *                                  *
*     *        nwpw_list_shift_range     *
*     *                                  *
*     ************************************
      subroutine nwpw_list_shift_range(nb,istart,iend)
      implicit none
      integer nb
      integer istart,iend

#include "nwpw_list_common.fh"

*     **** local variables ****
      integer i,tmp

      
!$OMP MASTER
      tmp = indx(istart,nb)
      do i=istart,iend-1
         indx(i,nb) = indx(i+1,nb)
      end do
      indx(iend,nb) = tmp
!$OMP END MASTER
      
      return
      end


*     ************************************
*     *                                  *
*     *        nwpw_list_remove          *
*     *                                  *
*     ************************************
      subroutine nwpw_list_remove(nb,k)
      implicit none
      integer nb
      integer k

#include "nwpw_list_common.fh"

*     **** local variables ****
      integer i,tmp
      
!$OMP MASTER
      tmp = indx(k,nb)
      do i=k,size_list(nb)-1
         indx(i,nb) = indx(i+1,nb)
      end do
      indx(size_list(nb),nb) = tmp
!$OMP END MASTER
      
      return
      end

*     ************************************
*     *                                  *
*     *        nwpw_list_swap            *
*     *                                  *
*     ************************************
      subroutine nwpw_list_swap(nb,j,k)
      implicit none
      integer nb
      integer j,k

#include "nwpw_list_common.fh"

*     **** local variables ****
      integer tmp
      
!$OMP MASTER
      tmp     = indx(k,nb)
      indx(k,nb) = indx(j,nb)
      indx(j,nb) = tmp
!$OMP END MASTER
      
      return
      end
