*
*     $Id: nwpw_list.F,v 1.3 2003-10-24 18:48:03 marat Exp $ 
*
      subroutine nwpw_list_start(size_list0,nsize0)
      implicit none
      integer size_list0,nsize0

#include "mafdecls.fh"
#include "nwpw_list_common.fh"

*     **** local variables ****
      logical value
      integer i

      size_list = size_list0
      nsize = nsize0

      
      value = .true.
      do i=1,size_list
         indx(i) = i
         value = value.and.
     >           MA_alloc_get(mt_dbl,nsize,
     >            'Grsm_list1',grsm_mem(2,i),grsm_mem(1,i))
      end do
      if (.not. value) call errquit('nwpw_list_start:get heap',0,0)
         
      return
      end

      subroutine nwpw_list_end()
      implicit none

#include "mafdecls.fh"
#include "nwpw_list_common.fh"

*     **** local variables ****
      logical value
      integer i

      
      value = .true.
      do i=1,size_list
         value = value.and.
     >           MA_free_heap(grsm_mem(2,i))
      end do
      if (.not. value) call errquit('nwpw_list_end:free heap',0,0)

      return
      end


*     ************************************
*     *                                  *
*     *        nwpw_list_ptr        *
*     *                                  *
*     ************************************
*
*  Warning - should not be used if Grassmann
*          list is stored on disk.
*
      subroutine nwpw_list_ptr(m,ptr)
      implicit none
      integer m
      integer ptr

#include "mafdecls.fh"
#include "nwpw_list_common.fh"

      ptr = grsm_mem(1,indx(m))
      return
      end


*     ************************************
*     *                                  *
*     *          nwpw_list_load          *
*     *                                  *
*     ************************************
      subroutine nwpw_list_load(m,A)
      implicit none
      integer m
      real*8 A(*)

#include "mafdecls.fh"
#include "nwpw_list_common.fh"

      call dcopy(nsize,dbl_mb(grsm_mem(1,indx(m))),1,A,1)
      
      return
      end


*     ************************************
*     *                                  *
*     *        nwpw_list_store      *
*     *                                  *
*     ************************************
      subroutine nwpw_list_store(m,A)
      implicit none
      integer m
      real*8 A(*)

#include "mafdecls.fh"
#include "nwpw_list_common.fh"

      call dcopy(nsize,A,1,dbl_mb(grsm_mem(1,indx(m))),1)
      return
      end


*     ************************************
*     *                                  *
*     *        nwpw_list_shift           *
*     *                                  *
*     ************************************
      subroutine nwpw_list_shift()
      implicit none

#include "nwpw_list_common.fh"

*     **** local variables ****
      integer i,tmp
      
      tmp = indx(1)
      do i=1,size_list-1
         indx(i) = indx(i+1)
      end do
      indx(size_list) = tmp
      
      return
      end

*     ************************************
*     *                                  *
*     *        nwpw_list_remove          *
*     *                                  *
*     ************************************
      subroutine nwpw_list_remove(k)
      implicit none
      integer k

#include "nwpw_list_common.fh"

*     **** local variables ****
      integer i,tmp
      
      tmp = indx(k)
      do i=k,size_list-1
         indx(i) = indx(i+1)
      end do
      indx(size_list) = tmp
      
      return
      end

*     ************************************
*     *                                  *
*     *        nwpw_list_swap            *
*     *                                  *
*     ************************************
      subroutine nwpw_list_swap(j,k)
      implicit none
      integer j,k

#include "nwpw_list_common.fh"

*     **** local variables ****
      integer tmp
      
      tmp     = indx(k)
      indx(k) = indx(j)
      indx(j) = tmp

      
      return
      end
