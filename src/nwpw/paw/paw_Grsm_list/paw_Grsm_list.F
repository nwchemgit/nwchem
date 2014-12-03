*
*     $Id$ 
*
      subroutine paw_Grsm_list_start()
      implicit none

#include "bafdecls.fh"
#include "paw_Grsm_list_common.fh"
#include "errquit.fh"

*     **** local variables ****
      logical value
      integer i

*     **** external functions ****
      logical  control_lmbfgs_ondisk
      integer  control_lmbfgs_size,paw_psi_ne
      external control_lmbfgs_ondisk
      external control_lmbfgs_size,paw_psi_ne

      size_list = 2*control_lmbfgs_size()
      call Pack_npack(1,nsize)
      nsize = 2*nsize
      nsize = nsize*(paw_psi_ne(1)+paw_psi_ne(2))
      ondisk = control_lmbfgs_ondisk()

      if (.not.ondisk) then
         value = .true.
         do i=1,size_list
            value = value.and.
     >               BA_alloc_get(mt_dbl,nsize,
     >                'Grsm_list1',grsm_mem(2,i),grsm_mem(1,i))
         end do
        if (.not. value) call errquit('pspw_Grsm_list_start:get heap',0,
     &       MA_ERR)
         
      end if

      return
      end

      subroutine paw_Grsm_list_end()
      implicit none

#include "bafdecls.fh"
#include "paw_Grsm_list_common.fh"
#include "errquit.fh"

*     **** local variables ****
      logical value
      integer i

      if (.not.ondisk) then
         value = .true.
         do i=1,size_list
            value = value.and.
     >               BA_free_heap(grsm_mem(2,i))
         end do
        if (.not. value) call errquit('paw_Grsm_list_end:free heap',0,
     &       MA_ERR)

      end if

      return
      end


      subroutine paw_Grsm_list_init(tag,size_list0,nsize0)
      implicit none
      character*(*) tag
      integer size_list0,nsize0

#include "paw_Grsm_list_common.fh"

*     **** local variables ****
      integer i,l

      size_list = size_list0
      nsize     = nsize0

      l = index(tag,' ') - 1
      do i=1,size_list
        indx(i) = i
        tag_list(i)  = tag//'1'//CHAR(ICHAR('a')+i-1)
      end do
 
      return
      end

*     ************************************
*     *                                  *
*     *        paw_Grsm_list_ptr        *
*     *                                  *
*     ************************************
*
*  Warning - should not be used if Grassmann
*          list is stored on disk.
*
      subroutine paw_Grsm_list_ptr(m,ptr)
      implicit none
      integer m
      integer ptr

#include "bafdecls.fh"
#include "paw_Grsm_list_common.fh"

      ptr = grsm_mem(1,indx(m))
      return
      end


*     ************************************
*     *                                  *
*     *        paw_Grsm_list_load       *
*     *                                  *
*     ************************************
      subroutine paw_Grsm_list_load(m,A)
      implicit none
      integer m
      real*8 A(*)

#include "bafdecls.fh"
#include "paw_Grsm_list_common.fh"

      if (ondisk) then
         call nwpw_scratch_read(tag_list(indx(m)), nsize,A)
      else
         call dcopy(nsize,dbl_mb(grsm_mem(1,indx(m))),1,A,1)
      end if
      return
      end


*     ************************************
*     *                                  *
*     *        paw_Grsm_list_store      *
*     *                                  *
*     ************************************
      subroutine paw_Grsm_list_store(m,A)
      implicit none
      integer m
      real*8 A(*)

#include "bafdecls.fh"
#include "paw_Grsm_list_common.fh"

      if (ondisk) then
        call nwpw_scratch_write(tag_list(indx(m)), nsize,A)
      else
         call dcopy(nsize,A,1,dbl_mb(grsm_mem(1,indx(m))),1)
      end if
      return
      end


*     ************************************
*     *                                  *
*     *        paw_Grsm_list_shift      *
*     *                                  *
*     ************************************
      subroutine paw_Grsm_list_shift()
      implicit none

#include "paw_Grsm_list_common.fh"

*     **** local variables ****
      integer i,tmp
      
      tmp = indx(1)
      do i=1,size_list-1
         indx(i) = indx(i+1)
      end do
      indx(size_list) = tmp
      
      return
      end
