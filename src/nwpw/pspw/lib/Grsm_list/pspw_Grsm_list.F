*
*     $Id: pspw_Grsm_list.F,v 1.4 2001-12-28 01:32:03 bylaska Exp $ 
*
      subroutine pspw_Grsm_list_start()
      implicit none

#include "mafdecls.fh"
#include "pspw_Grsm_list_common.fh"

*     **** local variables ****
      logical value
      integer i

*     **** external functions ****
      logical  control_lmbfgs_ondisk
      integer  control_lmbfgs_size,psi_ne
      external control_lmbfgs_ondisk
      external control_lmbfgs_size,psi_ne

      size_list = control_lmbfgs_size()
      call Pack_npack(1,nsize)
      nsize = 2*nsize
      nsize = nsize*(psi_ne(1)+psi_ne(2))
      ondisk = control_lmbfgs_ondisk()

      if (.not.ondisk) then
         value = .true.
         do i=1,size_list
            value = value.and.
     >               MA_alloc_get(mt_dbl,nsize,
     >                'Grsm_list1',grsm_mem(2,1,i),grsm_mem(1,1,i))
            value = value.and.
     >               MA_alloc_get(mt_dbl,nsize,
     >                'Grsm_list2',grsm_mem(2,2,i),grsm_mem(1,2,i))
         end do
        if (.not. value) call errquit('pspw_Grsm_list_start:get heap',0)
         
      end if

      return
      end

      subroutine pspw_Grsm_list_end()
      implicit none

#include "mafdecls.fh"
#include "pspw_Grsm_list_common.fh"

*     **** local variables ****
      logical value
      integer i

      if (.not.ondisk) then
         value = .true.
         do i=1,size_list
            value = value.and.
     >               MA_free_heap(grsm_mem(2,1,i))
            value = value.and.
     >               MA_free_heap(grsm_mem(2,2,i))
         end do
        if (.not. value) call errquit('pspw_Grsm_list_end:free heap',0)

      end if

      return
      end


      subroutine pspw_Grsm_list_init(tag,size_list0,nsize0)
      implicit none
      character*(*) tag
      integer size_list0,nsize0

#include "pspw_Grsm_list_common.fh"

*     **** local variables ****
      integer i,l

      size_list = size_list0
      nsize     = nsize0

      l = index(tag,' ') - 1
      do i=1,size_list
        indx(i) = i
        tag_list(i)  = tag//'1'//CHAR(ICHAR('a')+i-1)
        tag_list2(i) = tag//'2'//CHAR(ICHAR('a')+i-1)
      end do
 
      return
      end

      subroutine pspw_Grsm_list_load(m,A,B)
      implicit none
      integer m
      real*8 A(*), B(*)

#include "mafdecls.fh"
#include "pspw_Grsm_list_common.fh"

      if (ondisk) then
         call nwpw_scratch_read(tag_list(indx(m)), nsize,A)
         call nwpw_scratch_read(tag_list2(indx(m)),nsize,B)
      else
         call dcopy(nsize,dbl_mb(grsm_mem(1,1,indx(m))),1,A,1)
         call dcopy(nsize,dbl_mb(grsm_mem(1,2,indx(m))),1,B,1)
      end if
      return
      end


      subroutine pspw_Grsm_list_store(m,A,B)
      implicit none
      integer m
      real*8 A(*),B(*)

#include "mafdecls.fh"
#include "pspw_Grsm_list_common.fh"

      if (ondisk) then
        call nwpw_scratch_write(tag_list(indx(m)), nsize,A)
        call nwpw_scratch_write(tag_list2(indx(m)),nsize,B)
      else
         call dcopy(nsize,A,1,dbl_mb(grsm_mem(1,1,indx(m))),1)
         call dcopy(nsize,B,1,dbl_mb(grsm_mem(1,2,indx(m))),1)
      end if
      return
      end

      subroutine pspw_Grsm_list_shift()
      implicit none

#include "pspw_Grsm_list_common.fh"

*     **** local variables ****
      integer i,tmp
      
      tmp = indx(1)
      do i=1,size_list-1
         indx(i) = indx(i+1)
      end do
      indx(size_list) = tmp
      
      return
      end
