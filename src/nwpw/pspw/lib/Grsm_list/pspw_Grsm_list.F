*
*     $Id: pspw_Grsm_list.F,v 1.3 2001-12-25 22:54:10 bylaska Exp $                       
*

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

#include "pspw_Grsm_list_common.fh"

      call nwpw_scratch_read(tag_list(indx(m)), nsize,A)
      call nwpw_scratch_read(tag_list2(indx(m)),nsize,B)
      return
      end


      subroutine pspw_Grsm_list_store(m,A,B)
      implicit none
      integer m
      real*8 A(*),B(*)

#include "pspw_Grsm_list_common.fh"

      call nwpw_scratch_write(tag_list(indx(m)), nsize,A)
      call nwpw_scratch_write(tag_list2(indx(m)),nsize,B)
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
