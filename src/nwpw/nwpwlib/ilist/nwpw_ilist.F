*
* $Id$
*
#define NIL	-1

*     ***********************************************************
*     *															*
*     *   		   nwpw_ilist library							*
*     *															*
*     *   Author - Eric Bylaska									*
*     *   date   - 10/09/01										*
*     *															*
*     ***********************************************************

*	The nwpw_ilist  library is a linked list for integer values.
*
*

*     ***********************************
*     *									*
*     *	       nwpw_ilist_init			*	
*     *									*
*     ***********************************
      subroutine nwpw_ilist_init(mem)
      implicit none
      integer mem(2)


*     **** initialize list ****
!!$OMP SINGLE
      mem(1) = NIL
      mem(2) = NIL
!!OMP END SINGLE copyprivate(mem)
!OMP BARRIER

      return
      end


*     ***********************************
*     *									*
*     *	       nwpw_ilist_add			*	
*     *									*
*     ***********************************

      subroutine nwpw_ilist_add(mem,tag,ivalue)
      implicit none
#include "errquit.fh"
      integer mem(2)
      integer tag
      integer ivalue
     
#include "bafdecls.fh"


*     **** local variables ****
      logical value
      integer node(2),cur
      character*16 id
      
*     **** allocate q_map,p_map,k_map *****
      id='ilist'//CHAR(ICHAR('a')+tag)
      value = BA_alloc_get(mt_int,4,id,node(2),node(1))
      if (.not. value) then
        call errquit('nwpw_ilist_add: out of heap memory',0, MA_ERR)
      end if

      int_mb(node(1))   = NIL
      int_mb(node(1)+1) = NIL
      int_mb(node(1)+2) = tag
      int_mb(node(1)+3) = ivalue


*     *** add to the end of list ***
      if (mem(1).eq.NIL) then
         mem(1) = node(1)
         mem(2) = node(2)
      else
        cur = mem(1)
        do while (int_mb(cur).ne.NIL)
          cur = int_mb(cur)
        end do
        int_mb(cur)   = node(1)
        int_mb(cur+1) = node(2)
      end if

      return
      end

*     ***********************************
*     *									*
*     *	        nwpw_ilist_get   		*	
*     *									*
*     ***********************************
      subroutine nwpw_ilist_get(mem,tag,ivalue)
      implicit none
#include "errquit.fh"
      integer mem(2)
      integer tag
      integer ivalue

#include "bafdecls.fh"

*     **** local variables ****
      integer cur

      cur = mem(1)
      if (cur.eq.NIL) then
         call errquit('nwpw_ilist_get: empty list',0, MEM_ERR)
      end if

      do while ((int_mb(cur+2).ne.tag).and.
     >          (int_mb(cur).ne.NIL))
        cur = int_mb(cur)
      end do

*     **** error - tag not found ****
      if  (int_mb(cur+2).ne.tag) then
         write(*,*) 'tag,mem:',tag,mem
         call errquit('nwpw_ilist_get: tag not found',0, INT_ERR)
      end if

      ivalue = int_mb(cur+3)

      return
      end

         
*     ***********************************
*     *									*
*     *	       nwpw_ilist_delete		*
*     *									*
*     ***********************************

      subroutine nwpw_ilist_delete(mem,tag)      
      implicit none
#include "errquit.fh"
      integer mem(2)
      integer tag

#include "bafdecls.fh"

*     **** local variables ****
      integer cur,prev,del_node(2),next(2)

      prev = mem(1)
      cur  = mem(1)
      if (cur.eq.NIL) then
         write(*,*) 'warning empty ilist, tag,mem:',tag,mem
         return
      end if

*     *** delete first node ****
      if (int_mb(cur+2).eq.tag) then
         del_node(1) = mem(1)
         del_node(2) = mem(2)
         mem(1) = int_mb(cur)
         mem(2) = int_mb(cur+1)

*     *** delete after first node ****
      else
         do while ((int_mb(cur+2).ne.tag).and.
     >             (int_mb(cur).ne.NIL))
           prev = cur
           cur  = int_mb(cur)
         end do

*        **** error - tag not found ****
         if  (int_mb(cur+2).ne.tag) then
            write(*,*) 'tag,mem:',tag,mem
            call errquit('nwpw_ilist_get: tag not found',0, INT_ERR)
         end if

         del_node(1) = int_mb(prev)
         del_node(2) = int_mb(prev+1)
         next(1) = int_mb(cur)
         next(2) = int_mb(cur+1)

         int_mb(prev)   = next(1)
         int_mb(prev+1) = next(2)
      end if

*     **** remove node from heap ****
      if (.not.BA_free_heap(del_node(2))) then
         call errquit('nwpw_ilist_delete: cannot free heap',0, MA_ERR)
      end if
      
      return
      end



