*
* $Id$
*

*     ***********************************************************
*     *															*
*     *   		   		lcao_dnlist library						*
*     *															*
*     *   Author - Eric Bylaska									*
*     *   date   - 10/09/01										*
*     *															*
*     ***********************************************************

*
*

*     ***********************************
*     *									*
*     *	       lcao_dnlist_init			*	
*     *									*
*     ***********************************
      subroutine lcao_dnlist_init(maxsize0)
      implicit none
#include "errquit.fh"
      integer maxsize0

#include "bafdecls.fh"
#include "lcao_dnlist.fh"

*     **** local variables ****
      logical value

      maxsize = maxsize0
      size    = 0

*    ***** allocate memory from heap memory ****
      value = BA_alloc_get(mt_int,(maxsize),'dnlist',
     >                     dnlist(2),dnlist(1))
      value = value.and.
     >        BA_alloc_get(mt_int,(maxsize),'dnlist_tag',
     >                     dnlist_tag(2),dnlist_tag(1))
      if (.not. value)
     > call errquit('lcao_dnlist_init: out of heap memory',0, MA_ERR)

      return
      end


*     ***********************************
*     *									*
*     *	       lcao_dnlist_end			*	
*     *									*
*     ***********************************

      subroutine lcao_dnlist_end()
      implicit none
#include "errquit.fh"
     
#include "bafdecls.fh"
#include "lcao_dnlist.fh"


*     **** local variables ****
      logical value
      integer i,mindx

      value = .true.
      do i=1,size
         mindx = int_mb(dnlist_tag(1)+i-1)
         value = value.and.BA_free_heap(mindx)
      end do
      value = value.and.BA_free_heap(dnlist_tag(2))
      value = value.and.BA_free_heap(dnlist(2))
      if (.not. value)
     > call errquit('lcao_dnlist_init: freeing heap memory',0, MA_ERR)

      size    = 0
      maxsize = 0
      return
      end

*     ***********************************
*     *									*
*     *	        lcao_dnlist_new  		*	
*     *									*
*     ***********************************
      subroutine lcao_dnlist_new(n2ft3d,index)
      implicit none
#include "errquit.fh"
      integer n2ft3d
      integer index

#include "bafdecls.fh"
#include "lcao_dnlist.fh"

*     **** local variables ****
      logical value
      integer indx,mindx

      if (size.ge.maxsize) 
     > call errquit('lcao_dnlist_new: maxsize exceeded',0, INPUT_ERR)

      value = BA_alloc_get(mt_dbl,(n2ft3d),'dn',
     >                     mindx,indx)
      if (.not.value) 
     > call errquit('lcao_dnlist_new: out of heap memory',0, MA_ERR)

      size = size + 1
      int_mb(dnlist(1)    +size-1) = indx
      int_mb(dnlist_tag(1)+size-1) = mindx

      index  = indx
      return
      end
         


*     ***********************************
*     *									*
*     *	        lcao_dnlist_get  		*	
*     *									*
*     ***********************************
      subroutine lcao_dnlist_get(count,index)
      implicit none
#include "errquit.fh"
      integer count
      integer index

#include "bafdecls.fh"
#include "lcao_dnlist.fh"

      if ((count.gt.size).or.(count.lt.1)) 
     > call errquit('lcao_dnlist_get: bad count',0, INPUT_ERR)

      index   = int_mb(dnlist(1)+count-1)
      return
      end


*     ***********************************
*     *									*
*     *	        lcao_dnlist_size 		*	
*     *									*
*     ***********************************
      integer function lcao_dnlist_size()
      implicit none

#include "lcao_dnlist.fh"

      lcao_dnlist_size = size
      return
      end

*     ***********************************
*     *									*
*     *	        lcao_dnlist_maxsize 	*	
*     *									*
*     ***********************************
      integer function lcao_dnlist_maxsize()
      implicit none

#include "lcao_dnlist.fh"

      lcao_dnlist_maxsize = maxsize
      return
      end

*     ***********************************
*     *									*
*     *	        lcao_dnlist_shift 		*	
*     *									*
*     ***********************************
      subroutine lcao_dnlist_shift()
      implicit none

#include "bafdecls.fh"
#include "lcao_dnlist.fh"

*     **** local variables ****
      integer count,indx,indxp
      integer itmp

      itmp  = int_mb(dnlist(1))
      do count=1,size-1
         indx  = dnlist(1)+count-1
         indxp = dnlist(1)+count
         int_mb(indx) = int_mb(indxp)

      end do
      indx  = dnlist(1)+size-1
      int_mb(indx) = itmp

      return
      end


