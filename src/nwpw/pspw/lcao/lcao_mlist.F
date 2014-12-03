*
* $Id$
*

*     ***********************************************************
*     *															*
*     *   		   		lcao_mlist library						*
*     *															*
*     *   Author - Eric Bylaska									*
*     *   date   - 10/09/01										*
*     *															*
*     ***********************************************************

*
*

*     ***********************************
*     *									*
*     *	       lcao_mlist_init			*	
*     *									*
*     ***********************************
      subroutine lcao_mlist_init(maxsize0,ispin0,nbasis0)
      implicit none
#include "errquit.fh"
      integer maxsize0
      integer ispin0
      integer nbasis0

#include "bafdecls.fh"
#include "lcao_mlist.fh"

*     **** local variables ****
      logical value

      maxsize = maxsize0
      ispin   = ispin0
      nbasis  = nbasis0
      size    = 0

*    ***** allocate memory from heap memory ****
      value = BA_alloc_get(mt_int,(ispin*maxsize),'mlist',
     >                     mlist(2),mlist(1))
      if (.not. value)
     > call errquit('lcao_mlist_init: out of heap memory',0, MA_ERR)

      return
      end


*     ***********************************
*     *									*
*     *	       lcao_mlist_end			*	
*     *									*
*     ***********************************

      subroutine lcao_mlist_end()
      implicit none
#include "errquit.fh"
     
#include "bafdecls.fh"
#include "global.fh"
#include "lcao_mlist.fh"


*     **** local variables ****
      logical value
      integer i,mindx

      value = .true.
      do i=1,size
         mindx = int_mb(mlist(1)+i-1)
         value = value.and.ga_destroy(mindx)
      end do
      if (ispin.gt.1) then
        do i=1,size
         mindx = int_mb(mlist(1)+i-1+(ispin-1)*maxsize)
         value = value.and.ga_destroy(mindx)
        end do
      end if
      value = value.and.BA_free_heap(mlist(2))
      if (.not. value)
     > call errquit('lcao_mlist_init: freeing heap memory',0, MA_ERR)

      size    = 0
      maxsize = 0
      return
      end

*     ***********************************
*     *									*
*     *	        lcao_mlist_new  		*	
*     *									*
*     ***********************************
      subroutine lcao_mlist_new(index)
      implicit none
#include "errquit.fh"
      integer index(*)

#include "bafdecls.fh"
#include "global.fh"
#include "lcao_mlist.fh"

*     **** local variables ****
      logical value
      integer indx,indx2

      if (size.ge.maxsize) 
     > call errquit('lcao_mlist_new: maxsize exceeded',0, MA_ERR)

      value = ga_create(mt_dbl,nbasis,nbasis,'mlist',1,1,indx)
      if (ispin.gt.1)
     > value = value.and.
     >         ga_create(mt_dbl,nbasis,nbasis,'mlist',1,1,indx2)
      if (.not.value) 
     > call errquit('lcao_mlist_new: out of ga memory',0, GA_ERR)

      size = size + 1
      int_mb(mlist(1)+size-1) = indx
      if (ispin.gt.1)
     > int_mb(mlist(1)+size-1+(ispin-1)*maxsize) = indx2

      index(1)  = indx
      if (ispin.gt.1)
     > index(2)  = indx2
      return
      end
         


*     ***********************************
*     *									*
*     *	        lcao_mlist_get  		*	
*     *									*
*     ***********************************
      subroutine lcao_mlist_get(count,index)
      implicit none
#include "errquit.fh"
      integer count
      integer index(*)

#include "bafdecls.fh"
#include "lcao_mlist.fh"

      if ((count.gt.size).or.(count.lt.1)) 
     > call errquit('lcao_mlist_get: bad count',0, INPUT_ERR)

      index(1)   = int_mb(mlist(1)+count-1)
      if (ispin.gt.1)
     > index(2)   = int_mb(mlist(1)+count-1+(ispin-1)*maxsize)
      return
      end


*     ***********************************
*     *									*
*     *	        lcao_mlist_size 		*	
*     *									*
*     ***********************************
      integer function lcao_mlist_size()
      implicit none

#include "lcao_mlist.fh"

      lcao_mlist_size = size
      return
      end

*     ***********************************
*     *									*
*     *	        lcao_mlist_maxsize 	*	
*     *									*
*     ***********************************
      integer function lcao_mlist_maxsize()
      implicit none

#include "lcao_mlist.fh"

      lcao_mlist_maxsize = maxsize
      return
      end

*     ***********************************
*     *									*
*     *	        lcao_mlist_shift 		*	
*     *									*
*     ***********************************
      subroutine lcao_mlist_shift()
      implicit none

#include "bafdecls.fh"
#include "lcao_mlist.fh"

*     **** local variables ****
      integer count,indx,indxp
      integer itmp

      itmp  = int_mb(mlist(1))
      do count=1,size-1
         indx  = mlist(1)+count-1
         indxp = mlist(1)+count
         int_mb(indx) = int_mb(indxp)

      end do
      indx  = mlist(1)+size-1
      int_mb(indx) = itmp

      if (ispin.gt.1) then
        itmp  = int_mb(mlist(1)+(ispin-1)*maxsize)
        do count=1,size-1
           indx  = mlist(1)+count-1+(ispin-1)*maxsize
           indxp = mlist(1)+count  +(ispin-1)*maxsize
           int_mb(indx) = int_mb(indxp)

        end do
        indx  = mlist(1)+size-1+(ispin-1)*maxsize
        int_mb(indx) = itmp
      end if

      return
      end


