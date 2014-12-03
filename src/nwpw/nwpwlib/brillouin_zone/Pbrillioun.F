*
* $Id$
*

*     ***********************************
*     *					*
*     *	       Pbrillioun_init		*
*     *					*
*     ***********************************
      subroutine Pbrillioun_init()
      implicit none

#include "bafdecls.fh"
#include "errquit.fh"
#include "brillioun_common.fh"

*     **** local variables ****
      logical value
      integer nb,count
      real*8  kvec(3),ksvec(4)
      real*8 total_weight 

*     **** external functions ****
      printbrill = (nbrillioun.le.pbrillmax)

      if (printbrill) then

*     **** allocate from MA heap ****
      value =BA_alloc_get(mt_dbl,(nbrillioun),'pweight',
     >                    pweight(2),pweight(1))
      value =value.and.
     >       BA_alloc_get(mt_dbl,(3*nbrillioun),'pks',pks(2),pks(1))
      value =value.and.
     >       BA_alloc_get(mt_dbl,(3*nbrillioun),'pkv',pkv(2),pkv(1))
      if (.not. value)
     > call errquit('Pbrillioun_init:out of heap memory',0, MA_ERR)


      do nb=1,nbrillioun
         call control_ksvector(nb,ksvec)
         call control_kvector(nb,kvec)
         dbl_mb(pkv(1)+3*(nb-1))    = kvec(1)
         dbl_mb(pkv(1)+3*(nb-1)+1)  = kvec(2)
         dbl_mb(pkv(1)+3*(nb-1)+2)  = kvec(3)

         dbl_mb(pks(1)+3*(nb-1))    = ksvec(1)
         dbl_mb(pks(1)+3*(nb-1)+1)  = ksvec(2)
         dbl_mb(pks(1)+3*(nb-1)+2)  = ksvec(3)
         dbl_mb(pweight(1)+nb-1)    = ksvec(4)
      end do


*     *** set the default weights *****
      total_weight = 1.0d0
      count = 0
      do nb=1,nbrillioun
         if (dbl_mb(pweight(1)+nb-1).gt.0.0d0) then
           total_weight = total_weight - dbl_mb(pweight(1)+nb-1)
         else
           count = count + 1
         end if
      end do
     
      if (count.gt.0) then
         total_weight = total_weight/dble(count)
         do nb=1,nbrillioun
            if (dbl_mb(pweight(1)+nb-1).le.0.0d0) then
                dbl_mb(pweight(1)+nb-1) = total_weight
            end if
         end do
      end if

      end if

      return
      end


*     ***********************************
*     *					*
*     *		Pbrillioun_end		*
*     *					*
*     ***********************************
      subroutine Pbrillioun_end()
      implicit none

#include "bafdecls.fh"
#include "errquit.fh"
#include "brillioun_common.fh"

*     **** local variables ****
      logical value


      if (printbrill) then
      value =           BA_free_heap(pweight(2))
      value = value.and.BA_free_heap(pks(2))
      value = value.and.BA_free_heap(pkv(2))
      if (.not. value)
     > call errquit('brillioun_end:freeing heap memory',0, MA_ERR)
      end if

      return
      end 

*     ***********************************
*     *					*
*     *		Pbrillioun_weight	*
*     *					*
*     ***********************************
      real*8 function Pbrillioun_weight(nb)
      implicit none
      integer nb

#include "bafdecls.fh"
#include "brillioun_common.fh"

      if (printbrill) then
      Pbrillioun_weight = dbl_mb(pweight(1)+nb-1)
      else
      Pbrillioun_weight = 0.0d0
      end if
      return
      end

*     ***********************************
*     *					*
*     *		 Pbrillioun_ks		*
*     *					*
*     ***********************************
      real*8 function Pbrillioun_ks(i,nb)
      implicit none
      integer i,nb

#include "bafdecls.fh"
#include "brillioun_common.fh"

      if (printbrill) then
      Pbrillioun_ks = dbl_mb(pks(1)+3*(nb-1)+i-1)
      else
      Pbrillioun_ks = 0.0d0
      end if
      return
      end

*     ***********************************
*     *					*
*     *		  Pbrillioun_k	        *
*     *					*
*     ***********************************
      real*8 function Pbrillioun_k(i,nb)
      implicit none
      integer i,nb

#include "bafdecls.fh"
#include "brillioun_common.fh"

      if (printbrill) then
      Pbrillioun_k = dbl_mb(pkv(1)+3*(nb-1)+i-1)
      else
      Pbrillioun_k = 0.0d0
      end if

      return
      end


*     ***********************************
*     *                                 *
*     *           Pbrillioun_print      *
*     *                                 *
*     ***********************************
      logical function Pbrillioun_print()
      implicit none
      integer i,nb

#include "bafdecls.fh"
#include "brillioun_common.fh"

      Pbrillioun_print = printbrill
      return
      end


