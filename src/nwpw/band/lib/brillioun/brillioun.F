*
* $Id$
*

*     ***********************************
*     *					*
*     *	       brillioun_init		*
*     *					*
*     ***********************************
      subroutine brillioun_init()
      implicit none
#include "errquit.fh"

#include "bafdecls.fh"
#include "brillioun_common.fh"

*     **** local variables ****
      logical value
      integer nb,count
      real*8  kvec(3),ksvec(4)
      real*8 total_weight 

*     **** external functions ****
      integer  control_num_kvectors
      external control_num_kvectors

      nbrillioun = control_num_kvectors()

*     **** allocate from MA heap ****
      value = BA_alloc_get(mt_dbl,(nbrillioun),
     >                     'weight',weight(2),weight(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,(3*nbrillioun),'ks',ks(2),ks(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,(3*nbrillioun),'kv',kv(2),kv(1))
      if (.not. value)
     > call errquit('brillioun_init:out of heap memory',0, MA_ERR)


      do nb=1,nbrillioun
         call control_ksvector(nb,ksvec)
         call control_kvector(nb,kvec)
         dbl_mb(kv(1)+3*(nb-1))    = kvec(1)
         dbl_mb(kv(1)+3*(nb-1)+1)  = kvec(2)
         dbl_mb(kv(1)+3*(nb-1)+2)  = kvec(3)

         dbl_mb(ks(1)+3*(nb-1))    = ksvec(1)
         dbl_mb(ks(1)+3*(nb-1)+1)  = ksvec(2)
         dbl_mb(ks(1)+3*(nb-1)+2)  = ksvec(3)
         dbl_mb(weight(1)+nb-1)    = ksvec(4)
      end do


*     *** set the default weights *****
      total_weight = 1.0d0
      count = 0
      do nb=1,nbrillioun
         if (dbl_mb(weight(1)+nb-1).gt.0.0d0) then
           total_weight = total_weight - dbl_mb(weight(1)+nb-1)
         else
           count = count + 1
         end if
      end do
      if (count.gt.0) then
         total_weight = total_weight/dble(count)
         do nb=1,nbrillioun
            if (dbl_mb(weight(1)+nb-1).le.0.0d0) then
               dbl_mb(weight(1)+nb-1) = total_weight
            end if
         end do
      end if

      return
      end


*     ***********************************
*     *					*
*     *	   brillioun_structure_init	*
*     *					*
*     ***********************************
      subroutine brillioun_structure_init(nb)
      implicit none
      integer nb

#include "bafdecls.fh"
#include "errquit.fh"
#include "brillioun_common.fh"

*     **** local variables ****
      logical value
      real*8  kvec(3),ksvec(4)

      nbrillioun = 1

*     **** allocate from MA heap ****
      value = BA_alloc_get(mt_dbl,(nbrillioun),
     >                     'weight',weight(2),weight(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,(3*nbrillioun),'ks',ks(2),ks(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,(3*nbrillioun),'kv',kv(2),kv(1))
      if (.not. value)
     > call errquit('brillioun_structure_init:out of heap memory',
     >              0,MA_ERR)


      call control_ksvector_structure(nb,ksvec)
      call control_kvector_structure(nb,kvec)
      dbl_mb(kv(1))    = kvec(1)
      dbl_mb(kv(1)+1)  = kvec(2)
      dbl_mb(kv(1)+2)  = kvec(3)

      dbl_mb(ks(1))    = ksvec(1)
      dbl_mb(ks(1)+1)  = ksvec(2)
      dbl_mb(ks(1)+2)  = ksvec(3)
      dbl_mb(weight(1))= 1.0d0

      return
      end



*     ***********************************
*     *					*
*     *		brillioun_end		*
*     *					*
*     ***********************************
      subroutine brillioun_end()
      implicit none
#include "errquit.fh"

#include "bafdecls.fh"
#include "brillioun_common.fh"

*     **** local variables ****
      logical value

      value = BA_free_heap(weight(2))
      value = value.and.BA_free_heap(ks(2))
      value = value.and.BA_free_heap(kv(2))
      if (.not. value)
     > call errquit('brillioun_end:freeing heap memory',0, MA_ERR)

      return
      end 

*     ***********************************
*     *					*
*     *		brillioun_nbrillioun 	*
*     *					*
*     ***********************************
      integer function brillioun_nbrillioun()
      implicit none
!      integer nb

#include "bafdecls.fh"
#include "brillioun_common.fh"

      brillioun_nbrillioun = nbrillioun
      return
      end

*     ***********************************
*     *					*
*     *		brillioun_weight	*
*     *					*
*     ***********************************
      real*8 function brillioun_weight(nb)
      implicit none
      integer nb

#include "bafdecls.fh"
#include "brillioun_common.fh"

      brillioun_weight = dbl_mb(weight(1)+nb-1)
      return
      end

*     ***********************************
*     *					*
*     *		 brillioun_ks		*
*     *					*
*     ***********************************
      real*8 function brillioun_ks(i,nb)
      implicit none
      integer i,nb

#include "bafdecls.fh"
#include "brillioun_common.fh"

      brillioun_ks = dbl_mb(ks(1)+3*(nb-1)+i-1)
      return
      end

*     ***********************************
*     *					*
*     *		  brillioun_k	        *
*     *					*
*     ***********************************
      real*8 function brillioun_k(i,nb)
      implicit none
      integer i,nb

#include "bafdecls.fh"
#include "brillioun_common.fh"

      brillioun_k = dbl_mb(kv(1)+3*(nb-1)+i-1)
      return
      end

*     *************************************
*     *                                   *
*     *       brillioun_symmetrize        *
*     *                                   *
*     *************************************

      subroutine brillioun_symmetrize()
      implicit none

#include "bafdecls.fh"
#include "brillioun_common.fh"
#include "errquit.fh"

*     **** local variables ****
      logical value
      integer nops,tnbrillioun,tweight(2),tks(2)

*     **** external functions ****
      integer  ion_sym_number_ops
      external ion_sym_number_ops

      nops = ion_sym_number_ops()

      write(*,*) "nops=",nops

c*     **** allocate from MA stack ****
c      value = BA_alloc_get(mt_dbl,(nbrillioun),
c     >                     'tweight',tweight(2),tweight(1))
c      value = value.and.
c     >        BA_alloc_get(mt_dbl,(3*nbrillioun),'tks',tks(2),tks(1))
c      if (.not. value)
c     > call errquit('brillioun_symmetrize:out of stack memory',0,MA_ERR)
c
c
c      value = BA_pop_stack(ks(2))
c      value = value.and.BA_pop_stack(weight(2))
c      if (.not. value)
c     > call errquit('brillioun_symmetrize:pop stack',0, MA_ERR)

      return
      end 

