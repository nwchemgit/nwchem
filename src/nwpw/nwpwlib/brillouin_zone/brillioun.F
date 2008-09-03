*
* $Id: brillioun.F,v 1.5 2008-09-03 17:07:35 bylaska Exp $
*

*     ***********************************
*     *					*
*     *	       brillioun_init		*
*     *					*
*     ***********************************
      subroutine brillioun_init()
      implicit none
      integer map_in

#include "mafdecls.fh"
#include "errquit.fh"
#include "brillioun_common.fh"

*     **** local variables ****
      logical value
      integer nb,nbq,count,taskid_k,p
      real*8  kvec(3),ksvec(4)
      real*8 total_weight 

*     **** external functions ****
      integer  control_num_kvectors
      external control_num_kvectors

      call Parallel3d_taskid_k(taskid_k)

      nbrillioun = control_num_kvectors()
      printbrill = (nbrillioun.le.pbrillmax)

      map_in = 1
      call K1dB_init(nbrillioun,map_in)
      call K1dB_nq(nbrillq)

*     **** allocate from MA heap ****
      value =MA_alloc_get(mt_dbl,(nbrillq),'weight',weight(2),weight(1))
      value =value.and.
     >       MA_alloc_get(mt_dbl,(3*nbrillq),'ks',ks(2),ks(1))
      value =value.and.
     >       MA_alloc_get(mt_dbl,(3*nbrillq),'kv',kv(2),kv(1))
      value =value.and.
     >       MA_alloc_get(mt_dbl,(3*nbrillioun),'kv_all',
     >                    kv_all(2),kv_all(1))
      if (.not. value)
     > call errquit('brillioun_init:out of heap memory',0, MA_ERR)

      do nb=1,nbrillioun
         call control_ksvector(nb,ksvec)
         call control_kvector(nb,kvec)

         dbl_mb(kv_all(1)+3*(nb-1))    = kvec(1)
         dbl_mb(kv_all(1)+3*(nb-1)+1)  = kvec(2)
         dbl_mb(kv_all(1)+3*(nb-1)+2)  = kvec(3)

         call K1dB_ktoqp(nb,nbq,p)

         if (p.eq.taskid_k) then
            dbl_mb(kv(1)+3*(nbq-1))    = kvec(1)
            dbl_mb(kv(1)+3*(nbq-1)+1)  = kvec(2)
            dbl_mb(kv(1)+3*(nbq-1)+2)  = kvec(3)

            dbl_mb(ks(1)+3*(nbq-1))    = ksvec(1)
            dbl_mb(ks(1)+3*(nbq-1)+1)  = ksvec(2)
            dbl_mb(ks(1)+3*(nbq-1)+2)  = ksvec(3)
            dbl_mb(weight(1)+nbq-1)    = ksvec(4)
         end if
      end do

*     *** set the default weights *****
      count = 0
      total_weight = 0.0d0
      do nbq=1,nbrillq
         if (dbl_mb(weight(1)+nbq-1).gt.0.0d0) then
           total_weight = total_weight - dbl_mb(weight(1)+nbq-1)
         else
           count = count + 1
         end if
      end do
      call K1dB_ISumAll(count)
      call K1dB_SumAll(total_weight)
      total_weight = 1.0d0 - total_weight
     
      if (count.gt.0) then
         total_weight = total_weight/dble(count)
         do nbq=1,nbrillq
            if (dbl_mb(weight(1)+nbq-1).le.0.0d0) then
               dbl_mb(weight(1)+nbq-1) = total_weight
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
      subroutine brillioun_structure_init(nbstart,nbk)
      implicit none
      integer nbstart,nbk

#include "mafdecls.fh"
#include "errquit.fh"
#include "brillioun_common.fh"

*     **** local variables ****
      logical value
      real*8  kvec(3),ksvec(4)
      integer map_in,nb,nbq,taskid_k,p

      call Parallel3d_taskid_k(taskid_k)
      nbrillioun = nbk
      
      map_in = 1
      call K1dB_init(nbrillioun,map_in)
      call K1dB_nq(nbrillq)

*     **** allocate from MA heap ****
      value = MA_alloc_get(mt_dbl,(nbrillq),
     >                     'weight',weight(2),weight(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,(3*nbrillq),'ks',ks(2),ks(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,(3*nbrillq),'kv',kv(2),kv(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,(3*nbrillioun),
     >                     'kv_all',kv_all(2),kv_all(1))
      if (.not. value)
     > call errquit('brillioun_structure_init:out of heap memory',
     >              0,MA_ERR)

      do nb=1,nbrillioun
        call control_ksvector_structure((nb-1)+nbstart,ksvec)
        call control_kvector_structure((nb-1) +nbstart,kvec)

        dbl_mb(kv_all(1)+3*(nb-1))    = kvec(1)
        dbl_mb(kv_all(1)+3*(nb-1)+1)  = kvec(2)
        dbl_mb(kv_all(1)+3*(nb-1)+2)  = kvec(3)

        call K1dB_ktoqp(nb,nbq,p)
        if (p.eq.taskid_k) then
          dbl_mb(kv(1)+3*(nbq-1))    = kvec(1)
          dbl_mb(kv(1)+3*(nbq-1)+1)  = kvec(2)
          dbl_mb(kv(1)+3*(nbq-1)+2)  = kvec(3)

          dbl_mb(ks(1)+3*(nbq-1))    = ksvec(1)
          dbl_mb(ks(1)+3*(nbq-1)+1)  = ksvec(2)
          dbl_mb(ks(1)+3*(nbq-1)+2)  = ksvec(3)
          dbl_mb(weight(1)+(nbq-1))  = 1.0d0
        end if
      end do

      return
      end


*     ***********************************
*     *					*
*     *		brillioun_end		*
*     *					*
*     ***********************************
      subroutine brillioun_end()
      implicit none

#include "mafdecls.fh"
#include "errquit.fh"
#include "brillioun_common.fh"

*     **** local variables ****
      logical value

      call K1dB_end()

      value =           MA_free_heap(weight(2))
      value = value.and.MA_free_heap(ks(2))
      value = value.and.MA_free_heap(kv(2))
      value = value.and.MA_free_heap(kv_all(2))
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

#include "mafdecls.fh"
#include "brillioun_common.fh"

      brillioun_nbrillioun = nbrillioun
      return
      end

*     ***********************************
*     *					*
*     *		brillioun_nbrillq   	*
*     *					*
*     ***********************************
      integer function brillioun_nbrillq()
      implicit none
!      integer nb

#include "mafdecls.fh"
#include "brillioun_common.fh"

      brillioun_nbrillq = nbrillq
      return
      end

*     ***********************************
*     *					*
*     *		brillioun_weight	*
*     *					*
*     ***********************************
      real*8 function brillioun_weight(nbq)
      implicit none
      integer nbq

#include "mafdecls.fh"
#include "brillioun_common.fh"

      brillioun_weight = dbl_mb(weight(1)+nbq-1)
      return
      end

*     ***********************************
*     *					*
*     *		 brillioun_ks		*
*     *					*
*     ***********************************
      real*8 function brillioun_ks(i,nbq)
      implicit none
      integer i,nbq

#include "mafdecls.fh"
#include "brillioun_common.fh"

      brillioun_ks = dbl_mb(ks(1)+3*(nbq-1)+i-1)
      return
      end

*     ***********************************
*     *					*
*     *		  brillioun_k	        *
*     *					*
*     ***********************************
      real*8 function brillioun_k(i,nbq)
      implicit none
      integer i,nbq

#include "mafdecls.fh"
#include "brillioun_common.fh"

      brillioun_k = dbl_mb(kv(1)+3*(nbq-1)+i-1)
      return
      end



*     ***********************************
*     *                                 *
*     *           brillioun_all_k       *
*     *                                 *
*     ***********************************
      real*8 function brillioun_all_k(i,nb)
      implicit none
      integer i,nb

#include "mafdecls.fh"
#include "brillioun_common.fh"

      brillioun_all_k = dbl_mb(kv_all(1)+3*(nb-1)+i-1)
      return
      end



*     ***********************************
*     *					*
*     *		brillioun_weight_brdcst	*
*     *					*
*     ***********************************
      real*8 function brillioun_weight_brdcst(nb)
      implicit none
      integer nb

#include "mafdecls.fh"
#include "brillioun_common.fh"

      integer nbq,taskid_k,p
      real*8 w 

      call Parallel3d_taskid_k(taskid_k)
      call K1dB_ktoqp(nb,nbq,p)

      w = 0.0d0
      if (p.eq.taskid_k) w = dbl_mb(weight(1)+nbq-1)
      call K1dB_SumAll(w)

      brillioun_weight_brdcst = w
      return
      end

*     ***********************************
*     *					*
*     *		 brillioun_ks_brdcst	*
*     *					*
*     ***********************************
      real*8 function brillioun_ks_brdcst(i,nb)
      implicit none
      integer i,nb

#include "mafdecls.fh"
#include "brillioun_common.fh"

      integer nbq,taskid_k,p
      real*8 w 

      call Parallel3d_taskid_k(taskid_k)
      call K1dB_ktoqp(nb,nbq,p)

      w = 0.0d0
      if (p.eq.taskid_k) w = dbl_mb(ks(1)+3*(nbq-1)+i-1)
      call K1dB_SumAll(w)

      brillioun_ks_brdcst = w
      return
      end

*     ***********************************
*     *					*
*     *		  brillioun_k_brdcst    *
*     *					*
*     ***********************************
      real*8 function brillioun_k_brdcst(i,nb)
      implicit none
      integer i,nb

#include "mafdecls.fh"
#include "brillioun_common.fh"

      integer nbq,taskid_k,p
      real*8 w 

      call Parallel3d_taskid_k(taskid_k)
      call K1dB_ktoqp(nb,nbq,p)

      w = 0.0d0
      if (p.eq.taskid_k) w = dbl_mb(kv(1)+3*(nbq-1)+i-1)
      call K1dB_SumAll(w)

      brillioun_k_brdcst = w
      return
      end

*     ***********************************
*     *                                 *
*     *           brillioun_print      *
*     *                                 *
*     ***********************************
      logical function brillioun_print()
      implicit none
      integer i,nb

#include "mafdecls.fh"
#include "brillioun_common.fh"

      brillioun_print = printbrill
      return
      end


