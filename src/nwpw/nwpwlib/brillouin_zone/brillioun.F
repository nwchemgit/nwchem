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
      integer map_in

#include "stdio.fh"
#include "util.fh"
#include "bafdecls.fh"
#include "errquit.fh"
#include "brillioun_common.fh"

*     **** local variables ****
      integer taskid,MASTER
      parameter (MASTER=0)

      logical value,use_sym,oprint
      integer nb,nb0,nbq,count,taskid_k,p,nbrillioun_print
      integer kindx(2),kw(2)
      real*8  kvec(3),ksvec(4)
      real*8 total_weight 

*     **** external functions ****
      logical  control_brillioun_use_symmetry,control_print
      external control_brillioun_use_symmetry,control_print
      integer  control_num_kvectors,control_num_kvectors_print
      external control_num_kvectors,control_num_kvectors_print

      call Parallel_taskid(taskid)
      oprint = (taskid.eq.MASTER).and.control_print(print_high)

      call Parallel3d_taskid_k(taskid_k)

      use_sym = control_brillioun_use_symmetry()
      nbrillioun  = control_num_kvectors()
      nbrillioun0 = nbrillioun

*     **** set up indexing based on symmetry ****
      if (use_sym) then
*        **** allocate from MA stack ****
         value =BA_push_get(mt_int,nbrillioun0,'kindx',
     >                      kindx(2),kindx(1))
         value = value.and.BA_push_get(mt_dbl,nbrillioun0,'kw',
     >                                 kw(2),kw(1))
         if (.not. value)
     >    call errquit('brillioun_init:out of stack memory',0, MA_ERR)

*        **** reduced brillouin zone using symmetry ***
         call brillioun_symmetrize_index_weight(nbrillioun0,nbrillioun,
     >                                int_mb(kindx(1)),
     >                                dbl_mb(kw(1)))

         if (oprint) write(luout,109) nbrillioun0,nbrillioun
 109  format(/" - Symmetry reduction of Brillioun zone:",
     >       I6," to ",I6," zone points"/,
     >  " - To not use symmetry set the nwpw:brillioun_use_symmetry",
     >  " rtdb variable to false, i.e."/,
     >  "    (set nwpw:brillioun_use_symmetry .false.)")

      end if

      nbrillioun_print = control_num_kvectors_print()
      if (nbrillioun_print.lt.1) nbrillioun_print = pbrillmax
      printbrill = (nbrillioun.le.nbrillioun_print)

      map_in = 1
      call K1dB_init(nbrillioun,map_in)
      call K1dB_nq(nbrillq)

*     **** allocate from MA heap ****
      value =BA_alloc_get(mt_dbl,(nbrillq),'weight',weight(2),weight(1))
      value =value.and.
     >       BA_alloc_get(mt_dbl,(3*nbrillq),'ks',ks(2),ks(1))
      value =value.and.
     >       BA_alloc_get(mt_dbl,(3*nbrillq),'kv',kv(2),kv(1))
      value =value.and.
     >       BA_alloc_get(mt_dbl,(3*nbrillioun),'kv_all',
     >                    kv_all(2),kv_all(1))
      if (.not. value)
     > call errquit('brillioun_init:out of heap memory',0, MA_ERR)

      do nb=1,nbrillioun
         if (use_sym) then
            nb0 = int_mb(kindx(1)+nb-1)
         else
            nb0 = nb
         end if
         call control_ksvector(nb0,ksvec)
         call control_kvector(nb0,kvec)
         if (use_sym) ksvec(4) = dbl_mb(kw(1)+nb-1)

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

*     **** deallocate MA stack ****
      if (use_sym) then
         value =           BA_pop_stack(kw(2))
         value = value.and.BA_pop_stack(kindx(2))
         if (.not. value)
     >    call errquit('brillioun_init:stack memory',0, MA_ERR)
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

#include "bafdecls.fh"
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
      value = BA_alloc_get(mt_dbl,(nbrillq),
     >                     'weight',weight(2),weight(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,(3*nbrillq),'ks',ks(2),ks(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,(3*nbrillq),'kv',kv(2),kv(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,(3*nbrillioun),
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

#include "bafdecls.fh"
#include "errquit.fh"
#include "brillioun_common.fh"

*     **** local variables ****
      logical value

      call K1dB_end()

      value =           BA_free_heap(weight(2))
      value = value.and.BA_free_heap(ks(2))
      value = value.and.BA_free_heap(kv(2))
      value = value.and.BA_free_heap(kv_all(2))
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
*     *                                 *
*     *         brillioun_nbrillioun    *
*     *                                 *
*     ***********************************
      integer function brillioun_nbrillioun0()
      implicit none

#include "bafdecls.fh"
#include "brillioun_common.fh"

      brillioun_nbrillioun0 = nbrillioun0
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

#include "bafdecls.fh"
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

#include "bafdecls.fh"
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

#include "bafdecls.fh"
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

#include "bafdecls.fh"
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

#include "bafdecls.fh"
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

#include "bafdecls.fh"
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

#include "bafdecls.fh"
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

#include "bafdecls.fh"
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

#include "bafdecls.fh"
#include "brillioun_common.fh"

      brillioun_print = printbrill
      return
      end



*     *************************************
*     *                                   *
*     * brillioun_symmetrize_index_weight *
*     *                                   *
*     *************************************

      subroutine brillioun_symmetrize_index_weight(nbrillioun0,
     >                                             nbrillioun,
     >                                             kindex,kw)
      implicit none
      integer nbrillioun,nbrillioun0,kindex(*)
      real*8  kw(*)

#include "bafdecls.fh"
#include "errquit.fh"

*     **** local variables ****
      logical value
      integer op,nops,tnbrillioun,tweight(2),tks(2)
      integer nb0,nb1
      real*8  ksvec(4),ksvec1(3),xx

*     **** external functions ****
      integer  brillioun_unique_index
      external brillioun_unique_index

      do nb0=1,nbrillioun0
         kw(nb0) = 0.0d0
      end do

*     **** allocate from MA stack ****
      value = BA_push_get(mt_dbl,(4*nbrillioun0),'tks',tks(2),tks(1))
      if (.not. value)
     >call errquit('brillioun_symmetrize:out of stack',0,MA_ERR)

      nbrillioun = 0
      do nb0=1,nbrillioun0
         call control_ksvector(nb0,dbl_mb(tks(1)+4*(nbrillioun)))
         nb1 = brillioun_unique_index(nbrillioun,dbl_mb(tks(1)))
         kw(nb1) = kw(nb1) + dbl_mb(tks(1)+4*(nbrillioun)+3)
         !write(*,*) "nb0 --> nb1=",nb0,nb1
         if (nb1.eq.(nbrillioun+1)) then
            nbrillioun = nbrillioun + 1
            kindex(nbrillioun) = nb0
         end if
      end do
      
c      xx = 0.0d0
c      do nb0=1,nbrillioun
c         xx = xx + kw(nb0)
c      end do
c      write(*,*) "xx=",xx

      value = BA_pop_stack(tks(2))
      if (.not. value)
     >call errquit('brillioun_symmetrize:pop stack',0, MA_ERR)

      return
      end

*     *****************************************
*     *                                       *
*     *        brillion_unique_index          *
*     *                                       *
*     *****************************************
      integer function brillioun_unique_index(nbrillm,ksvec_all)
      implicit none
      integer nbrillm
      real*8 ksvec_all(4,*)
      
*     ****local variables ****
      logical found
      integer nb1
      real*8  kc(4)

*     **** external functions ****
      logical  brillioun_k_sym_compare
      external brillioun_k_sym_compare

      kc(1) = ksvec_all(1,nbrillm+1)
      kc(2) = ksvec_all(2,nbrillm+1)
      kc(3) = ksvec_all(3,nbrillm+1)

      nb1 = 1
      found = .false.
      do while ((nb1.le.nbrillm).and.(.not.found))
         found = brillioun_k_sym_compare(kc,ksvec_all(1,nb1))
         if (.not.found) nb1 = nb1 + 1
      end do

      brillioun_unique_index = nb1
      return
      end

*     *****************************************
*     *                                       *
*     *        brillion_k_sym_compare         *
*     *                                       *
*     *****************************************
      logical function brillioun_k_sym_compare(k1,k2)
      implicit none
      real*8 k1(4),k2(4)

      logical same
      integer op,nops,i
      real*8 sksvec(3),ksvec1(3),ksvec2(3),xx,yy

*     **** external functions ****
      integer  ion_sym_number_ops
      external ion_sym_number_ops

      do i=1,3
         ksvec1(i) = k1(i) + 0.5d0
         ksvec2(i) = k2(i) + 0.5d0
         do while (ksvec1(i).ge.1.0d0)
            ksvec1(i) = ksvec1(i)-1.0d0
         end do 
         do while (ksvec1(i).lt.0.0d0)
            ksvec1(i) = ksvec1(i)+1.0d0
         end do

         do while (ksvec2(i).ge.1.0d0)
            ksvec2(i) = ksvec2(i)-1.0d0
         end do 
         do while (ksvec2(i).lt.0.0d0)
            ksvec2(i) = ksvec2(i)+1.0d0
         end do
      end do
      nops = ion_sym_number_ops()
      op = 1
      same = .false.
      do while ((op.le.nops).and.(.not.same))
         call ion_sym_apply_op(op,ksvec1,sksvec)
         op = op + 1
c         do i=1,3
c            do while (sksvec(i).ge.1.0d0)
c               sksvec(i) = sksvec(i)-1.0d0
c            end do 
c            do while (sksvec(i).lt.0.0d0)
c               sksvec(i) = sksvec(i)+1.0d0
c            end do
c         end do
         xx =      (sksvec(1) - ksvec2(1))**2
         xx = xx + (sksvec(2) - ksvec2(2))**2
         xx = xx + (sksvec(3) - ksvec2(3))**2
         xx = dsqrt(xx)
         yy =      (sksvec(1) + ksvec2(1))**2
         yy = yy + (sksvec(2) + ksvec2(2))**2
         yy = yy + (sksvec(3) + ksvec2(3))**2
         yy = dsqrt(yy)
         if ((xx.lt.1.0d-6).or.(yy.lt.1.0d-6)) same = .true.
      end do

      brillioun_k_sym_compare = same
      return
      end
