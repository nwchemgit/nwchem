       subroutine create_mr_evl_sorted(k_f1_offsetm)
        implicit none
#include "tce.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "errquit.fh"
#include "sym.fh"
#include "tce_mrcc.fh"
#include "tce_main.fh"
#include "global.fh"



      integer k_f1_offsetm(maxref)
c      integer nref,iref
      integer iref
c      logical nodezero 
      integer g6b,g5b
      integer i,j,k,l
      integer l_f1,k_f1
      integer size
      integer offsetf

c      nodezero = (ga_nodeid().eq.0)

      do iref=1,nref


         k = 0
         k_sym = k_symm(iref)
         k_offset = k_offsetm(iref)
         k_range = k_rangem(iref)
         k_spin = k_spinm(iref)
         k_movecs_sorted = k_movecs_sortedm(iref)

         noa = nblcks(1,iref)
         nob = nblcks(2,iref)
         nva = nblcks(3,iref)
         nvb = nblcks(4,iref)

         noab = noa+nob
         nvab = nva+nvb

      DO g6b = 1,noab+nvab
      g5b = g6b
      IF (int_mb(k_spin+g6b-1) .eq. int_mb(k_spin+g5b-1)) THEN
      IF (ieor(int_mb(k_sym+g6b-1),int_mb(k_sym+g5b-1)) .eq. irrep_f) TH
     &EN
      IF ((.not.restricted).or.(int_mb(k_spin+g6b-1)+int_mb(k_spin+g5b-1
     &).ne.4)) THEN

        size = int_mb(k_range+g6b-1) * int_mb(k_range+g5b-1)

        if (.not.ma_push_get(mt_dbl,size,'f1mi',l_f1,k_f1))
     1   call errquit('tce_mrcc_iface_f1: MA problem get',0,MA_ERR)


        call get_hash_block(d_f1m(iref),dbl_mb(k_f1),size,
     1   int_mb(k_f1_offsetm(iref)),g6b-1+(noab+nvab)*(g5b-1))
        l = 0
        do i=1,int_mb(k_range+g6b-1)
        do j=1,int_mb(k_range+g5b-1)
        l = l + 1
           if(i.eq.j) then
             k = k + 1
             dbl_mb(k_evl_sortedm(iref)+k-1) = dbl_mb(k_f1-1+l)
           endif
        enddo
        enddo

        if (.not.ma_pop_stack(l_f1))
     1   call errquit('tce_mrcc_iface_f1: MA problem pop',1,MA_ERR)

      
      END IF
      END IF
      END IF
      END DO

      enddo

      return
      end

c $Id$
