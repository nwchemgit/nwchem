      subroutine trans_t1_to_t1dagger(d_t1,k_t1_offset,d_y1,k_y1_offset)
      IMPLICIT NONE
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "errquit.fh"
#include "tce.fh"
c
      integer d_t1,k_t1_offset
      integer d_y1,k_y1_offset
c
      integer p1b, h2b
      integer p1,  h2
      integer k_r2,l_r2,k_r2_sort,l_r2_sort,size
      integer i
c
      do p1b = noab+1,noab+nvab
      do h2b = 1,noab
         if (int_mb(k_spin+p1b-1) .eq. int_mb(k_spin+h2b-1)) then
         if (ieor(int_mb(k_sym+p1b-1),int_mb(k_sym+h2b-1)) 
     1       .eq. irrep_t) then
         if ((.not.restricted).or.
     1       (int_mb(k_spin+p1b-1)+int_mb(k_spin+h2b-1).ne.4)) then
c
             size = int_mb(k_range+p1b-1) * int_mb(k_range+h2b-1)
             if (.not.ma_push_get(mt_dbl,size,'r2',l_r2,k_r2))
     1       call errquit('trans_t1_to_t1dagger: ma problem',0,MA_ERR)
c
             if (.not.ma_push_get(mt_dbl,size,'r2_sort',l_r2_sort,
     1           k_r2_sort))
     1       call errquit('trans_t1_to_t1dagger: ma problem',0,MA_ERR)
c
             call get_hash_block(d_t1,dbl_mb(k_r2),size,
     1            int_mb(k_t1_offset), h2b-1+(p1b-noab-1)*noab)
c
             i=0
             do p1=1,int_mb(k_range+p1b-1)
             do h2=1,int_mb(k_range+h2b-1)
                dbl_mb(k_r2_sort+p1-1+(h2-1)*int_mb(k_range+p1b-1)) 
     1          = dbl_mb(k_r2+i)
                i=i+1
             enddo
             enddo
c
             call put_hash_block(d_y1,dbl_mb(k_r2_sort),size,
     1               int_mb(k_y1_offset),
     1               (p1b-noab-1)+(h2b-1)*(nvab))
c
            if (.not.ma_pop_stack(l_r2_sort))
     1      call errquit('trans_t1_to_t1dagger: ma problem',1,MA_ERR)
c
            if (.not.ma_pop_stack(l_r2))
     1      call errquit('trans_t1_to_t1dagger: ma problem',1,MA_ERR)
c
         endif
         endif
         endif
c
      enddo
      enddo
c
      end
c $Id$
