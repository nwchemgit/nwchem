      subroutine trans_t2_to_t2dagger(d_t2,k_t2_offset,d_y2,k_y2_offset)
      IMPLICIT NONE
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "errquit.fh"
#include "tce.fh"
c
      integer d_t2,k_t2_offset
      integer d_y2,k_y2_offset
c
      integer p1b, p2b, h3b, h4b
      integer p1,  p2,  h3,  h4
      integer k_r2,l_r2,k_r2_sort,l_r2_sort,size
      integer i
c
      do p1b = noab+1,noab+nvab
      do p2b = p1b,noab+nvab
      do h3b = 1,noab
      do h4b = h3b,noab
         if (int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1) .eq. 
     1       int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1)) then
         if (ieor(int_mb(k_sym+p1b-1),ieor(int_mb(k_sym+p2b-1),
     1       ieor(int_mb(k_sym+h3b-1),int_mb(k_sym+h4b-1)))) .eq. 
     1       irrep_t) then
         if ((.not.restricted).or.
     1       (int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1)
     1       +int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1).ne.8)) then
c
             size = int_mb(k_range+p1b-1) * int_mb(k_range+p2b-1)
     1               * int_mb(k_range+h3b-1) * int_mb(k_range+h4b-1)
             if (.not.ma_push_get(mt_dbl,size,'r2',l_r2,k_r2))
     1       call errquit('trans_t2_to_t2dagger: ma problem',0,ma_err)
c
             if (.not.ma_push_get(mt_dbl,size,'r2_sort',l_r2_sort,
     1           k_r2_sort))
     1       call errquit('trans_t2_to_t2dagger: ma problem',0,ma_err)
c
             call get_hash_block(d_t2,dbl_mb(k_r2),size,
     1            int_mb(k_t2_offset), h4b-1+(h3b-1)*noab
     1            +(p2b-noab-1)*noab*noab+(p1b-noab-1)*noab*noab*nvab)
c
             i=0
             do p1=1,int_mb(k_range+p1b-1)
             do p2=1,int_mb(k_range+p2b-1)
             do h3=1,int_mb(k_range+h3b-1)
             do h4=1,int_mb(k_range+h4b-1)
                dbl_mb(k_r2_sort+p2-1+(p1-1)*int_mb(k_range+p2b-1)+
     1          (h4-1)*int_mb(k_range+p2b-1)*int_mb(k_range+p1b-1)+
     1          (h3-1)*int_mb(k_range+p2b-1)*int_mb(k_range+p1b-1)*
     1          int_mb(k_range+h4b-1)) = dbl_mb(k_r2+i)
                i=i+1
             enddo
             enddo
             enddo
             enddo
c
             call put_hash_block(d_y2,dbl_mb(k_r2_sort),size,
     1               int_mb(k_y2_offset),
     1               (p2b-noab-1)+(p1b-noab-1)*(nvab)+
     1               (h4b-1)*(nvab)*(nvab)+
     1               (h3b-1)*(nvab)*(nvab)*(noab))
c
            if (.not.ma_pop_stack(l_r2_sort))
     1      call errquit('trans_t2_to_t2dagger: ma problem',1,ma_err)
c
            if (.not.ma_pop_stack(l_r2))
     1      call errquit('trans_t2_to_t2dagger: ma problem',1,ma_err)
c
         endif
         endif
         endif
c
      enddo
      enddo
      enddo
      enddo
c
      end
c $Id$
