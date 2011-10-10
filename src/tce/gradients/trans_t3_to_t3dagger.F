      subroutine trans_t3_to_t3dagger(d_t3,k_t3_offset,d_y3,k_y3_offset)
      IMPLICIT NONE
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "errquit.fh"
#include "tce.fh"
c
      integer d_t3,k_t3_offset
      integer d_y3,k_y3_offset
c
      integer p1b, p2b, p3b, h4b, h5b, h6b
      integer p1,  p2,  p3,  h4,  h5,  h6
      integer k_r2,l_r2,k_r2_sort,l_r2_sort,size
      integer i
c
      do p1b = noab+1,noab+nvab
      do p2b = p1b,noab+nvab
      do p3b = p2b,noab+nvab
      do h4b = 1,noab
      do h5b = h4b,noab
      do h6b = h5b,noab
         if (int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1)
     1       +int_mb(k_spin+p3b-1) .eq. 
     1       int_mb(k_spin+h4b-1)+int_mb(k_spin+h5b-1)
     1       +int_mb(k_spin+h6b-1)) then
         if (ieor(int_mb(k_sym+p1b-1),ieor(int_mb(k_sym+p2b-1),
     1       ieor(int_mb(k_sym+p3b-1),ieor(int_mb(k_sym+h4b-1),
     1       ieor(int_mb(k_sym+h5b-1),int_mb(k_sym+h6b-1)))))) .eq. 
     1       irrep_t) then
         if ((.not.restricted).or.(int_mb(k_spin+p1b-1)
     1       +int_mb(k_spin+p2b-1)+int_mb(k_spin+p3b-1)
     1       +int_mb(k_spin+h4b-1)+int_mb(k_spin+h5b-1)
     1       +int_mb(k_spin+h6b-1).ne.12)) then
c
             size = int_mb(k_range+p1b-1) * int_mb(k_range+p2b-1) 
     1            * int_mb(k_range+p3b-1) * int_mb(k_range+h4b-1) 
     1            * int_mb(k_range+h5b-1) * int_mb(k_range+h6b-1)
             if (.not.ma_push_get(mt_dbl,size,'r2',l_r2,k_r2))
     1       call errquit('trans_t3_to_t3dagger: ma problem',0,ma_err)
c
             if (.not.ma_push_get(mt_dbl,size,'r2_sort',l_r2_sort,
     1           k_r2_sort))
     1       call errquit('trans_t3_to_t3dagger: ma problem',0,ma_err)
c
             call get_hash_block(d_t3,dbl_mb(k_r2),size,
     1            int_mb(k_t3_offset), 
     1            h6b-1
     1          +(h5b-1)*noab
     1          +(h4b-1)*noab*noab
     1          +(p3b-noab-1)*noab*noab*noab
     1          +(p2b-noab-1)*noab*noab*noab*nvab
     1          +(p1b-noab-1)*noab*noab*noab*nvab*nvab)
c
             i=0
             do p1=1,int_mb(k_range+p1b-1)
             do p2=1,int_mb(k_range+p2b-1)
             do p3=1,int_mb(k_range+p3b-1)
             do h4=1,int_mb(k_range+h4b-1)
             do h5=1,int_mb(k_range+h5b-1)
             do h6=1,int_mb(k_range+h6b-1)
                dbl_mb(k_r2_sort
     1         +(p3-1)
     1         +(p2-1)*int_mb(k_range+p3b-1)
     1         +(p1-1)*int_mb(k_range+p3b-1)*int_mb(k_range+p2b-1)
     1         +(h6-1)*int_mb(k_range+p3b-1)*int_mb(k_range+p2b-1)
     1                *int_mb(k_range+p1b-1)
     1         +(h5-1)*int_mb(k_range+p3b-1)*int_mb(k_range+p2b-1)
     1                *int_mb(k_range+p1b-1)*int_mb(k_range+h6b-1)
     1         +(h4-1)*int_mb(k_range+p3b-1)*int_mb(k_range+p2b-1)
     1                *int_mb(k_range+p1b-1)*int_mb(k_range+h6b-1)
     1                *int_mb(k_range+h5b-1)) = dbl_mb(k_r2+i)
                i=i+1
             enddo
             enddo
             enddo
             enddo
             enddo
             enddo
c
             call put_hash_block(d_y3,dbl_mb(k_r2_sort),size,
     1                int_mb(k_y3_offset),
     1                (p3b-noab-1)
     1              + (p2b-noab-1)*nvab
     1              + (p1b-noab-1)*nvab*nvab
     1              + (h6b-1)*nvab*nvab*nvab
     1              + (h5b-1)*nvab*nvab*nvab*noab
     1              + (h4b-1)*nvab*nvab*nvab*noab*noab)
c
            if (.not.ma_pop_stack(l_r2_sort))
     1      call errquit('trans_t3_to_t3dagger: ma problem',1,ma_err)
c
            if (.not.ma_pop_stack(l_r2))
     1      call errquit('trans_t3_to_t3dagger: ma problem',1,ma_err)
c
         endif
         endif
         endif
c
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
c
      end
c $Id$
