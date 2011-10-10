      subroutine trans_t4_to_t4dagger(d_t4,k_t4_offset,d_y4,k_y4_offset)
      IMPLICIT NONE
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "errquit.fh"
#include "tce.fh"
c
      integer d_t4,k_t4_offset
      integer d_y4,k_y4_offset
c
      integer p1b, p2b, p3b, p4b, h5b, h6b, h7b, h8b
      integer p1,  p2,  p3,  p4,  h5,  h6,  h7,  h8
      integer k_r2,l_r2,k_r2_sort,l_r2_sort,size
      integer i
c
      do p1b = noab+1,noab+nvab
      do p2b = p1b,noab+nvab
      do p3b = p2b,noab+nvab
      do p4b = p3b,noab+nvab
      do h5b = 1,noab
      do h6b = h5b,noab
      do h7b = h6b,noab
      do h8b = h7b,noab
         if (int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1)
     1      +int_mb(k_spin+p3b-1)+int_mb(k_spin+p4b-1).eq. 
     1       int_mb(k_spin+h5b-1)+int_mb(k_spin+h6b-1)
     1      +int_mb(k_spin+h7b-1)+int_mb(k_spin+h8b-1)) then
         if (ieor(int_mb(k_sym+p1b-1),ieor(int_mb(k_sym+p2b-1),
     1       ieor(int_mb(k_sym+p3b-1),ieor(int_mb(k_sym+p4b-1),
     1       ieor(int_mb(k_sym+h5b-1),ieor(int_mb(k_sym+h6b-1),
     1       ieor(int_mb(k_sym+h7b-1),int_mb(k_sym+h8b-1)))))))).eq. 
     1       irrep_t) then
         if ((.not.restricted).or.(int_mb(k_spin+p1b-1)
     1       +int_mb(k_spin+p2b-1)+int_mb(k_spin+p3b-1)
     1       +int_mb(k_spin+p4b-1)+int_mb(k_spin+h5b-1)
     1       +int_mb(k_spin+h6b-1)+int_mb(k_spin+h7b-1)
     1       +int_mb(k_spin+h8b-1).ne.16)) then
c
             size = int_mb(k_range+p1b-1) * int_mb(k_range+p2b-1) 
     1            * int_mb(k_range+p3b-1) * int_mb(k_range+p4b-1) 
     1            * int_mb(k_range+h5b-1) * int_mb(k_range+h6b-1)
     1            * int_mb(k_range+h7b-1) * int_mb(k_range+h8b-1)
             if (.not.ma_push_get(mt_dbl,size,'r2',l_r2,k_r2))
     1       call errquit('trans_t4_to_t4dagger: ma problem',0,ma_err)
c
             if (.not.ma_push_get(mt_dbl,size,'r2_sort',l_r2_sort,
     1           k_r2_sort))
     1       call errquit('trans_t4_to_t4dagger: ma problem',0,ma_err)
c
             call get_hash_block(d_t4,dbl_mb(k_r2),size,
     1            int_mb(k_t4_offset), 
     1           (h8b-1)
     1          +(h7b-1)*noab
     1          +(h6b-1)*noab*noab
     1          +(h5b-1)*noab*noab*noab
     1          +(p4b-noab-1)*noab*noab*noab*noab
     1          +(p3b-noab-1)*noab*noab*noab*noab*nvab
     1          +(p2b-noab-1)*noab*noab*noab*noab*nvab*nvab
     1          +(p1b-noab-1)*noab*noab*noab*noab*nvab*nvab*nvab)
c
             i=0
             do p1=1,int_mb(k_range+p1b-1)
             do p2=1,int_mb(k_range+p2b-1)
             do p3=1,int_mb(k_range+p3b-1)
             do p4=1,int_mb(k_range+p4b-1)
             do h5=1,int_mb(k_range+h5b-1)
             do h6=1,int_mb(k_range+h6b-1)
             do h7=1,int_mb(k_range+h7b-1)
             do h8=1,int_mb(k_range+h8b-1)
                dbl_mb(k_r2_sort
     1         +(p4-1)
     1         +(p3-1)*int_mb(k_range+p4b-1)
     1         +(p2-1)*int_mb(k_range+p4b-1)*int_mb(k_range+p3b-1)
     1         +(p1-1)*int_mb(k_range+p4b-1)*int_mb(k_range+p3b-1)
     1                *int_mb(k_range+p2b-1)
     1         +(h8-1)*int_mb(k_range+p4b-1)*int_mb(k_range+p3b-1)
     1                *int_mb(k_range+p2b-1)*int_mb(k_range+p1b-1)
     1         +(h7-1)*int_mb(k_range+p4b-1)*int_mb(k_range+p3b-1)
     1                *int_mb(k_range+p2b-1)*int_mb(k_range+p1b-1)
     1                *int_mb(k_range+h8b-1)
     1         +(h6-1)*int_mb(k_range+p4b-1)*int_mb(k_range+p3b-1)
     1                *int_mb(k_range+p2b-1)*int_mb(k_range+p1b-1)
     1                *int_mb(k_range+h8b-1)*int_mb(k_range+h7b-1)
     1         +(h5-1)*int_mb(k_range+p4b-1)*int_mb(k_range+p3b-1)
     1                *int_mb(k_range+p2b-1)*int_mb(k_range+p1b-1)
     1                *int_mb(k_range+h8b-1)*int_mb(k_range+h7b-1)
     1                *int_mb(k_range+h6b-1)) = dbl_mb(k_r2+i)
                i=i+1
             enddo
             enddo
             enddo
             enddo
             enddo
             enddo
             enddo
             enddo
c
             call put_hash_block(d_y4,dbl_mb(k_r2_sort),size,
     1                int_mb(k_y4_offset),
     1                (p4b-noab-1)
     1              + (p3b-noab-1)*nvab
     1              + (p2b-noab-1)*nvab*nvab
     1              + (p1b-noab-1)*nvab*nvab*nvab
     1              + (h8b-1)*nvab*nvab*nvab*nvab
     1              + (h7b-1)*nvab*nvab*nvab*nvab*noab
     1              + (h6b-1)*nvab*nvab*nvab*nvab*noab*noab
     1              + (h5b-1)*nvab*nvab*nvab*nvab*noab*noab*noab)
c
            if (.not.ma_pop_stack(l_r2_sort))
     1      call errquit('trans_t4_to_t4dagger: ma problem',1,ma_err)
c
            if (.not.ma_pop_stack(l_r2))
     1      call errquit('trans_t4_to_t4dagger: ma problem',1,ma_err)
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
      enddo
      enddo
c
      end
c $Id$
