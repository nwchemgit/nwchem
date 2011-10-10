      subroutine put_2pdm_hphh(d_2pdm,k_2pdm_offset,
     1                         d_hphh,k_hphh_offset)
      IMPLICIT NONE
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "errquit.fh"
#include "tce.fh"
c
      integer d_2pdm,k_2pdm_offset
      integer d_hphh,k_hphh_offset
      integer i
c
      integer h1b,p2b,h3b,h4b
      integer h1,p2,h3,h4
      integer k_r2,l_r2,k_r2_sort,l_r2_sort,size
c
      integer next
      integer count
      integer nxtask
      integer nprocs
      external nxtask
c
      nprocs=ga_nnodes()
      count=0
      next=nxtask(nprocs,1)
c
      do h1b = 1,noab
      do p2b = noab+1,noab+nvab
      do h3b = 1,noab
      do h4b = 1,noab
         if (next.eq.count) then
c
         if (int_mb(k_spin+h1b-1)+int_mb(k_spin+p2b-1)
     1       .eq. int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1)) then
         if ((.not.restricted).or.
     1      (int_mb(k_spin+h1b-1)+int_mb(k_spin+p2b-1)+
     2       int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1).ne.8)) then
         if (ieor(int_mb(k_sym+h1b-1),ieor(int_mb(k_sym+p2b-1),
     1       ieor(int_mb(k_sym+h3b-1),int_mb(k_sym+h4b-1))))
     2       .eq. irrep_e) then
             size = int_mb(k_range+h1b-1) * int_mb(k_range+p2b-1)
     1            * int_mb(k_range+h3b-1) * int_mb(k_range+h4b-1)
             if (.not.ma_push_get(mt_dbl,size,'r2',l_r2,k_r2))
     1       call errquit('put_2pdm_hphh: MA problem',0,MA_ERR)
c
             if (.not.ma_push_get(mt_dbl,size,'r2_sort',l_r2_sort,
     1           k_r2_sort))
     1       call errquit('put_2pdm_hphh: MA problem',0,MA_ERR)

c
             if (h4b.ge.h3b) then
                call get_hash_block(d_hphh,dbl_mb(k_r2),size,
     1               int_mb(k_hphh_offset),((((h1b-1)*nvab+
     1               p2b-noab-1)*noab+h3b-1)*noab+h4b-1))
c               
                call put_hash_block(d_2pdm,dbl_mb(k_r2),size,
     1               int_mb(k_2pdm_offset),
     1               (h4b-1)+(h3b-1)*(noab+nvab)+
     1               (p2b-1)*(noab+nvab)*(noab+nvab)+
     1               (h1b-1)*(noab+nvab)*(noab+nvab)*(noab+nvab))
c
                i=0
                do h1=1,int_mb(k_range+h1b-1)
                do p2=1,int_mb(k_range+p2b-1)
                do h3=1,int_mb(k_range+h3b-1)
                do h4=1,int_mb(k_range+h4b-1)
                   dbl_mb(k_r2_sort+(h4-1)+(h3-1)*int_mb(k_range+h4b-1)+
     1             (h1-1)*int_mb(k_range+h4b-1)*int_mb(k_range+h3b-1)+
     1             (p2-1)*int_mb(k_range+h4b-1)*int_mb(k_range+h3b-1)*
     1             int_mb(k_range+h1b-1)) = - dbl_mb(k_r2+i)
                   i=i+1
                enddo
                enddo
                enddo
                enddo
                call put_hash_block(d_2pdm,dbl_mb(k_r2_sort),size,
     1               int_mb(k_2pdm_offset),
     1               (h4b-1)+(h3b-1)*(noab+nvab)+
     1               (h1b-1)*(noab+nvab)*(noab+nvab)+
     1               (p2b-1)*(noab+nvab)*(noab+nvab)*(noab+nvab))
             endif
c
             if (h4b.lt.h3b) then
                call get_hash_block(d_hphh,dbl_mb(k_r2),size,
     1               int_mb(k_hphh_offset),((((h1b-1)*nvab+
     1               p2b-noab-1)*noab+h4b-1)*noab+h3b-1))
c
                i=0
                do h1=1,int_mb(k_range+h1b-1)
                do p2=1,int_mb(k_range+p2b-1)
                do h4=1,int_mb(k_range+h4b-1)
                do h3=1,int_mb(k_range+h3b-1)
                   dbl_mb(k_r2_sort+h4-1+(h3-1)*int_mb(k_range+h4b-1)+
     1             (p2-1)*int_mb(k_range+h4b-1)*int_mb(k_range+h3b-1)+
     1             (h1-1)*int_mb(k_range+h4b-1)*int_mb(k_range+h3b-1)*
     1             int_mb(k_range+p2b-1)) = - dbl_mb(k_r2+i)
                   i=i+1
                enddo
                enddo
                enddo
                enddo
c
                call put_hash_block(d_2pdm,dbl_mb(k_r2_sort),size,
     1               int_mb(k_2pdm_offset),
     1               (h4b-1)+(h3b-1)*(noab+nvab)+
     1               (p2b-1)*(noab+nvab)*(noab+nvab)+
     1               (h1b-1)*(noab+nvab)*(noab+nvab)*(noab+nvab))
c
                i=0
                do h1=1,int_mb(k_range+h1b-1)
                do p2=1,int_mb(k_range+p2b-1)
                do h4=1,int_mb(k_range+h4b-1)
                do h3=1,int_mb(k_range+h3b-1)
                   dbl_mb(k_r2_sort+h4-1+(h3-1)*int_mb(k_range+h4b-1)+
     1             (h1-1)*int_mb(k_range+h4b-1)*int_mb(k_range+h3b-1)+
     1             (p2-1)*int_mb(k_range+h4b-1)*int_mb(k_range+h3b-1)*
     1             int_mb(k_range+h1b-1)) =  dbl_mb(k_r2+i)
                   i=i+1
                enddo
                enddo
                enddo
                enddo
c
                call put_hash_block(d_2pdm,dbl_mb(k_r2_sort),size,
     1               int_mb(k_2pdm_offset),
     1               (h4b-1)+(h3b-1)*(noab+nvab)+
     1               (h1b-1)*(noab+nvab)*(noab+nvab)+
     1               (p2b-1)*(noab+nvab)*(noab+nvab)*(noab+nvab))

             endif
c
            if (.not.ma_pop_stack(l_r2_sort))
     1      call errquit('put_1pdm_hh: MA problem',1,MA_ERR)
c
            if (.not.ma_pop_stack(l_r2))
     1      call errquit('put_1pdm_hh: MA problem',1,MA_ERR)
         endif
         endif
         endif
c
         next=nxtask(nprocs,1)
         endif
         count=count+1
c
      enddo
      enddo
      enddo
      enddo
      next=nxtask(-nprocs,1)
      call ga_sync( )
c
      end
c $Id$
