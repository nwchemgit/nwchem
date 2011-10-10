      subroutine put_1pdm_hh(d_1pdm,k_1pdm_offset,d_hh,k_hh_offset)
c
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "util.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "tce.fh"
c
      integer d_1pdm,k_1pdm_offset
      integer d_hh,k_hh_offset
      integer h1b,h2b
      integer h1,h2
      integer size
      integer l_r1,k_r1
      integer i
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
      do h2b = 1,noab
         if(next.eq.count) then
         if (int_mb(k_spin+h1b-1) .eq. int_mb(k_spin+h2b-1)) then 
         if ((.not.restricted).or.(int_mb(k_spin+h1b-1)
     1      +int_mb(k_spin+h2b-1).ne.4)) then
         if (ieor(int_mb(k_sym+h1b-1),int_mb(k_sym+h2b-1))
     1      .eq.irrep_e) then
            size = int_mb(k_range+h1b-1) * int_mb(k_range+h2b-1)
            if (.not.ma_push_get(mt_dbl,size,'r1',l_r1,k_r1))
     1        call errquit('put_1pdm_hh: ma problem',0,MA_ERR)
            call get_hash_block(d_hh,dbl_mb(k_r1),size,
     1        int_mb(k_hh_offset),((h1b-1)*noab+h2b-1))
            call put_hash_block(d_1pdm,dbl_mb(k_r1),size,
     1        int_mb(k_1pdm_offset),(h2b-1)+(h1b-1)*(noab+nvab))
c
            if (.not.ma_pop_stack(l_r1))
     1      call errquit('put_1pdm_hh: ma problem',1,MA_ERR)
c
         endif
         endif
         endif
         next=nxtask(nprocs,1)
         endif
         count=count+1
      enddo
      enddo
      next=nxtask(-nprocs,1)
      call ga_sync( )
c
      end
c $Id$
