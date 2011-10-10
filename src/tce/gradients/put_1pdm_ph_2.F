      subroutine put_1pdm_ph_2(d_1pdm,k_1pdm_offset,d_ph,k_ph_offset)
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
      integer d_ph,k_ph_offset
      integer p1b,h2b
      integer p1,h2
      integer size
      integer l_r1,k_r1,l_r1_sort,k_r1_sort
      integer i
c
      integer nprocs
      integer next
      integer count
      integer nxtask
      external nxtask
c
      nprocs=ga_nnodes( )
      count = 0
      next = nxtask(nprocs,1)
c
      do p1b = noab+1,noab+nvab
      do h2b = 1,noab
         if(next.eq.count) then
c
         if (int_mb(k_spin+p1b-1) .eq. int_mb(k_spin+h2b-1)) then
         if ((.not.restricted).or.(int_mb(k_spin+p1b-1)
     1      +int_mb(k_spin+h2b-1).ne.4)) then
         if (ieor(int_mb(k_sym+p1b-1),int_mb(k_sym+h2b-1))
     1      .eq.irrep_e) then
            size = int_mb(k_range+p1b-1) * int_mb(k_range+h2b-1)
            if (.not.ma_push_get(mt_dbl,size,'r1',l_r1,k_r1))
     1        call errquit('put_1pdm_ph: MA problem',0,MA_ERR)
c
            call get_hash_block(d_ph,dbl_mb(k_r1),size,
     1        int_mb(k_ph_offset),((p1b-noab-1)*noab+h2b-1))
            call put_hash_block(d_1pdm,dbl_mb(k_r1),size,
     1        int_mb(k_1pdm_offset),(h2b-1)+(p1b-1)*(noab+nvab))
c
            if (.not.ma_pop_stack(l_r1))
     1      call errquit('put_1pdm_ph: MA problem',1,MA_ERR)
c
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
c
      next=nxtask(-nprocs,1)
      call ga_sync( )
c
      end

c $Id$
