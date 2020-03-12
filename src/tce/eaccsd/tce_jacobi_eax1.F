ckbn Kiran Bhaskaran-Nair modified from tce_jacobi_x1
      subroutine tce_jacobi_eax1(d_r1,k_r1_offset)
c
c $Id$
c
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "util.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
#include "tce_diis.fh"
      integer d_r1
      integer p1b
      integer h2b
      integer p1
      integer h2
      integer k_r1_offset
      integer size
      integer l_r1,k_r1
      integer i
      integer nprocs
      integer count
      integer next
cc      integer nxtask
cc      external nxtask
      INTEGER NXTASK
      EXTERNAL NXTASK
      logical nodezero
      logical noloadbalance
c
      nodezero = (ga_nodeid().eq.0)
      noloadbalance = ((ioalg.eq.4).or. 
     1                ((ioalg.eq.6).and.(.not.fileisga(d_r1))))
      nprocs = ga_nnodes()
      count = 0
cc      next = nxtask(nprocs,1)
      next = NXTASK(nprocs, 1)
      do p1b = noab+1,noab+nvab
c#ifndef EACCSD
cc        write(LuOut,*) "I am here 1jx1"
cc        call util_flush(Luout)
c        do h2b = 1,noab
c#else
c        write(LuOut,*) "I am here 1jeax1"
c        call util_flush(Luout)
        do h2b = 1,1
c#endif
c#ifndef EACCSD
c          if (noloadbalance.or.(next.eq.count)) then
c            if (int_mb(k_spin+p1b-1) .eq. int_mb(k_spin+h2b-1)) then 
c            if ((.not.restricted).or.(int_mb(k_spin+p1b-1)
c     1        +int_mb(k_spin+h2b-1).ne.4)) then
c            if (ieor(int_mb(k_sym+p1b-1),int_mb(k_sym+h2b-1))
c     1        .eq.irrep_x) then
c#else
          if (noloadbalance.or.(next.eq.count)) then
            if (int_mb(k_spin+p1b-1) .eq. ea_unused_spin) then 
            if ((.not.restricted).or.(int_mb(k_spin+p1b-1)
     1        +ea_unused_spin.ne.4)) then
            if (ieor(int_mb(k_sym+p1b-1),ea_unused_sym)
     1        .eq.irrep_x) then
c#endif
c#ifndef EACCSD
c              size = int_mb(k_range+p1b-1) * int_mb(k_range+h2b-1)
c#else
              size = int_mb(k_range+p1b-1) * 1
c#endif
              if (.not.ma_push_get(mt_dbl,size,'r1',l_r1,k_r1))
     1          call errquit('tce_jacobi_x1: MA problem',0,MA_ERR)
              call get_hash_block(d_r1,dbl_mb(k_r1),size,
     1          int_mb(k_r1_offset),((p1b-noab-1)*noab+h2b-1))
              i = 0
              do p1 = 1,int_mb(k_range+p1b-1)
c#ifndef EACCSD
c                do h2 = 1,int_mb(k_range+h2b-1)
c#else
c                do h2 = 1,int_mb(k_range+h2b-1)
c#endif
                  i = i + 1
                  dbl_mb(k_r1+i-1) = dbl_mb(k_r1+i-1)
c#ifndef EACCSD
c     1           / (-dbl_mb(k_evl_sorted+int_mb(k_offset+p1b-1)+p1-1)
c     2              +dbl_mb(k_evl_sorted+int_mb(k_offset+h2b-1)+h2-1))
c#else
     1           / (-dbl_mb(k_evl_sorted+int_mb(k_offset+p1b-1)+p1-1))
c     2              +dbl_mb(k_evl_sorted+int_mb(k_offset+h2b-1)+h2-1))
c#endif
c#ifndef EACCSD
c                enddo
c#else
c                enddo
c#endif
              enddo
              call put_hash_block(d_r1,dbl_mb(k_r1),size,
     1          int_mb(k_r1_offset),((p1b-noab-1)*noab+h2b-1))
              if (.not.ma_pop_stack(l_r1))
     1          call errquit('tce_jacobi_x1: MA problem',1,MA_ERR)
            endif
            endif
            endif
cc            next = nxtask(nprocs,1)
      next = NXTASK(nprocs, 1)
          endif
          count = count + 1
        enddo
      enddo
cc      next = nxtask(-nprocs,1)
      next = NXTASK(-nprocs, 1)
      call ga_sync()
      return
      end
