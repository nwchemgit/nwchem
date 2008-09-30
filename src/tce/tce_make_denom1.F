      subroutine tce_make_denom1(d_denom1,k_denom1_offset,irrep,
     1                           denom_power,omega,shift)
c
c $Id: tce_make_denom1.F,v 1.2 2008-09-30 18:43:17 jhammond Exp $
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
      integer d_denom1
      integer p1b
      integer h2b
      integer p1
      integer h2
      integer k_denom1_offset
      integer size
      integer l_denom1,k_denom1
      integer i,irrep
      integer nprocs
      integer count
      integer next
      INTEGER NXTASK
      EXTERNAL NXTASK
      integer denom_power
      double precision shift,omega
      logical nodezero
      logical noloadbalance
c
c     ================
c     Loop over blocks
c     ================
c
      nodezero = (ga_nodeid().eq.0)
      noloadbalance = ((ioalg.eq.4).or.
     1                ((ioalg.eq.6).and.(.not.fileisga(d_denom1))))
      nprocs = ga_nnodes()
      count = 0
      next = NXTASK(nprocs, 1)
      do p1b = noab+1,noab+nvab
        do h2b = 1,noab
          if (noloadbalance.or.(next.eq.count)) then
            if (int_mb(k_spin+p1b-1) .eq. int_mb(k_spin+h2b-1)) then 
            if ((.not.restricted).or.(int_mb(k_spin+p1b-1)
     1        +int_mb(k_spin+h2b-1).ne.4)) then
            if (ieor(int_mb(k_sym+p1b-1),int_mb(k_sym+h2b-1))
     1                                   .eq.irrep) then
              size = int_mb(k_range+p1b-1) * int_mb(k_range+h2b-1)
              if (.not.ma_push_get(mt_dbl,size,'e1',l_denom1,k_denom1))
     1          call errquit('make_denom1: MA problem',0,MA_ERR)
              i = 0
              do p1 = 1,int_mb(k_range+p1b-1)
                do h2 = 1,int_mb(k_range+h2b-1)
                   i = i + 1
c                   dbl_mb(k_denom1+i-1) = dbl_mb(k_denom1+i-1)
c     1          /  (-dbl_mb(k_evl_sorted+int_mb(k_offset+p1b-1)+p1-1)
c     2              +dbl_mb(k_evl_sorted+int_mb(k_offset+h2b-1)+h2-1))
                dbl_mb(k_denom1+i-1) =
     1          (( -dbl_mb(k_evl_sorted+int_mb(k_offset+p1b-1)+p1-1)
     2             +dbl_mb(k_evl_sorted+int_mb(k_offset+h2b-1)+h2-1) )
     3           **denom_power-omega**denom_power+shift)
                enddo
              enddo
              call put_hash_block(d_denom1,dbl_mb(k_denom1),size,
     1          int_mb(k_denom1_offset),((p1b-noab-1)*noab+h2b-1))
              if (.not.ma_pop_stack(l_denom1))
     1          call errquit('make_denom1: MA problem',1,MA_ERR)
            endif
            endif
            endif
            next = NXTASK(nprocs, 1)
          endif
          count = count + 1
        enddo
      enddo
      next = NXTASK(-nprocs, 1)
      call ga_sync()
      return
      end
