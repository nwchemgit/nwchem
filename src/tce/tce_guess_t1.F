      subroutine tce_guess_t1(d_t1,k_t1_offset)
c
c $Id: tce_guess_t1.F,v 1.5 2003-10-19 00:52:07 sohirata Exp $
c
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "errquit.fh"
#include "tce.fh"
      integer d_t1
      integer p1b
      integer h2b
      integer k_t1_offset
      integer size
      integer l_t1,k_t1
      integer i
      integer nprocs
      integer count
      integer next
      integer nxtval
      external nxtval
c
c     ================
c     Loop over blocks
c     ================
c
      nprocs = ga_nnodes()
      count = 0
      next = nxtval(nprocs)
      do p1b = noab+1,noab+nvab
        do h2b = 1,noab
          if (next.eq.count) then
            if (int_mb(k_spin+p1b-1) .eq. int_mb(k_spin+h2b-1)) then
            if ((.not.restricted).or.(int_mb(k_spin+p1b-1)
     1        +int_mb(k_spin+h2b-1).ne.4)) then
            if (ieor(int_mb(k_sym+p1b-1),int_mb(k_sym+h2b-1)).eq.0) then
              size = int_mb(k_range+p1b-1) * int_mb(k_range+h2b-1)
              if (.not.ma_push_get(mt_dbl,size,'t1',l_t1,k_t1))
     1          call errquit('tce_guess_t1: MA problem',0,MA_ERR)
              do i = 1,size
                dbl_mb(k_t1+i-1) = 0.0d0
              enddo
              call put_block(d_t1,dbl_mb(k_t1),size,
     1          int_mb(k_t1_offset+(p1b-noab-1)*noab+h2b-1))
              if (.not.ma_pop_stack(l_t1))
     1          call errquit('tce_guess_t1: MA problem',1,MA_ERR)
            endif
            endif
            endif
            next = nxtval(nprocs)
          endif
          count = count + 1
        enddo
      enddo
      next = nxtval(-nprocs)
      call ga_sync()
      return
      end  
