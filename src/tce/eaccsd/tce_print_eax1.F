ckbn Kiran Bhaskaran-Nair modified from tce_print_x1.F
ckbn Has to be debugged
      subroutine tce_print_eax1(d_r1,k_r1_offset,r1_thresh,irrep)
c
c $Id$
c
c Prints the X1 tensor
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
      integer d_r1
      integer p1b
      integer h2b
      integer p1
      integer h2
      integer k_r1_offset
      integer size
      integer l_r1,k_r1
      integer i
      integer irrep
      logical nodezero
      double precision r1_thresh
      character*4 irrep_p,irrep_h
      character*5 spinnames(2),spin_p,spin_h
      data spinnames/'alpha','beta '/
c
      nodezero = (ga_nodeid().eq.0)
      if (nodezero) write(LuOut,9000) 'Singles contributions EA-EOMCCSD'
      do p1b = noab+1,noab+nvab
        do h2b = 1,1
          if (int_mb(k_spin+p1b-1) .eq. ea_unused_spin) then 
            if ((.not.restricted).or.(int_mb(k_spin+p1b-1)
     1        + ea_unused_spin .ne.4)) then
              if (ieor(int_mb(k_sym+p1b-1),ea_unused_sym)
     1          .eq.irrep) then
                size = int_mb(k_range+p1b-1) * 1
                if (.not.ma_push_get(mt_dbl,size,'r1',l_r1,k_r1))
     1            call errquit('tce_print_x1: MA problem',0,MA_ERR)
                call get_hash_block(d_r1,dbl_mb(k_r1),size,
     1            int_mb(k_r1_offset),((p1b-noab-1)*noab+h2b-1))
                i = 0
                do p1 = 1,int_mb(k_range+p1b-1)
ckbn                  do h2 = 1,int_mb(k_range+h2b-1)
                  do h2 = 1,1
                    i = i + 1
                    if (dabs(dbl_mb(k_r1+i-1)).gt.r1_thresh) then
                      call sym_irrepname(geom,int_mb(k_sym+p1b-1)+1,
     1                                   irrep_p)
                      spin_p = spinnames(int_mb(k_spin+p1b-1))
                      if (nodezero) write(LuOut,9010) 
     1                (int_mb(k_mo_index+int_mb(k_offset+p1b-1)+p1-1)+
     1                 2-int_mb(k_spin+p1b-1))/2,
     1                irrep_p,spin_p,
     2                dbl_mb(k_r1+i-1)
c#endif
                    endif
                  enddo
                enddo
                if (.not.ma_pop_stack(l_r1))
     1            call errquit('tce_print_x1: MA problem',1,MA_ERR)
              endif
            endif
          endif
        enddo
      enddo
      call ga_sync()
 9000 format(/,1x,a)
 9010 format(1x,i5,a4,'(',a5,')',f20.10)
      return
      end
