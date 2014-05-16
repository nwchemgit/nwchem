ckbn Kiran Bhaskaran-Nair modified from tce_print_x1.F
ckbn Has to be debugged
      subroutine tce_print_ipx1(d_r1,k_r1_offset,r1_thresh,irrep)
c
c $Id: tce_print_x1.F 20265 2011-05-05 00:09:26Z d3y133 $
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
      if (nodezero) write(LuOut,9000) 'Singles contributions'
ckbn      do p1b = noab+1,noab+nvab
      do p1b = 1,1
        do h2b = 1,noab
ckbn          if (int_mb(k_spin+p1b-1) .eq. int_mb(k_spin+h2b-1)) then 
          if (ip_unused_spin  .eq. int_mb(k_spin+h2b-1)) then 
ckbn            if ((.not.restricted).or.(int_mb(k_spin+p1b-1)
ckbn     1        +int_mb(k_spin+h2b-1).ne.4)) then
            if ((.not.restricted).or.(ip_unused_spin 
     1        +int_mb(k_spin+h2b-1).ne.4)) then
ckbn              if (ieor(int_mb(k_sym+p1b-1),int_mb(k_sym+h2b-1))
ckbn     1          .eq.irrep) then
              if (ieor(ip_unused_sym ,int_mb(k_sym+h2b-1))
     1          .eq.irrep) then
ckbn                size = int_mb(k_range+p1b-1) * int_mb(k_range+h2b-1)
                size = 1 * int_mb(k_range+h2b-1)
                if (.not.ma_push_get(mt_dbl,size,'r1',l_r1,k_r1))
     1            call errquit('tce_print_x1: MA problem',0,MA_ERR)
                call get_hash_block(d_r1,dbl_mb(k_r1),size,
     1            int_mb(k_r1_offset),((p1b-noab-1)*noab+h2b-1))
                i = 0
ckbn                do p1 = 1,int_mb(k_range+p1b-1)
                do p1 = 1,1
                  do h2 = 1,int_mb(k_range+h2b-1)
                    i = i + 1
                    if (dabs(dbl_mb(k_r1+i-1)).gt.r1_thresh) then
ckbn                      call sym_irrepname(geom,int_mb(k_sym+p1b-1)+1,
ckbn     1                                   irrep_p)
                      call sym_irrepname(geom,ip_unused_sym+1,
     1                                   irrep_p)
                      call sym_irrepname(geom,int_mb(k_sym+h2b-1)+1,
     1                                   irrep_h)
ckbn                      spin_p = spinnames(int_mb(k_spin+p1b-1))
                      spin_p = spinnames(ip_unused_spin)
                      spin_h = spinnames(int_mb(k_spin+h2b-1))
                      if (nodezero) write(LuOut,9010) 
c     1                int_mb(k_offset+p1b-1)+p1,irrep_p,spin_p,
c     2                int_mb(k_offset+h2b-1)+h2,irrep_h,spin_h,
c     3                dbl_mb(k_r1+i-1)
ckbn     1                (int_mb(k_mo_index+int_mb(k_offset+p1b-1)+p1-1)+
ckbn     1                 2-int_mb(k_spin+p1b-1))/2,
ckbn     1                irrep_p,spin_p,
ckbn     2                (int_mb(k_mo_index+int_mb(k_offset+h2b-1)+h2-1)+
ckbn     1                 2-int_mb(k_spin+h2b-1))/2,
ckbn     2                irrep_h,spin_h,
ckbn     2                dbl_mb(k_r1+i-1)
     1                (int_mb(k_mo_index+int_mb(k_offset+p1b-1)+p1-1)+
     1                 2-ip_unused_spin)/2,
     1                irrep_p,spin_p,
     2                (int_mb(k_mo_index+int_mb(k_offset+h2b-1)+h2-1)+
     1                 2-int_mb(k_spin+h2b-1))/2,
     2                irrep_h,spin_h,
     2                dbl_mb(k_r1+i-1)
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
 9010 format(1x,i5,a4,'(',a5,') --- ',i5,a4,'(',a5,')',f20.10)
      return
      end


