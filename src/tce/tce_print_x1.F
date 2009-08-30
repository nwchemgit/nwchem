      subroutine tce_print_x1(d_r1,k_r1_offset,r1_thresh,irrep)
c
c $Id: tce_print_x1.F,v 1.11 2009-02-25 00:46:20 d3p975 Exp $
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
      do p1b = noab+1,noab+nvab
        do h2b = 1,noab
          if (int_mb(k_spin+p1b-1) .eq. int_mb(k_spin+h2b-1)) then 
          if ((.not.restricted).or.(int_mb(k_spin+p1b-1)
     1      +int_mb(k_spin+h2b-1).ne.4)) then
          if (ieor(int_mb(k_sym+p1b-1),int_mb(k_sym+h2b-1))
     1      .eq.irrep) then
            size = int_mb(k_range+p1b-1) * int_mb(k_range+h2b-1)
            if (.not.ma_push_get(mt_dbl,size,'r1',l_r1,k_r1))
     1        call errquit('tce_print_x1: MA problem',0,MA_ERR)
            call get_hash_block(d_r1,dbl_mb(k_r1),size,
     1        int_mb(k_r1_offset),((p1b-noab-1)*noab+h2b-1))
            i = 0
            do p1 = 1,int_mb(k_range+p1b-1)
              do h2 = 1,int_mb(k_range+h2b-1)
                i = i + 1
                if (dabs(dbl_mb(k_r1+i-1)).gt.r1_thresh) then
                  call sym_irrepname(geom,int_mb(k_sym+p1b-1)+1,irrep_p)
                  call sym_irrepname(geom,int_mb(k_sym+h2b-1)+1,irrep_h)
                  spin_p = spinnames(int_mb(k_spin+p1b-1))
                  spin_h = spinnames(int_mb(k_spin+h2b-1))
                  if (nodezero) write(LuOut,9010) 
c     1              int_mb(k_offset+p1b-1)+p1,irrep_p,spin_p,
c     2              int_mb(k_offset+h2b-1)+h2,irrep_h,spin_h,
c     3              dbl_mb(k_r1+i-1)
     1              (int_mb(k_mo_index+int_mb(k_offset+p1b-1)+p1-1)+
     1               2-int_mb(k_spin+p1b-1))/2,
     1              irrep_p,spin_p,
     2              (int_mb(k_mo_index+int_mb(k_offset+h2b-1)+h2-1)+
     1               2-int_mb(k_spin+h2b-1))/2,
     2              irrep_h,spin_h,
     2              dbl_mb(k_r1+i-1)
                endif
              enddo
            enddo
            if (.not.ma_pop_stack(l_r1))
     1        call errquit('tce_print_x1: MA problem',1,MA_ERR)
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

      subroutine tce_print_y1(d_r1,k_r1_offset,r1_thresh,irrep)
c
c Prints the Y1 tensor
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
      do h2b = 1,noab
        do p1b = noab+1,noab+nvab
          if (int_mb(k_spin+p1b-1) .eq. int_mb(k_spin+h2b-1)) then 
          if ((.not.restricted).or.(int_mb(k_spin+p1b-1)
     1      +int_mb(k_spin+h2b-1).ne.4)) then
          if (ieor(int_mb(k_sym+p1b-1),int_mb(k_sym+h2b-1))
     1      .eq.irrep) then
            size = int_mb(k_range+p1b-1) * int_mb(k_range+h2b-1)
            if (.not.ma_push_get(mt_dbl,size,'r1',l_r1,k_r1))
     1        call errquit('tce_print_y1: MA problem',0,MA_ERR)
            call get_hash_block(d_r1,dbl_mb(k_r1),size,
     1        int_mb(k_r1_offset),((h2b-1)*nvab+p1b-noab-1))
            i = 0
            do p1 = 1,int_mb(k_range+p1b-1)
              do h2 = 1,int_mb(k_range+h2b-1)
                i = i + 1
                if (dabs(dbl_mb(k_r1+i-1)).gt.r1_thresh) then
                  call sym_irrepname(geom,int_mb(k_sym+p1b-1)+1,irrep_p)
                  call sym_irrepname(geom,int_mb(k_sym+h2b-1)+1,irrep_h)
                  spin_p = spinnames(int_mb(k_spin+p1b-1))
                  spin_h = spinnames(int_mb(k_spin+h2b-1))
                  if (nodezero) write(LuOut,9010) 
     1              int_mb(k_offset+p1b-1)+p1,irrep_p,spin_p,
     2              int_mb(k_offset+h2b-1)+h2,irrep_h,spin_h,
     3              dbl_mb(k_r1+i-1)
                endif
              enddo
            enddo
            if (.not.ma_pop_stack(l_r1))
     1        call errquit('tce_print_y1: MA problem',1,MA_ERR)
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
c
c
c
c
c
      subroutine tce_print_g1(rtdb,d_r1,k_r1_offset,r1_thresh,irrep_r1)
c
c Prints the G1 tensor
c
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "util.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "rtdb.fh"
#include "tce.fh"
#include "tce_main.fh"
      integer rtdb
      integer d_r1
      integer g1b
      integer g2b
      integer g1
      integer g2
      integer k_r1_offset
      integer size
      integer l_r1,k_r1
      integer i
      integer irrep_r1
      integer eiformat
      logical nodezero
      double precision r1_thresh
      character*4 irrep_g1,irrep_g2
      character*5 spinnames(2),spin_g1,spin_g2
      data spinnames/'alpha','beta '/
c
      if (.not.rtdb_get(rtdb,'tce:eiformat',mt_int,1,eiformat)) then
        eiformat = 1
        if (.not.rtdb_put(rtdb,'tce:eiformat',mt_int,1,eiformat))
     1  call errquit('tce_init: failed to write eiformat',0,RTDB_ERR)
      endif
c
      nodezero = (ga_nodeid().eq.0)
c      if (nodezero) write(LuOut,9000) 'Singles contributions'
      do g2b = 1,noab+nvab
        do g1b = 1,noab+nvab
          if (int_mb(k_spin+g1b-1) .eq. int_mb(k_spin+g2b-1)) then
          if ((.not.restricted).or.(int_mb(k_spin+g1b-1)
     1      +int_mb(k_spin+g2b-1).ne.4)) then
          if (ieor(int_mb(k_sym+g1b-1),int_mb(k_sym+g2b-1))
     1      .eq.irrep_r1) then
            size = int_mb(k_range+g1b-1) * int_mb(k_range+g2b-1)
            if (.not.ma_push_get(mt_dbl,size,'r1',l_r1,k_r1))
     1        call errquit('tce_print_y1: MA problem',0,MA_ERR)
            call get_hash_block(d_r1,dbl_mb(k_r1),size,
     1        int_mb(k_r1_offset),(g2b-1)*(noab+nvab)+g1b-1)
            i = 0
            do g1 = 1,int_mb(k_range+g1b-1)
             do g2 = 1,int_mb(k_range+g2b-1)
               i = i + 1
               if (dabs(dbl_mb(k_r1+i-1)).gt.r1_thresh) then
                 call sym_irrepname(geom,int_mb(k_sym+g1b-1)+1,irrep_g1)
                 call sym_irrepname(geom,int_mb(k_sym+g2b-1)+1,irrep_g2)
                 spin_g1 = spinnames(int_mb(k_spin+g1b-1))
                 spin_g2 = spinnames(int_mb(k_spin+g2b-1))
                 if ((nodezero).and.(eiformat.eq.1)) then
                   write(LuOut,9011)
     1             int_mb(k_offset+g1b-1)+g1,irrep_g1,spin_g1,
     2             int_mb(k_offset+g2b-1)+g2,irrep_g2,spin_g2,
     3             dbl_mb(k_r1+i-1)
                 elseif ((nodezero).and.(eiformat.eq.2)) then
                   write(LuOut,9012)
     1             int_mb(k_offset+g1b-1)+g1,irrep_g1,spin_g1,
     2             int_mb(k_offset+g2b-1)+g2,irrep_g2,spin_g2,
     3             dbl_mb(k_r1+i-1)
                 elseif ((nodezero).and.(eiformat.eq.3)) then
                   write(LuOut,9013)
     1             int_mb(k_offset+g1b-1)+g1,int_mb(k_sym+g1b-1)+1,
     2             int_mb(k_spin+g1b-1),
     3             int_mb(k_offset+g2b-1)+g2,int_mb(k_sym+g2b-1)+1,
     4             int_mb(k_spin+g2b-1),
     5             dbl_mb(k_r1+i-1)
                 elseif ((nodezero).and.(eiformat.eq.4)) then
                   write(LuOut,9014)
     1             int_mb(k_offset+g1b-1)+g1,
     2             int_mb(k_offset+g2b-1)+g2,
     3             dbl_mb(k_r1+i-1)          
                 endif
               endif
             enddo
            enddo
            if (.not.ma_pop_stack(l_r1))
     1        call errquit('tce_print_y1: MA problem',1,MA_ERR)
          endif
          endif
          endif
        enddo
      enddo
      call ga_sync()
 9000 format(/,1x,a)
 9011 format(1x,i3,a4,'(',a5,') --- ',i3,a4,'(',a5,')',f20.14)
 9012 format(1x,i3,',  ',a4,',',a5,',',i3,',  ',a4,',',a5,',',f20.14)
 9013 format(1x,i3,',',i3,',',i3,',',i3,',',i3,',',i3,',',f20.14)
 9014 format(1x,i3,', ',i3,', ',f20.14)
      return
      end
