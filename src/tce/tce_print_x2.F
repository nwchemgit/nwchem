      subroutine tce_print_x2(d_r2,k_r2_offset,r2_thresh,irrep)
c
c $Id$
c
c Prints the X2 tensor
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
      integer d_r2
      integer p1b
      integer p2b
      integer h3b
      integer h4b
      integer p1
      integer p2
      integer h3
      integer h4
      integer k_r2_offset
      integer size
      integer l_r2,k_r2
      integer i
      integer irrep
      logical nodezero
      double precision r2_thresh
      integer p1b_in,h3b_in,p2b_in
      character*4 irrep_p1,irrep_p2,irrep_h3,irrep_h4
      character*5 spinnames(2),spin_p1,spin_p2,spin_h3,spin_h4
      data spinnames/'alpha','beta '/
c
      nodezero = (ga_nodeid().eq.0)
      if (nodezero) then
      write(LuOut,9000) 'Doubles contributions'
      do p1b = noab+1,noab+nvab
        do p2b = p1b,noab+nvab
          do h3b = 1,noab
            do h4b = h3b,noab
              if (int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1) 
     1          .eq. int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1)) then
              if ((.not.restricted).or.
     1          (int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1)+
     2          int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1).ne.8)) then
              if (ieor(int_mb(k_sym+p1b-1),ieor(int_mb(k_sym+p2b-1),
     1          ieor(int_mb(k_sym+h3b-1),int_mb(k_sym+h4b-1))))
     2          .eq. irrep) then
                size = int_mb(k_range+p1b-1) * int_mb(k_range+p2b-1)
     1               * int_mb(k_range+h3b-1) * int_mb(k_range+h4b-1)
                if (.not.ma_push_get(mt_dbl,size,'r2',l_r2,k_r2))
     1            call errquit('tce_print_x2: MA problem',0,MA_ERR)
                call get_hash_block(d_r2,dbl_mb(k_r2),size,
     1            int_mb(k_r2_offset),((((p1b-noab-1)*nvab+p2b-noab-1)
     2            *noab+h3b-1)*noab+h4b-1))
                i = 0
                do p1 = 1,int_mb(k_range+p1b-1)
                  do p2 = 1,int_mb(k_range+p2b-1)
                    do h3 = 1,int_mb(k_range+h3b-1)
                      do h4 = 1,int_mb(k_range+h4b-1)
                        i = i + 1
              if (dabs(dbl_mb(k_r2+i-1)).gt.r2_thresh) then
                call sym_irrepname(geom,int_mb(k_sym+p1b-1)+1,irrep_p1)
                call sym_irrepname(geom,int_mb(k_sym+p2b-1)+1,irrep_p2)
                call sym_irrepname(geom,int_mb(k_sym+h3b-1)+1,irrep_h3)
                call sym_irrepname(geom,int_mb(k_sym+h4b-1)+1,irrep_h4)
                spin_p1 = spinnames(int_mb(k_spin+p1b-1))
                spin_p2 = spinnames(int_mb(k_spin+p2b-1))
                spin_h3 = spinnames(int_mb(k_spin+h3b-1))
                spin_h4 = spinnames(int_mb(k_spin+h4b-1))
                if (nodezero) write(LuOut,9010)
c     1            int_mb(k_offset+p1b-1)+p1,irrep_p1,spin_p1,
c     2            int_mb(k_offset+p2b-1)+p2,irrep_p2,spin_p2,
c     3            int_mb(k_offset+h3b-1)+h3,irrep_h3,spin_h3,
c     4            int_mb(k_offset+h4b-1)+h4,irrep_h4,spin_h4,
c     5            dbl_mb(k_r2+i-1)
     1            (int_mb(k_mo_index+int_mb(k_offset+p1b-1)+p1-1)+
     1             2-int_mb(k_spin+p1b-1))/2,
     1            irrep_p1,spin_p1,
     2            (int_mb(k_mo_index+int_mb(k_offset+p2b-1)+p2-1)+
     2             2-int_mb(k_spin+p2b-1))/2,
     2            irrep_p2,spin_p2,
     3            (int_mb(k_mo_index+int_mb(k_offset+h3b-1)+h3-1)+
     3             2-int_mb(k_spin+h3b-1))/2,
     3            irrep_h3,spin_h3,
     4            (int_mb(k_mo_index+int_mb(k_offset+h4b-1)+h4-1)+
     4             2-int_mb(k_spin+h4b-1))/2,
     4            irrep_h4,spin_h4,
     5            dbl_mb(k_r2+i-1)
              endif
                      enddo
                    enddo
                  enddo
                enddo
                if (.not.ma_pop_stack(l_r2))
     1            call errquit('tce_print_x2: MA problem',1,MA_ERR)
              endif
              endif
              endif
            enddo
          enddo
        enddo
      enddo
      endif
      call ga_sync()
 9000 format(/,1x,a)
 9010 format(1x,i5,a4,'(',a5,') ',i5,a4,'(',a5,') --- ',
     1          i5,a4,'(',a5,') ',i5,a4,'(',a5,')',f20.10)
      return
      end
c
c
      subroutine tce_print_x2a(d_r2,k_r2_offset,r2_thresh,irrep)
c
c $Id$
c
c Prints the X2 tensor
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
      integer d_r2
      integer p1b
      integer p2b
      integer h3b
      integer h4b
      integer p1
      integer p2
      integer h3
      integer h4
      integer k_r2_offset
      integer size
      integer l_r2,k_r2
      integer i
      integer irrep
      logical nodezero
      double precision r2_thresh
      LOGICAL is_active_1,is_active_2,is_active_3,is_active_4
      character*4 irrep_p1,irrep_p2,irrep_h3,irrep_h4
      character*5 spinnames(2),spin_p1,spin_p2,spin_h3,spin_h4
      data spinnames/'alpha','beta '/
c
      nodezero = (ga_nodeid().eq.0)
      if (nodezero) then
         write(LuOut,9000) 'Doubles contributions'
      do p1b = noab+1,noab+nvab
        do p2b = p1b,noab+nvab
          do h3b = 1,noab
            do h4b = h3b,noab
              IF(is_active_4(p1b,p2b,h3b,h4b)) THEN
              if (int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1) 
     1          .eq. int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1)) then
              if ((.not.restricted).or.
     1          (int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1)+
     2          int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1).ne.8)) then
              if (ieor(int_mb(k_sym+p1b-1),ieor(int_mb(k_sym+p2b-1),
     1          ieor(int_mb(k_sym+h3b-1),int_mb(k_sym+h4b-1))))
     2          .eq. irrep) then
                size = int_mb(k_range+p1b-1) * int_mb(k_range+p2b-1)
     1               * int_mb(k_range+h3b-1) * int_mb(k_range+h4b-1)
                if (.not.ma_push_get(mt_dbl,size,'r2',l_r2,k_r2))
     1            call errquit('tce_print_x2: MA problem',0,MA_ERR)
                call get_hash_block(d_r2,dbl_mb(k_r2),size,
     1            int_mb(k_r2_offset),((((p1b-noab-1)*nvab+p2b-noab-1)
     2            *noab+h3b-1)*noab+h4b-1))
                i = 0
                do p1 = 1,int_mb(k_range+p1b-1)
                  do p2 = 1,int_mb(k_range+p2b-1)
                    do h3 = 1,int_mb(k_range+h3b-1)
                      do h4 = 1,int_mb(k_range+h4b-1)
                        i = i + 1
              if (dabs(dbl_mb(k_r2+i-1)).gt.r2_thresh) then
                call sym_irrepname(geom,int_mb(k_sym+p1b-1)+1,irrep_p1)
                call sym_irrepname(geom,int_mb(k_sym+p2b-1)+1,irrep_p2)
                call sym_irrepname(geom,int_mb(k_sym+h3b-1)+1,irrep_h3)
                call sym_irrepname(geom,int_mb(k_sym+h4b-1)+1,irrep_h4)
                spin_p1 = spinnames(int_mb(k_spin+p1b-1))
                spin_p2 = spinnames(int_mb(k_spin+p2b-1))
                spin_h3 = spinnames(int_mb(k_spin+h3b-1))
                spin_h4 = spinnames(int_mb(k_spin+h4b-1))
                if (nodezero) write(LuOut,9010)
c     1            int_mb(k_offset+p1b-1)+p1,irrep_p1,spin_p1,
c     2            int_mb(k_offset+p2b-1)+p2,irrep_p2,spin_p2,
c     3            int_mb(k_offset+h3b-1)+h3,irrep_h3,spin_h3,
c     4            int_mb(k_offset+h4b-1)+h4,irrep_h4,spin_h4,
c     5            dbl_mb(k_r2+i-1)
     1            (int_mb(k_mo_index+int_mb(k_offset+p1b-1)+p1-1)+
     1             2-int_mb(k_spin+p1b-1))/2,
     1            irrep_p1,spin_p1,
     2            (int_mb(k_mo_index+int_mb(k_offset+p2b-1)+p2-1)+
     2             2-int_mb(k_spin+p2b-1))/2,
     2            irrep_p2,spin_p2,
     3            (int_mb(k_mo_index+int_mb(k_offset+h3b-1)+h3-1)+
     3             2-int_mb(k_spin+h3b-1))/2,
     3            irrep_h3,spin_h3,
     4            (int_mb(k_mo_index+int_mb(k_offset+h4b-1)+h4-1)+
     4             2-int_mb(k_spin+h4b-1))/2,
     4            irrep_h4,spin_h4,
     5            dbl_mb(k_r2+i-1)
              endif
                      enddo
                    enddo
                  enddo
                enddo
                if (.not.ma_pop_stack(l_r2))
     1            call errquit('tce_print_x2: MA problem',1,MA_ERR)
              endif
              endif
              endif
              endif
            enddo
          enddo
        enddo
      enddo
      endif
      call ga_sync()
 9000 format(/,1x,a)
 9010 format(1x,i5,a4,'(',a5,') ',i5,a4,'(',a5,') --- ',
     1          i5,a4,'(',a5,') ',i5,a4,'(',a5,')',f20.10)
      return
      end
c
c
      subroutine tce_print_y2(d_r2,k_r2_offset,r2_thresh,irrep)
c
c Prints the Y2 tensor
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
      integer d_r2
      integer p1b
      integer p2b
      integer h3b
      integer h4b
      integer p1
      integer p2
      integer h3
      integer h4
      integer k_r2_offset
      integer size
      integer l_r2,k_r2
      integer i
      integer irrep
      logical nodezero
      double precision r2_thresh
      character*4 irrep_p1,irrep_p2,irrep_h3,irrep_h4
      character*5 spinnames(2),spin_p1,spin_p2,spin_h3,spin_h4
      data spinnames/'alpha','beta '/
c
      nodezero = (ga_nodeid().eq.0)
      if (nodezero) then
         write(LuOut,9000) 'Doubles contributions'
      do h3b = 1,noab
        do h4b = h3b,noab
          do p1b = noab+1,noab+nvab
            do p2b = p1b,noab+nvab
              if (int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1) 
     1          .eq. int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1)) then
              if ((.not.restricted).or.
     1          (int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1)+
     2          int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1).ne.8)) then
              if (ieor(int_mb(k_sym+p1b-1),ieor(int_mb(k_sym+p2b-1),
     1          ieor(int_mb(k_sym+h3b-1),int_mb(k_sym+h4b-1))))
     2          .eq. irrep) then
                size = int_mb(k_range+p1b-1) * int_mb(k_range+p2b-1)
     1               * int_mb(k_range+h3b-1) * int_mb(k_range+h4b-1)
                if (.not.ma_push_get(mt_dbl,size,'r2',l_r2,k_r2))
     1            call errquit('tce_print_y2: MA problem',0,MA_ERR)
                call get_hash_block(d_r2,dbl_mb(k_r2),size,
     1            int_mb(k_r2_offset),((((h3b-1)*noab+h4b-1)
     2            *nvab+p1b-noab-1)*nvab+p2b-noab-1))
                i = 0
                do h3 = 1,int_mb(k_range+h3b-1)
                  do h4 = 1,int_mb(k_range+h4b-1)
                    do p1 = 1,int_mb(k_range+p1b-1)
                      do p2 = 1,int_mb(k_range+p2b-1)
                        i = i + 1
              if (dabs(dbl_mb(k_r2+i-1)).gt.r2_thresh) then
                call sym_irrepname(geom,int_mb(k_sym+p1b-1)+1,irrep_p1)
                call sym_irrepname(geom,int_mb(k_sym+p2b-1)+1,irrep_p2)
                call sym_irrepname(geom,int_mb(k_sym+h3b-1)+1,irrep_h3)
                call sym_irrepname(geom,int_mb(k_sym+h4b-1)+1,irrep_h4)
                spin_p1 = spinnames(int_mb(k_spin+p1b-1))
                spin_p2 = spinnames(int_mb(k_spin+p2b-1))
                spin_h3 = spinnames(int_mb(k_spin+h3b-1))
                spin_h4 = spinnames(int_mb(k_spin+h4b-1))
                if (nodezero) write(LuOut,9010)
     1            int_mb(k_offset+p1b-1)+p1,irrep_p1,spin_p1,
     2            int_mb(k_offset+p2b-1)+p2,irrep_p2,spin_p2,
     3            int_mb(k_offset+h3b-1)+h3,irrep_h3,spin_h3,
     4            int_mb(k_offset+h4b-1)+h4,irrep_h4,spin_h4,
     5            dbl_mb(k_r2+i-1)
              endif
                      enddo
                    enddo
                  enddo
                enddo
                if (.not.ma_pop_stack(l_r2))
     1            call errquit('tce_print_y2: MA problem',1,MA_ERR)
              endif
              endif
              endif
            enddo
          enddo
        enddo
      enddo
      endif
      call ga_sync()
 9000 format(/,1x,a)
 9010 format(1x,i5,a4,'(',a5,') ',i5,a4,'(',a5,') --- ',
     1          i5,a4,'(',a5,') ',i5,a4,'(',a5,')',f20.10)
      return
      end
c
c
c
c
c
      subroutine tce_print_g2(rtdb,d_r2,k_r2_offset,r2_thresh,irrep_r2)
c
c Prints the G2 tensor
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
      integer d_r2
      integer g1b,g2b,g3b,g4b
      integer g1,g2,g3,g4
      integer k_r2_offset
      integer size
      integer l_r2,k_r2
      integer i
      integer irrep_r2
      integer eiformat
      logical nodezero
      double precision r2_thresh
      character*4 irrep_g1,irrep_g2,irrep_g3,irrep_g4
      character*5 spinnames(2),spin_g1,spin_g2,spin_g3,spin_g4
      data spinnames/'alpha','beta '/
c
      if (.not.rtdb_get(rtdb,'tce:eiformat',mt_int,1,eiformat)) then
        eiformat = 1
        if (.not.rtdb_put(rtdb,'tce:eiformat',mt_int,1,eiformat))
     1  call errquit('tce_init: failed to write eiformat',0,RTDB_ERR)
      endif
c
      nodezero = (ga_nodeid().eq.0)
      if(nodezero) then
c      if (nodezero) write(LuOut,9000) 'Doubles contributions'
      do g3b = 1,noab+nvab
        do g4b = g3b,noab+nvab
          do g1b = 1,noab+nvab
            do g2b = g1b,noab+nvab
              if (int_mb(k_spin+g1b-1)+int_mb(k_spin+g2b-1)
     1          .eq. int_mb(k_spin+g3b-1)+int_mb(k_spin+g4b-1)) then
              if ((.not.restricted).or.
     1          (int_mb(k_spin+g1b-1)+int_mb(k_spin+g2b-1)+
     2           int_mb(k_spin+g3b-1)+int_mb(k_spin+g4b-1).ne.8)) then
              if (ieor(int_mb(k_sym+g1b-1),ieor(int_mb(k_sym+g2b-1),
     1            ieor(int_mb(k_sym+g3b-1),int_mb(k_sym+g4b-1))))
     2          .eq. irrep_r2) then
                size = int_mb(k_range+g1b-1) * int_mb(k_range+g2b-1)
     1               * int_mb(k_range+g3b-1) * int_mb(k_range+g4b-1)
                if (.not.ma_push_get(mt_dbl,size,'r2',l_r2,k_r2))
     1            call errquit('tce_print_x2: MA problem',0,MA_ERR)
                call get_hash_block(d_r2,dbl_mb(k_r2),size,
     1            int_mb(k_r2_offset),g2b-1+(noab+nvab)
     2              *(g1b-1+(noab+nvab)*(g4b-1+(noab+nvab)*(g3b-1))))
                i = 0
                do g3 = 1,int_mb(k_range+g3b-1)
                  do g4 = 1,int_mb(k_range+g4b-1)
                    do g1 = 1,int_mb(k_range+g1b-1)
                      do g2 = 1,int_mb(k_range+g2b-1)
                        i = i + 1
              if (dabs(dbl_mb(k_r2+i-1)).gt.r2_thresh) then
                call sym_irrepname(geom,int_mb(k_sym+g3b-1)+1,irrep_g3)
                call sym_irrepname(geom,int_mb(k_sym+g4b-1)+1,irrep_g4)
                call sym_irrepname(geom,int_mb(k_sym+g1b-1)+1,irrep_g1)
                call sym_irrepname(geom,int_mb(k_sym+g2b-1)+1,irrep_g2)
                spin_g3 = spinnames(int_mb(k_spin+g3b-1))
                spin_g4 = spinnames(int_mb(k_spin+g4b-1))
                spin_g1 = spinnames(int_mb(k_spin+g1b-1))
                spin_g2 = spinnames(int_mb(k_spin+g2b-1))
                if ((nodezero).and.(eiformat.eq.1)) then
                  write(LuOut,9011)
     3            int_mb(k_offset+g3b-1)+g3,irrep_g3,spin_g3,
     4            int_mb(k_offset+g4b-1)+g4,irrep_g4,spin_g4,
     1            int_mb(k_offset+g1b-1)+g1,irrep_g1,spin_g1,
     2            int_mb(k_offset+g2b-1)+g2,irrep_g2,spin_g2,
     5            dbl_mb(k_r2+i-1)
                elseif ((nodezero).and.(eiformat.eq.2)) then
                  write(LuOut,9012)
     3            int_mb(k_offset+g3b-1)+g3,irrep_g3,spin_g3,
     4            int_mb(k_offset+g4b-1)+g4,irrep_g4,spin_g4,
     1            int_mb(k_offset+g1b-1)+g1,irrep_g1,spin_g1,
     2            int_mb(k_offset+g2b-1)+g2,irrep_g2,spin_g2,
     5            dbl_mb(k_r2+i-1)
                elseif ((nodezero).and.(eiformat.eq.3)) then
                  write(LuOut,9013)
     1             int_mb(k_offset+g3b-1)+g3,int_mb(k_sym+g3b-1)+1,
     2             int_mb(k_spin+g3b-1),
     3             int_mb(k_offset+g4b-1)+g4,int_mb(k_sym+g4b-1)+1,
     4             int_mb(k_spin+g4b-1),
     1             int_mb(k_offset+g1b-1)+g1,int_mb(k_sym+g1b-1)+1,
     2             int_mb(k_spin+g1b-1),
     3             int_mb(k_offset+g2b-1)+g2,int_mb(k_sym+g2b-1)+1,
     4             int_mb(k_spin+g2b-1),
     5            dbl_mb(k_r2+i-1)
                elseif ((nodezero).and.(eiformat.eq.4)) then
                  write(LuOut,9014)
     3            int_mb(k_offset+g3b-1)+g3,
     4            int_mb(k_offset+g4b-1)+g4,
     1            int_mb(k_offset+g1b-1)+g1,
     2            int_mb(k_offset+g2b-1)+g2,
     5            dbl_mb(k_r2+i-1)
                endif
              endif
                      enddo
                    enddo
                  enddo
                enddo
                if (.not.ma_pop_stack(l_r2))
     1            call errquit('tce_print_x2: MA problem',1,MA_ERR)
              endif
              endif
              endif
            enddo
          enddo
        enddo
      enddo
      endif
      call ga_sync()
 9000 format(/,1x,a)
 9011 format(1x,i3,a4,'(',a5,') ',i3,a4,'(',a5,') --- ',
     1          i3,a4,'(',a5,') ',i3,a4,'(',a5,')',f20.14)
 9012 format(1x,i3,',  ',a4,',',a5,',',i3,',  ',a4,',',a5,',',
     1          i3,',  ',a4,',',a5,',',i3,',  ',a4,',',a5,',',f20.14)
 9013 format(1x,i3,',',i3,',',i3,',',i3,',',i3,',',i3,',',
     1          i3,',',i3,',',i3,',',i3,',',i3,',',i3,',',f20.14)
 9014 format(1x,i3,', ',i3,', ',i3,', ',i3,', ',f20.14)
      return
      end
c
c
c
      subroutine tce_string_x2(d_r2,k_r2_offset,r2_thresh,irrep,xnorm)
c
c $Id: tce_print_x2.F 26707 2015-01-17 01:06:43Z edo $
c
c Prints the X2 tensor
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
      integer d_r2
      integer p1b
      integer p2b
      integer h3b
      integer h4b
      integer p1
      integer p2
      integer h3
      integer h4
      integer k_r2_offset
      integer size
      integer l_r2,k_r2
      integer i
      integer irrep
c string
       integer spinsum,j
       integer iposh3,iposh4,iposp1,iposp2
       integer iposh3a,iposh4a,iposp1a,iposp2a
       integer nosa,nosb,nvsa,nvsb
c ------
      logical nodezero
      double precision r2_thresh,xnorm
      integer p1b_in,h3b_in,p2b_in
      character*4 irrep_p1,irrep_p2,irrep_h3,irrep_h4
      character*5 spinnames(2),spin_p1,spin_p2,spin_h3,spin_h4
c string
      character*6 refstring(204)
      character*1 bra,num1
      character*2 num2
      character*3 keta,ketb,kea,keb
      character*4 cre4
      character*5 cre5
      character*6 cre6
c ------
      data spinnames/'alpha','beta '/
c
      nodezero = (ga_nodeid().eq.0)
      if (nodezero) then
c        write(6,*)'check 3: xnorm =',xnorm
c        call util_flush(6)
c string
c
        bra='('
        keta='a)+'
        ketb='b)+'
        kea='a)'
        keb='b)'
c
        nosa=0
        nosb=0
        nvsa=0
        nvsb=0
        do i=1,noa
         nosa=nosa+int_mb(k_range+i-1)
        enddo
        do i=1,nob
         nosb=nosb+int_mb(k_range+noa+i-1)
        enddo
        do i=1,nva
         nvsa=nvsa+int_mb(k_range+noab+i-1)
        enddo
        do i=1,nvb
         nvsb=nvsb+int_mb(k_range+noab+nva+i-1)
        enddo
c reference string
        do i=1,nosa
          if(i.le.9) then
           j=i
           write(num1,301) j
           cre5=bra//num1//keta
           refstring(i+4)=cre5
          else
           j=i
           write(num2,302) j
           cre6=bra//num2//keta
           refstring(i+4)=cre6
          endif
        enddo
c
        do i=1,nosb
          if(i.le.9) then
           j=i
           write(num1,301) j
           cre5=bra//num1//ketb
           refstring(i+nosa+4)=cre5
          else
           j=i
           write(num2,302) j
           cre6=bra//num2//ketb
           refstring(i+nosa+4)=cre6
          endif
        enddo
c
c ------
      write(LuOut,9000) 'Doubles strings'
      do p1b = noab+1,noab+nvab
        do p2b = p1b,noab+nvab
          do h3b = 1,noab
            do h4b = h3b,noab
              if (int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1) 
     1          .eq. int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1)) then
c string
             spinsum=int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1)+
     1               int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1)
c
              if ((.not.restricted).or.
     1          (int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1)+
     2          int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1).ne.8)) then
              if (ieor(int_mb(k_sym+p1b-1),ieor(int_mb(k_sym+p2b-1),
     1          ieor(int_mb(k_sym+h3b-1),int_mb(k_sym+h4b-1))))
     2          .eq. irrep) then
                size = int_mb(k_range+p1b-1) * int_mb(k_range+p2b-1)
     1               * int_mb(k_range+h3b-1) * int_mb(k_range+h4b-1)
                if (.not.ma_push_get(mt_dbl,size,'r2',l_r2,k_r2))
     1            call errquit('tce_print_x2: MA problem',0,MA_ERR)
                call get_hash_block(d_r2,dbl_mb(k_r2),size,
     1            int_mb(k_r2_offset),((((p1b-noab-1)*nvab+p2b-noab-1)
     2            *noab+h3b-1)*noab+h4b-1))
                i = 0
                do p1 = 1,int_mb(k_range+p1b-1)
                  do p2 = 1,int_mb(k_range+p2b-1)
                    do h3 = 1,int_mb(k_range+h3b-1)
                      do h4 = 1,int_mb(k_range+h4b-1)
                        i = i + 1
              if (dabs(dbl_mb(k_r2+i-1)).gt.r2_thresh) then
c
c string : valid only for RHF reference
c
                 iposh3a=int_mb(k_offset+h3b-1)+h3
                 iposh4a=int_mb(k_offset+h4b-1)+h4
                 iposp1a=int_mb(k_offset+p1b-1)+p1
                 iposp2a=int_mb(k_offset+p2b-1)+p2
         if((iposp1a.le.iposp2a).and.(iposh3a.le.iposh4a))then
c p1+ p2+ h4 h3
            if(spinsum.ne.4) then ! a<b a<b ----------------------
                 iposh3=int_mb(k_offset+h3b-1)+h3
                 iposh4=int_mb(k_offset+h4b-1)+h4-nosa
                 iposp1=int_mb(k_offset+p1b-1)+p1-nosb
                 iposp2=int_mb(k_offset+p2b-1)+p2-nosb-nvsa
c p1+
                 if(iposp1.le.9) then
                  write(num1,301) iposp1
                  cre5=bra//num1//keta
                  refstring(1)=cre5
                 else
                  write(num2,302) iposp1
                  cre6=bra//num2//keta
                  refstring(1)=cre6
                 endif
c p2+
                 if(iposp2.le.9) then
                  write(num1,301) iposp2
                  cre5=bra//num1//ketb
                  refstring(2)=cre5
                 else
                  write(num2,302) iposp2
                  cre6=bra//num2//ketb
                  refstring(2)=cre6
                 endif
c  h4
                 if(iposh4.le.9) then
                  write(num1,301) iposh4
                  cre4=bra//num1//keb
                  refstring(3)=cre4
                 else
                  write(num2,302) iposh4
                  cre5=bra//num2//keb
                  refstring(3)=cre5
                 endif
c h3
                 if(iposh3.le.9) then
                  write(num1,301) iposh3
                  cre4=bra//num1//kea
                  refstring(4)=cre4
                 else
                  write(num2,302) iposh3
                  cre5=bra//num2//kea
                  refstring(4)=cre5
                 endif
c yaml initial_state_suggestion state
c yaml "Doubles contribution"
        write(6,*) xnorm*dbl_mb(k_r2+i-1),':',
     &             (refstring(j),j=1,nosa+nosb+4),'|0>'
            else  ! a<a a<a --------------------------------------
                 iposh3=int_mb(k_offset+h3b-1)+h3
                 iposh4=int_mb(k_offset+h4b-1)+h4
                 iposp1=int_mb(k_offset+p1b-1)+p1-nosb
                 iposp2=int_mb(k_offset+p2b-1)+p2-nosb
c p1+
                 if(iposp1.le.9) then
                  write(num1,301) iposp1
                  cre5=bra//num1//keta
                  refstring(1)=cre5
                 else
                  write(num2,302) iposp1
                  cre6=bra//num2//keta
                  refstring(1)=cre6
                 endif
c p2+
                 if(iposp2.le.9) then
                  write(num1,301) iposp2
                  cre5=bra//num1//keta
                  refstring(2)=cre5
                 else
                  write(num2,302) iposp2
                  cre6=bra//num2//keta
                  refstring(2)=cre6
                 endif
c  h4
                 if(iposh4.le.9) then
                  write(num1,301) iposh4
                  cre4=bra//num1//kea
                  refstring(3)=cre4
                 else
                  write(num2,302) iposh4
                  cre5=bra//num2//kea
                  refstring(3)=cre5
                 endif
c h3
                 if(iposh3.le.9) then
                  write(num1,301) iposh3
                  cre4=bra//num1//kea
                  refstring(4)=cre4
                 else
                  write(num2,302) iposh3
                  cre5=bra//num2//kea
                  refstring(4)=cre5
                 endif
c yaml initial_state_suggestion state
c yaml "Doubles contribution"
        write(6,*) xnorm*dbl_mb(k_r2+i-1),':',
     &            (refstring(j),j=1,nosa+nosb+4),'|0>'
               if(spinsum.eq.4) then  !b<b b<b
c p1+
                 if(iposp1.le.9) then
                  write(num1,301) iposp1
                  cre5=bra//num1//ketb
                  refstring(1)=cre5
                 else
                  write(num2,302) iposp1
                  cre6=bra//num2//ketb
                  refstring(1)=cre6
                 endif
c p2+
                 if(iposp2.le.9) then
                  write(num1,301) iposp2
                  cre5=bra//num1//ketb
                  refstring(2)=cre5
                 else
                  write(num2,302) iposp2
                  cre6=bra//num2//ketb
                  refstring(2)=cre6
                 endif
c  h4
                 if(iposh4.le.9) then
                  write(num1,301) iposh4
                  cre4=bra//num1//keb
                  refstring(3)=cre4
                 else
                  write(num2,302) iposh4
                  cre5=bra//num2//keb
                  refstring(3)=cre5
                 endif
c h3
                 if(iposh3.le.9) then
                  write(num1,301) iposh3
                  cre4=bra//num1//keb
                  refstring(4)=cre4
                 else
                  write(num2,302) iposh3
                  cre5=bra//num2//keb
                  refstring(4)=cre5
                 endif
c yaml initial_state_suggestion state
c yaml "Doubles contribution"
        write(6,*) xnorm*dbl_mb(k_r2+i-1),':',
     &             (refstring(j),j=1,nosa+nosb+4),'|0>'
               end if !b<b b<b
            endif ! a<a a<a
c
              endif ! iposp1<iposp2 iposh3<iposh4
              endif !thresh
                      enddo
                    enddo
                  enddo
                enddo
                if (.not.ma_pop_stack(l_r2))
     1            call errquit('tce_print_x2: MA problem',1,MA_ERR)
              endif
              endif
              endif
            enddo
          enddo
        enddo
      enddo
      endif
      call ga_sync()
  301 format(i1)
  302 format(i2)
 9000 format(/,1x,a)
 9010 format(1x,i5,a4,'(',a5,') ',i5,a4,'(',a5,') --- ',
     1          i5,a4,'(',a5,') ',i5,a4,'(',a5,')',f20.10)
      return
      end
c
c
c
      subroutine tce_string_x2_norm(d_r2,k_r2_offset,r2_thresh,irrep,
     &           xnorm2)
c
c $Id: tce_print_x2.F 26707 2015-01-17 01:06:43Z edo $
c
c Prints the X2 tensor
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
      integer d_r2
      integer p1b
      integer p2b
      integer h3b
      integer h4b
      integer p1
      integer p2
      integer h3
      integer h4
      integer k_r2_offset
      integer size
      integer l_r2,k_r2
      integer i
      integer irrep
c string
       integer spinsum,j
       integer iposh3,iposh4,iposp1,iposp2
c ------
      logical nodezero
      double precision r2_thresh,xnorm2
      integer p1b_in,h3b_in,p2b_in
      character*4 irrep_p1,irrep_p2,irrep_h3,irrep_h4
      character*5 spinnames(2),spin_p1,spin_p2,spin_h3,spin_h4
c ------
      data spinnames/'alpha','beta '/
c xnorm2 defined only on a nodezero
      xnorm2=0.0d0
      nodezero = (ga_nodeid().eq.0)
      if (nodezero) then
c ------
      write(LuOut,9000) 'doubles norm'
      do p1b = noab+1,noab+nvab
        do p2b = p1b,noab+nvab
          do h3b = 1,noab
            do h4b = h3b,noab
              if (int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1) 
     1          .eq. int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1)) then
c string
             spinsum=int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1)+
     1               int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1)
c
              if ((.not.restricted).or.
     1          (int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1)+
     2          int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1).ne.8)) then
              if (ieor(int_mb(k_sym+p1b-1),ieor(int_mb(k_sym+p2b-1),
     1          ieor(int_mb(k_sym+h3b-1),int_mb(k_sym+h4b-1))))
     2          .eq. irrep) then
                size = int_mb(k_range+p1b-1) * int_mb(k_range+p2b-1)
     1               * int_mb(k_range+h3b-1) * int_mb(k_range+h4b-1)
                if (.not.ma_push_get(mt_dbl,size,'r2',l_r2,k_r2))
     1            call errquit('tce_print_x2: MA problem',0,MA_ERR)
                call get_hash_block(d_r2,dbl_mb(k_r2),size,
     1            int_mb(k_r2_offset),((((p1b-noab-1)*nvab+p2b-noab-1)
     2            *noab+h3b-1)*noab+h4b-1))
                i = 0
                do p1 = 1,int_mb(k_range+p1b-1)
                  do p2 = 1,int_mb(k_range+p2b-1)
                    do h3 = 1,int_mb(k_range+h3b-1)
                      do h4 = 1,int_mb(k_range+h4b-1)
                        i = i + 1
              if (dabs(dbl_mb(k_r2+i-1)).gt.r2_thresh) then
c
                 iposh3=int_mb(k_offset+h3b-1)+h3
                 iposh4=int_mb(k_offset+h4b-1)+h4
                 iposp1=int_mb(k_offset+p1b-1)+p1
                 iposp2=int_mb(k_offset+p2b-1)+p2
         if((iposp1.le.iposp2).and.(iposh3.le.iposh4)) then
c
            if(spinsum.ne.4) then ! a<b a<b ----------------------
                xnorm2=xnorm2+dbl_mb(k_r2+i-1)*dbl_mb(k_r2+i-1)
            else  ! a<a a<a --------------------------------------
                xnorm2=xnorm2+dbl_mb(k_r2+i-1)*dbl_mb(k_r2+i-1)
               if(spinsum.eq.4) then  !b<b b<b
                xnorm2=xnorm2+dbl_mb(k_r2+i-1)*dbl_mb(k_r2+i-1)
               end if !b<b b<b
            endif ! a<a a<a
c
             endif !iposp1<iposp2 iposh3<iposh4 (particle int.n.shifted)
              endif !thresh
                      enddo
                    enddo
                  enddo
                enddo
                if (.not.ma_pop_stack(l_r2))
     1            call errquit('tce_print_x2: MA problem',1,MA_ERR)
              endif
              endif
              endif
            enddo
          enddo
        enddo
      enddo
      endif
      call ga_sync()
 9000 format(/,1x,a)
 9010 format(1x,i5,a4,'(',a5,') ',i5,a4,'(',a5,') --- ',
     1          i5,a4,'(',a5,') ',i5,a4,'(',a5,')',f20.10)
      return
      end
c
c
c
