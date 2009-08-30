      subroutine tce_print_x2(d_r2,k_r2_offset,r2_thresh,irrep)
c
c $Id: tce_print_x2.F,v 1.11 2009-02-25 00:46:45 d3p975 Exp $
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
      character*4 irrep_p1,irrep_p2,irrep_h3,irrep_h4
      character*5 spinnames(2),spin_p1,spin_p2,spin_h3,spin_h4
      data spinnames/'alpha','beta '/
c
      nodezero = (ga_nodeid().eq.0)
      if (nodezero) write(LuOut,9000) 'Doubles contributions'
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
      call ga_sync()
 9000 format(/,1x,a)
 9010 format(1x,i5,a4,'(',a5,') ',i5,a4,'(',a5,') --- ',
     1          i5,a4,'(',a5,') ',i5,a4,'(',a5,')',f20.10)
      return
      end








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
      if (nodezero) write(LuOut,9000) 'Doubles contributions'
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
