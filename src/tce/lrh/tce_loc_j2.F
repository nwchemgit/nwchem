      subroutine tce_loc_j2(d_r2,k_r2_offset,d_t2,k_t2_offset)
c
c $Id$
c
c TO DO - initial 1.0-ing of the d_t2 file
c
c
c
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
      integer nprocs
      integer count
      integer next
      integer nxtask
c --- t2 scan ---
      integer d_t2
      integer k_t2_offset
      integer l_t2,k_t2
      integer pp1b,pp2b
      integer hh3b,hh4b
      integer tsize
      integer ii
      integer pp1,pp2,hh3,hh4
      integer ipa1,ipa2,iha3,iha4
      integer ip1,ip2,ih3,ih4
      integer spin_sum,spin_suml
      integer p2alpha
      integer h4alpha
      integer ip2alpha
      integer ih4alpha
c ---------------
      external nxtask
      logical nodezero
      logical noloadbalance
c
      nodezero = (ga_nodeid().eq.0)
      noloadbalance = ((ioalg.eq.4).or. 
     1                ((ioalg.eq.6).and.(.not.fileisga(d_r2))))
      nprocs = ga_nnodes()
      count = 0
      next = nxtask(nprocs,1)
      do p1b = noab+1,noab+nvab
        do p2b = p1b,noab+nvab
          do h3b = 1,noab
            do h4b = h3b,noab
              if (noloadbalance.or.(next.eq.count)) then
                if (int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1) 
     1            .eq. int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1)) then
                if ((.not.restricted).or.
     1            (int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1)+
     2            int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1).ne.8)) then
                if (ieor(int_mb(k_sym+p1b-1),ieor(int_mb(k_sym+p2b-1),
     1            ieor(int_mb(k_sym+h3b-1),int_mb(k_sym+h4b-1))))
     2            .eq. irrep_x) then
                  size = int_mb(k_range+p1b-1) * int_mb(k_range+p2b-1)
     1                 * int_mb(k_range+h3b-1) * int_mb(k_range+h4b-1)
                  if (.not.ma_push_get(mt_dbl,size,'r2',l_r2,k_r2))
     1              call errquit('tce_jacobi_x2: MA problem',0,MA_ERR)
                  call get_hash_block(d_r2,dbl_mb(k_r2),size,
     1              int_mb(k_r2_offset),((((p1b-noab-1)*nvab+p2b-noab-1)
     2              *noab+h3b-1)*noab+h4b-1))
c
                  if(restricted) then
                   spin_sum=int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1)+
     1                      int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1)
                   if(spin_sum.eq.6) then
                    p2alpha=int_mb(k_alpha+p2b-1)
                    h4alpha=int_mb(k_alpha+h4b-1) 
                   end if
                  end if 
                  i = 0
                  do p1 = 1,int_mb(k_range+p1b-1)
                    do p2 = 1,int_mb(k_range+p2b-1)
                      do h3 = 1,int_mb(k_range+h3b-1)
                        do h4 = 1,int_mb(k_range+h4b-1)
                          i = i + 1
                         ip1=int_mb(k_offset+p1b-1)+p1
                         ip2=int_mb(k_offset+p2b-1)+p2
                         ih3=int_mb(k_offset+h3b-1)+h3
                         ih4=int_mb(k_offset+h4b-1)+h4
                         if(restricted.and.(spin_sum.eq.6)) then
                          ip2alpha=int_mb(k_offset+p2alpha-1)+p2
                          ih4alpha=int_mb(k_offset+h4alpha-1)+h4 
                         end if
c                          dbl_mb(k_r2+i-1) = dbl_mb(k_r2+i-1)
c     1           / (-dbl_mb(k_evl_sorted+int_mb(k_offset+p1b-1)+p1-1)
c     2              -dbl_mb(k_evl_sorted+int_mb(k_offset+p2b-1)+p2-1)
c     3              +dbl_mb(k_evl_sorted+int_mb(k_offset+h3b-1)+h3-1)
c     4              +dbl_mb(k_evl_sorted+int_mb(k_offset+h4b-1)+h4-1))
c --- loop over t2 -----------
      do pp1b = noab+1,noab+nvab
        do pp2b = pp1b,noab+nvab
          do hh3b = 1,noab
            do hh4b = hh3b,noab
                if (int_mb(k_spin+pp1b-1)+int_mb(k_spin+pp2b-1)
     1            .eq. int_mb(k_spin+hh3b-1)+int_mb(k_spin+hh4b-1)) then
                if ((.not.restricted).or.
     1           (int_mb(k_spin+pp1b-1)+int_mb(k_spin+pp2b-1)+
     2           int_mb(k_spin+hh3b-1)+int_mb(k_spin+hh4b-1).ne.8)) then
                if (ieor(int_mb(k_sym+pp1b-1),ieor(int_mb(k_sym+pp2b-1),
     1            ieor(int_mb(k_sym+hh3b-1),int_mb(k_sym+hh4b-1))))
     2            .eq. 0) then
                 tsize = int_mb(k_range+pp1b-1) * int_mb(k_range+pp2b-1)
     1                 * int_mb(k_range+hh3b-1) * int_mb(k_range+hh4b-1)
                  if (.not.ma_push_get(mt_dbl,tsize,'t2',l_t2,k_t2))
     1              call errquit('tce_jacobi_t2: MA problem',0,MA_ERR)
                  call get_hash_block(d_t2,dbl_mb(k_t2),tsize,
     1              int_mb(k_t2_offset),
     2              ((((pp1b-noab-1)*nvab+pp2b-noab-1)
     2              *noab+hh3b-1)*noab+hh4b-1))
                  if(restricted) then
                  spin_suml=int_mb(k_spin+pp1b-1)+int_mb(k_spin+pp2b-1)+
     1                      int_mb(k_spin+hh3b-1)+int_mb(k_spin+hh4b-1)
                  end if
                  ii = 0
                  do pp1 = 1,int_mb(k_range+pp1b-1)
                    do pp2 = 1,int_mb(k_range+pp2b-1)
                      do hh3 = 1,int_mb(k_range+hh3b-1)
                        do hh4 = 1,int_mb(k_range+hh4b-1)
                          ii = ii + 1
c
                         ipa1=int_mb(k_offset+pp1b-1)+pp1
                         ipa2=int_mb(k_offset+pp2b-1)+pp2
                         iha3=int_mb(k_offset+hh3b-1)+hh3
                         iha4=int_mb(k_offset+hh4b-1)+hh4
c
              if  (ipa1.lt.ipa2.AND.iha3.lt.iha4) THEN
               if ((ip1.eq.ipa1).or.(ip1.eq.ipa2).or.
     1            (ip2.eq.ipa1).or.(ip2.eq.ipa2).or.
     2            (ih3.eq.iha3).or.(ih3.eq.iha4).or.
     3            (ih4.eq.iha3).or.(ih4.eq.iha4)) then
                   dbl_mb(k_r2+i-1)=dbl_mb(k_r2+i-1)+
     1                              dbl_mb(k_t2+ii-1)*dbl_mb(k_t2+ii-1)
               end if 
               if (restricted.and.(spin_sum.eq.6).and.(spin_suml.eq.4)) 
     1            then
                if ((ip2alpha.eq.ipa1).or.(ip2alpha.eq.ipa2).or.
     1             (ih4alpha.eq.iha3).or.(ih4alpha.eq.iha4)) THEN
                   dbl_mb(k_r2+i-1)=dbl_mb(k_r2+i-1)+
     1                              dbl_mb(k_t2+ii-1)*dbl_mb(k_t2+ii-1)
                end if
               end if
              end if
c
c                          dbl_mb(k_t2+ii-1) = dbl_mb(k_t2+ii-1)
c     1           / (-dbl_mb(k_evl_sorted+int_mb(k_offset+pp1b-1)+pp1-1)
c     2              -dbl_mb(k_evl_sorted+int_mb(k_offset+pp2b-1)+pp2-1)
c     3              +dbl_mb(k_evl_sorted+int_mb(k_offset+hh3b-1)+hh3-1)
c     4              +dbl_mb(k_evl_sorted+int_mb(k_offset+hh4b-1)+hh4-1))
                        enddo
                      enddo
                    enddo
                  enddo
                  if (.not.ma_pop_stack(l_t2))
     1              call errquit('tce_jacobi_t2: MA problem',1,MA_ERR)
                endif
                endif
                endif
            enddo
          enddo
        enddo
      enddo
c --- end of t2 loop -----

                        enddo
                      enddo
                    enddo
                  enddo
                  call put_hash_block(d_r2,dbl_mb(k_r2),size,
     1              int_mb(k_r2_offset),((((p1b-noab-1)*nvab+p2b-noab-1)
     2              *noab+h3b-1)*noab+h4b-1))
                  if (.not.ma_pop_stack(l_r2))
     1              call errquit('tce_jacobi_x2: MA problem',1,MA_ERR)
                endif
                endif
                endif
                next = nxtask(nprocs,1)
              endif
              count = count + 1
            enddo
          enddo
        enddo
      enddo
      next = nxtask(-nprocs,1)
      call ga_sync()
      return
      end
c
c
c
c
c
      subroutine tce_loc_j2_one(d_r2,k_r2_offset)
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
      integer nprocs
      integer count
      integer next
      integer nxtask
      external nxtask
      logical nodezero
      logical noloadbalance
c
      nodezero = (ga_nodeid().eq.0)
      noloadbalance = ((ioalg.eq.4).or. 
     1                ((ioalg.eq.6).and.(.not.fileisga(d_r2))))
      nprocs = ga_nnodes()
      count = 0
      next = nxtask(nprocs,1)
      do p1b = noab+1,noab+nvab
        do p2b = p1b,noab+nvab
          do h3b = 1,noab
            do h4b = h3b,noab
              if (noloadbalance.or.(next.eq.count)) then
                if (int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1) 
     1            .eq. int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1)) then
                if ((.not.restricted).or.
     1            (int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1)+
     2            int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1).ne.8)) then
                if (ieor(int_mb(k_sym+p1b-1),ieor(int_mb(k_sym+p2b-1),
     1            ieor(int_mb(k_sym+h3b-1),int_mb(k_sym+h4b-1))))
     2            .eq. irrep_x) then
                  size = int_mb(k_range+p1b-1) * int_mb(k_range+p2b-1)
     1                 * int_mb(k_range+h3b-1) * int_mb(k_range+h4b-1)
                  if (.not.ma_push_get(mt_dbl,size,'r2',l_r2,k_r2))
     1              call errquit('tce_jacobi_x2: MA problem',0,MA_ERR)
                  call get_hash_block(d_r2,dbl_mb(k_r2),size,
     1              int_mb(k_r2_offset),((((p1b-noab-1)*nvab+p2b-noab-1)
     2              *noab+h3b-1)*noab+h4b-1))
                  i = 0
                  do p1 = 1,int_mb(k_range+p1b-1)
                    do p2 = 1,int_mb(k_range+p2b-1)
                      do h3 = 1,int_mb(k_range+h3b-1)
                        do h4 = 1,int_mb(k_range+h4b-1)
                          i = i + 1
                          dbl_mb(k_r2+i-1) = 1.0d0
                        enddo
                      enddo
                    enddo
                  enddo
                  call put_hash_block(d_r2,dbl_mb(k_r2),size,
     1              int_mb(k_r2_offset),((((p1b-noab-1)*nvab+p2b-noab-1)
     2              *noab+h3b-1)*noab+h4b-1))
                  if (.not.ma_pop_stack(l_r2))
     1              call errquit('tce_jacobi_x2: MA problem',1,MA_ERR)
                endif
                endif
                endif
                next = nxtask(nprocs,1)
              endif
              count = count + 1
            enddo
          enddo
        enddo
      enddo
      next = nxtask(-nprocs,1)
      call ga_sync()
      return
      end
c
c
c
c
c
      subroutine tce_loc_j2_inv(d_r2,k_r2_offset)
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
      integer nprocs
      integer count
      integer next
      integer nxtask
      external nxtask
      logical nodezero
      logical noloadbalance
c
      nodezero = (ga_nodeid().eq.0)
      noloadbalance = ((ioalg.eq.4).or. 
     1                ((ioalg.eq.6).and.(.not.fileisga(d_r2))))
      nprocs = ga_nnodes()
      count = 0
      next = nxtask(nprocs,1)
      do p1b = noab+1,noab+nvab
        do p2b = p1b,noab+nvab
          do h3b = 1,noab
            do h4b = h3b,noab
              if (noloadbalance.or.(next.eq.count)) then
                if (int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1) 
     1            .eq. int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1)) then
                if ((.not.restricted).or.
     1            (int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1)+
     2            int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1).ne.8)) then
                if (ieor(int_mb(k_sym+p1b-1),ieor(int_mb(k_sym+p2b-1),
     1            ieor(int_mb(k_sym+h3b-1),int_mb(k_sym+h4b-1))))
     2            .eq. irrep_x) then
                  size = int_mb(k_range+p1b-1) * int_mb(k_range+p2b-1)
     1                 * int_mb(k_range+h3b-1) * int_mb(k_range+h4b-1)
                  if (.not.ma_push_get(mt_dbl,size,'r2',l_r2,k_r2))
     1              call errquit('tce_jacobi_x2: MA problem',0,MA_ERR)
                  call get_hash_block(d_r2,dbl_mb(k_r2),size,
     1              int_mb(k_r2_offset),((((p1b-noab-1)*nvab+p2b-noab-1)
     2              *noab+h3b-1)*noab+h4b-1))
                  i = 0
                  do p1 = 1,int_mb(k_range+p1b-1)
                    do p2 = 1,int_mb(k_range+p2b-1)
                      do h3 = 1,int_mb(k_range+h3b-1)
                        do h4 = 1,int_mb(k_range+h4b-1)
                          i = i + 1
                          dbl_mb(k_r2+i-1) = 1.0d0/dbl_mb(k_r2+i-1)
                        enddo
                      enddo
                    enddo
                  enddo
                  call put_hash_block(d_r2,dbl_mb(k_r2),size,
     1              int_mb(k_r2_offset),((((p1b-noab-1)*nvab+p2b-noab-1)
     2              *noab+h3b-1)*noab+h4b-1))
                  if (.not.ma_pop_stack(l_r2))
     1              call errquit('tce_jacobi_x2: MA problem',1,MA_ERR)
                endif
                endif
                endif
                next = nxtask(nprocs,1)
              endif
              count = count + 1
            enddo
          enddo
        enddo
      enddo
      next = nxtask(-nprocs,1)
      call ga_sync()
      return
      end
c
c
c
c
c
c
      subroutine tce_j2_x2_scaling(d_r2,k_r2_offset,d_j2,k_j2_offset)
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
c ---
      integer d_j2
      integer k_j2_offset
      integer l_j2,k_j2
c ---
      integer nprocs
      integer count
      integer next
      integer nxtask
      external nxtask
      logical nodezero
      logical noloadbalance
c
      nodezero = (ga_nodeid().eq.0)
      noloadbalance = ((ioalg.eq.4).or. 
     1                ((ioalg.eq.6).and.(.not.fileisga(d_r2))))
      nprocs = ga_nnodes()
      count = 0
      next = nxtask(nprocs,1)
      do p1b = noab+1,noab+nvab
        do p2b = p1b,noab+nvab
          do h3b = 1,noab
            do h4b = h3b,noab
              if (noloadbalance.or.(next.eq.count)) then
                if (int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1) 
     1            .eq. int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1)) then
                if ((.not.restricted).or.
     1            (int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1)+
     2            int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1).ne.8)) then
                if (ieor(int_mb(k_sym+p1b-1),ieor(int_mb(k_sym+p2b-1),
     1            ieor(int_mb(k_sym+h3b-1),int_mb(k_sym+h4b-1))))
     2            .eq. irrep_x) then
                  size = int_mb(k_range+p1b-1) * int_mb(k_range+p2b-1)
     1                 * int_mb(k_range+h3b-1) * int_mb(k_range+h4b-1)
                  if (.not.ma_push_get(mt_dbl,size,'r2',l_r2,k_r2))
     1              call errquit('tce_jacobi_x2: MA problem',0,MA_ERR)
                  if (.not.ma_push_get(mt_dbl,size,'r2',l_j2,k_j2))
     1              call errquit('tce_jacobi_x2: MA problem',0,MA_ERR)
                  call get_hash_block(d_r2,dbl_mb(k_r2),size,
     1              int_mb(k_r2_offset),((((p1b-noab-1)*nvab+p2b-noab-1)
     2              *noab+h3b-1)*noab+h4b-1))
                  call get_hash_block(d_j2,dbl_mb(k_j2),size,
     1              int_mb(k_j2_offset),((((p1b-noab-1)*nvab+p2b-noab-1)
     2              *noab+h3b-1)*noab+h4b-1))
                  i = 0
                  do p1 = 1,int_mb(k_range+p1b-1)
                    do p2 = 1,int_mb(k_range+p2b-1)
                      do h3 = 1,int_mb(k_range+h3b-1)
                        do h4 = 1,int_mb(k_range+h4b-1)
                          i = i + 1
                    dbl_mb(k_r2+i-1) = dbl_mb(k_r2+i-1)*dbl_mb(k_j2+i-1)
                        enddo
                      enddo
                    enddo
                  enddo
                  call put_hash_block(d_r2,dbl_mb(k_r2),size,
     1              int_mb(k_r2_offset),((((p1b-noab-1)*nvab+p2b-noab-1)
     2              *noab+h3b-1)*noab+h4b-1))
                  if (.not.ma_pop_stack(l_j2))
     1              call errquit('tce_jacobi_x2: MA problem',1,MA_ERR)
                  if (.not.ma_pop_stack(l_r2))
     1              call errquit('tce_jacobi_x2: MA problem',1,MA_ERR)
                endif
                endif
                endif
                next = nxtask(nprocs,1)
              endif
              count = count + 1
            enddo
          enddo
        enddo
      enddo
      next = nxtask(-nprocs,1)
      call ga_sync()
      return
      end
