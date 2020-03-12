      subroutine tce_loc_j1(d_r1,k_r1_offset,d_t1,k_t1_offset)
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
      integer nxtask
c --- t2 scan ---
      integer d_t1
      integer k_t1_offset
      integer l_t1,k_t1
      integer pp1b
      integer hh2b
      integer tsize
      integer ii
      integer pp1,hh2
      integer ipa1,iha2
      integer ip1,ih2
c ---------------
      external nxtask
      logical nodezero
      logical noloadbalance
c
      nodezero = (ga_nodeid().eq.0)
      noloadbalance = ((ioalg.eq.4).or. 
     1                ((ioalg.eq.6).and.(.not.fileisga(d_r1))))
      nprocs = ga_nnodes()
      count = 0
      next = nxtask(nprocs,1)
      do p1b = noab+1,noab+nvab
        do h2b = 1,noab
          if (noloadbalance.or.(next.eq.count)) then
            if (int_mb(k_spin+p1b-1) .eq. int_mb(k_spin+h2b-1)) then 
            if ((.not.restricted).or.(int_mb(k_spin+p1b-1)
     1        +int_mb(k_spin+h2b-1).ne.4)) then
            if (ieor(int_mb(k_sym+p1b-1),int_mb(k_sym+h2b-1))
     1        .eq.irrep_x) then
              size = int_mb(k_range+p1b-1) * int_mb(k_range+h2b-1)
              if (.not.ma_push_get(mt_dbl,size,'r1',l_r1,k_r1))
     1          call errquit('tce_jacobi_x1: MA problem',0,MA_ERR)
              call get_hash_block(d_r1,dbl_mb(k_r1),size,
     1          int_mb(k_r1_offset),((p1b-noab-1)*noab+h2b-1))
              i = 0
              do p1 = 1,int_mb(k_range+p1b-1)
                do h2 = 1,int_mb(k_range+h2b-1)
                  i = i + 1
                    ip1=int_mb(k_offset+p1b-1)+p1
                    ih2=int_mb(k_offset+h2b-1)+h2
c ------ t1 loop --------------------------
      do pp1b = noab+1,noab+nvab
        do hh2b = 1,noab
          if (int_mb(k_spin+pp1b-1) .eq. int_mb(k_spin+hh2b-1)) then
          if ((.not.restricted).or.(int_mb(k_spin+pp1b-1)
     1        +int_mb(k_spin+hh2b-1).ne.4)) then
          if (ieor(int_mb(k_sym+pp1b-1),int_mb(k_sym+hh2b-1)).eq.0) then
            tsize = int_mb(k_range+pp1b-1) * int_mb(k_range+hh2b-1)
              if (.not.ma_push_get(mt_dbl,tsize,'t1',l_t1,k_t1))
     1          call errquit('tce_jacobi_t1: MA problem',0,MA_ERR)
              call get_hash_block(d_t1,dbl_mb(k_t1),tsize,
     1          int_mb(k_t1_offset),((pp1b-noab-1)*noab+hh2b-1))
              ii = 0
              do pp1 = 1,int_mb(k_range+pp1b-1)
                do hh2 = 1,int_mb(k_range+hh2b-1)
                  ii = ii + 1
                   ipa1=int_mb(k_offset+pp1b-1)+pp1
                   iha2=int_mb(k_offset+hh2b-1)+hh2
                   if ((ip1.eq.ipa1).or.(ih2.eq.iha2)) then
                     dbl_mb(k_r1+i-1)=dbl_mb(k_r1+i-1)+
     1               dbl_mb(k_t1+ii-1)*dbl_mb(k_t1+ii-1)
                   end if
                enddo
              enddo
              if (.not.ma_pop_stack(l_t1))
     1          call errquit('tce_jacobi_t1: MA problem',1,MA_ERR)
            endif
            endif
            endif
        enddo
      enddo
c ------ end of the t1 loop ---------------
                enddo
              enddo
              call put_hash_block(d_r1,dbl_mb(k_r1),size,
     1          int_mb(k_r1_offset),((p1b-noab-1)*noab+h2b-1))
              if (.not.ma_pop_stack(l_r1))
     1          call errquit('tce_jacobi_x1: MA problem',1,MA_ERR)
            endif
            endif
            endif
            next = nxtask(nprocs,1)
          endif
          count = count + 1
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
c
c
      subroutine tce_loc_j1_one(d_r1,k_r1_offset)
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
      integer nxtask
      external nxtask
      logical nodezero
      logical noloadbalance
c
      nodezero = (ga_nodeid().eq.0)
      noloadbalance = ((ioalg.eq.4).or. 
     1                ((ioalg.eq.6).and.(.not.fileisga(d_r1))))
      nprocs = ga_nnodes()
      count = 0
      next = nxtask(nprocs,1)
      do p1b = noab+1,noab+nvab
        do h2b = 1,noab
          if (noloadbalance.or.(next.eq.count)) then
            if (int_mb(k_spin+p1b-1) .eq. int_mb(k_spin+h2b-1)) then 
            if ((.not.restricted).or.(int_mb(k_spin+p1b-1)
     1        +int_mb(k_spin+h2b-1).ne.4)) then
            if (ieor(int_mb(k_sym+p1b-1),int_mb(k_sym+h2b-1))
     1        .eq.irrep_x) then
              size = int_mb(k_range+p1b-1) * int_mb(k_range+h2b-1)
              if (.not.ma_push_get(mt_dbl,size,'r1',l_r1,k_r1))
     1          call errquit('tce_jacobi_x1: MA problem',0,MA_ERR)
              call get_hash_block(d_r1,dbl_mb(k_r1),size,
     1          int_mb(k_r1_offset),((p1b-noab-1)*noab+h2b-1))
              i = 0
              do p1 = 1,int_mb(k_range+p1b-1)
                do h2 = 1,int_mb(k_range+h2b-1)
                  i = i + 1
                  dbl_mb(k_r1+i-1) = 1.0d0
                enddo
              enddo
              call put_hash_block(d_r1,dbl_mb(k_r1),size,
     1          int_mb(k_r1_offset),((p1b-noab-1)*noab+h2b-1))
              if (.not.ma_pop_stack(l_r1))
     1          call errquit('tce_jacobi_x1: MA problem',1,MA_ERR)
            endif
            endif
            endif
            next = nxtask(nprocs,1)
          endif
          count = count + 1
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
      subroutine tce_loc_j1_inv(d_r1,k_r1_offset)
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
      integer nxtask
      external nxtask
      logical nodezero
      logical noloadbalance
c 
      nodezero = (ga_nodeid().eq.0)
      noloadbalance = ((ioalg.eq.4).or. 
     1                ((ioalg.eq.6).and.(.not.fileisga(d_r1))))
      nprocs = ga_nnodes()
      count = 0
      next = nxtask(nprocs,1)
      do p1b = noab+1,noab+nvab
        do h2b = 1,noab
          if (noloadbalance.or.(next.eq.count)) then
            if (int_mb(k_spin+p1b-1) .eq. int_mb(k_spin+h2b-1)) then 
            if ((.not.restricted).or.(int_mb(k_spin+p1b-1)
     1        +int_mb(k_spin+h2b-1).ne.4)) then
            if (ieor(int_mb(k_sym+p1b-1),int_mb(k_sym+h2b-1))
     1        .eq.irrep_x) then
              size = int_mb(k_range+p1b-1) * int_mb(k_range+h2b-1)
              if (.not.ma_push_get(mt_dbl,size,'r1',l_r1,k_r1))
     1          call errquit('tce_jacobi_x1: MA problem',0,MA_ERR)
              call get_hash_block(d_r1,dbl_mb(k_r1),size,
     1          int_mb(k_r1_offset),((p1b-noab-1)*noab+h2b-1))
              i = 0
              do p1 = 1,int_mb(k_range+p1b-1)
                do h2 = 1,int_mb(k_range+h2b-1)
                  i = i + 1
                  dbl_mb(k_r1+i-1) = 1.0d0/dbl_mb(k_r1+i-1)
                enddo
              enddo
              call put_hash_block(d_r1,dbl_mb(k_r1),size,
     1          int_mb(k_r1_offset),((p1b-noab-1)*noab+h2b-1))
              if (.not.ma_pop_stack(l_r1))
     1          call errquit('tce_jacobi_x1: MA problem',1,MA_ERR)
            endif
            endif
            endif
            next = nxtask(nprocs,1)
          endif
          count = count + 1
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
c
      subroutine tce_j1_x1_scaling(d_r1,k_r1_offset,d_j1,k_j1_offset)
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
      integer d_r1
      integer p1b
      integer h2b
      integer p1
      integer h2
      integer k_r1_offset
      integer size
      integer l_r1,k_r1
      integer i
c ---
      integer d_j1
      integer l_j1,k_j1
      integer k_j1_offset
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
     1                ((ioalg.eq.6).and.(.not.fileisga(d_r1))))
      nprocs = ga_nnodes()
      count = 0
      next = nxtask(nprocs,1)
      do p1b = noab+1,noab+nvab
        do h2b = 1,noab
          if (noloadbalance.or.(next.eq.count)) then
            if (int_mb(k_spin+p1b-1) .eq. int_mb(k_spin+h2b-1)) then 
            if ((.not.restricted).or.(int_mb(k_spin+p1b-1)
     1        +int_mb(k_spin+h2b-1).ne.4)) then
            if (ieor(int_mb(k_sym+p1b-1),int_mb(k_sym+h2b-1))
     1        .eq.irrep_x) then
              size = int_mb(k_range+p1b-1) * int_mb(k_range+h2b-1)
              if (.not.ma_push_get(mt_dbl,size,'r1',l_r1,k_r1))
     1          call errquit('tce_jacobi_x1: MA problem',0,MA_ERR)
              if (.not.ma_push_get(mt_dbl,size,'r1',l_j1,k_j1))
     1          call errquit('tce_jacobi_x1: MA problem',0,MA_ERR)
              call get_hash_block(d_r1,dbl_mb(k_r1),size,
     1          int_mb(k_r1_offset),((p1b-noab-1)*noab+h2b-1))
              call get_hash_block(d_j1,dbl_mb(k_j1),size,
     1          int_mb(k_j1_offset),((p1b-noab-1)*noab+h2b-1))
              i = 0
              do p1 = 1,int_mb(k_range+p1b-1)
                do h2 = 1,int_mb(k_range+h2b-1)
                  i = i + 1
             dbl_mb(k_r1+i-1) = dbl_mb(k_r1+i-1)*dbl_mb(k_j1+i-1)
                enddo
              enddo
              call put_hash_block(d_r1,dbl_mb(k_r1),size,
     1          int_mb(k_r1_offset),((p1b-noab-1)*noab+h2b-1))
              if (.not.ma_pop_stack(l_j1))
     1          call errquit('tce_jacobi_x1: MA problem',1,MA_ERR)
              if (.not.ma_pop_stack(l_r1))
     1          call errquit('tce_jacobi_x1: MA problem',1,MA_ERR)
            endif
            endif
            endif
            next = nxtask(nprocs,1)
          endif
          count = count + 1
        enddo
      enddo
      next = nxtask(-nprocs,1)
      call ga_sync()
      return
      end
