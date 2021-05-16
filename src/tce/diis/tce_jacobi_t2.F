      subroutine tce_jacobi_t2(d_r2,d_t2,k_t2_offset,shift)
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
      integer d_t2
      integer p1b
      integer p2b
      integer h3b
      integer h4b
      integer p1
      integer p2
      integer h3
      integer h4
      integer k_t2_offset
      integer size
      integer l_r2,k_r2
      integer i
      integer nprocs
      integer count
      integer next
      integer iter 
      INTEGER NXTASK
      EXTERNAL NXTASK
      logical nodezero
      logical noloadbalance
      double precision shift
c
c      if(iter.le.100) then
c       shift=-(2.0d0)*zlshift
c      else
c       shift=0.0d0
c      end if
c
c     ================
c     Loop over blocks
c     ================
c
      nodezero = (ga_nodeid().eq.0)
      noloadbalance = ((ioalg.eq.4).or. 
     1                ((ioalg.eq.6).and.(.not.fileisga(d_r2))))
      nprocs = ga_nnodes()
      count = 0
cc      next = nxtask(nprocs,1)
      next = NXTASK(nprocs, 1)
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
     2            .eq. 0) then
                  size = int_mb(k_range+p1b-1) * int_mb(k_range+p2b-1)
     1                 * int_mb(k_range+h3b-1) * int_mb(k_range+h4b-1)
                  if (.not.ma_push_get(mt_dbl,size,'r2',l_r2,k_r2))
     1              call errquit('tce_jacobi_t2: MA problem',0,MA_ERR)
                  call get_hash_block(d_r2,dbl_mb(k_r2),size,
     1              int_mb(k_t2_offset),((((p1b-noab-1)*nvab+p2b-noab-1)
     2              *noab+h3b-1)*noab+h4b-1))
                  i = 0
                  do p1 = 1,int_mb(k_range+p1b-1)
                    do p2 = 1,int_mb(k_range+p2b-1)
                      do h3 = 1,int_mb(k_range+h3b-1)
                        do h4 = 1,int_mb(k_range+h4b-1)
                          i = i + 1
                          dbl_mb(k_r2+i-1) = dbl_mb(k_r2+i-1)
     1           / (-dbl_mb(k_evl_sorted+int_mb(k_offset+p1b-1)+p1-1)
     2              -dbl_mb(k_evl_sorted+int_mb(k_offset+p2b-1)+p2-1)
     3              +dbl_mb(k_evl_sorted+int_mb(k_offset+h3b-1)+h3-1)
     4              +dbl_mb(k_evl_sorted+int_mb(k_offset+h4b-1)+h4-1)
     5             +shift)
                        enddo
                      enddo
                    enddo
                  enddo
                  call add_hash_block(d_t2,dbl_mb(k_r2),size,
     1              int_mb(k_t2_offset),((((p1b-noab-1)*nvab+p2b-noab-1)
     2              *noab+h3b-1)*noab+h4b-1))
c update of the res.-double vector to the form of increment used in DIIS proc.
ccx                  call put_hash_block(d_r2,dbl_mb(k_r2),size,
ccx     1              int_mb(k_t2_offset),((((p1b-noab-1)*nvab+p2b-noab-1)
ccx     2              *noab+h3b-1)*noab+h4b-1))
c ----------------------------------------------------------------------------
                  if (nodezero.and.util_print('t2',print_debug)) 
     1              then
                    call get_hash_block(d_t2,dbl_mb(k_r2),size,
     1              int_mb(k_t2_offset),((((p1b-noab-1)*nvab+p2b-noab-1)
     2              *noab+h3b-1)*noab+h4b-1))
                    call ma_print_compact
     2              (dbl_mb(k_r2),size,1,'t2')
                  endif
                  if (.not.ma_pop_stack(l_r2))
     1              call errquit('tce_jacobi_t2: MA problem',1,MA_ERR)
                endif
                endif
                endif
cc                next = nxtask(nprocs,1)
      next = NXTASK(nprocs, 1)
              endif
              count = count + 1
            enddo
          enddo
        enddo
      enddo
cc      next = nxtask(-nprocs,1)
      next = NXTASK(-nprocs, 1)
      call ga_sync()
      return
      end
c
c
c
c
c
c
      subroutine tce_r2_divide(d_r2,k_t2_offset)
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
      integer k_t2_offset
      integer size
      integer l_r2,k_r2
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
c     ================
c     Loop over blocks
c     ================
c
      nodezero = (ga_nodeid().eq.0)
      noloadbalance = ((ioalg.eq.4).or. 
     1                ((ioalg.eq.6).and.(.not.fileisga(d_r2))))
      nprocs = ga_nnodes()
      count = 0
cc      next = nxtask(nprocs,1)
      next = NXTASK(nprocs, 1)
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
     2            .eq. 0) then
                  size = int_mb(k_range+p1b-1) * int_mb(k_range+p2b-1)
     1                 * int_mb(k_range+h3b-1) * int_mb(k_range+h4b-1)
                  if (.not.ma_push_get(mt_dbl,size,'r2',l_r2,k_r2))
     1              call errquit('tce_r2_divide: MA problem',0,MA_ERR)
                  call get_hash_block(d_r2,dbl_mb(k_r2),size,
     1              int_mb(k_t2_offset),((((p1b-noab-1)*nvab+p2b-noab-1)
     2              *noab+h3b-1)*noab+h4b-1))
                  i = 0
                  do p1 = 1,int_mb(k_range+p1b-1)
                    do p2 = 1,int_mb(k_range+p2b-1)
                      do h3 = 1,int_mb(k_range+h3b-1)
                        do h4 = 1,int_mb(k_range+h4b-1)
                          i = i + 1
                          dbl_mb(k_r2+i-1) = dbl_mb(k_r2+i-1)
     1           / (-dbl_mb(k_evl_sorted+int_mb(k_offset+p1b-1)+p1-1)
     2              -dbl_mb(k_evl_sorted+int_mb(k_offset+p2b-1)+p2-1)
     3              +dbl_mb(k_evl_sorted+int_mb(k_offset+h3b-1)+h3-1)
     4              +dbl_mb(k_evl_sorted+int_mb(k_offset+h4b-1)+h4-1))
                        enddo
                      enddo
                    enddo
                  enddo
c update of the res.-double vector to the form of increment used in DIIS proc.
                  call put_hash_block(d_r2,dbl_mb(k_r2),size,
     1              int_mb(k_t2_offset),((((p1b-noab-1)*nvab+p2b-noab-1)
     2              *noab+h3b-1)*noab+h4b-1))
c ----------------------------------------------------------------------------
                  if (.not.ma_pop_stack(l_r2))
     1              call errquit('tce_r2_divide: MA problem',1,MA_ERR)
                endif
                endif
                endif
cc                next = nxtask(nprocs,1)
      next = NXTASK(nprocs, 1)
              endif
              count = count + 1
            enddo
          enddo
        enddo
      enddo
cc      next = nxtask(-nprocs,1)
      next = NXTASK(-nprocs, 1)
      call ga_sync()
      return
      end
c
c
c
c
c
      subroutine tce_jacobi_t2a(d_r2,d_t2,k_t2_offset,iter)
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
      integer d_t2
      integer p1b
      integer p2b
      integer h3b
      integer h4b
      integer p1
      integer p2
      integer h3
      integer h4
      integer k_t2_offset
      integer size
      integer l_r2,k_r2
      integer i
      integer nprocs
      integer count
      integer next
      LOGICAL is_active_1,is_active_2,is_active_3,is_active_4
cc      integer nxtask
c
      integer iter 
c
cc      external nxtask
      INTEGER NXTASK
      EXTERNAL NXTASK
      logical nodezero
      logical noloadbalance
c --- shift ---
      double precision shift
c -------------
      if(iter.le.100) then
       shift=-(2.0d0)*zlshift
      else
       shift=0.0d0
      end if
c
c     ================
c     Loop over blocks
c     ================
c
      nodezero = (ga_nodeid().eq.0)
      noloadbalance = ((ioalg.eq.4).or. 
     1                ((ioalg.eq.6).and.(.not.fileisga(d_r2))))
      nprocs = ga_nnodes()
      count = 0
cc      next = nxtask(nprocs,1)
      next = NXTASK(nprocs, 1)
      do p1b = noab+1,noab+nvab
        do p2b = p1b,noab+nvab
          do h3b = 1,noab
            do h4b = h3b,noab
              if (noloadbalance.or.(next.eq.count)) then
                IF(is_active_4(p1b,p2b,h3b,h4b)) THEN
                if (int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1) 
     1            .eq. int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1)) then
                if ((.not.restricted).or.
     1            (int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1)+
     2            int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1).ne.8)) then
                if (ieor(int_mb(k_sym+p1b-1),ieor(int_mb(k_sym+p2b-1),
     1            ieor(int_mb(k_sym+h3b-1),int_mb(k_sym+h4b-1))))
     2            .eq. 0) then
                  size = int_mb(k_range+p1b-1) * int_mb(k_range+p2b-1)
     1                 * int_mb(k_range+h3b-1) * int_mb(k_range+h4b-1)
                  if (.not.ma_push_get(mt_dbl,size,'r2',l_r2,k_r2))
     1              call errquit('tce_jacobi_t2: MA problem',0,MA_ERR)
                  call get_hash_block(d_r2,dbl_mb(k_r2),size,
     1              int_mb(k_t2_offset),((((p1b-noab-1)*nvab+p2b-noab-1)
     2              *noab+h3b-1)*noab+h4b-1))
                  i = 0
                  do p1 = 1,int_mb(k_range+p1b-1)
                    do p2 = 1,int_mb(k_range+p2b-1)
                      do h3 = 1,int_mb(k_range+h3b-1)
                        do h4 = 1,int_mb(k_range+h4b-1)
                          i = i + 1
                          dbl_mb(k_r2+i-1) = dbl_mb(k_r2+i-1)
     1           / (-dbl_mb(k_evl_sorted+int_mb(k_offset+p1b-1)+p1-1)
     2              -dbl_mb(k_evl_sorted+int_mb(k_offset+p2b-1)+p2-1)
     3              +dbl_mb(k_evl_sorted+int_mb(k_offset+h3b-1)+h3-1)
     4              +dbl_mb(k_evl_sorted+int_mb(k_offset+h4b-1)+h4-1)
     5              +shift)
                        enddo
                      enddo
                    enddo
                  enddo
                  call add_hash_block(d_t2,dbl_mb(k_r2),size,
     1              int_mb(k_t2_offset),((((p1b-noab-1)*nvab+p2b-noab-1)
     2              *noab+h3b-1)*noab+h4b-1))
c update of the res.-double vector to the form of increment used in DIIS proc.
ccx                  call put_hash_block(d_r2,dbl_mb(k_r2),size,
ccx     1              int_mb(k_t2_offset),((((p1b-noab-1)*nvab+p2b-noab-1)
ccx     2              *noab+h3b-1)*noab+h4b-1))
c ----------------------------------------------------------------------------
                  if (nodezero.and.util_print('t2',print_debug)) 
     1              then
                    call get_hash_block(d_t2,dbl_mb(k_r2),size,
     1              int_mb(k_t2_offset),((((p1b-noab-1)*nvab+p2b-noab-1)
     2              *noab+h3b-1)*noab+h4b-1))
                    call ma_print_compact
     2              (dbl_mb(k_r2),size,1,'t2')
                  endif
                  if (.not.ma_pop_stack(l_r2))
     1              call errquit('tce_jacobi_t2: MA problem',1,MA_ERR)
                endif
                endif
                endif
                endif
cc                next = nxtask(nprocs,1)
      next = NXTASK(nprocs, 1)
              endif
              count = count + 1
            enddo
          enddo
        enddo
      enddo
cc      next = nxtask(-nprocs,1)
      next = NXTASK(-nprocs, 1)
      call ga_sync()
      return
      end
c

