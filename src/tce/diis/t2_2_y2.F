      subroutine t2_2_y2(d_t2,d_y2,k_t2_offset,k_y2_offset)
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
      integer d_y2
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
      integer k_y2_offset
      integer size
      integer l_t2,k_t2
      integer l_y2,k_y2
      integer dp1,dp2,dh3,dh4
      integer i,ipos
      integer nprocs
      integer count
      integer next
      integer iter 
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
     1                ((ioalg.eq.6).and.(.not.fileisga(d_t2))))
      nprocs = ga_nnodes()
      count = 0
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
                 dp1=int_mb(k_range+p1b-1)
                 dp2=int_mb(k_range+p2b-1)
                 dh3=int_mb(k_range+h3b-1)
                 dh4=int_mb(k_range+h4b-1)
c                  size = int_mb(k_range+p1b-1) * int_mb(k_range+p2b-1)
c     1                 * int_mb(k_range+h3b-1) * int_mb(k_range+h4b-1)
                  size = dp1*dp2*dh3*dh4
                  if (.not.ma_push_get(mt_dbl,size,'t2',l_t2,k_t2))
     1              call errquit('tce_jacobi_t2t: MA problem',0,MA_ERR)
                  if (.not.ma_push_get(mt_dbl,size,'y2',l_y2,k_y2))
     1              call errquit('tce_jacobi_t2y: MA problem',0,MA_ERR)
                  do i=1,size
                   dbl_mb(k_t2+i-1)=0.0d0
                   dbl_mb(k_y2+i-1)=0.0d0
                  enddo
c
c      ((((h3b-1)*noab+h4b-1)
c     2              *nvab+p1b-noab-1)*nvab+p2b-noab-1)
c
                  call get_hash_block(d_t2,dbl_mb(k_t2),size,
     1              int_mb(k_t2_offset),((((p1b-noab-1)*nvab+p2b-noab-1)
     2              *noab+h3b-1)*noab+h4b-1))
                  i = 0
                  do p1 = 1,int_mb(k_range+p1b-1)
                    do p2 = 1,int_mb(k_range+p2b-1)
                      do h3 = 1,int_mb(k_range+h3b-1)
                        do h4 = 1,int_mb(k_range+h4b-1)
                          i = i + 1
                          ipos=p2+dp2*(p1-1+dp1*(h4-1+dh4*(h3-1)))
                          dbl_mb(k_y2+ipos-1)=dbl_mb(k_t2+i-1)
                        enddo
                      enddo
                    enddo
                  enddo
                  call put_hash_block(d_y2,dbl_mb(k_y2),size,
     1              int_mb(k_y2_offset),((((h3b-1)*noab+h4b-1)
     2              *nvab+p1b-noab-1)*nvab+p2b-noab-1))
c ----------------------------------------------------------------------------
                  if (.not.ma_pop_stack(l_y2))
     1              call errquit('tce_jacobi_t2ye: MA problem',1,MA_ERR)
                  if (.not.ma_pop_stack(l_t2))
     1              call errquit('tce_jacobi_t2te: MA problem',1,MA_ERR)
                endif
                endif
                endif
      next = NXTASK(nprocs, 1)
              endif
              count = count + 1
            enddo
          enddo
        enddo
      enddo
      next = NXTASK(-nprocs, 1)
      call ga_sync()
      return
      end
c
c
c
c
c
