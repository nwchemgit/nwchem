      subroutine tce_mrcc_mul_1(alpha,d_a,beta,d_b,k_c_offset,d_c,iref)
        implicit none
#include "global.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "util.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_mrcc.fh"
#include "tce_main.fh"

      double precision alpha,beta
      integer d_a,d_b,d_c
      integer k_c_offset
      integer nprocs
      integer count
      integer next
      INTEGER NXTASK
      EXTERNAL NXTASK
      INTEGER NXTASKsub
      EXTERNAL NXTASKsub
      logical nodezero
      integer iref
      integer p1b,h2b,p1,h2
      integer inoab,invab
      integer size
      integer k_a,k_b,k_c
      integer l_a,l_b,l_c
      integer i
      logical noloadbalance

      nodezero = (ga_nodeid().eq.0)
      count = 0
      if(lusesub) then
      nprocs = ga_pgroup_nnodes(mypgid)
      call ga_pgroup_sync(mypgid)
      next = NXTASKsub(nprocs, 1,mypgid)
      else
      nprocs = ga_nnodes()
      call ga_sync()
      next = NXTASK(nprocs, 1)
      endif
      noloadbalance = ((ioalg.eq.4).or.
     1                ((ioalg.eq.6).and.(.not.fileisga(d_c))))

      inoab = nblcks(1,iref)+nblcks(2,iref)
      invab = nblcks(3,iref)+nblcks(4,iref)

      do p1b = inoab+1,inoab+invab
        do h2b = 1,inoab

          if (noloadbalance.or.(next.eq.count)) then

            if (int_mb(k_spinm(iref)+p1b-1) .eq.
     1  int_mb(k_spinm(iref)+h2b-1)) then
            if ((.not.restricted).or.(int_mb(k_spinm(iref)+p1b-1)
     1        +int_mb(k_spinm(iref)+h2b-1).ne.4)) then
            if (ieor(int_mb(k_symm(iref)+p1b-1),
     1 int_mb(k_symm(iref)+h2b-1)).eq.0) then

              size = int_mb(k_rangem(iref)+p1b-1) *
     1  int_mb(k_rangem(iref)+h2b-1)

              if (.not.ma_push_get(mt_dbl,size,'a',l_a,k_a))
     1          call errquit('tce_jacobi_t1: MA problem',0,MA_ERR)
              if (.not.ma_push_get(mt_dbl,size,'b',l_b,k_b))
     1          call errquit('tce_jacobi_t1: MA problem',0,MA_ERR)
              if (.not.ma_push_get(mt_dbl,size,'c',l_c,k_c))
     1          call errquit('tce_jacobi_t1: MA problem',0,MA_ERR)

              call get_hash_block(d_a,dbl_mb(k_a),size,
     1          int_mb(k_c_offset),((p1b-inoab-1)*inoab+h2b-1))

              call get_hash_block(d_b,dbl_mb(k_b),size,
     1          int_mb(k_c_offset),((p1b-inoab-1)*inoab+h2b-1))

              i = 0

              do p1 = 1,int_mb(k_rangem(iref)+p1b-1)
                do h2 = 1,int_mb(k_rangem(iref)+h2b-1)
                  i = i + 1
                  dbl_mb(k_c+i-1) = alpha*dbl_mb(k_a+i-1)+
     1                               beta*dbl_mb(k_b+i-1)
c        write(6,"('alpha/beta',2F16.12,2F16.12)")alpha,beta
c        write(6,"(I2,':',I4,I4,I4,I4,2F16.12,2F16.12,2F16.12)")
c     1ga_nodeid(),p1b,h2b,i,iref,
c     1 dbl_mb(k_a+i-1),dbl_mb(k_b+i-1),dbl_mb(k_c+i-1)
                enddo
              enddo

              call put_hash_block(d_c,dbl_mb(k_c),size,
     1          int_mb(k_c_offset),((p1b-inoab-1)*inoab+h2b-1))

              if (.not.ma_pop_stack(l_c))
     1          call errquit('tce_mrcc_mul: MA problem',1,MA_ERR)
              if (.not.ma_pop_stack(l_b))
     1          call errquit('tce_mrcc_mul: MA problem',1,MA_ERR)
              if (.not.ma_pop_stack(l_a))
     1          call errquit('tce_mrcc_mul: MA problem',1,MA_ERR)
            endif
            endif
            endif
cc            next = nxtask(nprocs,1)
      if(lusesub) then
       next = NXTASKsub(nprocs,1,mypgid)
      else
       next = NXTASK(nprocs, 1)
      endif
          endif
          count = count + 1
        enddo
      enddo
cc    next = nxtask(-nprocs,1)
      if(lusesub) then
       next = NXTASKsub(-nprocs,1,mypgid)
       call ga_pgroup_sync(mypgid)
      else
       next = NXTASK(-nprocs, 1)
       call ga_sync()
      endif

      return
      end

c $Id$
