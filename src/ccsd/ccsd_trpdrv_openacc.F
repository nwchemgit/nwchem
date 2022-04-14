      subroutine ccsd_trpdrv_openacc(t1,xeorb,
     &     g_objo,g_objv,g_coul,g_exch,
     &     ncor,nocc,nvir,iprt,emp4,emp5,
     &     oseg_lo,oseg_hi, kchunk)
      use iso_fortran_env
      use cudafor
      use cublas
      implicit none
!
#include "errquit.fh"
#include "global.fh"
#include "ccsd_len.fh"
#include "ccsdps.fh"
#include "util.fh"
#include "msgids.fh"
#include "yflop.fh"
!
      double precision, intent(inout) :: emp4,emp5
      double precision, intent(in) :: t1(*)
      double precision, intent(in) :: xeorb(*)
      integer, intent(in) :: ncor,nocc,nvir
      integer, intent(in) :: iprt
      integer, intent(in) :: g_objo,g_objv,g_coul,g_exch
      integer, intent(in) :: oseg_lo,oseg_hi, kchunk
      double precision, allocatable, device :: eorb(:)
      double precision, allocatable, device :: f1n(:,:)
      double precision, allocatable, device :: f2n(:,:)
      double precision, allocatable, device :: f3n(:,:)
      double precision, allocatable, device :: f4n(:,:)
      double precision, allocatable, device :: f1t(:,:)
      double precision, allocatable, device :: f2t(:,:)
      double precision, allocatable, device :: f3t(:,:)
      double precision, allocatable, device :: f4t(:,:)
      double precision, allocatable, pinned :: Tij(:), Tkj(:)
      double precision, allocatable, pinned :: Tia(:), Tka(:)
      double precision, allocatable, pinned :: Xia(:), Xka(:)
      double precision, allocatable, pinned :: Jia(:), Jka(:)
      double precision, allocatable, pinned :: Jij(:), Jkj(:)
      double precision, allocatable, pinned :: Kia(:), Kka(:)
      double precision, allocatable, pinned :: Kij(:), Kkj(:)
      double precision, allocatable, pinned :: Dja(:), Djka(:), Djia(:)
      double precision, allocatable, device :: xTij(:), xTkj(:)
      double precision, allocatable, device :: xTia(:), xTka(:)
      double precision, allocatable, device :: xXia(:), xXka(:)
      double precision, allocatable, device :: xJia(:), xJka(:)
      double precision, allocatable, device :: xJij(:), xJkj(:)
      double precision, allocatable, device :: xKia(:), xKka(:)
      double precision, allocatable, device :: xKij(:), xKkj(:)
! used to make inline threaded tengy correct - for now
! it is correct that dint[cx]1 are paired with t1v2 and vice versa
! in the inlined tengy loops.  see ccsd_tengy in ccsd_trpdrv.F for
! verification of the i-k and k-i pairing of these.
      double precision, allocatable, device :: dintc1(:),dintc2(:)
      double precision, allocatable, device :: dintx1(:),dintx2(:)
      double precision, allocatable, device :: t1v1(:),t1v2(:)
      integer :: alloc_error, err
!
      double precision :: emp4i,emp5i,emp4k,emp5k
      double precision, device :: eaijk
      double precision :: denom
      integer :: inode,next,nodes,me
      integer :: a,b,c,i,j,k,akold,av,nbf
      integer :: klo, khi
      integer nxtask
      external nxtask
      integer :: dgemm_flops, tengy_flops
      double precision agg_flops
!
!  Dependencies (global array, local array, handle):
!
!  These are waited on first
!
!      g_objv, Dja,  nbh_objv1
!      g_objv, Djka(1+(k-klo)*nvir), nbh_objv4(k)
!      g_objv, Djia, nbh_objv5
!
!  These are waited on later
!
!      g_objv, Tka,  nbh_objv2
!      g_objv, Xka,  nbh_objv3
!      g_objv, Tia,  nbh_objv6
!      g_objv, Xia,  nbh_objv7
!      g_objo, Tkj,  nbh_objo1
!      g_objo, Jkj,  nbh_objo2
!      g_objo, Kkj,  nbh_objo3
!      g_objo, Tij,  nbh_objo4
!      g_objo, Jij,  nbh_objo5
!      g_objo, Kij,  nbh_objo6
!      g_exch, Kka,  nbh_exch1
!      g_exch, Kia,  nbh_exch2
!      g_coul, Jka,  nbh_coul1
!      g_coul, Jia,  nbh_coul2
!
!  non-blocking handles
!
      integer nbh_objv1,nbh_objv2,nbh_objv3
      integer nbh_objv5,nbh_objv6,nbh_objv7
      integer nbh_objv4(nocc)
!
      integer nbh_objo1,nbh_objo2,nbh_objo3
      integer nbh_objo4,nbh_objo5,nbh_objo6
!
      integer nbh_exch1,nbh_exch2,nbh_coul1,nbh_coul2
      integer n_progr,pct_progr
      parameter(n_progr=20)
      logical i_progr(n_progr+1)
      logical got_ak
      ! timers
      double precision :: tt0, tt1, tc0, tc1
      integer(INT32) :: shi
      type(cublasHandle) :: handle(8)
      integer(kind=cuda_stream_kind) :: stream(8)
      integer(INT32) :: nv4, no4 ! cublasDgemm requires 32-bit integers
      integer(INT32), parameter :: cu_op_n = CUBLAS_OP_N
      integer(INT32), parameter :: cu_op_t = CUBLAS_OP_T
!
      nodes = ga_nnodes()
      me = ga_nodeid()
!
      if (me.eq.0) then
        write(6,99)
      endif
   99 format(2x,'Using Fortran OpenACC+CUBLAS in CCSD(T)')
      agg_flops = 0
!
      tt0 = util_wallsec()
      ! setup CUDA streams
      do shi=1,8
        err = cudaStreamCreate(stream(shi))
        if (err.ne.0) call errquit('cudaStreamCreate',err,UNKNOWN_ERR)
        err = cublasCreate(handle(shi))
        if (err.ne.0) call errquit('cublasCreate',err,UNKNOWN_ERR)
        err = cublasSetStream(handle(shi), stream(shi))
        if (err.ne.0) call errquit('cublasSetStream',err,UNKNOWN_ERR)
      end do
!
! device-only temp arrays
! produced by DGEMM, consumed by TENGY
!
      allocate( f1n(1:nvir,1:nvir), f1t(1:nvir,1:nvir),
     &          f2n(1:nvir,1:nvir), f2t(1:nvir,1:nvir),
     &          f3n(1:nvir,1:nvir), f3t(1:nvir,1:nvir),
     &          f4n(1:nvir,1:nvir), f4t(1:nvir,1:nvir),
     &          stat=alloc_error)
      if (alloc_error.ne.0) call errquit('f[1234][tn]',8,MA_ERR)
!
! device-only copy of input eorb
!
      nbf = ncor + nocc + nvir
      allocate( eorb(1:nbf), stat=alloc_error)
      if (alloc_error.ne.0) call errquit('eorb',10,MA_ERR)
      eorb(1:nbf) = xeorb(1:nbf)
!
! for TENGY
!
      allocate( dintc1(1:nvir), dintc2(1:nvir),
     &          dintx1(1:nvir), dintx2(1:nvir),
     &          t1v1(1:nvir), t1v2(1:nvir), stat=alloc_error)
      if (alloc_error.ne.0) call errquit('CXT1 temps',6,MA_ERR)
!
! H/D arrays, produced by GA Get, consumed by DGEMM
!
      allocate( Tij(1:lnvv),      Tkj(1:kchunk*lnvv),
     &          Tia(1:lnov*nocc), Tka(1:kchunk*lnov),
     &          Xia(1:lnov*nocc), Xka(1:kchunk*lnov),
     &          Jia(1:lnvv),      Jka(1:kchunk*lnvv),
     &          Kia(1:lnvv),      Kka(1:kchunk*lnvv),
     &          Jij(1:lnov*nocc), Jkj(1:kchunk*lnov),
     &          Kij(1:lnov*nocc), Kkj(1:kchunk*lnov),
     &          Dja(1:lnov),      Djka(1:nvir*kchunk),
     &          Djia(1:nvir*nocc), stat=alloc_error)
      if (alloc_error.ne.0) call errquit('TKJKD alloc',1,MA_ERR)
      allocate( xTij(1:lnvv),      xTkj(1:kchunk*lnvv),
     &          xTia(1:lnov*nocc), xTka(1:kchunk*lnov),
     &          xXia(1:lnov*nocc), xXka(1:kchunk*lnov),
     &          xJia(1:lnvv),      xJka(1:kchunk*lnvv),
     &          xKia(1:lnvv),      xKka(1:kchunk*lnvv),
     &          xJij(1:lnov*nocc), xJkj(1:kchunk*lnov),
     &          xKij(1:lnov*nocc), xKkj(1:kchunk*lnov),
     &          stat=alloc_error)
      if (alloc_error.ne.0) call errquit('TKJKD GPU alloc',1,MA_ERR)

      tt1 = util_wallsec()
      if (me.eq.0) then
        write(6,500) tt1-tt0
  500   format('CU+MEM init took ',e15.5,' seconds')
      endif
!
!      call ga_sync() ! ga_sync called just before trpdrv in aoccsd2
!
      if (occsdps) then
         call pstat_on(ps_trpdrv)
      else
         call qenter('trpdrv',0)
      endif
       do klo=1,n_progr+1
          i_progr(klo)=.true.
       enddo
      inode=-1
      next=nxtask(nodes, 1)
      do klo = 1, nocc, kchunk
         akold=0
         khi = min(nocc, klo+kchunk-1)
         do a=oseg_lo,oseg_hi
            av=a-ncor-nocc
            do j=1,nocc
               inode=inode+1
               if (inode.eq.next)then

                  call ga_nbget(g_objv,1+(j-1)*lnov,j*lnov,av,av,Dja,
     &                          lnov,nbh_objv1)
                  do k = klo, khi
                     call ga_nbget(g_objv,1+(j-1)*nvir+(k-1)*lnov,
     &                    j*nvir+(k-1)*lnov,av,av,
     &                    Djka(1+(k-klo)*nvir),nvir,nbh_objv4(k))
                  enddo
                  call ga_nbget(g_objo,(klo-1)*lnvv+1,khi*lnvv,j,j,Tkj,
     &                          (khi-klo+1)*lnvv,nbh_objo1)
                  call ga_nbget(g_objo,lnovv+(klo-1)*lnov+1,
     &                          lnovv+khi*lnov,j,j,Jkj,
     &                          (khi-klo+1)*lnov,nbh_objo2)
                  call ga_nbget(g_objo,lnovv+lnoov+(klo-1)*lnov+1,
     &                          lnovv+lnoov+khi*lnov,j,j,Kkj,
     &                          (khi-klo+1)*lnov,nbh_objo3)
                  if (akold .ne. a) then
                     akold = a
                     call ga_nbget(g_coul,1,lnvv,(a-oseg_lo)*nocc+klo,
     &                    (a-oseg_lo)*nocc+khi,Jka,lnvv,nbh_coul1)
                     call ga_nbget(g_exch,1,lnvv,(a-oseg_lo)*nocc+klo,
     &                    (a-oseg_lo)*nocc+khi,Kka,lnvv,nbh_exch1)
                     call ga_nbget(g_objv,1+lnoov+(klo-1)*lnov,
     &                    lnoov+khi*lnov,av,av,Tka,(khi-klo+1)*lnov,
     &                    nbh_objv2)
                     call ga_nbget(g_objv,1+2*lnoov+(klo-1)*lnov,
     &                    2*lnoov+khi*lnov,av,av,Xka,(khi-klo+1)*lnov,
     &                    nbh_objv3)
                     got_ak=.true.
                  else
                     got_ak=.false.
                  endif

                  emp4i = 0.0d0
                  emp5i = 0.0d0
                  emp4k = 0.0d0
                  emp5k = 0.0d0

                  do i=1,nocc

                     call ga_nbget(g_objv,1+(j-1)*nvir+(i-1)*lnov,
     &                    j*nvir+(i-1)*lnov,av,av,Djia,nvir,nbh_objv5)
                     call ga_nbget(g_objo,(i-1)*lnvv+1,i*lnvv,j,j,Tij,
     &                    lnvv,nbh_objo4)
                     call ga_nbget(g_objo,lnovv+(i-1)*lnov+1,
     &                    lnovv+i*lnov,j,j,Jij,lnov,nbh_objo5)
                     call ga_nbget(g_objo,lnovv+lnoov+(i-1)*lnov+1,
     &                    lnovv+lnoov+i*lnov,j,j,Kij,lnov,nbh_objo6)
                     call ga_nbget(g_coul,1,lnvv,(a-oseg_lo)*nocc+i,
     &                    (a-oseg_lo)*nocc+i,Jia,lnvv,nbh_coul2)
                     call ga_nbget(g_exch,1,lnvv,(a-oseg_lo)*nocc+i,
     &                    (a-oseg_lo)*nocc+i,Kia,lnvv,nbh_exch2)
                     call ga_nbget(g_objv,1+lnoov+(i-1)*lnov,
     &                    lnoov+i*lnov,av,av,Tia,lnov,nbh_objv6)
                     call ga_nbget(g_objv,1+2*lnoov+(i-1)*lnov,
     &                    2*lnoov+i*lnov,av,av,Xia,lnov,nbh_objv7)

                     if (occsdps) then
                        call pstat_on(ps_accwait)
                     else
                        call qenter('accwait',0)
                     endif
!$acc wait(9)
!$acc wait(10)
                     if (occsdps) then
                        call pstat_off(ps_accwait)
                     else
                        call qexit('accwait',0)
                     endif

                     t1v2(:) = t1((i-1)*nvir+1:i*nvir)
                     if(i.eq.1) then
                        call ga_nbwait(nbh_objv1) ! Dja
                     endif
                     dintc1(:) = Dja(1+(i-1)*nvir:i*nvir)
                     call ga_nbwait(nbh_objv5) ! Djia
                     dintx1(:) = Djia(1:nvir)

                     do k=klo,min(khi,i)
                        if (occsdps) then
                           call pstat_on(ps_accwait)
                        else
                           call qenter('accwait',0)
                        endif
!$acc wait(9)           
!$acc wait(10)          
                        if (occsdps) then
                           call pstat_off(ps_accwait)
                        else
                           call qexit('accwait',0)
                        endif

                        t1v1(:) = t1((k-1)*nvir+1:k*nvir)
                        dintc2(:) = Dja(1+(k-1)*nvir:k*nvir)
                        if(i.eq.1) then
                           call ga_nbwait(nbh_objv4(k)) ! Djka
                        endif
                        dintx2(:) = Djka(1+(k-klo)*nvir:(k-klo+1)*nvir)
!
!  These are the input dependencies for the DGEMM calls below.
!  We wait on all of them here because GA is not even remotely thread-safe.
!  All of these are independent of k, so we wait on them only
!  at the first trip of the loop.
!
                        if (k.eq.klo) then
                           if (occsdps) then
                              call pstat_on(ps_gpumove)
                           else
                              call qenter('gpumove',0)
                           endif
                           call ga_nbwait(nbh_coul2)
                           !xJia  = Jia
                           err = cudaMemcpyAsync(xJia,Jia,size(Jia),stream(1))
                           if (err.ne.0) then
                             call errquit('cudaMemcpyAsync',err,UNKNOWN_ERR)
                           endif
                           call ga_nbwait(nbh_objv6)
                           !xTia  = Tia
                           err = cudaMemcpyAsync(xTia,Tia,size(Tia),stream(1))
                           if (err.ne.0) then
                             call errquit('cudaMemcpyAsync',err,UNKNOWN_ERR)
                           endif
                           call ga_nbwait(nbh_exch2)
                           !xKia  = Kia
                           err = cudaMemcpyAsync(xKia,Kia,size(Kia),stream(2))
                           if (err.ne.0) then
                             call errquit('cudaMemcpyAsync',err,UNKNOWN_ERR)
                           endif
                           call ga_nbwait(nbh_objv7)
                           !xXia  = Xia
                           err = cudaMemcpyAsync(xXia,Xia,size(Xia),stream(2))
                           if (err.ne.0) then
                             call errquit('cudaMemcpyAsync',err,UNKNOWN_ERR)
                           endif
                           if (i.eq.1) then
                              call ga_nbwait(nbh_objo1)
                              !xTkj  = Tkj
                              err = cudaMemcpyAsync(xTkj,Tkj,size(Tkj),stream(1))
                              if (err.ne.0) then
                                call errquit('cudaMemcpyAsync',err,UNKNOWN_ERR)
                              endif
                              call ga_nbwait(nbh_objo3)
                              !xKkj  = Kkj
                              err = cudaMemcpyAsync(xKkj,Kkj,size(Kkj),stream(1))
                              if (err.ne.0) then
                                call errquit('cudaMemcpyAsync',err,UNKNOWN_ERR)
                              endif
                              call ga_nbwait(nbh_objo2)
                              !xJkj  = Jkj
                              err = cudaMemcpyAsync(xJkj,Jkj,size(Jkj),stream(3))
                              if (err.ne.0) then
                                call errquit('cudaMemcpyAsync',err,UNKNOWN_ERR)
                              endif
                              if (got_ak) then
                                 call ga_nbwait(nbh_coul1)
                                 !xJka  = Jka
                                 err = cudaMemcpyAsync(xJka,Jka,size(Jka),stream(5))
                                 if (err.ne.0) then
                                   call errquit('cudaMemcpyAsync',err,UNKNOWN_ERR)
                                 endif
                                 call ga_nbwait(nbh_objv2)
                                 !xTka  = Tka
                                 err = cudaMemcpyAsync(xTka,Tka,size(Tka),stream(5))
                                 if (err.ne.0) then
                                   call errquit('cudaMemcpyAsync',err,UNKNOWN_ERR)
                                 endif
                                 call ga_nbwait(nbh_exch1)
                                 !xKka  = Kka
                                 err = cudaMemcpyAsync(xKka,Kka,size(Kka),stream(6))
                                 if (err.ne.0) then
                                   call errquit('cudaMemcpyAsync',err,UNKNOWN_ERR)
                                 endif
                                 call ga_nbwait(nbh_objv3)
                                 !xXka  = Xka
                                 err = cudaMemcpyAsync(xXka,Xka,size(Xka),stream(6))
                                 if (err.ne.0) then
                                   call errquit('cudaMemcpyAsync',err,UNKNOWN_ERR)
                                 endif
                              endif ! got_ak
                           endif ! i==1
                           call ga_nbwait(nbh_objo4)
                           !xTij  = Tij
                           err = cudaMemcpyAsync(xTij,Tij,size(Tij),stream(5))
                           if (err.ne.0) then
                             call errquit('cudaMemcpyAsync',err,UNKNOWN_ERR)
                           endif
                           call ga_nbwait(nbh_objo6)
                           !xKij  = Kij
                           err = cudaMemcpyAsync(xKij,Kij,size(Kij),stream(5)) ! and 6
                           if (err.ne.0) then
                             call errquit('cudaMemcpyAsync',err,UNKNOWN_ERR)
                           endif
                           call ga_nbwait(nbh_objo5)
                           !xJij  = Jij
                           err = cudaMemcpyAsync(xJij,Jij,size(Jij),stream(7)) ! and 8
                           if (err.ne.0) then
                             call errquit('cudaMemcpyAsync',err,UNKNOWN_ERR)
                           endif

                           ! arrays and thus copies contribute to more than one CUBLAS call
                           ! but the copies on streams 1:4 and 5:8 are separable.
                           do shi=1,4
                            err = cudaStreamSynchronize(stream(shi))
                            if (err.ne.0) then
                             call errquit('cudaStreamSync',err,UNKNOWN_ERR)
                            endif
                           enddo
                           if (occsdps) then
                              call pstat_off(ps_gpumove)
                           else
                              call qexit('gpumove',0)
                           endif
                        endif ! k==klo

                        if (occsdps) then
                           call pstat_on(ps_doxxx)
                        else
                           call qenter('doxxx',0)
                        endif

                        tc0 = util_wallsec()

                        nv4 = nvir ! no possibility of overflow
                        no4 = nocc

                        if (occsdps) then
                           call pstat_on(ps_accwait)
                        else
                           call qenter('accwait',0)
                        endif
!$acc wait(9)           
!$acc wait(10)          
                        if (occsdps) then
                           call pstat_off(ps_accwait)
                        else
                           call qexit('accwait',0)
                        endif

                        err = cublasDgemm_v2(handle(1),
     &                        cu_op_n,cu_op_t,
     &                        nv4,nv4,nv4,1.0d0,
     &                        xJia,nv4,xTkj(1+(k-klo)*lnvv),nv4,0.0d0,
     &                        f1n,nv4)
                        if (err.ne.0) then
                          call errquit('cublasDgemm_v2',err,UNKNOWN_ERR)
                        endif
                        err = cublasDgemm_v2(handle(1),
     &                        cu_op_n,cu_op_n,
     &                        nv4,nv4,no4,-1.0d0,
     &                        xTia,nv4,xKkj(1+(k-klo)*lnov),no4,1.0d0,
     &                        f1n,nv4)
                        if (err.ne.0) then
                          call errquit('cublasDgemm_v2',err,UNKNOWN_ERR)
                        endif

                        err = cublasDgemm_v2(handle(2),
     &                        cu_op_n,cu_op_t,
     &                        nv4,nv4,nv4,1.0d0,
     &                        xKia,nv4,xTkj(1+(k-klo)*lnvv),nv4,0.0d0,
     &                        f2n,nv4)
                        if (err.ne.0) then
                          call errquit('cublasDgemm_v2',err,UNKNOWN_ERR)
                        endif
                        err = cublasDgemm_v2(handle(2),
     &                        cu_op_n,cu_op_n,
     &                        nv4,nv4,no4,-1.0d0,
     &                        xXia,nv4,xKkj(1+(k-klo)*lnov),no4,1.0d0,
     &                        f2n,nv4)
                        if (err.ne.0) then
                          call errquit('cublasDgemm_v2',err,UNKNOWN_ERR)
                        endif

                        err = cublasDgemm_v2(handle(3),
     &                        cu_op_n,cu_op_n,
     &                        nv4,nv4,nv4,1.0d0,
     &                        xJia,nv4,xTkj(1+(k-klo)*lnvv),nv4,0.0d0,
     &                        f3n,nv4)
                        if (err.ne.0) then
                          call errquit('cublasDgemm_v2',err,UNKNOWN_ERR)
                        endif
                        err = cublasDgemm_v2(handle(3),
     &                        cu_op_n,cu_op_n,
     &                        nv4,nv4,no4,-1.0d0,
     &                        xTia,nv4,xJkj(1+(k-klo)*lnov),no4,1.0d0,
     &                        f3n,nv4)
                        if (err.ne.0) then
                          call errquit('cublasDgemm_v2',err,UNKNOWN_ERR)
                        endif

                        err = cublasDgemm_v2(handle(4),
     &                        cu_op_n,cu_op_n,
     &                        nv4,nv4,nv4,1.0d0,
     &                        xKia,nv4,xTkj(1+(k-klo)*lnvv),nv4,0.0d0,
     &                        f4n,nv4)
                        if (err.ne.0) then
                          call errquit('cublasDgemm_v2',err,UNKNOWN_ERR)
                        endif
                        err = cublasDgemm_v2(handle(4),
     &                        cu_op_n,cu_op_n,
     &                        nv4,nv4,no4,-1.0d0,
     &                        xXia,nv4,xJkj(1+(k-klo)*lnov),no4,1.0d0,
     &                        f4n,nv4)
                        if (err.ne.0) then
                          call errquit('cublasDgemm_v2',err,UNKNOWN_ERR)
                        endif

                        if (k.eq.klo) then
                           do shi=5,8
                            err = cudaStreamSynchronize(stream(shi))
                            if (err.ne.0) then
                             call errquit('cudaStreamSync',err,UNKNOWN_ERR)
                            endif
                           enddo
                        endif ! k==klo

                        err = cublasDgemm_v2(handle(5),
     &                        cu_op_n,cu_op_t,
     &                        nv4,nv4,nv4,1.0d0,
     &                        xJka(1+(k-klo)*lnvv),nv4,xTij,nv4,0.0d0,
     &                        f1t,nv4)
                        if (err.ne.0) then
                          call errquit('cublasDgemm_v2',err,UNKNOWN_ERR)
                        endif
                        err = cublasDgemm_v2(handle(5),
     &                        cu_op_n,cu_op_n,
     &                        nv4,nv4,no4,-1.0d0,
     &                        xTka(1+(k-klo)*lnov),nv4,xKij,no4,1.0d0,
     &                        f1t,nv4)
                        if (err.ne.0) then
                          call errquit('cublasDgemm_v2',err,UNKNOWN_ERR)
                        endif

                        err = cublasDgemm_v2(handle(6),
     &                        cu_op_n,cu_op_t,
     &                        nv4,nv4,nv4,1.0d0,
     &                        xKka(1+(k-klo)*lnvv),nv4,xTij,nv4,0.0d0,
     &                        f2t,nv4)
                        if (err.ne.0) then
                          call errquit('cublasDgemm_v2',err,UNKNOWN_ERR)
                        endif
                        err = cublasDgemm_v2(handle(6),
     &                        cu_op_n,cu_op_n,
     &                        nv4,nv4,no4,-1.0d0,
     &                        xXka(1+(k-klo)*lnov),nv4,xKij,no4,1.0d0,
     &                        f2t,nv4)
                        if (err.ne.0) then
                          call errquit('cublasDgemm_v2',err,UNKNOWN_ERR)
                        endif

                        err = cublasDgemm_v2(handle(7),
     &                        cu_op_n,cu_op_n,
     &                        nv4,nv4,nv4,1.0d0,
     &                        xJka(1+(k-klo)*lnvv),nv4,xTij,nv4,0.0d0,
     &                        f3t,nv4)
                        if (err.ne.0) then
                          call errquit('cublasDgemm_v2',err,UNKNOWN_ERR)
                        endif
                        err = cublasDgemm_v2(handle(7),
     &                        cu_op_n,cu_op_n,
     &                        nv4,nv4,no4,-1.0d0,
     &                        xTka(1+(k-klo)*lnov),nv4,xJij,no4,1.0d0,
     &                        f3t,nv4)
                        if (err.ne.0) then
                          call errquit('cublasDgemm_v2',err,UNKNOWN_ERR)
                        endif

                        err = cublasDgemm_v2(handle(8),
     &                        cu_op_n,cu_op_n,
     &                        nv4,nv4,nv4,1.0d0,
     &                        xKka(1+(k-klo)*lnvv),nv4,xTij,nv4,0.0d0,
     &                        f4t,nv4)
                        if (err.ne.0) then
                          call errquit('cublasDgemm_v2',err,UNKNOWN_ERR)
                        endif
                        err = cublasDgemm_v2(handle(8),
     &                        cu_op_n,cu_op_n,
     &                        nv4,nv4,no4,-1.0d0,
     &                        xXka(1+(k-klo)*lnov),nv4,xJij,no4,1.0d0,
     &                        f4t,nv4)
                        if (err.ne.0) then
                          call errquit('cublasDgemm_v2',err,UNKNOWN_ERR)
                        endif

                        err = cudaDeviceSynchronize()
                        if (err.ne.0) then
                          call errquit('cudaDeviceSync',err,UNKNOWN_ERR)
                        endif

                        ! 8 pairs of DGEMM w/ VVV and VVO cost, 2 for FMA
                        dgemm_flops = 8*nvir*nvir*(nocc+nvir)*2
                        agg_flops = agg_flops + dgemm_flops

                        if (occsdps) then
                           call pstat_off(ps_doxxx)
                           call pstat_on(ps_tengy)
                        else
                           call qexit('doxxx',0)
                           call qenter('tengy',0)
                        endif

!$acc kernels
                        eaijk=eorb(a) - (  eorb(ncor+i)
     &                                    +eorb(ncor+j)
     &                                    +eorb(ncor+k) )
!$acc end kernels

#ifdef USE_YFLOP
      flops_ycount = flops_ycount + nvir*nvir*(
     &                       3 + 2*(
     &                       12 +
     &                       11 +
     &                       11 ) +
     &                       2*27 )
#endif

!$acc parallel loop tile(32,32) private(denom)
!$acc&         reduction(+:emp4i,emp5i)
!$acc&         async(9)
               do b=1,nvir
                 do c=1,nvir
                   denom=-1.0d0/( eorb(ncor+nocc+b)
     &                           +eorb(ncor+nocc+c)+eaijk )
                   emp4i=emp4i+denom*
     &                  (f1t(b,c)+f1n(c,b)+f2t(c,b)+f3n(b,c)+f4n(c,b))*
     &                  (f1t(b,c)-2*f2t(b,c)-2*f3t(b,c)+f4t(b,c))
     &                        -denom*
     &                  (f1n(b,c)+f1t(c,b)+f2n(c,b)+f3n(c,b))*
     &                  (2*f1t(b,c)-f2t(b,c)-f3t(b,c)+2*f4t(b,c))
     &                        +3*denom*(
     &                  f1n(b,c)*(f1n(b,c)+f3n(c,b)+2*f4t(c,b))+
     &                  f2n(b,c)*f2t(c,b)+f3n(b,c)*f4t(b,c))
                   emp5i=emp5i+denom*t1v1(b)*dintx1(c)*
     &                 (    f1t(b,c)+f2n(b,c)+f4n(c,b)
     &                  -2*(f3t(b,c)+f4n(b,c)+f2n(c,b)+
     &                      f1n(b,c)+f2t(b,c)+f3n(c,b))
     &                  +4*(f3n(b,c)+f4t(b,c)+f1n(c,b)))
     &                        +denom*t1v1(b)*dintc1(c)*
     &                 (     f1n(b,c)+f4n(b,c)+f1t(c,b)
     &                   -2*(f2n(b,c)+f3n(b,c)+f2t(c,b)))
                 end do
               end do
               if (i.ne.k) then
!$acc parallel loop tile(32,32) private(denom)
!$acc&         reduction(+:emp4k,emp5k)
!$acc&         async(10)
               do b=1,nvir
                 do c=1,nvir
                   denom=-1.0d0/( eorb(ncor+nocc+b)
     &                           +eorb(ncor+nocc+c)+eaijk )
                   emp4k=emp4k+denom*
     &                  (f1n(b,c)+f1t(c,b)+f2n(c,b)+f3t(b,c)+f4t(c,b))*
     &                  (f1n(b,c)-2*f2n(b,c)-2*f3n(b,c)+f4n(b,c))
     &                        -denom*
     &                  (f1t(b,c)+f1n(c,b)+f2t(c,b)+f3t(c,b))*
     &                  (2*f1n(b,c)-f2n(b,c)-f3n(b,c)+2*f4n(b,c))
     &                        +3*denom*(
     &                  f1t(b,c)*(f1t(b,c)+f3t(c,b)+2*f4n(c,b))+
     &                  f2t(b,c)*f2n(c,b)+f3t(b,c)*f4n(b,c))
                   emp5k=emp5k+denom*t1v2(b)*dintx2(c)*
     &                 (    f1n(b,c)+f2t(b,c)+f4t(c,b)
     &                  -2*(f3n(b,c)+f4t(b,c)+f2t(c,b)+
     &                      f1t(b,c)+f2n(b,c)+f3t(c,b))
     &                  +4*(f3t(b,c)+f4n(b,c)+f1t(c,b)))
     &                        +denom*t1v2(b)*dintc2(c)*
     &                 (     f1t(b,c)+f4t(b,c)+f1n(c,b)
     &                   -2*(f2t(b,c)+f3t(b,c)+f2n(c,b)))
                 end do
               end do
               end if ! (i.ne.k)
               tengy_flops = nvir*nvir*( 3 + 2*( 12 + 11 + 11 ) + 2*27 )
               agg_flops = agg_flops + tengy_flops

                         tc1 = util_wallsec()

                         if (occsdps) then
                            call pstat_off(ps_tengy)
                         else
                            call qexit('tengy',0)
                         endif

                     end do    ! k
                  end do       ! i

                  if (occsdps) then
                     call pstat_on(ps_accwait)
                  else
                     call qenter('accwait',0)
                  endif
!$acc wait(9)     
!$acc wait(10)    
                  if (occsdps) then
                     call pstat_off(ps_accwait)
                  else
                     call qexit('accwait',0)
                  endif

                  emp4 = emp4 + emp4i
                  emp5 = emp5 + emp5i
                  emp4 = emp4 + emp4k
                  emp5 = emp5 + emp5k
                  if (iprt.gt.50)then
                     write(6,1234)me,a,j,emp4,emp5
 1234                format(' me aijk',3i5,2e15.5)
                  end if
                  next=nxtask(nodes, 1)
            if(me.eq.0) then
               pct_progr=(a-(ncor+nocc)+((klo-1)/kchunk)*nvir)*n_progr/
     &                   ((nocc/kchunk)*nvir)+1
               if(i_progr(pct_progr)) then
                  i_progr(pct_progr)=.false.
               write(6,4321) ' ccsd(t): done ',
     &              a-(ncor+nocc)+((klo-1)/kchunk)*nvir,
     &              ' out of ',(nocc/kchunk)*nvir,
     &              ' progress: ',
     &              ((a-(ncor+nocc)+((klo-1)/kchunk)*nvir)*100)/
     &              ((nocc/kchunk)*nvir),
     &            '%, Gflops=',1e-9*(dgemm_flops+tengy_flops)/(tc1-tc0),
     &                 ' at ',(util_wallsec()-tt0),' secs'
               call util_flush(6)
 4321          format(a,i8,a,i8,a,i3,a,1pg11.4,a,0pf10.1,a)
               endif
            endif
               end if
            end do
         end do
      end do
      call ga_sync()
      next=nxtask(-nodes, 1)
      tt1=util_wallsec()
      call ga_dgop(msg_cc_diis1,agg_flops,1, '+')
      if(me.eq.0) then
         write(6,4322) ' ccsd(t): 100% done, Aggregate Gflops=',
     &        1e-9*agg_flops/(tt1-tt0),' in ',(tt1-tt0),' secs'
 4322    format(a,1pg11.4,a,0pf10.1,a)
         call util_flush(6)
      endif
      call ga_sync()
      if (occsdps) then
         call pstat_off(ps_trpdrv)
      else
         call qexit('trpdrv',0)
      endif
!
      tt0 = util_wallsec()
      deallocate( f1n, f1t, f2n, f2t, f3n, f3t, f4n, f4t,
     &            stat=alloc_error)
      if (alloc_error.ne.0) call errquit('free f[1234][tn]',8,MA_ERR)

      deallocate( eorb, dintc1, dintx1, t1v1, dintc2, dintx2, t1v2,
     &            stat=alloc_error)
      if (alloc_error.ne.0) call errquit('free CXT1 temps',6,MA_ERR)

      deallocate( Tij, Tkj, Tia, Tka, Xia, Xka,
     &            Jia, Jka, Kia, Kka, Jij, Jkj, Kij, Kkj,
     &            Dja, Djka, Djia, stat=alloc_error)
      if (alloc_error.ne.0) call errquit('free TKJKD',1,MA_ERR)
      deallocate( xTij, xTkj, xTia, xTka, xXia, xXka,
     &            xJia, xJka, xKia, xKka, xJij, xJkj, xKij, xKkj,
     &            stat=alloc_error)
      if (alloc_error.ne.0) call errquit('free TKJKD GPU',1,MA_ERR)
!
! CUDA stuff
!
      do shi=1,8
        err = cublasDestroy(handle(shi))
        if (err.ne.0) call errquit('cublasDestroy',err,UNKNOWN_ERR)
        err = cudaStreamDestroy(stream(shi))
        if (err.ne.0) call errquit('cudaStreamDestroy',err,UNKNOWN_ERR)
      end do
!
      tt1 = util_wallsec()
      if (me.eq.0) then
        write(6,501) tt1-tt0
  501   format('CU+MEM free took ',e15.5,' seconds')
      endif
!
      end
