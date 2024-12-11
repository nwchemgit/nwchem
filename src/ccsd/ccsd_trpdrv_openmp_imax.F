!#define USE_CPU_REDUCTION

      module reduction_data
        real(8) :: emp4i = 0.0
        real(8) :: emp4k = 0.0
        real(8) :: emp5i = 0.0
        real(8) :: emp5k = 0.0
        !$omp declare target to(emp4i,emp4k,emp5i,emp5ki)
      end module reduction_data


      subroutine ccsd_trpdrv_offload_xe(t1,xeorb,
     &     g_objo,g_objv,g_coul,g_exch,
     &     ncor,nocc,nvir,iprt,emp4,emp5,
     &     oseg_lo,oseg_hi, kchunk)
      use iso_fortran_env
!      use cudafor
!      use cublas
      use omp_lib
      use onemkl_blas_omp_offload_ilp64
      use reduction_data
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

!      double precision, allocatable, device :: eorb(:)
!      double precision, allocatable, device :: f1n(:,:)
!      double precision, allocatable, device :: f2n(:,:)
!      double precision, allocatable, device :: f3n(:,:)
!      double precision, allocatable, device :: f4n(:,:)
!      double precision, allocatable, device :: f1t(:,:)
!      double precision, allocatable, device :: f2t(:,:)
!      double precision, allocatable, device :: f3t(:,:)
!      double precision, allocatable, device :: f4t(:,:)
      double precision, allocatable :: eorb(:)
      double precision, allocatable :: f1n(:,:)
      double precision, allocatable :: f2n(:,:)
      double precision, allocatable :: f3n(:,:)
      double precision, allocatable :: f4n(:,:)
      double precision, allocatable :: f1t(:,:)
      double precision, allocatable :: f2t(:,:)
      double precision, allocatable :: f3t(:,:)
      double precision, allocatable :: f4t(:,:)

!      double precision, allocatable, pinned :: Tij(:), Tkj(:)
!      double precision, allocatable, pinned :: Tia(:), Tka(:)
!      double precision, allocatable, pinned :: Xia(:), Xka(:)
!      double precision, allocatable, pinned :: Jia(:), Jka(:)
!      double precision, allocatable, pinned :: Jij(:), Jkj(:)
!      double precision, allocatable, pinned :: Kia(:), Kka(:)
!      double precision, allocatable, pinned :: Kij(:), Kkj(:)
!      double precision, allocatable, pinned :: Dja(:), Djka(:), Djia(:)
      double precision, allocatable :: Tij(:), Tkj(:)
      double precision, allocatable :: Tia(:), Tka(:)
      double precision, allocatable :: Xia(:), Xka(:)
      double precision, allocatable :: Jia(:), Jka(:)
      double precision, allocatable :: Jij(:), Jkj(:)
      double precision, allocatable :: Kia(:), Kka(:)
      double precision, allocatable :: Kij(:), Kkj(:)
      double precision, allocatable :: Dja(:), Djka(:), Djia(:)

!      double precision, allocatable, device :: xTij(:), xTkj(:)
!      double precision, allocatable, device :: xTia(:), xTka(:)
!      double precision, allocatable, device :: xXia(:), xXka(:)
!      double precision, allocatable, device :: xJia(:), xJka(:)
!      double precision, allocatable, device :: xJij(:), xJkj(:)
!      double precision, allocatable, device :: xKia(:), xKka(:)
!      double precision, allocatable, device :: xKij(:), xKkj(:)
! used to make inline threaded tengy correct - for now
! it is correct that dint[cx]1 are paired with t1v2 and vice versa
! in the inlined tengy loops.  see ccsd_tengy in ccsd_trpdrv.F for
! verification of the i-k and k-i pairing of these.
!      double precision, allocatable, device :: dintc1(:),dintc2(:)
!      double precision, allocatable, device :: dintx1(:),dintx2(:)
!      double precision, allocatable, device :: t1v1(:),t1v2(:)
      double precision, allocatable :: dintc1(:),dintc2(:)
      double precision, allocatable :: dintx1(:),dintx2(:)
      double precision, allocatable :: t1v1(:),t1v2(:)
      integer :: alloc_error, err
!
      !double precision :: emp4i,emp5i,emp4k,emp5k
      !double precision, device :: eaijk
      double precision :: eaijk
      double precision :: denom
      integer :: inode,next,nodes,me
      integer :: a,b,c,i,j,k,akold,av,nbf
      integer :: klo, khi
      integer nxtask
      external nxtask
      integer :: dgemm_flops, tengy_flops
      double precision agg_flops

#ifdef USE_BATCHDGEMM_TRPDRV
      integer, parameter :: batch_size=8
      integer(KIND=C_SIZE_T) :: a_array1(batch_size)
      integer(KIND=C_SIZE_T) :: b_array1(batch_size)
      integer(KIND=C_SIZE_T) :: c_array1(batch_size)
      integer(KIND=C_SIZE_T) :: a_array2(batch_size)
      integer(KIND=C_SIZE_T) :: b_array2(batch_size)
      integer(KIND=C_SIZE_T) :: c_array2(batch_size)
      
      integer, parameter :: grp_count=2
      integer            :: group_size(grp_count)
      character*1        :: transa_array(grp_count)
      character*1        :: transb_array(grp_count)
      integer            :: m_array(grp_count)
      integer            :: n_array(grp_count)
      integer            :: k_array(grp_count)
      integer            :: lda_array(grp_count)
      integer            :: ldb_array(grp_count)
      integer            :: ldc_array(grp_count)
      double precision   :: alpha_array(grp_count)
      double precision   :: beta_array(grp_count)
      integer            :: group_count
#endif
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
      integer thread_num

!     timers
      double precision :: tt0, tt1, tc0, tc1
      double precision :: t_dgemm0, t_dgemm1, t_dgemm_total
      double precision :: t_red0, t_red1, t_red_total

#ifndef USE_BATCHDGEMM_TRPDRV
!     OpenMP interop objects
      integer(kind = omp_interop_kind) :: obj0 = omp_interop_none
      integer(kind = omp_interop_kind) :: obj1 = omp_interop_none
      integer(kind = omp_interop_kind) :: obj2 = omp_interop_none
      integer(kind = omp_interop_kind) :: obj3 = omp_interop_none
      integer(kind = omp_interop_kind) :: obj4 = omp_interop_none
      integer(kind = omp_interop_kind) :: obj5 = omp_interop_none
      integer(kind = omp_interop_kind) :: obj6 = omp_interop_none
      integer(kind = omp_interop_kind) :: obj7 = omp_interop_none
      integer(kind = omp_interop_kind) :: obj_lev0 = omp_interop_none

      !$omp interop init(prefer_type("sycl"),targetsync: obj0)
      !$omp interop init(prefer_type("sycl"),targetsync: obj1)
      !$omp interop init(prefer_type("sycl"),targetsync: obj2)
      !$omp interop init(prefer_type("sycl"),targetsync: obj3)
      !$omp interop init(prefer_type("sycl"),targetsync: obj4)
      !$omp interop init(prefer_type("sycl"),targetsync: obj5)
      !$omp interop init(prefer_type("sycl"),targetsync: obj6)
      !$omp interop init(prefer_type("sycl"),targetsync: obj7)

      !$omp interop init(prefer_type("level_zero"),targetsync: obj_lev0)
#endif

      t_dgemm_total = 0.0
      t_red_total = 0.0

!
      nodes = ga_nnodes()
      me = ga_nodeid()
!
      if (me.eq.0) then
        write(6,99)
      endif
!  99 format(2x,'Using Fortran OpenACC+CUBLAS in CCSD(T)')
   99 format(2x,'Using Fortran OpenMP+MKL in CCSD(T)')
      agg_flops = 0
!
      tt0 = util_wallsec()
!      ! setup CUDA streams
!      do shi=1,8
!        err = cudaStreamCreate(stream(shi))
!        if (err.ne.0) call errquit('cudaStreamCreate',err,UNKNOWN_ERR)
!        err = cublasCreate(handle(shi))
!        if (err.ne.0) call errquit('cublasCreate',err,UNKNOWN_ERR)
!        err = cublasSetStream(handle(shi), stream(shi))
!        if (err.ne.0) call errquit('cublasSetStream',err,UNKNOWN_ERR)
!      end do
!
! device-only temp arrays
! produced by DGEMM, consumed by TENGY
!
!      allocate( f1n(1:nvir,1:nvir), f1t(1:nvir,1:nvir),
!     &          f2n(1:nvir,1:nvir), f2t(1:nvir,1:nvir),
!     &          f3n(1:nvir,1:nvir), f3t(1:nvir,1:nvir),
!     &          f4n(1:nvir,1:nvir), f4t(1:nvir,1:nvir),
!     &          stat=alloc_error)
!      if (alloc_error.ne.0) call errquit('f[1234][tn]',8,MA_ERR)

      !$omp allocate allocator(omp_target_device_mem_alloc)
      allocate( f1n(1:nvir,1:nvir) )
      !$omp allocate allocator(omp_target_device_mem_alloc)
      allocate( f1t(1:nvir,1:nvir) )
      !$omp allocate allocator(omp_target_device_mem_alloc)
      allocate( f2n(1:nvir,1:nvir) )
      !$omp allocate allocator(omp_target_device_mem_alloc)
      allocate( f2t(1:nvir,1:nvir) )
      !$omp allocate allocator(omp_target_device_mem_alloc)
      allocate( f3n(1:nvir,1:nvir) )
      !$omp allocate allocator(omp_target_device_mem_alloc)
      allocate( f3t(1:nvir,1:nvir) )
      !$omp allocate allocator(omp_target_device_mem_alloc)
      allocate( f4n(1:nvir,1:nvir) )
      !$omp allocate allocator(omp_target_device_mem_alloc)
      allocate( f4t(1:nvir,1:nvir) )
!
! device-only copy of input eorb
!
      nbf = ncor + nocc + nvir
      !$omp allocate allocator(omp_target_host_mem_alloc)
      allocate( eorb(1:nbf) )
      ! allocate( eorb(1:nbf), stat=alloc_error)
      ! if (alloc_error.ne.0) call errquit('eorb',10,MA_ERR)
      eorb(1:nbf) = xeorb(1:nbf)
!
! for TENGY
!
!      allocate( dintc1(1:nvir), dintc2(1:nvir),
!     &          dintx1(1:nvir), dintx2(1:nvir),
!     &          t1v1(1:nvir), t1v2(1:nvir), stat=alloc_error)
!      if (alloc_error.ne.0) call errquit('CXT1 temps',6,MA_ERR)

      !$omp allocate allocator(omp_target_host_mem_alloc)
      allocate( dintc1(1:nvir) )
      !$omp allocate allocator(omp_target_host_mem_alloc)
      allocate( dintc2(1:nvir) )
      !$omp allocate allocator(omp_target_host_mem_alloc)
      allocate( dintx1(1:nvir) )
      !$omp allocate allocator(omp_target_host_mem_alloc)
      allocate( dintx2(1:nvir) )
      !$omp allocate allocator(omp_target_host_mem_alloc)
      allocate( t1v1(1:nvir) )
      !$omp allocate allocator(omp_target_host_mem_alloc)
      allocate( t1v2(1:nvir) )

!
! H/D arrays, produced by GA Get, consumed by DGEMM
!
!      allocate( Tij(1:lnvv),      Tkj(1:kchunk*lnvv),
!     &          Tia(1:lnov*nocc), Tka(1:kchunk*lnov),
!     &          Xia(1:lnov*nocc), Xka(1:kchunk*lnov),
!     &          Jia(1:lnvv),      Jka(1:kchunk*lnvv),
!     &          Kia(1:lnvv),      Kka(1:kchunk*lnvv),
!     &          Jij(1:lnov*nocc), Jkj(1:kchunk*lnov),
!     &          Kij(1:lnov*nocc), Kkj(1:kchunk*lnov),
!     &          Dja(1:lnov),      Djka(1:nvir*kchunk),
!     &          Djia(1:nvir*nocc), stat=alloc_error)
!      if (alloc_error.ne.0) call errquit('TKJKD alloc',1,MA_ERR)

      !$omp allocate allocator(omp_target_host_mem_alloc)
      allocate( Tij(1:lnvv) )
      !$omp allocate allocator(omp_target_host_mem_alloc)
      allocate( Tkj(1:kchunk*lnvv) )
      !$omp allocate allocator(omp_target_host_mem_alloc)
      allocate( Tia(1:lnov*nocc) )
      !$omp allocate allocator(omp_target_host_mem_alloc)
      allocate( Tka(1:kchunk*lnov) )
      !$omp allocate allocator(omp_target_host_mem_alloc)
      allocate( Xia(1:lnov*nocc) )
      !$omp allocate allocator(omp_target_host_mem_alloc)
      allocate( Xka(1:kchunk*lnov) )
      !$omp allocate allocator(omp_target_host_mem_alloc)
      allocate( Jia(1:lnvv) )
      !$omp allocate allocator(omp_target_host_mem_alloc)
      allocate( Jka(1:kchunk*lnvv) )
      !$omp allocate allocator(omp_target_host_mem_alloc)
      allocate( Kia(1:lnvv) )
      !$omp allocate allocator(omp_target_host_mem_alloc)
      allocate( Kka(1:kchunk*lnvv) )
      !$omp allocate allocator(omp_target_host_mem_alloc)
      allocate( Jij(1:lnov*nocc) )
      !$omp allocate allocator(omp_target_host_mem_alloc)
      allocate( Jkj(1:kchunk*lnov) )
      !$omp allocate allocator(omp_target_host_mem_alloc)
      allocate( Kij(1:lnov*nocc) )
      !$omp allocate allocator(omp_target_host_mem_alloc)
      allocate( Kkj(1:kchunk*lnov) )
      !$omp allocate allocator(omp_target_host_mem_alloc)
      allocate( Dja(1:lnov) )
      !$omp allocate allocator(omp_target_host_mem_alloc)
      allocate( Djka(1:nvir*kchunk) )
      !$omp allocate allocator(omp_target_host_mem_alloc)
      allocate(Djia(1:nvir*nocc) )

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

!$omp target data
     &  map(to:eorb)
     &  map(to:dintc1,dintx1,t1v1,dintc2,dintx2,t1v2)
     &  map(to:Jia,Tkj,Tia,Kkj,Kia,Tka,Kij)
     &  map(to:Xka,Jij,Xia,Jkj,Jka,Tij,Kka)
#ifdef USE_BATCHDGEMM_TRPDRV
     &  map(to:a_array1, b_array1, c_array1)
     &  map(to:a_array2, b_array2, c_array2)
#endif
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
!$omp target update to (t1v2, dintc1, dintx1)

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
!$omp target update to (t1v1, dintc2, dintx2)
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
                           !!xJia  = Jia
                           !err = cudaMemcpyAsync(xJia,Jia,size(Jia),stream(1))
                           !if (err.ne.0) then
                           !  call errquit('cudaMemcpyAsync',err,UNKNOWN_ERR)
                           !endif
!!$omp target update to (Jia)
                           call ga_nbwait(nbh_objv6)
                           !!xTia  = Tia
                           !err = cudaMemcpyAsync(xTia,Tia,size(Tia),stream(1))
                           !if (err.ne.0) then
                           !  call errquit('cudaMemcpyAsync',err,UNKNOWN_ERR)
                           !endif
!!$omp target update to (Tia)
                           call ga_nbwait(nbh_exch2)
                           !!xKia  = Kia
                           !err = cudaMemcpyAsync(xKia,Kia,size(Kia),stream(2))
                           !if (err.ne.0) then
                           !  call errquit('cudaMemcpyAsync',err,UNKNOWN_ERR)
                           !endif
!!$omp target update to (Kia)
                           call ga_nbwait(nbh_objv7)
                           !!xXia  = Xia
                           !err = cudaMemcpyAsync(xXia,Xia,size(Xia),stream(2))
                           !if (err.ne.0) then
                           !  call errquit('cudaMemcpyAsync',err,UNKNOWN_ERR)
                           !endif
!!$omp target update to (Xia)
!$omp target update to (Jia,Tia,Kia,Xia)
                           if (i.eq.1) then
                              call ga_nbwait(nbh_objo1)
                              !!xTkj  = Tkj
                              !err = cudaMemcpyAsync(xTkj,Tkj,size(Tkj),stream(1))
                              !if (err.ne.0) then
                              !  call errquit('cudaMemcpyAsync',err,UNKNOWN_ERR)
                              !endif
!!$omp target update to (Tkj)
                              call ga_nbwait(nbh_objo3)
                              !!xKkj  = Kkj
                              !err = cudaMemcpyAsync(xKkj,Kkj,size(Kkj),stream(1))
                              !if (err.ne.0) then
                              !  call errquit('cudaMemcpyAsync',err,UNKNOWN_ERR)
                              !endif
!!$omp target update to (Kkj)
                              call ga_nbwait(nbh_objo2)
                              !!xJkj  = Jkj
                              !err = cudaMemcpyAsync(xJkj,Jkj,size(Jkj),stream(3))
                              !if (err.ne.0) then
                              !  call errquit('cudaMemcpyAsync',err,UNKNOWN_ERR)
                              !endif
!!$omp target update to (Jkj)
!$omp target update to (Tkj,Kkj,Jkj)

                              if (got_ak) then
                                 call ga_nbwait(nbh_coul1)
                                 !!xJka  = Jka
                                 !err = cudaMemcpyAsync(xJka,Jka,size(Jka),stream(5))
                                 !if (err.ne.0) then
                                 !  call errquit('cudaMemcpyAsync',err,UNKNOWN_ERR)
                                 !endif
!!$omp target update to (Jka)
                                 call ga_nbwait(nbh_objv2)
                                 !!xTka  = Tka
                                 !err = cudaMemcpyAsync(xTka,Tka,size(Tka),stream(5))
                                 !if (err.ne.0) then
                                 !  call errquit('cudaMemcpyAsync',err,UNKNOWN_ERR)
                                 !endif
!!$omp target update to (Tka)
                                 call ga_nbwait(nbh_exch1)
                                 !!xKka  = Kka
                                 !err = cudaMemcpyAsync(xKka,Kka,size(Kka),stream(6))
                                 !if (err.ne.0) then
                                 !  call errquit('cudaMemcpyAsync',err,UNKNOWN_ERR)
                                 !endif
!!$omp target update to (Kka)
                                 call ga_nbwait(nbh_objv3)
                                 !!xXka  = Xka
                                 !err = cudaMemcpyAsync(xXka,Xka,size(Xka),stream(6))
                                 !if (err.ne.0) then
                                 !  call errquit('cudaMemcpyAsync',err,UNKNOWN_ERR)
                                 !endif
!!$omp target update to (Xka)
!$omp target update to (Jka,Tka,Kka,Xka)
                              endif ! got_ak
                           endif ! i==1
                           call ga_nbwait(nbh_objo4)
                           !!xTij  = Tij
                           !err = cudaMemcpyAsync(xTij,Tij,size(Tij),stream(5))
                           !if (err.ne.0) then
                           !  call errquit('cudaMemcpyAsync',err,UNKNOWN_ERR)
                           !endif
!!$omp target update to (Tij)
                           call ga_nbwait(nbh_objo6)
                           !!xKij  = Kij
                           !err = cudaMemcpyAsync(xKij,Kij,size(Kij),stream(5)) ! and 6
                           !if (err.ne.0) then
                           !  call errquit('cudaMemcpyAsync',err,UNKNOWN_ERR)
                           !endif
!!$omp target update to (Kij)
                           call ga_nbwait(nbh_objo5)
                           !!xJij  = Jij
                           !err = cudaMemcpyAsync(xJij,Jij,size(Jij),stream(7)) ! and 8
                           !if (err.ne.0) then
                           !  call errquit('cudaMemcpyAsync',err,UNKNOWN_ERR)
                           !endif
!!$omp target update to (Jij)
!$omp target update to (Tij,Kij,Jij)

                           ! arrays and thus copies contribute to more than one CUBLAS call
                           ! but the copies on streams 1:4 and 5:8 are separable.
                           !do shi=1,4
                           ! err = cudaStreamSynchronize(stream(shi))
                           ! if (err.ne.0) then
                           !  call errquit('cudaStreamSync',err,UNKNOWN_ERR)
                           ! endif
                           !enddo
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


              tc0 = util_wallsec()
              t_dgemm0 = util_wallsec()

#ifdef USE_BATCHDGEMM_TRPDRV
              ! Update the arrays of pointers to arrays (a_array1, b_array1,
              ! c_array1, a_array2, b_array2, and c_array2) on the device.
              !
              !$omp target
                a_array1(1) = LOC(Jia(1))
                a_array1(2) = LOC(Kia(1))
                a_array1(3) = LOC(Jka(1+(k-klo)*lnvv))
                a_array1(4) = LOC(Kka(1+(k-klo)*lnvv))
                a_array1(5) = LOC(Jia(1))
                a_array1(6) = LOC(Kia(1))
                a_array1(7) = LOC(Jka(1+(k-klo)*lnvv))
                a_array1(8) = LOC(Kka(1+(k-klo)*lnvv))

                b_array1(1) = LOC(Tkj(1+(k-klo)*lnvv))
                b_array1(2) = LOC(Tkj(1+(k-klo)*lnvv))
                b_array1(3) = LOC(Tij(1))
                b_array1(4) = LOC(Tij(1))
                b_array1(5) = LOC(Tkj(1+(k-klo)*lnvv))
                b_array1(6) = LOC(Tkj(1+(k-klo)*lnvv))
                b_array1(7) = LOC(Tij(1))
                b_array1(8) = LOC(Tij(1))

                c_array1(1) = LOC(f1n(1,1))
                c_array1(2) = LOC(f2n(1,1))
                c_array1(3) = LOC(f1t(1,1))
                c_array1(4) = LOC(f2t(1,1))
                c_array1(5) = LOC(f3n(1,1))
                c_array1(6) = LOC(f4n(1,1))
                c_array1(7) = LOC(f3t(1,1))
                c_array1(8) = LOC(f4t(1,1))

                a_array2(1) = LOC(Tia(1))
                a_array2(2) = LOC(Xia(1))
                a_array2(3) = LOC(Tia(1))
                a_array2(4) = LOC(Xia(1))
                a_array2(5) = LOC(Tka(1+(k-klo)*lnov))
                a_array2(6) = LOC(Xka(1+(k-klo)*lnov))
                a_array2(7) = LOC(Tka(1+(k-klo)*lnov))
                a_array2(8) = LOC(Xka(1+(k-klo)*lnov))

                b_array2(1) = LOC(Kkj(1+(k-klo)*lnov))
                b_array2(2) = LOC(Kkj(1+(k-klo)*lnov))
                b_array2(3) = LOC(Jkj(1+(k-klo)*lnov))
                b_array2(4) = LOC(Jkj(1+(k-klo)*lnov))
                b_array2(5) = LOC(Kij(1))
                b_array2(6) = LOC(Kij(1))
                b_array2(7) = LOC(Jij(1))
                b_array2(8) = LOC(Jij(1))

                c_array2(1) = LOC(f1n(1,1))
                c_array2(2) = LOC(f2n(1,1))
                c_array2(3) = LOC(f3n(1,1))
                c_array2(4) = LOC(f4n(1,1))
                c_array2(5) = LOC(f1t(1,1))
                c_array2(6) = LOC(f2t(1,1))
                c_array2(7) = LOC(f3t(1,1))
                c_array2(8) = LOC(f4t(1,1))
              !$omp end target

              !
              ! First call to dgemm_batch (2 groups, each of size 4)
              !

              ! First group in first batch call.
              ! 4 matrix multiplications. So group size is 4.
              group_count = 2
              group_size(1) = 4

              transa_array(1) = 'n'
              transb_array(1) = 't'

              m_array(1) = nvir
              n_array(1) = nvir
              k_array(1) = nvir

              lda_array(1) = nvir
              ldb_array(1) = nvir
              ldc_array(1) = nvir

              alpha_array(1) = 1.0d0
              beta_array(1)  = 0.0d0

             ! Second group in first batch call.
             ! 4 matrix multiplications. So group size is 4.
             group_size(2) = 4

             transa_array(2) = 'n'
             transb_array(2) = 'n'

             m_array(2) = nvir
             n_array(2) = nvir
             k_array(2) = nvir

             lda_array(2) = nvir
             ldb_array(2) = nvir
             ldc_array(2) = nvir

             alpha_array(2) = 1.0d0
             beta_array(2)  = 0.0d0

             !$omp dispatch
               call dgemm_batch(transa_array,
     &                          transb_array,
     &                          m_array,
     &                          n_array,
     &                          k_array,
     &                          alpha_array,
     &                          a_array1,
     &                          lda_array,
     &                          b_array1,
     &                          ldb_array,
     &                          beta_array,
     &                          c_array1,
     &                          ldc_array,
     &                          group_count,
     &                          group_size)

             !
             ! Second call to dgemm_batch (1 group of size 8)
             !

             ! 8 matrix multiplications. So group size is 8.
             group_count = 1
             group_size(1) = 8

             transa_array(1) = 'n'
             transb_array(1) = 'n'

             m_array(1) = nvir
             n_array(1) = nvir
             k_array(1) = nocc

             lda_array(1) = nvir
             ldb_array(1) = nocc
             ldc_array(1) = nvir

             alpha_array(1) = -1.0d0
             beta_array(1) = 1.0d0

             !$omp dispatch
               call dgemm_batch(transa_array,
     &                          transb_array,
     &                          m_array,
     &                          n_array,
     &                          k_array,
     &                          alpha_array,
     &                          a_array2,
     &                          lda_array,
     &                          b_array2,
     &                          ldb_array,
     &                          beta_array,
     &                          c_array2,
     &                          ldc_array,
     &                          group_count,
     &                          group_size)
#else
              !$omp parallel sections num_threads(8)

                !$omp section
                    !$omp dispatch interop(obj0) nowait
                        call dgemm('n','t',nvir,nvir,nvir,1.0d0,
     1                      Jia,nvir,Tkj(1+(k-klo)*lnvv),nvir,0.0d0,
     2                      f1n,nvir)
                  
                    !$omp interop use(obj0) nowait

                    !$omp dispatch interop(obj0) nowait
                        call dgemm('n','n',nvir,nvir,nocc,-1.0d0,
     1                      Tia,nvir,Kkj(1+(k-klo)*lnov),nocc,1.0d0,
     2                      f1n,nvir)

                !$omp section
                    !$omp dispatch interop(obj1) nowait
                        call dgemm('n','t',nvir,nvir,nvir,1.0d0,
     1                      Kia,nvir,Tkj(1+(k-klo)*lnvv),nvir,0.0d0,
     2                      f2n,nvir)

                    !$omp interop use(obj1) nowait

                    !$omp dispatch interop(obj1) nowait
                        call dgemm('n','n',nvir,nvir,nocc,-1.0d0,
     1                      Xia,nvir,Kkj(1+(k-klo)*lnov),nocc,1.0d0,
     2                      f2n,nvir)
                
                !$omp section
                    !$omp dispatch interop(obj2) nowait
                        call dgemm('n','n',nvir,nvir,nvir,1.0d0,
     1                      Jia,nvir,Tkj(1+(k-klo)*lnvv),nvir,0.0d0,
     2                      f3n,nvir)

                    !$omp interop use(obj2) nowait

                    !$omp dispatch interop(obj2) nowait
                        call dgemm('n','n',nvir,nvir,nocc,-1.0d0,
     1                      Tia,nvir,Jkj(1+(k-klo)*lnov),nocc,1.0d0,
     2                      f3n,nvir)

                !$omp section
                    !$omp dispatch interop(obj3) nowait
                        call dgemm('n','n',nvir,nvir,nvir,1.0d0,
     1                      Kia,nvir,Tkj(1+(k-klo)*lnvv),nvir,0.0d0,
     2                      f4n,nvir)

                    !$omp interop use(obj3) nowait

                    !$omp dispatch interop(obj3) nowait
                        call dgemm('n','n',nvir,nvir,nocc,-1.0d0,
     1                      Xia,nvir,Jkj(1+(k-klo)*lnov),nocc,1.0d0,
     2                      f4n,nvir)

                !$omp section
                    !$omp dispatch interop(obj4) nowait
                        call dgemm('n','t',nvir,nvir,nvir,1.0d0,
     1                      Jka(1+(k-klo)*lnvv),nvir,Tij,nvir,0.0d0,
     2                      f1t,nvir)

                    !$omp interop use(obj4) nowait

                    !$omp dispatch interop(obj4) nowait
                        call dgemm('n','n',nvir,nvir,nocc,-1.0d0,
     1                      Tka(1+(k-klo)*lnov),nvir,Kij,nocc,1.0d0,
     2                      f1t,nvir)

                !$omp section
                    !$omp dispatch interop(obj5) nowait
                        call dgemm('n','t',nvir,nvir,nvir,1.0d0,
     1                      Kka(1+(k-klo)*lnvv),nvir,Tij,nvir,0.0d0,
     2                      f2t,nvir)

                    !$omp interop use(obj5) nowait

                    !$omp dispatch interop(obj5) nowait
                        call dgemm('n','n',nvir,nvir,nocc,-1.0d0,
     1                      Xka(1+(k-klo)*lnov),nvir,Kij,nocc,1.0d0,
     2                      f2t,nvir)

                !$omp section
                    !$omp dispatch interop(obj6) nowait
                        call dgemm('n','n',nvir,nvir,nvir,1.0d0,
     1                      Jka(1+(k-klo)*lnvv),nvir,Tij,nvir,0.0d0,
     2                      f3t,nvir)

                    !$omp interop use(obj6) nowait

                    !$omp dispatch interop(obj6) nowait
                        call dgemm('n','n',nvir,nvir,nocc,-1.0d0,
     1                      Tka(1+(k-klo)*lnov),nvir,Jij,nocc,1.0d0,
     2                      f3t,nvir)
      
                !$omp section
                    !$omp dispatch interop(obj7) nowait
                        call dgemm('n','n',nvir,nvir,nvir,1.0d0,
     1                      Kka(1+(k-klo)*lnvv),nvir,Tij,nvir,0.0d0,
     2                      f4t,nvir)

                    !$omp interop use(obj7) nowait

                    !$omp dispatch interop(obj7) nowait
                        call dgemm('n','n',nvir,nvir,nocc,-1.0d0,
     1                      Xka(1+(k-klo)*lnov),nvir,Jij,nocc,1.0d0,
     2                      f4t,nvir)
              !$omp end parallel sections

              ! "omp interop use(obj_lev0)" generates a device barrier
              ! on level0 queue. This behavior is specific to Intel and
              ! works because of the mapping of level0 queues to
              ! hardware
              ! queues -- only one hardware queue is used.
              !
              ! For both Sycl and Level0 queues:
              ! interop use(obj) nowait == device barrier
              ! interop use(obj)        == device barrier + EventHostSynchronize on host
               
              !$omp interop use(obj_lev0) nowait
#endif
              t_dgemm1 = util_wallsec()
              t_dgemm_total = t_dgemm_total + (t_dgemm1 - t_dgemm0)
                        
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

                        eaijk=eorb(a) - (  eorb(ncor+i)
     &                                    +eorb(ncor+j)
     &                                    +eorb(ncor+k) )

#ifdef USE_YFLOP
      flops_ycount = flops_ycount + nvir*nvir*(
     &                       3 + 2*(
     &                       12 +
     &                       11 +
     &                       11 ) +
     &                       2*27 )
#endif

              t_red0 = util_wallsec()
      
#ifdef USE_CPU_REDUCTION

!$omp target data map(to:f1n,f1t,f2n,f2t,f3n,f3t,f4n,f4t)
!$omp target data map(tofrom: emp4i,emp5i) map(to: eaijk)
!$omp target teams distribute parallel do
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
!$omp end target teams distribute parallel do
!$omp end target data
               
               if (i.ne.k) then
!$omp target data map(tofrom: emp4k,emp5k) map(to: eaijk)
!$omp target teams distribute parallel do
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
!$omp end target teams distribute parallel do
!$omp end target data
               end if ! (i.ne.k)
!$omp end target data

               emp4 = emp4 + emp4i
               emp5 = emp5 + emp5i
               emp4 = emp4 + emp4k
               emp5 = emp5 + emp5k
#else

           
          call ccsd_trpdrv_omp_fbody_reduce_new (f1n, f1t, f2n, f2t,
     &                                 f3n, f3t, f4n, f4t,
     &                                 eorb,
     &                                 ncor, nocc, nvir,
     &                                 emp4, emp5,
     &                                 i, k,
     &                                 eaijk,
     &                                 dintc1, dintx1, t1v1,
     &                                 dintc2, dintx2, t1v2)

#endif

              t_red1 = util_wallsec()
              t_red_total = t_red_total + (t_red1 - t_red0)
               
              tc1 = util_wallsec()

              tengy_flops = nvir*nvir*( 3 + 2*( 12 + 11 + 11 ) + 2*27 )
              agg_flops = agg_flops + tengy_flops

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
     &                 ' at ',(util_wallsec()-tt0),' secs',
     &            ', Time for Dgemms =', (t_dgemm1 - t_dgemm0), ' secs',
     &            ', Time for Reduction =', (t_red1 - t_red0), ' secs'

               call util_flush(6)
 4321          format(a,i8,a,i8,a,i3,a,1pg11.4,a,0pf10.1,a,
     &                a, 0pf10.5, a, a, 0pf10.5, a)

               endif
            endif
               end if
            end do
         end do
      end do

      print *, " "
      print *, "- - - - - - - - - - - - - - - - - - - -  - - - - - - -"
      write(6,4331)
     & 'TOTAL Time for Dgemms =', t_dgemm_total, ' secs;  ',
     & 'TOTAL Time for Reduction =', t_red_total, ' secs'
      print *, "- - - - - - - - - - - - - - - - - -- - - - - - - - - -"
      print *, " "

      call util_flush(6)
 4331 format(a, 0pf10.5, a, a, 0pf10.5, a)
               
! end mapping of all data below
!$omp end target data

#ifndef USE_BATCHDGEMM_TRPDRV
      !$omp interop destroy(obj0)
      !$omp interop destroy(obj1)
      !$omp interop destroy(obj2)
      !$omp interop destroy(obj3)
      !$omp interop destroy(obj4)
      !$omp interop destroy(obj5)
      !$omp interop destroy(obj6)
      !$omp interop destroy(obj7)
      !$omp interop destroy(obj_lev0)
#endif
      
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
!      deallocate( xTij, xTkj, xTia, xTka, xXia, xXka,
!     &            xJia, xJka, xKia, xKka, xJij, xJkj, xKij, xKkj,
!     &            stat=alloc_error)
!      if (alloc_error.ne.0) call errquit('free TKJKD GPU',1,MA_ERR)
!
! CUDA stuff
!
      !do shi=1,8
      !  err = cublasDestroy(handle(shi))
      !  if (err.ne.0) call errquit('cublasDestroy',err,UNKNOWN_ERR)
      !  err = cudaStreamDestroy(stream(shi))
      !  if (err.ne.0) call errquit('cudaStreamDestroy',err,UNKNOWN_ERR)
      !end do
!
      tt1 = util_wallsec()
      if (me.eq.0) then
        write(6,501) tt1-tt0
  501   format('CU+MEM free took ',e15.5,' seconds')
      endif
!
      end
