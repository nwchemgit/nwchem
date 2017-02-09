#ifndef TRPMIXP_OFF
      subroutine ccsd_trpdrv_omp_mixprec(t1,
     &     eorb,g_objo,g_objv,g_coul,g_exch,
     &     ncor,nocc,nvir,iprt,emp4,emp5,
     &     oseg_lo,oseg_hi, kchunk,
     &     Tij, Tkj, Tia, Tka, Xia, Xka, Jia, Jka, Kia, Kka,
     &     Jij, Jkj, Kij, Kkj, Dja, Djka, Djia)
!     $Id: ccsd_trpdrv_omp_mixprec.F 28839 2016-11-22 21:26:48Z jhammond $
      implicit none
!
#include "errquit.fh"
#include "global.fh"
#include "ccsd_len.fh"
#include "ccsdps.fh"
!
      integer, parameter :: sp = selected_real_kind(6, 37)
      integer, parameter :: dp = selected_real_kind(15, 307)
!
      real(kind=dp), intent(inout) :: emp4,emp5
      real(kind=dp), intent(in) :: t1(*)
      integer, intent(in) :: ncor,nocc,nvir
      integer, intent(in) :: iprt
      integer, intent(in) :: g_objo,g_objv,g_coul,g_exch
      integer, intent(in) :: oseg_lo,oseg_hi, kchunk
! TODO: allocate a copy of eorb on heap with fastmem for KNL
      real(kind=dp), intent(in) :: eorb(*)
! TODO: determine if these should be allocated w/ fastmem in parent subroutine
      real(kind=dp), intent(in) :: Xia(*), Xka(*)
      real(kind=dp), intent(in) :: Tij(*), Tkj(*), Tia(*), Tka(*)
      real(kind=dp), intent(in) :: Jia(*), Jka(*), Jij(*), Jkj(*)
      real(kind=dp), intent(in) :: Kia(*), Kka(*), Kij(*), Kkj(*)
      real(kind=dp), intent(in) :: Dja(*), Djka(*), Djia(*)
!
! These arrays are touched by GEMM:
!
!    Jia Tia Kia Xia
!    Jka Tka Kka Xka
!    Jij Tij Kij
!    Jkj Tkj Kkj
!
      real(kind=sp), allocatable :: sXia(:), sXka(:)
      real(kind=sp), allocatable :: sTij(:), sTkj(:), sTia(:), sTka(:)
      real(kind=sp), allocatable :: sJia(:), sJka(:), sJij(:), sJkj(:)
      real(kind=sp), allocatable :: sKia(:), sKka(:), sKij(:), sKkj(:)
      real(kind=sp), allocatable :: f1n(:,:), f1t(:,:)
      real(kind=sp), allocatable :: f2n(:,:), f2t(:,:)
      real(kind=sp), allocatable :: f3n(:,:), f3t(:,:)
      real(kind=sp), allocatable :: f4n(:,:), f4t(:,:)
      integer :: alloc_error_sp(22)

! used to make inline threaded tengy correct - for now
#ifdef USE_F90_ALLOCATABLE
      real(kind=sp), allocatable :: dintc1(:),dintc2(:)
      real(kind=sp), allocatable :: dintx1(:),dintx2(:)
      real(kind=sp), allocatable :: t1v1(:),t1v2(:)
#ifdef USE_FASTMEM
!dec$ attributes fastmem :: dintc1,dintc2,dintx1,dintx2,t1v1,t1v2
#endif
      integer :: alloc_error(6)
#else
      real(kind=sp) :: dintc1(nvir),dintx1(nvir),t1v1(nvir)
      real(kind=sp) :: dintc2(nvir),dintx2(nvir),t1v2(nvir)
#endif
!
      real(kind=sp) :: emp4i,emp5i,emp4k,emp5k
      real(kind=dp) :: eaijk
      real(kind=sp) :: denom
      integer :: inode,next,nodes,iam
      integer :: a,b,c,i,j,k,akold,av
      ! chunking is the loop blocking size in the loop nest
      ! formerly associated with the tengy routine.
      ! we have not explored this paramater space but 32 is
      ! optimal for TLB blocking in matrix transpose on most
      ! architectures (especially x86).
      integer, parameter :: chunking = 32
      integer :: bb,cc
      integer :: klo, khi
      integer :: offo, offv
      integer nxtask
      external nxtask
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
!
#if defined(USE_OPENMP)
      integer  omp_get_thread_num
      external omp_get_thread_num
      integer  omp_get_num_threads
      external omp_get_num_threads
      integer  omp_get_max_threads
      external omp_get_max_threads
      if (ga_nodeid().eq.0) write(6,99) omp_get_max_threads()
   99 format(2x,'Using ',i2,' OpenMP threads in CCSD(T)')
#else
      if (ga_nodeid().eq.0) then
        write(6,99) 1
        write(6,999)
      endif
   99 format(2x,'Using ',i2,' OpenMP thread in CCSD(T)')
  999 format(2x,'Recompile w/ USE_OPENMP=1 to use threads in CCSD(T)')
#endif
      allocate( sXia(lnvv),        stat=alloc_error_sp( 1) )
      allocate( sXka(kchunk*lnov), stat=alloc_error_sp( 2) )
      allocate( sTij(lnvv),        stat=alloc_error_sp( 3) )
      allocate( sTkj(kchunk*lnvv), stat=alloc_error_sp( 4) )
      allocate( sTia(lnov*nocc),   stat=alloc_error_sp( 5) )
      allocate( sTka(kchunk*lnov), stat=alloc_error_sp( 6) )
      allocate( sJia(lnvv),        stat=alloc_error_sp( 7) )
      allocate( sJka(kchunk*lnvv), stat=alloc_error_sp( 8) )
      allocate( sJij(lnov*nocc),   stat=alloc_error_sp( 9) )
      allocate( sJkj(kchunk*lnov), stat=alloc_error_sp(10) )
      allocate( sKia(lnvv),        stat=alloc_error_sp(11) )
      allocate( sKka(kchunk*lnvv), stat=alloc_error_sp(12) )
      allocate( sKij(lnov*nocc),   stat=alloc_error_sp(13) )
      allocate( sKkj(kchunk*lnov), stat=alloc_error_sp(14) )
      allocate( f1n(nvir,nvir),    stat=alloc_error_sp(15) )
      allocate( f2n(nvir,nvir),    stat=alloc_error_sp(16) )
      allocate( f3n(nvir,nvir),    stat=alloc_error_sp(17) )
      allocate( f4n(nvir,nvir),    stat=alloc_error_sp(18) )
      allocate( f1t(nvir,nvir),    stat=alloc_error_sp(19) )
      allocate( f2t(nvir,nvir),    stat=alloc_error_sp(20) )
      allocate( f3t(nvir,nvir),    stat=alloc_error_sp(21) )
      allocate( f4t(nvir,nvir),    stat=alloc_error_sp(22) )
      ! see https://gcc.gnu.org/onlinedocs/gfortran/ALL.html
      if (.not.all(alloc_error_sp.eq.0)) then
          call errquit('allocate sp',1,MA_ERR)
      endif
#ifdef USE_F90_ALLOCATABLE
      allocate( dintc1(1:nvir), stat=alloc_error(1))
      allocate( dintx1(1:nvir), stat=alloc_error(2))
      allocate( t1v1(1:nvir),   stat=alloc_error(3))
      allocate( dintc2(1:nvir), stat=alloc_error(4))
      allocate( dintx2(1:nvir), stat=alloc_error(5))
      allocate( t1v2(1:nvir),   stat=alloc_error(6))
      if (.not.all(alloc_error.eq.0)) then
          call errquit('allocate {dint,t1}[cx][12]',1,MA_ERR)
      endif
#endif
!
      nodes = ga_nnodes()
      iam = ga_nodeid()
!
!      call ga_sync() ! ga_sync called just before trpdrv in aoccsd2
!
      if (occsdps) then
         call pstat_on(ps_trpdrv)
      else
         call qenter('trpdrv',0)
      endif
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
                  endif

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

                     !call dcopy(nvir,t1((i-1)*nvir+1),1,t1v2,1)
                     call copy_64_to_32(nvir, t1((i-1)*nvir+1),t1v2)
                     call ga_nbwait(nbh_objv1) ! Dja
                     !call dcopy(nvir,Dja(1+(i-1)*nvir),1,dintc1,1)
                     call copy_64_to_32(nvir,Dja(1+(i-1)*nvir),dintc1)
                     call ga_nbwait(nbh_objv5) ! Djia
                     !call dcopy(nvir,Djia,1,dintx1,1)
                     call copy_64_to_32(nvir,Djia,dintx1)

                     do k=klo,min(khi,i)

                        !call dcopy(nvir,t1((k-1)*nvir+1),1,t1v1,1)
                        call copy_64_to_32(nvir,t1((k-1)*nvir+1),
     &                                          t1v1)
                        !call dcopy(nvir,Dja(1+(k-1)*nvir),1,dintc2,1)
                        call copy_64_to_32(nvir,Dja(1+(k-1)*nvir),
     &                                          dintc2)
                        call ga_nbwait(nbh_objv4(k)) ! Djka
                        !call dcopy(nvir,Djka(1+(k-klo)*nvir),1,dintx2,1)
                        call copy_64_to_32(nvir,Djka(1+(k-klo)*nvir),
     &                                          dintx2)
                        if (occsdps) then
                           call pstat_on(ps_doxxx)
                        else
                           call qenter('doxxx',0)
                        endif
!
!  These are the input dependencies for the DGEMM calls below.
!  We wait on all of them here because GA is not even remotely thread-safe.
!  All of these are independent of k, so we wait on them only
!  at the first trip of the loop.
!
                        if (k.eq.klo) then
                            call ga_nbwait(nbh_objv2)
                            call ga_nbwait(nbh_objv3)
                            call ga_nbwait(nbh_objv6)
                            call ga_nbwait(nbh_objv7)
                            call ga_nbwait(nbh_objo1)
                            call ga_nbwait(nbh_objo2)
                            call ga_nbwait(nbh_objo3)
                            call ga_nbwait(nbh_objo4)
                            call ga_nbwait(nbh_objo5)
                            call ga_nbwait(nbh_objo6)
                            call ga_nbwait(nbh_exch1)
                            call ga_nbwait(nbh_exch2)
                            call ga_nbwait(nbh_coul1)
                            call ga_nbwait(nbh_coul2)
                        endif
!$omp parallel
!$omp& shared(eorb)
!$omp& shared(f1n,f2n,f3n,f4n,f1t,f2t,f3t,f4t)
!$omp& shared(t1v1,dintc1,dintx1)
!$omp& shared(t1v2,dintc2,dintx2)
!$omp& private(eaijk,denom)
!$omp& firstprivate(ncor,nocc,nvir,lnov,lnvv,i,j,k,klo)

           offo = 1+(k-klo)*lnov
           offv = 1+(k-klo)*lnvv
!
! Performance Note:
!
! By definition, the following does not scale to more than 8 threads
! unless nested parallelism (i.e. inside of DGEMM) is used.
! It may be prudent to write a manually threaded wrapper for the
! cases where single-threaded BLAS is used.
!
!$omp sections
!$omp section
           call copy_64_to_32(nvir*nvir, Jia, sJia)
           call copy_64_to_32(nvir*nvir, Tkj(offv), sTkj(offv))
           call sgemm('n','t',nvir,nvir,nvir,1.0,sJia,nvir,
     1                sTkj(offv),nvir,0.0,f1n,nvir)

           call copy_64_to_32(nvir*nocc, Tia, sTia)
           call copy_64_to_32(nvir*nocc, Kkj(offo), sKkj(offo))
           call sgemm('n','n',nvir,nvir,nocc,-1.0,sTia,nvir,
     1                sKkj(offo),nocc,1.0,f1n,nvir)
!$omp section
           call copy_64_to_32(nvir*nvir, Kia, sKia)
           call copy_64_to_32(nvir*nvir, Tkj(offv), sTkj(offv))
           call sgemm('n','t',nvir,nvir,nvir,1.0,sKia,nvir,
     1                sTkj(offv),nvir,0.0,f2n,nvir)

           call copy_64_to_32(nvir*nocc, Xia, sXia)
           call copy_64_to_32(nvir*nocc, Kkj(offo), sKkj(offo))
           call sgemm('n','n',nvir,nvir,nocc,-1.0,sXia,nvir,
     1                sKkj(offo),nocc,1.0,f2n,nvir)
!$omp section
           call copy_64_to_32(nvir*nvir, Jia, sJia)
           call copy_64_to_32(nvir*nvir, Tkj(offv), sTkj(offv))
           call sgemm('n','n',nvir,nvir,nvir,1.0,sJia,nvir,
     1                sTkj(offv),nvir,0.0,f3n,nvir)

           call copy_64_to_32(nvir*nocc, Tia, sTia)
           call copy_64_to_32(nvir*nocc, Jkj(offo), sJkj(offo))
           call sgemm('n','n',nvir,nvir,nocc,-1.0,sTia,nvir,
     1                sJkj(offo),nocc,1.0,f3n,nvir)
!$omp section
           call copy_64_to_32(nvir*nvir, Kia, sKia)
           call copy_64_to_32(nvir*nvir, Tkj(offv), sTkj(offv))
           call sgemm('n','n',nvir,nvir,nvir,1.0,sKia,nvir,
     1                sTkj(offv),nvir,0.0,f4n,nvir)

           call copy_64_to_32(nvir*nocc, Xia, sXia)
           call copy_64_to_32(nvir*nocc, Jkj(offo), sJkj(offo))
           call sgemm('n','n',nvir,nvir,nocc,-1.0,sXia,nvir,
     1                sJkj(offo),nocc,1.0,f4n,nvir)
!$omp section
           call copy_64_to_32(nvir*nvir, Jka(offv), sJka(offv))
           call copy_64_to_32(nvir*nvir, Tij, sTij)
           call sgemm('n','t',nvir,nvir,nvir,1.0,sJka(offv),
     1                nvir,sTij,nvir,0.0,f1t,nvir)

           call copy_64_to_32(nvir*nocc, Tka(offo), sTka(offo))
           call copy_64_to_32(nvir*nocc, Kij, sKij)
           call sgemm('n','n',nvir,nvir,nocc,-1.0,sTka(offo),
     1                nvir,sKij,nocc,1.0,f1t,nvir)
!$omp section
           call copy_64_to_32(nvir*nvir, Kka(offv), sKka(offv))
           call copy_64_to_32(nvir*nvir, Tij, sTij)
           call sgemm('n','t',nvir,nvir,nvir,1.0,sKka(offv),
     1                nvir,sTij,nvir,0.0,f2t,nvir)

           call copy_64_to_32(nvir*nocc, Xka(offo), sXka(offo))
           call copy_64_to_32(nvir*nocc, Kij, sKij)
           call sgemm('n','n',nvir,nvir,nocc,-1.0,sXka(offo),
     1                nvir,sKij,nocc,1.0,f2t,nvir)
!$omp section
           call copy_64_to_32(nvir*nvir, Jka(offv), sJka(offv))
           call copy_64_to_32(nvir*nvir, Tij, sTij)
           call sgemm('n','n',nvir,nvir,nvir,1.0,sJka(offv),
     1                nvir,sTij,nvir,0.0,f3t,nvir)

           call copy_64_to_32(nvir*nocc, Tka(offo), sTka(offo))
           call copy_64_to_32(nvir*nocc, Jij, sJij)
           call sgemm('n','n',nvir,nvir,nocc,-1.0,sTka(offo),
     1                nvir,sJij,nocc,1.0,f3t,nvir)
!$omp section
           call copy_64_to_32(nvir*nvir, Kka(offv), sKka(offv))
           call copy_64_to_32(nvir*nvir, Tij, sTij)
           call sgemm('n','n',nvir,nvir,nvir,1.0,sKka(offv),
     1                nvir,sTij,nvir,0.0,f4t,nvir)

           call copy_64_to_32(nvir*nocc, Xka(offo), sXka(offo))
           call copy_64_to_32(nvir*nocc, Jij, sJij)
           call sgemm('n','n',nvir,nvir,nocc,-1.0,sXka(offo),
     1                nvir,sJij,nocc,1.0,f4t,nvir)
!$omp end sections

!$omp master
           if (occsdps) then
              call pstat_off(ps_doxxx)
              call pstat_on(ps_tengy)
           else
              call qexit('doxxx',0)
              call qenter('tengy',0)
           endif

           emp4i = 0.
           emp5i = 0.
           emp4k = 0.
           emp5k = 0.
!$omp end master

           eaijk = eorb(a) - eorb(ncor+i) - eorb(ncor+j) - eorb(ncor+k)
!$omp do collapse(2)
!$omp& schedule(static)
!$omp& reduction(+:emp5i,emp4i)
!$omp& reduction(+:emp5k,emp4k)
! WARNING: Do not add IVDEP here.  Code will be incorrect.
           do bb=1,nvir,chunking
             do cc=1,nvir,chunking
               do b=bb,min(bb+chunking-1,nvir)
                 do c=cc,min(cc+chunking-1,nvir)
! less accurate
!                   denom=real(-1.0d0/(eorb(ncor+nocc+b)
!     &                               +eorb(ncor+nocc+c)+eaijk),kind=sp)
! more accurate
!                   denom=-1./(eorb(ncor+nocc+b)+eorb(ncor+nocc+c)+eaijk)
! most accurate
                   denom=-1.0d0/(eorb(ncor+nocc+b)
     &                          +eorb(ncor+nocc+c)+eaijk)
! fusing emp[45][ki] accumulates may help vectorization...
                   emp4i=emp4i+denom*
     &                  (f1t(b,c)+f1n(c,b)+f2t(c,b)
     &                           +f3n(b,c)+f4n(c,b))*
     &                  (f1t(b,c)-2*f2t(b,c)-2*f3t(b,c)+f4t(b,c))
     &                        -denom*
     &                  (f1n(b,c)+f1t(c,b)+f2n(c,b)+f3n(c,b))*
     &                  (2*f1t(b,c)-f2t(b,c)-f3t(b,c)+2*f4t(b,c))
     &                        +3*denom*(
     &                  f1n(b,c)*(f1n(b,c)+f3n(c,b)+2*f4t(c,b))+
     &                  f2n(b,c)*f2t(c,b)+f3n(b,c)*f4t(b,c))
                   emp4k=emp4k+denom*
     &                  (f1n(b,c)+f1t(c,b)+f2n(c,b)
     &                           +f3t(b,c)+f4t(c,b))*
     &                  (f1n(b,c)-2*f2n(b,c)-2*f3n(b,c)+f4n(b,c))
     &                        -denom*
     &                  (f1t(b,c)+f1n(c,b)+f2t(c,b)+f3t(c,b))*
     &                  (2*f1n(b,c)-f2n(b,c)-f3n(b,c)+2*f4n(b,c))
     &                        +3*denom*(
     &                  f1t(b,c)*(f1t(b,c)+f3t(c,b)+2*f4n(c,b))+
     &                  f2t(b,c)*f2n(c,b)+f3t(b,c)*f4n(b,c))
                   emp5i=emp5i+denom*t1v1(b)*dintx1(c)*
     &                 (    f1t(b,c)+f2n(b,c)+f4n(c,b)
     &                  -2*(f3t(b,c)+f4n(b,c)+f2n(c,b)+
     &                      f1n(b,c)+f2t(b,c)+f3n(c,b))
     &                  +4*(f3n(b,c)+f4t(b,c)+f1n(c,b)))
     &                        +denom*t1v1(b)*dintc1(c)*
     &                 (     f1n(b,c)+f4n(b,c)+f1t(c,b)
     &                   -2*(f2n(b,c)+f3n(b,c)+f2t(c,b)))
                   emp5k=emp5k+denom*t1v2(b)*dintx2(c)*
     &                 (    f1n(b,c)+f2t(b,c)+f4t(c,b)
     &                  -2*(f3n(b,c)+f4t(b,c)+f2t(c,b)+
     &                      f1t(b,c)+f2n(b,c)+f3t(c,b))
     &                  +4*(f3t(b,c)+f4n(b,c)+f1t(c,b)))
     &                        +denom*t1v2(b)*dintc2(c)*
     &                 (     f1t(b,c)+f4t(b,c)+f1n(c,b)
     &                   -2*(f2t(b,c)+f3t(b,c)+f2n(c,b)))
                 enddo
               enddo
             enddo
           enddo
!$omp end do
!$omp end parallel
                         if (occsdps) then
                            call pstat_off(ps_tengy)
                         else
                            call qexit('tengy',0)
                         endif

                         emp4 = emp4 + real(emp4i,kind=dp)
                         emp5 = emp5 + real(emp5i,kind=dp)
                         if (i.ne.k) then
                             emp4 = emp4 + real(emp4k,kind=dp)
                             emp5 = emp5 + real(emp5k,kind=dp)
                         end if ! (i.ne.k)
                     end do    ! k
                  end do       ! i
                  if (iprt.gt.50)then
                     write(6,1234)iam,a,j,emp4,emp5
 1234                format(' iam aijk',3i5,2e15.5)
                  end if
                  next=nxtask(nodes, 1)
               end if
            end do
            if(ga_nodeid().eq.0) then
               write(6,4321) ' ccsd(t): done ',
     A              a-(ncor+nocc)+((klo-1)/kchunk)*nvir,
     O              ' out of ',(nocc/kchunk)*nvir,
     O              ' progress: ',
     O              ((a-(ncor+nocc)+((klo-1)/kchunk)*nvir)*100d0)/
     D              ((nocc/kchunk)*nvir),
     P              '%'
 4321          format(a,i8,a,i8,a,f6.1,a1)
            endif
         end do
      end do
      call ga_sync()
      next=nxtask(-nodes, 1)
      call ga_sync()
      if (occsdps) then
         call pstat_off(ps_trpdrv)
      else
         call qexit('trpdrv',0)
      endif
!
#ifdef USE_F90_ALLOCATABLE
      deallocate( dintc1, stat=alloc_error(1))
      deallocate( dintx1, stat=alloc_error(2))
      deallocate( t1v1,   stat=alloc_error(3))
      deallocate( dintc2, stat=alloc_error(4))
      deallocate( dintx2, stat=alloc_error(5))
      deallocate( t1v2,   stat=alloc_error(6))
      if (.not.all(alloc_error.eq.0)) then
          call errquit('deallocate {dint,t1}[cx][12]',2,MA_ERR)
      endif
#endif
      deallocate( sXia, stat=alloc_error_sp( 1) )
      deallocate( sXka, stat=alloc_error_sp( 2) )
      deallocate( sTij, stat=alloc_error_sp( 3) )
      deallocate( sTkj, stat=alloc_error_sp( 4) )
      deallocate( sTia, stat=alloc_error_sp( 5) )
      deallocate( sTka, stat=alloc_error_sp( 6) )
      deallocate( sJia, stat=alloc_error_sp( 7) )
      deallocate( sJka, stat=alloc_error_sp( 8) )
      deallocate( sJij, stat=alloc_error_sp( 9) )
      deallocate( sJkj, stat=alloc_error_sp(10) )
      deallocate( sKia, stat=alloc_error_sp(11) )
      deallocate( sKka, stat=alloc_error_sp(12) )
      deallocate( sKij, stat=alloc_error_sp(13) )
      deallocate( sKkj, stat=alloc_error_sp(14) )
      deallocate( f1n,  stat=alloc_error_sp(15) )
      deallocate( f2n,  stat=alloc_error_sp(16) )
      deallocate( f3n,  stat=alloc_error_sp(17) )
      deallocate( f4n,  stat=alloc_error_sp(18) )
      deallocate( f1t,  stat=alloc_error_sp(19) )
      deallocate( f2t,  stat=alloc_error_sp(20) )
      deallocate( f3t,  stat=alloc_error_sp(21) )
      deallocate( f4t,  stat=alloc_error_sp(22) )
      if (.not.all(alloc_error_sp.eq.0)) then
          call errquit('deallocate sp',2,MA_ERR)
      endif
!
      end
#endif
