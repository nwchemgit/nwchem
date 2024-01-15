#define USE_OMP_SECTIONS 1
c#define BLAS_THREAD 1
      subroutine ccsd_trpdrv_omp(t1,
     &     f1n,f1t,f2n,f2t,f3n,f3t,f4n,f4t,eorb,
     &     g_objo,g_objv,g_coul,g_exch,
     &     ncor,nocc,nvir,iprt,emp4,emp5,
     &     oseg_lo,oseg_hi, kchunk,
     &     Tij, Tkj, Tia, Tka, Xia, Xka, Jia, Jka, Kia, Kka,
     &     Jij, Jkj, Kij, Kkj, Dja, Djka, Djia)
!     $Id$
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
      integer, intent(in) :: ncor,nocc,nvir
      integer, intent(in) :: iprt
      integer, intent(in) :: g_objo,g_objv,g_coul,g_exch
      integer, intent(in) :: oseg_lo,oseg_hi, kchunk
! FYI: f[1234][nt] parent routine allocates these on the heap (+fastmem)
      double precision, intent(in) :: f1n(nvir,nvir),f1t(nvir,nvir)
      double precision, intent(in) :: f2n(nvir,nvir),f2t(nvir,nvir)
      double precision, intent(in) :: f3n(nvir,nvir),f3t(nvir,nvir)
      double precision, intent(in) :: f4n(nvir,nvir),f4t(nvir,nvir)
! TODO: allocate a copy of eorb on heap with fastmem for KNL
      double precision, intent(in) :: eorb(*)
! TODO: determine if these should be allocated w/ fastmem in parent subroutine
      double precision, intent(in) :: Tij(*), Tkj(*), Tia(*), Tka(*)
      double precision, intent(in) :: Xia(*), Xka(*)
      double precision, intent(in) :: Jia(*), Jka(*), Jij(*), Jkj(*)
      double precision, intent(in) :: Kia(*), Kka(*), Kij(*), Kkj(*)
      double precision, intent(in) :: Dja(*), Djka(*), Djia(*)
! used to make inline threaded tengy correct - for now
! it is correct that dint[cx]1 are paired with t1v2 and vice versa
! in the inlined tengy loops.  see ccsd_tengy in ccsd_trpdrv.F for
! verification of the i-k and k-i pairing of these.
#ifdef USE_F90_ALLOCATABLE
      double precision, allocatable :: dintc1(:),dintc2(:)
      double precision, allocatable :: dintx1(:),dintx2(:)
      double precision, allocatable :: t1v1(:),t1v2(:)
#ifdef USE_FASTMEM
!dec$ attributes fastmem :: dintc1,dintc2,dintx1,dintx2,t1v1,t1v2
#endif
      integer :: alloc_error
#else
      double precision :: dintc1(nvir),dintx1(nvir),t1v1(nvir)
      double precision :: dintc2(nvir),dintx2(nvir),t1v2(nvir)
#endif
!
      double precision :: emp4i,emp5i,emp4k,emp5k
      double precision :: eaijk,denom
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
      integer nxtask
      external nxtask
      double precision perfm_flop,tzero,flopzero,t_flops,agg_flops
      external perfm_flop

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
      logical ldyn_org
      logical got_ak
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
      tzero=util_wallsec()
      flopzero=perfm_flop()
#ifdef USE_F90_ALLOCATABLE
      allocate( dintc1(1:nvir), dintc2(1:nvir),
     &          dintx1(1:nvir), dintx2(1:nvir),
     &          t1v1(1:nvir), t1v2(1:nvir), stat=alloc_error)
      if (alloc_error.ne.0) call errquit('CXT1 temps',6,MA_ERR)
#endif
!
#ifdef BLAS_THREAD
      call omp_set_nested(1)
      call omp_get_dynamic(ldyn_org)
      call omp_set_dynamic(.true.)
      call  util_setthreads_fromenv()
c      call util_blas_set_num_threads(4)
#endif
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

                     call dcopy(nvir,t1((i-1)*nvir+1),1,t1v2,1)
                     if(i.eq.1) then
                        call ga_nbwait(nbh_objv1) ! Dja
                     endif
                     call dcopy(nvir,Dja(1+(i-1)*nvir),1,dintc1,1)
                     call ga_nbwait(nbh_objv5) ! Djia
                     call dcopy(nvir,Djia,1,dintx1,1)

                     do k=klo,min(khi,i)

                        call dcopy(nvir,t1((k-1)*nvir+1),1,t1v1,1)
                        call dcopy(nvir,Dja(1+(k-1)*nvir),1,dintc2,1)
                        if(i.eq.1) then
                           call ga_nbwait(nbh_objv4(k)) ! Djka
                        endif
                        call dcopy(nvir,Djka(1+(k-klo)*nvir),1,dintx2,1)
                        emp4i = 0.0d0
                        emp5i = 0.0d0
                        emp4k = 0.0d0
                        emp5k = 0.0d0
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
                           if (i.eq.1) then
                              if (got_ak) then
                                 call ga_nbwait(nbh_exch1)
                                 call ga_nbwait(nbh_coul1)
                                 call ga_nbwait(nbh_objv2)
                                 call ga_nbwait(nbh_objv3)
                              endif ! got_ak
                              call ga_nbwait(nbh_objo1)
                              call ga_nbwait(nbh_objo2)
                              call ga_nbwait(nbh_objo3)
                           endif ! i==1
                           call ga_nbwait(nbh_objo4)
                           call ga_nbwait(nbh_objo5)
                           call ga_nbwait(nbh_objo6)
                           call ga_nbwait(nbh_coul2)
                           call ga_nbwait(nbh_exch2)
                           call ga_nbwait(nbh_objv6)
                           call ga_nbwait(nbh_objv7)
                        endif ! k==klo

#if USE_OMP_SECTIONS
!$omp parallel
!$omp& shared(eorb)
!$omp& shared(f1n,f2n,f3n,f4n,f1t,f2t,f3t,f4t)
!$omp& shared(t1v1,dintc1,dintx1)
!$omp& shared(t1v2,dintc2,dintx2)
!$omp& private(eaijk,denom)
!$omp& firstprivate(ncor,nocc,nvir,lnov,lnvv,i,j,k,klo)
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
#endif
                        call dgemm('n','t',nvir,nvir,nvir,1.0d0,
     1                       Jia,nvir,Tkj(1+(k-klo)*lnvv),nvir,0.0d0,
     2                       f1n,nvir)
                        call dgemm('n','n',nvir,nvir,nocc,-1.0d0,
     1                       Tia,nvir,Kkj(1+(k-klo)*lnov),nocc,1.0d0,
     2                       f1n,nvir)
#if USE_OMP_SECTIONS
!$omp section
#endif
                        call dgemm('n','t',nvir,nvir,nvir,1.0d0,
     1                       Kia,nvir,Tkj(1+(k-klo)*lnvv),nvir,0.0d0,
     2                       f2n,nvir)
                        call dgemm('n','n',nvir,nvir,nocc,-1.0d0,
     1                       Xia,nvir,Kkj(1+(k-klo)*lnov),nocc,1.0d0,
     2                       f2n,nvir)
#if USE_OMP_SECTIONS
!$omp section
#endif
                        call dgemm('n','n',nvir,nvir,nvir,1.0d0,
     1                       Jia,nvir,Tkj(1+(k-klo)*lnvv),nvir,0.0d0,
     2                       f3n,nvir)
                        call dgemm('n','n',nvir,nvir,nocc,-1.0d0,
     1                       Tia,nvir,Jkj(1+(k-klo)*lnov),nocc,1.0d0,
     2                       f3n,nvir)
#if USE_OMP_SECTIONS
!$omp section
#endif
                        call dgemm('n','n',nvir,nvir,nvir,1.0d0,
     1                       Kia,nvir,Tkj(1+(k-klo)*lnvv),nvir,0.0d0,
     2                       f4n,nvir)
                        call dgemm('n','n',nvir,nvir,nocc,-1.0d0,
     1                       Xia,nvir,Jkj(1+(k-klo)*lnov),nocc,1.0d0,
     2                       f4n,nvir)
#if USE_OMP_SECTIONS
!$omp section
#endif
                        call dgemm('n','t',nvir,nvir,nvir,1.0d0,
     1                       Jka(1+(k-klo)*lnvv),nvir,Tij,nvir,0.0d0,
     2                       f1t,nvir)
                        call dgemm('n','n',nvir,nvir,nocc,-1.0d0,
     1                       Tka(1+(k-klo)*lnov),nvir,Kij,nocc,1.0d0,
     2                       f1t,nvir)
#if USE_OMP_SECTIONS
!$omp section
#endif
                        call dgemm('n','t',nvir,nvir,nvir,1.0d0,
     1                       Kka(1+(k-klo)*lnvv),nvir,Tij,nvir,0.0d0,
     2                       f2t,nvir)
                        call dgemm('n','n',nvir,nvir,nocc,-1.0d0,
     1                       Xka(1+(k-klo)*lnov),nvir,Kij,nocc,1.0d0,
     2                       f2t,nvir)
#if USE_OMP_SECTIONS
!$omp section
#endif
                        call dgemm('n','n',nvir,nvir,nvir,1.0d0,
     1                       Jka(1+(k-klo)*lnvv),nvir,Tij,nvir,0.0d0,
     2                       f3t,nvir)
                        call dgemm('n','n',nvir,nvir,nocc,-1.0d0,
     1                       Tka(1+(k-klo)*lnov),nvir,Jij,nocc,1.0d0,
     2                       f3t,nvir)
#if USE_OMP_SECTIONS
!$omp section
#endif
                        call dgemm('n','n',nvir,nvir,nvir,1.0d0,
     1                       Kka(1+(k-klo)*lnvv),nvir,Tij,nvir,0.0d0,
     2                       f4t,nvir)
                        call dgemm('n','n',nvir,nvir,nocc,-1.0d0,
     1                       Xka(1+(k-klo)*lnov),nvir,Jij,nocc,1.0d0,
     2                       f4t,nvir)
#if USE_OMP_SECTIONS
!$omp end sections
!$omp master
#endif
                        if (occsdps) then
                           call pstat_off(ps_doxxx)
                           call pstat_on(ps_tengy)
                        else
                           call qexit('doxxx',0)
                           call qenter('tengy',0)
                        endif
#if USE_OMP_SECTIONS
!$omp end master
#endif

                        eaijk=eorb(a) - (  eorb(ncor+i)
     &                                    +eorb(ncor+j)
     &                                    +eorb(ncor+k) )
#ifdef USE_YFLOP
      flops_ycount = flops_ycount + nvir*nvir*(
     D                       3 + 2*(
     E                       12 +
     E                       11 +
     E                       11 ) +
     5                       2*27 )
#endif

#if USE_OMP_SECTIONS
!$omp do collapse(2)
#else
!$omp parallel do collapse(2)
!$omp& private(denom)
!$omp& firstprivate(ncor,nocc,nvir,eaijk)
#endif
!$omp& schedule(static)
!$omp& reduction(+:emp5i,emp4i)
!$omp& reduction(+:emp5k,emp4k)
           do bb=1,nvir,chunking
             do cc=1,nvir,chunking
               do b=bb,min(bb+chunking-1,nvir)
                 do c=cc,min(cc+chunking-1,nvir)
                   denom=-1.0d0/( eorb(ncor+nocc+b)
     &                           +eorb(ncor+nocc+c)+eaijk )
! fusing emp[45][ki] accumulates may help vectorization...
                   emp4i=emp4i+denom*
     &                  (f1t(b,c)+f1n(c,b)+f2t(c,b)+f3n(b,c)+f4n(c,b))*
     &                  (f1t(b,c)-2*f2t(b,c)-2*f3t(b,c)+f4t(b,c))
     &                        -denom*
     &                  (f1n(b,c)+f1t(c,b)+f2n(c,b)+f3n(c,b))*
     &                  (2*f1t(b,c)-f2t(b,c)-f3t(b,c)+2*f4t(b,c))
     &                        +3*denom*(
     &                  f1n(b,c)*(f1n(b,c)+f3n(c,b)+2*f4t(c,b))+
     &                  f2n(b,c)*f2t(c,b)+f3n(b,c)*f4t(b,c))
                   emp4k=emp4k+denom*
     &                  (f1n(b,c)+f1t(c,b)+f2n(c,b)+f3t(b,c)+f4t(c,b))*
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
#if USE_OMP_SECTIONS
!$omp end do
!$omp end parallel
#endif
                         if (occsdps) then
                            call pstat_off(ps_tengy)
                         else
                            call qexit('tengy',0)
                         endif

                         emp4 = emp4 + emp4i
                         emp5 = emp5 + emp5i
                         if (i.ne.k) then
                             emp4 = emp4 + emp4k
                             emp5 = emp5 + emp5k
                         end if ! (i.ne.k)
                     end do    ! k
                  end do       ! i
                  if (iprt.gt.50)then
                     write(6,1234)iam,a,j,emp4,emp5
 1234                format(' iam aijk',3i5,2e15.5)
                  end if
                  next=nxtask(nodes, 1)
            if(ga_nodeid().eq.0) then
               pct_progr=(a-(ncor+nocc)+((klo-1)/kchunk)*nvir)*n_progr/
     /          ((nocc/kchunk)*nvir)+1
               if(i_progr(pct_progr)) then
                  i_progr(pct_progr)=.false.
               write(6,4321) ' ccsd(t): done ',
     A              a-(ncor+nocc)+((klo-1)/kchunk)*nvir,
     O              ' out of ',(nocc/kchunk)*nvir,
     O              ' progress: ',
     O              ((a-(ncor+nocc)+((klo-1)/kchunk)*nvir)*100)/
     D              ((nocc/kchunk)*nvir),
     P            '%, Gflops=',(perfm_flop()-flopzero)/
     D              (util_wallsec()-tzero),
     P                 ' at ',(util_wallsec()-tzero),' secs'
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
      t_flops=util_wallsec()-tzero
      agg_flops=perfm_flop()-flopzero
      call ga_dgop(msg_cc_diis1,agg_flops,1, '+')
      if(ga_nodeid().eq.0) then
         write(6,4322) ' ccsd(t): 100% done, Aggregate Gflops=',
     P        agg_flops/t_flops,
     P                 ' in ',t_flops,' secs'
 4322    format(a,1pg11.4,a,0pf10.1,a)
         call util_flush(6)
      endif
      call ga_sync()
      if (occsdps) then
         call pstat_off(ps_trpdrv)
      else
         call qexit('trpdrv',0)
      endif
#ifdef BLAS_THREAD
      call omp_set_nested(0)
c     call mkl_set_dynamic(0)
      call omp_set_dynamic(ldyn_org)
      call util_blas_set_num_threads(1)
#endif
!
#ifdef USE_F90_ALLOCATABLE
      deallocate( dintc1, dintx1, t1v1, dintc2, dintx2, t1v2,
     &            stat=alloc_error)
      if (alloc_error.ne.0) call errquit('free CXT1 temps',6,MA_ERR)
#endif
!
      end
