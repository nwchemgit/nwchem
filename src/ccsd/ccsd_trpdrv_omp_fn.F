      subroutine ccsd_trpdrv_omp_fn(t1,
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
      integer :: inode,next,nodes,iam
      integer :: a,i,j,k,akold,av
      integer :: klo, khi
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
   99 format(2x,'Using ',i2,' OpenMP threads in CCSD(T) - C version')
#else
      if (ga_nodeid().eq.0) then
        write(6,99) 1
        write(6,999)
      endif
   99 format(2x,'Using ',i2,' OpenMP thread in CCSD(T)')
  999 format(2x,'Recompile w/ USE_OPENMP=1 to use threads in CCSD(T)')
#endif

#ifdef USE_F90_ALLOCATABLE
      allocate( dintc1(1:nvir), stat=alloc_error)
      if (alloc_error.ne.0) call errquit('dintc1',1,MA_ERR)
      allocate( dintx1(1:nvir), stat=alloc_error)
      if (alloc_error.ne.0) call errquit('dintx1',2,MA_ERR)
      allocate( t1v1(1:nvir), stat=alloc_error)
      if (alloc_error.ne.0) call errquit('t1v1',3,MA_ERR)
      allocate( dintc2(1:nvir), stat=alloc_error)
      if (alloc_error.ne.0) call errquit('dintc2',4,MA_ERR)
      allocate( dintx2(1:nvir), stat=alloc_error)
      if (alloc_error.ne.0) call errquit('dintx2',5,MA_ERR)
      allocate( t1v2(1:nvir), stat=alloc_error)
      if (alloc_error.ne.0) call errquit('t1v2',6,MA_ERR)
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

                     call dcopy(nvir,t1((i-1)*nvir+1),1,t1v2,1)
                     call ga_nbwait(nbh_objv1) ! Dja
                     call dcopy(nvir,Dja(1+(i-1)*nvir),1,dintc1,1)
                     call ga_nbwait(nbh_objv5) ! Djia
                     call dcopy(nvir,Djia,1,dintx1,1)

                     do k=klo,min(khi,i)

                        call dcopy(nvir,t1((k-1)*nvir+1),1,t1v1,1)
                        call dcopy(nvir,Dja(1+(k-1)*nvir),1,dintc2,1)
                        call ga_nbwait(nbh_objv4(k)) ! Djka
                        call dcopy(nvir,Djka(1+(k-klo)*nvir),1,dintx2,1)
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
                        call ccsd_trpdrv_omp_body(
     &                       f1n,f1t,f2n,f2t,f3n,f3t,f4n,f4t,eorb,
     &                       ncor,nocc,nvir, emp4,emp5,
     &                       a,i,j,k,klo,
     &                       Tij, Tkj, Tia, Tka, Xia, Xka,
     &                       Jia, Jka, Kia, Kka,
     &                       Jij, Jkj, Kij, Kkj,
     &                       dintc1,dintx1,t1v1,
     &                       dintc2,dintx2,t1v2)
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
      deallocate( dintc1, stat=alloc_error)
      if (alloc_error.ne.0) call errquit('dintc1',11,MA_ERR)
      deallocate( dintx1, stat=alloc_error)
      if (alloc_error.ne.0) call errquit('dintx1',12,MA_ERR)
      deallocate( t1v1, stat=alloc_error)
      if (alloc_error.ne.0) call errquit('t1v1',13,MA_ERR)
      deallocate( dintc2, stat=alloc_error)
      if (alloc_error.ne.0) call errquit('dintc2',14,MA_ERR)
      deallocate( dintx2, stat=alloc_error)
      if (alloc_error.ne.0) call errquit('dintx2',15,MA_ERR)
      deallocate( t1v2, stat=alloc_error)
      if (alloc_error.ne.0) call errquit('t1v2',16,MA_ERR)
#endif
!
      end
