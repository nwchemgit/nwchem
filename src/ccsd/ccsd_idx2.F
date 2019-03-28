      subroutine ccsd_idx2_omp(snsi,sisn,cmo,lscr,
     &                         nfi,nfk,ncor,nocc,nact,nbf,
     &                         ilo,ihi,klo,khi,offik,offki,
     &                         g_x,g_c)
      implicit none
#include "errquit.fh"
#include "ccsdps.fh"
      ! lscr=max(maxbfsh*nocc*nbf,maxbfsh*maxbfsh*nbf)
      integer, intent(in) :: lscr
      integer, intent(in) :: nfi,nfk,ncor,nocc,nact,nbf
      integer, intent(in) :: ilo,klo
      integer, intent(in) :: ihi,khi ! unused
      integer, intent(in) :: offik,offki
      integer, intent(in) :: g_x,g_c
      double precision, intent(in) :: snsi(nbf,nfi,nfk,nocc)
      double precision, intent(in) :: sisn(nbf,nfk,nfi,nocc)
      double precision, intent(in) :: cmo(nbf,nbf)
      ! if stack buffers, may segfault if large.
      ! allocating these hits IFRTL global lock issue.
      ! allocate num_max_thrd*nbf*nbf array outside and
      ! index to private portion using thread id.
      double precision, allocatable :: bufa(:,:)
      double precision, allocatable :: bufb(:,:)
      integer :: alloc_error
      !
      integer :: nfik,ii,jj,ik,ad1,ad2
      !
      integer :: nt, tid
#if defined(USE_OPENMP)
      integer  :: omp_get_thread_num
      external :: omp_get_thread_num
      integer  :: omp_get_num_threads
      external :: omp_get_num_threads
      integer  :: omp_get_max_threads
      external :: omp_get_max_threads
#endif
      !
      if (occsdps) then
         call pstat_on(ps_idx2)
      else
         call qenter('idx2',0)
      endif
      !
      nfik=nfi*nfk
      !
#if defined(USE_OPENMP)
      nt = omp_get_max_threads()
#else
      nt = 1
#endif
      !
      allocate( bufa(lscr,nt), stat=alloc_error)
      if (alloc_error.ne.0) call errquit('ccsd_idx2_omp: bufa',1,MA_ERR)
      allocate( bufb(lscr,nt), stat=alloc_error)
      if (alloc_error.ne.0) call errquit('ccsd_idx2_omp: bufb',2,MA_ERR)
      !
!$omp parallel
!$omp& default(none)
!$omp& private(ii,jj,ik,ad1,ad2,tid)
!$omp& shared(nfik,nfi,nfk,offik,offki,ilo,klo)
!$omp& shared(bufa,bufb)
!$omp& shared(g_x,g_c,cmo,snsi,sisn,ncor,nocc,nact,nbf,lscr)
#if defined(USE_OPENMP)
      ! OpenMP thread ids run from 0 to nt-1
      tid = omp_get_thread_num() + 1
#else
      tid = 1
#endif
!$omp do schedule(guided)
      do ii=1,nocc
        call dgemm('t','n',nact,nfik,nbf,1.0d0,cmo(1,ncor+1),nbf,
     &             snsi(1,1,1,ii),nbf,0.0d0,bufa(1,tid),nact)
        ! matrix transpose
        do jj=1,nocc
          do ik=1,nfik
            bufb((jj-1)*nfik+ik,tid)=bufa((ik-1)*nact+jj,tid)
          enddo
        enddo
        ad1=offki
        ad2=(ii-1)*nocc
!$omp critical
        call ga_put(g_x,ad1+1,ad1+nfik,ad2+1,ad2+nocc,bufb(1,tid),nfik)
!$omp end critical
        call dgemm('n','n',nact*nfi,nocc,nfk,1.0d0,bufa(1,tid),nact*nfi,
     &               cmo(klo,ncor+1),nbf,0.0d0,bufb(1,tid),nact*nfi)
        ad1=(ilo-1)*nact
        ad2=(ii-1)*nocc
!$omp critical
#ifdef ACC_STRIPS
        do jj=1,nfi
           call ga_acc(g_c,ad1+1+(jj-1)*nact,ad1+nact*jj,ad2+1,ad2+nocc,
     S          bufb((jj-1)*nact+1,tid),
     &          nact*nfi,1.0d0)
        enddo
#else
        call ga_acc(g_c,ad1+1,ad1+nact*nfi,ad2+1,ad2+nocc,bufb(1,tid),
     &              nact*nfi,1.0d0)
#endif
!$omp end critical
        if (ilo.ne.klo)then
          call dgemm('t','n',nact,nfik,nbf,1.0d0,cmo(1,ncor+1),nbf,
     &               sisn(1,1,1,ii),nbf,0.0d0,bufa(1,tid),nact)
          ! matrix transpose
          do jj=1,nocc
            do ik=1,nfik
              bufb((jj-1)*nfik+ik,tid)=bufa((ik-1)*nact+jj,tid)
            enddo
          enddo
          ad1=offik
          ad2=(ii-1)*nocc
!$omp critical
          call ga_put(g_x,ad1+1,ad1+nfik,ad2+1,ad2+nocc,
     &                bufb(1,tid),nfik)
!$omp end critical
          call dgemm('n','n',nact*nfk,nocc,nfi,1.0d0,bufa(1,tid),
     &               nact*nfk,cmo(ilo,ncor+1),nbf,0.0d0,bufb(1,tid),
     &               nact*nfk)
          ad1=(klo-1)*nact
          ad2=(ii-1)*nocc
!$omp critical
#ifdef ACC_STRIPS
        do jj=1,nfk
           call ga_acc(g_c,ad1+1+(jj-1)*nact,ad1+nact*jj,ad2+1,ad2+nocc,
     S          bufb((jj-1)*nact+1,tid),
     &          nact*nfk,1.0d0)
        enddo
#else
          call ga_acc(g_c,ad1+1,ad1+nact*nfk,ad2+1,ad2+nocc,bufb(1,tid),
     &                nact*nfk,1.0d0)
#endif
!$omp end critical
        endif
      enddo
!$omp end do
!$omp end parallel
      !
      deallocate( bufa, stat=alloc_error)
      if (alloc_error.ne.0) call errquit('ccsd_idx2_omp: bufa',3,MA_ERR)
      deallocate( bufb, stat=alloc_error)
      if (alloc_error.ne.0) call errquit('ccsd_idx2_omp: bufb',4,MA_ERR)
      !
      if (occsdps) then
         call pstat_off(ps_idx2)
      else
         call qexit('idx2',0)
      endif
      return
      end
      subroutine ccsd_idx2(snsi,sisn,cmo,scra,scrb,lscr,
     &                     nfi,nfk,ncor,nocc,nact,nbf,
     &                     ilo,ihi,klo,khi,offik,offki,
     &                     g_x,g_c)
      implicit none
#include "ccsdps.fh"
      integer lscr,nfi,nfk,ncor,nocc,nact,nbf,ilo,ihi,klo,khi,
     &        offik,offki,g_x,g_c
      double precision snsi(nbf,nfi,nfk,nocc),
     &                 sisn(nbf,nfk,nfi,nocc),
     &                 cmo(nbf,nbf),
     &                 scra(lscr),scrb(lscr)
! lscr=max(maxbfsh*nocc*nbf,maxbfsh*maxbfsh*nbf)
      integer nfik,ii,jj,ik,ad1,ad2
      if (occsdps) then
         call pstat_on(ps_idx2)
      else
         call qenter('idx2',0)
      endif
      nfik=nfi*nfk
      do ii=1,nocc
        call dgemm('t','n',nact,nfik,nbf,1.0d0,cmo(1,ncor+1),nbf,
     &             snsi(1,1,1,ii),nbf,0.0d0,scra,nact)
        ! matrix transpose
        do jj=1,nocc
          do ik=1,nfik
            scrb((jj-1)*nfik+ik)=scra((ik-1)*nact+jj)
          enddo
        enddo
        ad1=offki
        ad2=(ii-1)*nocc
        call ga_put(g_x,ad1+1,ad1+nfik,ad2+1,ad2+nocc,scrb,nfik)
        call dgemm('n','n',nact*nfi,nocc,nfk,1.0d0,scra,nact*nfi,
     &               cmo(klo,ncor+1),nbf,0.0d0,scrb,nact*nfi)
        ad1=(ilo-1)*nact
        ad2=(ii-1)*nocc
#ifdef ACC_STRIPS
        call util_gacc_strips(g_c,ad1+1,ad1+nact*nfi,ad2+1,ad2+nocc,
     c       scrb,nact*nfi,1d0,5*nact)
#else
        call ga_acc(g_c,ad1+1,ad1+nact*nfi,ad2+1,ad2+nocc,scrb,
     &              nact*nfi,1.0d0)
#endif
        if (ilo.ne.klo)then
          call dgemm('t','n',nact,nfik,nbf,1.0d0,cmo(1,ncor+1),nbf,
     &               sisn(1,1,1,ii),nbf,0.0d0,scra,nact)
          ! matrix transpose
          do jj=1,nocc
            do ik=1,nfik
              scrb((jj-1)*nfik+ik)=scra((ik-1)*nact+jj)
            enddo
          enddo
          ad1=offik
          ad2=(ii-1)*nocc
          call ga_put(g_x,ad1+1,ad1+nfik,ad2+1,ad2+nocc,scrb,nfik)
          call dgemm('n','n',nact*nfk,nocc,nfi,1.0d0,scra,nact*nfk,
     &                 cmo(ilo,ncor+1),nbf,0.0d0,scrb,nact*nfk)
          ad1=(klo-1)*nact
          ad2=(ii-1)*nocc
#ifdef ACC_STRIPS
        call util_gacc_strips(g_c,ad1+1,ad1+nact*nfk,ad2+1,ad2+nocc,
     c       scrb,nfk*nact,1d0,5*nact)
#else
          call ga_acc(g_c,ad1+1,ad1+nact*nfk,ad2+1,ad2+nocc,scrb,
     &                nact*nfk,1.0d0)
#endif
        endif
      enddo
      if (occsdps) then
         call pstat_off(ps_idx2)
      else
         call qexit('idx2',0)
      endif
      return
      end     
