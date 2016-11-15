#if defined(USE_OPENMP)
      subroutine ccsd_idx2(snsi,sisn,cmo,scra,scrb,lscr,
     &                     nfi,nfk,ncor,nocc,nact,nbf,
     &                     ilo,ihi,klo,khi,offik,offki,
     &                     g_x,g_c)
      implicit none
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
      ! allocated externally, challenge for thread-safety
      ! unused in the OpenMP version of this function
      double precision, intent(inout) :: scra(lscr),scrb(lscr)
      ! stack buffers, may segfault if large
      ! allocating these hits IFRTL global lock issue
      ! allocate num_max_thrd*nbf*nbf array outside and
      ! have TLS pointer to chunk for each thread
      double precision :: bufa(lscr),bufb(lscr)
      integer nfik,ii,jj,ik,ad1,ad2
      if (occsdps) then
         call pstat_on(ps_idx2)
      else
         call qenter('idx2',0)
      endif
      nfik=nfi*nfk
!$omp parallel do schedule(guided)
!$omp& default(none)
!$omp& firstprivate(nfik,nfi,nfk,offik,offki,ilo,klo)
!$omp& private(ii,jj,ik,ad1,ad2,bufa,bufb)
!$omp& shared(g_x,g_c,cmo,snsi,sisn,ncor,nocc,nact,nbf,lscr)
      do ii=1,nocc
        call dgemm('t','n',nact,nfik,nbf,1.0d0,cmo(1,ncor+1),nbf,
     &             snsi(1,1,1,ii),nbf,0.0d0,bufa,nact)
        ! matrix transpose
        do jj=1,nocc
          do ik=1,nfik
            bufb((jj-1)*nfik+ik)=bufa((ik-1)*nact+jj)
          enddo
        enddo
        ad1=offki
        ad2=(ii-1)*nocc
!$omp critical
        call ga_put(g_x,ad1+1,ad1+nfik,ad2+1,ad2+nocc,bufb,nfik)
!$omp end critical
        call dgemm('n','n',nact*nfi,nocc,nfk,1.0d0,bufa,nact*nfi,
     &               cmo(klo,ncor+1),nbf,0.0d0,bufb,nact*nfi)
        ad1=(ilo-1)*nact
        ad2=(ii-1)*nocc
!$omp critical
        call ga_acc(g_c,ad1+1,ad1+nact*nfi,ad2+1,ad2+nocc,bufb,
     &              nact*nfi,1.0d0)
!$omp end critical
        if (ilo.ne.klo)then
          call dgemm('t','n',nact,nfik,nbf,1.0d0,cmo(1,ncor+1),nbf,
     &               sisn(1,1,1,ii),nbf,0.0d0,bufa,nact)
          ! matrix transpose
          do jj=1,nocc
            do ik=1,nfik
              bufb((jj-1)*nfik+ik)=bufa((ik-1)*nact+jj)
            enddo
          enddo
          ad1=offik
          ad2=(ii-1)*nocc
!$omp critical
          call ga_put(g_x,ad1+1,ad1+nfik,ad2+1,ad2+nocc,bufb,nfik)
!$omp end critical
          call dgemm('n','n',nact*nfk,nocc,nfi,1.0d0,bufa,nact*nfk,
     &                 cmo(ilo,ncor+1),nbf,0.0d0,bufb,nact*nfk)
          ad1=(klo-1)*nact
          ad2=(ii-1)*nocc
!$omp critical
          call ga_acc(g_c,ad1+1,ad1+nact*nfk,ad2+1,ad2+nocc,bufb,
     &                nact*nfk,1.0d0)
!$omp end critical
        endif
      enddo
!$omp end parallel do
      if (occsdps) then
         call pstat_off(ps_idx2)
      else
         call qexit('idx2',0)
      endif
      return
      end     
#else
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
        call ga_acc(g_c,ad1+1,ad1+nact*nfi,ad2+1,ad2+nocc,scrb,
     &              nact*nfi,1.0d0)
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
          call ga_acc(g_c,ad1+1,ad1+nact*nfk,ad2+1,ad2+nocc,scrb,
     &                nact*nfk,1.0d0)
        endif
      enddo
      if (occsdps) then
         call pstat_off(ps_idx2)
      else
         call qexit('idx2',0)
      endif
      return
      end     
#endif
