      subroutine ccsd_idx34_omp(basis,cmo,
     &                          nsh,ncor,nocc,nact,nbf,
     &                          g_x,g_c)
      implicit none
#include "errquit.fh"
      integer, intent(in) ::  basis,nsh,ncor,nocc,nact,nbf,g_x,g_c
      double precision, intent(in) :: cmo(nbf,nbf)
      ! if stack buffers, may segfault if large.
      ! allocating these hits IFRTL global lock issue.
      ! allocate num_max_thrd*nbf*nbf array outside and
      ! index to private portion using thread id.
      double precision, allocatable :: bufa(:,:)
      double precision, allocatable :: bufb(:,:)
      integer :: alloc_error
      !
#include "bas.fh"
#include "ccsdps.fh"
      integer :: ii,jj,ish,jsh,i,j,ilo,ihi,jlo,jhi,iam,ga_nodeid
      integer :: g_jlo,g_jhi,g_ilo,g_ihi,ad1,ad2
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
      if (occsdps) then
         call pstat_on(ps_idx34)
      else
         call qenter('idx34',0)
      endif
      !
#if defined(USE_OPENMP)
      nt = omp_get_max_threads()
#else
      nt = 1
#endif
      !
      allocate( bufa(nbf*nbf,nt), stat=alloc_error)
      if (alloc_error.ne.0)
     &       call errquit('ccsd_idx34_omp: bufa',1,MA_ERR)
      allocate( bufb(nbf*nbf,nt), stat=alloc_error)
      if (alloc_error.ne.0)
     &       call errquit('ccsd_idx34_omp: bufb',2,MA_ERR)
      !
      iam=ga_nodeid()
      call ga_distribution(g_x,iam,g_jlo,g_jhi,g_ilo,g_ihi)
!$omp parallel
!$omp& default(none)
!$omp& private(ilo,ihi,jlo,jhi,ad1,ad2,tid)
!$omp& shared(bufa,bufb)
!$omp& shared(g_x,g_c,cmo,nbf,nsh,ncor,nact,nocc,basis)
!$omp& shared(g_ilo,g_ihi)
#if defined(USE_OPENMP)
      ! OpenMP thread ids run from 0 to nt-1
      tid = omp_get_thread_num() + 1
#else
      tid = 1
#endif
!$omp do collapse(2) schedule(guided)
      do ii=1,nocc
        do jj=1,nocc
          ad2=(ii-1)*nocc+jj
          if (ad2.ge.g_ilo.and.ad2.le.g_ihi)then
!$omp critical
            call ga_get(g_x,1,nbf*nbf,ad2,ad2,bufa(1,tid),nbf*nbf)
!$omp end critical
            ! bas_cn2bfr is thread-safe (accessor function)
            ad1=0
            do ish=1,nsh
              if (.not. bas_cn2bfr(basis,ish,ilo,ihi))
     &            call errquit('vvvv: bas_cn2bfr',ish, BASIS_ERR)
              do jsh=1,nsh
                if (.not. bas_cn2bfr(basis,jsh,jlo,jhi))
     &              call errquit('vvvv: bas_cn2bfr',jsh, BASIS_ERR)
                ! matrix transpose
                do i=ilo,ihi 
                  do j=jlo,jhi 
                    ad1=ad1+1
                    bufb((i-1)*nbf+j,tid)=bufa(ad1,tid)
                  enddo
                enddo
              enddo
            enddo
            call dgemm('t','n',nact,nbf,nbf,1.0d00,cmo(1,ncor+1),nbf,
     &                 bufb(1,tid),nbf,0.0d00,bufa(1,tid),nact)
            call dgemm('n','n',nact,nact,nbf,1.0d00,bufa(1,tid),nact,
     &                 cmo(1,ncor+1),nbf,0.0d00,bufb(1,tid),nact)
!$omp critical
            call ga_put(g_x,1,nact*nact,ad2,ad2,bufb(1,tid),nact*nact)
            call ga_get(g_c,1,nact*nbf,ad2,ad2,bufa(1,tid),nact*nbf)
!$omp end critical
            call dgemm('n','n',nact,nact,nbf,1.0d00,bufa(1,tid),nact,
     &                 cmo(1,ncor+1),nbf,0.0d00,bufb(1,tid),nact)
!$omp critical
            call ga_put(g_c,1,nact*nact,ad2,ad2,bufb(1,tid),nact*nact)
!$omp end critical
          endif
        enddo
      enddo
!$omp end do
!$omp end parallel
      !
      deallocate( bufa, stat=alloc_error)
      if (alloc_error.ne.0) call errquit('ccsd_idx34: bufa',3,MA_ERR)
      deallocate( bufb, stat=alloc_error)
      if (alloc_error.ne.0) call errquit('ccsd_idx34: bufb',4,MA_ERR)
      !
      if (occsdps) then
         call pstat_off(ps_idx34)
      else
         call qexit('idx34',0)
      endif
      return
      end     
      subroutine ccsd_idx34(basis,cmo,scra,scrb,nsh,
     &                      ncor,nocc,nact,nbf,
     &                      g_x,g_c)
      implicit none
#include "errquit.fh"
      integer, intent(in) ::  basis,nsh,ncor,nocc,nact,nbf,g_x,g_c
      double precision, intent(in) :: cmo(nbf,nbf)
      double precision, intent(inout) :: scra(nbf*nbf),scrb(nbf*nbf)
#include "bas.fh"
#include "ccsdps.fh"
      integer :: ii,jj,ish,jsh,i,j,ilo,ihi,jlo,jhi,iam,ga_nodeid
      integer :: g_jlo,g_jhi,g_ilo,g_ihi,ad1,ad2
      if (occsdps) then
         call pstat_on(ps_idx34)
      else
         call qenter('idx34',0)
      endif
      iam=ga_nodeid()
      call ga_distribution(g_x,iam,g_jlo,g_jhi,g_ilo,g_ihi)
      do ii=1,nocc
        do jj=1,nocc
          ad2=(ii-1)*nocc+jj
          if (ad2.ge.g_ilo.and.ad2.le.g_ihi)then
            call ga_get(g_x,1,nbf*nbf,ad2,ad2,scra,nbf*nbf)
            ad1=0
            do ish=1,nsh
              if (.not. bas_cn2bfr(basis,ish,ilo,ihi))
     &            call errquit('vvvv: bas_cn2bfr',ish, BASIS_ERR)
              do jsh=1,nsh
                if (.not. bas_cn2bfr(basis,jsh,jlo,jhi))
     &              call errquit('vvvv: bas_cn2bfr',jsh, BASIS_ERR)
                do i=ilo,ihi 
                  do j=jlo,jhi 
                    ad1=ad1+1
                    scrb((i-1)*nbf+j)=scra(ad1)
                  enddo
                enddo
              enddo
            enddo
            call dgemm('t','n',nact,nbf,nbf,1.0d00,cmo(1,ncor+1),nbf,
     &                 scrb,nbf,0.0d00,scra,nact)
            call dgemm('n','n',nact,nact,nbf,1.0d00,scra,nact,
     &                 cmo(1,ncor+1),nbf,0.0d00,scrb,nact)
            call ga_put(g_x,1,nact*nact,ad2,ad2,scrb,nact*nact)
            call ga_get(g_c,1,nact*nbf,ad2,ad2,scra,nact*nbf)
            call dgemm('n','n',nact,nact,nbf,1.0d00,scra,nact,
     &                 cmo(1,ncor+1),nbf,0.0d00,scrb,nact)
            call ga_put(g_c,1,nact*nact,ad2,ad2,scrb,nact*nact)
          endif
        enddo
      enddo
      if (occsdps) then
         call pstat_off(ps_idx34)
      else
         call qexit('idx34',0)
      endif
      return
      end     
