#if defined(USE_OPENMP)
      subroutine ccsd_idx34(basis,cmo,
     &                      scra,scrb, ! unused
     &                      nsh,ncor,nocc,nact,nbf,
     &                      g_x,g_c)
      implicit none
#include "errquit.fh"
      integer, intent(in) ::  basis,nsh,ncor,nocc,nact,nbf,g_x,g_c
      double precision, intent(in) :: cmo(nbf,nbf)
      ! allocated externally, challenge for thread-safety
      ! unused in the OpenMP version of this function
      double precision, intent(inout) :: scra(nbf*nbf),scrb(nbf*nbf)
      ! stack buffers, may segfault if large
      ! allocating these hits IFRTL global lock issue
      ! allocate num_max_thrd*nbf*nbf array outside and
      ! have TLS pointer to chunk for each thread
      double precision :: bufa(nbf*nbf),bufb(nbf*nbf)
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
!$omp parallel do collapse(2) schedule(dynamic,1) !!!schedule(guided)
!$omp& default(none)
!$omp& firstprivate(nbf,nsh,ncor,nact,nocc,basis,g_ilo,g_ihi)
!$omp& private(ilo,ihi,jlo,jhi,ad1,ad2,bufa,bufb)
!$omp& shared(g_x,g_c,cmo)
      do ii=1,nocc
        do jj=1,nocc
          ad2=(ii-1)*nocc+jj
          if (ad2.ge.g_ilo.and.ad2.le.g_ihi)then
!$omp critical
            call ga_get(g_x,1,nbf*nbf,ad2,ad2,bufa,nbf*nbf)
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
                    bufb((i-1)*nbf+j)=bufa(ad1)
                  enddo
                enddo
              enddo
            enddo
            call dgemm('t','n',nact,nbf,nbf,1.0d00,cmo(1,ncor+1),nbf,
     &                 bufb,nbf,0.0d00,bufa,nact)
            call dgemm('n','n',nact,nact,nbf,1.0d00,bufa,nact,
     &                 cmo(1,ncor+1),nbf,0.0d00,bufb,nact)
!$omp critical
            call ga_put(g_x,1,nact*nact,ad2,ad2,bufb,nact*nact)
!$omp end critical
!$omp critical
            call ga_get(g_c,1,nact*nbf,ad2,ad2,bufa,nact*nbf)
!$omp end critical
            call dgemm('n','n',nact,nact,nbf,1.0d00,bufa,nact,
     &                 cmo(1,ncor+1),nbf,0.0d00,bufb,nact)
!$omp critical
            call ga_put(g_c,1,nact*nact,ad2,ad2,bufb,nact*nact)
!$omp end critical
          endif
        enddo
      enddo
!$omp end parallel do
      if (occsdps) then
         call pstat_off(ps_idx34)
      else
         call qexit('idx34',0)
      endif
      return
      end     
#else
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
#endif
