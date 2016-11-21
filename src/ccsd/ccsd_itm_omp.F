      subroutine ccsd_itm_omp(nocc,nvir,nact,hiu,giu,habe,gabe,
     &                        hia,t1,ht1,
     &                        scra,scrb,scrc,scrd,scre,scrf,
     &                        g_nt2,g_nht2,g_ncoul,g_nexch,
     &                        g_x,g_c,g_sht2,
     &                        iprt)
!     $Id: ccsd_itm.F 27896 2015-12-04 22:29:52Z edo $
      implicit none
#include "errquit.fh"
      integer, intent(in)             :: nocc,nvir,nact
      integer, intent(in)             :: g_nt2,g_nht2,g_ncoul,g_nexch
      integer, intent(in)             :: g_x,g_c,g_sht2
      integer, intent(in)             :: iprt
      double precision, intent(in)    :: t1(nocc,nvir)
      double precision, intent(inout) :: hiu(nocc,nocc)
      double precision, intent(inout) :: giu(nocc,nocc)
      double precision, intent(inout) :: hia(nocc,nvir)
      double precision, intent(inout) :: habe(nvir,nvir)
      double precision, intent(inout) :: gabe(nvir,nvir)
      double precision, intent(inout) :: ht1(nocc,nvir)
      ! allocated externally, challenge for thread-safety
      ! unused in the OpenMP version of this function
      double precision, intent(in) :: scra(nact*nact)
      double precision, intent(in) :: scrb(nact*nact)
      double precision, intent(in) :: scrc(nact,nact)
      double precision, intent(in) :: scrd(nact,nact)
      double precision, intent(inout) :: scre(nact*nact)
      double precision, intent(inout) :: scrf(nact*nact)
      ! if stack buffers, may segfault if large.
      ! allocating these hits IFRTL global lock issue.
      ! allocate num_max_thrd*nbf*nbf array outside and
      ! index to private portion using thread id.
      double precision, allocatable   :: bufa(:,:)
      double precision, allocatable   :: bufb(:,:)
      double precision, allocatable   :: bufc(:,:,:)
      double precision, allocatable   :: bufd(:,:,:)
      double precision, allocatable   :: bufe(:,:)
      double precision, allocatable   :: buff(:,:)
      integer :: alloc_error
!
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "rtdb.fh"
#include "msgids.fh"
#include "ccsd_debug.fh"
#include "ccsdps.fh"
!
      integer :: g_aitm,g_tmt1,g_tmt2,g_tmt,g_tmp,g_tmp2
      integer :: g_jlo_nexch,g_jhi_nexch,g_ilo_nexch,g_ihi_nexch
      integer :: g_jlo_sht2,g_jhi_sht2,g_ilo_sht2,g_ihi_sht2
      integer :: g_jlo_x,g_jhi_x,g_ilo_x,g_ihi_x
      integer :: g_jlo_c,g_jhi_c,g_ilo_c,g_ihi_c
      integer :: g_jlo_nt2,g_jhi_nt2,g_ilo_nt2,g_ihi_nt2
      integer :: g_jlo_nht2,g_jhi_nht2,g_ilo_nht2,g_ihi_nht2
      integer :: g_jlo_tmp2,g_jhi_tmp2,g_ilo_tmp2,g_ihi_tmp2
      integer :: i,j,k,l,ij,a,b,c
      integer :: ad1,ad2,ad3,ad4,ad5,ad6
      integer :: lnoo,lnov,lnvv
      integer :: inode,next
      integer :: nodes, iam
      double precision :: temp
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
!
      integer :: nxtask
      external nxtask
!
      nodes = ga_nnodes()
      iam = ga_nodeid()
!
      if (occsdps) then
         call pstat_on(ps_itm)
      else
         call qenter('itm',0)
      endif
!
      lnoo=nocc*nocc
      lnov=nocc*nvir
      lnvv=nvir*nvir
!
! - test x and c by doing its contribution to (2j-k)*(2t-t)
      if (.not.ga_create(MT_DBL,lnoo,lnoo,'aitm',nocc,nocc,g_aitm))
     &     call errquit('ccsd_itm: ga_create g_aitm failed',0, GA_ERR)
      if (.not.ga_create(MT_DBL,lnvv,lnoo,'tmt1',lnvv,0,g_tmt1))
     &     call errquit('ccsd_itm: ga_create g_tmt1 failed',0, GA_ERR)
      if (.not.ga_create(MT_DBL,lnvv,lnoo,'tmt2',lnvv,0,g_tmt2))
     &     call errquit('ccsd_itm: ga_create g_tmt2 failed',0, GA_ERR)
      !
#if defined(USE_OPENMP)
      nt = omp_get_max_threads()
#else
      nt = 1
#endif
      !
      allocate( bufa(nact*nact,nt), stat=alloc_error)
      if (alloc_error.ne.0) call errquit('ccsd_itm_omp: bufa',1,MA_ERR)
      allocate( bufb(nact*nact,nt), stat=alloc_error)
      if (alloc_error.ne.0) call errquit('ccsd_itm_omp: bufb',2,MA_ERR)
      allocate( bufc(nact,nact,nt), stat=alloc_error)
      if (alloc_error.ne.0) call errquit('ccsd_itm_omp: bufc',3,MA_ERR)
      allocate( bufd(nact,nact,nt), stat=alloc_error)
      if (alloc_error.ne.0) call errquit('ccsd_itm_omp: bufd',4,MA_ERR)
      allocate( bufe(nact*nact,nt), stat=alloc_error)
      if (alloc_error.ne.0) call errquit('ccsd_itm_omp: bufe',5,MA_ERR)
      allocate( buff(nact*nact,nt), stat=alloc_error)
      if (alloc_error.ne.0) call errquit('ccsd_itm_omp: buff',6,MA_ERR)
      !
      call ga_distribution(g_nexch,iam,g_jlo_nexch,g_jhi_nexch,
     &                                 g_ilo_nexch,g_ihi_nexch)
      call ga_distribution(g_sht2,iam,g_jlo_sht2,g_jhi_sht2,
     &                                g_ilo_sht2,g_ihi_sht2)
      call ga_distribution(g_x,iam,g_jlo_x,g_jhi_x,
     &                             g_ilo_x,g_ihi_x)
      call ga_distribution(g_nt2,iam,g_jlo_nt2,g_jhi_nt2,
     &                               g_ilo_nt2,g_ihi_nt2)
      call ga_distribution(g_nht2,iam,g_jlo_nht2,g_jhi_nht2,
     &                                g_ilo_nht2,g_ihi_nht2)
      call ga_distribution(g_c,iam,g_jlo_c,g_jhi_c,g_ilo_c,g_ihi_c)
!
!$omp parallel 
!$omp& default(none)
!$omp& shared(t1)
!$omp& shared(g_nexch,g_nt2,g_sht2,g_tmt1,g_tmt2)
!$omp& shared(g_c,g_x,g_nht2,g_aitm)
!$omp& shared(g_jlo_x,g_jhi_x,g_ilo_x,g_ihi_x)
!$omp& shared(g_jlo_sht2,g_jhi_sht2,g_ilo_sht2,g_ihi_sht2)
!$omp& shared(g_jlo_nt2,g_jhi_nt2,g_ilo_nt2,g_ihi_nt2)
!$omp& shared(g_jlo_nht2,g_jhi_nht2,g_ilo_nht2,g_ihi_nht2)
!$omp& shared(g_jlo_nexch,g_jhi_nexch,g_ilo_nexch,g_ihi_nexch)
!$omp& shared(nact,nocc,nvir,iam,iprt)
!$omp& shared(lnoo,lnov,lnvv)
!$omp& private(ij,tid)
!$omp& private(ad1,ad2,ad3,ad4,ad5,ad6)
!$omp& private(scre,scrf)
!$omp& shared(bufa,bufb,bufc,bufd,bufe,buff)
!$omp& shared(ht1,gabe,habe,giu,hiu,hia)
!$omp& shared(doa,dob)
!
#if defined(USE_OPENMP)
      ! OpenMP thread ids run from 0 to nt-1
      tid = omp_get_thread_num() + 1
#else
      tid = 1
#endif
!
!$omp do collapse(2) schedule(static)
      do k=1,nvir
        do l=1,nvir
          habe(l,k) = 0.0d0
        enddo
      enddo
!$omp end do nowait
!$omp do collapse(2) schedule(static)
      do k=1,nocc
        do l=1,nocc
          hiu(l,k) = 0.0d0
        enddo
      enddo
!$omp end do nowait
!$omp do collapse(2) schedule(static)
      do k=1,nvir
        do l=1,nocc
          hia(l,k) = 0.0d0
        enddo
      enddo
!$omp end do nowait
!
!$omp do collapse(2) schedule(guided)
!$omp& reduction(+:ht1)
!$omp& reduction(+:habe)
      do i=1,nocc
        do j=1,nocc
          ad1=(i-1)*nact
          ad2=(j-1)*nact
          ad3=(i-1)*nvir
          ad4=(j-1)*nvir
          if ((ad1+1.ge.g_ilo_nexch).and.(ad1+1.le.g_ihi_nexch)) then
          if ((ad2+1.ge.g_jlo_nexch).and.(ad2+1.le.g_jhi_nexch)) then
!-----------------------------------------------------------------------
!  S7
!-----------------------------------------------------------------------
!$omp critical
            call ga_get(g_nexch,ad2+1,ad2+nact,ad1+1,ad1+nact,
     &                  bufc(1,1,tid),nact)
            call ga_get(g_nt2,ad4+1,ad4+nvir,ad3+1,ad3+nvir,
     &                  scre,nvir)
!$omp end critical
            do k=1,nocc
            do a=1,nvir
            do b=1,nvir
            ht1(k,a)=ht1(k,a)-(2.0d0*bufc(nocc+b,k,tid)
     &                              -bufc(k,nocc+b,tid))*
     &                        (scre((a-1)*nvir+b)+t1(i,a)*t1(j,b))
            enddo
            enddo
            enddo
            do c=1,nvir
            do a=1,nvir
            do b=1,nvir
            habe(a,b)=habe(a,b)-(2.0d0*bufc(nocc+b,nocc+c,tid)
     &                                -bufc(nocc+c,nocc+b,tid))*
     &                          (scre((c-1)*nvir+a)+t1(i,c)*t1(j,a))
            enddo
            enddo
            enddo
          endif
          endif
        enddo
      enddo
!$omp end do
!
!$omp do collapse(2) schedule(guided)
!$omp& reduction(+:hiu)
      do i=1,nocc
        do j=1,nocc
          ij=(i-1)*nocc+j
          if (ij.ge.g_ilo_sht2.and.ij.le.g_ihi_sht2)then
!$omp critical
            call ga_get(g_sht2,1,nact*nact,ij,ij,scre,nact*nact)
!$omp end critical
            !call dscal(nact*nact,2.0d0,scre,1)
            do k=1,nact*nact
                scre(k) = scre(k) * 2.0d0
            enddo
!
! H2
            do k=1,nocc
                hiu(k,i)=hiu(k,i)+scre((k-1)*nact+j)
     &                           +scre((k-1)*nact+j)
     &                           -scre((j-1)*nact+k)
            enddo
          endif
        enddo
      enddo
!$omp end do
!
!$omp single
      if (iprt.gt.5.and.iam.eq.0)print *,'forming aitm'
!$omp end single
!
!$omp do collapse(2) schedule(static)
      do k=1,nocc
        do l=1,nocc
          giu(l,k) = hiu(l,k)
        enddo
      enddo
!$omp end do nowait
!$omp do collapse(2) schedule(static)
      do k=1,nvir
        do l=1,nvir
          gabe(l,k) = habe(l,k)
        enddo
      enddo
!$omp end do
!
!$omp do collapse(2) schedule(guided)
!$omp& reduction(+:hia)
!$omp& reduction(+:ht1)
!$omp& reduction(+:giu)
!$omp& reduction(+:gabe)
      do i=1,nocc
        do j=1,nocc
          ij=(i-1)*nocc+j
          if (ij.ge.g_ilo_x.and.ij.le.g_ihi_x)then
            ad1=(i-1)*nact
            ad2=(j-1)*nact
!
!-----------------------------------------------------------------------
!  A1, A2, A3 and A4 intermediates
!-----------------------------------------------------------------------
!
! A1
            do k=1,nact*nact
                bufa(k,tid) = 0.0d0
            enddo
            if (doa(1).eq.2)then
!$omp critical
              call ga_get(g_nexch,ad2+1,ad2+nocc,ad1+1,ad1+nocc,
     &                    bufb(1,tid),nocc)
!$omp end critical
              do k=1,nocc*nocc
                  bufa(k,tid) = bufa(k,tid) + 0.5d0 * bufb(k,tid)
              enddo
            endif
!---------------------------------------------------------
!  S5, D1, D6 and D7 
!---------------------------------------------------------
!$omp critical
            call ga_get(g_nexch,ad2+nocc+1,ad2+nact,ad1+1,ad1+nact,
     &                  scre,nvir)
!$omp end critical
            do k=1,nact*nact
                scrf(k) = 0.0d0
            enddo
!
! D1,D7
            ad1=0
            do a=1,nvir
              do b=nocc+1,nact
                ad1=ad1+1
                scrf(ad1)=scrf(ad1)+0.5d0*scre((b-1)*nvir+a)
                do k=1,nocc
                  scrf(ad1)=scrf(ad1)-scre((k-1)*nvir+a)*t1(k,b-nocc)
                enddo
              enddo
            enddo
!
            ad1=(i-1)*nocc+j
!$omp critical
            call ga_get(g_c,1,nact*nact,ad1,ad1,bufb(1,tid),nact*nact)
!$omp end critical
            do k=1,nact*nact
                scre(k) = bufb(k,tid)
            enddo
!
! H8
            do a=nocc+1,nact
              hia(j,a-nocc)=hia(j,a-nocc)-scre((a-1)*nact+i)
            enddo
!
! G4
            if (i.eq.j) then
             do b=nocc+1,nact
              do a=nocc+1,nact
                gabe(a-nocc,b-nocc)=gabe(a-nocc,b-nocc)
     &                             +scre((a-1)*nact+b)
     &                             +scre((a-1)*nact+b)
              enddo
             enddo
            endif
!
            ad1=(i-1)*nocc+j
!$omp critical
            call ga_get(g_x,1,nact*nact,ad1,ad1,bufb(1,tid),nact*nact)
!$omp end critical
            do k=1,nact*nact
                scre(k) = bufb(k,tid)
            enddo
!
! H8
            do a=nocc+1,nact
              hia(j,a-nocc)=hia(j,a-nocc)
     &                     +scre((i-1)*nact+a)
     &                     +scre((i-1)*nact+a)
            enddo
!
! S5
            do a=nocc+1,nact
              ht1(j,a-nocc)=ht1(j,a-nocc)
     &                     +scre((i-1)*nact+a)
     &                     +scre((i-1)*nact+a)
     &                     -scre((a-1)*nact+i)
            enddo
!
! G2
            do k=1,nocc
              giu(k,j)=giu(k,j)
     &                +scre((i-1)*nact+k)
     &                +scre((i-1)*nact+k)
     &                -scre((k-1)*nact+i)
            enddo
!
! G4
            if (i.eq.j) then
             do b=nocc+1,nact
              do a=nocc+1,nact
                gabe(a-nocc,b-nocc)=gabe(a-nocc,b-nocc)
     &                             -scre((a-1)*nact+b)
              enddo
             enddo
            endif
!
! D6,D7
            ad1=0
            do a=nocc+1,nact
              do b=nocc+1,nact
                ad1=ad1+1
                scrf(ad1)=scrf(ad1)+scre((b-1)*nact+a)
                do k=1,nocc
                  scrf(ad1)=scrf(ad1)-scre((b-1)*nact+k)*t1(k,a-nocc)
                  scrf(ad1)=scrf(ad1)-scre((k-1)*nact+a)*t1(k,b-nocc)
                enddo
              enddo
            enddo
            ad1=(i-1)*nvir
            ad2=(j-1)*nvir
!$omp critical
            call ga_acc(g_nht2,ad1+1,ad1+nvir,ad2+1,ad2+nvir,
     &                  scrf,nvir,1.0d0)
!$omp end critical
!
            if (doa(2).eq.2.and.doa(3).eq.2) then
            ad1=0
            do k=1,nocc
              do l=1,nocc
                ad1=ad1+1
                bufa(ad1,tid)=bufa(ad1,tid)+0.5d0*bufb((k-1)*nact+l,tid)
              enddo
            enddo
            ad1=(j-1)*nocc+i
!$omp critical
            call ga_get(g_x,1,nact*nact,ad1,ad1,bufb(1,tid),nact*nact)
!$omp end critical
            ad1=0
            do k=1,nocc
              do l=1,nocc
                ad1=ad1+1
                bufa(ad1,tid)=bufa(ad1,tid)+0.5d0*bufb((l-1)*nact+k,tid)
              enddo
            enddo
            endif
!
!$omp critical
            call ga_get(g_sht2,1,nact*nact,ij,ij,bufb(1,tid),nact*nact)
!$omp end critical
            do k=1,nact*nact
                scre(k) = 0.0d0
            enddo
            do k=1,nact*nact
                scre(k) = scre(k) + 2.0d0 * bufb(k,tid)
            enddo
!
! S6
            do a=nocc+1,nact
              ht1(j,a-nocc)=ht1(j,a-nocc)
     &                     +scre((i-1)*nact+a)
     &                     +scre((i-1)*nact+a)
     &                     -scre((a-1)*nact+i)
            enddo
!
! A4
            if (doa(4).eq.2)then
             ad1=0
             do k=1,nocc
               do l=1,nocc
                 ad1=ad1+1
                 bufa(ad1,tid)=bufa(ad1,tid)+bufb((k-1)*nact+l,tid)
               enddo
             enddo
            endif
!$omp critical
            call ga_put(g_aitm,1,lnoo,ij,ij,bufa(1,tid),lnoo)
!$omp end critical
!-----------------------------------------------------------------------
!  B2 and B3 intermediate and their contribution to D3
!-----------------------------------------------------------------------
            if (dob(2).eq.2)then
            ad1=nocc+1
            call dgemm('n','n',nvir,nvir,nocc,1.0d0,bufb(ad1,tid),
     &                 nact,t1,nocc,0.0d0,bufa(1,tid),nvir)
            ad1=nocc*nact+1
            call dgemm('t','n',nvir,nvir,nocc,1.0d0,t1,nocc,
     &                 bufb(ad1,tid),nact,1.0d0,bufa(1,tid),nvir)
            ad1=(i-1)*nvir
            ad2=(j-1)*nvir
!$omp critical
            call ga_acc(g_nht2,ad2+1,ad2+nvir,ad1+1,ad1+nvir,
     &                  bufa(1,tid),nvir,-1.0d0)
!$omp end critical
            endif
!
          endif
        enddo
      enddo
!$omp end do
!
!$omp single
      call ga_dgop(msg_cc_term1,hiu,lnoo, '+')
      call ga_dgop(msg_cc_term2,habe,lnvv, '+')
      call ga_dgop(msg_cc_term3,giu,lnoo, '+')
      call ga_dgop(msg_cc_term4,hia,lnov, '+')
      call ga_dgop(msg_cc_term5,gabe,lnvv, '+')
!$omp end single
!
!-----------------------------------------------------------------------
!  D2 (product of A intermediate and Tau)
!-----------------------------------------------------------------------
! construct tau
!
!$omp do collapse(2) schedule(guided)
      do i=1,nocc
        do j=1,nocc
          ad1=(i-1)*nvir
          ad2=(j-1)*nvir
          if (ad1+1.ge.g_ilo_nt2.and.ad1+1.le.g_ihi_nt2)then
          if (ad2+1.ge.g_jlo_nt2.and.ad2+1.le.g_jhi_nt2)then
            ij=(i-1)*nocc+j
!$omp critical
            call ga_get(g_nt2,ad2+1,ad2+nvir,ad1+1,ad1+nvir,
     &                  bufa(1,tid),nvir)
!$omp end critical
            ad3=0
            do a=1,nvir
              do b=1,nvir
                ad3=ad3+1
                bufa(ad3,tid)=bufa(ad3,tid)+t1(i,a)*t1(j,b)
              enddo
            enddo
!$omp critical
            call ga_put(g_tmt1,1,lnvv,ij,ij,bufa(1,tid),lnvv)
!$omp end critical
          endif
          endif
        enddo
      enddo
!$omp end do
!
!$omp single
      call ga_dgemm('n','n',lnvv,lnoo,lnoo,
     &              1.0d0,g_tmt1,g_aitm,0.0d0,g_tmt2)
!$omp end single
!
!$omp do collapse(2)
      do i=1,nocc
        do j=1,nocc
          ad1=(i-1)*nvir
          ad2=(j-1)*nvir
          if (ad1+1.ge.g_ilo_nht2.and.ad1+1.le.g_ihi_nht2)then
          if (ad2+1.ge.g_jlo_nht2.and.ad2+1.le.g_jhi_nht2)then
            ij=(i-1)*nocc+j
!$omp critical
            call ga_get(g_tmt2,1,lnvv,ij,ij,bufa(1,tid),lnvv)
            call ga_acc(g_nht2,ad2+1,ad2+nvir,ad1+1,ad1+nvir,
     &                  bufa(1,tid),nvir,1.0d0)
!$omp end critical
          endif
          endif
        enddo
      enddo
!$omp end do
!
!-----------------------------------------------------------------------
!  S4
!-----------------------------------------------------------------------
!
! construct (2*t2_ij^ab - t2_ij^ba + t_i^a t_j^b)
!
!$omp do collapse(2) schedule(guided)
!$omp& reduction(+:ht1)
      do i=1,nocc
        do j=1,nocc
          ad1=(i-1)*nvir
          ad2=(j-1)*nvir
          if (ad1+1.ge.g_ilo_nt2.and.ad1+1.le.g_ihi_nt2)then
          if (ad2+1.ge.g_jlo_nt2.and.ad2+1.le.g_jhi_nt2)then
!$omp critical
            call ga_get(g_nt2,ad2+1,ad2+nvir,ad1+1,ad1+nvir,
     &                  bufa(1,tid),nvir)
!$omp end critical
            do a=1,nvir
              do b=1,nvir
              ad3=(b-1)*nvir+a
              ad4=(a-1)*nvir+b
              ht1(i,a)=ht1(i,a)+(bufa(ad4,tid)+bufa(ad4,tid)
     &                     -bufa(ad3,tid)+t1(i,b)*t1(j,a))*hia(j,b)
              enddo
            enddo
          endif
          endif
        enddo
      enddo
!$omp end do
!
!
!-----------------------------------------------------------------------
!  D4
!-----------------------------------------------------------------------
!
!$omp do collapse(2) schedule(guided)
      do i=1,nocc
        do j=1,nocc
          ad1=(i-1)*nvir
          ad2=(j-1)*nvir
          if ((ad1+1).ge.g_ilo_nht2.and.(ad1+1).le.g_ihi_nht2)then
          if ((ad2+1).ge.g_jlo_nht2.and.(ad2+1).le.g_jhi_nht2)then
!$omp critical
            call ga_get(g_nt2,ad2+1,ad2+nvir,ad1+1,ad1+nvir,scre,nvir)
!$omp end critical
            call dgemm('n','n',nvir,nvir,nvir,1.0d0,gabe,
     &                 nvir,scre,nvir,0.0d0,bufd(1,1,tid),nvir)
!$omp critical
            call ga_acc(g_nht2,ad2+1,ad2+nvir,ad1+1,ad1+nvir,
     &                  bufd(1,1,tid),nvir,1.0d0)
!$omp end critical
          endif
          endif
        enddo
      enddo
!$omp end do
!
!$omp single
      call ga_sync()
!$omp end single
!
!$omp end parallel
!
!-----------------------------------------------------------------------
!  D5
!-----------------------------------------------------------------------
!
      inode=-1
      do a=1,nvir
        do i=1,nocc
          ad1=(i-1)*nvir+a
          inode=inode+1
          if (mod(inode, nodes) .eq. iam) then
!$omp critical
            call ga_get(g_nt2,1,lnov,ad1,ad1,scre,lnov)
!$omp end critical
            call dgemm('n','n',nvir,nocc,nocc,-1.0d0,scre,
     &                 nvir,giu,nocc,0.0d0,scrf,nvir)
!$omp critical
            call ga_acc(g_nht2,1,lnov,ad1,ad1,scrf,lnov,1.0d0)
!$omp end critical
          endif
        enddo
      enddo
!
!$omp single
      call ga_sync()
      if (.not.ga_destroy(g_tmt2))
     &    call errquit('ga_dest g_tmt2 fail',0, GA_ERR)
      if (.not.ga_destroy(g_tmt1))
     &    call errquit('ga_dest g_tmt1 fail',0, GA_ERR)
      if (.not.ga_destroy(g_aitm))
     &    call errquit('ga_dest g_aitm fail',0, GA_ERR)
!$omp end single
!
!-----------------------------------------------------------------------
!  J1, J2, J3 and K1, K2, K3 intermediates
!-----------------------------------------------------------------------
!
!$omp parallel 
!$omp& default(none)
!$omp& shared(t1)
!$omp& shared(g_nexch,g_nt2,g_sht2,g_tmt1,g_tmt2)
!$omp& shared(g_c,g_x,g_nht2,g_aitm,g_ncoul,g_tmt,g_tmp,g_tmp2)
!$omp& shared(g_jlo_c,g_jhi_c,g_ilo_c,g_ihi_c)
!$omp& shared(g_jlo_x,g_jhi_x,g_ilo_x,g_ihi_x)
!$omp& shared(g_jlo_nt2,g_jhi_nt2,g_ilo_nt2,g_ihi_nt2)
!$omp& shared(g_jlo_nht2,g_jhi_nht2,g_ilo_nht2,g_ihi_nht2)
!$omp& shared(g_jlo_tmp2,g_jhi_tmp2,g_ilo_tmp2,g_ihi_tmp2)
!$omp& shared(g_jlo_nexch,g_jhi_nexch,g_ilo_nexch,g_ihi_nexch)
!$omp& shared(nact,nocc,nvir,nodes,iam,iprt)
!$omp& shared(lnoo,lnov,lnvv)
!$omp& shared(doa,dob)
!$omp& shared(hia,ht1,giu,gabe,habe)
!$omp& private(ad1,ad2,ad3,ad4,ad5,ad6)
!$omp& private(ij,tid)
!$omp& private(scre,scrf)
!$omp& shared(bufa,bufb,bufc,bufd,bufe,buff)
!
#if defined(USE_OPENMP)
      ! OpenMP thread ids run from 0 to nt-1
      tid = omp_get_thread_num() + 1
#else
      tid = 1
#endif
!
!$omp do collapse(2) schedule(guided)
      do i=1,nocc
        do j=1,nocc
          ij=(i-1)*nocc+j
          if (ij.ge.g_ilo_x.and.ij.le.g_ihi_x)then
            ad1=(i-1)*nact
            ad2=(j-1)*nact
!$omp critical
            call ga_get(g_nexch,ad2+1,ad2+nact,ad1+1,ad1+nact,
     &                  bufd(1,1,tid),nact)
!$omp end critical
            call dgemm('n','n',nvir,nvir,nocc,-1.0d0,bufd(nocc+1,1,tid),
     &                 nact,t1,nocc,0.0d0,bufc(nocc+1,nocc+1,tid),nact)
            do k=1,nact
              do l=1,nact
                bufd(l,k,tid) = bufd(l,k,tid) + bufc(l,k,tid)
              enddo
            enddo
!$omp critical
            call ga_acc(g_x,1,nact*nact,ij,ij,bufd(1,1,tid),
     &                  nact*nact,1.0d0)
!$omp end critical
          endif
        enddo
      enddo
!$omp end do
!
!$omp do collapse(2) schedule(guided)
      do i=1,nocc
        do j=1,nocc
          ij=(i-1)*nocc+j
          if (ij.ge.g_ilo_c.and.ij.le.g_ihi_c)then
            ad1=(i-1)*nact
            ad2=(j-1)*nact
!$omp critical
            call ga_get(g_ncoul,ad2+1,ad2+nact,ad1+1,ad1+nact,
     &                  bufd(1,1,tid),nact)
!$omp end critical
            call dgemm('n','n',nvir,nvir,nocc,-1.0d0,bufd(nocc+1,1,tid),
     &                 nact,t1,nocc,0.0d0,bufc(nocc+1,nocc+1,tid),nact)
            do k=1,nact
              do l=1,nact
                bufd(l,k,tid) = bufd(l,k,tid) + bufc(l,k,tid)
              enddo
            enddo
!$omp critical
            call ga_acc(g_c,1,nact*nact,ij,ij,bufd(1,1,tid),
     &                  nact*nact,1.0d0)
!$omp end critical
          endif
        enddo
      enddo
!$omp end do
!
!$omp single
      if (.not.ga_create(MT_DBL,lnov,lnov,'tmt',nvir,nvir,g_tmt))
     &     call errquit('ccsd_itm: ga_create g_tmt failed',0, GA_ERR)
      if (.not.ga_create(MT_DBL,lnov,lnov,'tmp',nvir,nvir,g_tmp))
     &     call errquit('ccsd_itm: ga_create g_tmp failed',0, GA_ERR)
      if (.not.ga_create(MT_DBL,lnov,lnov,'tmp2',nvir,nvir,g_tmp2))
     &     call errquit('ccsd_itm: ga_create g_tmp2 failed',0, GA_ERR)
      call ga_distribution(g_tmp2,iam,g_jlo_tmp2,g_jhi_tmp2,
     &                                g_ilo_tmp2,g_ihi_tmp2)
!$omp end single
!
!-----------------------------------------------------------------------
!  K4 intermediate
!-----------------------------------------------------------------------
!
! construct (t2_ij^ab + 2*t_i^a*t_j^b)
!
!$omp do collapse(2) schedule(guided)
      do i=1,nocc
        do j=1,nocc
          ad1=(i-1)*nvir
          ad2=(j-1)*nvir
          if (ad1+1.ge.g_ilo_nt2.and.ad1+1.le.g_ihi_nt2)then
          if (ad2+1.ge.g_jlo_nt2.and.ad2+1.le.g_jhi_nt2)then
!$omp critical
            call ga_get(g_nt2,ad2+1,ad2+nvir,ad1+1,ad1+nvir,
     &                  bufa(1,tid),nvir)
!$omp end critical
            do a=1,nvir
              do b=1,nvir
              ad3=(a-1)*nvir+b
              bufa(ad3,tid)=bufa(ad3,tid)+2.0d0*t1(i,a)*t1(j,b)
              enddo
            enddo
!$omp critical
            call ga_put(g_tmt,ad1+1,ad1+nvir,ad2+1,ad2+nvir,
     &                  bufa(1,tid),nvir)
!$omp end critical
          endif
          endif
        enddo
      enddo
!$omp end do
!
!$omp do collapse(2) schedule(guided)
      do i=1,nocc
        do j=1,nocc
          ad2=(j-1)*nact
          ad1=(i-1)*nact
          if (ad1+1.ge.g_ilo_nexch.and.ad1+1.le.g_ihi_nexch)then
          if (ad2+1.ge.g_jlo_nexch.and.ad2+1.le.g_jhi_nexch)then
            ad3=(i-1)*nvir
            ad4=(j-1)*nvir
!$omp critical
            call ga_get(g_nexch,ad2+nocc+1,ad2+nact,ad1+nocc+1,ad1+nact,
     &                  bufa(1,tid),nvir)
            call ga_put(g_tmp,ad3+1,ad3+nvir,ad4+1,ad4+nvir,
     &                  bufa(1,tid),nvir)
!$omp end critical
          endif
          endif
        enddo
      enddo
!$omp end do
!
!$omp single
      call ga_dgemm('n','n',lnov,lnov,lnov,
     &              -0.5d0,g_tmp,g_tmt,0.0d0,g_tmp2)
!-----------------------------------------------------------------------
!  D9 and D10 (product of K intermediate and T2)
!-----------------------------------------------------------------------
!$omp end single
!
!$omp do collapse(2) schedule(guided)
      do i=1,nocc
        do j=1,nocc
          ad1=(i-1)*nvir
          ad2=(j-1)*nvir
          if (ad1+1.ge.g_ilo_tmp2.and.ad1+1.le.g_ihi_tmp2)then
          if (ad2+1.ge.g_jlo_tmp2.and.ad2+1.le.g_jhi_tmp2)then
            ij=(i-1)*nocc+j
!$omp critical
            call ga_get(g_c,1,nact*nact,ij,ij,bufb(1,tid),nact*nact)
!$omp end critical
            ad3=0
            do a=nocc+1,nact
              do b=nocc+1,nact
                ad3=ad3+1
                bufa(ad3,tid)=bufb((a-1)*nact+b,tid)
              enddo
            enddo
!$omp critical
            call ga_acc(g_tmp2,ad2+1,ad2+nvir,ad1+1,ad1+nvir,
     &                  bufa(1,tid),nvir,1.0d0)
!$omp end critical
          endif
          endif
        enddo
      enddo
!$omp end do
!
! construct (t2_ij^ba)
!
!$omp do collapse(2) schedule(guided)
      do i=1,nocc
        do j=1,nocc
          ad1=(i-1)*nvir
          ad2=(j-1)*nvir
          if (ad1+1.ge.g_ilo_nt2.and.ad1+1.le.g_ihi_nt2)then
          if (ad2+1.ge.g_jlo_nt2.and.ad2+1.le.g_jhi_nt2)then
!$omp critical
            call ga_get(g_nt2,ad1+1,ad1+nvir,ad2+1,ad2+nvir,
     &                  bufa(1,tid),nvir)
            call ga_put(g_tmt,ad2+1,ad2+nvir,ad1+1,ad1+nvir,
     &                  bufa(1,tid),nvir)
!$omp end critical
          endif
          endif
        enddo
      enddo
!$omp end do
!
!$omp single
      call ga_dgemm('t','n',lnov,lnov,lnov,
     &              1.0d0,g_tmp2,g_tmt,0.0d0,g_tmp)
!$omp end single
!
!     d9 and d10
!
!$omp do collapse(2) schedule(guided)
      do i=1,nocc
        do j=1,nocc
          ad1=(i-1)*nvir
          ad2=(j-1)*nvir
          if (ad1+1.ge.g_ilo_nht2.and.ad1+1.le.g_ihi_nht2)then
          if (ad2+1.ge.g_jlo_nht2.and.ad2+1.le.g_jhi_nht2)then
!$omp critical
            call ga_get(g_tmp,ad2+1,ad2+nvir,ad1+1,ad1+nvir,
     &                  bufa(1,tid),nvir)
!$omp end critical
            do a=1,nvir
              do b=1,nvir
                ad3=(a-1)*nvir+b
                ad4=(b-1)*nvir+a
                bufb(ad3,tid)=-0.5d0*bufa(ad3,tid)-bufa(ad4,tid)
              enddo
            enddo
!$omp critical
            call ga_acc(g_nht2,ad2+1,ad2+nvir,ad1+1,ad1+nvir,
     &                  bufb(1,tid),nvir,1.0d0)
!$omp end critical
          endif
          endif
        enddo
      enddo
!$omp end do
!
!$omp single
      call ga_sync()
      if (iprt.gt.5.and.iam.eq.0)print *,'forming j - k'
!$omp end single
!
!$omp do collapse(2) schedule(guided)
      do i=1,nocc
        do j=1,nocc
          ij=(i-1)*nocc+j
          if (ij.ge.g_ilo_x.and.ij.le.g_ihi_x)then
            ad1=(i-1)*nact
            ad2=(j-1)*nact
            do k=1,nact*nact
                bufa(k,tid) = 0.0d0
            enddo
!$omp critical
            call ga_get(g_x,1,nact*nact,ij,ij,bufb(1,tid),nact*nact)
!$omp end critical
            do k=1,nact*nact
                bufa(k,tid) = bufa(k,tid) + bufb(k,tid)
            enddo
!$omp critical
            call ga_get(g_c,1,nact*nact,ij,ij,bufb(1,tid),nact*nact)
!$omp end critical
            do k=1,nact*nact
                bufa(k,tid) = bufa(k,tid) - 0.5d0 * bufb(k,tid)
            enddo
            ad1=0
            do a=nocc+1,nact
              do b=nocc+1,nact
                ad1=ad1+1
                bufb(ad1,tid)=bufa((a-1)*nact+b,tid)
              enddo
            enddo
            ad1=(j-1)*nvir
            ad2=(i-1)*nvir
!$omp critical
            call ga_put(g_tmp,ad1+1,ad1+nvir,ad2+1,ad2+nvir,
     &                  bufb(1,tid),nvir)
!$omp end critical
          endif
        enddo
      enddo
!$omp end do
!-----------------------------------------------------------------------
!  J4 + J5 - K4 intermediate 
!-----------------------------------------------------------------------
!
! construct (2*t2_ij^ab - t2_ij^ba - 2*t_i^a t_j^b)
!
!$omp do collapse(2) schedule(guided)
      do i=1,nocc
        do j=1,nocc
          ad1=(i-1)*nvir
          ad2=(j-1)*nvir
          if (ad1+1.ge.g_ilo_nt2.and.ad1+1.le.g_ihi_nt2)then
          if (ad2+1.ge.g_jlo_nt2.and.ad2+1.le.g_jhi_nt2)then
!$omp critical
            call ga_get(g_nt2,ad2+1,ad2+nvir,ad1+1,ad1+nvir,
     &                  bufa(1,tid),nvir)
!$omp end critical
            do a=1,nvir
              do b=1,nvir
                ad3=(b-1)*nvir+a
                ad4=(a-1)*nvir+b
                bufb(ad3,tid)=bufa(ad3,tid)+bufa(ad3,tid)-bufa(ad4,tid)
     &                    -2.0d0*t1(i,a)*t1(j,b)
              enddo
            enddo
!$omp critical
            call ga_put(g_tmt,ad2+1,ad2+nvir,ad1+1,ad1+nvir,
     &                  bufb(1,tid),nvir)
!$omp end critical
          endif
          endif
        enddo
      enddo
!$omp end do
!
!$omp do collapse(2) schedule(guided)
      do i=1,nocc
        do j=1,nocc
          ad1=(i-1)*nact
          ad2=(j-1)*nact
          if (ad1+1.ge.g_ilo_nexch.and.ad1+1.le.g_ihi_nexch)then
          if (ad2+1.ge.g_jlo_nexch.and.ad2+1.le.g_jhi_nexch)then
            ad3=(i-1)*nvir
            ad4=(j-1)*nvir
!$omp critical
            call ga_get(g_nexch,ad2+nocc+1,ad2+nact,ad1+nocc+1,ad1+nact,
     &                  bufa(1,tid),nvir)
!$omp end critical
              do a=1,nvir
                do b=1,nvir
                  ad5=(a-1)*nvir+b
                  ad6=(b-1)*nvir+a
                  bufb(ad5,tid)=bufa(ad5,tid)+bufa(ad5,tid)
     &                         -bufa(ad6,tid)
                enddo
              enddo
!$omp critical
            call ga_put(g_tmp2,ad4+1,ad4+nvir,ad3+1,ad3+nvir,
     &                  bufb(1,tid),nvir)
!$omp end critical
          endif
          endif
        enddo
      enddo
!$omp end do
!
!$omp single
      call ga_dgemm('n','t',lnov,lnov,lnov,
     &              0.25d0,g_tmp2,g_tmt,1.0d0,g_tmp)
!$omp end single
!
! construct (2*t2_ij^ab - t2_ij^ba)
!
!$omp do collapse(2) schedule(guided)
      do i=1,nocc
        do j=1,nocc
          ad1=(i-1)*nvir
          ad2=(j-1)*nvir
          if (ad1+1.ge.g_ilo_nt2.and.ad1+1.le.g_ihi_nt2)then
          if (ad2+1.ge.g_jlo_nt2.and.ad2+1.le.g_jhi_nt2)then
!$omp critical
            call ga_get(g_nt2,ad2+1,ad2+nvir,ad1+1,ad1+nvir,
     &                  bufa(1,tid),nvir)
!$omp end critical
            do a=1,nvir
              do b=1,nvir
                ad3=(a-1)*nvir+b
                ad4=(b-1)*nvir+a
                bufb(ad3,tid)=bufa(ad3,tid)+bufa(ad3,tid)-bufa(ad4,tid)
              enddo
            enddo
!$omp critical
            call ga_put(g_tmt,ad2+1,ad2+nvir,ad1+1,ad1+nvir,
     &                  bufb(1,tid),nvir)
!$omp end critical
          endif
          endif
        enddo
      enddo
!$omp end do
!$omp end parallel
!
!$omp single
      call ga_dgemm('t','n',lnov,lnov,lnov,
     &              1.0d0,g_tmp,g_tmt,1.0d0,g_nht2)
      if (.not.ga_destroy(g_tmp2))
     &    call errquit('ga_dest g_tmp2 fail',0, GA_ERR)
      if (.not.ga_destroy(g_tmp))
     &    call errquit('ga_dest g_tmp fail',0, GA_ERR)
      if (.not.ga_destroy(g_tmt))
     &    call errquit('ga_dest g_tmt fail',0, GA_ERR)
      call ga_sync()
!$omp end single
!
#if 0
!
! There is absolutely no good reason to use NXTVAL to
! load-balance such a trivial loop...
!
! term s1 n**3
      inode=-1
      next=nxtask(nodes, 1)
      do a=1,nvir
       do b=1,nvir
        inode=inode+1
        if (inode.eq.next)then
         if (dos(1).ne.0)then
          do k=1,nocc
            ht1(k,b) = ht1(k,b) + habe(b,a) * t1(k,a)
          enddo
         endif
!$omp critical
         next=nxtask(nodes, 1)
!$omp end critical
        endif
       enddo
      enddo
      call ga_sync()
!
      inode=-1
      next=nxtask(-nodes, 1)
      call ga_sync()
! -- add some trivial terms to ht1 (nothing to ht2)
! term s2 n**3
      do a=1,nvir
       do i=1,nocc
        inode=inode+1
        if (mod(inode,nodes).eq.iam)then
         if (dos(2).ne.0)then
          temp = 0.0d0
          do k=1,nocc
            temp = temp + hiu(k,i) * t1(k,a)
          enddo
          ht1(i,a) = ht1(i,a) - temp
         endif
        endif
       enddo
      enddo
#else
! -- add some trivial terms to ht1 (nothing to ht2)
! term s1 n**3
      if (iam.eq.0) then
        if (dos(1).ne.0)then
#if 0
          do k=1,nocc
            do b=1,nvir
              do a=1,nvir
                ht1(k,b) = ht1(k,b) + t1(k,a) * habe(b,a)
              enddo
            enddo
          enddo
#else
          ! DGEMM: A * B^T
          ! alpha = 1.0d0
          ! A = t1
          ! B = habe
          ! beta = 1.0d0
          ! C = ht1
          call dgemm('n','t',nocc,nvir,nvir,
     &               1.0d0,t1,nocc,habe,nvir,
     &               1.0d0,ht1,nocc)
#endif
        endif ! dos(1)
! term s2 n**3
        if (dos(2).ne.0)then
#if 0
          do i=1,nocc
            do a=1,nvir
              do k=1,nocc
                ht1(i,a) = ht1(i,a) - hiu(k,i) * t1(k,a)
              enddo
            enddo
          enddo
#else
          ! DGEMM: A^T * B
          ! alpha = -1.0d0
          ! A = hiu
          ! B = t1
          ! beta = 1.0d0
          ! C = ht1
          call dgemm('t','n',nocc,nvir,nocc,
     &               -1.0d0,hiu,nocc,t1,nocc,
     &               1.0d0,ht1,nocc)
#endif
        endif ! dos(2)
      endif ! iam.eq.0
      call ga_sync()
#endif
      call ga_dgop(msg_cc_term6,ht1,lnov, '+')
      !
      deallocate( bufa, stat=alloc_error)
      deallocate( bufb, stat=alloc_error)
      deallocate( bufc, stat=alloc_error)
      deallocate( bufd, stat=alloc_error)
      deallocate( bufe, stat=alloc_error)
      deallocate( buff, stat=alloc_error)
      if (alloc_error.ne.0) call errquit('ccsd_itm_omp: bufa', 7,MA_ERR)
      if (alloc_error.ne.0) call errquit('ccsd_itm_omp: bufb', 8,MA_ERR)
      if (alloc_error.ne.0) call errquit('ccsd_itm_omp: bufc', 9,MA_ERR)
      if (alloc_error.ne.0) call errquit('ccsd_itm_omp: bufd',10,MA_ERR)
      if (alloc_error.ne.0) call errquit('ccsd_itm_omp: bufe',11,MA_ERR)
      if (alloc_error.ne.0) call errquit('ccsd_itm_omp: buff',12,MA_ERR)
      !
      if (occsdps) then
         call pstat_off(ps_itm)
      else
         call qexit('itm',0)
      endif
      return
      end
