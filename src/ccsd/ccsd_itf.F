      subroutine ccsd_itf(basis,ncor,nocc,nvir,ndel,nbf,nsh,
     &                    hiu,giu,habe,gabe,hia,ht1,
     &                    idiis,cmo,eorb,iprt,
     &                    g_t2,g_ht2,g_ncoul,g_nexch,tklst)
      implicit none
      integer basis,ncor,nocc,nvir,ndel,nbf,nsh,
     &        idiis,iprt,g_t2,g_ht2
      real*8 cmo(*),eorb(*),hiu(*),giu(*),habe(*),gabe(*),ht1(*),hia(*)
#include "mafdecls.fh"
#include "bas.fh"
#include "rtdb.fh"
#include "global.fh"
c
      integer doa,dob,dog,doh,dojk,dos,dod
      common/debug/doa(15),dob(15),dog(15),doh(15),dojk(15),
     &             dos(15),dod(15)
c
      integer lnov,g_nt2,g_nht2,offt2,lnoov,iam,g_jlo,g_jhi,
     &        g_ilo,g_ihi,a,i,j,b,ad1,ad2,lnoo,lnvv,
     &        l_t1,k_t1,l_sa,k_sa,l_sb,k_sb,lsab,
     &        ij,lnof,g_ncoul,g_nexch,tklst(nsh*(nsh+1)/2,2)
      real*8 tol2e
      logical stat
c
      iam=ga_nodeid()
      lnoo=nocc*nocc
      lnov=nocc*nvir
      lnof=nocc*nbf
      lnoov=nocc*nocc*nvir
      lnvv=nvir*nvir
      tol2e=1.0d-12
      offt2=(idiis-1)*lnoov
c
      lsab=max(lnoov,nbf*nbf)
      stat=.true.
      stat=stat.and.ma_push_get(MT_DBL,lnov,'t1',l_t1, k_t1)
      stat=stat.and.ma_push_get(MT_DBL,lsab,'sa',l_sa, k_sa)
      stat=stat.and.ma_push_get(MT_DBL,lsab,'sb',l_sb, k_sb)
      if (.not.stat)call errquit('ma_push i_t_f ',0)
      call ga_get(g_t2,offt2+1,offt2+lnov,nvir+1,nvir+1,
     &            dbl_mb(k_t1),lnov)
c
c ------------------------------------------------------------
c create new GAs with proposed final ordering
c ------------------------------------------------------------
      if (.not.ga_create(MT_DBL,lnov,lnov,'nt2',
     &                   nvir,nvir,g_nt2))
     &     call errquit('ga_create g_nt2 failed',0)
      if (.not.ga_create(MT_DBL,lnov,lnov,'nht2',
     &                   nvir,nvir,g_nht2))
     &     call errquit('ga_create g_nht2 failed',0)
CRK   if (.not.ga_create(MT_DBL,lnof,lnof,'ncoul',
CRK  &                   nbf,nbf,g_ncoul))
CRK  &     call errquit('ga_create g_ncoul failed',0)
CRK   if (.not.ga_create(MT_DBL,lnof,lnof,'nexch',
CRK  &                   nbf,nbf,g_nexch))
CRK  &     call errquit('ga_create g_nexch failed',0)

c ------------------------------------------------------------
c fill new GAs
c ------------------------------------------------------------
      call ga_distribution(g_t2,iam,g_jlo,g_jhi,g_ilo,g_ihi)
      do a=1,nvir
        if (a.ge.g_ilo.and.a.le.g_ihi)then
          call ga_get(g_t2,offt2+1,offt2+lnoov,a,a,
     &                 dbl_mb(k_sa),lnoov)
          do i=1,nocc
            do j=1,nocc
              do b=1,nvir
                ad1=k_sa+(b-1)*lnoo+(i-1)*nocc+j-1
                ad2=k_sb+(j-1)*nvir+b-1
                dbl_mb(ad2)=dbl_mb(ad1)
              enddo
            enddo
            ad1=(i-1)*nvir+a
            call ga_put(g_nt2,1,lnov,ad1,ad1,dbl_mb(k_sb),lnov)
          enddo
        endif
      enddo
c
      call ga_zero(g_nht2)
c
CRK   call ga_distribution(g_coul,iam,g_jlo,g_jhi,g_ilo,g_ihi)
CRK   ij=0
CRK   do i=1,nocc
CRK     do j=1,i
CRK       ij=ij+1
CRK       if (ij.ge.g_ilo.and.ij.le.g_ihi)then
CRK         ad1=(i-1)*nbf
CRK         ad2=(j-1)*nbf
CRK         call ga_get(g_coul,1,nbf*nbf,ij,ij,dbl_mb(k_sa),nbf*nbf)
CRK         call ga_put(g_ncoul,ad2+1,ad2+nbf,ad1+1,ad1+nbf,
CRK  &                  dbl_mb(k_sa),nbf)
CRK         call mxtrn(dbl_mb(k_sa),dbl_mb(k_sb),nbf,nbf)
CRK         call ga_put(g_ncoul,ad1+1,ad1+nbf,ad2+1,ad2+nbf,
CRK  &                  dbl_mb(k_sb),nbf)
c
CRK         call ga_get(g_exch,1,nbf*nbf,ij,ij,dbl_mb(k_sa),nbf*nbf)
CRK         call ga_put(g_nexch,ad2+1,ad2+nbf,ad1+1,ad1+nbf,
CRK  &                  dbl_mb(k_sa),nbf)
CRK         call mxtrn(dbl_mb(k_sa),dbl_mb(k_sb),nbf,nbf)
CRK         call ga_put(g_nexch,ad1+1,ad1+nbf,ad2+1,ad2+nbf,
CRK  &                  dbl_mb(k_sb),nbf)
CRK       endif
CRK     enddo
CRK   enddo
c
      call ccsd_pampt3(basis,ncor,nocc,nvir,ndel,nbf,nsh,
     &                 tol2e,cmo,eorb,iprt,hiu,giu,habe,gabe,hia,
     &                 dbl_mb(k_t1),ht1,idiis,g_nt2,g_nht2,
     &                 g_ncoul,g_nexch,tklst)
c
      call ga_distribution(g_ht2,iam,g_jlo,g_jhi,g_ilo,g_ihi)
      do a=1,nvir
        if (a.ge.g_ilo.and.a.le.g_ihi)then
          do i=1,nocc
            ad1=(i-1)*nvir+a
            call ga_get(g_nht2,1,lnov,ad1,ad1,dbl_mb(k_sb),lnov)
            do j=1,nocc
              do b=1,nvir
                ad1=k_sa+(b-1)*lnoo+(i-1)*nocc+j-1
                ad2=k_sb+(j-1)*nvir+b-1
                dbl_mb(ad1)=dbl_mb(ad2)
              enddo
            enddo
          enddo
          call ga_acc(g_ht2,offt2+1,offt2+lnoov,a,a,
     &                 dbl_mb(k_sa),lnoov,1.0d00)
        endif
      enddo
CRK   if (.not.ga_destroy(g_nexch))
CRK  &    call errquit('ga_dest g_nexch fail',0)
CRK   if (.not.ga_destroy(g_ncoul))
CRK  &    call errquit('ga_dest g_ncoul fail',0)
      if (.not.ga_destroy(g_nht2))
     &    call errquit('ga_dest g_nht2 fail',0)
      if (.not.ga_destroy(g_nt2))
     &    call errquit('ga_dest g_nt2 fail',0)
      stat=.true.
      stat=stat.and.ma_pop_stack(l_sb)
      stat=stat.and.ma_pop_stack(l_sa)
      stat=stat.and.ma_pop_stack(l_t1)
      if (.not.stat)call errquit('ma_pop i_t_f ',0)
c
      return
      end
