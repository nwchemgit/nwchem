      subroutine ccsd_zitf(basis,nsh,ncor,nocc,nvir,nact,nbf,
     &                     t1,z1,hiu,giu,habe,gabe,hia,hz1,
     &                    idiis,cmo,eorb,iprt,
     &                    g_t2,g_z2,g_hz2,g_ncoul,g_nexch,tklst)
C     $Id: ccsd_zitf.F,v 2.4 2003-10-17 22:54:25 carlfahl Exp $
      implicit none
#include "errquit.fh"
      integer basis,nsh,ncor,nocc,nvir,nact,nbf,
     &        idiis,iprt,g_t2,g_z2,g_hz2
      double precision cmo(*),eorb(*),t1(*),hiu(*),giu(*),habe(*),
     &                 gabe(*),hz1(*),hia(*),z1(*)
#include "mafdecls.fh"
#include "bas.fh"
#include "global.fh"
#include "ccsd_debug.fh"
c
      integer lnov,g_nt2,g_nz2,g_nhz2,offt2,lnoov,g_jlo,g_jhi,
     &        g_ilo,g_ihi,a,i,j,b,ad1,ad2,lnoo,
     &        l_sa,k_sa,l_sb,k_sb,lsab,
     &        g_ncoul,g_nexch,tklst(nsh*(nsh+1)/2,2)
      double precision tol2e
      logical stat
      Integer IAm
c
      IAM = GA_NodeID()
C
      lnoo=nocc*nocc
      lnov=nocc*nvir
      lnoov=nocc*nocc*nvir
      tol2e=1.0d-12
      offt2=(idiis-1)*lnoov
c
      lsab=max(lnoov,nbf*nbf)
      stat=.true.
      stat=stat.and.ma_push_get(MT_DBL,lsab,'sa',l_sa, k_sa)
      stat=stat.and.ma_push_get(MT_DBL,lsab,'sb',l_sb, k_sb)
      if (.not.stat)call errquit('ma_push zitf ',0, MA_ERR)
c
c ------------------------------------------------------------
c create new GAs with proposed final ordering
c ------------------------------------------------------------
      if (.not.ga_create(MT_DBL,lnov,lnov,'nt2',
     &                   nvir,nvir,g_nt2))
     &     call errquit('ga_create g_nt2 failed',0, GA_ERR)
      if (.not.ga_create(MT_DBL,lnov,lnov,'nz2',
     &                   nvir,nvir,g_nz2))
     &     call errquit('ga_create g_nz2 failed',0, GA_ERR)
      if (.not.ga_create(MT_DBL,lnov,lnov,'nhz2',
     &                   nvir,nvir,g_nhz2))
     &     call errquit('ga_create g_nhz2 failed',0, GA_ERR)

c ------------------------------------------------------------
c fill new GAs
c ------------------------------------------------------------
      call ga_distribution(g_t2,iam,g_jlo,g_jhi,g_ilo,g_ihi)
      do a=1,nvir
        if (a.ge.g_ilo.and.a.le.g_ihi)then
          call ga_get(g_t2,1,lnoov,a,a,
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
      call ga_distribution(g_z2,iam,g_jlo,g_jhi,g_ilo,g_ihi)
      do a=1,nvir
        if (a.ge.g_ilo.and.a.le.g_ihi)then
          call ga_get(g_z2,offt2+1,offt2+lnoov,a,a,
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
            call ga_put(g_nz2,1,lnov,ad1,ad1,dbl_mb(k_sb),lnov)
          enddo
        endif
      enddo
c
      call ga_zero(g_nhz2)
c
      call ccsd_pzamp(basis,nsh,ncor,nocc,nvir,nact,nbf,
     &                 tol2e,cmo,eorb,iprt,hiu,giu,habe,gabe,hia,
     &                 t1,z1,hz1,idiis,g_nt2,g_nz2,g_nhz2,
     &                 g_ncoul,g_nexch,tklst)
c
      call ga_distribution(g_hz2,iam,g_jlo,g_jhi,g_ilo,g_ihi)
      do a=1,nvir
        if (a.ge.g_ilo.and.a.le.g_ihi)then
          do i=1,nocc
            ad1=(i-1)*nvir+a
            call ga_get(g_nhz2,1,lnov,ad1,ad1,dbl_mb(k_sb),lnov)
            do j=1,nocc
              do b=1,nvir
                ad1=k_sa+(b-1)*lnoo+(i-1)*nocc+j-1
                ad2=k_sb+(j-1)*nvir+b-1
                dbl_mb(ad1)=dbl_mb(ad2)
              enddo
            enddo
          enddo
          call ga_acc(g_hz2,offt2+1,offt2+lnoov,a,a,
     &                 dbl_mb(k_sa),lnoov,1.0d00)
        endif
      enddo
c
      if (.not.ga_destroy(g_nhz2))
     &    call errquit('ga_dest g_nhz2 fail',0, GA_ERR)
      if (.not.ga_destroy(g_nz2))
     &    call errquit('ga_dest g_nz2 fail',0, GA_ERR)
      if (.not.ga_destroy(g_nt2))
     &    call errquit('ga_dest g_nt2 fail',0, GA_ERR)
      stat=.true.
      stat=stat.and.ma_pop_stack(l_sb)
      stat=stat.and.ma_pop_stack(l_sa)
      if (.not.stat)call errquit('ma_pop zitf ',0, MA_ERR)
c
      return
      end
