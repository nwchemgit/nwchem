      subroutine ccsd_zsig2(nocc,nvir,nbf,hiu,giu,habe,gabe,hia,t1,
     &                     z1,scra,scrb,scrc,scrd,scre,scrf,
     &                     g_nt2,g_nz2,g_nhz2,g_ncoul,g_nexch,g_x,g_c,
     &                     g_sht2,g_zx,g_zc,g_aitm,g_tz3,iprt)
      implicit none
      integer nocc,nvir,nbf,g_nt2,g_nz2,g_nhz2,g_ncoul,
     &        g_nexch,g_x,g_c,g_sht2,g_zx,g_zc,g_tz3,iprt
      double precision hiu(nocc,nocc),giu(nocc,nocc),hia(nocc,nvir)
      double precision habe(nvir,nvir),gabe(nvir,nvir)
      double precision t1(nocc,nvir),z1(nocc,nvir)
      double precision scra(nbf*nbf),scrb(nbf*nbf),
     &                 scrc(nbf,nbf),scrd(nbf,nbf),
     &                 scre(nbf*nbf),scrf(nbf*nbf)
c
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "bas.fh"
#include "ccsd_debug.fh"
c
      integer g_aitm,g_tmt1,g_tmt2,g_tmt,g_tmp,g_tmp2,
     &        g_jlo,g_jhi,g_ilo,g_ihi,
     &        i,j,k,l,ij,a,b,c,d,ad1,ad2,ad3,ad4,ad5,ad6,lnoo,lnov,lnvv
      integer inode,next,aold
      double precision xdot,ddot
c
      Integer Nodes, IAm
      Nodes = GA_NNodes()
      IAM = GA_NodeID()
C
      lnoo=nocc*nocc
      lnov=nocc*nvir
      lnvv=nvir*nvir
c
      if (.not.ga_create(MT_DBL,lnvv,lnoo,'tmp',
     &                   lnvv,1,g_tmp))
     &     call errquit('ga_create g_tmp failed',0)
      if (.not.ga_create(MT_DBL,lnvv,lnoo,'tmp2',
     &                   lnvv,1,g_tmp2))
     &     call errquit('ga_create g_tmp2 failed',0)
c
      call ga_distribution(g_nexch,iam,g_jlo,g_jhi,g_ilo,g_ihi)
      do i=1,nocc
        ad1=(i-1)*nbf
        if (ad1+1.ge.g_ilo.and.ad1+1.le.g_ihi)then
          do j=1,nocc
            ad2=(j-1)*nbf
            if (ad2+1.ge.g_jlo.and.ad2+1.le.g_jhi)then
            call ga_get(g_nexch,ad2+1,ad2+nbf,ad1+1,ad1+nbf,
     &                  scrc,nbf)

            ad3=0
            do a=1,nvir
            do b=1,nvir
              ad3=ad3+1
              scra(ad3)=-0.5d0*(2.0d0*scrc(nocc+b,nocc+a)
     &                               -scrc(nocc+a,nocc+b))
            enddo
            enddo
            ad3=(i-1)*nvir
            ad4=(j-1)*nvir
            call ga_acc(g_nhz2,ad4+1,ad4+nvir,ad3+1,ad3+nvir,
     &                    scra,nvir,1.0d0)
            endif
          enddo
        endif
      enddo
c
      call ga_distribution(g_nhz2,iam,g_jlo,g_jhi,g_ilo,g_ihi)
      do i=1,nocc
        ad1=(i-1)*nvir
        if (ad1+1.ge.g_ilo.and.ad1+1.le.g_ihi)then
          do j=1,nocc
            ad2=(j-1)*nvir
            if (ad2+1.ge.g_jlo.and.ad2+1.le.g_jhi)then
              call ga_get(g_nz2,ad2+1,ad2+nvir,ad1+1,ad1+nvir,
     &                  scre,nvir)
               ad3=0
               do a=1,nvir
               do b=1,nvir
                 ad3=ad3+1
                 scra(ad3)=-hia(i,a)*z1(j,b)+0.5d0*hia(i,b)*z1(j,a)
               enddo
               enddo
              call dgemm('t','n',nvir,nvir,nvir,1.0d0,gabe,nvir,
     &                   scre,nvir,1.0d0,scra,nvir)
            call ga_acc(g_nhz2,ad2+1,ad2+nvir,ad1+1,ad1+nvir,
     &                  scra,nvir,1.0d00)
            endif
          enddo
        endif
      enddo
c
      call ga_distribution(g_zx,iam,g_jlo,g_jhi,g_ilo,g_ihi)
      do i=1,nocc
      do j=1,nocc
        ij=(i-1)*nocc+j
        if (ij.ge.g_ilo.and.ij.le.g_ihi)then
c
          ad1=(i-1)*nocc+j
          call ga_get(g_zx,1,nbf*nbf,ad1,ad1,scra,nbf*nbf)
            ad1=0
            do a=nocc+1,nbf
            do b=nocc+1,nbf
               ad1=ad1+1
               scrb(ad1)=-0.5d0*(2.0d0*scra((b-1)*nbf+a)
     &                   -scra((a-1)*nbf+b))
            enddo
            enddo
            ad3=(i-1)*nvir
            ad4=(j-1)*nvir
            call ga_acc(g_nhz2,ad3+1,ad3+nvir,ad4+1,ad4+nvir,
     &                    scrb,nvir,1.0d00)
        endif
      enddo
      enddo
c
      call ga_distribution(g_nz2,iam,g_jlo,g_jhi,g_ilo,g_ihi)
      do i=1,nocc
        ad1=(i-1)*nvir
        if (ad1+1.ge.g_ilo.and.ad1+1.le.g_ihi)then
          do j=1,nocc
            ad2=(j-1)*nvir
            if (ad2+1.ge.g_jlo.and.ad2+1.le.g_jhi)then
              ij=(i-1)*nocc+j
              call ga_get(g_nz2,ad2+1,ad2+nvir,ad1+1,ad1+nvir,
     &                    scra,nvir)
              call ga_put(g_tmp,1,lnvv,ij,ij,scra,lnvv)
            endif
          enddo
        endif
      enddo
c
      call ga_dgemm('n','t',lnvv,lnoo,lnoo,
     &              1.0d00,g_tmp,g_aitm,0.0d00,g_tmp2)
c
      call ga_distribution(g_nhz2,iam,g_jlo,g_jhi,g_ilo,g_ihi)
      do i=1,nocc
        ad1=(i-1)*nvir
        if (ad1+1.ge.g_ilo.and.ad1+1.le.g_ihi)then
          do j=1,nocc
            ad2=(j-1)*nvir
            if (ad2+1.ge.g_jlo.and.ad2+1.le.g_jhi)then
              ij=(i-1)*nocc+j
              call ga_get(g_tmp2,1,lnvv,ij,ij,scra,lnvv)
              call ga_acc(g_nhz2,ad2+1,ad2+nvir,ad1+1,ad1+nvir,
     &                    scra,nvir,1.0d00)
            endif
          enddo
        endif
      enddo
c
      inode=-1
      next=nxtval(-nodes)
      next=nxtval(nodes)
      do a=1,nvir
        do i=1,nocc
          ad1=(i-1)*nvir+a
          inode=inode+1
          if (inode.eq.next)then
            call ga_get(g_nz2,1,lnov,ad1,ad1,scre,1)
            call dgemm('n','t',nvir,nocc,nocc,-1.0d00,scre,
     &                 nvir,giu,nocc,0.0d00,scrf,nvir)
            call ga_acc(g_nhz2,1,lnov,ad1,ad1,scrf,lnov,1.0d00)
            next=nxtval(nodes)
          endif
        enddo
      enddo
c
      if (.not.ga_destroy(g_tmp2))
     &    call errquit('ga_dest g_tmp2 fail',0)
      if (.not.ga_destroy(g_tmp))
     &    call errquit('ga_dest g_tmp fail',0)
c
c-----------------------------------------------------------------------
c  J1, J2, J3 and K1, K2, K3 intermediates
c-----------------------------------------------------------------------
c
      if (.not.ga_create(MT_DBL,lnov,lnov,'tmt',
     &                   nvir,nvir,g_tmt))
     &     call errquit('ga_create g_tmt failed',0)
      if (.not.ga_create(MT_DBL,lnov,lnov,'tmp',
     &                   nvir,nvir,g_tmp))
     &     call errquit('ga_create g_tmp failed',0)
      if (.not.ga_create(MT_DBL,lnov,lnov,'tmp2',
     &                   nvir,nvir,g_tmp2))
     &     call errquit('ga_create g_tmp2 failed',0)
c
      call ga_distribution(g_x,iam,g_jlo,g_jhi,g_ilo,g_ihi)
      do i=1,nocc
        do j=1,nocc
          ij=(i-1)*nocc+j
          if (ij.ge.g_ilo.and.ij.le.g_ihi)then
            ad1=(i-1)*nbf
            ad2=(j-1)*nbf
c
            call ga_get(g_nexch,ad2+1,ad2+nbf,ad1+1,ad1+nbf,
     &                  scrd,nbf)
            call dgemm('n','n',nvir,nvir,nocc,-1.0d00,scrd(nocc+1,1),
     &                 nbf,t1,nocc,0.0d00,scrc(nocc+1,nocc+1),nbf)
            call daxpy(nbf*nbf,1.0d00,scrc,1,scrd,1)
            call ga_acc(g_x,1,nbf*nbf,ij,ij,scrd,nbf*nbf,1.0d00)
          endif
        enddo
      enddo
c
      call ga_distribution(g_c,iam,g_jlo,g_jhi,g_ilo,g_ihi)
      do i=1,nocc
        do j=1,nocc
          ij=(i-1)*nocc+j
          if (ij.ge.g_ilo.and.ij.le.g_ihi)then
            ad1=(i-1)*nbf
            ad2=(j-1)*nbf
            call ga_get(g_ncoul,ad2+1,ad2+nbf,ad1+1,ad1+nbf,
     &                    scrd,nbf)
            call dgemm('n','n',nvir,nvir,nocc,-1.0d00,scrd(nocc+1,1),
     &                 nbf,t1,nocc,0.0d00,scrc(nocc+1,nocc+1),nbf)
            call daxpy(nbf*nbf,1.0d00,scrc,1,scrd,1)
            call ga_acc(g_c,1,nbf*nbf,ij,ij,scrd,nbf*nbf,1.0d00)
          endif
        enddo
      enddo
c
c-----------------------------------------------------------------------
c  K4 intermediate
c-----------------------------------------------------------------------
c
c construct (t2_ij^ab + 2*t_i^a*t_j^b)
c
      call ga_distribution(g_nt2,iam,g_jlo,g_jhi,g_ilo,g_ihi)
      do i=1,nocc
        ad1=(i-1)*nvir
        if (ad1+1.ge.g_ilo.and.ad1+1.le.g_ihi)then
          do j=1,nocc
            ad2=(j-1)*nvir
            if (ad2+1.ge.g_jlo.and.ad2+1.le.g_jhi)then
              call ga_get(g_nt2,ad2+1,ad2+nvir,ad1+1,ad1+nvir,
     &                    scra,nvir)
              do a=1,nvir
                do b=1,nvir
                ad3=(a-1)*nvir+b
                scra(ad3)=scra(ad3)+2.0d0*t1(i,a)*t1(j,b)
                enddo
              enddo
              call ga_put(g_tmt,ad1+1,ad1+nvir,ad2+1,ad2+nvir,
     &                     scra,nvir)
            endif
          enddo
        endif
      enddo
c
      call ga_distribution(g_nexch,iam,g_jlo,g_jhi,g_ilo,g_ihi)
      do i=1,nocc
        do j=1,nocc
          ij=(i-1)*nocc+j
          if (ij.ge.g_ilo.and.ij.le.g_ihi)then
            ad1=(i-1)*nbf
            ad2=(j-1)*nbf
            ad3=(i-1)*nvir
            ad4=(j-1)*nvir
            call ga_get(g_nexch,ad2+nocc+1,ad2+nbf,ad1+nocc+1,ad1+nbf,
     &                  scra,nvir)
            call ga_put(g_tmp,ad3+1,ad3+nvir,ad4+1,ad4+nvir,
     &                  scra,nvir)
          endif
        enddo
      enddo
c
      call ga_dgemm('n','n',lnov,lnov,lnov,
     &              -0.5d00,g_tmp,g_tmt,0.0d00,g_tmp2)
c
c-----------------------------------------------------------------------
c product of K intermediate and T2
c-----------------------------------------------------------------------
c
      call ga_distribution(g_c,iam,g_jlo,g_jhi,g_ilo,g_ihi)
      do i=1,nocc
        do j=1,nocc
          ij=(i-1)*nocc+j
          if (ij.ge.g_ilo.and.ij.le.g_ihi)then
            ad1=(i-1)*nbf
            ad2=(j-1)*nbf
            call ga_get(g_c,1,nbf*nbf,ij,ij,scrb,nbf*nbf)
            ad1=0
            do a=nocc+1,nbf
              do b=nocc+1,nbf
                ad1=ad1+1
                scra(ad1)=scrb((a-1)*nbf+b)
              enddo
            enddo
            ad1=(j-1)*nvir
            ad2=(i-1)*nvir
            call ga_acc(g_tmp2,ad1+1,ad1+nvir,ad2+1,ad2+nvir,
     &                  scra,nvir,1.0d00)
          endif
        enddo
      enddo
c
      call ga_dgemm('n','n',lnov,lnov,lnov,
     &              1.0d00,g_tmp2,g_nz2,0.0d00,g_tmp)
c
      call ga_distribution(g_tmp,iam,g_jlo,g_jhi,g_ilo,g_ihi)
      do i=1,nocc
        ad1=(i-1)*nvir
        if (ad1+1.ge.g_ilo.and.ad1+1.le.g_ihi)then
          do j=1,nocc
            ad2=(j-1)*nvir
            if (ad2+1.ge.g_jlo.and.ad2+1.le.g_jhi)then
              call ga_get(g_tmp,ad2+1,ad2+nvir,ad1+1,ad1+nvir,
     &                    scra,nvir)
              do a=1,nvir
                do b=1,nvir
                ad3=(a-1)*nvir+b
                ad4=(b-1)*nvir+a
crk             scrb(ad3)=-0.5d0*scra(ad3)-scra(ad4)
                scrb(ad3)=-0.5d0*scra(ad3)
                enddo
              enddo
c             call ga_acc(g_nhz2,ad2+1,ad2+nvir,ad1+1,ad1+nvir,
c    &                    scrb,nvir,1.0d00)
            endif
          enddo
        endif
      enddo
c
      if (.not.ga_destroy(g_tmp2))
     &    call errquit('ga_dest g_tmp2 fail',0)
      if (.not.ga_destroy(g_tmp))
     &    call errquit('ga_dest g_tmp fail',0)
      if (.not.ga_destroy(g_tmt))
     &    call errquit('ga_dest g_tmt fail',0)
      if (.not.ga_destroy(g_aitm))
     &    call errquit('ga_dest g_aitm fail',0)
      return
      end
