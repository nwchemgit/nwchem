      subroutine ccsd_zsig2(nocc,nvir,nbf,hiu,giu,habe,gabe,hia,t1,
     &                     z1,tzvv,tzoo,
     &                     scra,scrb,scrc,scrd,scre,scrf,
     &                     g_nt2,g_nz2,g_nhz2,g_ncoul,g_nexch,g_x,g_c,
     &                     g_sht2,g_zx,g_zc,g_aitm,g_tz1,g_tz2,g_tz3,
     &                     iprt)
C     $Id: ccsd_zsig2.F,v 2.4 1997-10-22 16:40:56 d3g681 Exp $
      implicit none
      integer nocc,nvir,nbf,g_nt2,g_nz2,g_nhz2,g_ncoul,g_nexch,g_x,g_c,
     &        g_sht2,g_zx,g_zc,g_tz1,g_tz2,g_tz3,iprt
      double precision hiu(nocc,nocc),giu(nocc,nocc),hia(nocc,nvir)
      double precision habe(nvir,nvir),gabe(nvir,nvir)
      double precision t1(nocc,nvir),z1(nocc,nvir)
      double precision tzoo(nocc,nocc),tzvv(nvir,nvir)
      double precision scra(nbf*nbf),scrb(nbf*nbf),
     &                 scrc(nbf,nbf),scrd(nbf,nbf),
     &                 scre(nbf*nbf),scrf(nbf*nbf)
c
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "ccsd_debug.fh"
c
      integer g_aitm,g_tmp,g_tmp2,g_tmp3,
     &        g_jlo,g_jhi,g_ilo,g_ihi,
     &        i,j,k,ij,a,b,c,ad1,ad2,ad3,ad4,ad5,ad6,lnoo,lnov,lnvv
      integer inode,next
      integer nxtval
      external_nxtval
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
      do a=1,nvir
      do c=1,nvir
        scrc(c,a)=tzvv(c,a)
        do k=1,nocc
          scrc(c,a)=scrc(c,a)-0.5d0*t1(k,c)*z1(k,a)
        enddo
      enddo
      enddo
c
      call ga_distribution(g_nexch,iam,g_jlo,g_jhi,g_ilo,g_ihi)
      do i=1,nocc
        ad1=(i-1)*nbf
        if (ad1+1.ge.g_ilo.and.ad1+1.le.g_ihi)then
          do j=1,nocc
            ad2=(j-1)*nbf
            if (ad2+1.ge.g_jlo.and.ad2+1.le.g_jhi)then
            call ga_get(g_nexch,ad2+1,ad2+nbf,ad1+1,ad1+nbf,
     &                  scrd,nbf)
            ad3=0
            do a=1,nvir
            do b=1,nvir
              ad3=ad3+1
              scre(ad3)=-0.5d0*(2.0d0*scrd(nocc+b,nocc+a)
     &                               -scrd(nocc+a,nocc+b))
            do c=1,nvir
              scre(ad3)=scre(ad3)-(2.0d0*scrd(nocc+b,nocc+c)
     &                               -scrd(nocc+c,nocc+b))*scrc(c,a)
            enddo
            do k=1,nocc
              scre(ad3)=scre(ad3)+0.5d0*(2.0d0*scrd(nocc+b,k)
     &                               -scrd(k,nocc+b))*z1(k,a)
            enddo
            enddo
            enddo
            ad3=(i-1)*nvir
            ad4=(j-1)*nvir
            call ga_acc(g_nhz2,ad4+1,ad4+nvir,ad3+1,ad3+nvir,
     &                    scre,nvir,1.0d0)
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
      call ga_distribution(g_nexch,iam,g_jlo,g_jhi,g_ilo,g_ihi)
      do i=1,nocc
        ad1=(i-1)*nbf
        if (ad1+1.ge.g_ilo.and.ad1+1.le.g_ihi)then
          do j=1,nocc
            ad2=(j-1)*nbf
            if (ad2+1.ge.g_jlo.and.ad2+1.le.g_jhi)then
              ij=(i-1)*nocc+j
              call ga_get(g_nexch,ad2+1,ad2+nbf,ad1+1,ad1+nbf,
     &                    scra,nbf)
              ad3=0
              do a=nocc+1,nbf
              do b=nocc+1,nbf
                ad3=ad3+1
                scrb(ad3)=scra((a-1)*nbf+b)
              enddo
              enddo
              call ga_put(g_tmp,1,lnvv,ij,ij,scrb,lnvv)
            endif
          enddo
        endif
      enddo
c  
      call ga_dgemm('n','t',lnvv,lnoo,lnoo,
     &              0.5d00,g_tmp,g_tz3,0.0d00,g_tmp2)
c
      call ga_distribution(g_tmp2,iam,g_jlo,g_jhi,g_ilo,g_ihi)
      do i=1,nocc
      do j=1,nocc
        ij=(i-1)*nocc+j
        if (ij.ge.g_ilo.and.ij.le.g_ihi)then
          ad1=(i-1)*nvir
          ad2=(j-1)*nvir
          call ga_get(g_tmp2,1,lnvv,ij,ij,scra,nvir)
          call ga_acc(g_nhz2,ad2+1,ad2+nvir,ad1+1,ad1+nvir,
     &                scra,nvir,1.0d00)
        endif
      enddo
      enddo
c
      do j=1,nocc
      do k=1,nocc
        scrc(k,j)=tzoo(k,j)
        do c=1,nvir
          scrc(k,j)=scrc(k,j)-0.5d0*t1(k,c)*z1(j,c)
        enddo
      enddo
      enddo
c
      call ga_distribution(g_nexch,iam,g_jlo,g_jhi,g_ilo,g_ihi)
      do i=1,nocc
        ad1=(i-1)*nbf
        if (ad1+1.ge.g_ilo.and.ad1+1.le.g_ihi)then
          do j=1,nocc
            ad2=(j-1)*nbf
            if (ad2+1.ge.g_jlo.and.ad2+1.le.g_jhi)then
              call dfill(nvir*nvir,0.0d00,scrf,1)
              do k=1,nocc
                ad3=(k-1)*nbf
                call ga_get(g_nexch,ad3+1,ad3+nbf,ad1+1,ad1+nbf,
     &                      scre,nbf)
                do a=1,nvir
                do b=1,nvir
                  ad4=(a-1)*nvir+b
                  ad5=(nocc+a-1)*nbf+nocc+b
                  ad6=(nocc+b-1)*nbf+nocc+a
                  scrf(ad4)=scrf(ad4)-(2.0d0*scre(ad5)-scre(ad6))
     &                                *scrc(k,j)
                enddo
                enddo
              enddo
              ad5=(i-1)*nvir
              ad6=(j-1)*nvir
              call ga_acc(g_nhz2,ad6+1,ad6+nvir,ad5+1,ad5+nvir,
     &                    scrf,nvir,1.0d00)
            endif
          enddo
        endif
      enddo
c
      inode=-1
      call ga_sync
      next=nxtval(-nodes)
      call ga_sync
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
      if (.not.ga_create(MT_DBL,lnov,lnov,'tmp',
     &                   nvir,nvir,g_tmp))
     &     call errquit('ga_create g_tmp failed',0)
      if (.not.ga_create(MT_DBL,lnov,lnov,'tmp2',
     &                   nvir,nvir,g_tmp2))
     &     call errquit('ga_create g_tmp2 failed',0)
      if (.not.ga_create(MT_DBL,lnov,lnov,'tmp3',
     &                   nvir,nvir,g_tmp3))
     &     call errquit('ga_create g_tmp3 failed',0)
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
c
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
              call ga_put(g_tmp,ad1+1,ad1+nvir,ad2+1,ad2+nvir,
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
            call ga_put(g_tmp2,ad3+1,ad3+nvir,ad4+1,ad4+nvir,
     &                  scra,nvir)
          endif
        enddo
      enddo
c
      call ga_dgemm('n','n',lnov,lnov,lnov,
     &              -0.5d00,g_tmp2,g_tmp,0.0d00,g_tmp3)
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
            call ga_acc(g_tmp3,ad1+1,ad1+nvir,ad2+1,ad2+nvir,
     &                  scra,nvir,1.0d00)
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
              call ga_get(g_nz2,ad2+1,ad2+nvir,ad1+1,ad1+nvir,
     &                    scra,nvir)
              do a=1,nvir
              do b=1,nvir
                 ad3=(b-1)*nvir+a
                 ad4=(a-1)*nvir+b
                 scrb(ad4)=0.5d0*scra(ad4)+scra(ad3)
              enddo
              enddo
              call ga_put(g_tmp,ad2+1,ad2+nvir,ad1+1,ad1+nvir,
     &                    scrb,nvir)
            endif
          enddo
        endif
      enddo
      call ga_dgemm('n','n',lnov,lnov,lnov,
     &             -1.0d00,g_tmp3,g_tmp,0.0d00,g_tmp2)
c
      call ga_distribution(g_tmp2,iam,g_ilo,g_ihi,g_jlo,g_jhi)
      do i=1,nocc
        ad1=(i-1)*nvir
        if (ad1+1.ge.g_ilo.and.ad1+1.le.g_ihi)then
          do j=1,nocc
            ad2=(j-1)*nvir
            if (ad2+1.ge.g_jlo.and.ad2+1.le.g_jhi)then
              call ga_get(g_tmp2,ad1+1,ad1+nvir,ad2+1,ad2+nvir,
     &                    scra,nvir)
              do a=1,nvir
                do b=1,nvir
                ad3=(a-1)*nvir+b
                ad4=(b-1)*nvir+a
                scrb(ad3)=scra(ad4)
                enddo
              enddo
              call ga_acc(g_nhz2,ad1+1,ad1+nvir,ad2+1,ad2+nvir,
     &                    scrb,nvir,1.0d0)
            endif
          enddo
        endif
      enddo
c
      if (iprt.gt.5.and.iam.eq.0)print *,'forming jmk'
      call ga_distribution(g_x,iam,g_jlo,g_jhi,g_ilo,g_ihi)
      do i=1,nocc
        do j=1,nocc
          ij=(i-1)*nocc+j
          if (ij.ge.g_ilo.and.ij.le.g_ihi)then
            ad1=(i-1)*nbf
            ad2=(j-1)*nbf
            call dfill(nbf*nbf,0.0d00,scra,1)
              call ga_get(g_x,1,nbf*nbf,ij,ij,scrb,nbf*nbf)
              call daxpy(nbf*nbf,1.0d00,scrb,1,scra,1)
              call ga_get(g_c,1,nbf*nbf,ij,ij,scrb,nbf*nbf)
              call daxpy(nbf*nbf,-0.5d00,scrb,1,scra,1)
            ad1=0
            do a=nocc+1,nbf
              do b=nocc+1,nbf
                ad1=ad1+1
                scrb(ad1)=scra((a-1)*nbf+b)
              enddo
            enddo
            ad1=(j-1)*nvir
            ad2=(i-1)*nvir
            call ga_put(g_tmp,ad1+1,ad1+nvir,ad2+1,ad2+nvir,
     &                  scrb,nvir)
          endif
        enddo
      enddo
c-----------------------------------------------------------------------
c  J4 + J5 - K4 intermediate
c-----------------------------------------------------------------------
c
c construct (2*t2_ij^ab - t2_ij^ba - 2*t_i^a t_j^b)
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
                ad3=(b-1)*nvir+a
                ad4=(a-1)*nvir+b
                scrb(ad3)=scra(ad3)+scra(ad3)-scra(ad4)
     &                    -2.0d00*t1(i,a)*t1(j,b)
                enddo
              enddo
              call ga_put(g_tmp3,ad2+1,ad2+nvir,ad1+1,ad1+nvir,
     &                     scrb,nvir)
            endif
          enddo
        endif
      enddo
c
      call ga_distribution(g_nexch,iam,g_jlo,g_jhi,g_ilo,g_ihi)
      do i=1,nocc
        ad1=(i-1)*nbf
        if (ad1+1.ge.g_ilo.and.ad1+1.le.g_ihi)then
          do j=1,nocc
            ad2=(j-1)*nbf
            if (ad2+1.ge.g_jlo.and.ad2+1.le.g_jhi)then
            ad3=(i-1)*nvir
            ad4=(j-1)*nvir
            call ga_get(g_nexch,ad2+nocc+1,ad2+nbf,ad1+nocc+1,ad1+nbf,
     &                  scra,nvir)
              do a=1,nvir
                do b=1,nvir
                ad5=(a-1)*nvir+b
                ad6=(b-1)*nvir+a
                scrb(ad5)=scra(ad5)+scra(ad5)-scra(ad6)
                enddo
              enddo
            call ga_put(g_tmp2,ad4+1,ad4+nvir,ad3+1,ad3+nvir,
     &                  scrb,nvir)
          endif
          enddo
        endif
      enddo
c
      call ga_dgemm('n','n',lnov,lnov,lnov,
     &              0.25d00,g_tmp2,g_tmp3,1.0d00,g_tmp)
      call ga_dgemm('n','n',lnov,lnov,lnov,
     &              1.0d00,g_tmp,g_nz2,0.0d00,g_tmp2)
c
      call ga_distribution(g_tmp2,iam,g_ilo,g_ihi,g_jlo,g_jhi)
      do i=1,nocc
        ad1=(i-1)*nvir
        if (ad1+1.ge.g_ilo.and.ad1+1.le.g_ihi)then
          do j=1,nocc
            ad2=(j-1)*nvir
            if (ad2+1.ge.g_jlo.and.ad2+1.le.g_jhi)then
              call ga_get(g_tmp2,ad1+1,ad1+nvir,ad2+1,ad2+nvir,
     &                    scra,nvir)
              do a=1,nvir
                do b=1,nvir
                ad3=(a-1)*nvir+b
                ad4=(b-1)*nvir+a
                scrb(ad3)=scra(ad3)+scra(ad3)-scra(ad4)
                enddo
              enddo
              call ga_acc(g_nhz2,ad1+1,ad1+nvir,ad2+1,ad2+nvir,
     &                    scrb,nvir,1.0d00)
            endif
          enddo
        endif
      enddo
c
      call ga_distribution(g_nexch,iam,g_jlo,g_jhi,g_ilo,g_ihi)
      do i=1,nocc
        ad1=(i-1)*nbf
        if (ad1+1.ge.g_ilo.and.ad1+1.le.g_ihi)then
          do j=1,nocc
            ad2=(j-1)*nbf
            if (ad2+1.ge.g_jlo.and.ad2+1.le.g_jhi)then
              call ga_get(g_nexch,ad2+1,ad2+nbf,ad1+1,ad1+nbf,
     &                    scra,nbf)
              do a=1,nvir
              do b=1,nvir
                ad3=(nocc+b-1)*nbf+nocc+a
                ad4=(nocc+a-1)*nbf+nocc+b
                scrb((b-1)*nvir+a)=2.0d0*scra(ad3)-scra(ad4)
                scre((b-1)*nvir+a)=scra(ad3)
                scrf((b-1)*nvir+a)=scra(ad4)
              enddo
              enddo
              ad5=(i-1)*nvir
              ad6=(j-1)*nvir
              call ga_put(g_tmp,ad6+1,ad6+nvir,ad5+1,ad5+nvir,
     &                    scrb,nvir)
              call ga_put(g_tmp2,ad6+1,ad6+nvir,ad5+1,ad5+nvir,
     &                    scre,nvir)
              call ga_put(g_tmp3,ad6+1,ad6+nvir,ad5+1,ad5+nvir,
     &                    scrf,nvir)
            endif
          enddo
        endif
      enddo
c
      call ga_dgemm('n','n',lnov,lnov,lnov,
     &              0.5d00,g_tmp,g_tz1,1.0d00,g_nhz2)
      call ga_dgemm('n','n',lnov,lnov,lnov,
     &              -0.5d00,g_tmp2,g_tz1,0.0d00,g_tmp)
      call ga_dgemm('n','n',lnov,lnov,lnov,
     &              0.5d00,g_tmp3,g_tz2,0.0d00,g_tmp2)
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
              call ga_get(g_tmp2,ad2+1,ad2+nvir,ad1+1,ad1+nvir,
     &                    scrb,nvir)
              do a=1,nvir
              do b=1,nvir
                scrb((b-1)*nvir+a)=scrb((b-1)*nvir+a)+scra((b-1)*nvir+a)
              enddo
              enddo
              call ga_acc(g_nhz2,ad1+1,ad1+nvir,ad2+1,ad2+nvir,
     &                    scrb,nvir,1.0d00)
            endif
          enddo
        endif
      enddo
c
      if (.not.ga_destroy(g_tmp3))
     &    call errquit('ga_dest g_tmp3 fail',0)
      if (.not.ga_destroy(g_tmp2))
     &    call errquit('ga_dest g_tmp2 fail',0)
      if (.not.ga_destroy(g_tmp))
     &    call errquit('ga_dest g_tmp fail',0)
c
      return
      end
