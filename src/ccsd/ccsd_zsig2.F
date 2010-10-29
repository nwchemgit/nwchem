      subroutine ccsd_zsig2(nocc,nvir,nbf,hiu,giu,habe,gabe,hia,t1,
     &                     z1,tzvv,tzoo,
     &                     scra,scrb,scrc,scrd,scre,scrf,
     &                     g_nt2,g_nz2,g_nhz2,g_ncoul,g_nexch,g_x,g_c,
     &                     g_sht2,g_zx,g_zc,g_aitm,g_qitm,g_tz3,
     &                     iprt)
C     $Id$
      implicit none
#include "errquit.fh"
      integer nocc,nvir,nbf,g_nt2,g_nz2,g_nhz2,g_ncoul,g_nexch,g_x,g_c,
     &        g_sht2,g_zx,g_zc,g_tz3,iprt
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
      integer g_aitm,g_qitm,g_ritm,g_tmp,g_tmp2,
     &        g_jlo,g_jhi,g_ilo,g_ihi,
     &        i,j,k,l,ij,a,b,c,d,ad1,ad2,ad3,ad4,ad5,ad6,lnoo,lnov,lnvv
      integer inode,next
      integer nxtask
      external nxtask
c
      Integer Nodes, IAm
      Nodes = GA_NNodes()
      IAM = GA_NodeID()
C
      lnoo=nocc*nocc
      lnov=nocc*nvir
      lnvv=nvir*nvir
C
*ga:1:0
      if (.not.ga_create(MT_DBL,lnvv,lnoo,'tmp',
     &                   lnvv,0,g_tmp))
     &     call errquit('ga_create g_tmp failed',0, GA_ERR)
*ga:1:0
      if (.not.ga_create(MT_DBL,lnvv,lnoo,'tmp2',
     &                   lnvv,0,g_tmp2))
     &     call errquit('ga_create g_tmp2 failed',0, GA_ERR)
c
      do a=1,nvir
      do c=1,nvir
        scrc(c,a)=tzvv(c,a)
        do k=1,nocc
          scrc(c,a)=scrc(c,a)+t1(k,c)*z1(k,a)
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
            do a=1,nvir
            do b=1,nvir
              ad3=(a-1)*nvir+b
              scre(ad3)=-(2.0d0*scrd(nocc+b,nocc+a)
     &                               -scrd(nocc+a,nocc+b))
            do c=1,nvir
              scre(ad3)=scre(ad3)+(2.0d0*scrd(nocc+c,nocc+a)
     &                               -scrd(nocc+a,nocc+c))*scrc(c,b)
            enddo
            do k=1,nocc
              scre(ad3)=scre(ad3)+(2.0d0*scrd(k,nocc+a)
     &                               -scrd(nocc+a,k))*z1(k,b)
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
               do a=1,nvir
               do b=1,nvir
               ad3=(a-1)*nvir+b
               scra(ad3)=-2.0d0*hia(i,a)*z1(j,b)+hia(i,b)*z1(j,a)
               enddo
               enddo
              call dgemm('n','n',nvir,nvir,nvir,-1.0d0,scre,nvir,
     &                   gabe,nvir,1.0d0,scra,nvir)
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
               scrb(ad1)=-(2.0d0*scra((b-1)*nbf+a)
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
     &              -1.0d00,g_tmp,g_aitm,0.0d00,g_tmp2)
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
     &              -0.5d00,g_tmp,g_tz3,0.0d00,g_tmp2)
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
          scrc(k,j)=scrc(k,j)+t1(k,c)*z1(j,c)
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
              call dfill(nbf*nbf,0.0d00,scrf,1)
              do k=1,nocc
                ad3=(k-1)*nbf
                call ga_get(g_nexch,ad3+1,ad3+nbf,ad1+1,ad1+nbf,
     &                      scre,nbf)
                do a=1,nvir
                do b=1,nvir
                  ad4=(a-1)*nvir+b
                  ad5=(nocc+a-1)*nbf+nocc+b
                  ad6=(nocc+b-1)*nbf+nocc+a
                  scrf(ad4)=scrf(ad4)+(2.0d0*scre(ad5)-scre(ad6))
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
      next=nxtask(-nodes, 1)
      call ga_sync
      next=nxtask(nodes, 1)
      do a=1,nvir
        do i=1,nocc
          ad1=(i-1)*nvir+a
          inode=inode+1
          if (inode.eq.next)then
            call ga_get(g_nz2,1,lnov,ad1,ad1,scre,1)
            call dgemm('n','t',nvir,nocc,nocc,1.0d00,scre,
     &                 nvir,giu,nocc,0.0d00,scrf,nvir)
            call ga_acc(g_nhz2,1,lnov,ad1,ad1,scrf,lnov,1.0d00)
            next=nxtask(nodes, 1)
          endif
        enddo
      enddo
c
      if (.not.ga_destroy(g_tmp2))
     &    call errquit('ga_dest g_tmp2 fail',0, GA_ERR)
      if (.not.ga_destroy(g_tmp))
     &    call errquit('ga_dest g_tmp fail',0, GA_ERR)
c
      if (.not.ga_create(MT_DBL,lnov,lnov,'ritm',
     &                   nvir,nvir,g_ritm))
     &     call errquit('ga_create g_ritm failed',0, GA_ERR)
      if (.not.ga_create(MT_DBL,lnov,lnov,'tmp',
     &                   nvir,nvir,g_tmp))
     &     call errquit('ga_create g_tmp failed',0, GA_ERR)
      if (.not.ga_create(MT_DBL,lnov,lnov,'tmp2',
     &                   nvir,nvir,g_tmp2))
     &     call errquit('ga_create g_tmp2 failed',0, GA_ERR)
c
      call ga_dgemm('n','n',lnov,lnov,lnov,
     &              1.0d00,g_nz2,g_qitm,0.0d00,g_tmp)
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
              call ga_get(g_tmp,ad1+1,ad1+nvir,ad2+1,ad2+nvir,
     &                    scrb,nvir)
              do a=1,nvir
              do b=1,nvir
                scra((a-1)*nvir+b)=scra((a-1)*nvir+b)
     &                            -0.5d0*scrb((a-1)*nvir+b)
              enddo
              enddo
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
          do k=1,nocc
            ad2=(k-1)*nbf
            if (ad2+1.ge.g_jlo.and.ad2+1.le.g_jhi)then
              call ga_get(g_ncoul,ad2+1,ad2+nbf,ad1+1,ad1+nbf,
     &                    scra,nbf)
              call dfill(lnvv,0.0d0,scre,1)
              do c=1,nvir
              do a=1,nvir
                scre((c-1)*nvir+a)=scre((c-1)*nvir+a)
     &                             +scra((nocc+c-1)*nbf+nocc+a)
              enddo
              enddo
              do c=1,nvir
              do a=1,nvir
              do l=1,nocc
                scre((c-1)*nvir+a)=scre((c-1)*nvir+a)
     &                             -scra((nocc+c-1)*nbf+l)*t1(l,a)
              enddo
              enddo
              enddo
              ad3=(i-1)*nvir
              ad4=(k-1)*nvir
              call ga_put(g_ritm,ad4+1,ad4+nvir,ad3+1,ad3+nvir,
     &                    scre,nvir)
            endif
          enddo
        endif
      enddo
      call ga_sync()
c
      call ga_distribution(g_c,iam,g_jlo,g_jhi,g_ilo,g_ihi)
      do i=1,nocc
      do k=1,nocc
        ad1=(i-1)*nocc+k
        if (ad1.ge.g_ilo.and.ad1.le.g_ihi)then
          call ga_get(g_c,1,nbf*nbf,ad1,ad1,scra,nbf*nbf)
              do c=1,nvir
              do a=1,nvir
                scre((a-1)*nvir+c)=scra((nocc+c-1)*nbf+nocc+a)
              enddo
              enddo
              ad3=(i-1)*nvir
              ad4=(k-1)*nvir
              call ga_acc(g_ritm,ad3+1,ad3+nvir,ad4+1,ad4+nvir,
     &                    scre,nvir,1.0d0)
            endif
          enddo
      enddo
c
      call ga_distribution(g_nexch,iam,g_jlo,g_jhi,g_ilo,g_ihi)
      do k=1,nocc
        ad1=(k-1)*nbf
        if (ad1+1.ge.g_ilo.and.ad1+1.le.g_ihi)then
          do i=1,nocc
            ad2=(i-1)*nvir
            if (ad2+1.ge.g_jlo.and.ad2+1.le.g_jhi)then
              call dfill(lnvv,0.0d0,scre,1)
              do l=1,nocc
              ad3=(l-1)*nbf
              ad4=(l-1)*nvir
              call ga_get(g_nexch,ad3+1,ad3+nbf,ad1+1,ad1+nbf,
     &                    scra,nbf)
              call ga_get(g_nt2,ad4+1,ad4+nvir,ad2+1,ad2+nvir,
     &                    scrb,nvir)
              do c=1,nvir
              do d=1,nvir
              do a=1,nvir
           scre((c-1)*nvir+a)=scre((c-1)*nvir+a)
     &                    -scra((nocc+d-1)*nbf+nocc+c)
     &               *(scrb((d-1)*nvir+a)+t1(i,d)*t1(l,a))
              enddo
              enddo
              enddo
              enddo
              ad3=(k-1)*nvir
              call ga_acc(g_ritm,ad2+1,ad2+nvir,ad3+1,ad3+nvir,
     &                    scre,nvir,1.0d0)
            endif
          enddo
        endif
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
              do a=1,nvir
              do b=1,nvir
              scrb((b-1)*nvir+a)=0.5d0*scra((b-1)*nvir+a)
     &                          +scra((a-1)*nvir+b)
              enddo
              enddo
              call ga_put(g_tmp,ad2+1,ad2+nvir,ad1+1,ad1+nvir,
     &                    scrb,nvir)
            endif
          enddo
        endif
      enddo
c
      call ga_dgemm('n','n',lnov,lnov,lnov,
     &              1.0d00,g_tmp,g_ritm,0.0d00,g_tmp2)
c
      call ga_distribution(g_tmp2,iam,g_jlo,g_jhi,g_ilo,g_ihi)
      do i=1,nocc
        ad1=(i-1)*nvir
        if (ad1+1.ge.g_ilo.and.ad1+1.le.g_ihi)then
          do j=1,nocc
            ad2=(j-1)*nvir
            if (ad2+1.ge.g_jlo.and.ad2+1.le.g_jhi)then
              ij=(i-1)*nocc+j
              call ga_get(g_tmp2,ad2+1,ad2+nvir,ad1+1,ad1+nvir,
     &                    scra,nvir)
              do a=1,nvir
              do b=1,nvir
              scrb((b-1)*nvir+a)=scra((a-1)*nvir+b)
              enddo
              enddo
              call ga_acc(g_nhz2,ad2+1,ad2+nvir,ad1+1,ad1+nvir,
     &                    scrb,nvir,1.0d0)
            endif
          enddo
        endif
      enddo
c
      if (.not.ga_destroy(g_tmp2))
     &    call errquit('ga_dest g_tmp2 fail',0, GA_ERR)
      if (.not.ga_destroy(g_tmp))
     &    call errquit('ga_dest g_tmp fail',0, GA_ERR)
      if (.not.ga_destroy(g_ritm))
     &    call errquit('ga_dest g_ritm fail',0, GA_ERR)
      if (.not.ga_destroy(g_tz3))
     &    call errquit('ga_dest g_tz3 fail',0, GA_ERR)
      if (.not.ga_destroy(g_qitm))
     &    call errquit('ga_dest g_qitm fail',0, GA_ERR)
      if (.not.ga_destroy(g_aitm))
     &    call errquit('ga_dest g_aitm fail',0, GA_ERR)
c
      return
      end
