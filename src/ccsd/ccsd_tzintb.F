      subroutine ccsd_tzintb(nocc,nvir,nbf,t1,tzoo,scra,scrb,scrc,
     &                       g_nt2,g_nz2,g_ncoul,g_nexch,g_c,g_x,
     &                       g_sht2,g_aitm,g_qitm,g_tz3,iprt)
* $Id$
      implicit none
#include "errquit.fh"
      integer nocc,nvir,nbf,lnoo,lnov,lnvv,g_nt2,g_nz2,g_ncoul,g_nexch,
     &        g_c,g_x,g_sht2,g_aitm,g_qitm,g_tz3,g_tmp,iprt
      integer i,j,k,l,a,b,c,d,ij,kl,ad1,ad2,ad3,ad4,
     &        g_jlo,g_jhi,g_ilo,g_ihi
      double precision t1(nocc,nvir),tzoo(nocc,nocc),scra(nbf*nbf),
     &                 scrb(nbf*nbf),scrc(nbf*nbf)
      Integer IAm
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
c
      IAM = GA_NodeID()
C
      lnoo=nocc*nocc
      lnov=nocc*nvir
      lnvv=nvir*nvir
c
      if (.not.ga_create(MT_DBL,lnov,lnov,'tmp',
     &                   nvir,nvir,g_tmp))
     &     call errquit('ga_create g_tmp failed',0, GA_ERR)
      call ga_dgemm('n','t',lnov,lnov,lnov,
     &              1.0d00,g_nt2,g_nz2,0.0d00,g_tmp)
c
      call dfill(lnoo,0.0d00,tzoo,1)
      call ga_distribution(g_tmp,iam,g_jlo,g_jhi,g_ilo,g_ihi)
      do i=1,nocc
        ad1=(i-1)*nvir
        if (ad1+1.ge.g_ilo.and.ad1+1.le.g_ihi)then
          do j=1,nocc
            ad2=(j-1)*nvir
            if (ad2+1.ge.g_jlo.and.ad2+1.le.g_jhi)then
              call ga_get(g_tmp,ad1+1,ad1+nvir,ad2+1,ad2+nvir,
     &                    scra,nvir)
              do a=1,nvir
                tzoo(i,j)=tzoo(i,j)+scra((a-1)*nvir+a)
              enddo
            endif
          enddo
        endif
      enddo
c
      call ga_dgop(912,tzoo,lnoo, '+')
c
      if (.not.ga_destroy(g_tmp))
     &    call errquit('ga_dest g_tmp fail',0, GA_ERR)
c
      if (.not.ga_create(MT_DBL,lnoo,lnoo,'aitm',
     &                   nocc,nocc,g_aitm))
     &     call errquit('ga_create g_aitm failed',0, GA_ERR)
      if (.not.ga_create(MT_DBL,lnov,lnov,'qitm',
     &                   nvir,nvir,g_qitm))
     &     call errquit('ga_create g_qitm failed',0, GA_ERR)
      if (.not.ga_create(MT_DBL,lnoo,lnoo,'tz3',
     &                   nocc,nocc,g_tz3))
     &     call errquit('ga_create g_tz3 failed',0, GA_ERR)
c
      call ga_distribution(g_x,iam,g_jlo,g_jhi,g_ilo,g_ihi)
      do i=1,nocc
        do j=1,nocc
          ij=(i-1)*nocc+j
          if (ij.ge.g_ilo.and.ij.le.g_ihi)then
            ad1=(i-1)*nbf
            ad2=(j-1)*nbf
c
c-----------------------------------------------------------------------
c  A1, A2, A3 and A4 intermediates
c-----------------------------------------------------------------------
c
c A1
            call dfill(nbf*nbf,0.0d00,scra,1)
            call ga_get(g_nexch,ad2+1,ad2+nocc,ad1+1,ad1+nocc,
     &                  scrb,nocc)
            call daxpy(nocc*nocc,0.5d00,scrb,1,scra,1)
            ad1=(i-1)*nocc+j
            call ga_get(g_x,1,nbf*nbf,ad1,ad1,scrb,nbf*nbf)
            ad1=0
            do k=1,nocc
              do l=1,nocc
                ad1=ad1+1
                scra(ad1)=scra(ad1)+0.5d00*scrb((k-1)*nbf+l)
              enddo
            enddo
            ad1=(j-1)*nocc+i
            call ga_get(g_x,1,nbf*nbf,ad1,ad1,scrb,nbf*nbf)
            ad1=0
            do k=1,nocc
              do l=1,nocc
                ad1=ad1+1
                scra(ad1)=scra(ad1)+0.5d00*scrb((l-1)*nbf+k)
              enddo
            enddo
            call ga_get(g_sht2,1,nbf*nbf,ij,ij,scrb,nbf*nbf)
            ad1=0
            do k=1,nocc
              do l=1,nocc
                ad1=ad1+1
                scra(ad1)=scra(ad1)+scrb((k-1)*nbf+l)
              enddo
            enddo
            call ga_put(g_aitm,1,lnoo,ij,ij,scra,lnoo)
          endif
        enddo
      enddo
c
c-----------------------------------------------------------------------
c  Q intermediate
c-----------------------------------------------------------------------
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
              call ga_get(g_nexch,ad2+1,ad2+nbf,ad1+1,ad1+nbf,
     &                    scrb,nbf)
              call dfill(lnvv,0.0d0,scrc,1)
              do c=1,nvir
              do a=1,nvir
                scrc((c-1)*nvir+a)=scrc((c-1)*nvir+a)
     &          -(2.0d0*scrb((nocc+c-1)*nbf+nocc+a)
     &                             -scra((nocc+c-1)*nbf+nocc+a))
              enddo
              enddo
              do c=1,nvir
              do a=1,nvir
              do l=1,nocc
                scrc((c-1)*nvir+a)=scrc((c-1)*nvir+a)
     &          +(2.0d0*scrb((nocc+c-1)*nbf+l)
     &                             -scra((nocc+c-1)*nbf+l))*t1(l,a)
              enddo
              enddo
              enddo
              ad3=(i-1)*nvir
              ad4=(k-1)*nvir
              call ga_put(g_qitm,ad4+1,ad4+nvir,ad3+1,ad3+nvir,
     &                    scrc,nvir)
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
c
          call ga_get(g_c,1,nbf*nbf,ad1,ad1,scra,nbf*nbf)
          call ga_get(g_x,1,nbf*nbf,ad1,ad1,scrb,nbf*nbf)
c
              do c=1,nvir
              do a=1,nvir
                scrc((a-1)*nvir+c)=-(2.0d0*scrb((nocc+c-1)*nbf+nocc+a)
     &                             -scra((nocc+c-1)*nbf+nocc+a))
              enddo
              enddo
              ad3=(i-1)*nvir
              ad4=(k-1)*nvir
              call ga_acc(g_qitm,ad3+1,ad3+nvir,ad4+1,ad4+nvir,
     &                    scrc,nvir,1.0d0)
            endif
          enddo
      enddo
c
      call ga_distribution(g_qitm,iam,g_jlo,g_jhi,g_ilo,g_ihi)
      do k=1,nocc
        ad1=(k-1)*nvir
        if (ad1+1.ge.g_ilo.and.ad1+1.le.g_ihi)then
          do i=1,nocc
            ad2=(i-1)*nvir
            if (ad2+1.ge.g_jlo.and.ad2+1.le.g_jhi)then
              call dfill(lnvv,0.0d0,scrc,1)
              do l=1,nocc
              ad3=(k-1)*nbf
              ad4=(l-1)*nbf
              call ga_get(g_nexch,ad4+1,ad4+nbf,ad3+1,ad3+nbf,
     &                    scra,nbf)
              ad4=(l-1)*nvir
              call ga_get(g_nt2,ad4+1,ad4+nvir,ad2+1,ad2+nvir,
     &                    scrb,nvir)
              do c=1,nvir
              do d=1,nvir
              do a=1,nvir
           scrc((c-1)*nvir+a)=scrc((c-1)*nvir+a)
     &                  +(2.0d0*scra((nocc+c-1)*nbf+nocc+d)
     &                    -scra((nocc+d-1)*nbf+nocc+c))
     &               *(scrb((d-1)*nvir+a)+t1(i,d)*t1(l,a)
     &               -2.0d0*scrb((a-1)*nvir+d))
              enddo
              enddo
              enddo
              enddo
              call ga_acc(g_qitm,ad2+1,ad2+nvir,ad1+1,ad1+nvir,
     &                    scrc,nvir,1.0d0)
            endif
          enddo
        endif
      enddo
c
      call ga_distribution(g_nt2,iam,g_jlo,g_jhi,g_ilo,g_ihi)
      do i=1,nocc
        ad1=(i-1)*nvir
        if (ad1+1.ge.g_ilo.and.ad1+1.le.g_ihi)then
          do j=1,nocc
            ad2=(j-1)*nvir
            if (ad2+1.ge.g_jlo.and.ad2+1.le.g_jhi)then
              ij=(i-1)*nocc+j
              call ga_get(g_nt2,ad2+1,ad2+nvir,ad1+1,ad1+nvir,
     &                    scra,nvir)
              call dfill(nocc*nocc,0.0d00,scrc,1)
              do k=1,nocc
              do l=1,nocc
                kl=(k-1)*nocc+l
                ad3=(k-1)*nvir
                ad4=(l-1)*nvir
                call ga_get(g_nz2,ad4+1,ad4+nvir,ad3+1,ad3+nvir,
     &                      scrb,nvir)
                do a=1,nvir
                do b=1,nvir
                  scrc((k-1)*nocc+l)=scrc((k-1)*nocc+l)
     &                           +(scra((a-1)*nvir+b)+t1(i,a)*t1(j,b))
     &                           *scrb((a-1)*nvir+b)
                enddo
                enddo
              enddo
              enddo
              call ga_put(g_tz3,1,lnoo,ij,ij,scrc,lnoo)
            endif
          enddo
        endif
      enddo
c
      call ga_sync()
      return
      end
