      subroutine ccsd_tzintb(nocc,nvir,nbf,t1,tzoo,scra,scrb,scrc,
     &                       g_nt2,g_nz2,g_nexch,g_x,g_sht2,g_aitm,
     &                       g_tz3,iprt)
C     $Id: ccsd_tzintb.F,v 2.2 1997-04-14 06:35:03 gg502 Exp $
      implicit none
      integer nocc,nvir,nbf,lnoo,lnov,g_nt2,g_nz2,g_nexch,g_x,g_sht2,
     &        g_aitm,g_tz3,g_tmp,iprt
      integer i,j,k,l,a,b,ij,kl,ad1,ad2,ad3,ad4,g_jlo,g_jhi,g_ilo,g_ihi
      double precision t1(nocc,nvir),tzoo(nocc,nocc),scra(nbf*nbf),
     &                 scrb(nbf*nbf),scrc(nocc*nocc),x,ddot
      Integer IAm
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
c
      IAM = GA_NodeID()
C
      lnoo=nocc*nocc
      lnov=nocc*nvir
c
      call dfill(lnoo,0.0d00,tzoo,1)
      if (.not.ga_create(MT_DBL,lnov,lnov,'tmp',
     &                   nvir,nvir,g_tmp))
     &     call errquit('ga_create g_tmp failed',0)
      call ga_dgemm('n','t',lnov,lnov,lnov,
     &              1.0d00,g_nt2,g_nz2,0.0d00,g_tmp)
c
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
      if (.not.ga_destroy(g_tmp))
     &    call errquit('ga_dest g_tmp fail',0)
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
      return
      end
