      subroutine ccsd_zsig1(nocc,nvir,nbf,hiu,giu,habe,gabe,hia,t1,hz1,
     &                     z1,tzvv,tzoo,
     &                     scra,scrb,scrc,scrd,scre,scrf,
     &                     g_nt2,g_nz2,g_ncoul,g_nexch,g_x,g_c,
     &                     g_sht2,g_zx,g_zc,g_shz2,g_aitm,g_tz1,g_tz2,
     &                     g_tz3,iprt)
C     $Id: ccsd_zsig1.F,v 2.2 1997-04-14 06:35:05 gg502 Exp $
      implicit none
      integer nocc,nvir,nbf,g_nt2,g_nz2,g_ncoul,g_nexch,g_x,g_c,g_sht2,
     &        g_zx,g_zc,g_shz2,g_tz1,g_tz2,g_tz3,iprt
      double precision hiu(nocc,nocc),giu(nocc,nocc),hia(nocc,nvir)
      double precision habe(nvir,nvir),gabe(nvir,nvir)
      double precision t1(nocc,nvir),hz1(nocc,nvir),z1(nocc,nvir)
      double precision tzvv(nvir,nvir),tzoo(nocc,nocc)
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
      integer g_aitm,
     &        g_jlo,g_jhi,g_ilo,g_ihi,
     &        i,j,k,l,ij,a,b,c,d,ad1,ad2,ad3,ad4,ad5,lnoo,lnov
c
      Integer IAm
      IAM = GA_NodeID()
C
      lnoo=nocc*nocc
      lnov=nocc*nvir
c
      call daxpy(lnov,-2.0d00,hia,1,hz1,1)
c
      call ga_distribution(g_c,iam,g_jlo,g_jhi,g_ilo,g_ihi)
      do i=1,nocc
      do j=1,nocc
        ij=(i-1)*nocc+j
        if (ij.ge.g_ilo.and.ij.le.g_ihi)then
c
          ad1=(i-1)*nocc+j
          call ga_get(g_c,1,nbf*nbf,ad1,ad1,scra,nbf*nbf)
          call ga_get(g_x,1,nbf*nbf,ad1,ad1,scrb,nbf*nbf)
          if (i.eq.j) then
            do k=1,nocc
            do l=1,nocc
              ad1=(k-1)*nbf+l
              scre((l-1)*nocc+k)=2.0d0*scra(ad1)-scrb(ad1)
            enddo
            enddo
            call dgemm('n','n',nocc,nvir,nocc,1.0d00,scre,nocc,z1,
     &                 nocc,1.0d00,hz1,nocc)
            do a=1,nvir
            do c=1,nvir
              ad1=(nocc+a-1)*nbf+nocc+c 
              scre((c-1)*nvir+a)=2.0d0*scra(ad1)-scrb(ad1)
            enddo
            enddo
            call dgemm('n','n',nocc,nvir,nvir,-1.0d00,z1,nocc,scre,
     &                 nvir,1.0d00,hz1,nocc)
            do k=1,nocc
            do a=1,nvir
              ad1=(k-1)*nbf+nocc+a
              scre((a-1)*nocc+k)=2.0d0*scra(ad1)-scrb(ad1)
            enddo
            enddo
            call dgemm('n','n',nocc,nvir,nvir,-2.0d00,scre,nocc,tzvv,
     &                 nvir,1.0d00,hz1,nocc)
            call dgemm('t','n',nocc,nvir,nocc,-2.0d00,tzoo,nocc,scre,
     &                 nocc,1.0d00,hz1,nocc)
          endif
          do a=1,nvir
          do c=1,nvir
            hz1(j,c)=hz1(j,c)-(2.0d0*scrb((nocc+a-1)*nbf+nocc+c)
     &                   -scra((nocc+a-1)*nbf+nocc+c))*z1(i,a)
          enddo
          enddo
          do k=1,nocc
          do c=1,nvir
            hz1(k,c)=hz1(k,c)-2.0d0*(2.0d0*scra((k-1)*nbf+nocc+c)
     &                 -scrb((k-1)*nbf+nocc+c))*tzoo(j,i)
          enddo
          enddo
          ad1=(i-1)*nvir
          ad2=(j-1)*nvir
          call ga_get(g_nz2,ad2+1,ad2+nvir,ad1+1,ad1+nvir,
     &                scre,nvir)
          ad3=(i-1)*nbf
          ad4=(j-1)*nbf
          call ga_get(g_nexch,ad4+1,ad4+nbf,ad3+1,ad3+nbf,
     &                scrf,nbf)
          call dgemm('n','t',nocc,nvir,nvir,-2.0d00,scrd(1,nocc+1),nbf,
     &               scre,nvir,1.0d00,hz1,nocc)
          do a=1,nvir
          do k=1,nocc
            ad1=(k-1)*nbf+nocc+a
            scra((a-1)*nocc+k)=scrf(ad1)+scrb(ad1)
          enddo
          enddo
          call dgemm('n','n',nocc,nvir,nvir,-2.0d00,scra,nocc,scre,
     &               nvir,1.0d00,hz1,nocc)
        endif
      enddo
      enddo
c
      call dgemm('n','n',nocc,nvir,nocc,1.0d00,hiu,nocc,z1,nocc,1.0d00,
     &           hz1,nocc)
      call dgemm('n','t',nocc,nvir,nvir,-1.0d00,z1,nocc,habe,nvir,
     &            1.0d00,hz1,nocc)
c
      call ga_distribution(g_zc,iam,g_jlo,g_jhi,g_ilo,g_ihi)
      do i=1,nocc
      do j=1,nocc
        ij=(i-1)*nocc+j
        if (ij.ge.g_ilo.and.ij.le.g_ihi)then
c
          ad1=(i-1)*nocc+j
          call ga_get(g_zc,1,nbf*nbf,ad1,ad1,scra,nbf*nbf)
          call ga_get(g_zx,1,nbf*nbf,ad1,ad1,scrb,nbf*nbf)
          if(i.eq.j)then
            do c=1,nvir
            do k=1,nocc
               hz1(k,c)=hz1(k,c)-(2.0d0*scra((nocc+c-1)*nbf+k)
     &                   -scrb((nocc+c-1)*nbf+k))
            enddo
            enddo
          endif
        endif
      enddo
      enddo
c
      call dfill(lnoo,0.0d0,scrc,1)
      do i=1,nocc
      do l=1,nocc
        scrc(l,i)=-2.0d0*tzoo(l,i)
      do a=1,nvir
        scrc(l,i)=scrc(l,i)+t1(l,a)*z1(i,a)
      enddo
      enddo
      enddo
c
      call dgemm('n','t',nocc,nvir,nvir,-2.0d00,t1,nocc,tzvv,nvir,
     &            0.0d00,scrd,nbf)
c
      call ga_distribution(g_nt2,iam,g_jlo,g_jhi,g_ilo,g_ihi)
      do i=1,nocc
        ad1=(i-1)*nvir
        if (ad1+1.ge.g_ilo.and.ad1+1.le.g_ihi)then
          do j=1,nocc
            ad2=(j-1)*nvir
            if (ad2+1.ge.g_jlo.and.ad2+1.le.g_jhi)then
              call ga_get(g_nt2,ad2+1,ad2+nvir,ad1+1,ad1+nvir,
     &                    scre,nvir)
              do a=1,nvir
                do b=1,nvir
                ad3=(a-1)*nvir+b
                ad4=(b-1)*nvir+a
                scrd(i,b)=scrd(i,b)+(scre(ad3)+t1(i,a)*t1(j,b)
     &                              -2.0d0*scre(ad4))*z1(j,a)
                enddo
              enddo
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
              call ga_get(g_ncoul,ad2+1,ad2+nbf,ad1+1,ad1+nbf,
     &                    scre,nbf)
              call ga_get(g_nexch,ad2+1,ad2+nbf,ad1+1,ad1+nbf,
     &                    scrf,nbf)
              do c=1,nvir
              do k=1,nocc
                hz1(k,c)=hz1(k,c)+(2.0d0*scre((nocc+c-1)*nbf+k)
     &                             -scrf((nocc+c-1)*nbf+k))*scrc(i,j)
              enddo
              enddo
              do c=1,nvir
              do d=1,nvir
                hz1(i,c)=hz1(i,c)+(2.0d0*scrf((nocc+c-1)*nbf+nocc+d)
     &                         -scrf((nocc+d-1)*nbf+nocc+c))*scrd(j,d)
              enddo
              enddo
            endif
          enddo
        endif
      enddo
c
      call ga_distribution(g_aitm,iam,g_jlo,g_jhi,g_ilo,g_ihi)
      do i=1,nocc
      do j=1,nocc
        ij=(i-1)*nocc+j
        if (ij.ge.g_ilo.and.ij.le.g_ihi)then
          call ga_get(g_aitm,1,lnoo,ij,ij,scra,lnoo)
          call dgemm('n','n',nocc,nvir,nocc,1.0d00,scra,nocc,t1,nocc,
     &               0.0d00,scrd,nocc)
          ad1=(i-1)*nvir
          ad2=(j-1)*nvir
          call ga_get(g_nz2,ad2+1,ad2+nvir,ad1+1,ad1+nvir,scra,nvir)
          call dgemm('n','t',nocc,nvir,nvir,4.0d00,scrd,nocc,scra,nvir,
     &               1.0d00,hz1,nocc)
        endif
      enddo
      enddo
c
      call ga_distribution(g_tz3,iam,g_jlo,g_jhi,g_ilo,g_ihi)
      do i=1,nocc
      do j=1,nocc
        ij=(i-1)*nocc+j
        if (ij.ge.g_ilo.and.ij.le.g_ihi)then
          call ga_get(g_tz3,1,lnoo,ij,ij,scra,lnoo)
          ad1=(i-1)*nbf
          ad2=(j-1)*nbf
          call ga_get(g_nexch,ad2+1,ad2+nbf,ad1+1,ad1+nbf,scrb,nbf)
          do c=1,nvir
          do k=1,nocc
          do l=1,nocc
            hz1(l,c)=hz1(l,c)+2.0d0*scrb((nocc+c-1)*nbf+k)
     &                         *scra((l-1)*nocc+k)
          enddo
          enddo
          enddo
        endif
      enddo
      enddo
c
      call ga_distribution(g_tz3,iam,g_ilo,g_ihi,g_jlo,g_jhi)
      do i=1,nocc
        ad1=(i-1)*nocc
        if (ad1+1.ge.g_ilo.and.ad1+1.le.g_ihi)then
          do j=1,nocc
            ad2=(j-1)*nocc
            if (ad2+1.ge.g_jlo.and.ad2+1.le.g_jhi)then
              call ga_get(g_tz3,ad1+1,ad1+nocc,ad2+1,ad2+nocc,
     &                    scra,nocc)
              ij=(i-1)*nocc+j
              call ga_get(g_c,1,nbf*nbf,ij,ij,scrb,nbf*nbf)
              do k=1,nocc
              do l=1,nocc
              do c=1,nvir
                hz1(k,c)=hz1(k,c)+2.0d0*scrb((nocc+c-1)*nbf+l)
     &                              *scra((l-1)*nocc+k)
              enddo
              enddo
              enddo
            endif
          enddo
        endif
      enddo
c
      call ga_distribution(g_shz2,iam,g_jlo,g_jhi,g_ilo,g_ihi)
      do i=1,nocc
        do j=1,nocc
          ij=(i-1)*nocc+j
          if (ij.ge.g_ilo.and.ij.le.g_ihi)then
             call ga_get(g_shz2,1,nbf*nbf,ij,ij,scra,nbf*nbf)
             do c=1,nvir
               hz1(j,c)=hz1(j,c)+4.0d0*scra((i-1)*nbf+nocc+c)
             enddo
             do c=1,nvir
             do d=1,nvir
              hz1(i,c)=hz1(i,c)+4.0d0*scra((nocc+c-1)*nbf+nocc+d)
     &                           *t1(j,d)
             enddo
             enddo
          endif
        enddo
      enddo
c
      call ga_distribution(g_sht2,iam,g_jlo,g_jhi,g_ilo,g_ihi)
      do i=1,nocc
        do j=1,nocc
          ij=(i-1)*nocc+j
          if (ij.ge.g_ilo.and.ij.le.g_ihi)then
             call ga_get(g_sht2,1,nbf*nbf,ij,ij,scra,nbf*nbf)
             ad1=(i-1)*nvir
             ad2=(j-1)*nvir
             call ga_get(g_nz2,ad2+1,ad2+nvir,ad1+1,ad1+nvir,
     &                   scrb,nvir)
             do k=1,nocc
             do a=1,nvir
             do c=1,nvir
              hz1(k,c)=hz1(k,c)-4.0d0*scra((k-1)*nbf+nocc+a)
     &                               *scrb((c-1)*nvir+a)
             enddo
             enddo
             enddo
          endif
        enddo
      enddo
c
 105  continue
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
              call dfill(nbf*nbf,0.0d0,scrd,1)
              do a=1,nvir
              do b=1,nvir
              do l=1,nocc
                scrd(l,a)=scrd(l,a)+t1(l,b)*scra((a-1)*nvir+b)
              enddo
              enddo
              enddo
              do l=1,nocc
                ad3=(i-1)*nbf
                ad4=(l-1)*nbf
                ad5=(i-1)*nocc+l
                call ga_get(g_x,1,nbf*nbf,ad5,ad5,scre,nbf*nbf)
                call ga_get(g_nexch,ad4+1,ad4+nbf,ad3+1,ad3+nbf,scrf,
     &                     nbf)
                do a=1,nvir
                do c=1,nvir
                  ad3=(nocc+a-1)*nbf+nocc+c
                  scre(ad3)=scre(ad3)+scrf(ad3)
                enddo
                enddo
                do a=1,nvir
                do c=1,nvir
                  hz1(j,c)=hz1(j,c)-2.0d0*scre((nocc+a-1)*nbf+nocc+c)
     &                              *scrd(l,a)
                enddo
                enddo
              enddo
              call dfill(nbf*nbf,0.0d0,scrd,1)
              do a=1,nvir
              do b=1,nvir
              do l=1,nocc
                scrd(l,a)=scrd(l,a)+t1(l,b)*scra((b-1)*nvir+a)
              enddo
              enddo
              enddo
              do l=1,nocc
                ad3=(i-1)*nbf
                ad4=(l-1)*nbf
                ad5=(i-1)*nocc+l
                call ga_get(g_c,1,nbf*nbf,ad5,ad5,scre,nbf*nbf)
                call ga_get(g_ncoul,ad4+1,ad4+nbf,ad3+1,ad3+nbf,scrf,
     &                      nbf)
                do a=1,nvir
                do c=1,nvir
                  ad3=(nocc+a-1)*nbf+nocc+c
                  scre(ad3)=scre(ad3)+scrf(ad3)
                enddo
                enddo
                do a=1,nvir
                do c=1,nvir
                  hz1(j,c)=hz1(j,c)-2.0d0*scre((nocc+a-1)*nbf+nocc+c)
     &                              *scrd(l,a)
                enddo
                enddo
              enddo
            endif
          enddo
        endif
      enddo
c
      call ga_distribution(g_tz1,iam,g_jlo,g_jhi,g_ilo,g_ihi)
      do i=1,nocc
        ad1=(i-1)*nvir
        if (ad1+1.ge.g_ilo.and.ad1+1.le.g_ihi)then
          do j=1,nocc
            ad2=(j-1)*nvir
            if (ad2+1.ge.g_jlo.and.ad2+1.le.g_jhi)then
              call ga_get(g_tz1,ad2+1,ad2+nvir,ad1+1,ad1+nvir,
     &                    scra,nvir)
              ad3=(i-1)*nbf
              ad4=(j-1)*nbf
              call ga_get(g_nexch,ad4+1,ad4+nbf,ad3+1,ad3+nbf,scrb,nbf)
              ij=(i-1)*nocc+j
              call ga_get(g_x,1,nbf*nbf,ij,ij,scre,nbf*nbf)
              do a=1,nvir
              do k=1,nocc
                ad5=(k-1)*nbf+nocc+a
                scrb(ad5)=scrb(ad5)+scre(ad5)
              enddo
              enddo
              do a=1,nvir
              do c=1,nvir
              do k=1,nocc
                 hz1(k,c)=hz1(k,c)-2.0d0*scrb((k-1)*nbf+nocc+a)
     &                    *scra((c-1)*nvir+a)
              enddo
              enddo
              enddo
            endif
          enddo
        endif
      enddo
c
      call ga_distribution(g_tz2,iam,g_jlo,g_jhi,g_ilo,g_ihi)
      do i=1,nocc
        ad1=(i-1)*nvir
        if (ad1+1.ge.g_ilo.and.ad1+1.le.g_ihi)then
          do j=1,nocc
            ad2=(j-1)*nvir
            if (ad2+1.ge.g_jlo.and.ad2+1.le.g_jhi)then
              call ga_get(g_tz2,ad2+1,ad2+nvir,ad1+1,ad1+nvir,
     &                    scra,nvir)
              ad3=(i-1)*nbf
              ad4=(j-1)*nbf
              call ga_get(g_ncoul,ad4+1,ad4+nbf,ad3+1,ad3+nbf,scrb,nbf)
              ij=(i-1)*nocc+j
              call ga_get(g_c,1,nbf*nbf,ij,ij,scre,nbf*nbf)
              do a=1,nvir
              do k=1,nocc
                ad5=(k-1)*nbf+nocc+a
                scrb(ad5)=scrb(ad5)+scre(ad5)
              enddo
              enddo
              do a=1,nvir
              do c=1,nvir
              do k=1,nocc
                 hz1(k,c)=hz1(k,c)+2.0d0*scrb((k-1)*nbf+nocc+a)
     &                    *scra((c-1)*nvir+a)
              enddo
              enddo
              enddo
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
     &                    scrc,nbf)
              do k=1,nocc
                ad3=(i-1)*nvir
                ad4=(k-1)*nvir
                call ga_get(g_tz1,ad3+1,ad3+nvir,ad4+1,ad4+nvir,
     &                    scre,nvir)
                call ga_get(g_tz2,ad3+1,ad3+nvir,ad4+1,ad4+nvir,
     &                    scrf,nvir)
                call dgemm('n','n',nvir,nvir,nvir,-2.0d00,
     &                     scrc(nocc+1,nocc+1),nbf,scre,nvir,0.0d00,
     &                     scrd,nbf)
                call dgemm('t','n',nvir,nvir,nvir,2.0d00,
     &                     scrc(nocc+1,nocc+1),nbf,scrf,nvir,1.0d00,
     &                     scrd,nbf)
                do a=1,nvir
                do c=1,nvir
                  hz1(k,c)=hz1(k,c)+scrd(c,a)*t1(j,a)
                enddo
                enddo
              enddo
            endif
          enddo
        endif
      enddo
c
      call ga_dgop(913,hz1,lnov, '+')
c
      return
      end
