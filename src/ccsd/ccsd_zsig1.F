      subroutine ccsd_zsig1(nocc,nvir,nbf,hiu,giu,habe,gabe,hia,t1,hz1,
     &                     z1,tzvv,tzoo,scra,scrb,scrc,scrd,scre,scrf,
     &                     g_nt2,g_nz2,g_ncoul,g_nexch,g_x,g_c,
     &                     g_sht2,g_zx,g_zc,g_shz2,g_aitm,g_tz3,iprt)
      implicit none
      integer nocc,nvir,nbf,g_nt2,g_nz2,g_ncoul,
     &        g_nexch,g_x,g_c,g_sht2,g_zx,g_zc,g_shz2,g_tz3,iprt
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
      if (.not.ga_create(MT_DBL,lnoo,lnoo,'aitm',
     &                   nocc,nocc,g_aitm))
     &     call errquit('ga_create g_aitm failed',0)
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
     &                    scrb,nocc)
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
            do c=1,nvir
              hz1(k,c)=hz1(k,c)+
     &             (2.0d0*scra((k-1)*nbf+l)-scrb((k-1)*nbf+l))*z1(l,c)
            enddo
            enddo
            enddo
            do k=1,nocc
            do a=1,nvir
            do c=1,nvir
              hz1(k,c)=hz1(k,c)-(2.0d0*scra((nocc+a-1)*nbf+nocc+c)
     &                   -scrb((nocc+a-1)*nbf+nocc+c))*z1(k,a)
            enddo
            enddo
            enddo
          endif
          do a=1,nvir
          do c=1,nvir
            hz1(j,c)=hz1(j,c)-(2.0d0*scrb((nocc+a-1)*nbf+nocc+c)
     &                   -scra((nocc+a-1)*nbf+nocc+c))*z1(i,a)
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
          do a=1,nvir
          do c=1,nvir
          do k=1,nocc
            hz1(k,c)=hz1(k,c)
     &                 -2.0d0*scrb((nocc+a-1)*nbf+k)*scre((a-1)*nvir+c)
            hz1(k,c)=hz1(k,c)
     &           -2.0d0*(scrf((k-1)*nbf+nocc+a)+scrb((k-1)*nbf+nocc+a))
     &                       *scre((c-1)*nvir+a)
          enddo
          enddo
          enddo
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
      call dfill(nbf*nbf,0.0d00,scrc,1)
      do i=1,nocc
      do l=1,nocc
      do a=1,nvir
        scrc(l,i)=scrc(l,i)+t1(l,a)*z1(i,a)
      enddo
      enddo
      enddo
c
      call dfill(nbf*nbf,0.0d00,scrd,1)
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
        do j=1,nocc
          ij=(i-1)*nocc+j
          if (ij.ge.g_ilo.and.ij.le.g_ihi)then
            ad1=(i-1)*nbf
            ad2=(j-1)*nbf
            call ga_get(g_ncoul,ad2+1,ad2+nbf,ad1+1,ad1+nbf,
     &                  scre,nbf)
            call ga_get(g_nexch,ad2+1,ad2+nbf,ad1+1,ad1+nbf,
     &                  scrf,nbf)
            do c=1,nvir
            do k=1,nocc
              hz1(k,c)=hz1(k,c)+(2.0d0*scre((nocc+c-1)*nbf+k)
     &                              -scrf((nocc+c-1)*nbf+k))*scrc(i,j)
            enddo
            enddo
            do c=1,nvir
            do d=1,nvir
              hz1(i,c)=hz1(i,c)+(2.0d0*scrf((nocc+c-1)*nbf+nocc+d)
     &                   -scrf((nocc+d-1)*nbf+nocc+c))*scrd(j,d)
            enddo
            enddo
          endif
        enddo
      enddo
c
      call ga_distribution(g_aitm,iam,g_jlo,g_jhi,g_ilo,g_ihi)
      do i=1,nocc
      do j=1,nocc
        ij=(i-1)*nocc+j
        if (ij.ge.g_ilo.and.ij.le.g_ihi)then
          call ga_get(g_aitm,1,lnoo,ij,ij,scra,lnoo)
          call dfill(nbf*nbf,0.0d00,scrd,1)
          do k=1,nocc
          do l=1,nocc
          do b=1,nvir
            scrd(k,b)=scrd(k,b)+scra((l-1)*nocc+k)*t1(l,b)
          enddo
          enddo
          enddo
          ad1=(i-1)*nvir
          ad2=(j-1)*nvir
          call ga_get(g_nz2,ad2+1,ad2+nvir,ad1+1,ad1+nvir,scra,nvir)
          do b=1,nvir
          do c=1,nvir
          do k=1,nocc
            hz1(k,c)=hz1(k,c)+4.0d0*scrd(k,b)*scra((b-1)*nvir+c)
          enddo
          enddo
          enddo
        endif
      enddo
      enddo
c
      call ga_distribution(g_shz2,iam,g_jlo,g_jhi,g_ilo,g_ihi)
      do i=1,nocc
        do j=1,nocc
          ij=(i-1)*nocc+j
          if (ij.ge.g_ilo.and.ij.le.g_ihi)then
             call ga_get(g_shz2,1,nbf*nbf,ij,ij,scra,nbf*nbf)
             do c=1,nvir
crk           hz1(j,c)=hz1(j,c)+4.0d0*scra((i-1)*nbf+nocc+c)
             enddo
             do c=1,nvir
             do d=1,nvir
crk           hz1(i,c)=hz1(i,c)+4.0d0*scra((nocc+c-1)*nbf+nocc+d)
crk  &                           *t1(j,d)
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
      do 95 a=1,nvir
      do 95 c=1,nvir
      do 95 d=1,nvir
      do 95 k=1,nocc
      do 95 l=1,nocc
CRK   hz1(k,c)=hz1(k,c)-2.0d0*(2.0d0*moints(k,nocc+a,l,nocc+d)
CRK  &        -moints(k,nocc+d,l,nocc+a))*t1(l,d)*tzvv(a,c)
 95   continue
      do 105 a=1,nvir
      do 105 c=1,nvir
      do 105 d=1,nvir
      do 105 k=1,nocc
      do 105 l=1,nocc
CRK   hz1(k,c)=hz1(k,c)-2.0d0*(2.0d0*moints(k,nocc+c,l,nocc+d)
CRK  &        -moints(k,nocc+d,l,nocc+c))*t1(l,a)*tzvv(d,a)
 105  continue
      do 115 a=1,nvir
      do 115 c=1,nvir
      do 115 d=1,nvir
      do 115 k=1,nocc
CRK   hz1(k,c)=hz1(k,c)+2.0d0*(2.0d0*moints(nocc+a,nocc+d,nocc+c,k)
CRK  &        -moints(nocc+a,nocc+c,nocc+d,k))*tzvv(d,a)
 115  continue
      do 125 c=1,nvir
      do 125 d=1,nvir
      do 125 i=1,nocc
      do 125 k=1,nocc
      do 125 l=1,nocc
CRK   hz1(k,c)=hz1(k,c)-2.0d0*(2.0d0*moints(i,nocc+c,l,nocc+d)
CRK  &        -moints(i,nocc+d,l,nocc+c))*t1(l,d)*tzoo(i,k)
CRK   hz1(k,c)=hz1(k,c)-2.0d0*(2.0d0*moints(k,nocc+c,l,nocc+d)
CRK  &        -moints(k,nocc+d,l,nocc+c))*t1(i,d)*tzoo(l,i)
 125  continue
      do 135 c=1,nvir
      do 135 i=1,nocc
      do 135 k=1,nocc
      do 135 l=1,nocc
CRK   hz1(k,c)=hz1(k,c)-2.0d0*(2.0d0*moints(i,l,nocc+c,k)
CRK  &        -moints(i,k,nocc+c,l))*tzoo(l,i)
 135  continue
c
      call ga_dgop(913,hz1,lnov, '+')
c
      return
      end
