      subroutine ccsd_itm(ncor,ndel,nocc,nvir,nbf,hiu,giu,habe,gabe,
     &                    hia,t1,ht1,
     &                    scra,scrb,scrc,scrd,scre,scrf,
     &                    g_nt2,g_nht2,g_ncoul,g_nexch,g_x,g_c,g_sht2,
     &                    iprt)
      implicit none
      integer ncor,ndel,nocc,nvir,nbf,g_nt2,g_nht2,g_ncoul,g_nexch,
     &        g_x,g_c,g_sht2,iprt
      double precision hiu(nocc,nocc),giu(nocc,nocc),hia(nocc,nvir)
      double precision habe(nvir,nvir),gabe(nvir,nvir)
      double precision t1(nocc,nvir),ht1(nocc,nvir)
      double precision scra(nbf*nbf),scrb(nbf*nbf),
     &                 scrc(nbf,nbf),scrd(nbf,nbf),
     &                 scre(nbf*nbf),scrf(nbf*nbf)
c
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "bas.fh"
#include "rtdb.fh"
c
      integer doa,dob,dog,doh,dojk,dos,dod
      common/debug/doa(15),dob(15),dog(15),doh(15),dojk(15),
     &             dos(15),dod(15)
c
      integer g_aitm,g_tmt1,g_tmt2,g_tmt,g_tmp,g_tmp2,
     &        g_jlo,g_jhi,g_ilo,g_ihi,
     &        i,j,k,l,ij,a,b,c,ad1,ad2,ad3,ad4,ad5,ad6,lnoo,lnov,lnvv
      integer inode,next,aold
c
      integer nodes,iam
      common/para1/nodes,iam
c
      call qenter('itm',0)
c
      lnoo=nocc*nocc
      lnov=nocc*nvir
      lnvv=nvir*nvir
c
c - test x and c by doing its contribution to (2j-k)*(2t-t)
      if (.not.ga_create(MT_DBL,lnoo,lnoo,'aitm',
     &                   nocc,nocc,g_aitm))
     &     call errquit('ga_create g_aitm failed',0)
      if (.not.ga_create(MT_DBL,lnvv,lnoo,'tmt1',
     &                   lnvv,1,g_tmt1))
     &     call errquit('ga_create g_tmt1 failed',0)
      if (.not.ga_create(MT_DBL,lnvv,lnoo,'tmt2',
     &                   lnvv,1,g_tmt2))
     &     call errquit('ga_create g_tmt2 failed',0)
c
      call dfill(lnvv,0.0d00,habe,1)
      call ga_distribution(g_nexch,iam,g_jlo,g_jhi,g_ilo,g_ihi)
      do i=1,nocc
        ad1=(i-1)*nbf
        if (ad1+1.ge.g_ilo.and.ad1+1.le.g_ihi)then
          do j=1,nocc
            ad2=(j-1)*nbf
            if (ad2+1.ge.g_jlo.and.ad2+1.le.g_jhi)then
c-----------------------------------------------------------------------
c  S7
c-----------------------------------------------------------------------
            call ga_get(g_nexch,ad2+1,ad2+nbf,ad1+1,ad1+nbf,
     &                  scrc,nbf)
            ad3=(i-1)*nvir
            ad4=(j-1)*nvir
            call ga_get(g_nt2,ad4+1,ad4+nvir,ad3+1,ad3+nvir,
     &                    scre,nvir)
            do k=1,nocc
            do a=1,nvir
            do b=1,nvir
            ht1(k,a)=ht1(k,a)-(2.0d00*scrc(nocc+b,k)-scrc(k,nocc+b))*
     &         (scre((a-1)*nvir+b)+t1(i,a)*t1(j,b))
            enddo
            enddo
            enddo
c
            do c=1,nvir
            do a=1,nvir
            do b=1,nvir
            habe(a,b)=habe(a,b)-(2.0d00*scrc(nocc+b,nocc+c)
     &         -scrc(nocc+c,nocc+b))*
     &         (scre((c-1)*nvir+a)+t1(i,c)*t1(j,a))
            enddo
            enddo
            enddo
            endif
          enddo
        endif
      enddo
c
      call dfill(lnoo,0.0d00,hiu,1)
      call ga_distribution(g_sht2,iam,g_jlo,g_jhi,g_ilo,g_ihi)
      do i=1,nocc
        do j=1,nocc
          ij=(i-1)*nocc+j
          if (ij.ge.g_ilo.and.ij.le.g_ihi)then
            call ga_get(g_sht2,1,nbf*nbf,ij,ij,scre,nbf*nbf)
            call sscal(nbf*nbf,2.0d00,scre,1)
c
c H2
            do k=1,nocc
                hiu(k,i)=hiu(k,i)+scre((k-1)*nbf+j)+scre((k-1)*nbf+j)
     &                          -scre((j-1)*nbf+k)
            enddo
          endif
        enddo
      enddo
c
      call scopy(lnoo,hiu,1,giu,1)
      call scopy(lnvv,habe,1,gabe,1)
      call dfill(lnov,0.0d00,hia,1)
c
      if (iprt.gt.5.and.iam.eq.0)print *,'forming aitm'
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
            if (doa(1).eq.2)then
              call ga_get(g_nexch,ad2+1,ad2+nocc,ad1+1,ad1+nocc,
     &                    scrb,nocc)
              call saxpy(nocc*nocc,0.5d00,scrb,1,scra,1)
            endif
c---------------------------------------------------------
c  S5, D1, D6 and D7 
c---------------------------------------------------------
            call ga_get(g_nexch,ad2+nocc+1,ad2+nbf,ad1+1,ad1+nbf,
     &                  scre,nvir)
            call dfill(nbf*nbf,0.0d00,scrf,1)
c
c D1,D7
            ad1=0
            do a=1,nvir
              do b=nocc+1,nbf
                ad1=ad1+1
                scrf(ad1)=scrf(ad1)+0.5d00*scre((b-1)*nvir+a)
                do k=1,nocc
                  scrf(ad1)=scrf(ad1)-scre((k-1)*nvir+a)*t1(k,b-nocc)
                enddo
              enddo
            enddo
c
            ad1=(i-1)*nocc+j
            call ga_get(g_c,1,nbf*nbf,ad1,ad1,scrb,nbf*nbf)
            call scopy(nbf*nbf,scrb,1,scre,1)
c
c H8
            do a=nocc+1,nbf
              hia(j,a-nocc)=hia(j,a-nocc)-scre((a-1)*nbf+i)
            enddo
c
c G4
            if (i.eq.j) then
            do b=nocc+1,nbf
            do a=nocc+1,nbf
              gabe(a-nocc,b-nocc)=gabe(a-nocc,b-nocc)+scre((a-1)*nbf+b)
     &                            +scre((a-1)*nbf+b)
            enddo
            enddo
            endif
c
            ad1=(i-1)*nocc+j
            call ga_get(g_x,1,nbf*nbf,ad1,ad1,scrb,nbf*nbf)
            call scopy(nbf*nbf,scrb,1,scre,1)
c
c H8
            do a=nocc+1,nbf
              hia(j,a-nocc)=hia(j,a-nocc)+scre((i-1)*nbf+a)
     &                      +scre((i-1)*nbf+a)
            enddo
c
c S5
            do a=nocc+1,nbf
              ht1(j,a-nocc)=ht1(j,a-nocc)+scre((i-1)*nbf+a)
     &                      +scre((i-1)*nbf+a)-scre((a-1)*nbf+i)
            enddo
c
c G2
            do k=1,nocc
              giu(k,j)=giu(k,j)+scre((i-1)*nbf+k)
     &                      +scre((i-1)*nbf+k)-scre((k-1)*nbf+i)
            enddo
c
c G4
            if (i.eq.j) then
            do b=nocc+1,nbf
            do a=nocc+1,nbf
              gabe(a-nocc,b-nocc)=gabe(a-nocc,b-nocc)-scre((a-1)*nbf+b)
            enddo
            enddo
            endif
c
c D6,D7
            ad1=0
            do a=nocc+1,nbf
              do b=nocc+1,nbf
                ad1=ad1+1
                scrf(ad1)=scrf(ad1)+scre((b-1)*nbf+a)
                do k=1,nocc
                  scrf(ad1)=scrf(ad1)-scre((b-1)*nbf+k)*t1(k,a-nocc)
                  scrf(ad1)=scrf(ad1)-scre((k-1)*nbf+a)*t1(k,b-nocc)
                enddo
              enddo
            enddo
            ad1=(i-1)*nvir
            ad2=(j-1)*nvir
            call ga_acc(g_nht2,ad1+1,ad1+nvir,ad2+1,ad2+nvir,
     &                  scrf,nvir,1.0d00)
c
            if (doa(2).eq.2.and.doa(3).eq.2)then
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
            endif
c
            call ga_get(g_sht2,1,nbf*nbf,ij,ij,scrb,nbf*nbf)
            call dfill(nbf*nbf,0.0d00,scre,1)
            call saxpy(nbf*nbf,2.0d00,scrb,1,scre,1)
c
c S6
            do a=nocc+1,nbf
              ht1(j,a-nocc)=ht1(j,a-nocc)+scre((i-1)*nbf+a)
     &                      +scre((i-1)*nbf+a)-scre((a-1)*nbf+i)
            enddo
c
c A4
            if (doa(4).eq.2)then
            ad1=0
            do k=1,nocc
              do l=1,nocc
                ad1=ad1+1
                scra(ad1)=scra(ad1)+scrb((k-1)*nbf+l)
              enddo
            enddo
            endif
            call ga_put(g_aitm,1,lnoo,ij,ij,scra,lnoo)
c-----------------------------------------------------------------------
c  B2 and B3 intermediate and their contribution to D3
c-----------------------------------------------------------------------
            if (dob(2).eq.2)then
            ad1=nocc+1
            call sgemm('n','n',nvir,nvir,nocc,1.0d00,scrb(ad1),
     &                 nbf,t1,nocc,0.0d00,scra,nvir)
            ad1=nocc*nbf+1
            call sgemm('t','n',nvir,nvir,nocc,1.0d00,t1,nocc,
     &                 scrb(ad1),nbf,1.0d00,scra,nvir)
            ad1=(i-1)*nvir
            ad2=(j-1)*nvir
            call ga_acc(g_nht2,ad2+1,ad2+nvir,ad1+1,ad1+nvir,
     &                  scra,nvir,-1.0d00)
            endif
c
          endif
        enddo
      enddo
c
      call ga_dgop(900,hiu,lnoo, '+')
      call ga_dgop(901,habe,lnvv, '+')
      call ga_dgop(902,giu,lnoo, '+')
      call ga_dgop(903,hia,lnov, '+')
      call ga_dgop(904,gabe,lnvv, '+')
c
c-----------------------------------------------------------------------
c  D2 (product of A intermediate and Tau)
c-----------------------------------------------------------------------
c construct tau
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
              ad3=0
              do a=1,nvir
                do b=1,nvir
                  ad3=ad3+1
                  scra(ad3)=scra(ad3)+t1(i,a)*t1(j,b)
                enddo
              enddo
              call ga_put(g_tmt1,1,lnvv,ij,ij,scra,lnvv)
            endif
          enddo
        endif
      enddo
c
      call ga_dgemm('n','n',lnvv,lnoo,lnoo,
     &              1.0d00,g_tmt1,g_aitm,0.0d00,g_tmt2)
c
      call ga_distribution(g_nht2,iam,g_jlo,g_jhi,g_ilo,g_ihi)
      do i=1,nocc
        ad1=(i-1)*nvir
        if (ad1+1.ge.g_ilo.and.ad1+1.le.g_ihi)then
          do j=1,nocc
            ad2=(j-1)*nvir
            if (ad2+1.ge.g_jlo.and.ad2+1.le.g_jhi)then
              ij=(i-1)*nocc+j
              call ga_get(g_tmt2,1,lnvv,ij,ij,scra,lnvv)
              call ga_acc(g_nht2,ad2+1,ad2+nvir,ad1+1,ad1+nvir,
     &                    scra,nvir,1.0d00)
            endif
          enddo
        endif
      enddo
c
c-----------------------------------------------------------------------
c  S4
c-----------------------------------------------------------------------
c
c construct (2*t2_ij^ab - t2_ij^ba + t_i^a t_j^b)
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
                ht1(i,a)=ht1(i,a)+(scra(ad4)+scra(ad4)-scra(ad3)
     &                    +t1(i,b)*t1(j,a))*hia(j,b)
                enddo
              enddo
            endif
          enddo
        endif
      enddo
c
c-----------------------------------------------------------------------
c  D4
c-----------------------------------------------------------------------
c
      call ga_distribution(g_nht2,iam,g_jlo,g_jhi,g_ilo,g_ihi)
      do i=1,nocc
        ad1=(i-1)*nvir
        if (ad1+1.ge.g_ilo.and.ad1+1.le.g_ihi)then
          do j=1,nocc
            ad2=(j-1)*nvir
            if (ad2+1.ge.g_jlo.and.ad2+1.le.g_jhi)then
              call ga_get(g_nt2,ad2+1,ad2+nvir,ad1+1,ad1+nvir,
     &                    scre,nvir)
              call sgemm('n','n',nvir,nvir,nvir,1.0d00,gabe,
     &                   nvir,scre,nvir,0.0d00,scrd,nvir)
              call ga_acc(g_nht2,ad2+1,ad2+nvir,ad1+1,ad1+nvir,
     &                    scrd,nvir,1.0d00)
            endif
          enddo
        endif
      enddo
c
c-----------------------------------------------------------------------
c  D5
c-----------------------------------------------------------------------
c
c
      inode=-1
      next=nxtval(-nodes)
      next=nxtval(nodes)
      do a=1,nvir
          do i=1,nocc
            ad1=(i-1)*nvir+a
        inode=inode+1
        if (inode.eq.next)then
            call ga_get(g_nt2,1,lnov,ad1,ad1,scre,1)
            call sgemm('n','n',nvir,nocc,nocc,-1.0d00,scre,
     &                 nvir,giu,nocc,0.0d00,scrf,nvir)
            call ga_acc(g_nht2,1,lnov,ad1,ad1,scrf,lnov,1.0d00)
      next=nxtval(nodes)
          endif
        enddo
      enddo
c
      if (.not.ga_destroy(g_tmt2))
     &    call errquit('ga_dest g_tmt2 fail',0)
      if (.not.ga_destroy(g_tmt1))
     &    call errquit('ga_dest g_tmt1 fail',0)
      if (.not.ga_destroy(g_aitm))
     &    call errquit('ga_dest g_aitm fail',0)
c
c-----------------------------------------------------------------------
c  J1, J2, J3 and K1, K2, K3 intermediates
c-----------------------------------------------------------------------
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
            call sgemm('n','n',nvir,nvir,nocc,-1.0d00,scrd(nocc+1,1),
     &                 nbf,t1,nocc,0.0d00,scrc(nocc+1,nocc+1),nbf)
            call saxpy(nbf*nbf,1.0d00,scrc,1,scrd,1)
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
            call sgemm('n','n',nvir,nvir,nocc,-1.0d00,scrd(nocc+1,1),
     &                 nbf,t1,nocc,0.0d00,scrc(nocc+1,nocc+1),nbf)
            call saxpy(nbf*nbf,1.0d00,scrc,1,scrd,1)
            call ga_acc(g_c,1,nbf*nbf,ij,ij,scrd,nbf*nbf,1.0d00)
c
          endif
        enddo
      enddo
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
c-----------------------------------------------------------------------
c  D9 and D10 (product of K intermediate and T2)
c-----------------------------------------------------------------------
      call ga_distribution(g_tmp2,iam,g_jlo,g_jhi,g_ilo,g_ihi)
      do i=1,nocc
        ad1=(i-1)*nvir
        if (ad1+1.ge.g_ilo.and.ad1+1.le.g_ihi)then
          do j=1,nocc
            ad2=(j-1)*nvir
            if (ad2+1.ge.g_jlo.and.ad2+1.le.g_jhi)then
              ij=(i-1)*nocc+j
              call ga_get(g_c,1,nbf*nbf,ij,ij,scrb,nbf*nbf)
              ad3=0
              do a=nocc+1,nbf
                do b=nocc+1,nbf
                  ad3=ad3+1
                  scra(ad3)=scrb((a-1)*nbf+b)
                enddo
              enddo
              call ga_acc(g_tmp2,ad2+1,ad2+nvir,ad1+1,ad1+nvir,
     &                    scra,nvir,1.0d00)
            endif
          enddo
        endif
      enddo
c
c construct (t2_ij^ba)
c
      call ga_distribution(g_nt2,iam,g_jlo,g_jhi,g_ilo,g_ihi)
      do i=1,nocc
        ad1=(i-1)*nvir
        if (ad1+1.ge.g_ilo.and.ad1+1.le.g_ihi)then
          do j=1,nocc
            ad2=(j-1)*nvir
            if (ad2+1.ge.g_jlo.and.ad2+1.le.g_jhi)then
              call ga_get(g_nt2,ad1+1,ad1+nvir,ad2+1,ad2+nvir,
     &                    scra,nvir)
              call ga_put(g_tmt,ad2+1,ad2+nvir,ad1+1,ad1+nvir,
     &                     scra,nvir)
            endif
          enddo
        endif
      enddo
c
      call ga_dgemm('t','n',lnov,lnov,lnov,
     &              1.0d00,g_tmp2,g_tmt,0.0d00,g_tmp)
c
c     d9 and d10
c
      call ga_distribution(g_nht2,iam,g_jlo,g_jhi,g_ilo,g_ihi)
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
                scrb(ad3)=-0.5d0*scra(ad3)-scra(ad4)
                enddo
              enddo
              call ga_acc(g_nht2,ad2+1,ad2+nvir,ad1+1,ad1+nvir,
     &                    scrb,nvir,1.0d00)
            endif
          enddo
        endif
      enddo
c
      call ga_sync()
c
      if (iprt.gt.5.and.iam.eq.0)print *,'forming j - k'
      call ga_distribution(g_x,iam,g_jlo,g_jhi,g_ilo,g_ihi)
      do i=1,nocc
        do j=1,nocc
          ij=(i-1)*nocc+j
          if (ij.ge.g_ilo.and.ij.le.g_ihi)then
            ad1=(i-1)*nbf
            ad2=(j-1)*nbf
            call dfill(nbf*nbf,0.0d00,scra,1)
              call ga_get(g_x,1,nbf*nbf,ij,ij,scrb,nbf*nbf)
              call saxpy(nbf*nbf,1.0d00,scrb,1,scra,1)
              call ga_get(g_c,1,nbf*nbf,ij,ij,scrb,nbf*nbf)
              call saxpy(nbf*nbf,-0.5d00,scrb,1,scra,1)
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
              call ga_put(g_tmt,ad2+1,ad2+nvir,ad1+1,ad1+nvir,
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
      call ga_dgemm('n','t',lnov,lnov,lnov,
     &              0.25d00,g_tmp2,g_tmt,1.0d00,g_tmp)
c
c construct (2*t2_ij^ab - t2_ij^ba)
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
                ad4=(b-1)*nvir+a
                scrb(ad3)=scra(ad3)+scra(ad3)-scra(ad4)
                enddo
              enddo
              call ga_put(g_tmt,ad2+1,ad2+nvir,ad1+1,ad1+nvir,
     &                     scrb,nvir)
            endif
          enddo
        endif
      enddo
c
      call ga_dgemm('t','n',lnov,lnov,lnov,
     &              1.0d00,g_tmp,g_tmt,1.0d00,g_nht2)
c
      if (.not.ga_destroy(g_tmp2))
     &    call errquit('ga_dest g_tmp2 fail',0)
      if (.not.ga_destroy(g_tmp))
     &    call errquit('ga_dest g_tmp fail',0)
      if (.not.ga_destroy(g_tmt))
     &    call errquit('ga_dest g_tmt fail',0)
c
      aold=0
      inode=-1
      next=nxtval(-nodes)
      next=nxtval(nodes)
      do a=1,nvir
       do b=1,nvir
c
c -- parallel task --
        inode=inode+1
        if (inode.eq.next)then
c
         call ccsd_fsig1(a,b,ncor,nocc,nvir,ndel,nbf,iprt,
     $                   t1,ht1,habe)
         next=nxtval(nodes)
        endif
c -- end parallel task --
c
       enddo
      enddo
c
c -- add some trivial terms to ht1 (nothing to ht2)
      inode=-1
      do a=1,nvir
       do i=1,nocc
c
c -- parallel task --
        inode=inode+1
        if (mod(inode,nodes).eq.iam)then
         call fsig2(nocc,nvir,a,i,t1,ht1,hiu)
        endif
c -- end parallel task --
c
       enddo
      enddo
      call ga_dgop(913,ht1,lnov, '+')
      next=nxtval(-nodes)
c
      call qexit('itm',0)
      return
      end
