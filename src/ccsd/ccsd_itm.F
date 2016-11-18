#ifdef CRAYXT
c#warning undeffing nxtval
#undef NXT
#endif      
      subroutine ccsd_itm(nocc,nvir,nact,hiu,giu,habe,gabe,
     &                    hia,t1,ht1,
     &                    scra,scrb,scrc,scrd,scre,scrf,
     &                    g_nt2,g_nht2,g_ncoul,g_nexch,g_x,g_c,g_sht2,
     &                    iprt)
C     $Id$
      implicit none
#include "errquit.fh"
      integer nocc,nvir,nact,g_nt2,g_nht2,g_ncoul,g_nexch,
     &        g_x,g_c,g_sht2,iprt
      double precision hiu(nocc,nocc),giu(nocc,nocc),hia(nocc,nvir)
      double precision habe(nvir,nvir),gabe(nvir,nvir)
      double precision t1(nocc,nvir),ht1(nocc,nvir)
      double precision scra(nact*nact),scrb(nact*nact),
     &                 scrc(nact,nact),scrd(nact,nact),
     &                 scre(nact*nact),scrf(nact*nact)
c
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "rtdb.fh"
#include "msgids.fh"
#include "ccsd_debug.fh"
#include "ccsdps.fh"
c
      integer g_aitm,g_tmt1,g_tmt2,g_tmt,g_tmp,g_tmp2,
     &        g_jlo,g_jhi,g_ilo,g_ihi,
     &        i,j,k,l,ij,a,b,c,ad1,ad2,ad3,ad4,ad5,ad6,lnoo,lnov,lnvv
      integer inode,next
      Integer Nodes, IAm
c
      integer nxtask
      external nxtask
c
      Nodes = GA_NNodes()
      IAm = GA_NodeID()
c
      if (occsdps) then
         call pstat_on(ps_itm)
      else
         call qenter('itm',0)
      endif
c
      lnoo=nocc*nocc
      lnov=nocc*nvir
      lnvv=nvir*nvir
c
c - test x and c by doing its contribution to (2j-k)*(2t-t)
      if (.not.ga_create(MT_DBL,lnoo,lnoo,'aitm',
     &                   nocc,nocc,g_aitm))
     &     call errquit('ccsd_itm: ga_create g_aitm failed',0, GA_ERR)
*ga:1:0
      if (.not.ga_create(MT_DBL,lnvv,lnoo,'tmt1',
     &                   lnvv,0,g_tmt1))
     &     call errquit('ccsd_itm: ga_create g_tmt1 failed',0, GA_ERR)
*ga:1:0
      if (.not.ga_create(MT_DBL,lnvv,lnoo,'tmt2',
     &                   lnvv,0,g_tmt2))
     &     call errquit('ccsd_itm: ga_create g_tmt2 failed',0, GA_ERR)
c
      call dfill(lnvv,0.0d00,habe,1)
      call ga_distribution(g_nexch,iam,g_jlo,g_jhi,g_ilo,g_ihi)
      do i=1,nocc
        ad1=(i-1)*nact
        if (ad1+1.ge.g_ilo.and.ad1+1.le.g_ihi)then
          do j=1,nocc
            ad2=(j-1)*nact
            if (ad2+1.ge.g_jlo.and.ad2+1.le.g_jhi)then
c-----------------------------------------------------------------------
c  S7
c-----------------------------------------------------------------------
            call ga_get(g_nexch,ad2+1,ad2+nact,ad1+1,ad1+nact,
     &                  scrc,nact)
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
            call ga_get(g_sht2,1,nact*nact,ij,ij,scre,nact*nact)
            call dscal(nact*nact,2.0d00,scre,1)
c
c H2
            do k=1,nocc
                hiu(k,i)=hiu(k,i)+scre((k-1)*nact+j)+scre((k-1)*nact+j)
     &                          -scre((j-1)*nact+k)
            enddo
          endif
        enddo
      enddo
c
      call dcopy(lnoo,hiu,1,giu,1)
      call dcopy(lnvv,habe,1,gabe,1)
      call dfill(lnov,0.0d00,hia,1)
c
      if (iprt.gt.5.and.iam.eq.0)print *,'forming aitm'
      call ga_distribution(g_x,iam,g_jlo,g_jhi,g_ilo,g_ihi)
      do i=1,nocc
        do j=1,nocc
          ij=(i-1)*nocc+j
          if (ij.ge.g_ilo.and.ij.le.g_ihi)then
            ad1=(i-1)*nact
            ad2=(j-1)*nact
c
c-----------------------------------------------------------------------
c  A1, A2, A3 and A4 intermediates
c-----------------------------------------------------------------------
c
c A1
            call dfill(nact*nact,0.0d00,scra,1)
            if (doa(1).eq.2)then
              call ga_get(g_nexch,ad2+1,ad2+nocc,ad1+1,ad1+nocc,
     &                    scrb,nocc)
              call daxpy(nocc*nocc,0.5d00,scrb,1,scra,1)
            endif
c---------------------------------------------------------
c  S5, D1, D6 and D7 
c---------------------------------------------------------
            call ga_get(g_nexch,ad2+nocc+1,ad2+nact,ad1+1,ad1+nact,
     &                  scre,nvir)
            call dfill(nact*nact,0.0d00,scrf,1)
c
c D1,D7
            ad1=0
            do a=1,nvir
              do b=nocc+1,nact
                ad1=ad1+1
                scrf(ad1)=scrf(ad1)+0.5d00*scre((b-1)*nvir+a)
                do k=1,nocc
                  scrf(ad1)=scrf(ad1)-scre((k-1)*nvir+a)*t1(k,b-nocc)
                enddo
              enddo
            enddo
c
            ad1=(i-1)*nocc+j
            call ga_get(g_c,1,nact*nact,ad1,ad1,scrb,nact*nact)
            call dcopy(nact*nact,scrb,1,scre,1)
c
c H8
            do a=nocc+1,nact
              hia(j,a-nocc)=hia(j,a-nocc)-scre((a-1)*nact+i)
            enddo
c
c G4
            if (i.eq.j) then
            do b=nocc+1,nact
            do a=nocc+1,nact
              gabe(a-nocc,b-nocc)=gabe(a-nocc,b-nocc)+scre((a-1)*nact+b)
     &                            +scre((a-1)*nact+b)
            enddo
            enddo
            endif
c
            ad1=(i-1)*nocc+j
            call ga_get(g_x,1,nact*nact,ad1,ad1,scrb,nact*nact)
            call dcopy(nact*nact,scrb,1,scre,1)
c
c H8
            do a=nocc+1,nact
              hia(j,a-nocc)=hia(j,a-nocc)+scre((i-1)*nact+a)
     &                      +scre((i-1)*nact+a)
            enddo
c
c S5
            do a=nocc+1,nact
              ht1(j,a-nocc)=ht1(j,a-nocc)+scre((i-1)*nact+a)
     &                      +scre((i-1)*nact+a)-scre((a-1)*nact+i)
            enddo
c
c G2
            do k=1,nocc
              giu(k,j)=giu(k,j)+scre((i-1)*nact+k)
     &                      +scre((i-1)*nact+k)-scre((k-1)*nact+i)
            enddo
c
c G4
            if (i.eq.j) then
            do b=nocc+1,nact
            do a=nocc+1,nact
              gabe(a-nocc,b-nocc)=gabe(a-nocc,b-nocc)-scre((a-1)*nact+b)
            enddo
            enddo
            endif
c
c D6,D7
            ad1=0
            do a=nocc+1,nact
              do b=nocc+1,nact
                ad1=ad1+1
                scrf(ad1)=scrf(ad1)+scre((b-1)*nact+a)
                do k=1,nocc
                  scrf(ad1)=scrf(ad1)-scre((b-1)*nact+k)*t1(k,a-nocc)
                  scrf(ad1)=scrf(ad1)-scre((k-1)*nact+a)*t1(k,b-nocc)
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
                scra(ad1)=scra(ad1)+0.5d00*scrb((k-1)*nact+l)
              enddo
            enddo
            ad1=(j-1)*nocc+i
            call ga_get(g_x,1,nact*nact,ad1,ad1,scrb,nact*nact)
            ad1=0
            do k=1,nocc
              do l=1,nocc
                ad1=ad1+1
                scra(ad1)=scra(ad1)+0.5d00*scrb((l-1)*nact+k)
              enddo
            enddo
            endif
c
            call ga_get(g_sht2,1,nact*nact,ij,ij,scrb,nact*nact)
            call dfill(nact*nact,0.0d00,scre,1)
            call daxpy(nact*nact,2.0d00,scrb,1,scre,1)
c
c S6
            do a=nocc+1,nact
              ht1(j,a-nocc)=ht1(j,a-nocc)+scre((i-1)*nact+a)
     &                      +scre((i-1)*nact+a)-scre((a-1)*nact+i)
            enddo
c
c A4
            if (doa(4).eq.2)then
            ad1=0
            do k=1,nocc
              do l=1,nocc
                ad1=ad1+1
                scra(ad1)=scra(ad1)+scrb((k-1)*nact+l)
              enddo
            enddo
            endif
            call ga_put(g_aitm,1,lnoo,ij,ij,scra,lnoo)
c-----------------------------------------------------------------------
c  B2 and B3 intermediate and their contribution to D3
c-----------------------------------------------------------------------
            if (dob(2).eq.2)then
            ad1=nocc+1
            call dgemm('n','n',nvir,nvir,nocc,1.0d00,scrb(ad1),
     &                 nact,t1,nocc,0.0d00,scra,nvir)
            ad1=nocc*nact+1
            call dgemm('t','n',nvir,nvir,nocc,1.0d00,t1,nocc,
     &                 scrb(ad1),nact,1.0d00,scra,nvir)
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
      call ga_dgop(msg_cc_term1,hiu,lnoo, '+')
      call ga_dgop(msg_cc_term2,habe,lnvv, '+')
      call ga_dgop(msg_cc_term3,giu,lnoo, '+')
      call ga_dgop(msg_cc_term4,hia,lnov, '+')
      call ga_dgop(msg_cc_term5,gabe,lnvv, '+')
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
              call dgemm('n','n',nvir,nvir,nvir,1.0d00,gabe,
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
      call ga_sync
#ifdef NXT
      next=nxtask(nodes, 1)
#endif
      do a=1,nvir
          do i=1,nocc
            ad1=(i-1)*nvir+a
        inode=inode+1
#ifdef NXT
        if (inode.eq.next)then
#else
            if (mod(inode, nodes) .eq. iam) then
#endif
            call ga_get(g_nt2,1,lnov,ad1,ad1,scre,lnov)
            call dgemm('n','n',nvir,nocc,nocc,-1.0d00,scre,
     &                 nvir,giu,nocc,0.0d00,scrf,nvir)
            call ga_acc(g_nht2,1,lnov,ad1,ad1,scrf,lnov,1.0d00)
#ifdef NXT
      next=nxtask(nodes, 1)
#endif
          endif
        enddo
      enddo
#ifdef NXT
      call ga_sync
      next=nxtask(-nodes, 1)
#endif
      call ga_sync

c
      if (.not.ga_destroy(g_tmt2))
     &    call errquit('ga_dest g_tmt2 fail',0, GA_ERR)
      if (.not.ga_destroy(g_tmt1))
     &    call errquit('ga_dest g_tmt1 fail',0, GA_ERR)
      if (.not.ga_destroy(g_aitm))
     &    call errquit('ga_dest g_aitm fail',0, GA_ERR)
c
c-----------------------------------------------------------------------
c  J1, J2, J3 and K1, K2, K3 intermediates
c-----------------------------------------------------------------------
      call ga_distribution(g_x,iam,g_jlo,g_jhi,g_ilo,g_ihi)
      do i=1,nocc
        do j=1,nocc
          ij=(i-1)*nocc+j
          if (ij.ge.g_ilo.and.ij.le.g_ihi)then
            ad1=(i-1)*nact
            ad2=(j-1)*nact
c
            call ga_get(g_nexch,ad2+1,ad2+nact,ad1+1,ad1+nact,
     &                  scrd,nact)
            call dgemm('n','n',nvir,nvir,nocc,-1.0d00,scrd(nocc+1,1),
     &                 nact,t1,nocc,0.0d00,scrc(nocc+1,nocc+1),nact)
            call daxpy(nact*nact,1.0d00,scrc,1,scrd,1)
            call ga_acc(g_x,1,nact*nact,ij,ij,scrd,nact*nact,1.0d00)
          endif
        enddo
      enddo
c
      call ga_distribution(g_c,iam,g_jlo,g_jhi,g_ilo,g_ihi)
      do i=1,nocc
        do j=1,nocc
          ij=(i-1)*nocc+j
          if (ij.ge.g_ilo.and.ij.le.g_ihi)then
            ad1=(i-1)*nact
            ad2=(j-1)*nact
            call ga_get(g_ncoul,ad2+1,ad2+nact,ad1+1,ad1+nact,
     &                    scrd,nact)
            call dgemm('n','n',nvir,nvir,nocc,-1.0d00,scrd(nocc+1,1),
     &                 nact,t1,nocc,0.0d00,scrc(nocc+1,nocc+1),nact)
            call daxpy(nact*nact,1.0d00,scrc,1,scrd,1)
            call ga_acc(g_c,1,nact*nact,ij,ij,scrd,nact*nact,1.0d00)
c
          endif
        enddo
      enddo
c
      if (.not.ga_create(MT_DBL,lnov,lnov,'tmt',
     &                   nvir,nvir,g_tmt))
     &     call errquit('ccsd_itm: ga_create g_tmt failed',0, GA_ERR)
      if (.not.ga_create(MT_DBL,lnov,lnov,'tmp',
     &                   nvir,nvir,g_tmp))
     &     call errquit('ccsd_itm: ga_create g_tmp failed',0, GA_ERR)
      if (.not.ga_create(MT_DBL,lnov,lnov,'tmp2',
     &                   nvir,nvir,g_tmp2))
     &     call errquit('ccsd_itm: ga_create g_tmp2 failed',0, GA_ERR)
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
        ad1=(i-1)*nact
        if (ad1+1.ge.g_ilo.and.ad1+1.le.g_ihi)then
          do j=1,nocc
            ad2=(j-1)*nact
            if (ad2+1.ge.g_jlo.and.ad2+1.le.g_jhi)then
            ad3=(i-1)*nvir
            ad4=(j-1)*nvir
            call ga_get(g_nexch,ad2+nocc+1,ad2+nact,ad1+nocc+1,ad1+nact,
     &                  scra,nvir)
            call ga_put(g_tmp,ad3+1,ad3+nvir,ad4+1,ad4+nvir,
     &                  scra,nvir)
            endif
          enddo
        endif
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
              call ga_get(g_c,1,nact*nact,ij,ij,scrb,nact*nact)
              ad3=0
              do a=nocc+1,nact
                do b=nocc+1,nact
                  ad3=ad3+1
                  scra(ad3)=scrb((a-1)*nact+b)
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
            ad1=(i-1)*nact
            ad2=(j-1)*nact
            call dfill(nact*nact,0.0d00,scra,1)
              call ga_get(g_x,1,nact*nact,ij,ij,scrb,nact*nact)
              call daxpy(nact*nact,1.0d00,scrb,1,scra,1)
              call ga_get(g_c,1,nact*nact,ij,ij,scrb,nact*nact)
              call daxpy(nact*nact,-0.5d00,scrb,1,scra,1)
            ad1=0
            do a=nocc+1,nact
              do b=nocc+1,nact
                ad1=ad1+1
                scrb(ad1)=scra((a-1)*nact+b)
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
        ad1=(i-1)*nact
        if (ad1+1.ge.g_ilo.and.ad1+1.le.g_ihi)then
          do j=1,nocc
            ad2=(j-1)*nact
            if (ad2+1.ge.g_jlo.and.ad2+1.le.g_jhi)then
            ad3=(i-1)*nvir
            ad4=(j-1)*nvir
            call ga_get(g_nexch,ad2+nocc+1,ad2+nact,ad1+nocc+1,ad1+nact,
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
     &    call errquit('ga_dest g_tmp2 fail',0, GA_ERR)
      if (.not.ga_destroy(g_tmp))
     &    call errquit('ga_dest g_tmp fail',0, GA_ERR)
      if (.not.ga_destroy(g_tmt))
     &    call errquit('ga_dest g_tmt fail',0, GA_ERR)
c
      inode=-1
      call ga_sync
      next=nxtask(nodes, 1)
      do a=1,nvir
       do b=1,nvir
c
c -- parallel task --
        inode=inode+1
        if (inode.eq.next)then
c
         call ccsd_fsig1(a,b,nocc,nvir,nact,iprt,
     $                   t1,ht1,habe)
         next=nxtask(nodes, 1)
        endif
c -- end parallel task --
c
       enddo
      enddo
      call ga_sync
      next=nxtask(-nodes, 1)
      call ga_sync
c
c -- add some trivial terms to ht1 (nothing to ht2)
      inode=-1
      do a=1,nvir
       do i=1,nocc
c
c -- parallel task --
        inode=inode+1
        if (mod(inode,nodes).eq.iam)then
         call ccsd_fsig2(nocc,nvir,a,i,t1,ht1,hiu)
        endif
c -- end parallel task --
c
       enddo
      enddo
      call ga_dgop(msg_cc_term6,ht1,lnov, '+')
c
      if (occsdps) then
         call pstat_off(ps_itm)
      else
         call qexit('itm',0)
      endif
      return
      end
      subroutine ccsd_fsig1(av,bv,nocc,nvir,nbf,iprt,t1,ht1,habe)
      implicit none
      integer av,bv,nocc,nvir,nbf,iprt
      double precision t1(nocc,nvir),ht1(nocc,nvir),habe(nvir,nvir)
#include "ccsd_debug.fh"
#include "ccsd_len.fh"
! term s1 n**3
      if (dos(1).ne.0)then
       call daxpy(nocc,habe(bv,av),t1(1,av),1,ht1(1,bv),1)
      endif
      return
      end
      subroutine ccsd_fsig2(nocc,nvir,bv,u,t1,ht1,hiu)
      implicit none
      integer nocc,nvir,bv,u
      double precision t1(nocc,nvir),ht1(nocc,nvir),hiu(nocc,nvir)
      double precision ddot
      external ddot
#include "ccsd_debug.fh"
#include "ccsd_len.fh"
! term s2 n**3
      if (dos(2).ne.0)then
       ht1(u,bv)=ht1(u,bv)-ddot(nocc,hiu(1,u),1,t1(1,bv),1)
      endif
      return
      end
