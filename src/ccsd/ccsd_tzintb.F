      subroutine ccsd_tzintb(nocc,nvir,nbf,t1,tzoo,scra,scrb,scrc,
     &                       g_nt2,g_nz2,g_tz3,iprt)
      implicit none
      integer nocc,nvir,nbf,lnoo,lnov,g_nt2,g_nz2,g_tz3,g_tmp,iprt
      integer i,j,k,l,a,b,ij,kl,ad1,ad2,ad3,ad4,g_jlo,g_jhi,g_ilo,g_ihi
      double precision t1(nocc,nvir),tzoo(nocc,nocc),scra(nvir,nvir),
     &                 scrb(nvir,nvir),scrc(nocc,nocc),x,ddot
      Integer Nodes, IAm
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
c
      Nodes = GA_NNodes()
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
                tzoo(i,j)=tzoo(i,j)+scra(a,a)
              enddo
            endif
          enddo
        endif
      enddo
c
      if (.not.ga_destroy(g_tmp))
     &    call errquit('ga_dest g_tmp fail',0)
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
              call dfill(lnoo,0.0d00,scrc,1)
              do k=1,nocc
              do l=1,nocc
              kl=(k-1)*nocc+l
                ad3=(k-1)*nvir
                ad4=(l-1)*nvir
                call ga_get(g_nz2,ad4+1,ad4+nvir,ad3+1,ad3+nvir,
     &                      scrb,nvir)
c               do a=1,nvir
c               do b=1,nvir
c                 scrc(k,l)=scrc(k,l)+(scra(a,b)+t1(i,a)*t1(j,b))
c    &                               *scrb(a,b)
c               enddo
c               enddo
              x=ddot(nvir*nvir,scra,1,scrb,1) 
              call ga_put(g_tz3,kl,kl,ij,ij,x,1)
              enddo
              enddo
c             call ga_put(g_tz3,1,lnoo,ij,ij,scrc,lnoo)
            endif
          enddo
        endif
      enddo
c
c     call ga_print(g_tz3)
c
      return
      end
