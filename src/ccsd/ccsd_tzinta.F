      subroutine ccsd_tzinta(nocc,nvir,nbf,t1,tzvv,scra,scrb,g_nt2,
     &                       g_nz2,g_tz1,g_tz2,iprt)
      implicit none
      integer nocc,nvir,nbf,lnov,lnvv,g_nt2,g_nz2,g_tz1,g_tz2,g_tmp,
     &        iprt
      integer i,j,a,b,ad1,ad2,ad3,ad4,g_jlo,g_jhi,g_ilo,g_ihi
      double precision t1(nocc,nvir),tzvv(nvir,nvir),scra(nvir,nvir),
     &                 scrb(nvir,nvir)
      Integer Nodes, IAm
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
c
      Nodes = GA_NNodes()
      IAM = GA_NodeID()
C
      lnov=nocc*nvir
      lnvv=nvir*nvir
c
      call dfill(lnvv,0.0d0,tzvv,1)
      call ga_distribution(g_nt2,iam,g_jlo,g_jhi,g_ilo,g_ihi)
      do i=1,nocc
        ad1=(i-1)*nvir
        if (ad1+1.ge.g_ilo.and.ad1+1.le.g_ihi)then
          do j=1,nocc
            ad2=(j-1)*nvir
            if (ad2+1.ge.g_jlo.and.ad2+1.le.g_jhi)then
              call ga_get(g_nt2,ad1+1,ad1+nvir,ad2+1,ad2+nvir,
     &                    scra,nvir)
              call ga_get(g_nz2,ad1+1,ad1+nvir,ad2+1,ad2+nvir,
     &                    scrb,nvir)
              call dgemm('n','t',nvir,nvir,nvir,1.0d0,scra,nvir,
     &                   scrb,nvir,1.0d0,tzvv,nvir)
            endif
          enddo
        endif
      enddo
c
      if (.not.ga_create(MT_DBL,lnov,lnov,'tmp',
     &                   nvir,nvir,g_tmp))
     &     call errquit('ga_create g_tmp failed',0)
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
                scrb(a,b)=2.0d0*scra(a,b)-scra(b,a)
              enddo
              enddo
              call ga_put(g_tmp,ad2+1,ad2+nvir,ad1+1,ad1+nvir,
     &                    scrb,nvir)
            endif
          enddo
        endif
      enddo
c
      call ga_dgemm('n','t',lnov,lnov,lnov,
     &              1.0d00,g_tmp,g_nz2,0.0d00,g_tz1)
crk   call ga_print(g_tz1)
c
      if (.not.ga_destroy(g_tmp))
     &    call errquit('ga_dest g_tmp fail',0)
c
      return
      end
