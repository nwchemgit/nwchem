      subroutine ccsd_aosig(basis,nsh,nocc,nvir,nbf,cmo,scre,eri1,eri2,
     &                      t1,hz1,tzvv,tzao,hzao,t1ao,scra,scrb,
     &                      mem2,max2e,tol2e,g_nt2,g_nz2,g_nhz2,g_tz1,
     &                      g_tz2)
      implicit none
#include "errquit.fh"
C     $Id$
      integer nsh,nocc,nvir,nbf,next,icnt,basis,max2e,mem2
      integer ish,jsh,ksh,lsh,ilo,ihi,jlo,jhi,klo,khi,llo,lhi
      integer g_ilo,g_ihi,g_jlo,g_jhi
      integer nfi,nfj,nfk,nfl,lnijkl,lnof
      integer p,q,a,b,c,i,j,k,ad1,ad2,ad3,ad4
      double precision cmo(nbf,nbf),scra(nbf,nbf),scrb(nbf*nbf)
      double precision t1(nocc,nvir),hz1(nocc,nvir),tzvv(nvir,nvir),
     &                 tzao(nbf,nbf),hzao(nbf,nbf),t1ao(nbf,nvir)
      double precision tol2e,eri1(max2e),scre(mem2),schwarz1,
     &                 eri2(max2e),schwarz2
      integer g_nt2,g_nz2,g_nhz2,g_tz1,g_tz2,g_tmp,g_tmp2,g_tmp3,g_tmp4
      integer nxtask
      external nxtask
c
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "schwarz.fh"
c
      Integer Nodes, IAm
      Nodes = GA_NNodes()
      IAM = GA_NodeID()
C
      lnof=nocc*nbf
c
      if (.not.ga_create(MT_DBL,lnof,lnof,'tmp',
     &                   nbf,nbf,g_tmp))
     &     call errquit('ga_create g_tmp failed',0, GA_ERR)
      if (.not.ga_create(MT_DBL,lnof,lnof,'tmp2',
     &                   nbf,nbf,g_tmp2))
     &     call errquit('ga_create g_tmp2 failed',0, GA_ERR)
      if (.not.ga_create(MT_DBL,lnof,lnof,'tmp3',
     &                   nbf,nbf,g_tmp3))
     &     call errquit('ga_create g_tmp3 failed',0, GA_ERR)
      if (.not.ga_create(MT_DBL,lnof,lnof,'tmp4',
     &                   nbf,nbf,g_tmp4))
     &     call errquit('ga_create g_tmp4 failed',0, GA_ERR)
c
c     Back transform the intermediates
c
      call dfill(nbf*nbf,0.0d00,tzao,1)
      do p=1,nbf
      do q=1,nbf
      do a=1,nvir
      do b=1,nvir
        tzao(p,q)=tzao(p,q)+tzvv(a,b)*cmo(p,nocc+a)*cmo(q,nocc+b)
      enddo
      enddo
      enddo
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
     &                    scrb,nvir)
              call dfill(nbf*nbf,0.0d00,scra,1)
              do p=1,nbf
              do q=1,nbf
              do a=1,nvir
              do b=1,nvir
                scra(p,q)=scra(p,q)+scrb((b-1)*nvir+a)*cmo(p,nocc+a)
     &                    *cmo(q,nocc+b)
              enddo
              enddo
              enddo
              enddo
              ad3=(i-1)*nbf
              ad4=(j-1)*nbf
              call ga_put(g_tmp,ad4+1,ad4+nbf,ad3+1,ad3+nbf,
     &                    scra,nbf)
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
     &                    scrb,nvir)
              call dfill(nbf*nbf,0.0d00,scra,1)
              do p=1,nbf
              do q=1,nbf
              do a=1,nvir
              do b=1,nvir
                scra(p,q)=scra(p,q)+scrb((b-1)*nvir+a)*cmo(p,nocc+a)
     &                    *cmo(q,nocc+b)
              enddo
              enddo
              enddo
              enddo
              ad3=(i-1)*nbf
              ad4=(j-1)*nbf
              call ga_put(g_tmp2,ad4+1,ad4+nbf,ad3+1,ad3+nbf,
     &                    scra,nbf)
            endif
          enddo
        endif
      enddo
c
      call dfill(nbf*nvir,0.0d00,t1ao,1)
      do p=1,nbf
      do c=1,nvir
      do k=1,nocc
        t1ao(p,c)=t1ao(p,c)+t1(k,c)*cmo(p,k)
      enddo
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
     &                    scrb,nvir)
              call dfill(nbf*nbf,0.0d00,scra,1)
              do p=1,nbf
              do q=1,nbf
              do a=1,nvir
              do b=1,nvir
                scra(p,q)=scra(p,q)+scrb((b-1)*nvir+a)*cmo(p,nocc+a)
     &                    *t1ao(q,b)
              enddo
              enddo
              enddo
              enddo
              ad3=(i-1)*nbf
              ad4=(j-1)*nbf
              call ga_put(g_tmp3,ad4+1,ad4+nbf,ad3+1,ad3+nbf,
     &                    scra,nbf)
            endif
          enddo
        endif
      enddo
c
c     Generate integrals
c
      call ga_zero(g_tmp4)
      nodes=ga_nnodes()
      call ga_sync()
      next=nxtask(nodes, 1)
      icnt=-1
      call dfill(nbf*nbf,0.0d00,hzao,1)
      do ish=1,nsh
        if (.not. bas_cn2bfr(basis,ish,ilo,ihi))
     $      call errquit('aosig: bas_cn2bfr',ish, BASIS_ERR)
        nfi=ihi-ilo+1
        do jsh=1,nsh
          if (.not. bas_cn2bfr(basis,jsh,jlo,jhi))
     $        call errquit('aosig: bas_cn2bfr',jsh, BASIS_ERR)
          nfj=jhi-jlo+1
c
c // parallel task here
          icnt=icnt+1
          if (icnt.eq.next)then
            do ksh=1,nsh
              if (.not. bas_cn2bfr(basis,ksh,klo,khi))
     $            call errquit('aosig: bas_cn2bfr',ksh, BASIS_ERR)
              nfk=khi-klo+1
              do lsh=1,nsh
                If (.not. bas_cn2bfr(basis,lsh,llo,lhi))
     $              call errquit('aosig: bas_cn2bfr',lsh, BASIS_ERR)
                nfl=lhi-llo+1
                lnijkl=nfi*nfj*nfk*nfl
                schwarz1=schwarz_shell(ish,jsh)*schwarz_shell(ksh,lsh)
                schwarz2=schwarz_shell(ish,lsh)*schwarz_shell(ksh,jsh)
                if (schwarz1.ge.tol2e)then
                    call int_2e4c(basis, ish, jsh, basis, ksh, lsh,
     $                            mem2, scre, max2e, eri1 )
                  else
                    call dfill(lnijkl,0.0d00,eri1,1)
                  endif
                  if (schwarz2.ge.tol2e)then
                    call int_2e4c(basis, ish, lsh, basis, ksh, jsh,
     $                            mem2, scre, max2e, eri2 )
                  else
                    call dfill(lnijkl,0.0d00,eri2,1)
                  endif
                  call z2eri(ilo,ihi,jlo,jhi,klo,khi,llo,lhi,eri1,eri2,
     &                       tzao,hzao,nbf,hz1,scra,scrb,cmo,nocc,nvir,
     &                       g_tmp,g_tmp2,g_tmp3,g_tmp4)
              enddo
            enddo
          next=nxtask(nodes, 1)
          endif
        enddo
      enddo
c
      call ga_sync()
      next=nxtask(-nodes, 1)
      call ga_sync
      call ga_dgop(910,hzao,nbf*nbf, '+')
c
      if(iam.eq.0)then
        do p=1,nbf
        do q=1,nbf
        do k=1,nocc
        do c=1,nvir
          hz1(k,c)=hz1(k,c)+hzao(p,q)*cmo(p,k)*cmo(q,nocc+c)
        enddo
        enddo
        enddo
        enddo
      endif
c
      call ga_zero(g_nhz2)
      call ga_distribution(g_tmp4,iam,g_jlo,g_jhi,g_ilo,g_ihi)
      do i=1,nocc
        ad1=(i-1)*nbf
        if (ad1+1.ge.g_ilo.and.ad1+1.le.g_ihi)then
          do j=1,nocc
            ad2=(j-1)*nbf
            if (ad2+1.ge.g_jlo.and.ad2+1.le.g_jhi)then
              call ga_get(g_tmp4,ad2+1,ad2+nbf,ad1+1,ad1+nbf,
     &                    scra,nbf)
              call dfill(nvir*nvir,0.0d00,scrb,1)
              do p=1,nbf
              do q=1,nbf
              do a=1,nvir
              do b=1,nvir
                scrb((b-1)*nvir+a)=scrb((b-1)*nvir+a)+scra(p,q)
     &                             *cmo(p,nocc+b)*cmo(q,nocc+a)
              enddo
              enddo
              enddo
              enddo
              ad3=(i-1)*nvir
              ad4=(j-1)*nvir
              call ga_acc(g_nhz2,ad4+1,ad4+nvir,ad3+1,ad3+nvir,
     &                    scrb,nvir,1.0d0)
            endif
          enddo
        endif
      enddo
c
      if (.not.ga_destroy(g_tmp4))
     &    call errquit('ga_dest g_tmp4 fail',0, GA_ERR)
      if (.not.ga_destroy(g_tmp3))
     &    call errquit('ga_dest g_tmp3 fail',0, GA_ERR)
      if (.not.ga_destroy(g_tmp2))
     &    call errquit('ga_dest g_tmp2 fail',0, GA_ERR)
      if (.not.ga_destroy(g_tmp))
     &    call errquit('ga_dest g_tmp fail',0, GA_ERR)
      return
      end
      subroutine z2eri(ilo,ihi,jlo,jhi,klo,khi,llo,lhi,eri1,eri2,
     &                 tzao,hzao,nbf,hz1,scra,scrb,cmo,nocc,nvir,
     &                       g_tmp,g_tmp2,g_tmp3,g_tmp4)
      implicit none
      integer ilo,ihi,jlo,jhi,klo,khi,llo,lhi,i,j,nbf,nocc,nvir
      integer p,q,r,s,c,ad1,ad2
      integer g_tmp,g_tmp2,g_tmp3,g_tmp4
      double precision eri1(llo:lhi,klo:khi,jlo:jhi,ilo:ihi),
     &                 eri2(jlo:jhi,klo:khi,llo:lhi,ilo:ihi),
     &                 tzao(nbf,nbf),hzao(nbf,nbf),scra(nbf,nbf),
     &                 scrb(nbf,nbf)
      double precision hz1(nocc,nvir),cmo(nbf,nbf)
c
#include "global.fh"
c
      Integer Nodes, IAm
      Nodes = GA_NNodes()
      IAM = GA_NodeID()
C
      do p=ilo,ihi
      do q=jlo,jhi
      do r=klo,khi
      do s=llo,lhi
      hzao(r,s)=hzao(r,s)-(2.0d0*eri1(s,r,q,p)-eri2(q,r,s,p))
     &   *tzao(q,p)
      enddo
      enddo
      enddo
      enddo
c
      do i=1,nocc
        ad1=(i-1)*nbf
          do j=1,nocc
            ad2=(j-1)*nbf
              call ga_get(g_tmp,ad1+1,ad1+nbf,ad2+1,ad2+nbf,
     &                    scra,nbf)
              do p=ilo,ihi
              do q=jlo,jhi
              do r=klo,khi
              do s=llo,lhi
              do c=1,nvir
                hz1(j,c)=hz1(j,c)-eri2(q,r,s,p)*cmo(r,i)
     &                            *cmo(s,nocc+c)*scra(q,p)
              enddo
              enddo
              enddo
              enddo
              enddo
          enddo
      enddo
c
      do i=1,nocc
        ad1=(i-1)*nbf
          do j=1,nocc
            ad2=(j-1)*nbf
              call ga_get(g_tmp2,ad1+1,ad1+nbf,ad2+1,ad2+nbf,
     &                    scra,nbf)
              do p=ilo,ihi
              do q=jlo,jhi
              do r=klo,khi
              do s=llo,lhi
              do c=1,nvir
                hz1(j,c)=hz1(j,c)+eri1(s,r,q,p)*cmo(r,i)
     &                            *cmo(s,nocc+c)*scra(q,p)
              enddo
              enddo
              enddo
              enddo
              enddo
          enddo
      enddo
c
      do i=1,nocc
        ad1=(i-1)*nbf
          do j=1,nocc
            ad2=(j-1)*nbf
              call ga_get(g_tmp3,ad2+1,ad2+nbf,ad1+1,ad1+nbf,
     &                    scra,nbf)
              call dfill(nbf*nbf,0.0d00,scrb,1)
              do p=ilo,ihi
              do q=jlo,jhi
              do r=klo,khi
              do s=llo,lhi
                scrb(r,s)=scrb(r,s)+eri2(q,r,s,p)*scra(p,q)
              enddo
              enddo
              enddo
              enddo
              call ga_acc(g_tmp4,ad2+1,ad2+nbf,ad1+1,ad1+nbf,
     &                    scrb,nbf,1.0d0)
          enddo
      enddo
c
      return
      end
