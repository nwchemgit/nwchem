      subroutine ccsd_z2pm(basis,nsh,ncor,nocc,nvir,nbf,g_sz2,g_nz2,cmo,
     &                     scra,scrb)
C     $Id$
      implicit none
#include "errquit.fh"
      integer basis,nsh,ncor,nocc,nvir,nbf,g_sz2,g_nz2
      double precision cmo(nbf,nbf),scra(nbf*nbf),scrb(nbf*nbf)
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "bas.fh"
#include "rtdb.fh"
#include "ccsdps.fh"
c
      integer g_sz2t,lnoo
      integer i,j,k,l,ad1,ad2,ad3,g_jlo,g_jhi,g_ilo,g_ihi,
     &        ish,ilo,ihi,jsh,jlo,jhi,ksh,klo,khi,lsh,llo,lhi,
     &        xy,x,y,ipp,imm,ii,jj,iijj,jjii
c
      Integer Nodes, IAm
      Nodes = GA_NNodes()
      IAM = GA_NodeID()
c
      if (occsdps) then
         call pstat_on(ps_z2pm)
      else
         call qenter('z2pm',0)
      endif
C
      lnoo=nocc*nocc
c
*ga:1:0
      if (.not.ga_create(MT_DBL,nbf*nbf,lnoo,'sz2t',
     &                   nbf*nbf,0,g_sz2t))
     &     call errquit('ccsd_z2pm: ga_create g_sz2t failed',0, GA_ERR)
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
              call dgemm('n','t',nvir,nbf,nvir,1.0d00,scra,nvir,
     &                   cmo(1,ncor+nocc+1),nbf,0.0d00,scrb,nvir)
              call dgemm('n','n',nbf,nbf,nvir,1.0d00,cmo(1,ncor+nocc+1),
     &                   nbf,scrb,nvir,0.0d00,scra,nbf)
              ad3=0
              do ksh=1,nsh
                if (.not. bas_cn2bfr(basis,ksh,klo,khi))
     &              call errquit('vvvv: bas_cn2bfr',ksh, BASIS_ERR)
                do lsh=1,nsh
                  if (.not. bas_cn2bfr(basis,lsh,llo,lhi))
     &                call errquit('vvvv: bas_cn2bfr',lsh, BASIS_ERR)
                  do k=klo,khi
                    do l=llo,lhi
                      ad3=ad3+1
                      scrb(ad3)=scra((k-1)*nbf+l)
                    enddo
                  enddo
                enddo
              enddo
              ad3=(i-1)*nocc+j
              call ga_put(g_sz2t,1,nbf*nbf,ad3,ad3,
     &                    scrb,nbf*nbf)
            endif
          enddo
        endif
      enddo
c
      call ga_transpose(g_sz2t,g_sz2)
      if (.not.ga_destroy(g_sz2t))
     &    call errquit('ga_dest g_sz2t fail',0, GA_ERR)
c
c ------------------------------------------------------------
c - form z2+/- (see gustavos paper)
c - z2+ = t_ij^ab + t_ij^ba 
c - z2- = t_ij^ab - t_ij^ba 
c - some prefactors absorbed
c ------------------------------------------------------------
      call ga_distribution(g_sz2,iam,g_jlo,g_jhi,g_ilo,g_ihi)
      xy=0
      do ish=1,nsh
        if (.not. bas_cn2bfr(basis,ish,ilo,ihi))
     $      call errquit('vvvv: bas_cn2bfr',ish, BASIS_ERR)
        do jsh=1,nsh
          if (.not. bas_cn2bfr(basis,jsh,jlo,jhi))
     $        call errquit('vvvv: bas_cn2bfr',jsh, BASIS_ERR)
          do x=ilo,ihi
            do y=jlo,jhi
              xy=xy+1
              if (xy.ge.g_ilo.and.xy.le.g_ihi)then
                call ga_get(g_sz2,1,lnoo,xy,xy,scra,lnoo)
                ipp=0
                imm=nocc*(nocc+1)/2
                do ii=1,nocc
                  do jj=1,ii-1
                    iijj=(ii-1)*nocc+jj
                    jjii=(jj-1)*nocc+ii
                    ipp=ipp+1
                    imm=imm+1
                    scrb(ipp)=(scra(iijj)+scra(jjii))
                    scrb(imm)=(scra(iijj)-scra(jjii))
                  enddo
                  iijj=(ii-1)*nocc+ii
                  ipp=ipp+1
                  scrb(ipp)=scra(iijj)+scra(iijj)
                enddo
                if (x.eq.y)then
                  call dscal(lnoo,0.25d00,scrb,1)
                else
                  call dscal(lnoo,0.5d00,scrb,1)
                endif
                call ga_put(g_sz2,1,lnoo,xy,xy,scrb,lnoo)
              endif
            enddo
          enddo
        enddo
      enddo
c
      if (occsdps) then
         call pstat_off(ps_z2pm)
      else
         call qexit('z2pm',0)
      endif
c
      return
      end
      subroutine ccsd_hz2pm(basis,nsh,ncor,nocc,nvir,nbf,g_shz2,
     &                      g_nhz2,cmo,
     &                      scra,scrb,offsh)
      implicit none
#include "errquit.fh"
      integer basis,nsh,ncor,nocc,nvir,nbf,g_shz2,g_nhz2,
     &        offsh(nsh,nsh)
      double precision cmo(nbf,nbf),scra(nbf*nbf),scrb(nbf*nbf)

#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "bas.fh"
#include "rtdb.fh"
#include "ccsd_debug.fh"
#include "ccsdps.fh"
c
      integer g_shz2t,g_jlo,g_jhi,g_ilo,g_ihi,
     &        ish,ilo,ihi,jsh,jlo,jhi,ksh,klo,khi,lsh,llo,lhi,
     &        x,y,xy,ymax,nfi,nfj,ipp,imm,ii,jj,iijj,
     &        jjii,i,j,k,l,ad1,ad2,ad3,lnoo,a,b
c
      Integer Nodes, IAm
      Nodes = GA_NNodes()
      IAM = GA_NodeID()
c
      if (occsdps) then
         call pstat_on(ps_hz2pm)
      else
         call qenter('hz2pm',0)
      endif
c
      lnoo=nocc*nocc
c
      call ga_distribution(g_shz2,iam,g_jlo,g_jhi,g_ilo,g_ihi)
      do ish=1,nsh
        if (.not. bas_cn2bfr(basis,ish,ilo,ihi))
     $      call errquit('vvvv: bas_cn2bfr',ish, BASIS_ERR)
        nfi=ihi-ilo+1
        do jsh=1,ish
          if (.not. bas_cn2bfr(basis,jsh,jlo,jhi))
     $        call errquit('vvvv: bas_cn2bfr',jsh, BASIS_ERR)
          nfj=jhi-jlo+1
          do x=1,nfi
            ymax=nfj
            if (ish.eq.jsh)ymax=x
            do y=1,ymax
              xy=offsh(ish,jsh)+(x-1)*nfj+y
              if (xy.ge.g_ilo.and.xy.le.g_ihi)then
                call ga_get(g_shz2,1,lnoo,xy,xy,scra,lnoo)
                ipp=0
                imm=nocc*(nocc+1)/2
                do ii=1,nocc
                  do jj=1,ii-1
                    ipp=ipp+1
                    imm=imm+1
                    iijj=(ii-1)*nocc+jj
                    jjii=(jj-1)*nocc+ii
                    scrb(iijj)=scra(ipp)+scra(imm)
                    scrb(jjii)=scra(ipp)-scra(imm)
                  enddo
                  iijj=(ii-1)*nocc+ii
                  ipp=ipp+1
                  scrb(iijj)=scra(ipp)
                enddo
                call ga_put(g_shz2,1,lnoo,xy,xy,scrb,lnoo)
                if (ish.ne.jsh.or.x.ne.y)then
                  xy=offsh(jsh,ish)+(y-1)*nfi+x
                  do ii=1,nocc
                    do jj=1,nocc
                     iijj=(ii-1)*nocc+jj
                     jjii=(jj-1)*nocc+ii
                     scra(iijj)=scrb(jjii)
                    enddo
                  enddo
                  call ga_put(g_shz2,1,lnoo,xy,xy,scra,lnoo)
                endif
              endif
            enddo
          enddo
        enddo
      enddo
      call ga_sync()
c
c ------------------------------------------------------------
c - transform ao indices of hz2 array into the mo basis
c ------------------------------------------------------------

*ga:1:0
      if (.not.ga_create(MT_DBL,nbf*nbf,lnoo,'shz2t',
     &                   nbf*nbf,0,g_shz2t))
     &     call errquit('ccsd_hz2pm: ga_create g_shz2t failed',0,
     &       GA_ERR)

      call ga_transpose(g_shz2,g_shz2t)
c - redefine g_shz2
      if (.not.ga_destroy(g_shz2))
     &    call errquit('ga_dest g_shz2 fail',0, GA_ERR)
*ga:1:0
      if (.not.ga_create(MT_DBL,nbf*nbf,lnoo,'shz2',
     &                   nbf*nbf,0,g_shz2))
     &     call errquit('ccsd_hz2pm: ga_create g_shz2 failed',0, GA_ERR)

      call ga_distribution(g_nhz2,iam,g_jlo,g_jhi,g_ilo,g_ihi)
      do i=1,nocc
        ad1=(i-1)*nvir
        if (ad1+1.ge.g_ilo.and.ad1+1.le.g_ihi)then
          do j=1,nocc
            ad2=(j-1)*nvir
            if (ad2+1.ge.g_jlo.and.ad2+1.le.g_jhi)then
            ad3=(i-1)*nocc+j
              call ga_get(g_shz2t,1,nbf*nbf,ad3,ad3,
     &                    scrb,nbf*nbf)
              ad3=0
              do ksh=1,nsh
                if (.not. bas_cn2bfr(basis,ksh,klo,khi))
     &              call errquit('vvvv: bas_cn2bfr',ksh, BASIS_ERR)
                do lsh=1,nsh
                  if (.not. bas_cn2bfr(basis,lsh,llo,lhi))
     &                call errquit('vvvv: bas_cn2bfr',lsh, BASIS_ERR)
                  do k=klo,khi
                    do l=llo,lhi
                      ad3=ad3+1
                      scra((k-1)*nbf+l)=scrb(ad3)
                    enddo
                  enddo
                enddo
              enddo
              call dgemm('n','n',nbf,nbf,nbf,1.0d00,scra,nbf,
     &                   cmo(1,ncor+1),nbf,0.0d00,scrb,nbf)
              call dgemm('t','n',nbf,nbf,nbf,1.0d00,cmo(1,ncor+1),nbf,
     &                   scrb,nbf,0.0d00,scra,nbf)
              ad3=(i-1)*nocc+j
              call ga_put(g_shz2,1,nbf*nbf,ad3,ad3,
     &                    scra,nbf*nbf)
c
              if (dob(1).eq.2)then
                ad3=0
                do a=nocc+1,nbf
                  do b=nocc+1,nbf
                    ad3=ad3+1
                    scrb(ad3)=-scra((a-1)*nbf+b)
                  enddo
                enddo
                call ga_acc(g_nhz2,ad2+1,ad2+nvir,ad1+1,ad1+nvir,
     &                      scrb,nvir,1.0d00)
              endif
c
            endif
          enddo
        endif
      enddo
      if (.not.ga_destroy(g_shz2t))
     &    call errquit('ga_dest g_shz2t fail',0, GA_ERR)
c
      if (occsdps) then
         call pstat_off(ps_hz2pm)
      else
         call qexit('hz2pm',0)
      endif
c
      return
      end
