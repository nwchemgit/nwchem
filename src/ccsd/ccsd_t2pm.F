      subroutine ccsd_ht2pm(basis,nsh,ncor,nocc,nvir,nact,nbf,g_sht2,
     &                      g_nht2,cmo,scra,scrb,offsh)
C     $Id: ccsd_t2pm.F,v 2.14 1999-06-16 00:10:12 d3e129 Exp $
      implicit none
      integer basis,nsh,ncor,nocc,nvir,nact,nbf,g_sht2,g_nht2,
     &        offsh(nsh,nsh)
      double precision cmo(nbf,nbf),scra(nbf*nbf),scrb(nbf*nbf)

#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "rtdb.fh"
#include "ccsd_debug.fh"
c
      integer g_sht2t,g_jlo,g_jhi,g_ilo,g_ihi,
     &        ish,ilo,ihi,jsh,jlo,jhi,ksh,klo,khi,lsh,llo,lhi,
     &        x,y,xy,ymax,nfi,nfj,ipp,imm,ii,jj,iijj,
     &        jjii,i,j,k,l,ad1,ad2,ad3,iam,lnoo,a,b
      integer nxtval
      external nxtval
c
      call qenter('ht2pm',0)
c
      iam=ga_nodeid()
      lnoo=nocc*nocc
c
      call ga_distribution(g_sht2,iam,g_jlo,g_jhi,g_ilo,g_ihi)
      do ish=1,nsh
        if (.not. bas_cn2bfr(basis,ish,ilo,ihi))
     $      call errquit('vvvv: bas_cn2bfr',ish)
        nfi=ihi-ilo+1
        do jsh=1,ish
          if (.not. bas_cn2bfr(basis,jsh,jlo,jhi))
     $        call errquit('vvvv: bas_cn2bfr',jsh)
          nfj=jhi-jlo+1
          do x=1,nfi
            ymax=nfj
            if (ish.eq.jsh)ymax=x
            do y=1,ymax
              xy=offsh(ish,jsh)+(x-1)*nfj+y
              if (xy.ge.g_ilo.and.xy.le.g_ihi)then
                call ga_get(g_sht2,1,lnoo,xy,xy,scra,lnoo)
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
                call ga_put(g_sht2,1,lnoo,xy,xy,scrb,lnoo)
                if (ish.ne.jsh.or.x.ne.y)then
                  xy=offsh(jsh,ish)+(y-1)*nfi+x
                  do ii=1,nocc
                    do jj=1,nocc
                     iijj=(ii-1)*nocc+jj
                     jjii=(jj-1)*nocc+ii
                     scra(iijj)=scrb(jjii)
                    enddo
                  enddo
                  call ga_put(g_sht2,1,lnoo,xy,xy,scra,lnoo)
                endif
              endif
            enddo
          enddo
        enddo
      enddo
      call ga_sync()
c
c ------------------------------------------------------------
c - transform ao indices of ht2 array into the mo basis
c ------------------------------------------------------------

*ga:1:0
      if (.not.ga_create(MT_DBL,nbf*nbf,lnoo,'sht2t',
     &                   nbf*nbf,0,g_sht2t))
     &     call errquit('ga_create g_sht2t failed',0)

      call ga_transpose(g_sht2,g_sht2t)
c - redefine g_sht2
      if (.not.ga_destroy(g_sht2))
     &    call errquit('ga_dest g_sht2 fail',0)
*ga:1:0
      if (.not.ga_create(MT_DBL,nact*nact,lnoo,'sht2',
     &                   nact*nact,0,g_sht2))
     &     call errquit('ga_create g_sht2 failed',0)

      call ga_distribution(g_nht2,iam,g_jlo,g_jhi,g_ilo,g_ihi)
      do i=1,nocc
        ad1=(i-1)*nvir
        if (ad1+1.ge.g_ilo.and.ad1+1.le.g_ihi)then
          do j=1,nocc
            ad2=(j-1)*nvir
            if (ad2+1.ge.g_jlo.and.ad2+1.le.g_jhi)then
            ad3=(i-1)*nocc+j
              call ga_get(g_sht2t,1,nbf*nbf,ad3,ad3,
     &                    scrb,nbf*nbf)
              ad3=0
              do ksh=1,nsh
                if (.not. bas_cn2bfr(basis,ksh,klo,khi))
     &              call errquit('vvvv: bas_cn2bfr',ksh)
                do lsh=1,nsh
                  if (.not. bas_cn2bfr(basis,lsh,llo,lhi))
     &                call errquit('vvvv: bas_cn2bfr',lsh)
                  do k=klo,khi
                    do l=llo,lhi
                      ad3=ad3+1
                      scra((k-1)*nbf+l)=scrb(ad3)
                    enddo
                  enddo
                enddo
              enddo
              call dgemm('n','n',nbf,nact,nbf,1.0d00,scra,nbf,
     &                   cmo(1,ncor+1),nbf,0.0d00,scrb,nbf)
              call dgemm('t','n',nact,nact,nbf,1.0d00,cmo(1,ncor+1),nbf,
     &                   scrb,nbf,0.0d00,scra,nact)
              ad3=(i-1)*nocc+j
              call ga_put(g_sht2,1,nact*nact,ad3,ad3,
     &                    scra,nact*nact)
c
              if (dob(1).eq.2)then
                ad3=0
                do a=nocc+1,nact
                  do b=nocc+1,nact
                    ad3=ad3+1
                    scrb(ad3)=scra((a-1)*nact+b)
                  enddo
                enddo
                call ga_acc(g_nht2,ad2+1,ad2+nvir,ad1+1,ad1+nvir,
     &                      scrb,nvir,1.0d00)
              endif
c
            endif
          enddo
        endif
      enddo
      if (.not.ga_destroy(g_sht2t))
     &    call errquit('ga_dest g_sht2t fail',0)
c
      call qexit('ht2pm',0)
c
      return
      end
      subroutine ccsd_t2pm(basis,nsh,ncor,nocc,nvir,nbf,g_st2,g_nt2,cmo,
     &                     t1,scra,scrb)
      implicit none
      integer basis,nsh,ncor,nocc,nvir,nbf,g_st2,g_nt2
      double precision cmo(nbf,nbf),scra(nbf*nbf),scrb(nbf*nbf),
     &                 t1(nocc,nvir)
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "bas.fh"
#include "rtdb.fh"
c
      integer g_st2t,lnoo,iam
      integer i,j,k,l,a,b,ad1,ad2,ad3,g_jlo,g_jhi,g_ilo,g_ihi,
     &        ish,ilo,ihi,jsh,jlo,jhi,ksh,klo,khi,lsh,llo,lhi,
     &        xy,x,y,ipp,imm,ii,jj,iijj,jjii
c
      call qenter('t2pm',0)
c
      iam=ga_nodeid()
      lnoo=nocc*nocc
c
*ga:1:0
      if (.not.ga_create(MT_DBL,nbf*nbf,lnoo,'st2t',
     &                   nbf*nbf,0,g_st2t))
     &     call errquit('ga_create g_st2t failed',0)
c
      call ga_distribution(g_st2t,iam,g_jlo,g_jhi,g_ilo,g_ihi)
      do i=1,nocc
          do j=1,nocc
            ad3=(i-1)*nocc+j
            if (ad3.ge.g_ilo.and.ad3.le.g_ihi)then
              ad1=(i-1)*nvir
              ad2=(j-1)*nvir
              call ga_get(g_nt2,ad2+1,ad2+nvir,ad1+1,ad1+nvir,
     &                    scra,nvir)
              ad3=0
              do a=1,nvir
                do b=1,nvir
                  ad3=ad3+1
                  scra(ad3)=scra(ad3)+t1(i,a)*t1(j,b)
                enddo
              enddo
              call dgemm('n','t',nvir,nbf,nvir,1.0d00,scra,nvir,
     &                   cmo(1,ncor+nocc+1),nbf,0.0d00,scrb,nvir)
              call dgemm('n','n',nbf,nbf,nvir,1.0d00,cmo(1,ncor+nocc+1),
     &                   nbf,scrb,nvir,0.0d00,scra,nbf)
              ad3=0
              do ksh=1,nsh
                if (.not. bas_cn2bfr(basis,ksh,klo,khi))
     &              call errquit('vvvv: bas_cn2bfr',ksh)
                do lsh=1,nsh
                  if (.not. bas_cn2bfr(basis,lsh,llo,lhi))
     &                call errquit('vvvv: bas_cn2bfr',lsh)
                  do k=klo,khi
                    do l=llo,lhi
                      ad3=ad3+1
                      scrb(ad3)=scra((k-1)*nbf+l)
                    enddo
                  enddo
                enddo
              enddo
              ad3=(i-1)*nocc+j
              call ga_put(g_st2t,1,nbf*nbf,ad3,ad3,
     &                    scrb,nbf*nbf)
            endif
          enddo
      enddo
c
      call ga_transpose(g_st2t,g_st2)
      if (.not.ga_destroy(g_st2t))
     &    call errquit('ga_dest g_st2t fail',0)
c
c ------------------------------------------------------------
c - form t2+/- (see gustavos paper)
c - t2+ = t_ij^ab + t_ij^ba 
c - t2- = t_ij^ab - t_ij^ba 
c - some prefactors absorbed
c ------------------------------------------------------------
      call ga_distribution(g_st2,iam,g_jlo,g_jhi,g_ilo,g_ihi)
      xy=0
      do ish=1,nsh
        if (.not. bas_cn2bfr(basis,ish,ilo,ihi))
     $      call errquit('vvvv: bas_cn2bfr',ish)
        do jsh=1,nsh
          if (.not. bas_cn2bfr(basis,jsh,jlo,jhi))
     $        call errquit('vvvv: bas_cn2bfr',jsh)
          do x=ilo,ihi
            do y=jlo,jhi
              xy=xy+1
              if (xy.ge.g_ilo.and.xy.le.g_ihi)then
                call ga_get(g_st2,1,lnoo,xy,xy,scra,lnoo)
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
                call ga_put(g_st2,1,lnoo,xy,xy,scrb,lnoo)
              endif
            enddo
          enddo
        enddo
      enddo
c
      call qexit('t2pm',0)
c
      return
      end
      subroutine ccsd_sxy(basis,nsh,ncor,nocc,nvir,nact,nbf,g_st2,
     &                    g_sht2,g_c,g_x,offsh,snsi,sisn,lssni,scre,
     &                    mem2,max2e,eri1,eri2,t1,cmo,t1ao,scra,
     &                    scrb,lscr,tol2e,iprt,tklst)
      implicit none
      integer basis,nsh,ncor,nocc,nvir,nact,nbf,g_st2,g_sht2,
     &        g_c,g_x,lssni,
     &        mem2,max2e,lscr,iprt
      integer offsh(nsh,nsh),tklst(nsh*(nsh+1)/2,2)
      double precision tol2e,eri1(max2e),eri2(max2e),scre(mem2),
     &                 t1(nocc*nvir),cmo(nbf,nbf),scra(lscr),
     &                 scrb(lscr),snsi(lssni),sisn(lssni),
     &                 t1ao(nocc*nbf),tx(2)

#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "schwarz.fh"
*rak-s
#include "ccsd_time.fh"
*rak-e
      integer ad1,ish,ilo,ihi,jsh,jlo,jhi,ksh,klo,khi,lsh,llo,lhi,next,
     &        icnt,nfi,nfk,lnijkl,offjl,offlj,lnoo,nfj,nfl,nodes,iam
      double precision schwarz1,schwarz2
c
      call qenter('sxy',0)
c
      lnoo=nocc*nocc
      iam=ga_nodeid()
      nodes=ga_nnodes()
c
c ------------------------------------------------------------
c - work out an index array for shell offset
c ------------------------------------------------------------
      call ga_sync()
      ad1=0
      do ish=1,nsh
        if (.not. bas_cn2bfr(basis,ish,ilo,ihi))
     $      call errquit('vvvv: bas_cn2bfr',ish)
        nfi=ihi-ilo+1
        do jsh=1,nsh
          if (.not. bas_cn2bfr(basis,jsh,jlo,jhi))
     $        call errquit('vvvv: bas_cn2bfr',jsh)
          nfj=jhi-jlo+1
          offsh(ish,jsh)=ad1
          ad1=ad1+nfi*nfj
        enddo
      enddo
      call ga_sync()
c ------------------------------------------------------------
c create t1ao
c ------------------------------------------------------------
      call dgemm('n','t',nbf,nocc,nvir,1.0d00,cmo(1,ncor+nocc+1),
     &           nbf,t1,nocc,0.0d00,t1ao,nbf)
c
c ------------------------------------------------------------
c - loop over the integral generation this is what needs to be
c - efficient, so gather timings for this
c - at the moment we are shell blocked, we could go to atom
c - blocking if this is inefficient
c - note integrals computed 4 times minimal list
c ------------------------------------------------------------
c
      call ga_sync()
      tx(1)=tcgtime()
      if (iam.eq.0.and.iprt.gt.5)
     &    print *,' begin parallel integral generation'
      call ga_zero(g_sht2)
      call ga_zero(g_c)
      call ga_sync
      next=nxtval(nodes)
      do icnt=0,nsh*(nsh+1)/2-1
c // parallel task here
          if (icnt.eq.next)then
            jsh=(tklst(icnt+1,2)-1)/nsh+1
            lsh=tklst(icnt+1,2)-(jsh-1)*nsh
c
            if (.not. bas_cn2bfr(basis,jsh,jlo,jhi))
     $          call errquit('vvvv: bas_cn2bfr',jsh)
            nfj=jhi-jlo+1
            if (.not. bas_cn2bfr(basis,lsh,llo,lhi))
     $         call errquit('directbf: bas_cn2bfr',lsh)
            nfl=lhi-llo+1
c
            call dfill(lnoo*nfj*nfl,0.0d00,scrb,1)
            call dfill(nfj*nfl*nbf*nocc,0.0d00,snsi,1)
            call dfill(nfj*nfl*nbf*nocc,0.0d00,sisn,1)
            do ish=1,nsh
              if (.not. bas_cn2bfr(basis,ish,ilo,ihi))
     $            call errquit('vvvv: bas_cn2bfr',ish)
              nfi=ihi-ilo+1
              do ksh=1,ish
                if (.not. bas_cn2bfr(basis,ksh,klo,khi))
     $              call errquit('vvvv: bas_cn2bfr',ksh)
                nfk=khi-klo+1
                ad1=offsh(ish,ksh)
                call ga_get(g_st2,1,lnoo,ad1+1,ad1+nfi*nfk,scra,lnoo)
                schwarz1=schwarz_shell(ish,jsh)*schwarz_shell(ksh,lsh)
                schwarz2=schwarz_shell(ish,lsh)*schwarz_shell(ksh,jsh)
                if (schwarz1.ge.tol2e.or.schwarz2.ge.tol2e)then
                  lnijkl=nfi*nfj*nfk*nfl
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
                  call t2eri(ilo,ihi,jlo,jhi,klo,khi,llo,lhi,eri1,
     &                       eri2,scra,scrb,lnoo,nocc)
                  if (ish.eq.ksh)then
                    call dscal(lnijkl,0.5d00,eri1,1)
                    call dscal(lnijkl,0.5d00,eri2,1)
                  endif
                  call ccsd_idx1(snsi,sisn,eri1,eri2,t1ao,nfj,nfl,nocc,
     &                           nbf,ilo,ihi,klo,khi,jlo,llo)
                endif
              enddo
            enddo
            ad1=offsh(jsh,lsh)
            call ga_acc(g_sht2,1,lnoo,ad1+1,ad1+nfj*nfl,
     &                  scrb,lnoo,0.5d00)
            offjl=offsh(jsh,lsh)
            offlj=offsh(lsh,jsh)
            call ccsd_idx2(snsi,sisn,cmo,scra,scrb,lscr,
     &                     nfj,nfl,ncor,nocc,nact,nbf,jlo,jhi,llo,lhi,
     &                     offjl,offlj,g_x,g_c)
            next=nxtval(nodes)
          endif
c // end parallel task
c
      enddo
      call ga_sync
      next=nxtval(-nodes)
      call ga_sync
c sync before 3 and 4 index transformation
      call ccsd_idx34(basis,cmo,scra,scrb,nsh,ncor,nocc,nact,nbf,
     &                g_x,g_c)
c
      call ga_sync()
      tx(2)=tcgtime()
      if (iam.eq.0) then
*rak     write(6,*)'Time around main block',tx(2)-tx(1)
         main_block_time = tx(2)-tx(1)
      endif
c
      call qexit('sxy',0)
c
      return
      end
