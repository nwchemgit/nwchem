      subroutine ccsd_pampt3(basis,nsh,ncor,nocc,nvir,nact,nbf,
     &                       tol2e,cmo,eorb,iprt,hiu,giu,habe,gabe,
     &                       hia,t1,ht1,idiis,g_nt2,g_nht2,
     &                       g_ncoul,g_nexch,tklst)
C     $Id: ccsd_pampt3.F,v 2.11 1999-06-16 00:10:12 d3e129 Exp $
      implicit none
      integer basis,nsh,ncor,nocc,nvir,nact,nbf,
     &        iprt,iam,nodes,idiis,g_nt2,g_nht2,
     &        g_ncoul,g_nexch,tklst(nsh*(nsh+1)/2,2)

      double precision tol2e,cmo(nbf,nbf),eorb(nbf),
     &                 hiu(*),giu(*),habe(*),gabe(*),hia(*),
     &                 t1(nocc*nvir),ht1(nocc*nvir)

#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "schwarz.fh"
#include "ccsd_debug.fh"
C
      integer l_snsi,k_snsi,l_sisn,k_sisn,l_eri1,k_eri1,
     &        l_eri2,k_eri2,l_scre,k_scre,l_t1ao,k_t1ao,
     &        l_offsh,k_offsh,l_sc1,k_sc1,l_sc2,k_sc2,
     &        l_sc3,k_sc3,l_sc4,k_sc4,
     &        maxbfsh,max1e,mem1,max2e,mem2,lscr,lssni,
     &        g_x,g_c,g_st2,g_sht2
      logical stat
c
      Nodes = GA_NNodes()
      IAm = GA_NodeID()
c
      if (iam.eq.0.and.iprt.gt.5)print *,' ccsd_vvvv '
      stat = bas_numcont(basis,nsh)
      stat=stat.and.bas_nbf_cn_max(basis,maxbfsh)
      if (.not. stat)
     $   call errquit('ccsd_vvvv: failed to get info', 0)
c
c ------------------------------------------------------------
c create global arrays
c ------------------------------------------------------------
*ga:1:0
      if (.not.ga_create(MT_DBL,nbf*nbf,nocc*nocc,'x',
     &                   nbf*nbf,0,g_x))
     &     call errquit('ccsd_pampt3: ga_create g_x failed',0)
c
*ga:1:0
      if (.not.ga_create(MT_DBL,nbf*nbf,nocc*nocc,'c',
     &                   nbf*nbf,0,g_c))
     &     call errquit('ccsd_pampt3: ga_create g_c failed',0)
c -- global array for ht2 with 2 ao indices
*ga:1:0
      if (.not.ga_create(MT_DBL,nocc*nocc,nbf*nbf,'sht2',
     &                   nocc*nocc,0,g_sht2))
     &     call errquit('ccsd_pampt3: ga_create g_sht2 failed',0)
c
c -- global array for t2 with 2 ao indices
*ga:1:0
      if (.not.ga_create(MT_DBL,nocc*nocc,nbf*nbf,'st2',
     &                   nocc*nocc,0,g_st2))
     &     call errquit('ccsd_pampt3: ga_create g_st2 failed',0)
c
c ------------------------------------------------------------
c - Make t2+/- with 2 ao and 2 occ indices
c ------------------------------------------------------------
      lscr=nbf*nbf
      stat=.true.
      stat=stat.and.ma_push_get(MT_DBL,lscr,'scra1',l_sc1, k_sc1)
      stat=stat.and.ma_push_get(MT_DBL,lscr,'scrb',l_sc2, k_sc2)
      if (.not.stat)call errquit('ma_get ccsd_t2pm ',0)
      call ccsd_t2pm(basis,nsh,ncor,nocc,nvir,nbf,g_st2,g_nt2,cmo,
     &               t1,dbl_mb(k_sc1),dbl_mb(k_sc2))
      stat=stat.and.ma_pop_stack(l_sc2)
      stat=stat.and.ma_pop_stack(l_sc1)
      if (.not.stat)call errquit('ma_pop ccsd_t2pm ',0)
c ------------------------------------------------------------
c generate integrals
c ------------------------------------------------------------
      call int_mem_1e(max1e, mem1)
      call int_mem_2e4c(max2e, mem2)
      mem2=max(mem2,nbf*nbf)
      mem2=max(mem2,mem1)
      lssni=nocc*nbf*maxbfsh*maxbfsh
      lscr=max(nbf*nbf,nocc*nocc*maxbfsh*maxbfsh)
      lscr=max(lscr,maxbfsh*nocc*nbf)
      lscr=max(lscr,maxbfsh*maxbfsh*nbf)
      stat=.true.
      stat=stat.and.ma_push_get(MT_INT,nsh*nsh,'offsh',l_offsh,k_offsh)
      stat=stat.and.ma_push_get(MT_DBL,lssni,'snsi',l_snsi, k_snsi)
      stat=stat.and.ma_push_get(MT_DBL,lssni,'sisn',l_sisn, k_sisn)
      stat=stat.and.ma_push_get(MT_DBL,max2e,'eri1',l_eri1, k_eri1)
      stat=stat.and.ma_push_get(MT_DBL,max2e,'eri2',l_eri2, k_eri2)
      stat=stat.and.ma_push_get(MT_DBL,mem2,'scre',l_scre,k_scre)
      stat=stat.and.ma_push_get(MT_DBL,nocc*nbf,'t1ao',l_t1ao, k_t1ao)
      stat=stat.and.ma_push_get(MT_DBL,lscr,'scra2',l_sc1, k_sc1)
      stat=stat.and.ma_push_get(MT_DBL,lscr,'scrb',l_sc2, k_sc2)
      if (.not.stat)call errquit('ma_get ccsd_sxy ',0)
      call ccsd_sxy(basis,nsh,ncor,nocc,nvir,nact,nbf,g_st2,g_sht2,
     &              g_c,g_x,int_mb(k_offsh),dbl_mb(k_snsi),
     &              dbl_mb(k_sisn),lssni,dbl_mb(k_scre),
     &              mem2,max2e,dbl_mb(k_eri1),
     &              dbl_mb(k_eri2),t1,cmo,dbl_mb(k_t1ao),dbl_mb(k_sc1),
     &              dbl_mb(k_sc2),lscr,tol2e,iprt,tklst)
      stat=stat.and.ma_pop_stack(l_sc2)
      stat=stat.and.ma_pop_stack(l_sc1)
      stat=stat.and.ma_pop_stack(l_t1ao)
      stat=stat.and.ma_pop_stack(l_scre)
      stat=stat.and.ma_pop_stack(l_eri2)
      stat=stat.and.ma_pop_stack(l_eri1)
      stat=stat.and.ma_pop_stack(l_sisn)
      stat=stat.and.ma_pop_stack(l_snsi)
      if (.not.stat)call errquit('ma_pop ccsd_sxy ',0)
c
      if (.not.ga_destroy(g_st2))
     &    call errquit('ga_dest g_st2 fail',0)
c ------------------------------------------------------------
c - Change ht2+/- with 2 ao and 2 occ indices into normal ht2
c ------------------------------------------------------------
      stat=.true.
      lscr=nbf*nbf
      stat=stat.and.ma_push_get(MT_DBL,lscr,'scra3',l_sc1, k_sc1)
      stat=stat.and.ma_push_get(MT_DBL,lscr,'scrb',l_sc2, k_sc2)
      if (.not.stat)call errquit('ma_get ccsd_ht2pm ',0)
      call ccsd_ht2pm(basis,nsh,ncor,nocc,nvir,nact,nbf,g_sht2,g_nht2,
     &                cmo,dbl_mb(k_sc1),dbl_mb(k_sc2),
     &                int_mb(k_offsh))
      stat=stat.and.ma_pop_stack(l_sc2)
      stat=stat.and.ma_pop_stack(l_sc1)
      stat=stat.and.ma_pop_stack(l_offsh)
c ------------------------------------------------------------
c other terms using B intermediate
c ------------------------------------------------------------
      lscr=nact*nact
      stat=stat.and.ma_push_get(MT_DBL,lscr,'scra4',l_sc1, k_sc1)
      stat=stat.and.ma_push_get(MT_DBL,lscr,'scrb',l_sc2, k_sc2)
      stat=stat.and.ma_push_get(MT_DBL,lscr,'scre',l_sc3, k_sc3)
      stat=stat.and.ma_push_get(MT_DBL,lscr,'scrf',l_sc4, k_sc4)
      call ccsd_itm(nocc,nvir,nact,hiu,giu,habe,gabe,hia,
     &              t1,ht1,dbl_mb(k_sc1),dbl_mb(k_sc2),
     &              dbl_mb(k_sc1),dbl_mb(k_sc2),
     &              dbl_mb(k_sc3),dbl_mb(k_sc4),
     &              g_nt2,g_nht2,g_ncoul,g_nexch,g_x,g_c,g_sht2,
     &              iprt)
      stat=stat.and.ma_pop_stack(l_sc4)
      stat=stat.and.ma_pop_stack(l_sc3)
      stat=stat.and.ma_pop_stack(l_sc2)
      stat=stat.and.ma_pop_stack(l_sc1)
c ------------------------------------------------------------
c remove global arrays associated with B intermediate
c ------------------------------------------------------------
      if (.not.ga_destroy(g_sht2))
     &    call errquit('ga_dest g_sht2 fail',0)
c ------------------------------------------------------------
c remove global arrays associated with vvvo integrals (g_x and g_c)
c ------------------------------------------------------------
      if (.not.ga_destroy(g_c))
     &    call errquit('ga_dest g_c fail',0)
      if (.not.ga_destroy(g_x))
     &    call errquit('ga_dest g_x fail',0)
c
      return
      end
      subroutine ccsd_idx1(snsi,sisn,eri1,eri2,t1ao,nfj,nfl,nocc,
     &                     nbf,ilo,ihi,klo,khi,jlo,llo)
      implicit none
      integer nfj,nfl,nocc,nbf,ilo,ihi,klo,khi,jlo,llo
      double precision eri1(nfl,klo:khi,nfj,ilo:ihi),
     &                 eri2(nfj,klo:khi,nfl,ilo:ihi),
     &                 snsi(nbf,nfj,nfl,nocc),
     &                 sisn(nbf,nfl,nfj,nocc),
     &                 t1ao(nbf,nocc)
c
      integer i,j,k,l,ii
c
      call qenter('idx1',0)
      do j=1,nfj
        do i=ilo,ihi
          do l=1,nfl
            do k=klo,khi
              do ii=1,nocc
                snsi(i,j,l,ii)=snsi(i,j,l,ii)+t1ao(k,ii)*eri1(l,k,j,i)
                snsi(k,j,l,ii)=snsi(k,j,l,ii)+t1ao(i,ii)*eri2(j,k,l,i)
c
                sisn(k,l,j,ii)=sisn(k,l,j,ii)+t1ao(i,ii)*eri1(l,k,j,i)
                sisn(i,l,j,ii)=sisn(i,l,j,ii)+t1ao(k,ii)*eri2(j,k,l,i)
              enddo
            enddo
          enddo
        enddo
      enddo
      call qexit('idx1',0)
      return
      end
      subroutine ccsd_idx2(snsi,sisn,cmo,scra,scrb,lscr,
     &                     nfi,nfk,ncor,nocc,nact,nbf,
     &                     ilo,ihi,klo,khi,offik,offki,
     &                     g_x,g_c)
      implicit none
      integer lscr,nfi,nfk,ncor,nocc,nact,nbf,ilo,ihi,klo,khi,
     &        offik,offki,g_x,g_c
      double precision snsi(nbf,nfi,nfk,nocc),
     &                 sisn(nbf,nfk,nfi,nocc),
     &                 cmo(nbf,nbf),
     &                 scra(lscr),scrb(lscr)
c lscr=max(maxbfsh*nocc*nbf,maxbfsh*maxbfsh*nbf)
c
      integer nfik,ii,jj,ik,ad1,ad2
c
      call qenter('idx2',0)
c
      nfik=nfi*nfk
      do ii=1,nocc
        call dgemm('t','n',nact,nfik,nbf,1.0d00,cmo(1,ncor+1),nbf,
     &             snsi(1,1,1,ii),nbf,0.0d00,scra,nact)
        do jj=1,nocc
          do ik=1,nfik
            scrb((jj-1)*nfik+ik)=scra((ik-1)*nact+jj)
          enddo
        enddo
        ad1=offki
        ad2=(ii-1)*nocc
        call ga_put(g_x,ad1+1,ad1+nfik,ad2+1,ad2+nocc,scrb,nfik)
        call dgemm('n','n',nact*nfi,nocc,nfk,1.0d00,scra,nact*nfi,
     &               cmo(klo,ncor+1),nbf,0.0d00,scrb,nact*nfi)
        ad1=(ilo-1)*nact
        ad2=(ii-1)*nocc
        call ga_acc(g_c,ad1+1,ad1+nact*nfi,ad2+1,ad2+nocc,scrb,
     &              nact*nfi,1.0d00)
c
        if (ilo.ne.klo)then
          call dgemm('t','n',nact,nfik,nbf,1.0d00,cmo(1,ncor+1),nbf,
     &               sisn(1,1,1,ii),nbf,0.0d00,scra,nact)
          do jj=1,nocc
            do ik=1,nfik
              scrb((jj-1)*nfik+ik)=scra((ik-1)*nact+jj)
            enddo
          enddo
          ad1=offik
          ad2=(ii-1)*nocc
          call ga_put(g_x,ad1+1,ad1+nfik,ad2+1,ad2+nocc,scrb,nfik)
          call dgemm('n','n',nact*nfk,nocc,nfi,1.0d00,scra,nact*nfk,
     &                 cmo(ilo,ncor+1),nbf,0.0d00,scrb,nact*nfk)
          ad1=(klo-1)*nact
          ad2=(ii-1)*nocc
          call ga_acc(g_c,ad1+1,ad1+nact*nfk,ad2+1,ad2+nocc,scrb,
     &                nact*nfk,1.0d00)
        endif
      enddo
c
      call qexit('idx2',0)
c
      return
      end     
      subroutine ccsd_idx34(basis,cmo,scra,scrb,nsh,
     &                      ncor,nocc,nact,nbf,
     &                      g_x,g_c)
      implicit none
      integer basis,nsh,ncor,nocc,nact,nbf,g_x,g_c
      double precision cmo(nbf,nbf),
     &                 scra(nbf*nbf),scrb(nbf*nbf)
c
#include "bas.fh"
c
      integer ii,jj,ish,jsh,i,j,ilo,ihi,jlo,jhi,iam,ga_nodeid,
     &        g_jlo,g_jhi,g_ilo,g_ihi,ad1,ad2
c
      call qenter('idx34',0)
      iam=ga_nodeid()
      call ga_distribution(g_x,iam,g_jlo,g_jhi,g_ilo,g_ihi)
      do ii=1,nocc
        do jj=1,nocc
          ad2=(ii-1)*nocc+jj
c // here
          if (ad2.ge.g_ilo.and.ad2.le.g_ihi)then
            call ga_get(g_x,1,nbf*nbf,ad2,ad2,scra,nbf*nbf)
            ad1=0
            do ish=1,nsh
              if (.not. bas_cn2bfr(basis,ish,ilo,ihi))
     $            call errquit('vvvv: bas_cn2bfr',ish)
              do jsh=1,nsh
                if (.not. bas_cn2bfr(basis,jsh,jlo,jhi))
     $              call errquit('vvvv: bas_cn2bfr',jsh)
                do i=ilo,ihi 
                  do j=jlo,jhi 
                    ad1=ad1+1
                    scrb((i-1)*nbf+j)=scra(ad1)
                  enddo
                enddo
              enddo
            enddo
            call dgemm('t','n',nact,nbf,nbf,1.0d00,cmo(1,ncor+1),nbf,
     &                 scrb,nbf,0.0d00,scra,nact)
            call dgemm('n','n',nact,nact,nbf,1.0d00,scra,nact,
     &                 cmo(1,ncor+1),nbf,0.0d00,scrb,nact)
            call ga_put(g_x,1,nact*nact,ad2,ad2,scrb,nact*nact)
c
            call ga_get(g_c,1,nact*nbf,ad2,ad2,scra,nact*nbf)
            call dgemm('n','n',nact,nact,nbf,1.0d00,scra,nact,
     &                 cmo(1,ncor+1),nbf,0.0d00,scrb,nact)
            call ga_put(g_c,1,nact*nact,ad2,ad2,scrb,nact*nact)
          endif
c end //
        enddo
      enddo
c
      call qexit('idx34',0)
c
      return
      end     
      subroutine t2eri(ilo,ihi,jlo,jhi,klo,khi,llo,lhi,eri1,eri2,
     &                 t2,ht2,lnoo,nocc)
      implicit none
      integer ilo,ihi,jlo,jhi,klo,khi,llo,lhi,lnoo,nocc
      double precision eri1(llo:lhi,klo:khi,jlo:jhi,ilo:ihi),
     &                 eri2(jlo:jhi,klo:khi,llo:lhi,ilo:ihi),
     &                 t2(lnoo,klo:khi,ilo:ihi),
     &                 ht2(lnoo,llo:lhi,jlo:jhi)
c
      integer i,j,k,l,kmax,lmax
      integer ipp,imm
      double precision gp,gm
c
      call qenter('t2eri',0)
c
c - we could block the ints together and do an mxm if
c - this performs poorly
c
      do i=ilo,ihi
        do j=jlo,jhi
          kmax=khi
          if (ilo.eq.klo)kmax=i
          do k=klo,kmax
            lmax=lhi
            if (jlo.eq.llo)lmax=j
            do l=llo,lmax
              gp=(eri1(l,k,j,i)+eri2(j,k,l,i))
              gm=(eri1(l,k,j,i)-eri2(j,k,l,i))
              ipp=nocc*(nocc+1)/2
              imm=nocc*(nocc-1)/2
c
c - here is the n**6, note inside a 4 fold loop over nbf which is
c - restricted twice, giving a factor of 1/4
c - the two daxpy give nocc*nocc ops
c
              call daxpy(ipp,gp,t2(1,k,i),1,ht2(1,l,j),1)
              call daxpy(imm,gm,t2(ipp+1,k,i),1,ht2(ipp+1,l,j),1)
            enddo
          enddo
        enddo
      enddo
c
      call qexit('t2eri',0)
c
      return
      end

