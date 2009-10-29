      subroutine ccsd_pampt3(basis,nsh,ncor,nocc,nvir,nact,nbf,
     &                       tol2e,cmo,eorb,iprt,hiu,giu,habe,gabe,
     &                       hia,t1,ht1,idiis,g_nt2,g_nht2,
     &                       g_ncoul,g_nexch,tklst)
C     $Id: ccsd_pampt3.F,v 2.16 2006-11-03 14:36:10 edo Exp $
      implicit none
#include "errquit.fh"
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
#include "ccsdps.fh"
C
      integer l_snsi,k_snsi,l_sisn,k_sisn,l_eri1,k_eri1,
     &        l_eri2,k_eri2,l_scre,k_scre,l_t1ao,k_t1ao,
     &        l_offsh,k_offsh,l_sc1,k_sc1,l_sc2,k_sc2,
     &        l_sc3,k_sc3,l_sc4,k_sc4,intblk,kscr,
     &        maxbfsh,max1e,mem1,max2e,mem2,lscr,lssni,
     &        g_x,g_c,g_st2,g_sht2,l_sch,k_sch,l_ish,k_ish,nrints,
     &        l_shinf,k_shinf
      logical stat
      integer max_sht2_blk
      integer ilo,ihi,jlo,jhi,nbfdim,ish,jsh,nfi,nfj
c
      Nodes = GA_NNodes()
      IAm = GA_NodeID()
c
      if (iam.eq.0.and.iprt.gt.5)print *,' ccsd_vvvv '
      stat = bas_numcont(basis,nsh)
      stat=stat.and.bas_nbf_cn_max(basis,maxbfsh)
      if (.not. stat)
     $   call errquit('ccsd_vvvv: failed to get info', 0, BASIS_ERR)
c
c ------------------------------------------------------------
c create global arrays
c ------------------------------------------------------------
*ga:1:0
      if (.not.ga_create(MT_DBL,nbf*nbf,nocc*nocc,'x',
     &                   nbf*nbf,0,g_x))
     &     call errquit('ccsd_pampt3: ga_create g_x failed',0, GA_ERR)
c
*ga:1:0
      if (.not.ga_create(MT_DBL,nbf*nbf,nocc*nocc,'c',
     &                   nbf*nbf,0,g_c))
     &     call errquit('ccsd_pampt3: ga_create g_c failed',0, GA_ERR)
c -- global array for ht2 with 2 ao indices
*ga:1:0
      if (.not.ga_create(MT_DBL,nocc*nocc,nbf*nbf,'sht2',
     &                   nocc*nocc,0,g_sht2))
     &     call errquit('ccsd_pampt3: ga_create g_sht2 failed',0,
     &       GA_ERR)
c
c -- global array for t2 with 2 ao indices
*ga:1:0
      nbfdim=0
      do ish=1,nsh
        if (.not. bas_cn2bfr(basis,ish,ilo,ihi))
     $      call errquit('vvvv: bas_cn2bfr',ish, BASIS_ERR)
        nfi=ihi-ilo+1
        do jsh=1,ish
          if (.not. bas_cn2bfr(basis,jsh,jlo,jhi))
     $        call errquit('vvvv: bas_cn2bfr',jsh, BASIS_ERR)
          nfj=jhi-jlo+1
          nbfdim=nbfdim+nfi*nfj
        enddo
      enddo
      if (.not.ga_create(MT_DBL,nocc*nocc,nbfdim,'st2',
     &                   nocc*nocc,0,g_st2))
     &     call errquit('ccsd_pampt3: ga_create g_st2 failed',0, GA_ERR)
c
c ------------------------------------------------------------
c - Make t2+/- with 2 ao and 2 occ indices
c ------------------------------------------------------------
      lscr=nbf*nbf
      stat=.true.
      stat=stat.and.ma_push_get(MT_DBL,lscr,'scra1',l_sc1, k_sc1)
      stat=stat.and.ma_push_get(MT_DBL,lscr,'scrb',l_sc2, k_sc2)
      if (.not.stat)call errquit('ma_get ccsd_t2pm ',0, MA_ERR)
      call ccsd_t2pm(basis,nsh,ncor,nocc,nvir,nbf,g_st2,g_nt2,cmo,
     &               t1,dbl_mb(k_sc1),dbl_mb(k_sc2),nbfdim)
      stat=stat.and.ma_pop_stack(l_sc2)
      stat=stat.and.ma_pop_stack(l_sc1)
      if (.not.stat)call errquit('ma_pop ccsd_t2pm ',0, MA_ERR)
c ------------------------------------------------------------
c generate integrals
c ------------------------------------------------------------
      call int_mem_1e(max1e, mem1)
      call int_mem_2e4c(max2e, mem2)
      mem2=max(mem2,nbf*nbf)
      mem2=max(mem2,mem1)
      lssni=nocc*nbf*maxbfsh*maxbfsh
      max_sht2_blk=nocc*nocc*maxbfsh*maxbfsh
      lscr=max(nbf*nbf,max_sht2_blk)
      lscr=max(lscr,maxbfsh*nocc*nbf)
      lscr=max(lscr,maxbfsh*maxbfsh*nbf)
      kscr=max(lscr,2*max_sht2_blk)
      stat=.true.
      stat=stat.and.ma_push_get(MT_INT,nsh*nsh*2,'off',l_offsh,k_offsh)
      stat=stat.and.ma_push_get(MT_INT,nsh*3,'shinf',l_shinf,k_shinf)
      stat=stat.and.ma_push_get(MT_DBL,lssni,'snsi',l_snsi, k_snsi)
      stat=stat.and.ma_push_get(MT_DBL,lssni,'sisn',l_sisn, k_sisn)
      stat=stat.and.ma_push_get(MT_DBL,mem2,'scre',l_scre,k_scre)
      stat=stat.and.ma_push_get(MT_DBL,nocc*nbf,'t1ao',l_t1ao, k_t1ao)
      stat=stat.and.ma_push_get(MT_DBL,lscr,'scrb',l_sc2, k_sc2)
      stat=stat.and.ma_push_get(MT_DBL,kscr,'scra2',l_sc1, k_sc1)
      stat=stat.and.ma_push_get(MT_INT,nsh*nsh+2,'ish',l_ish,k_ish)
c 
c Get memory availability and allocate multiple blocks if possible
c Also need an index array for blocks of ish of length max2e*intblk
c (in case we have a whole series of <ss|ss> integrals
c
c Maximum memory needed for <ij|kl> block is max2e * 2
c Use intblk variable to define multiple blocks
c
      intblk=ma_inquire_avail(MT_DBL)/(2*max2e)
      intblk=max(1,intblk)
      nrints=max2e*intblk
      stat=stat.and.ma_push_get(MT_DBL,nrints,'eri1',l_eri1,k_eri1)
      stat=stat.and.ma_push_get(MT_DBL,nrints,'eri2',l_eri2,k_eri2)
      call ccsd_sxy(basis,nsh,ncor,nocc,nvir,nact,nbf,g_st2,g_sht2,
     &              g_c,g_x,int_mb(k_offsh),dbl_mb(k_snsi),
     &              dbl_mb(k_sisn),lssni,dbl_mb(k_scre),
     &              mem2,max2e,dbl_mb(k_eri1),
     &              dbl_mb(k_eri2),t1,cmo,dbl_mb(k_t1ao),dbl_mb(k_sc1),
     &              dbl_mb(k_sc2),lscr,kscr,tol2e,iprt,tklst,nrints,
     &              int_mb(k_ish),int_mb(k_shinf),max_sht2_blk)
      stat=stat.and.ma_pop_stack(l_eri2)
      stat=stat.and.ma_pop_stack(l_eri1)
      stat=stat.and.ma_pop_stack(l_ish)
      stat=stat.and.ma_pop_stack(l_sc1)
      stat=stat.and.ma_pop_stack(l_sc2)
      stat=stat.and.ma_pop_stack(l_t1ao)
      stat=stat.and.ma_pop_stack(l_scre)
      stat=stat.and.ma_pop_stack(l_sisn)
      stat=stat.and.ma_pop_stack(l_snsi)
      if (.not.stat)
     &   call errquit('ccsd_pampt3: ma_pop_stack ccsd_sxy ',0, MA_ERR)
c
      if (.not.ga_destroy(g_st2))
     &    call errquit('ga_dest g_st2 fail',0, GA_ERR)
c ------------------------------------------------------------
c - Change ht2+/- with 2 ao and 2 occ indices into normal ht2
c ------------------------------------------------------------
      stat=.true.
      lscr=nbf*nbf
      stat=stat.and.ma_push_get(MT_DBL,lscr,'scra3',l_sc1, k_sc1)
      stat=stat.and.ma_push_get(MT_DBL,lscr,'scrb',l_sc2, k_sc2)
      if (.not.stat)call errquit('ma_get ccsd_ht2pm ',0, MA_ERR)
      call ccsd_ht2pm(basis,nsh,ncor,nocc,nvir,nact,nbf,g_sht2,g_nht2,
     &                cmo,dbl_mb(k_sc1),dbl_mb(k_sc2),
     &                int_mb(k_offsh))
      stat=stat.and.ma_pop_stack(l_sc2)
      stat=stat.and.ma_pop_stack(l_sc1)
      stat=stat.and.ma_pop_stack(l_shinf)
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
     &    call errquit('ga_dest g_sht2 fail',0, GA_ERR)
c ------------------------------------------------------------
c remove global arrays associated with vvvo integrals (g_x and g_c)
c ------------------------------------------------------------
      if (.not.ga_destroy(g_c))
     &    call errquit('ga_dest g_c fail',0, GA_ERR)
      if (.not.ga_destroy(g_x))
     &    call errquit('ga_dest g_x fail',0, GA_ERR)
c
      return
      end
      subroutine ccsd_idx1(snsi,sisn,eri1,eri2,t1ao,nfj,nfl,nocc,
     &                     nbf,idx,idim,nsh,shinf)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "schwarz.fh"
#include "eaf.fh"
#include "ccsdps.fh"
      integer nfj,nfl,nocc,nbf,idim,nsh
      integer idx(idim),shinf(nsh,3)
      double precision eri1(*),
     &                 eri2(*),
     &                 snsi(nbf,nfj,nfl,nocc),
     &                 sisn(nbf,nfl,nfj,nocc),
     &                 t1ao(nbf,nocc)
c
      integer i,j,k,l,ii,ncnt,ia,ka
      integer nfi,nfk,ilo,ihi,klo,khi,ish,ksh
      integer ict,off1,off2
      integer mlkj,mkl,mkj
     
      double precision int1, int2, scale
c
      if (occsdps) then
         call pstat_on(ps_idx1)
      else
         call qenter('idx1',0)
      endif
c
      ict=1
      do ncnt=3,idim
        ksh=(idx(ncnt)-1)/nsh+1
        ish=idx(ncnt)-(ksh-1)*nsh
        nfk=shinf(ksh,1)
        klo=shinf(ksh,2)
        khi=shinf(ksh,3)
        nfi=shinf(ish,1)
        ilo=shinf(ish,2)
        ihi=shinf(ish,3)
        if (ish.eq.ksh) then
          call ccsd_idx1_wrk2(snsi,sisn,eri1(ict),eri2(ict),t1ao,nfj,
     &                        nfl,nocc,nbf,ilo,ihi,klo,khi)
        else
          call ccsd_idx1_wrk1(snsi,sisn,eri1(ict),eri2(ict),t1ao,nfj,
     &                        nfl,nocc,nbf,ilo,ihi,klo,khi)
        endif
        ict=ict+nfk*nfl*nfi*nfj
      enddo
      if (occsdps) then
         call pstat_off(ps_idx1)
      else
         call qexit('idx1',0)
      endif
      return
      end
      subroutine ccsd_idx1_wrk1(snsi,sisn,eri1,eri2,t1ao,nfj,nfl,nocc,
     &                     nbf,ilo,ihi,klo,khi)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "ccsdps.fh"
      integer nfj,nfl,nocc,nbf,ilo,ihi,klo,khi
      integer i,j,k,l,ii
      double precision eri1(nfl,klo:khi,nfj,ilo:ihi),
     &                 eri2(nfj,klo:khi,nfl,ilo:ihi),
     &                 snsi(nbf,nfj,nfl,nocc),
     &                 sisn(nbf,nfl,nfj,nocc),
     &                 t1ao(nbf,nocc)
      double precision int1,int2
      call pstat_on(ps_idx11)
      do l=1,nfl
        do j=1,nfj
          do i=ilo,ihi
            do k=klo,khi
              int1=eri1(l,k,j,i)
              int2=eri2(j,k,l,i)
              do ii=1,nocc
               snsi(i,j,l,ii)=snsi(i,j,l,ii)+t1ao(k,ii)*int1
               snsi(k,j,l,ii)=snsi(k,j,l,ii)+t1ao(i,ii)*int2
               sisn(k,l,j,ii)=sisn(k,l,j,ii)+t1ao(i,ii)*int1
               sisn(i,l,j,ii)=sisn(i,l,j,ii)+t1ao(k,ii)*int2
              enddo
            enddo
          enddo
        enddo
      enddo
      call pstat_off(ps_idx11)
      return
      end
      subroutine ccsd_idx1_wrk2(snsi,sisn,eri1,eri2,t1ao,nfj,nfl,nocc,
     &                     nbf,ilo,ihi,klo,khi)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "ccsdps.fh"
      integer nfj,nfl,nocc,nbf,ilo,ihi,klo,khi
      integer i,j,k,l,ii
      double precision eri1(nfl,klo:khi,nfj,ilo:ihi),
     &                 eri2(nfj,klo:khi,nfl,ilo:ihi),
     &                 snsi(nbf,nfj,nfl,nocc),
     &                 sisn(nbf,nfl,nfj,nocc),
     &                 t1ao(nbf,nocc)
      double precision int1,int2
      call pstat_on(ps_idx12)
      do l=1,nfl
        do j=1,nfj
          do i=ilo,ihi
            do k=klo,khi
              int1=eri1(l,k,j,i)*0.5d0
              int2=eri2(j,k,l,i)*0.5d0
              do ii=1,nocc
               snsi(i,j,l,ii)=snsi(i,j,l,ii)+t1ao(k,ii)*int1
               snsi(k,j,l,ii)=snsi(k,j,l,ii)+t1ao(i,ii)*int2
               sisn(k,l,j,ii)=sisn(k,l,j,ii)+t1ao(i,ii)*int1
               sisn(i,l,j,ii)=sisn(i,l,j,ii)+t1ao(k,ii)*int2
              enddo
            enddo
          enddo
        enddo
      enddo
      call pstat_off(ps_idx12)
      return
      end
      subroutine ccsd_idx2(snsi,sisn,cmo,scra,scrb,lscr,
     &                     nfi,nfk,ncor,nocc,nact,nbf,
     &                     ilo,ihi,klo,khi,offik,offki,
     &                     g_x,g_c)
      implicit none
#include "ccsdps.fh"
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
      if (occsdps) then
         call pstat_on(ps_idx2)
      else
         call qenter('idx2',0)
      endif
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
      if (occsdps) then
         call pstat_off(ps_idx2)
      else
         call qexit('idx2',0)
      endif
c
      return
      end     
      subroutine ccsd_idx34(basis,cmo,scra,scrb,nsh,
     &                      ncor,nocc,nact,nbf,
     &                      g_x,g_c)
      implicit none
#include "errquit.fh"
      integer basis,nsh,ncor,nocc,nact,nbf,g_x,g_c
      double precision cmo(nbf,nbf),
     &                 scra(nbf*nbf),scrb(nbf*nbf)
c
#include "bas.fh"
#include "ccsdps.fh"
c
      integer ii,jj,ish,jsh,i,j,ilo,ihi,jlo,jhi,iam,ga_nodeid,
     &        g_jlo,g_jhi,g_ilo,g_ihi,ad1,ad2
c
      if (occsdps) then
         call pstat_on(ps_idx34)
      else
         call qenter('idx34',0)
      endif
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
     $            call errquit('vvvv: bas_cn2bfr',ish, BASIS_ERR)
              do jsh=1,nsh
                if (.not. bas_cn2bfr(basis,jsh,jlo,jhi))
     $              call errquit('vvvv: bas_cn2bfr',jsh, BASIS_ERR)
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
      if (occsdps) then
         call pstat_off(ps_idx34)
      else
         call qexit('idx34',0)
      endif
c
      return
      end     
      subroutine t2eri(idim,idx,jlo,jhi,nfj,llo,lhi,nfl,nsh,eri1,eri2,
     &                 t2,ht2,lnoo,nocc,offsh,nbf,g_st2,
     &                 shinf,max_st2_blk,snsi,sisn,t1ao)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "schwarz.fh"
#include "eaf.fh"
#include "ccsdps.fh"
      integer idim,jlo,jhi,nsh,llo,lhi,lnoo,nocc
      integer idx(idim)
      integer nbf,g_st2,max_st2_blk
      integer offsh(nsh,nsh,2),shinf(nsh,3)
      integer nfj,nfl,ksh,ish,ilo,ihi,klo,khi,nfk,nfi
      double precision eri1(*),
     &                 eri2(*),
     &                 t2(0:(2*max_st2_blk-1)),
     &                 ht2(lnoo,(lhi-llo+1),(jhi-jlo+1)),
     &                 snsi(nbf,nfj,nfl,nocc),
     &                 sisn(nbf,nfl,nfj,nocc),
     &                 t1ao(nbf,nocc)
c
      integer i,j,k,l,kmax,lmax,ik,ad1
      integer ipp,imm,ippp,ict,off1,off2,off3
      integer mlkj,mkl,mkj,ii,jj
      integer kmsh,imsh,nmfi,nmfk,nbhandle,blk1,blk2
      integer kmlo,kmhi,imlo,imhi
      double precision gp,gm
      logical keqi,jeql
c
      nbhandle=4321+ga_nodeid()
c
      if (occsdps) then
         call pstat_on(ps_t2eri)
      else
         call qenter('t2eri',0)
      endif
c
      jeql=jlo.eq.llo
      ipp=nocc*(nocc+1)/2
      ippp=ipp+1
      imm=nocc*(nocc-1)/2
      ict=1
c
c Get first block ik
c
      ksh=(idx(3)-1)/nsh+1
      ish=idx(3)-(ksh-1)*nsh
      nfk=shinf(ksh,1)
      klo=shinf(ksh,2)
      khi=shinf(ksh,3)
      nfi=shinf(ish,1)
      ilo=shinf(ish,2)
      ihi=shinf(ish,3)
      ad1=offsh(ish,ksh,2)
      call ga_get(g_st2,1,lnoo,ad1+1,ad1+nfi*nfk,t2(0),lnoo)
      blk1=0
      blk2=max_st2_blk
c
c Now loop over remaining blocks
c
      do ik=4,idim
c
c Prefetch block ik, start using previous block
c
        kmsh=(idx(ik)-1)/nsh+1
        imsh=idx(ik)-(kmsh-1)*nsh
        nmfk=shinf(kmsh,1)
        kmlo=shinf(kmsh,2)
        kmhi=shinf(kmsh,3)
        nmfi=shinf(imsh,1)
        imlo=shinf(imsh,2)
        imhi=shinf(imsh,3)
        ad1=offsh(imsh,kmsh,2)
c       call pstat_on(ps_t2erin)
        call ga_nbget(g_st2,1,lnoo,ad1+1,ad1+nmfi*nmfk,t2(blk2),lnoo,
     &                nbhandle)
c       call pstat_off(ps_t2erin)
c
c  Do the actual work
c
        keqi=ish.eq.ksh
        if (keqi) then
          call ccsd_idx1_wrk2(snsi,sisn,eri1(ict),eri2(ict),t1ao,nfj,
     &                        nfl,nocc,nbf,ilo,ihi,klo,khi)
        else
          call ccsd_idx1_wrk1(snsi,sisn,eri1(ict),eri2(ict),t1ao,nfj,
     &                        nfl,nocc,nbf,ilo,ihi,klo,khi)
        endif
        if (keqi.and.jeql) then
           call t2eri_wrk_aleq(eri1(ict),eri2(ict),t2(blk1),ht2,
     &                         nfi,nfj,nfk,nfl,nocc,lnoo)
        elseif (keqi) then
           call t2eri_wrk_ikeq(eri1(ict),eri2(ict),t2(blk1),ht2,
     &                         nfi,nfj,nfk,nfl,nocc,lnoo)
        elseif (jeql) then
           call t2eri_wrk_jleq(eri1(ict),eri2(ict),t2(blk1),ht2,
     &                         nfi,nfj,nfk,nfl,nocc,lnoo)
        else
           call t2eri_wrk_noeq(eri1(ict),eri2(ict),t2(blk1),ht2,
     &                         nfi,nfj,nfk,nfl,nocc,lnoo)
        endif
c
c Increase offset of integral block to next ik block
c
        ict=ict+nfi*nfj*nfk*nfl
c
c Replace ish and ksh data with new ones from imsh and kmsh
c
        ish=imsh
        ksh=kmsh
        nfi=nmfi
        nfk=nmfk
        ilo=imlo
        klo=kmlo
        ihi=imhi
        khi=kmhi
c
c Wait for non-blocking get to finish and replace blk with correct pointer
c
c       call pstat_on(ps_t2erih)
        call ga_nbwait(nbhandle)
c       call pstat_off(ps_t2erih)
        blk1=-blk1+max_st2_blk
        blk2=-blk2+max_st2_blk
      enddo
c
c  Do the actual work on final block
c
      keqi=ish.eq.ksh
      if (keqi) then
         call ccsd_idx1_wrk2(snsi,sisn,eri1(ict),eri2(ict),t1ao,nfj,
     &                       nfl,nocc,nbf,ilo,ihi,klo,khi)
      else
         call ccsd_idx1_wrk1(snsi,sisn,eri1(ict),eri2(ict),t1ao,nfj,
     &                       nfl,nocc,nbf,ilo,ihi,klo,khi)
      endif
      if (keqi.and.jeql) then
         call t2eri_wrk_aleq(eri1(ict),eri2(ict),t2(blk1),ht2,
     &                       nfi,nfj,nfk,nfl,nocc,lnoo)
      elseif (keqi) then
         call t2eri_wrk_ikeq(eri1(ict),eri2(ict),t2(blk1),ht2,
     &                       nfi,nfj,nfk,nfl,nocc,lnoo)
      elseif (jeql) then
         call t2eri_wrk_jleq(eri1(ict),eri2(ict),t2(blk1),ht2,
     &                       nfi,nfj,nfk,nfl,nocc,lnoo)
      else
         call t2eri_wrk_noeq(eri1(ict),eri2(ict),t2(blk1),ht2,
     &                       nfi,nfj,nfk,nfl,nocc,lnoo)
      endif
c
      if (occsdps) then
         call pstat_off(ps_t2eri)
      else
         call qexit('t2eri',0)
      endif
c
      return
      end
      subroutine t2eri_wrk_noeq(eri1,eri2,t2,ht2,ni,nj,nk,nl,nocc,lnoo)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
      integer ni,nj,nk,nl,nocc,lnoo
      double precision eri1(nl,nk,nj,ni), eri2(nj,nk,nl,ni)
      double precision t2(lnoo,nk,ni), ht2(lnoo,nl,nj)
      integer ipp,ippp,imm,i,j,k,l
      double precision gp, gm
      ipp=nocc*(nocc+1)/2
      ippp=ipp+1
      imm=nocc*(nocc-1)/2
      do i=1,ni
        do j=1,nj
          do k=1,nk
            do l=1,nl
              gp=(eri1(l,k,j,i)+eri2(j,k,l,i))*0.5d0
              gm=(eri1(l,k,j,i)-eri2(j,k,l,i))*0.5d0
              call daxpy(ipp,gp,t2(1,k,i),1,ht2(1,l,j),1)
              call daxpy(imm,gm,t2(ippp,k,i),1,ht2(ippp,l,j),1)
            enddo
          enddo
        enddo
      enddo
      return
      end
      subroutine t2eri_wrk_jleq(eri1,eri2,t2,ht2,ni,nj,nk,nl,nocc,lnoo)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
      integer ni,nj,nk,nl,nocc,lnoo
      double precision eri1(nl,nk,nj,ni), eri2(nj,nk,nl,ni)
      double precision t2(lnoo,nk,ni), ht2(lnoo,nl,nj)
      integer ipp,ippp,imm,i,j,k,l
      double precision gp, gm
      ipp=nocc*(nocc+1)/2
      ippp=ipp+1
      imm=nocc*(nocc-1)/2
      do i=1,ni
        do j=1,nj
          do k=1,nk
            do l=1,j
              gp=(eri1(l,k,j,i)+eri2(j,k,l,i))*0.5d0
              gm=(eri1(l,k,j,i)-eri2(j,k,l,i))*0.5d0
              call daxpy(ipp,gp,t2(1,k,i),1,ht2(1,l,j),1)
              call daxpy(imm,gm,t2(ippp,k,i),1,ht2(ippp,l,j),1)
            enddo
          enddo
        enddo
      enddo
      return
      end
      subroutine t2eri_wrk_ikeq(eri1,eri2,t2,ht2,ni,nj,nk,nl,nocc,lnoo)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
      integer ni,nj,nk,nl,nocc,lnoo
      double precision eri1(nl,nk,nj,ni), eri2(nj,nk,nl,ni)
      double precision t2(lnoo,nk,ni), ht2(lnoo,nl,nj)
      integer ipp,ippp,imm,i,j,k,l
      double precision gp, gm
      ipp=nocc*(nocc+1)/2
      ippp=ipp+1
      imm=nocc*(nocc-1)/2
      do i=1,ni
        do j=1,nj
          do k=1,i
            do l=1,nl
              gp=(eri1(l,k,j,i)+eri2(j,k,l,i))*0.5d0
              gm=(eri1(l,k,j,i)-eri2(j,k,l,i))*0.5d0
              call daxpy(ipp,gp,t2(1,k,i),1,ht2(1,l,j),1)
              call daxpy(imm,gm,t2(ippp,k,i),1,ht2(ippp,l,j),1)
            enddo
          enddo
        enddo
      enddo
      return
      end
      subroutine t2eri_wrk_aleq(eri1,eri2,t2,ht2,ni,nj,nk,nl,nocc,lnoo)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
      integer ni,nj,nk,nl,nocc,lnoo
      double precision eri1(nl,nk,nj,ni), eri2(nj,nk,nl,ni)
      double precision t2(lnoo,nk,ni), ht2(lnoo,nl,nj)
      integer ipp,ippp,imm,i,j,k,l
      double precision gp, gm
      ipp=nocc*(nocc+1)/2
      ippp=ipp+1
      imm=nocc*(nocc-1)/2
      do i=1,ni
        do j=1,nj
          do k=1,i
            do l=1,j
              gp=(eri1(l,k,j,i)+eri2(j,k,l,i))*0.5d0
              gm=(eri1(l,k,j,i)-eri2(j,k,l,i))*0.5d0
              call daxpy(ipp,gp,t2(1,k,i),1,ht2(1,l,j),1)
              call daxpy(imm,gm,t2(ippp,k,i),1,ht2(ippp,l,j),1)
            enddo
          enddo
        enddo
      enddo
      return
      end
