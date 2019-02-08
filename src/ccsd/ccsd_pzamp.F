      subroutine ccsd_pzamp(basis,nsh,ncor,nocc,nvir,nact,nbf,
     &                       tol2e,cmo,eorb,iprt,hiu,giu,habe,gabe,
     &                       hia,t1,z1,hz1,idiis,g_nt2,g_nz2,
     &                       g_nhz2,g_ncoul,g_nexch,tklst,use_ccsd_omp)
      implicit none
#include "errquit.fh"
C     $Id$
      integer basis,nsh,ncor,nocc,nvir,nact,nbf,
     &        iprt,idiis,g_nt2,g_nz2,g_nhz2,g_aitm,
     &        g_ncoul,g_nexch,tklst(nsh*(nsh+1)/2,2)

      double precision tol2e,cmo(nbf,nbf),eorb(nbf),
     &                 hiu(*),giu(*),habe(*),gabe(*),hia(*),
     &                 t1(nocc*nvir),hz1(nocc*nvir),z1(nocc,nvir)
c
c
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "bas.fh"
#include "ccsd_debug.fh"
c

      integer l_snsi,k_snsi,l_sisn,k_sisn,l_eri1,k_eri1,
     &        l_eri2,k_eri2,l_scre,k_scre,l_t1ao,k_t1ao,
     &        l_offsh,k_offsh,l_sc1,k_sc1,l_sc2,k_sc2,
     &        l_sc3,k_sc3,l_sc4,k_sc4,
     &        l_tz1,k_tz1,l_tz2,k_tz2
      integer maxbfsh,max1e,mem1,max2e,mem2,lscr,lssni,
     &        g_x,g_c,g_st2,g_sht2,g_zx,g_zc,g_sz2,g_shz2,
     &        g_tz1,g_tz2,g_tz3,g_qitm
      integer lnoo,lnov,lnvv
      logical stat
c
      Integer Nodes, IAm
      logical, optional, intent(in) :: use_ccsd_omp
      call errquit('ccsd_pzamp: broken ccsd_sxy interface',0,0)
      if (.not.present(use_ccsd_omp)) then
          call errquit('ccsd_pzamp: use_ccsd_omp not present!',0,0)
      endif
      Nodes = GA_NNodes()
      IAM = GA_NodeID()
c
      lnoo=nocc*nocc
      lnov=nocc*nvir
      lnvv=nvir*nvir
c
c     write(6,*)'entering pzamp'
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
     &     call errquit('ga_create g_x failed',0, GA_ERR)
c
*ga:1:0
      if (.not.ga_create(MT_DBL,nbf*nbf,nocc*nocc,'c',
     &                   nbf*nbf,0,g_c))
     &     call errquit('ga_create g_c failed',0, GA_ERR)
c
*ga:1:0
      if (.not.ga_create(MT_DBL,nbf*nbf,nocc*nocc,'zx',
     &                   nbf*nbf,0,g_zx))
     &     call errquit('ga_create g_zx failed',0, GA_ERR)
c
*ga:1:0
      if (.not.ga_create(MT_DBL,nbf*nbf,nocc*nocc,'zc',
     &                   nbf*nbf,0,g_zc))
     &     call errquit('ga_create g_zc failed',0, GA_ERR)
c -- global array for t2 with 2 ao indices
*ga:1:0
      if (.not.ga_create(MT_DBL,nocc*nocc,nbf*nbf,'st2',
     &                   nocc*nocc,0,g_st2))
     &     call errquit('ga_create g_st2 failed',0, GA_ERR)
*ga:1:0
      if (.not.ga_create(MT_DBL,nocc*nocc,nbf*nbf,'sz2',
     &                   nocc*nocc,0,g_sz2))
     &     call errquit('ga_create g_sz2 failed',0, GA_ERR)
c
c -- global array for ht2 with 2 ao indices
*ga:1:0
      if (.not.ga_create(MT_DBL,nocc*nocc,nbf*nbf,'sht2',
     &                   nocc*nocc,0,g_sht2))
     &     call errquit('ga_create g_sht2 failed',0, GA_ERR)
*ga:1:0
      if (.not.ga_create(MT_DBL,nocc*nocc,nbf*nbf,'shz2',
     &                   nocc*nocc,0,g_shz2))
     &     call errquit('ga_create g_shz2 failed',0, GA_ERR)
c
c
c ------------------------------------------------------------
c - Make t2+/- with 2 ao and 2 occ indices
c ------------------------------------------------------------
      lscr=nbf*nbf
      stat=.true.
      stat=stat.and.ma_push_get(MT_DBL,lscr,'scra',l_sc1, k_sc1)
      stat=stat.and.ma_push_get(MT_DBL,lscr,'scrb',l_sc2, k_sc2)
      if (.not.stat)call errquit('ma_get ccsd_t2pm ',0, MA_ERR)
      call ccsd_t2pm(basis,nsh,ncor,nocc,nvir,nbf,g_st2,g_nt2,cmo,
     &               t1,dbl_mb(k_sc1),dbl_mb(k_sc2))
      call ccsd_z2pm(basis,nsh,ncor,nocc,nvir,nbf,g_sz2,g_nz2,cmo,
     &               dbl_mb(k_sc1),dbl_mb(k_sc2))
      stat=stat.and.ma_pop_stack(l_sc2)
      stat=stat.and.ma_pop_stack(l_sc1)
      if (.not.stat)call errquit('ma_pop ccsd_t2pm ',0, MA_ERR)
c
c ------------------------------------------------------------
c  Form TZ intermediates
c ------------------------------------------------------------
c
      if (.not.ga_create(MT_DBL,lnov,lnov,'tz1',
     &                   nvir,nvir,g_tz1))
     &     call errquit('ga_create g_tz1 failed',0, GA_ERR)
      if (.not.ga_create(MT_DBL,lnov,lnov,'tz2',
     &                   nvir,nvir,g_tz2))
     &     call errquit('ga_create g_tz2 failed',0, GA_ERR)
      lscr=nvir*nvir
      stat=.true.
      stat=stat.and.ma_push_get(MT_DBL,nvir*nvir,'tzvv',l_tz1, k_tz1)
      stat=stat.and.ma_push_get(MT_DBL,lscr,'scra',l_sc1, k_sc1)
      stat=stat.and.ma_push_get(MT_DBL,lscr,'scrb',l_sc2, k_sc2)
      stat=stat.and.ma_push_get(MT_DBL,lscr,'scrc',l_sc3, k_sc3)
      call ccsd_tzinta(nocc,nvir,nbf,t1,dbl_mb(k_tz1),dbl_mb(k_sc1),
     &                 dbl_mb(k_sc2),dbl_mb(k_sc3),g_nt2,g_nz2,g_tz1,
     &                 g_tz2,iprt)
      stat=stat.and.ma_pop_stack(l_sc3)
      stat=stat.and.ma_pop_stack(l_sc2)
      stat=stat.and.ma_pop_stack(l_sc1)
c
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
      stat=stat.and.ma_push_get(MT_DBL,lscr,'scra',l_sc1, k_sc1)
      stat=stat.and.ma_push_get(MT_DBL,lscr,'scrb',l_sc2, k_sc2)
      stat=stat.and.ma_push_get(MT_DBL,nbf*nbf,'scrc',l_sc3, k_sc3)
      stat=stat.and.ma_push_get(MT_DBL,nbf*nbf,'scrd',l_sc4, k_sc4)
      stat=stat.and.ma_push_get(MT_DBL,nocc*nbf,'t1ao',l_t1ao, k_t1ao)
      if (.not.stat)
     &   call errquit('ccsd_pzamp: ma_push_get ccsd_sxy ',0, MA_ERR)
      call ccsd_sxy(basis,nsh,ncor,nocc,nvir,nact,nbf,g_st2,g_sht2,
     &              g_c,g_x,int_mb(k_offsh),dbl_mb(k_snsi),
     &              dbl_mb(k_sisn),lssni,dbl_mb(k_scre),
     &              mem2,max2e,dbl_mb(k_eri1),
     &              dbl_mb(k_eri2),t1,cmo,dbl_mb(k_t1ao),dbl_mb(k_sc1),
     &              dbl_mb(k_sc2),lscr,tol2e,iprt,tklst,use_ccsd_omp)
      stat=stat.and.ma_pop_stack(l_t1ao)
      stat=stat.and.ma_pop_stack(l_sc4)
      stat=stat.and.ma_pop_stack(l_sc3)
      stat=stat.and.ma_pop_stack(l_sc2)
      stat=stat.and.ma_pop_stack(l_sc1)
      stat=stat.and.ma_pop_stack(l_scre)
      stat=stat.and.ma_pop_stack(l_eri2)
      stat=stat.and.ma_pop_stack(l_eri1)
      stat=stat.and.ma_pop_stack(l_sisn)
      stat=stat.and.ma_pop_stack(l_snsi)
      if (.not.stat)
     &   call errquit('ccsd_pzamp: ma_pop_stack ccsd_sxy ',0, MA_ERR)
c
c ------------------------------------------------------------
c - Change ht2+/- with 2 ao and 2 occ indices into normal ht2
c ------------------------------------------------------------
      stat=.true.
      lscr=nbf*nbf
      stat=stat.and.ma_push_get(MT_DBL,lscr,'scra',l_sc1, k_sc1)
      stat=stat.and.ma_push_get(MT_DBL,lscr,'scrb',l_sc2, k_sc2)
      if (.not.stat)call errquit('ma_get ccsd_ht2pm ',0, MA_ERR)
      call ccsd_ht2pm(basis,nsh,ncor,nocc,nvir,nact,nbf,g_sht2,g_nhz2,
     &                cmo,
     &                dbl_mb(k_sc1),dbl_mb(k_sc2),int_mb(k_offsh))
      stat=stat.and.ma_pop_stack(l_sc2)
      stat=stat.and.ma_pop_stack(l_sc1)
      stat=stat.and.ma_pop_stack(l_offsh)
c
      call ga_zero(g_nhz2)
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
      stat=stat.and.ma_push_get(MT_DBL,lscr,'scra',l_sc1, k_sc1)
      stat=stat.and.ma_push_get(MT_DBL,lscr,'scrb',l_sc2, k_sc2)
      stat=stat.and.ma_push_get(MT_DBL,nbf*nbf,'scrc',l_sc3, k_sc3)
      stat=stat.and.ma_push_get(MT_DBL,nbf*nbf,'scrd',l_sc4, k_sc4)
      stat=stat.and.ma_push_get(MT_DBL,nocc*nbf,'t1ao',l_t1ao, k_t1ao)
      if (.not.stat)
     &   call errquit('ccsd_pzamp: ma_push_get ccsd_sxy ',0, MA_ERR)
      call ccsd_sxy(basis,nsh,ncor,nocc,nvir,nact,nbf,g_sz2,g_shz2,
     &              g_zc,g_zx,int_mb(k_offsh),dbl_mb(k_snsi),
     &              dbl_mb(k_sisn),lssni,dbl_mb(k_scre),
     &              mem2,max2e,dbl_mb(k_eri1),
     &              dbl_mb(k_eri2),z1,cmo,dbl_mb(k_t1ao),dbl_mb(k_sc1),
     &              dbl_mb(k_sc2),lscr,tol2e,iprt,tklst,use_ccsd_omp)
      stat=stat.and.ma_pop_stack(l_t1ao)
      stat=stat.and.ma_push_get(MT_DBL,nvir*nbf,'t1ao',l_t1ao, k_t1ao)
      call ccsd_aosig(basis,nsh,nocc,nvir,nbf,cmo,dbl_mb(k_scre),
     &               dbl_mb(k_eri1),dbl_mb(k_eri2),t1,hz1,
     &               dbl_mb(k_tz1),dbl_mb(k_sc1),dbl_mb(k_sc2),
     &               dbl_mb(k_t1ao),dbl_mb(k_sc3),dbl_mb(k_sc4),mem2,
     &               max2e,tol2e,g_nt2,g_nz2,g_nhz2,g_tz1,g_tz2)
      stat=stat.and.ma_pop_stack(l_t1ao)
      stat=stat.and.ma_pop_stack(l_sc4)
      stat=stat.and.ma_pop_stack(l_sc3)
      stat=stat.and.ma_pop_stack(l_sc2)
      stat=stat.and.ma_pop_stack(l_sc1)
      stat=stat.and.ma_pop_stack(l_scre)
      stat=stat.and.ma_pop_stack(l_eri2)
      stat=stat.and.ma_pop_stack(l_eri1)
      stat=stat.and.ma_pop_stack(l_sisn)
      stat=stat.and.ma_pop_stack(l_snsi)
      if (.not.stat) 
     &   call errquit('ccsd_pzamp: ma_pop_stack ccsd_sxy ',0, MA_ERR)
      stat=.true.
      lscr=nbf*nbf
      stat=stat.and.ma_push_get(MT_DBL,lscr,'scra',l_sc1, k_sc1)
      stat=stat.and.ma_push_get(MT_DBL,lscr,'scrb',l_sc2, k_sc2)
      call ccsd_hz2pm(basis,nsh,ncor,nocc,nvir,nbf,g_shz2,g_nhz2,
     &                cmo,
     &                dbl_mb(k_sc1),dbl_mb(k_sc2),int_mb(k_offsh))
      stat=stat.and.ma_pop_stack(l_sc2)
      stat=stat.and.ma_pop_stack(l_sc1)
      stat=stat.and.ma_pop_stack(l_offsh)
c
c ------------------------------------------------------------
c Z vector routines
c ------------------------------------------------------------
c
      lscr=nbf*nbf
      stat=.true.
      stat=stat.and.ma_push_get(MT_DBL,lnoo,'tzoo',l_tz2, k_tz2)
      stat=stat.and.ma_push_get(MT_DBL,lscr,'scra',l_sc1, k_sc1)
      stat=stat.and.ma_push_get(MT_DBL,lscr,'scrb',l_sc2, k_sc2)
      stat=stat.and.ma_push_get(MT_DBL,lscr,'scre',l_sc3, k_sc3)
      stat=stat.and.ma_push_get(MT_DBL,lscr,'scrf',l_sc4, k_sc4)
      call ccsd_zsig1(nocc,nvir,nbf,hiu,giu,habe,gabe,hia,t1,hz1,
     &                z1,dbl_mb(k_tz1),dbl_mb(k_tz2),
     &                dbl_mb(k_sc1),dbl_mb(k_sc2),dbl_mb(k_sc1),
     &                dbl_mb(k_sc2),dbl_mb(k_sc3),dbl_mb(k_sc4),
     &                g_nt2,g_nz2,g_ncoul,g_nexch,g_x,g_c,
     &                g_sht2,g_zx,g_zc,g_shz2,g_aitm,g_qitm,
     &                g_tz1,g_tz2,g_tz3,iprt)
      call ccsd_zsig2(nocc,nvir,nbf,hiu,giu,habe,gabe,hia,t1,
     &                z1,dbl_mb(k_tz1),dbl_mb(k_tz2),
     &                dbl_mb(k_sc1),dbl_mb(k_sc2),dbl_mb(k_sc1),
     &                dbl_mb(k_sc2),dbl_mb(k_sc3),dbl_mb(k_sc4),
     &                g_nt2,g_nz2,g_nhz2,g_ncoul,g_nexch,g_x,g_c,
     &                g_sht2,g_zx,g_zc,g_aitm,g_qitm,g_tz3,iprt)
      stat=.true.
      stat=stat.and.ma_pop_stack(l_sc4)
      stat=stat.and.ma_pop_stack(l_sc3)
      stat=stat.and.ma_pop_stack(l_sc2)
      stat=stat.and.ma_pop_stack(l_sc1)
      stat=stat.and.ma_pop_stack(l_tz2)
      stat=stat.and.ma_pop_stack(l_tz1)
c
c ------------------------------------------------------------
c remove global arrays associated with B intermediate
c ------------------------------------------------------------
      if (.not.ga_destroy(g_shz2))
     &    call errquit('ga_dest g_shz2 fail',0, GA_ERR)
      if (.not.ga_destroy(g_sht2))
     &    call errquit('ga_dest g_sht2 fail',0, GA_ERR)
      if (.not.ga_destroy(g_sz2))
     &    call errquit('ga_dest g_sz2 fail',0, GA_ERR)
      if (.not.ga_destroy(g_st2))
     &    call errquit('ga_dest g_st2 fail',0, GA_ERR)
c ------------------------------------------------------------
c remove global arrays associated with vvvo integrals (g_x and g_c)
c ------------------------------------------------------------
      if (.not.ga_destroy(g_zc))
     &    call errquit('ga_dest g_zc fail',0, GA_ERR)
      if (.not.ga_destroy(g_zx))
     &    call errquit('ga_dest g_zx fail',0, GA_ERR)
      if (.not.ga_destroy(g_c))
     &    call errquit('ga_dest g_c fail',0, GA_ERR)
      if (.not.ga_destroy(g_x))
     &    call errquit('ga_dest g_x fail',0, GA_ERR)
c
      return
      end
