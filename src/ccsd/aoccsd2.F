      subroutine aoccsd(basis,ncor,nocc,nvir,ndel,nact,nbf,maxit,
     &                  convi,iprt,cmo,eorb,blen,
     &                  g_ncoul, g_nexch, RefEner,
     $   CC_Theory, rtdb, mxvec, geom, Tol2e, occd, oconverged)
      implicit none
C     $Id: aoccsd2.F,v 2.25 1999-06-16 00:10:12 d3e129 Exp $
c
      integer basis,ncor,nocc,nvir,ndel,nact,nbf,maxit,iprt,
     &        blen,g_objv,g_coul,g_exch,g_ncoul,g_nexch, mxvec
      double precision cmo(nbf,nbf),eorb(nbf)
      Double Precision RefEner, Convi, Tol2E
      Character*(*) CC_Theory
      Integer RTDB, geom
      logical oconverged, occd
c
#include "ccsd_len.fh"
#include "tcgmsg.fh"
#include "bas.fh"
#include "rtdb.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "stdio.fh"
#include "msgids.fh"
#include "util.fh"
c
c ccsd
      Integer i
      logical stat
      integer nsh,maxbfsh,max2e,mem2
      integer l_t1,l_ht1,l_scra,l_scrb,
     &        l_hiu,l_hia,l_habe,l_giu,l_gabe,
     &        l_bbkp,l_bb
      integer k_t1,k_ht1,k_scra,k_scrb,
     &        k_hiu,k_hia,k_habe,k_giu,k_gabe,
     &        k_bbkp,k_bb
      integer g_t2,g_ht2
      integer ad1,ad2,ad3,ad4
c triples
      integer lbfo,lbfv,mem_avail
      integer l_f1n,l_f1t,l_f2n,l_f2t,l_f3n,l_f3t,l_f4n,l_f4t,
     &        l_buf1,l_buf2,l_tk
      integer k_f1n,k_f1t,k_f2n,k_f2t,k_f3n,k_f3t,k_f4n,k_f4t,
     &        k_buf1,k_buf2,k_tk
      integer g_objo
      integer vchunk,noper,occ_lo,vir_lo,vir_hi,oseg_lo,oseg_hi,
     &        nvpass
      integer ga_create_jkblocked
      logical oexch,ocoul,ohalf
      double precision eccsdt,eccsd,empt(2),tx(3)
      Integer Nodes, IAm, memdrv, memtrn
c
      integer kchunk,nkpass,    ! For reduced-memory triples
     $        k_trp_Tij, k_trp_Tkj, k_trp_Tia, 
     $        k_trp_Tka, k_trp_Xia, k_trp_Xka, 
     $        k_trp_Jia, k_trp_Jka, k_trp_Kia, 
     $        k_trp_Kka, k_trp_Jij, k_trp_Jkj, 
     $        k_trp_Kij, k_trp_Kkj, k_trp_Dja,
     $        k_trp_Djka, k_trp_Djia,
     $        l_trp_Tij, l_trp_Tkj, l_trp_Tia, 
     $        l_trp_Tka, l_trp_Xia, l_trp_Xka, 
     $        l_trp_Jia, l_trp_Jka, l_trp_Kia, 
     $        l_trp_Kka, l_trp_Jij, l_trp_Jkj, 
     $        l_trp_Kij, l_trp_Kkj, l_trp_Dja,
     $        l_trp_Djka, l_trp_Djia
c
      
      integer moints_trp_mem
      external moints_trp_mem
C
      double precision zip
*     double precision trace1, trace2, trace3, trace4
      integer lhoo,lhoov
C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
C     Temporary stuff added to play density matrix games
C
      Integer A, J, B, LSab
      Integer L_Sb, K_Sb, L_Sa, K_Sa
      Integer H_DenZ1, I_DenZ1, H_DenT1, I_DenT1
      Integer G_DenT2, G_DenZ2, G_1PDM
      Logical Do_CCDen
      Integer GA_Create_Atom_Blocked
      External GA_Create_Atom_Blocked
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
c gradient
      integer l_z1
      integer k_z1
      integer g_nt2
      integer ilo,ihi,jlo,jhi,offt2,nvp1
      Logical Do_Z
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
* DATA MUST FOLLOW ALL DELCARATIONS 
      Data Do_Z /.FALSE./
      Data Do_CCDen /.FALSE./
C
      data oexch,ocoul/.true.,.true./
      data ohalf/.false./
      data zip/0.0d00/
c$$$      double precision rtc
c$$$      external rtc
c
      call ga_sync()
      call qenter('aoccsd',0)
c
      nodes = ga_nnodes()
      iam = ga_nodeid()
      stat=bas_numcont(basis,nsh)
      stat=stat.and.bas_nbf_cn_max(basis,maxbfsh)
      if (.not.stat) call errquit('moints: cannot get basis info',0)
c
c      print *,'basis set handle ',basis
c      print *,' ga allocated start of aoccsd '
c      call ga_summarize(.true.)
c
      if (iam.eq.0)then
      write(6,1001)nodes
 1001 format(//,1x,76(1h*),/,14x,
     &       'the segmented parallel titan ccsd program:',
     &       i3,' nodes',/,1x,76(1h*),//)
      endif
      if (iam.eq.0)then
      write (*,1318) CC_Theory( :Inp_StrLen(CC_Theory) ),
     $      ncor,nocc,nvir,ndel,nbf,nsh,basis
 1318 format(//,' level of theory    ',A,/
     $          ' number of core     ',i5,/,
     &          ' number of occupied ',i5,/,
     &          ' number of virtual  ',i5,/,
     &          ' number of deleted  ',i5,/,
     &          ' total functions    ',i5,/,
     &          ' number of shells   ',i5,/,
     &          ' basis label        ',i5,/)
      write(*,6003) iprt,convi,maxit, mxvec
 6003 format(//,2x,' ***** ccsd parameters *****',
     &        /,2x,' iprt   = ',i5,
     &        /,2x,' convi  = ',e9.3,
     &        /,2x,' maxit  = ',i5,
     $        /,2x,' mxvec  = ',i5)
      endif
C
C Initialize GA_TRACING
C
CRK   call trace_init(1000)

c
c----------------------------------------------------------------
c Begin parallel ccsd code
c----------------------------------------------------------------
c
      mem_avail = ma_inquire_avail(MT_DBL)
      if(iam.eq.0)write(6,*)'memory',mem_avail
      call int_mem_2e4c(max2e, mem2)
      mem2=max(mem2,nbf*nbf)
      if (iprt.gt.10.and.iam.eq.0)then
         write(6,*)'lnoo,lnov,lnvv,lnoov,lnooo,lnovv',
     &              lnoo,lnov,lnvv,lnoov,lnooo,lnovv
      endif
      lnobj=max(lnoov,nbf*nbf)
      stat=.true.
      stat=stat.and.ma_push_get(MT_DBL,lnov,'t1',l_t1,k_t1)
      stat=stat.and.ma_push_get(MT_DBL,lnov,'ht1',l_ht1,k_ht1)
      stat=stat.and.ma_push_get(MT_DBL,lnoo,'hiu',l_hiu,k_hiu)
      stat=stat.and.ma_push_get(MT_DBL,lnov,'hia',l_hia,k_hia)
      stat=stat.and.ma_push_get(MT_DBL,lnvv,'habe',l_habe,k_habe)
      stat=stat.and.ma_push_get(MT_DBL,lnoo,'giu',l_giu,k_giu)
      stat=stat.and.ma_push_get(MT_DBL,lnvv,'gabe',l_gabe,k_gabe)
      stat=stat.and.ma_push_get(MT_DBL,mxvec*(mxvec+1)/2,'bbkp',
     &                          l_bbkp,k_bbkp)
      stat=stat.and.ma_push_get(MT_DBL,(mxvec+1)*(mxvec+2),'bb',
     &                          l_bb,k_bb)
      stat=stat.and.ma_push_get(MT_INT,nsh*(nsh+1),'tklst',l_tk,k_tk)
      stat=stat.and.ma_push_get(MT_DBL,lnobj,'1scra',l_scra,k_scra)
      stat=stat.and.ma_push_get(MT_DBL,lnobj,'scrb',l_scrb,k_scrb)
      if (.not.stat)call errquit('ccsd: failed to allocate ',0)
c
        lhoo=nocc*(nocc+1)/2
        lhoov=nvir*nocc*(nocc+1)/2
*ga:1:0
      if (.not.ga_create(MT_DBL,mxvec*lhoov,nvir+1,'g_t2', 
     &                   mxvec*lhoov,0,g_t2))
     &     call errquit('ga_create failed for g_t2',0)
*ga:1:0
      if (.not.ga_create(MT_DBL,mxvec*lhoov,nvir+1,'g_ht2', 
     &                   mxvec*lhoov,0,g_ht2))
     &     call errquit('ga_create failed for g_ht2',0)
c      
      call ccsd_iterdrv2(rtdb,basis,nsh,ncor,nocc,nvir,nact,nbf,
     &                  maxit,convi,iprt,cmo,eorb,
     &                  dbl_mb(k_t1),dbl_mb(k_ht1),
     &                  dbl_mb(k_scra),dbl_mb(k_scrb),
     &                  dbl_mb(k_hiu),dbl_mb(k_hia),dbl_mb(k_habe),
     &                  dbl_mb(k_giu),dbl_mb(k_gabe),
     &                  dbl_mb(k_bbkp),dbl_mb(k_bb),
     &                  g_t2,g_ht2,
     &                  mxvec,eccsd,max2e,mem2,g_ncoul,g_nexch,
     &                  int_mb(k_tk),offt2,nvp1, refener, Tol2E, 
     $     occd, oconverged)
c     
      call ga_sync
c
         if (.not.ga_create(MT_DBL,lnov,lnov,'g_nt2',
     &      nvir,nvir,g_nt2))
     &      call errquit('ga_create failed g_nt2',0)
c
      call ga_get(g_t2,offt2+1,offt2+lnov,nvp1,nvp1,
     $            dbl_mb(k_t1),lnov)
c
      call ga_distribution(g_t2,iam,jlo,jhi,ilo,ihi)
      do a=1,nvir
        if (a.ge.ilo.and.a.le.ihi)then
          call ga_get(g_t2,offt2+1,offt2+lhoov,a,a,
     &                 dbl_mb(k_scra),lhoov)
          do i=1,nocc
            call dfill(lnov,0.0d00,dbl_mb(k_scrb),1)
            do j=1,i
              do b=1,nvir
                ad1=(b-1)*lhoo+i*(i-1)/2+j
                ad2=(j-1)*nvir+b
                dbl_mb(k_scrb+ad2-1)=dbl_mb(k_scra+ad1-1)
              enddo
            enddo
            ad1=(i-1)*nvir+a
            call ga_put(g_nt2,1,i*nvir,ad1,ad1,dbl_mb(k_scrb),i*nvir)
          enddo
        endif
      enddo
      call ga_sync
c
      call ga_distribution(g_nt2,iam,jlo,jhi,ilo,ihi)
      do i=1,nocc
        ad1=(i-1)*nvir
        if (ad1+1.ge.ilo.and.ad1+1.le.ihi)then
          do j=1,i-1
            ad2=(j-1)*nvir
            if (ad2+1.ge.jlo.and.ad2+1.le.jhi)then
              call ga_get(g_nt2,ad2+1,ad2+nvir,ad1+1,ad1+nvir,
     &                    dbl_mb(k_scra),nvir)
              do a=1,nvir
              do b=1,nvir
                ad3=(b-1)*nvir+a
                ad4=(a-1)*nvir+b
                dbl_mb(k_scrb+ad4-1)=dbl_mb(k_scra+ad3-1)
              enddo
              enddo
              call ga_put(g_nt2,ad1+1,ad1+nvir,ad2+1,ad2+nvir,
     &                    dbl_mb(k_scrb),nvir)
            endif
          enddo
        endif
      enddo
      call ga_sync
c
         if (iam.eq.0.and.iprt.gt.5)then
            print *,'starting t1 vector'
            write(6,139)(dbl_mb(k_t1-1+i),i=1,lnov)
 139        format(4e14.4)
            call ga_print(g_nt2)
         endif
c
      stat=stat.and.ma_pop_stack(l_scrb)
      stat=stat.and.ma_pop_stack(l_scra)
c
      lnobj=2*lnoov+max(lnoov,lnooo)
*ga:1:0
      if (.not.ga_create(MT_DBL,lnobj,nvir+1,'objv',0,0,
     &                   g_objv))
     &     call errquit('ga_create g_objv failed',0)
      call ga_zero(g_objv)
c
      lnobj=lnoov+max(lnooo,lnoov)
      stat=stat.and.ma_push_get(MT_DBL,lnobj,'1scra',l_scra,k_scra)
      stat=stat.and.ma_push_get(MT_DBL,nocc*nact,'scrb',l_scrb,k_scrb)
      call ccsd_trpmo(ncor,nocc,nvir,nact,nbf,dbl_mb(k_scra),
     $            dbl_mb(k_scrb),dbl_mb(k_t1),g_nt2,g_objv,g_ncoul,
     $            g_nexch,iprt)
      stat=stat.and.ma_pop_stack(l_scrb)
      stat=stat.and.ma_pop_stack(l_scra)
c
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
CRK
CRK   Gradient driver to go in here for development purposes 
CRK 
c
      If ( Do_Z) then
        lnobj=max(lnoov,nbf*nbf)
        stat=stat.and.ma_push_get(MT_DBL,lnobj,'3scra',l_scra,k_scra)
        stat=stat.and.ma_push_get(MT_DBL,lnobj,'scrb',l_scrb,k_scrb)
         stat=stat.and.ma_push_get(MT_DBL,lnov,'z1',l_z1,k_z1)
         call ccsd_graddrv(basis,nsh,ncor,nocc,nvir,nact,nbf,
     &      maxit,convi,iprt,cmo,eorb,dbl_mb(k_t1),dbl_mb(k_ht1),
     &      dbl_mb(k_scra),dbl_mb(k_scrb),dbl_mb(k_z1),
     &      dbl_mb(k_hiu),dbl_mb(k_hia),dbl_mb(k_habe),
     &      dbl_mb(k_giu),dbl_mb(k_gabe),
     &      dbl_mb(k_bbkp),dbl_mb(k_bb),
     &      g_nt2,g_t2,g_ht2,
     &      mxvec,eccsd,max2e,mem2,g_ncoul,g_nexch,
     &      int_mb(k_tk), Tol2E)
      EndIf
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
c
      stat=.true.
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      If ( Do_Z ) then
         stat=stat.and.ma_pop_stack(l_z1)
         stat=stat.and.ma_pop_stack(l_scrb)
         stat=stat.and.ma_pop_stack(l_scra)
      EndIf
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      stat=stat.and.ma_pop_stack(l_tk)
      stat=stat.and.ma_pop_stack(l_bb)
      stat=stat.and.ma_pop_stack(l_bbkp)
      stat=stat.and.ma_pop_stack(l_gabe)
      stat=stat.and.ma_pop_stack(l_giu)
      stat=stat.and.ma_pop_stack(l_habe)
      stat=stat.and.ma_pop_stack(l_hia)
      stat=stat.and.ma_pop_stack(l_hiu)
      stat=stat.and.ma_pop_stack(l_ht1)
      stat=stat.and.ma_pop_stack(l_t1)
      if (.not.stat)call errquit('ccsd: cannot pop stack?',0)
      call ga_sync()
      if (.not. ga_destroy(g_nt2))
     &     call errquit('ga_destroy failed',0)
      if (.not. ga_destroy(g_ht2))
     &     call errquit('ga_destroy failed',0)
      if (.not. ga_destroy(g_t2))
     &     call errquit('ga_destroy failed',0)
c
      call qexit('aoccsd',0)
c
      If ( GA_NodeID() .eq. 0 .AND.
     $   ( Util_Print('byproduct energies', Print_Default)
     $   .OR. Inp_Compare(.FALSE., 'ccsd', CC_Theory) ) ) then
         Call Banner(LuOut,
     $      'CCSD Energy', '-',
     $      .TRUE., .TRUE., .FALSE.)
         Write (LuOut, '(1X, A, F30.15)') 'Reference energy:  ',
     $      RefEner
         Write (LuOut, '(1X, A, F30.15)') 'CCSD corr. energy: ',
     $      eccsd
         Write (LuOut, '(1X, A, F30.15)') 'Total CCSD energy: ',
     $      RefEner + eccsd
         call ecce_print1('ccsd correlation energy', mt_dbl, eccsd, 1)
         call ecce_print1('ccsd total energy', mt_dbl, refener+eccsd, 1)
         call ecce_print1('total ccsd energy', mt_dbl, refener+eccsd, 1)
      EndIf
C
C     Write our energy out to the database
C
      If ( .NOT. RTDB_Put( RTDB, 'ccsd:ccsd correlation energy',
     $   MT_Dbl, 1, Eccsd) )
     $   Call ErrQuit('CCSD: Can''t write to RTDB', 1)
      If ( .NOT. RTDB_Put( RTDB, 'ccsd:ccsd energy', MT_Dbl, 1,
     $   RefEner+eccsd) ) 
     $   Call ErrQuit('CCSD: Can''t write to RTDB', 2)
C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
C     Temporarily a good place to mess around with some density
C     matrix stuff
C
      If (Do_CCDen) then
C
C     Start by creating the data structures we need
C
         If ( .NOT. MA_Push_Get(MT_Dbl, NOcc*NVir, 'CCDens T1',
     $      H_DenT1, I_DenT1) ) Call ErrQuit(
     $      'CCDen (prep): unable to allocate T1', NOcc*NVir)
C
         If ( .NOT. MA_Push_Get(MT_Dbl, NOcc*NVir, 'CCDens Z1',
     $      H_DenZ1, I_DenZ1) ) Call ErrQuit(
     $      'CCDen (prep): unable to allocate Z1', NOcc*NVir)
C
         If ( .NOT. GA_Create(MT_Dbl, NVir*NOcc, NVir*NOcc,
     $      'CCDen T2', NVir, NVir, G_DenT2) ) Call ErrQuit(
     $      'CCDen (prep): unable to allocate T2', NVir*NOcc*NVir*NOcc)
C
         If ( .NOT. GA_Create(MT_Dbl, NVir*NOcc, NVir*NOcc,
     $      'CCDen Z2', NVir, NVir, G_DenZ2) ) Call ErrQuit(
     $      'CCDen (prep): unable to allocate Z2', NVir*NOcc*NVir*NOcc)
C
         G_1PDM = GA_Create_Atom_Blocked(Geom, Basis, 'CCDen 1-PDM')
C
C        Now rearrange the existing data into our new structures
C
C        First, T1 (and for now fill Z1 with T1 too)
C        (T1 is replicated for density, so all nodes do this)
C
         LnObj = lnoov + lnooo
         Call GA_Get(g_ObjV, LnObj+1, LnObj+Lnov, NVir+1, NVir+1,
     &      Dbl_MB(I_DenT1), NOcc*NVir)
         Call GA_Get(g_ObjV, LnObj+1, LnObj+Lnov, NVir+1, NVir+1,
     &      Dbl_MB(I_DenZ1), NOcc*NVir)
C
C        Now T2 (then copy to Z2 for the time being)
C        (T2 is a GA for density, so have node 0 do it for simplicity)
C
         If ( GA_NodeID() .eq. 0) then
            lsab=max(lnoov,nbf*nbf)
            If ( .NOT. ma_push_get(MT_DBL,lsab,'sa',l_sa, k_sa) )
     $         Call ErrQuit('CCDen (prep): can''t allocate sa', lsab)
            If (.NOT. ma_push_get(MT_DBL,lsab,'sb',l_sb, k_sb) )
     $         Call ErrQuit('CCDen (prep): can''t allocate sb', lsab)
            do a=1,nvir
               call ga_get(g_ObjV,lnobj+1, lnobj+lnoov,a,a,
     &            dbl_mb(k_sa),lnoov)
               do i=1,nocc
                  do j=1,nocc
                     do b=1,nvir
                        ad1=k_sa+(b-1)*lnoo+(i-1)*nocc+j-1
                        ad2=k_sb+(j-1)*nvir+b-1
                        dbl_mb(ad2)=dbl_mb(ad1)
                     enddo
                  enddo
                  ad1=(i-1)*nvir+a
                  call ga_put(g_DenT2,1,lnov,ad1,ad1,dbl_mb(k_sb),lnov)
               enddo
            enddo
            If ( .NOT. MA_Pop_Stack(L_Sb)) Call ErrQuit(
     $         'CCDen (prep): freeing sb', 0)
            If ( .NOT. MA_Pop_Stack(L_Sa)) Call ErrQuit(
     $         'CCDen (prep): freeing sa', 0)
         EndIf
         Call GA_Sync
         Call GA_Copy(G_DenT2, G_DenZ2)
C
C        Now call the Density driver routine
C
         Call CCDen_Driver(NOcc, NVir, Dbl_MB(I_DenT1), NOcc, G_DenT2,
     $      Dbl_MB(I_DenZ1), NOcc, G_DenZ2, G_1PDM)
C
         If ( .NOT. GA_Destroy(G_1PDM) ) Call ErrQuit(
     $      'CCDen (prep): unable to free 1PDM', 0)
         If ( .NOT. GA_Destroy(G_DenZ2) ) Call ErrQuit(
     $      'CCDen (prep): unable to free Z2', 0)
         If ( .NOT. GA_Destroy(G_DenT2) ) Call ErrQuit(
     $      'CCDen (prep): unable to free T2', 0)
         If ( .NOT. MA_Pop_Stack(H_DenZ1) ) Call ErrQuit(
     $      'CCDen (prep): unable to free Z1', 0)
         If ( .NOT. MA_Pop_Stack(H_DenT1) ) Call ErrQuit(
     $      'CCDen (prep): unable to free T1', 0)
      EndIf
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
c
      if (.not. oconverged) goto 999
c
C     Do we continue to the triples or not?
C     Note: task_energy expects ccsd:energy as a marker of successful
C     completion.  Therefore, we don't want to write it unless we
C     are certain we have completed the expected task.
C
      If ( .NOT. ( Inp_Compare(.FALSE., 'ccsd(t)', CC_Theory)
     $   .OR. Inp_Compare(.FALSE., 'ccsd+t(ccsd)', CC_Theory) ) ) then
         call ecce_print1('total energy', mt_dbl, refener+eccsd, 1)
         If ( .NOT. RTDB_Put( RTDB, 'ccsd:energy', MT_Dbl, 1,
     $      RefEner+eccsd) ) 
     $      Call ErrQuit('CCSD: Can''t write to RTDB', 3)
         Goto 999
      EndIf
c
c ----------------------------
c the new parallel triples
c ----------------------------
c
      mem_avail = ma_inquire_avail(MT_DBL)
      if(iam.eq.0)write(6,*)'memory',mem_avail
      call qenter('triples',0)
      if (iam.eq.0)write(6,293)
 293  format(//,'*********triples calculation*********',/)
      lbfo=lnovv+lnoov+lnoov
      lbfv=lnoov+lnoov+lnoov
      lbfv=max(lbfv,nbf*nbf)
      lnobj=max(2*lnoov+lnooo,3*lnoov)
      stat=.true.
      stat=stat.and.ma_push_get(MT_DBL,lnov,'t1',l_t1,k_t1)
      stat=stat.and.ma_push_get(MT_DBL,lnobj,'buf1',l_buf1,k_buf1)
      stat=stat.and.ma_push_get(MT_DBL,lnobj,'buf2',l_buf2,k_buf2)
c
*ga:1:0
*      if (.not.ga_create(MT_DBL,lbfo,nocc,'objo',lbfo,0,g_objo))
*     &     call errquit('ga_create failed for g_objo',0)
*ga:1:0
      if (.not.ga_create(MT_DBL,lbfo,nocc,'objo',0,0,g_objo))
     &     call errquit('ga_create failed for g_objo',0)
c
      call ga_zero(g_objo)
      call ccsd_rdtrpo
     $     (dbl_mb(k_t1),dbl_mb(k_buf1),dbl_mb(k_buf2),
     $                 g_objo,g_objv,nocc,nvir,iprt)
      stat=.true.
      stat=stat.and.ma_pop_stack(l_buf2)
      stat=stat.and.ma_pop_stack(l_buf1)
      if (.not.stat)call errquit('rdtrpo: cannot pop stack?',0)
c
c      print *,' ga allocated before trpdrv '
c      call ga_summarize(.true.)
c
c     Iteratively determine nvpass (vchunk) and 
c     internal chunking for the tripes (kchunk)
c
      mem_avail = ma_inquire_avail(MT_DBL)
      call ga_igop(msg_cc_mem, mem_avail, 1, 'min')
c
      nvpass = 1
c
 7651 nkpass = 1
 7652 kchunk = (nocc - 1)/nkpass + 1
      vchunk = (nvir - 1)/nvpass + 1
c
      memdrv = lnvv*(11+3*kchunk) + lnov*(5+4*kchunk) + nvir*(1+kchunk)
      memtrn = moints_trp_mem(basis, vchunk, blen)
c
      if (iam .eq. 0) then
         write(6,1101) nkpass, nvpass,memdrv,memtrn,mem_avail
 1101    format('nkpass=',i2,'; nvpass=',i3,'; memdrv=',i8,
     $        '; memtrn=',i8,'; memavail=',i8)
         call util_flush(6)
      endif
c
      noper = nodes * 
     $     ((mem_avail - max(memdrv, memtrn))/(nbf*nbf + nbf*nvir))
      if (vchunk .gt. max(min(nvir,noper/nocc),1)) then
c
c     If the triples part is using more memory than the 4-index
c     then increase nkpass, otherwise increase nvpass
c
         if (memdrv .gt. memtrn) then
            nkpass = nkpass + 1
            if (nkpass .gt. nocc) call errquit
     $           ('ccsd triples: too many k passes', nkpass)
            goto 7652
         else
            nvpass = nvpass + 1
            if (nvpass .gt. nvir) call errquit
     $           ('ccsd triples: too many v passes ', nvpass)
            goto 7651
         endif
      endif
c
      noper=nocc*vchunk
c
      if (iam.eq.0)write(*,845)mem_avail,nvir,vchunk,nvpass,
     $     kchunk,nkpass
  845 format(' memory available/node                ',i8,/,
     &       ' total number of virtual orbitals     ',i5,/,
     &       ' number of virtuals per integral pass ',i5,/,
     &       ' number of integral evaluations       ',i5,/,
     $       ' number of occupied per triples pass  ',i5,/,
     $       ' number of triples passes             ',i5/)
      g_coul = ga_create_JKblocked(noper,nbf,nbf,'Coulomb oper')
      g_exch = ga_create_JKblocked(noper,nbf,nvir,'X oper')
      nvpass=0
      occ_lo = ncor+1
      vir_lo = ncor+nocc+1
      vir_hi = ncor+nocc+nvir
      oseg_hi = vir_lo - 1
      empt(1)=zip
      empt(2)=zip
      do while (oseg_hi.lt.vir_hi)
         nvpass = nvpass + 1
         oseg_lo = oseg_hi + 1
         oseg_hi = min(oseg_hi+vchunk,vir_hi)
         if (iam.eq.0)write(*,846)nvpass
 846     format(' commencing integral evaluation ',i8)
         tx(1)=tcgtime()
         call moints_trp(basis,  ohalf, occ_lo, oseg_lo, oseg_hi,
     &      vir_lo, vir_hi,
     &      g_coul, ocoul, g_exch, oexch,
     &      blen, cmo, ncor, nocc, nvir, Tol2E)
c        print *,'call trpdrv ',nvpass

*          if (iam .eq. 0) then
*            write(6,*) ' Out of MOINTS '
*            call util_flush(6)
*          endif

         if (.not. ma_push_get(MT_DBL,lnvv,'f1n',l_f1n,k_f1n))
     $        call errquit('aoccsd2:triples f1n',lnvv)
         if (.not. ma_push_get(MT_DBL,lnvv,'f1t',l_f1t,k_f1t))
     $        call errquit('aoccsd2:triples f1t',lnvv)
         if (.not. ma_push_get(MT_DBL,lnvv,'f2n',l_f2n,k_f2n))
     $        call errquit('aoccsd2:triples f2n',lnvv)
         if (.not. ma_push_get(MT_DBL,lnvv,'f2t',l_f2t,k_f2t))
     $        call errquit('aoccsd2:triples f2t',lnvv)
         if (.not. ma_push_get(MT_DBL,lnvv,'f3n',l_f3n,k_f3n))
     $        call errquit('aoccsd2:triples f3n',lnvv)
         if (.not. ma_push_get(MT_DBL,lnvv,'f3t',l_f3t,k_f3t))
     $        call errquit('aoccsd2:triples f3t',lnvv)
         if (.not. ma_push_get(MT_DBL,lnvv,'f4n',l_f4n,k_f4n))
     $        call errquit('aoccsd2:triples f4n',lnvv)
         if (.not. ma_push_get(MT_DBL,lnvv,'f4t',l_f4t,k_f4t))
     $        call errquit('aoccsd2:triples f4t',lnvv)
c
      if (.not. ma_push_get(MT_DBL, lnvv,'k_trp_Tij',
     $        l_trp_Tij,k_trp_Tij))
     $     call errquit('aoccsd2: triples: k_trp_Tij ', lnvv)
      if (.not. ma_push_get(MT_DBL, kchunk*lnvv,'k_trp_Tkj',
     $     l_trp_Tkj,k_trp_Tkj))
     $     call errquit('aoccsd2: triples: k_trp_Tkj ', kchunk*lnvv)
      if (.not. ma_push_get(MT_DBL, lnov,'k_trp_Tia',
     $     l_trp_Tia,k_trp_Tia))
     $     call errquit('aoccsd2: triples: k_trp_Tia ', lnov)
      if (.not. ma_push_get(MT_DBL, kchunk*lnov,'k_trp_Tka',
     $     l_trp_Tka,k_trp_Tka))
     $     call errquit('aoccsd2: triples: k_trp_Tka ', kchunk*lnov)
      if (.not. ma_push_get(MT_DBL, lnov,'k_trp_Xia',
     $     l_trp_Xia,k_trp_Xia))
     $     call errquit('aoccsd2: triples: k_trp_Xia ', lnov)
      if (.not. ma_push_get(MT_DBL, kchunk*lnov,'k_trp_Xka',
     $     l_trp_Xka,k_trp_Xka))
     $     call errquit('aoccsd2: triples: k_trp_Xka ', kchunk*lnov)
      if (.not. ma_push_get(MT_DBL, lnvv,'k_trp_Jia',
     $     l_trp_Jia,k_trp_Jia))
     $     call errquit('aoccsd2: triples: k_trp_Jia ', lnvv)
      if (.not. ma_push_get(MT_DBL, kchunk*lnvv,'k_trp_Jka',
     $     l_trp_Jka,k_trp_Jka))
     $     call errquit('aoccsd2: triples: k_trp_Jka ', kchunk*lnvv)
      if (.not. ma_push_get(MT_DBL, lnvv,'k_trp_Kia',
     $     l_trp_Kia,k_trp_Kia))
     $     call errquit('aoccsd2: triples: k_trp_Kia ', lnvv)
      if (.not. ma_push_get(MT_DBL, kchunk*lnvv,'k_trp_Kka',
     $     l_trp_Kka,k_trp_Kka))
     $     call errquit('aoccsd2: triples: k_trp_Kka ', kchunk*lnvv)
      if (.not. ma_push_get(MT_DBL, lnov,'k_trp_Jij',
     $     l_trp_Jij,k_trp_Jij))
     $     call errquit('aoccsd2: triples: k_trp_Jij ', lnov)
      if (.not. ma_push_get(MT_DBL, kchunk*lnov,'k_trp_Jkj',
     $     l_trp_Jkj,k_trp_Jkj))
     $     call errquit('aoccsd2: triples: k_trp_Jkj ', kchunk*lnov)
      if (.not. ma_push_get(MT_DBL, lnov,'k_trp_Kij',
     $     l_trp_Kij,k_trp_Kij))
     $     call errquit('aoccsd2: triples: k_trp_Kij ', lnov)
      if (.not. ma_push_get(MT_DBL, kchunk*lnov,'k_trp_Kkj',
     $     l_trp_Kkj,k_trp_Kkj))
     $     call errquit('aoccsd2: triples: k_trp_Kkj ', kchunk*lnov)
      if (.not. ma_push_get(MT_DBL, lnov,'k_trp_Dja',
     $     l_trp_Dja,k_trp_Dja))
     $     call errquit('aoccsd2: triples: k_trp_Dja ', lnov)
      if (.not. ma_push_get(MT_DBL, nvir*kchunk,'k_trp_Djka',
     $     l_trp_Djka,k_trp_Djka))
     $     call errquit('aoccsd2: triples: k_trp_Djka ', nvir*kchunk)
      if (.not. ma_push_get(MT_DBL, nvir,'k_trp_Djia',
     $     l_trp_Djia,k_trp_Djia))
     $     call errquit('aoccsd2: triples: k_trp_Djia ', nvir)
c

         call ga_sync
         tx(2)=tcgtime()
c
         call ccsd_trpdrv(dbl_mb(k_t1),dbl_mb(k_f1n),
     $        dbl_mb(k_f1t),dbl_mb(k_f2n),dbl_mb(k_f2t),
     $        dbl_mb(k_f3n),dbl_mb(k_f3t),dbl_mb(k_f4n),
     $        dbl_mb(k_f4t),eorb,eccsdt,g_objo,g_objv,
     $        g_coul,g_exch,ncor,nocc,nvir,iprt,
     $        empt(1),empt(2),oseg_lo,oseg_hi,kchunk,
     $        dbl_mb(k_trp_Tij), dbl_mb(k_trp_Tkj), dbl_mb(k_trp_Tia), 
     $        dbl_mb(k_trp_Tka), dbl_mb(k_trp_Xia), dbl_mb(k_trp_Xka), 
     $        dbl_mb(k_trp_Jia), dbl_mb(k_trp_Jka), dbl_mb(k_trp_Kia), 
     $        dbl_mb(k_trp_Kka), dbl_mb(k_trp_Jij), dbl_mb(k_trp_Jkj), 
     $        dbl_mb(k_trp_Kij), dbl_mb(k_trp_Kkj), dbl_mb(k_trp_Dja),
     $        dbl_mb(k_trp_Djka), dbl_mb(k_trp_Djia))
c
         tx(3)=tcgtime()

         if (.not. ma_chop_stack(l_f1n)) call errquit
     $        ('ccsd triples: stack corrupt (f1n)',0)

         if (iam.eq.0)write(6,'(1x,a,i5,f9.2)')
     &       'Time for integral evaluation pass ',
     $       nvpass,tx(2)-tx(1)
         if (iam.eq.0)write(6,'(1x,a,i5,f9.2)')
     &       'Time for triples evaluation pass  ',
     $      nvpass,tx(3)-tx(2)
      enddo
      empt(1)=empt(1)+empt(1)
      empt(2)=empt(2)+empt(2)
      call ga_dgop(msg_cc_emp2,empt,2, '+')
      eccsdt=empt(1)+empt(2)
      if (iam.eq.0)then
       write(6,100)empt(1),empt(2),eccsdt
 100   format(/,
     &     ' pseudo-e(mp4) ',e22.14,/,
     &     ' pseudo-e(mp5) ',e22.14,/,
     &     '        e(t)   ',e22.14,/)
      endif

      stat=ga_destroy(g_coul)
      stat=stat.and.ga_destroy(g_exch)
c
      If ( IAm .eq. 0) then
         call ecce_print1('total ccsd+t(ccsd) energy', 
     $        mt_dbl, refener+eccsd+empt(1), 1)
         call ecce_print1('total ccsd(t) energy', 
     $        mt_dbl, RefEner + eccsd + eccsdt, 1)
         Call Banner(LuOut,
     $      'CCSD(T) Energy', '-',
     $      .TRUE., .TRUE., .FALSE.)
         Write (LuOut, '(1X, A, F30.15)') 'Reference energy:          ',
     $      RefEner
         If ( Util_Print('byproduct energies', Print_Default)
     $      .OR. Inp_Compare(.FALSE., 'ccsd+t(ccsd)', CC_Theory) ) then
            Write (LuOut, '(1X)')
            Write (LuOut, '(1X, A, F30.15)')
     $         'CCSD corr. energy:         ', eccsd
            Write (LuOut, '(1X, A, F30.15)')
     $         'T(CCSD) corr. energy:      ', empt(1)
            write (LuOut, '(1X, A, F30.15)')
     $         'Total CCSD+T(CCSD) energy: ', RefEner + eccsd + empt(1)
         EndIf
         If ( Util_Print('byproduct energies', Print_Default)
     $      .OR. Inp_Compare(.FALSE., 'ccsd(t)', CC_Theory) ) then
            Write (LuOut, '(1X)')
            Write (LuOut, '(1X, A, F30.15)')
     $         'CCSD corr. energy:         ', eccsd
            Write (LuOut, '(1X, A, F30.15)')
     $         '(T) corr. energy:          ', eccsdt
            Write (LuOut, '(1X, A, F30.15)')
     $         'Total CCSD(T) energy:      ', RefEner + eccsd + eccsdt
         EndIf
         call util_flush(6)
      EndIf
C
C     Write our energy out to the database
C
      If ( .NOT. RTDB_Put( RTDB, 'ccsd:ccsd+t(ccsd) correlation energy',
     $   MT_Dbl, 1, Eccsd+empt(1) ) )
     $   Call ErrQuit('CCSD: Can''t write to RTDB', 10)
      If ( .NOT. RTDB_Put( RTDB, 'ccsd:ccsd+t(ccsd) energy', MT_Dbl, 1,
     $   RefEner+eccsd+empt(1) ) ) 
     $   Call ErrQuit('CCSD: Can''t write to RTDB', 11)
      If ( .NOT. RTDB_Put( RTDB, 'ccsd:ccsd(t) correlation energy',
     $   MT_Dbl, 1, Eccsd+eccsdt) )
     $   Call ErrQuit('CCSD: Can''t write to RTDB', 12)
      If ( .NOT. RTDB_Put( RTDB, 'ccsd:ccsd(t) energy', MT_Dbl, 1,
     $   RefEner+eccsd+eccsdt) ) 
     $   Call ErrQuit('CCSD: Can''t write to RTDB', 13)
C
C     This reports back to task_energy that the requested task has been
C     completed successfully.
C
      If ( Inp_Compare(.FALSE., 'ccsd(t)', CC_Theory) ) then
            call ecce_print1('total energy', mt_dbl, 
     $        RefEner+eccsd+eccsdt, 1)
         If ( .NOT. RTDB_Put( RTDB, 'ccsd:energy', MT_Dbl, 1,
     $      RefEner+eccsd+eccsdt) ) 
     $      Call ErrQuit('CCSD: Can''t write to RTDB', 14)
      ElseIf ( Inp_Compare(.FALSE., 'ccsd+t(ccsd)', CC_Theory) ) then
         call ecce_print1('total energy', mt_dbl, 
     $        RefEner+eccsd+empt(1), 1)
         If ( .NOT. RTDB_Put( RTDB, 'ccsd:energy', MT_Dbl, 1,
     $      RefEner+eccsd+empt(1) ) ) 
     $      Call ErrQuit('CCSD: Can''t write to RTDB', 14)
      EndIf
C
C
C
      stat=.true.
      stat=stat.and.ma_pop_stack(l_t1)
      if (.not.stat)call errquit('trpdrv: cannot pop stack?',0)
c
      if (.not. ga_destroy(g_objo))
     &     call errquit('g_objo ga_destroy failed',0)
c
c -- end of triples calculation
c
c      call ga_sync()
c      if (iam.eq.0)then
c       print *,' ga allocated end of aoccsd '
c       call ga_summarize(.false.)
c      endif
      call ga_sync()
      call qexit('triples',0)
c
 999  if (.not. ga_destroy(g_objv))
     &     call errquit('g_objv ga_destroy failed',0)
c
      call qstat()
      return
      end
