      subroutine aoccsd(basis,ncor,nocc,nvir,ndel,nbf,nact,maxit,
     &                  convi,iprt,cmo,eorb,blen,g_objv,
     &                  g_coul,g_exch,g_ncoul, g_nexch, RefEner,
     $   CC_Theory, rtdb, mxvec, geom)
      implicit none
c
      integer basis,ncor,nocc,nvir,ndel,nbf,nact,mxov,maxit,convi,iprt,
     &        blen,g_objv,g_coul,g_exch,g_ncoul,g_nexch, mxvec
      double precision cmo(nbf,nbf),eorb(nbf)
      Double Precision RefEner
      Character*(*) CC_Theory
      Integer RTDB, geom
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
      logical stat
      integer nsh,maxbfsh,max2e,mem2
      integer l_t1,l_t2,l_ht1,l_ht2,l_scra,l_scrb,
     &        l_hiu,l_hia,l_habe,l_giu,l_gabe,
     &        l_bbkp,l_bb
      integer k_t1,k_t2,k_ht1,k_ht2,k_scra,k_scrb,
     &        k_hiu,k_hia,k_habe,k_giu,k_gabe,
     &        k_bbkp,k_bb
      integer g_t2,g_ht2
c triples
      integer lbfo,lbfv,lbfscr,mem_avail
      integer l_f1n,l_f1t,l_f2n,l_f2t,l_f3n,l_f3t,l_f4n,l_f4t,
     &        l_buf1,l_buf2,l_buff,l_bufo,l_bufv,l_tk
      integer k_f1n,k_f1t,k_f2n,k_f2t,k_f3n,k_f3t,k_f4n,k_f4t,
     &        k_buf1,k_buf2,k_buff,k_bufo,k_bufv,k_tk
      integer g_objo
      integer nvpass,noper,occ_lo,occ_hi,vir_lo,vir_hi,oseg_lo,oseg_hi,
     &        npass
      integer ga_create_jkblocked
      logical oexch,ocoul,ohalf
      double precision eccsdt,escf,eccsd,empt(2),tx(3),rx(2)
      Integer Nodes, IAm
C
      double precision zip
      integer ii,jj
C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
C     Temporary stuff added to play density matrix games
C
      Integer A, I, Ad1, Ad2, J, B, LSab
      Integer L_Sb, K_Sb, L_Sa, K_Sa
      Integer H_DenZ1, I_DenZ1, H_DenT1, I_DenT1
      Integer G_DenT2, G_DenZ2, G_1PDM
      Logical Do_CCDen
      Data Do_CCDen /.FALSE./
      Integer GA_Create_Atom_Blocked
      External GA_Create_Atom_Blocked
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
* data must follow all delcarations 

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
      escf=zip
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
      convi=-convi
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
     &        /,2x,' convi  = ',i5,
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
      lnobj=lnoo*nbf
      if (iprt.gt.10)then
         write(6,*)'lnoo,lnov,lnvv,lnoov,lnooo,lnobj,lnovv',
     &              lnoo,lnov,lnvv,lnoov,lnooo,lnobj,lnovv
      endif
      stat=.true.
      stat=stat.and.ma_push_get(MT_DBL,lnov,'t1',l_t1,k_t1)
      stat=stat.and.ma_push_get(MT_DBL,2*lnoov,'t2',l_t2,k_t2)
      stat=stat.and.ma_push_get(MT_DBL,lnov,'ht1',l_ht1,k_ht1)
      stat=stat.and.ma_push_get(MT_DBL,lnoov,'ht2',l_ht2,k_ht2)
      stat=stat.and.ma_push_get(MT_DBL,lnobj,'scra',l_scra,k_scra)
      stat=stat.and.ma_push_get(MT_DBL,lnobj,'scrb',l_scrb,k_scrb)
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
      if (.not.stat)call errquit('ccsd: failed to allocate ',0)
c
      if (.not.ga_create(MT_DBL,mxvec*lnoov,nvir+1,'g_t2', 
     &                   mxvec*lnoov,1,g_t2))
     &     call errquit('ga_create failed',0)
      if (.not.ga_create(MT_DBL,mxvec*lnoov,nvir+1,'g_ht2', 
     &                   mxvec*lnoov,1,g_ht2))
     &     call errquit('ga_create failed',0)
c      
      call ccsd_iterdrv2(basis,ncor,nocc,nvir,ndel,nbf,nact,nsh,
     &                  maxit,convi,iprt,cmo,eorb,
     &                  dbl_mb(k_t1),dbl_mb(k_t2),
     &                  dbl_mb(k_ht1),dbl_mb(k_ht2),
     &                  dbl_mb(k_scra),dbl_mb(k_scrb),
     &                  dbl_mb(k_hiu),dbl_mb(k_hia),dbl_mb(k_habe),
     &                  dbl_mb(k_giu),dbl_mb(k_gabe),
     &                  dbl_mb(k_bbkp),dbl_mb(k_bb),
     &                  g_t2,g_ht2,g_objv,
     &                  mxvec,eccsd,max2e,mem2,g_ncoul,g_nexch,
     &                  int_mb(k_tk))
      if (.not. ga_destroy(g_ht2))
     &     call errquit('ga_destroy failed',0)
      if (.not. ga_destroy(g_t2))
     &     call errquit('ga_destroy failed',0)
      stat=.true.
      stat=stat.and.ma_pop_stack(l_tk)
      stat=stat.and.ma_pop_stack(l_bb)
      stat=stat.and.ma_pop_stack(l_bbkp)
      stat=stat.and.ma_pop_stack(l_gabe)
      stat=stat.and.ma_pop_stack(l_giu)
      stat=stat.and.ma_pop_stack(l_habe)
      stat=stat.and.ma_pop_stack(l_hia)
      stat=stat.and.ma_pop_stack(l_hiu)
      stat=stat.and.ma_pop_stack(l_scrb)
      stat=stat.and.ma_pop_stack(l_scra)
      stat=stat.and.ma_pop_stack(l_ht2)
      stat=stat.and.ma_pop_stack(l_ht1)
      stat=stat.and.ma_pop_stack(l_t2)
      stat=stat.and.ma_pop_stack(l_t1)
      if (.not.stat)call errquit('ccsd: cannot pop stack?',0)
      call ga_sync()
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
C
C     Do we continue to the triples or not?
C     Note: task_energy expects ccsd:energy as a marker of successful
C     completion.  Therefore, we don't want to write it unless we
C     are certain we have completed the expected task.
C
      If ( .NOT. ( Inp_Compare(.FALSE., 'ccsd(t)', CC_Theory)
     $   .OR. Inp_Compare(.FALSE., 'ccsd+t(ccsd)', CC_Theory) ) ) then
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
c$$$      rx(1)=rtc()
      call qenter('triples',0)
      if (iam.eq.0)write(6,293)
 293  format(//,'*********triples calculation*********',/)
      lbfo=lnovv+lnoov+lnoov
      lbfv=lnoov+lnoov+lnoov
      lbfscr=lnovv+lnovv
      lbfv=max(lbfv,nbf*nbf)
      lbfscr=max(lbfscr,3)
      lnobj=max(lbfv,lbfo)
      lnobj=max(lnobj,2*lnoov+lnooo)
      stat=.true.
      stat=stat.and.ma_push_get(MT_DBL,lnov,'t1',l_t1,k_t1)
      stat=stat.and.ma_push_get(MT_DBL,lnobj,'buf1',l_buf1,k_buf1)
      stat=stat.and.ma_push_get(MT_DBL,lnobj,'buf2',l_buf2,k_buf2)
      if (.not.ga_create(MT_DBL,lbfo,nocc,'objo',lbfo,1,g_objo))
     &     call errquit('ga_create failed',0)
      call ccsd_rdtrpo(dbl_mb(k_t1),dbl_mb(k_buf1),dbl_mb(k_buf2),
     $                 g_objo,g_objv,nocc,nvir,iprt)
      stat=.true.
      stat=stat.and.ma_pop_stack(l_buf2)
      stat=stat.and.ma_pop_stack(l_buf1)
      if (.not.stat)call errquit('rdtrpo: cannot pop stack?',0)
c
      stat=.true.
      stat=stat.and.ma_push_get(MT_DBL,lnvv,'f1n',l_f1n,k_f1n)
      stat=stat.and.ma_push_get(MT_DBL,lnvv,'f1t',l_f1t,k_f1t)
      stat=stat.and.ma_push_get(MT_DBL,lnvv,'f2n',l_f2n,k_f2n)
      stat=stat.and.ma_push_get(MT_DBL,lnvv,'f2t',l_f2t,k_f2t)
      stat=stat.and.ma_push_get(MT_DBL,lnvv,'f3n',l_f3n,k_f3n)
      stat=stat.and.ma_push_get(MT_DBL,lnvv,'f3t',l_f3t,k_f3t)
      stat=stat.and.ma_push_get(MT_DBL,lnvv,'f4n',l_f4n,k_f4n)
      stat=stat.and.ma_push_get(MT_DBL,lnvv,'f4t',l_f4t,k_f4t)
      stat=stat.and.ma_push_get(MT_DBL,lbfscr,'buff',l_buff,k_buff)
      stat=stat.and.ma_push_get(MT_DBL,lbfo,'bufo',l_bufo,k_bufo)
      stat=stat.and.ma_push_get(MT_DBL,lbfv,'bufv',l_bufv,k_bufv)
      if (.not.stat)call errquit('trpdrv: cannot allocate stack?',0)
c      print *,' ga allocated before trpdrv '
c      call ga_summarize(.true.)
c
      mem_avail = ma_inquire_avail(MT_DBL)
      call ga_sync()
      call ga_igop(msg_cc_mem,mem_avail,1,'min')
      mem_avail=nodes*mem_avail/2
      noper=mem_avail/(nbf*nbf)
      nvpass=min(nvir,noper/nocc)
      nvpass=max(nvpass,1)
      noper=nocc*nvpass
      npass=(nvir+nvpass-1)/nvpass
      if (iam.eq.0)write(*,845)mem_avail,nvir,nvpass,npass
  845 format(' memory available                     ',i8,/,
     &       ' total number of virtual orbitals     ',i5,/,
     &       ' number of virtuals per integral pass ',i5,/,
     &       ' number of integral evaluations       ',i5,/)
      g_coul = ga_create_JKblocked(noper,nbf,nbf,'Coulomb oper')
      g_exch = ga_create_JKblocked(noper,nbf,nvir,'X oper')
      npass=0
      occ_lo = ncor+1
      occ_hi = ncor+nocc
      vir_lo = ncor+nocc+1
      vir_hi = ncor+nocc+nvir
      oseg_hi = vir_lo - 1
      empt(1)=zip
      empt(2)=zip
      do while (oseg_hi.lt.vir_hi)
         npass = npass + 1
         oseg_lo = oseg_hi + 1
         oseg_hi = min(oseg_hi+nvpass,vir_hi)
         if (iam.eq.0)write(*,846)npass
 846     format(' commencing integral evaluation ',i5)
         tx(1)=tcgtime()
         call moints_trp(basis,  ohalf, occ_lo, oseg_lo, oseg_hi,
     &      vir_lo, vir_hi,
     &      g_coul, ocoul, g_exch, oexch,
     &      blen, cmo, ncor, ndel, nocc, nvir)
c        print *,'call trpdrv ',npass
c
         tx(2)=tcgtime()
         call ccsd_trpdrv(dbl_mb(k_t1),dbl_mb(k_buff),dbl_mb(k_bufo),
     $      lbfo,dbl_mb(k_bufv),lbfv,dbl_mb(k_f1n),
     $      dbl_mb(k_f1t),dbl_mb(k_f2n),dbl_mb(k_f2t),
     $      dbl_mb(k_f3n),dbl_mb(k_f3t),dbl_mb(k_f4n),
     $      dbl_mb(k_f4t),eorb,eccsdt,g_objo,g_objv,
     $      g_coul,g_exch,ncor,nocc,nvir,ndel,iprt,
     $      empt(1),empt(2),oseg_lo,oseg_hi)
         tx(3)=tcgtime()
         if (iam.eq.0)write(6,*)'Time for integral evaluation pass ',
     $      npass,tx(2)-tx(1)
         if (iam.eq.0)write(6,*)'Time for triples evaluation pass ',
     $      npass,tx(3)-tx(2)
      enddo
      empt(1)=empt(1)+empt(1)
      empt(2)=empt(2)+empt(2)
      call ga_dgop(msg_cc_emp2,empt,2, '+')
      eccsdt=empt(1)+empt(2)
      if (iam.eq.0)then
       write(6,100)empt(1),empt(2),eccsdt
 100   format(/,' pseudo-e(mp4) ',e22.14,/,
     &        ' pseudo-e(mp5) ',e22.14,/,' e(t)   ',e22.14,/)
      endif

      stat=ga_destroy(g_coul)
      stat=stat.and.ga_destroy(g_exch)
c
      If ( IAm .eq. 0) then
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
         If ( .NOT. RTDB_Put( RTDB, 'ccsd:energy', MT_Dbl, 1,
     $      RefEner+eccsd+eccsdt) ) 
     $      Call ErrQuit('CCSD: Can''t write to RTDB', 14)
      ElseIf ( Inp_Compare(.FALSE., 'ccsd+t(ccsd)', CC_Theory) ) then
         If ( .NOT. RTDB_Put( RTDB, 'ccsd:energy', MT_Dbl, 1,
     $      RefEner+eccsd+empt(1) ) ) 
     $      Call ErrQuit('CCSD: Can''t write to RTDB', 14)
      EndIf
C
C
C
      stat=.true.
      stat=stat.and.ma_pop_stack(l_bufv)
      stat=stat.and.ma_pop_stack(l_bufo)
      stat=stat.and.ma_pop_stack(l_buff)
      stat=stat.and.ma_pop_stack(l_f4t)
      stat=stat.and.ma_pop_stack(l_f4n)
      stat=stat.and.ma_pop_stack(l_f3t)
      stat=stat.and.ma_pop_stack(l_f3n)
      stat=stat.and.ma_pop_stack(l_f2t)
      stat=stat.and.ma_pop_stack(l_f2n)
      stat=stat.and.ma_pop_stack(l_f1t)
      stat=stat.and.ma_pop_stack(l_f1n)
      stat=stat.and.ma_pop_stack(l_t1)
      if (.not.stat)call errquit('trpdrv: cannot pop stack?',0)
c
      if (.not. ga_destroy(g_objo))
     &     call errquit('g_obj ga_destroy failed',0)
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
c$$$      rx(2)=rtc()
c$$$      if (iam.eq.0)write(6,*)'Wall-clock time for triples ',
c$$$     $                        (rx(2)-rx(1))*3.33e-09
  999 continue
      call qstat()
      return
      end
