c#define DOGRAD 
      subroutine aoccsd(basis,ncor,nocc,nvir,ndel,nact,nbf,maxit,
     &                  convi,iprt,cmo,eorb,blen,
     &                  g_ncoul, g_nexch, RefEner,
     $   CC_Theory, rtdb, mxvec, geom, Tol2e, occd, oconverged)
      implicit none
#include "errquit.fh"
C     $Id$
c
      integer basis,ncor,nocc,nvir,ndel,nact,nbf,maxit,iprt,
     &        blen,g_objv,g_coul,g_exch,g_ncoul,g_nexch, mxvec
      double precision cmo(nbf,nbf),eorb(nbf)
      Double Precision RefEner, Convi, Tol2E
      Character*(*) CC_Theory
      Integer RTDB, geom
      logical oconverged, occd, use_trpdrv_nb
      logical use_trpdrv_omp, use_trpdrv_bgp2
      logical use_trpdrv_omp_mp
      logical use_trpdrv_offload
c
#include "ccsd_len.fh"
#include "tcgmsg.fh"
#include "bas.fh"
#include "rtdb.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"
#include "stdio.fh"
#include "msgids.fh"
#include "ccsdps.fh"
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
#ifdef USE_F90_ALLOCATABLE
      double precision, allocatable :: f1n(:,:),f1t(:,:)
      double precision, allocatable :: f2n(:,:),f2t(:,:)
      double precision, allocatable :: f3n(:,:),f3t(:,:)
      double precision, allocatable :: f4n(:,:),f4t(:,:)
#ifdef USE_FASTMEM
!dec$ attributes fastmem :: f1n,f2n,f3n,f4n,f1t,f2t,f3t,f4t
#endif
      integer :: alloc_error
#else
      integer l_f1n,l_f1t,l_f2n,l_f2t,l_f3n,l_f3t,l_f4n,l_f4t
      integer k_f1n,k_f1t,k_f2n,k_f2t,k_f3n,k_f3t,k_f4n,k_f4t
#endif
      integer l_buf1,l_buf2,l_tk,k_buf1,k_buf2,k_tk
      integer g_objo
      integer vchunk,noper,occ_lo,vir_lo,vir_hi,oseg_lo,oseg_hi,
     &        nvpass
      integer ga_create_jkblocked
      logical oexch,ocoul,ohalf,oprint,oprinthigh
      double precision eccsdt,eccsd,empt(2),tx(3)
      double precision eccsdos,eccsdss,fos,fss
      Integer Nodes, IAm, memdrv, memtrn

      logical use_ccsd_omp
#if defined(USE_OPENMP)
      integer  omp_get_thread_num
      external omp_get_thread_num
      integer  omp_get_num_threads
      external omp_get_num_threads
      integer  omp_get_max_threads
      external omp_get_max_threads
#endif
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
      integer gmem_needed,gmem_avail,nvnnodes,ff,ff2,ff3,iin,jin
      integer memlimit
      integer MSG_MIN_LOCAL_MEM
      parameter(MSG_MIN_LOCAL_MEM=1971)
      integer nkpass0
      external ccsd_createg
      Logical Do_Z,ccsd_createg
      integer nthrd_org,ccsd_omp_threads
#ifdef USE_F90INTERFACE
      interface
      subroutine ccsd_iterdrv2(rtdb,basis,nsh,ncor,nocc,nvir,nact,nbf,
     &                  maxit,convi,iprt,cmo,eorb,
     &                  t1,ht1,
     &                  scra,scrb,
     &                  hiu,hia,habe,
     &                  giu,gabe,
     &                  bbkp,bb,
     &                  g_t2,g_ht2,
     &                  mxvec,eccsd,max2e,mem2,g_ncoul,g_nexch,
     &                  tk,offt2,nvp1, refener, Tol2E, 
     &                  occd, oconverged, eccsdos, use_ccsd_omp)
      integer rtdb,basis,nsh,ncor,nocc,nvir,nact,nbf,maxit,iprt,
     & g_t2,g_ht2,g_nt2,g_nht2,mxvec,g_ncoul,g_nexch,
     &     tklst
      Double Precision Thre, Tol2E
      integer offt2
      double precision time_iter_start,convi
      integer time_prev_iter
      integer tk
      logical oconverged        
      logical occd
      logical, optional, intent(in) :: use_ccsd_omp
      double precision
     &                  t1,ht1,
     &                  scra,scrb,
     &                  hiu,hia,habe,
     &                  giu,gabe,
     &                  bbkp,bb
      double precision cmo(nbf,*),eorb(*)
      double precision eccsd,eccsdos,refener
      end subroutine
      end interface
#endif      
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
      if (occsdps) then
         call pstat_on(ps_aoccsd)
      else
         call qenter('aoccsd',0)
      endif
      oprint = util_print('information',print_low)
      oprinthigh = util_print('ccsd extra',print_high)
c
      nodes = ga_nnodes()
      iam = ga_nodeid()
      stat=bas_numcont(basis,nsh)
      stat=stat.and.bas_nbf_cn_max(basis,maxbfsh)
      if (.not.stat) call errquit('moints: cannot get basis info',0,
     &       BASIS_ERR)
c
c      print *,'basis set handle ',basis
c      print *,' ga allocated start of aoccsd '
c      call ga_summarize(.true.)
c
      if (iam.eq.0.and.oprint)then
      write(6,1001)nodes
 1001 format(//,1x,76(1h*),/,14x,
     &       'the segmented parallel ccsd program:',
     &       i5,' nodes',/,1x,76(1h*),//)
      endif
      if (iam.eq.0.and.oprint)then
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
     &        /,2x,' convi  = ',e10.3,
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
      if(iam.eq.0.and.oprint)write(6,*)'memory',mem_avail
      call int_mem_2e4c(max2e, mem2)
      mem2=max(mem2,nbf*nbf)
      if (iprt.gt.10.and.iam.eq.0.and.oprint)then
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
      if (.not.stat)call errquit('ccsd: failed to allocate ',0, MEM_ERR)
c
        lhoo=nocc*(nocc+1)/2
        lhoov=nvir*nocc*(nocc+1)/2
c
c     check if enough GA is there for the next 3 gacreate
c
       gmem_avail=(ga_memory_avail()/MA_sizeof(MT_DBL,1,MT_BYTE))
       nvnnodes=max(1,(nvir+1)/ga_nnodes()+1)
       call ga_igop(MSG_MIN_LOCAL_MEM,gmem_avail,1,'min')
c     mem for t2
       gmem_needed=2*mxvec*lhoov*nvnnodes
c     mem for nt2 (see ccsd_iterdrv2)
       ff=lnov/nvir/int(dsqrt(dble(ga_nnodes())))+2
       gmem_needed = gmem_needed + 2*nvir*nvir*ff*ff
c     mem needed for x c sht2 (see pampt3)
       ff2=max(1,nocc*nocc/ga_nnodes()+1)
       ff3=max(1,nbf*nbf/ga_nnodes()+1)
       gmem_needed = gmem_needed + 2*ff2*nbf*nbf+ff3*nocc*nocc
       
       if(gmem_needed.gt.gmem_avail) then
                mxvec=min(10,((gmem_avail-ff3*nocc*nocc)/2
     -         -(nvir*nvir*ff*ff)-
     -         ff2*nbf*nbf)/
     /         (lhoov*nvnnodes))
                if(ga_nodeid().eq.0) then
                   write(luout,*) ' ff = ',ff
                   write(luout,*) ' mxvec decreased to ',mxvec
                endif
                if(mxvec.eq.0) then
                if(ga_nodeid().eq.0) then
                   write(luout,'(A,2I16)')
     W               ' gmems (dbls): needed and avail ',
     W                  gmem_needed,gmem_avail
                   write(luout,'(A,2I10)') ' lhoov nv ',lhoov,nvir
                   write(luout,'(A,I16)') ' gneede for mxvec=1 ',
     W                  nvnnodes*(2*lhoov+lhoov)
                   call util_flush(luout)
                endif
               
                 call errquit(' bailing out',0,0)
                endif
       endif
       

*ga:1:0
      if (.not.ga_create(MT_DBL,mxvec*lhoov,nvir+1,'g_t2', 
     &                   mxvec*lhoov,0,g_t2))
     &     call errquit('ga_create failed for g_t2',0, GA_ERR)
*ga:1:0
      g_ht2=123456789
#ifdef DOGRAD
      if (.not.ga_create(MT_DBL,mxvec*lhoov,nvir+1,'g_ht2', 
     &                   mxvec*lhoov,0,g_ht2))
     &     call errquit('ga_create failed for g_ht2',0, GA_ERR)
#endif
!
! Determine if OpenMP should be used in CCSD code
!
      if (.not.rtdb_get(rtdb,'ccsd:use_ccsd_omp',mt_log,1,use_ccsd_omp))
     &         use_ccsd_omp = .false. 
      if (.not.rtdb_get(rtdb,'ccsd:ccsd_omp_threads',mt_int,1,
     A     ccsd_omp_threads)) ccsd_omp_threads=0
#if defined(USE_OPENMP)
c hack to avoid problems with ompi 4 ucx on tahoma
        if (use_ccsd_omp.and.ccsd_omp_threads.ne.0) then
          nthrd_org= omp_get_max_threads()
          call omp_set_num_threads(ccsd_omp_threads)
       endif
#endif
      if (ga_nodeid().eq.0) then
#if defined(USE_OPENMP)
        if (use_ccsd_omp) then
           write(6,*) ' nthrd_org ',nthrd_org
          write(6,33) omp_get_max_threads()
        else
          write(6,33) 1
        endif
#else
        write(6,33) 1
        if (use_ccsd_omp) then
          write(6,333)
        endif
#endif
      endif
   33 format(2x,'Using ',i2,' OpenMP thread(s) in CCSD')
  333 format(2x,'Recompile w/ USE_OPENMP=1 to use threads in CCSD')
!
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
     &                  occd, oconverged, eccsdos, use_ccsd_omp)
!     
      call ga_sync()
#ifdef USE_OPENMP
c hack to avoid problems with ompi 4 ucx on tahoma
        if (use_ccsd_omp.and.ccsd_omp_threads.ne.0) then
          call omp_set_num_threads(nthrd_org)
          if (ga_nodeid().eq.0) then
             write(6,*) 'OMP threads from now on is',
     C            omp_get_max_threads()
          endif
       endif
#endif
!
         if (.not.ga_create(MT_DBL,lnov,lnov,'g_nt2',
     &      nvir,nvir,g_nt2))
     &      call errquit('ga_create failed g_nt2',0, GA_ERR)
         call ga_zero(g_nt2)
      if(ga_nodeid().eq.0.and.oprinthigh) then
         call ga_distribution(g_nt2, 0, ilo, ihi, jlo, jhi)
         write(luout,'(A,4I12)') ' distr for nt2 ',
     I        ilo, ihi, jlo, jhi
         call util_flush(luout)
      endif
c
#if 0
      call ga_get(g_t2,offt2+1,offt2+lnov,nvp1,nvp1,
     $            dbl_mb(k_t1),lnov)
#else
      if (ga_nodeid().eq.0) then
          call ga_get(g_t2,offt2+1,offt2+lnov,nvp1,nvp1,
     &                dbl_mb(k_t1),lnov)
      endif
      call ga_brdcst(65537, dbl_mb(k_t1),
     &               ma_sizeof(MT_DBL,1,MT_BYTE)*lnov, 0)
#endif
! TODO use OpenMP here
      call ga_distribution(g_t2,iam,jlo,jhi,ilo,ihi)
      do a=1,nvir
        if (a.ge.ilo.and.a.le.ihi)then
          call ga_get(g_t2,offt2+1,offt2+lhoov,a,a,
     &                 dbl_mb(k_scra),lhoov)
! staggering
          do iin=ga_nodeid()+1,ga_nodeid()+nocc
            i=mod(iin,nocc)+1
            call dcopy(lnov,0.0d00,0,dbl_mb(k_scrb),1)
            do jin=1+ga_nodeid(),i+ga_nodeid()
              j=mod(jin,i)+1
              do b=1,nvir
                ad1=(b-1)*lhoo+i*(i-1)/2+j
                dbl_mb(k_scrb+b-1)=dbl_mb(k_scra+ad1-1)
              enddo
            ad1=(i-1)*nvir+a
            call ga_put(g_nt2,(j-1)*nvir+1,j*nvir,ad1,ad1,
     &                  dbl_mb(k_scrb),nvir)
            enddo
          enddo
        endif
      enddo
      call ga_sync()
! TODO use OpenMP here
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
      call ga_sync()
c
         if (iam.eq.0.and.iprt.gt.5.and.oprint)then
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
      if(ga_nnodes().lt.3) then
         if (.not.ga_create(MT_DBL,lnobj,nvir+1,'objv',lnobj,1,
     &        g_objv))
     &     call errquit('ga_create g_objv failed',0, GA_ERR)
      else
         if(.not.ccsd_createg(g_objv, lnobj, nvir+1, 
     R        lnoov, max(lnoov,lnooo), lnoov,min(nocc,nvir),
     O        oprinthigh))
     &     call errquit('ga_create g_objv failed',0, GA_ERR)
      endif
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
#ifdef DOGRAD
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
#else
         call errquit('ccsd gradients not enabled',0,CAPMIS_ERR)
#endif
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
      if (.not.stat)call errquit('ccsd: cannot pop stack?',0,MA_ERR)
      call ga_sync()
      if (.not. ga_destroy(g_nt2))
     &     call errquit('ga_destroy failed',0, GA_ERR)
      if(g_ht2.ne.123456789) then
      if (.not. ga_destroy(g_ht2))
     &     call errquit('ga_destroy failed',0, GA_ERR)
      g_ht2=123456789
      endif
      if (.not. ga_destroy(g_t2))
     &     call errquit('ga_destroy failed',0, GA_ERR)
c
      if (occsdps) then
         call pstat_off(ps_aoccsd)
      else
         call qexit('aoccsd',0)
      endif

c    get scaling factors for SCS-CCSD
      if (.not. rtdb_get(rtdb, 'ccsd:fss', mt_dbl, 1, fss))
     $  fss = 1.13d0        ! same spin scaling factor
      if (.not. rtdb_get(rtdb, 'ccsd:fos', mt_dbl, 1, fos))
     $  fos = 1.27d0        ! opposite spin scaling factor

c
      If ( GA_NodeID() .eq. 0 .AND.
     $     Util_Print('information', Print_Low) .AND.
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

c       Spin Componenst Scaled CCSD

         write(LuOut,'(/)')
         eccsdss=eccsd-eccsdos
         Call Banner(LuOut,
     $      'Spin Component Scaled (SCS) CCSD', '-',
     $      .TRUE., .TRUE., .FALSE.)
         Write (LuOut, '(1X, A, F30.15)')'Same spin contribution:     ',
     $      eccsdss
         Write (LuOut, '(1X, A, F30.15)')'Same spin scaling factor:   ',
     $      fss
         Write (LuOut, '(1X, A, F30.15)')'Opposite spin contribution: ',
     $      eccsdos
         Write (LuOut, '(1X, A, F30.15)')'Opposite spin scaling fact.:',
     $      fos
         Write (LuOut, '(1X, A, F30.15)')'SCS-CCSD correlation energy:',
     $      eccsdos*fos+eccsdss*fss
         Write (LuOut, '(1X, A, F30.15)')'Total SCS-CCSD energy:      ',
     $      RefEner+eccsdos*fos+eccsdss*fss

      EndIf
C
C     Write our energy out to the database
C
      If ( .NOT. RTDB_Put( RTDB, 'ccsd:ccsd correlation energy',
     $   MT_Dbl, 1, Eccsd) )
     $   Call ErrQuit('CCSD: Can''t write to RTDB', 1,0)
      If ( .NOT. RTDB_Put( RTDB, 'ccsd:ccsd energy', MT_Dbl, 1,
     $   RefEner+eccsd) ) 
     $   Call ErrQuit('CCSD: Can''t write to RTDB', 2,0)
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
     $      'CCDen (prep): unable to allocate T1', NOcc*NVir,0)
C
         If ( .NOT. MA_Push_Get(MT_Dbl, NOcc*NVir, 'CCDens Z1',
     $      H_DenZ1, I_DenZ1) ) Call ErrQuit(
     $      'CCDen (prep): unable to allocate Z1', NOcc*NVir,0)
C
         If ( .NOT. GA_Create(MT_Dbl, NVir*NOcc, NVir*NOcc,
     $      'CCDen T2', NVir, NVir, G_DenT2) ) Call ErrQuit(
     $      'CCDen(prep): unable to alloc T2', NVir*NOcc*NVir*NOcc,0)
C
         If ( .NOT. GA_Create(MT_Dbl, NVir*NOcc, NVir*NOcc,
     $      'CCDen Z2', NVir, NVir, G_DenZ2) ) Call ErrQuit(
     $      'CCDen (prep): unable to alloc Z2',NVir*NOcc*NVir*NOcc,0)
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
     $         Call ErrQuit('CCDen (prep): can''t allocate sa', lsab,0)
            If (.NOT. ma_push_get(MT_DBL,lsab,'sb',l_sb, k_sb) )
     $         Call ErrQuit('CCDen (prep): can''t allocate sb', lsab,0)
! TODO use OpenMP here
            do a=1,nvir
               call ga_get(g_ObjV,lnobj+1, lnobj+lnoov,a,a,
     &                     dbl_mb(k_sa),lnoov)
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
     $         'CCDen (prep): freeing sb', 0,0)
            If ( .NOT. MA_Pop_Stack(L_Sa)) Call ErrQuit(
     $         'CCDen (prep): freeing sa', 0,0)
         EndIf
         Call GA_Sync()
         Call GA_Copy(G_DenT2, G_DenZ2)
C
C        Now call the Density driver routine
C
         Call CCDen_Driver(NOcc, NVir, Dbl_MB(I_DenT1), NOcc, G_DenT2,
     $      Dbl_MB(I_DenZ1), NOcc, G_DenZ2, G_1PDM)
C
         If ( .NOT. GA_Destroy(G_1PDM) ) Call ErrQuit(
     $      'CCDen (prep): unable to free 1PDM', 0,0)
         If ( .NOT. GA_Destroy(G_DenZ2) ) Call ErrQuit(
     $      'CCDen (prep): unable to free Z2', 0,0)
         If ( .NOT. GA_Destroy(G_DenT2) ) Call ErrQuit(
     $      'CCDen (prep): unable to free T2', 0,0)
         If ( .NOT. MA_Pop_Stack(H_DenZ1) ) Call ErrQuit(
     $      'CCDen (prep): unable to free Z1', 0,0)
         If ( .NOT. MA_Pop_Stack(H_DenT1) ) Call ErrQuit(
     $      'CCDen (prep): unable to free T1', 0,0)
      EndIf
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
c
      if (.not. oconverged) goto 999
c
C     Do we continue to the triples or not?
C     Note: task_energy expects ccsd:energy as a marker of successful
C     completion.  Therefore, we do not want to write it unless we
C     are certain we have completed the expected task.
C
      If ( .NOT. ( Inp_Compare(.FALSE., 'ccsd(t)', CC_Theory)
     $   .OR. Inp_Compare(.FALSE., 'ccsd+t(ccsd)', CC_Theory) ) ) then
         call ecce_print1('total energy', mt_dbl, refener+eccsd, 1)
         If ( .NOT. RTDB_Put( RTDB, 'ccsd:energy', MT_Dbl, 1,
     $      RefEner+eccsd) ) 
     $      Call ErrQuit('CCSD: Cannot write to RTDB', 3,0)
         Goto 999
      EndIf
!
! ----------------------------
! the new parallel triples
! ----------------------------
!
! Determine if OpenMP, etc. should be used in (T) code
!
      if (.not. rtdb_get(rtdb, 'ccsd:use_trpdrv_nb', mt_log, 1,
     1                   use_trpdrv_nb))
     2    use_trpdrv_nb=.false.
      if (.not. rtdb_get(rtdb, 'ccsd:use_trpdrv_omp', mt_log, 1,
     1                   use_trpdrv_omp))
     2    use_trpdrv_omp=.false.
      if (.not. rtdb_get(rtdb, 'ccsd:use_trpdrv_omp_mp', mt_log, 1,
     1                   use_trpdrv_omp_mp))
     2    use_trpdrv_omp_mp=.false.
      if (.not. rtdb_get(rtdb, 'ccsd:use_trpdrv_offload', mt_log, 1,
     1                   use_trpdrv_offload))
     2    use_trpdrv_offload=.false.
      if (.not. rtdb_get(rtdb, 'ccsd:use_trpdrv_bgp2', mt_log, 1,
     1                   use_trpdrv_bgp2))
     2    use_trpdrv_bgp2=.false.
!
      mem_avail = ma_inquire_avail(MT_DBL)
      if(iam.eq.0.and.oprint)write(luout,*)'memory',mem_avail
      if (occsdps) then
         call pstat_on(ps_tripls)
      else
         call qenter('triples',0)
      endif
      if (iam.eq.0.and.oprint)write(luout,293)
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
      if(ga_nnodes().lt.3) then
      if (.not.ga_create(MT_DBL,lbfo,nocc,'objo',lbfo,1,g_objo))
     &     call errquit('ga_create failed for g_objo',0, GA_ERR)
      else
         if(.not.ccsd_createg(g_objo, lbfo, nocc, lnovv, lnoov, lnoov,
     &                        min(nocc,nvir), oprinthigh))
     &     call errquit('ga_create g_objo failed',0, GA_ERR)
      endif
c
c
      call ga_zero(g_objo)
      call ccsd_rdtrpo(dbl_mb(k_t1),dbl_mb(k_buf1),dbl_mb(k_buf2),
     &                 g_objo,g_objv,nocc,nvir,iprt)
      stat=.true.
      stat=stat.and.ma_pop_stack(l_buf2)
      stat=stat.and.ma_pop_stack(l_buf1)
      if (.not.stat)call errquit('rdtrpo: cannot pop stack?',0, MA_ERR)
c
c      print *,' ga allocated before trpdrv '
c      call ga_summarize(.true.)
c
c     Iteratively determine nvpass (vchunk) and 
c     internal chunking for the tripes (kchunk)
c
      mem_avail = ma_inquire_avail(MT_DBL)
c     we can hack the mem available (in MB)
      if ( rtdb_get(rtdb, 'ccsdt:memlimit', mt_int, 1, memlimit))
     $  mem_avail=min(mem_avail,memlimit*1024*1024/
     /     MA_sizeof(MT_DBL,1,MT_BYTE))

      call ga_igop(msg_cc_mem, mem_avail, 1, 'min')
c
      nvpass = 1
      if (.not.rtdb_get(rtdb,'ccsd:nkpass',mt_int,1,nkpass0))
     &         nkpass0=1
      if (.not.rtdb_get(rtdb,'ccsd:nvpass',mt_int,1,nvpass))
     &         nvpass=1
c
 7651 nkpass = nkpass0
 7652 kchunk = (nocc - 1)/nkpass + 1
      vchunk = (nvir - 1)/nvpass + 1
c
cold      memdrv = lnvv*(11+3*kchunk) + lnov*(5+4*kchunk) + nvir*(1+kchunk)
      memdrv = lnvv*(11+3*kchunk) + lnov*(1+4*nocc+4*kchunk) +
     +     nvir*(1+kchunk)
      memtrn = moints_trp_mem(basis, vchunk, blen)
c
      if (iam .eq. 0.and.oprint) then
         write(luout,1101) nkpass, nvpass,memdrv,memtrn,mem_avail
 1101    format('nkpass=',i5,'; nvpass=',i5,'; memdrv=',i15,
     $        '; memtrn=',i15,'; memavail=',i15)
         call util_flush(luout)
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
     $           ('ccsd triples: too many k passes', nkpass,
     &       UNKNOWN_ERR)
            goto 7652
         else
            nvpass = nvpass + 1
            if (nvpass .gt. nvir) call errquit
     $           ('ccsd triples: too many v passes ', nvpass,
     &       UNKNOWN_ERR)
            goto 7651
         endif
      endif
c
      noper=nocc*vchunk
c
      if (iam.eq.0.and.oprint)write(*,845)mem_avail,nvir,vchunk,nvpass,
     $     kchunk,nkpass
  845 format(' memory available/node                ',i15,/,
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

#ifdef USE_F90_ALLOCATABLE
      allocate( f1n(1:nvir,1:nvir), stat=alloc_error)
      if (alloc_error.ne.0) call errquit('f1n',1,MA_ERR)
      allocate( f2n(1:nvir,1:nvir), stat=alloc_error)
      if (alloc_error.ne.0) call errquit('f2n',2,MA_ERR)
      allocate( f3n(1:nvir,1:nvir), stat=alloc_error)
      if (alloc_error.ne.0) call errquit('f3n',3,MA_ERR)
      allocate( f4n(1:nvir,1:nvir), stat=alloc_error)
      if (alloc_error.ne.0) call errquit('f4n',4,MA_ERR)
      allocate( f1t(1:nvir,1:nvir), stat=alloc_error)
      if (alloc_error.ne.0) call errquit('f1t',5,MA_ERR)
      allocate( f2t(1:nvir,1:nvir), stat=alloc_error)
      if (alloc_error.ne.0) call errquit('f2t',6,MA_ERR)
      allocate( f3t(1:nvir,1:nvir), stat=alloc_error)
      if (alloc_error.ne.0) call errquit('f3t',7,MA_ERR)
      allocate( f4t(1:nvir,1:nvir), stat=alloc_error)
      if (alloc_error.ne.0) call errquit('f4t',8,MA_ERR)
#else
      if (.not. ma_push_get(MT_DBL,lnvv,'f1n',l_f1n,k_f1n))
     $     call errquit('aoccsd2:triples f1n',lnvv, MA_ERR)
      if (.not. ma_push_get(MT_DBL,lnvv,'f1t',l_f1t,k_f1t))
     $     call errquit('aoccsd2:triples f1t',lnvv, MA_ERR)
      if (.not. ma_push_get(MT_DBL,lnvv,'f2n',l_f2n,k_f2n))
     $     call errquit('aoccsd2:triples f2n',lnvv, MA_ERR)
      if (.not. ma_push_get(MT_DBL,lnvv,'f2t',l_f2t,k_f2t))
     $     call errquit('aoccsd2:triples f2t',lnvv, MA_ERR)
      if (.not. ma_push_get(MT_DBL,lnvv,'f3n',l_f3n,k_f3n))
     $     call errquit('aoccsd2:triples f3n',lnvv, MA_ERR)
      if (.not. ma_push_get(MT_DBL,lnvv,'f3t',l_f3t,k_f3t))
     $     call errquit('aoccsd2:triples f3t',lnvv, MA_ERR)
      if (.not. ma_push_get(MT_DBL,lnvv,'f4n',l_f4n,k_f4n))
     $     call errquit('aoccsd2:triples f4n',lnvv, MA_ERR)
      if (.not. ma_push_get(MT_DBL,lnvv,'f4t',l_f4t,k_f4t))
     $     call errquit('aoccsd2:triples f4t',lnvv, MA_ERR)
#endif

      if (.not. ma_push_get(MT_DBL, lnvv,'k_trp_Tij',
     $        l_trp_Tij,k_trp_Tij))
     $     call errquit('aoccsd2: triples: k_trp_Tij ', lnvv, MA_ERR)
      if (.not. ma_push_get(MT_DBL, kchunk*lnvv,'k_trp_Tkj',
     $     l_trp_Tkj,k_trp_Tkj))
     $     call errquit('aoccsd2: triples: k_trp_Tkj ', kchunk*lnvv,
     &       MA_ERR)
      if (.not. ma_push_get(MT_DBL, lnov*nocc,'k_trp_Tia',
     $     l_trp_Tia,k_trp_Tia))
     $     call errquit('aoccsd2: triples: k_trp_Tia ', lnov, MA_ERR)
      if (.not. ma_push_get(MT_DBL, kchunk*lnov,'k_trp_Tka',
     $     l_trp_Tka,k_trp_Tka))
     $     call errquit('aoccsd2: triples: k_trp_Tka ', kchunk*lnov,
     &       MA_ERR)
      if (.not. ma_push_get(MT_DBL, lnov*nocc,'k_trp_Xia',
     $     l_trp_Xia,k_trp_Xia))
     $     call errquit('aoccsd2: triples: k_trp_Xia ', lnov, MA_ERR)
      if (.not. ma_push_get(MT_DBL, kchunk*lnov,'k_trp_Xka',
     $     l_trp_Xka,k_trp_Xka))
     $     call errquit('aoccsd2: triples: k_trp_Xka ', kchunk*lnov,
     &       MA_ERR)
      if (.not. ma_push_get(MT_DBL, lnvv,'k_trp_Jia',
     $     l_trp_Jia,k_trp_Jia))
     $     call errquit('aoccsd2: triples: k_trp_Jia ', lnvv, MA_ERR)
      if (.not. ma_push_get(MT_DBL, kchunk*lnvv,'k_trp_Jka',
     $     l_trp_Jka,k_trp_Jka))
     $     call errquit('aoccsd2: triples: k_trp_Jka ', kchunk*lnvv,
     &       MA_ERR)
      if (.not. ma_push_get(MT_DBL, lnvv,'k_trp_Kia',
     $     l_trp_Kia,k_trp_Kia))
     $     call errquit('aoccsd2: triples: k_trp_Kia ', lnvv, MA_ERR)
      if (.not. ma_push_get(MT_DBL, kchunk*lnvv,'k_trp_Kka',
     $     l_trp_Kka,k_trp_Kka))
     $     call errquit('aoccsd2: triples: k_trp_Kka ', kchunk*lnvv,
     &       MA_ERR)
      if (.not. ma_push_get(MT_DBL, lnov*nocc,'k_trp_Jij',
     $     l_trp_Jij,k_trp_Jij))
     $     call errquit('aoccsd2: triples: k_trp_Jij ', lnov, MA_ERR)
      if (.not. ma_push_get(MT_DBL, kchunk*lnov,'k_trp_Jkj',
     $     l_trp_Jkj,k_trp_Jkj))
     $     call errquit('aoccsd2: triples: k_trp_Jkj ', kchunk*lnov,
     &       MA_ERR)
      if (.not. ma_push_get(MT_DBL, lnov*nocc,'k_trp_Kij',
     $     l_trp_Kij,k_trp_Kij))
     $     call errquit('aoccsd2: triples: k_trp_Kij ', lnov, MA_ERR)
      if (.not. ma_push_get(MT_DBL, kchunk*lnov,'k_trp_Kkj',
     $     l_trp_Kkj,k_trp_Kkj))
     $     call errquit('aoccsd2: triples: k_trp_Kkj ', kchunk*lnov,
     &       MA_ERR)
      if (.not. ma_push_get(MT_DBL, lnov,'k_trp_Dja',
     $     l_trp_Dja,k_trp_Dja))
     $     call errquit('aoccsd2: triples: k_trp_Dja ', lnov, MA_ERR)
      if (.not. ma_push_get(MT_DBL, nvir*kchunk,'k_trp_Djka',
     $     l_trp_Djka,k_trp_Djka))
     $     call errquit('aoccsd2: triples: k_trp_Djka ', nvir*kchunk,
     &       MA_ERR)
      if (.not. ma_push_get(MT_DBL, nvir*nocc,'k_trp_Djia',
     $     l_trp_Djia,k_trp_Djia))
     $     call errquit('aoccsd2: triples: k_trp_Djia ', nvir, MA_ERR)
c
      do while (oseg_hi.lt.vir_hi)
         nvpass = nvpass + 1
         oseg_lo = oseg_hi + 1
         oseg_hi = min(oseg_hi+vchunk,vir_hi)
         tx(1)=tcgtime()
         if (iam.eq.0)write(*,846)nvpass,tx(1)
 846     format(' commencing integral evaluation ',i8,' at ',f14.2)
         call moints_trp(rtdb,basis,  ohalf, occ_lo, oseg_lo, oseg_hi,
     &      vir_lo, vir_hi,
     &      g_coul, ocoul, g_exch, oexch,
     &      blen, cmo, ncor, nocc, nvir, ndel, Tol2E)
c        print *,'call trpdrv ',nvpass
         call ga_sync()
         tx(2)=tcgtime()
c
         if (use_trpdrv_omp) then
c
         if (iam.eq.0.and.oprint) then
            write(luout,1848) nvpass,util_wallsec()
            call util_flush(luout)
         endif
 1848    format(' commencing triples evaluation - OpenMP version',i8,
     I        ' at ',f20.2,' secs')
         call ccsd_trpdrv_omp(dbl_mb(k_t1),
#ifdef USE_F90_ALLOCATABLE
     $        f1n,f1t,f2n,f2t,f3n,f3t,f4n,f4t,
#else
     $        dbl_mb(k_f1n),dbl_mb(k_f1t),dbl_mb(k_f2n),dbl_mb(k_f2t),
     $        dbl_mb(k_f3n),dbl_mb(k_f3t),dbl_mb(k_f4n),dbl_mb(k_f4t),
#endif
     $        eorb,g_objo,g_objv,g_coul,g_exch,ncor,nocc,nvir,iprt,
     $        empt(1),empt(2),oseg_lo,oseg_hi,kchunk,
     $        dbl_mb(k_trp_Tij), dbl_mb(k_trp_Tkj), dbl_mb(k_trp_Tia),
     $        dbl_mb(k_trp_Tka), dbl_mb(k_trp_Xia), dbl_mb(k_trp_Xka),
     $        dbl_mb(k_trp_Jia), dbl_mb(k_trp_Jka), dbl_mb(k_trp_Kia),
     $        dbl_mb(k_trp_Kka), dbl_mb(k_trp_Jij), dbl_mb(k_trp_Jkj),
     $        dbl_mb(k_trp_Kij), dbl_mb(k_trp_Kkj), dbl_mb(k_trp_Dja),
     $        dbl_mb(k_trp_Djka), dbl_mb(k_trp_Djia))
c
         elseif (use_trpdrv_omp_mp) then
#ifndef TRPMIXP_OFF
c
         if (iam.eq.0.and.oprint) then
            write(luout,1850) nvpass
            call util_flush(luout)
         endif
 1850    format(' commencing triples evaluation - OpenMP',
     $          ' mixed-precision version',i8)
         call ccsd_trpdrv_omp_mixprec(dbl_mb(k_t1),
     $        eorb,g_objo,g_objv,g_coul,g_exch,ncor,nocc,nvir,iprt,
     $        empt(1),empt(2),oseg_lo,oseg_hi,kchunk,
     $        dbl_mb(k_trp_Tij), dbl_mb(k_trp_Tkj), dbl_mb(k_trp_Tia),
     $        dbl_mb(k_trp_Tka), dbl_mb(k_trp_Xia), dbl_mb(k_trp_Xka),
     $        dbl_mb(k_trp_Jia), dbl_mb(k_trp_Jka), dbl_mb(k_trp_Kia),
     $        dbl_mb(k_trp_Kka), dbl_mb(k_trp_Jij), dbl_mb(k_trp_Jkj),
     $        dbl_mb(k_trp_Kij), dbl_mb(k_trp_Kkj), dbl_mb(k_trp_Dja),
     $        dbl_mb(k_trp_Djka), dbl_mb(k_trp_Djia))
#else
         call errquit('aoccsd: trpdrv_omp_mixprec disabled ',0,0)
#endif
c
         elseif (use_trpdrv_nb) then
c
         if (iam.eq.0.and.oprint) then
            write(luout,1847) nvpass
            call util_flush(luout)
         endif
 1847    format(' commencing triples evaluation - non-blocking',i8)
         call ccsd_trpdrv_nb(dbl_mb(k_t1),
#ifdef USE_F90_ALLOCATABLE
     $        f1n,f1t,f2n,f2t,f3n,f3t,f4n,f4t,
#else
     $        dbl_mb(k_f1n),dbl_mb(k_f1t),dbl_mb(k_f2n),dbl_mb(k_f2t),
     $        dbl_mb(k_f3n),dbl_mb(k_f3t),dbl_mb(k_f4n),dbl_mb(k_f4t),
#endif
     $        eorb,g_objo,g_objv,g_coul,g_exch,ncor,nocc,nvir,iprt,
     $        empt(1),empt(2),oseg_lo,oseg_hi,kchunk,
     $        dbl_mb(k_trp_Tij), dbl_mb(k_trp_Tkj), dbl_mb(k_trp_Tia),
     $        dbl_mb(k_trp_Tka), dbl_mb(k_trp_Xia), dbl_mb(k_trp_Xka),
     $        dbl_mb(k_trp_Jia), dbl_mb(k_trp_Jka), dbl_mb(k_trp_Kia),
     $        dbl_mb(k_trp_Kka), dbl_mb(k_trp_Jij), dbl_mb(k_trp_Jkj),
     $        dbl_mb(k_trp_Kij), dbl_mb(k_trp_Kkj), dbl_mb(k_trp_Dja),
     $        dbl_mb(k_trp_Djka), dbl_mb(k_trp_Djia))
c
#ifdef USE_MIC_TRPDRV
         elseif (use_trpdrv_offload) then
c
         if (iam.eq.0.and.oprint) then
            write(luout,1849) nvpass
            call util_flush(luout)
         endif
 1849    format(' commencing triples evaluation - offload version',i8)
         call ccsd_trpdrv_offload(dbl_mb(k_t1),
#ifdef USE_F90_ALLOCATABLE
     $        f1n,f1t,f2n,f2t,f3n,f3t,f4n,f4t,
#else
     $        dbl_mb(k_f1n),dbl_mb(k_f1t),dbl_mb(k_f2n),dbl_mb(k_f2t),
     $        dbl_mb(k_f3n),dbl_mb(k_f3t),dbl_mb(k_f4n),dbl_mb(k_f4t),
#endif
     $        eorb,g_objo,g_objv,g_coul,g_exch,ncor,nocc,nvir,
     $        empt(1),empt(2),oseg_lo,oseg_hi,kchunk,
     $        dbl_mb(k_trp_Tij), dbl_mb(k_trp_Tkj), dbl_mb(k_trp_Tia),
     $        dbl_mb(k_trp_Tka), dbl_mb(k_trp_Xia), dbl_mb(k_trp_Xka),
     $        dbl_mb(k_trp_Jia), dbl_mb(k_trp_Jka), dbl_mb(k_trp_Kia),
     $        dbl_mb(k_trp_Kka), dbl_mb(k_trp_Jij), dbl_mb(k_trp_Jkj),
     $        dbl_mb(k_trp_Kij), dbl_mb(k_trp_Kkj), dbl_mb(k_trp_Dja),
     $        dbl_mb(k_trp_Djka), dbl_mb(k_trp_Djia))
c
#endif
#ifdef BGP
         elseif (use_trpdrv_bgp2) then
c
         if (iam.eq.0.and.oprint) then
            write(luout,1849)nvpass
            call util_flush(luout)
         endif
 1849    format(' commencing triples evaluation - Blue Gene v2',i8)
         call ccsd_trpdrv_bgp2(dbl_mb(k_t1),
#ifdef USE_F90_ALLOCATABLE
     $        f1n,f1t,f2n,f2t,f3n,f3t,f4n,f4t,
#else
     $        dbl_mb(k_f1n),dbl_mb(k_f1t),dbl_mb(k_f2n),dbl_mb(k_f2t),
     $        dbl_mb(k_f3n),dbl_mb(k_f3t),dbl_mb(k_f4n),dbl_mb(k_f4t),
#endif
     $        eorb,g_objo,g_objv,g_coul,g_exch,ncor,nocc,nvir,iprt,
     $        empt(1),empt(2),oseg_lo,oseg_hi,kchunk,
     $        dbl_mb(k_trp_Tij), dbl_mb(k_trp_Tkj), dbl_mb(k_trp_Tia),
     $        dbl_mb(k_trp_Tka), dbl_mb(k_trp_Xia), dbl_mb(k_trp_Xka),
     $        dbl_mb(k_trp_Jia), dbl_mb(k_trp_Jka), dbl_mb(k_trp_Kia),
     $        dbl_mb(k_trp_Kka), dbl_mb(k_trp_Jij), dbl_mb(k_trp_Jkj),
     $        dbl_mb(k_trp_Kij), dbl_mb(k_trp_Kkj), dbl_mb(k_trp_Dja),
     $        dbl_mb(k_trp_Djka), dbl_mb(k_trp_Djia))
c
#endif
         else
c
         if (iam.eq.0.and.oprint) then
            write(luout,1846)nvpass
            call util_flush(luout)
         endif
 1846    format(' commencing triples evaluation - blocking',i8)
         call ccsd_trpdrv(dbl_mb(k_t1),
#ifdef USE_F90_ALLOCATABLE
     $        f1n,f1t,f2n,f2t,f3n,f3t,f4n,f4t,
#else
     $        dbl_mb(k_f1n),dbl_mb(k_f1t),dbl_mb(k_f2n),dbl_mb(k_f2t),
     $        dbl_mb(k_f3n),dbl_mb(k_f3t),dbl_mb(k_f4n),dbl_mb(k_f4t),
#endif
     $        eorb,g_objo,g_objv,g_coul,g_exch,ncor,nocc,nvir,iprt,
     $        empt(1),empt(2),oseg_lo,oseg_hi,kchunk,
     $        dbl_mb(k_trp_Tij), dbl_mb(k_trp_Tkj), dbl_mb(k_trp_Tia), 
     $        dbl_mb(k_trp_Tka), dbl_mb(k_trp_Xia), dbl_mb(k_trp_Xka), 
     $        dbl_mb(k_trp_Jia), dbl_mb(k_trp_Jka), dbl_mb(k_trp_Kia), 
     $        dbl_mb(k_trp_Kka), dbl_mb(k_trp_Jij), dbl_mb(k_trp_Jkj), 
     $        dbl_mb(k_trp_Kij), dbl_mb(k_trp_Kkj), dbl_mb(k_trp_Dja),
     $        dbl_mb(k_trp_Djka), dbl_mb(k_trp_Djia))
c
         endif
c
         tx(3)=tcgtime()

         if (iam.eq.0.and.oprint)then
            write(luout,'(1x,a,i5,f12.2)')
     &       'Time for integral evaluation pass ',
     $       nvpass,tx(2)-tx(1)
            write(luout,'(1x,a,i5,f12.2)')
     &       'Time for triples evaluation pass  ',
     $      nvpass,tx(3)-tx(2)
            call util_flush(luout)
         endif

      enddo

      if (.not. ma_chop_stack(l_trp_Tij)) call errquit
     $     ('ccsd triples: stack corrupt (l_trp_Tij)',0, MA_ERR)

#ifdef USE_F90_ALLOCATABLE
      deallocate( f1n, stat=alloc_error)
      if (alloc_error.ne.0) call errquit('f1n',11,MA_ERR)
      deallocate( f2n, stat=alloc_error)
      if (alloc_error.ne.0) call errquit('f2n',12,MA_ERR)
      deallocate( f3n, stat=alloc_error)
      if (alloc_error.ne.0) call errquit('f3n',13,MA_ERR)
      deallocate( f4n, stat=alloc_error)
      if (alloc_error.ne.0) call errquit('f4n',14,MA_ERR)
      deallocate( f1t, stat=alloc_error)
      if (alloc_error.ne.0) call errquit('f1t',15,MA_ERR)
      deallocate( f2t, stat=alloc_error)
      if (alloc_error.ne.0) call errquit('f2t',16,MA_ERR)
      deallocate( f3t, stat=alloc_error)
      if (alloc_error.ne.0) call errquit('f3t',17,MA_ERR)
      deallocate( f4t, stat=alloc_error)
      if (alloc_error.ne.0) call errquit('f4t',18,MA_ERR)
#else
      if (.not. ma_chop_stack(l_f1n)) call errquit
     $     ('ccsd triples: stack corrupt (f1n)',0, MA_ERR)
#endif

      empt(1)=empt(1)+empt(1)
      empt(2)=empt(2)+empt(2)
      call ga_dgop(msg_cc_emp2,empt,2, '+')
      eccsdt=empt(1)+empt(2)
      if (iam.eq.0.and.oprint)then
       write(luout,100)empt(1),empt(2),eccsdt
 100   format(/,
     &     ' pseudo-e(mp4) ',e22.14,/,
     &     ' pseudo-e(mp5) ',e22.14,/,
     &     '        e(t)   ',e22.14,/)
       call util_flush(luout)
      endif

      stat=ga_destroy(g_coul)
      stat=stat.and.ga_destroy(g_exch)
c
      If ( IAm .eq. 0.and.oprint) then
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
         call util_flush(luout)
      EndIf
C
C     Write our energy out to the database
C
      If ( .NOT. RTDB_Put( RTDB, 'ccsd:ccsd+t(ccsd) correlation energy',
     $   MT_Dbl, 1, Eccsd+empt(1) ) )
     $   Call ErrQuit('CCSD: Can''t write to RTDB', 10,0)
      If ( .NOT. RTDB_Put( RTDB, 'ccsd:ccsd+t(ccsd) energy', MT_Dbl, 1,
     $   RefEner+eccsd+empt(1) ) ) 
     $   Call ErrQuit('CCSD: Can''t write to RTDB', 11,0)
      If ( .NOT. RTDB_Put( RTDB, 'ccsd:ccsd(t) correlation energy',
     $   MT_Dbl, 1, Eccsd+eccsdt) )
     $   Call ErrQuit('CCSD: Can''t write to RTDB', 12,0)
      If ( .NOT. RTDB_Put( RTDB, 'ccsd:ccsd(t) energy', MT_Dbl, 1,
     $   RefEner+eccsd+eccsdt) ) 
     $   Call ErrQuit('CCSD: Can''t write to RTDB', 13,0)
C
C     This reports back to task_energy that the requested task has been
C     completed successfully.
C
      If ( Inp_Compare(.FALSE., 'ccsd(t)', CC_Theory) ) then
            call ecce_print1('total energy', mt_dbl, 
     $        RefEner+eccsd+eccsdt, 1)
         If ( .NOT. RTDB_Put( RTDB, 'ccsd:energy', MT_Dbl, 1,
     $      RefEner+eccsd+eccsdt) ) 
     $      Call ErrQuit('CCSD: Can''t write to RTDB', 14,0)
      ElseIf ( Inp_Compare(.FALSE., 'ccsd+t(ccsd)', CC_Theory) ) then
         call ecce_print1('total energy', mt_dbl, 
     $        RefEner+eccsd+empt(1), 1)
         If ( .NOT. RTDB_Put( RTDB, 'ccsd:energy', MT_Dbl, 1,
     $      RefEner+eccsd+empt(1) ) ) 
     $      Call ErrQuit('CCSD: Can''t write to RTDB', 14,0)
      EndIf
C
C
C
      stat=.true.
      stat=stat.and.ma_pop_stack(l_t1)
      if (.not.stat)call errquit('trpdrv: cannot pop stack?',0, MA_ERR)
c
      if (.not. ga_destroy(g_objo))
     &     call errquit('g_objo ga_destroy failed',0, MA_ERR)
c
c -- end of triples calculation
c
c      call ga_sync()
c      if (iam.eq.0)then
c       print *,' ga allocated end of aoccsd '
c       call ga_summarize(.false.)
c      endif
      call ga_sync()
      if (occsdps) then
         call pstat_off(ps_tripls)
      else
         call qexit('triples',0)
      endif
c
 999  if (.not. ga_destroy(g_objv))
     &     call errquit('g_objv ga_destroy failed',0, GA_ERR)
c
      if(.not.occsdps) call qstat()
      return
      end
      logical function ccsd_createg(g, rows, cols, r1, r2, r3, mcm,
     O     oprint)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
c
      integer g
      integer rows, cols, r1, r2, r3 ,mcm
      logical oprint
c
      integer mapr(3),l_mapc,k_mapc,i,nproc,
     C     nblocks,blksize,nblocksr,
     R     l_mapr,k_mapr,j
      integer blksize_0,i_switch
c
      nproc=ga_nnodes()
c
      nblocks=min(nproc/3,cols)
      blksize_0=max(cols/nblocks,1)
c
      if(ga_nodeid().eq.0.and.oprint) then
         write(6,1) 'nprocs ',nproc
         write(6,1) 'rows cols nprocs ',rows,cols,nproc
         write(6,1) 'r1 r2 r2 mcm',r1,r2,r3,mcm
      endif
      nblocksr=max(r1/mcm+r2/mcm+r3/mcm,1)
 11    continue
      if(nblocks*nblocksr.gt.nproc) then
         nblocksr=nproc/nblocks
         mcm=(r1+r2+r3)/nblocksr
         if(ga_nodeid().eq.0.and.oprint) then
            write(6,*) ' nblocksr ',nblocksr,
     c           nblocks*nblocksr
            write(6,*) ' np/nb ',
     D        dble(nproc)/dble(nblocks),nproc,nblocks
            write(6,*) ' new nblocksr ',nblocksr
            write(6,*) ' mcm increased to ',mcm
            write(6,1) ' segs ',r1/mcm,r2/mcm,r3/mcm
            write(6,1) ' nsegs*nbc ',
     C           (r1/mcm+r2/mcm+r3/mcm)*nblocks
         endif
 12      nblocksr=max(r1/mcm,1)+max(r2/mcm,1)+max(r3/mcm,1)
         if(nblocks*nblocksr.gt.nproc) then
         if(ga_nodeid().eq.0.and.oprint) 
     W           write(6,1) ' ouch ',nblocksr,nblocks,nproc
            mcm=mcm*2
            goto 12
         endif
      endif
      mapr(1)=1
      mapr(2)=r1+1
      mapr(3)=r1+r2+1
      if(.not.ma_push_get(MT_int,nblocks,
     .        'ga block maps',l_mapc,k_mapc))
     .     call errquit('ccsd_createg: not enough local mem',0,0)
      if(.not.ma_push_get(MT_int,nblocksr,
     .        'ga block mapr',l_mapr,k_mapr))
     .     call errquit('ccsd_createg: not enough local mem',0,0)
      nblocksr=0
      do i=1,max(r1/mcm,1)
         int_mb(k_mapr+i-1)=1+(i-1)*mcm
         nblocksr=nblocksr+1
      enddo
       if(ga_nodeid().eq.0.and.oprint) then
          write(6,*) ' 1st seg nblocksr ',nblocksr
          write(6,*) ' last map ',int_mb(k_mapr+nblocksr-1)
       endif
      do i=1,max(r2/mcm,1)
         int_mb(k_mapr+nblocksr)=r1+(i-1)*mcm
         nblocksr=nblocksr+1
      enddo
      if(ga_nodeid().eq.0.and.oprint) then
         write(6,*) ' 2nd seg nblocksr ',nblocksr
         write(6,*) ' last map ',int_mb(k_mapr+nblocksr-1)
      endif
      do i=1,max(r3/mcm,1)
         int_mb(k_mapr+nblocksr)=r1+r2+(i-1)*mcm
         nblocksr=nblocksr+1
      enddo
      if(ga_nodeid().eq.0.and.oprint) then
         write(6,*) ' 3rd seg nblocksr ',nblocksr
         write(6,*) ' last map ',int_mb(k_mapr+nblocksr-1)
         do i=1,nblocksr
            write(6,'(" MM ",5I13)')  i,int_mb(k_mapr+i-1)
         enddo
      endif
#if 0
      blksize=blksize_0
      do i=1,nblocks
         int_mb(k_mapc+i-1)=(i-1)*blksize+1
      enddo
#else
c     finer load distribution since we might get very large blocks for last proc
c     with the previous distribution
      blksize=blksize_0+1
      i_switch=cols-blksize_0*nblocks
      int_mb(k_mapc)=1
      do i=2,nblocks
         if(i.ge.i_switch) blksize=blksize_0
         int_mb(k_mapc+i-1)=int_mb(k_mapc+i-2)+blksize
!         if(ga_nodeid().eq.0)write(6,'(5i6)') i,blksize,
!     I        int_mb(k_mapc+i-1),cols-((i-1)*blksize+1),nblocks-i
      enddo
#endif
      if(ga_nodeid().eq.0.and.oprint) then
         write(6,*) ' after maps nblocksr ',nblocksr
         write(6,1) 'rows cols nprocs ',rows,cols,nproc
         write(6,1) 'r1 r2 r2 ',r1,r2,r3
         write(6,1) 'mapr ',mapr
         write(6,1) ' blksize nblocks ',blksize,nblocks
         write(6,1) ' map1 maplast',int_mb(k_mapc),
     I        int_mb(k_mapc+nblocks-1)
     
 1     format(A,7I10)
       call util_flush(6)
      endif
      call ga_sync()
      ccsd_createg =  ga_create_irreg(MT_DBL, rows, cols, 'dunno',
     $     int_mb(k_mapr), nblocksr,
     $     int_mb(k_mapc), nblocks, g)

      if (.not.ma_chop_stack(l_mapc))
     &     call errquit('ccsd_createg: cannot chop stack',0,0)
      
      return
      end
