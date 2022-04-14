      subroutine ccsd_pampt3(basis,nsh,ncor,nocc,nvir,nact,nbf,
     &                       tol2e,cmo,eorb,iprt,hiu,giu,habe,gabe,
     &                       hia,t1,ht1,idiis,g_nt2,g_nht2,
     &                       g_ncoul,g_nexch,tklst,use_ccsd_omp)
C     $Id$
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
#include "schwarz.fh"
#include "ccsd_debug.fh"
#include "ccsdps.fh"
#include "ccsd_data.fh"
#include "util.fh"
C
      integer l_snsi,k_snsi,l_sisn,k_sisn,l_eri1,k_eri1,
     &        l_eri2,k_eri2,l_scre,k_scre,l_t1ao,k_t1ao,
     &        l_offsh,k_offsh,l_sc1,k_sc1,l_sc2,k_sc2,
     &        l_sc3,k_sc3,l_sc4,k_sc4,intblk,kscr,
     &        maxbfsh,max1e,mem1,max2e,mem2,lscr,lssni,
     &        g_x,g_c,g_st2,g_sht2,l_sch,k_sch,l_ish,k_ish,nrints,
     &        l_shinf,k_shinf,nsz_st2
#ifdef USE_F90INTERFACE
      interface
      subroutine ccsd_sxy(basis,nsh,ncor,nocc,nvir,nact,nbf,g_st2,
     &                    g_sht2,g_c,g_x,offsh,snsi,sisn,lssni,scre,
     &                    mem2,max2e,eri1,eri2,t1,cmo,t1ao,scra,
     &                    scrb,lscr,kscr,tol2e,iprt,tklst,maxints,
     &                    ish_idx,shinf,max_sht2_blk,nbfdim,a_st2,
     &                    use_ccsd_omp)
      integer basis,nsh,ncor,nocc,nvir,nact,nbf,g_st2,g_sht2,
     &        g_c,g_x,lssni,
     &        mem2,max2e,lscr,kscr,iprt,max_sht2_blk,nbfdim
      integer offsh(nsh,nsh,*),tklst(*),maxints
      integer ish_idx(*),shinf(nsh,*)
      double precision tol2e,eri1(*),eri2(*),scre(*),
     &                 t1(nocc*nvir),cmo(nbf,*),scra(*),
     &                 scrb(*),snsi(*),sisn(*),
     &     t1ao(*),a_st2(nocc*nocc,*)
      logical, optional, intent(in) :: use_ccsd_omp
      end subroutine
      end interface
#endif
      logical schwarz1,schwarz2
      logical stat
      integer max_sht2_blk
      integer ilo,ihi,jlo,jhi,nbfdim,ish,jsh,nfi,nfj
      integer mem_ava,mem_left
      integer k_ast2,l_ast2
      integer lnoo, i
      character*255 st2file
      double precision tstart, tused
      integer g_st2_rep,blk(2)
      logical util_mirrmat_chunk
      external util_mirrmat_chunk
      logical, optional, intent(in) :: use_ccsd_omp
      if (.not.present(use_ccsd_omp)) then
          call errquit('ccsd_pampt2: use_ccsd_omp not present!',0,0)
      endif
      Nodes = GA_NNodes()
      IAm = GA_NodeID()
      lnoo = nocc * nocc
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
      call ga_distribution(g_st2,
     .     ga_nodeid(), ilo, ihi, jlo, jhi)
c      write(6,'(i3,a,4i5)') ga_nodeid(), ' old distr ',
c     L           ilo, ihi, jlo, jhi
#ifdef HAVE_SET_GA_PROPERTY
      if(use_rdonly_st2) then
         call ga_set_property(g_st2,'read_only')
      call ga_distribution(g_st2,
     .     ga_nodeid(), ilo, ihi, jlo, jhi)
c      write(6,'(i3,a,4i5)') ga_nodeid(), ' new distr ',
c     L           ilo, ihi, jlo, jhi
         if(ga_nodeid().eq.0) write(6,*) ' g_st2 read_only'
c     disable old replicated data
         use_inmem_st2=.false.
      endif
#endif
      if(use_mirr_st2) then
         blk(1)=nocc*nocc
         blk(2)=0
         if(util_mirrmat_chunk(1,g_st2,g_st2_rep,.true.,.false.,
     B        blk)) then
            use_inmem_st2=.false.
            if (.not.ga_destroy(g_st2))
     &           call errquit('ga_dest g_st2 fail',1, GA_ERR)
            g_st2=g_st2_rep
            if(ga_nodeid().eq.0) write(6,*) ' g_st2 mirrored'
         else
            use_mirr_st2=.false.
         endif
      endif
      call ga_distribution(g_st2,
     .     ga_nodeid(), ilo, ihi, jlo, jhi)
c      write(6,'(i3,a,4i5)') ga_nodeid(), ' new distr ',
c     L           ilo, ihi, jlo, jhi
      if (use_inmem_st2) then
c     show how much memory g_st2 (a_st2) needs
         nsz_st2 = lnoo * nbfdim*ma_sizeof(MT_DBL,1,MT_BYTE)
         mem_ava=ma_inquire_avail(MT_BYTE)
         if(iam.eq.0) then
          if (util_print('ccsd extra', print_high)) then
            write(*,'(a,i8,a)') " g_st2 size:",
     &           max(nsz_st2/1024/1024,1), " MB",
     &           " mem. avail ",max(mem_ava/1024/1024,1)," MB"
          endif
         endif
         tstart = util_wallsec()
c       use node-replicated a_st2 instead of g_st2
c     lets grab 80% of available memory
         mem_left=mem_ava*80/100-nsz_st2
         if(mem_left.gt.0) then
            use_inmem_st2 = ma_push_get(MT_DBL,lnoo*nbfdim,
     A           'ast2',l_ast2,k_ast2)
         else
            if(iam.eq.0) write(*,'(a,I10)')
     &           " Not enough mem for st2. Increase MA by",
     M           -mem_left/1024/1024," MB"
            use_inmem_st2 = .false.
         endif
         if(use_inmem_st2) then
           if(iam.eq.0) then
            if (util_print('ccsd extra', print_high)) then
             write(*,'(a)') " Memory based method: ST2 is allocated"
            endif
           endif
c     
c     compute ST2 in parallel: expect network congestion
           call util_mygabcast(g_st2,lnoo,nbfdim,dbl_mb(k_ast2),lnoo)
         else
          if(iam.eq.0) then
           write(*,'(a)') " Not Enough Memory to keep ST2 in local",
     &                    " memory - expect network congestion."
          endif
         endif

         
         tused = util_wallsec() - tstart
         if (iam.eq.0.and.use_inmem_st2) then
           if (util_print('ccsd extra', print_high)) then
             write(6,'(a,f9.2,a)') " ST2 array is replicated ",tused,"s"
           endif
           call util_flush(6)
         endif
         
      endif

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
      intblk=(ma_inquire_avail(MT_DBL)-2*ma_sizeof_overhead(MT_DBL))/
     &       (2*max2e)
c     leave 10 pct available
      intblk=intblk*9/10
      intblk=max(1,intblk)
      nrints=max2e*intblk
      stat=stat.and.ma_push_get(MT_DBL,nrints,'eri1',l_eri1,k_eri1)
      stat=stat.and.ma_push_get(MT_DBL,nrints,'eri2',l_eri2,k_eri2)
      if (.not.stat)
     &   call errquit('ccsd_pampt3: ma_push_get ccsd_sxy ',0, MA_ERR)
      call ccsd_sxy(basis,nsh,ncor,nocc,nvir,nact,nbf,g_st2,g_sht2,
     &              g_c,g_x,int_mb(k_offsh),dbl_mb(k_snsi),
     &              dbl_mb(k_sisn),lssni,dbl_mb(k_scre),
     &              mem2,max2e,dbl_mb(k_eri1),
     &              dbl_mb(k_eri2),t1,cmo,dbl_mb(k_t1ao),dbl_mb(k_sc1),
     &              dbl_mb(k_sc2),lscr,kscr,tol2e,iprt,tklst,nrints,
     &              int_mb(k_ish),int_mb(k_shinf),max_sht2_blk,
     &              nbfdim,dbl_mb(k_ast2),use_ccsd_omp)
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

#ifdef HAVE_SET_GA_PROPERTY
      if(use_rdonly_st2) call ga_unset_property(g_st2)
#endif
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
      if (use_inmem_st2 ) then
      if(.not.ma_pop_stack(l_ast2)) then
         call ma_summarize_allocated_blocks()
         call errquit('ccsd_pampt3: ma pop stack fails ',11,MA_ERR)
         endif
      endif
c ------------------------------------------------------------
c other terms using B intermediate
c ------------------------------------------------------------
      lscr=nact*nact
      stat=stat.and.ma_push_get(MT_DBL,lscr,'scra4',l_sc1, k_sc1)
      stat=stat.and.ma_push_get(MT_DBL,lscr,'scrb',l_sc2, k_sc2)
      stat=stat.and.ma_push_get(MT_DBL,lscr,'scre',l_sc3, k_sc3)
      stat=stat.and.ma_push_get(MT_DBL,lscr,'scrf',l_sc4, k_sc4)
!
! WARNING: dbl_mb(k_sc[12]) are aliased below, which violates Fortran...
!
      if (use_ccsd_omp) then
      call ccsd_itm_omp(nocc,nvir,nact,hiu,giu,habe,gabe,hia,t1,ht1,
     &                  g_nt2,g_nht2,g_ncoul,g_nexch,g_x,g_c,g_sht2,
     &                  iprt)
      else
      call ccsd_itm(nocc,nvir,nact,hiu,giu,habe,gabe,hia,t1,ht1,
     &              dbl_mb(k_sc1),dbl_mb(k_sc2),
     &              dbl_mb(k_sc1),dbl_mb(k_sc2),
     &              dbl_mb(k_sc3),dbl_mb(k_sc4),
     &              g_nt2,g_nht2,g_ncoul,g_nexch,g_x,g_c,g_sht2,
     &              iprt)
      endif
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
      subroutine t2eri(idim,idx,jlo,jhi,nfj,llo,lhi,nfl,nsh,eri1,eri2,
     &                 t2,ht2,lnoo,nocc,offsh,nbf,g_st2,
     &                 shinf,max_st2_blk,snsi,sisn,t1ao,nbfdim,a_st2,
     &                 use_ccsd_omp)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "bas.fh"
#include "geom.fh"
#include "schwarz.fh"
#include "eaf.fh"
#include "ccsdps.fh"
#include "ccsd_data.fh"
#define NBGET 1
c
#ifdef CRAYXT
c#warning undefining NBGET
#undef NBGET
#endif
      integer idim,jlo,jhi,nsh,llo,lhi,lnoo,nocc,nbfdim
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
     &                 t1ao(nbf,nocc), a_st2(lnoo,nbfdim)
      logical use_ccsd_omp
c
      integer i,j,k,l,kmax,lmax,ik,ad1,ptr
      integer ipp,imm,ict,off1,off2,off3
      integer mlkj,mkl,mkj,ii,jj
      integer blk1,blk2,ikin,ik0
#ifdef NBGET
      integer kmsh,imsh,nmfi,nmfk,nbhandle
      integer kmlo,kmhi,imlo,imhi
#endif
      double precision gp,gm,factor
      logical keqi,jeql
      integer szt
      szt = 2 * max_st2_blk - 1
c
      if (occsdps) then
         call pstat_on(ps_t2eri)
      else
         call qenter('t2eri',0)
      endif
c
      jeql=jlo.eq.llo
      ipp=nocc*(nocc+1)/2
      imm=nocc*(nocc-1)/2
      ict=1
cstag
c     figure out initial ict
      ikin=ga_nodeid()
      ik0=mod(ikin,(idim-3+1))+3
      if(ik0.gt.3) then
         do ik=3,ik0-1
            ksh=(idx(ik)-1)/nsh+1
            ish=idx(ik)-(ksh-1)*nsh
            nfk=shinf(ksh,1)
            nfi=shinf(ish,1)
            ict=ict+nfi*nfk*nfl*nfj
      if(nfi.lt.0) then
         write(6,1964) ga_nodeid(),' : nfi ',nfi, ' ish ',ish,
     A        ' nsh ',nsh, ' ik ',ik
 1964    format(I5,A,I24,A,I10,A,I10,A,I10)
         call util_flush(6)
         call errquit(' t2eri: bug1 ',0,0)
      endif
      if(nfk.lt.0) then
         write(6,*) ga_nodeid(),' : nfk ',nfk, ' ksh ',ksh,' nsh ',nsh
         call util_flush(6)
         call errquit(' t2eri: bug2 ',0,0)
      endif
         enddo
      endif
c     
c Get first block ik=ik0
c
      ksh=(idx(ik0)-1)/nsh+1
      ish=idx(ik0)-(ksh-1)*nsh
      nfk=shinf(ksh,1)
      klo=shinf(ksh,2)
      khi=shinf(ksh,3)
      nfi=shinf(ish,1)
      ilo=shinf(ish,2)
      ihi=shinf(ish,3)
      ad1=offsh(ish,ksh,2)
      ptr = ad1 + 1
      blk1=0
      blk2=max_st2_blk
#ifdef NBGET
      nbhandle=4321+ga_nodeid()
      if (.not.use_inmem_st2) then
          call ga_get(g_st2,1,lnoo,ad1+1,ad1+nfi*nfk,t2(0),lnoo)
      endif
c
c Now loop over remaining blocks
c
cold       do ik=4,idim
         do ikin=ga_nodeid()+1,(idim-3)+ga_nodeid()
            ik=mod(ikin,(idim-3+1))+3
            if(ik.eq.4) ict=1
#else
cold       do ik=3,idim
         do ikin=ga_nodeid(),(idim-3)+ga_nodeid()
            ik=mod(ikin,(idim-3+1))+3
            if(ik.eq.3) ict=1
#endif
#ifdef NBGET
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
        if (.not.use_inmem_st2) then
            call ga_nbget(g_st2,1,lnoo,ad1+1,ad1+nmfi*nmfk,t2(blk2),
     &                lnoo,nbhandle)
        endif
#else
        ksh=(idx(ik)-1)/nsh+1
        ish=idx(ik)-(ksh-1)*nsh
        nfk=shinf(ksh,1)
        klo=shinf(ksh,2)
        khi=shinf(ksh,3)
        nfi=shinf(ish,1)
        ilo=shinf(ish,2)
        ihi=shinf(ish,3)
        ad1=offsh(ish,ksh,2)
        if (use_inmem_st2) then
          ptr = ad1 + 1
        else
            call ga_get(g_st2,1,lnoo,ad1+1,ad1+nfi*nfk,t2(0),lnoo)
        endif
#endif
c       call pstat_off(ps_t2erin)
c
c  Do the actual work
c
        keqi=ish.eq.ksh
        if (keqi) then
          factor = 0.5d0
        else
          factor = 1.0d0
        endif
        call ccsd_idx1_wrk(snsi,sisn,eri1(ict),eri2(ict),t1ao,nfj,
     &                     nfl,nocc,nbf,ilo,ihi,klo,khi,factor,
     &                     use_ccsd_omp)
        if (use_inmem_st2) then
           call t2eri_wrk_all(eri1(ict),eri2(ict),a_st2(1,ptr),ht2,
     &                        nfi,nfj,nfk,nfl,nocc,lnoo,keqi,jeql,
     &                        use_ccsd_omp)
        else
           call t2eri_wrk_all(eri1(ict),eri2(ict),t2(blk1),ht2,
     &                        nfi,nfj,nfk,nfl,nocc,lnoo,keqi,jeql,
     &                        use_ccsd_omp)
        endif
c
c Increase offset of integral block to next ik block
c
        ict=ict+nfi*nfj*nfk*nfl
c
c Wait for non-blocking get to finish and replace blk with correct pointer
c
c       call pstat_on(ps_t2erih)
#ifdef NBGET
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
        if (use_inmem_st2) then
          ptr = ad1 + 1
        else
          call ga_nbwait(nbhandle)
        endif
        blk1=-blk1+max_st2_blk
        blk2=-blk2+max_st2_blk
#else
        blk1=0
        blk2=0
#endif
c       call pstat_off(ps_t2erih)
      enddo
c
c  Do the actual work on final block
c
#ifdef NBGET
      keqi=ish.eq.ksh
      if(ik.eq.3) ict=1
      if (keqi) then
         factor = 0.5d0
      else
         factor = 1.0d0
      endif
      call ccsd_idx1_wrk(snsi,sisn,eri1(ict),eri2(ict),t1ao,nfj,
     &                   nfl,nocc,nbf,ilo,ihi,klo,khi,factor,
     &                   use_ccsd_omp)
      if (use_inmem_st2) then
         call t2eri_wrk_all(eri1(ict),eri2(ict),a_st2(1,ptr),ht2,
     &                      nfi,nfj,nfk,nfl,nocc,lnoo,keqi,jeql,
     &                      use_ccsd_omp)
      else
         call t2eri_wrk_all(eri1(ict),eri2(ict),t2(blk1),ht2,
     &                      nfi,nfj,nfk,nfl,nocc,lnoo,keqi,jeql,
     &                      use_ccsd_omp)
      endif
#endif
c
      if (occsdps) then
         call pstat_off(ps_t2eri)
      else
         call qexit('t2eri',0)
      endif
c
      return
      end


      subroutine t2eri_wrk_all(eri1,eri2,t2,ht2,
     &                         ni,nj,nk,nl,nocc,lnoo,
     &                         keqi,jeql,use_ccsd_omp)
      implicit none
      integer, intent(in)             :: ni,nj,nk,nl,nocc,lnoo
      double precision, intent(in)    :: eri1(nl,nk,nj,ni)
      double precision, intent(in)    :: eri2(nj,nk,nl,ni)
      double precision, intent(in)    :: t2(lnoo,nk,ni)
      double precision, intent(inout) :: ht2(lnoo,nl,nj)
      logical, intent(in)             :: keqi,jeql,use_ccsd_omp
      double precision                :: gp, gm
      integer                         :: ipp,imm
      integer                         :: i,j,k,l,q
      integer                         :: kmax, lmax
      ipp=nocc*(nocc+1)/2
      imm=nocc*(nocc-1)/2
      kmax = nk
      lmax = nl
      if (use_ccsd_omp) then
!
! WARNING: You cannot use collapse on this loop because that will
!          cause concurrent updates to ht2, which is a race-condition
!          and leads to incorrectness.  In short, the J loop must be
!          the one and only that is threaded.  Okay, we can thread over
!          the L loop as well, but that may screw up memory access.
!
!          It may be appropriate to block these loops in any case.
!
!$omp parallel do schedule(static)
!$omp& default(none)
!$omp& private(gp,gm)
!$omp& firstprivate(kmax,lmax)
!$omp& shared(ni,nj,nk,nl,ipp,imm,keqi,jeql)
!$omp& shared(eri1,eri2,t2,ht2)
        do j=1,nj
          do i=1,ni
            ! branch predictor should do a good job here since it will
            ! be constant for any given invocation of the subroutine.
            if (keqi) kmax = i
            if (jeql) lmax = j
            do k=1,kmax
              do l=1,lmax
                gp=(eri1(l,k,j,i)+eri2(j,k,l,i))*0.5d0
                gm=(eri1(l,k,j,i)-eri2(j,k,l,i))*0.5d0
                do q=1,ipp
                  ht2(q,l,j) = ht2(q,l,j) + gp * t2(q,k,i)
                enddo
                do q=1,imm
                  ht2(q+ipp,l,j) = ht2(q+ipp,l,j) + gm * t2(q+ipp,k,i)
                enddo
              enddo
            enddo
          enddo
        enddo
!$omp end parallel do
      else
        do j=1,nj
          do i=1,ni
            ! branch predictor should do a good job here since it will
            ! be constant for any given invocation of the subroutine.
            if (keqi) kmax = i
            if (jeql) lmax = j
            do k=1,kmax
              do l=1,lmax
                gp=(eri1(l,k,j,i)+eri2(j,k,l,i))*0.5d0
                gm=(eri1(l,k,j,i)-eri2(j,k,l,i))*0.5d0
                do q=1,ipp
                  ht2(q,l,j) = ht2(q,l,j) + gp * t2(q,k,i)
                enddo
                do q=1,imm
                  ht2(q+ipp,l,j) = ht2(q+ipp,l,j) + gm * t2(q+ipp,k,i)
                enddo
              enddo
            enddo
          enddo
        enddo
      endif
      return
      end
