      subroutine ccsd_iterdrv2(rtdb,basis,nsh,ncor,nocc,nvir,nact,nbf,
     &                        maxit,Thre,iprt,cmo,eorb,t1,
     &                        ht1,scra,scrb,hiu,hia,habe,giu,gabe,
     &                        bbkp,bb,g_t2,g_ht2,
     &                        mxvec,eccsd,max2e,mem2,
     &                        g_ncoul,g_nexch,tklst,offt2,nvp1,
     &                        refener, Tol2E, occd, oconverged,
     &                        eccsdos, use_ccsd_omp)
C     $Id$
      implicit none
#include "errquit.fh"
#include "util.fh"
#include "stdio.fh"
c
      integer rtdb,basis,nsh,ncor,nocc,nvir,nact,nbf,maxit,iprt,
     & g_t2,g_ht2,g_nt2,g_nht2,mxvec,g_ncoul,g_nexch,
     &        tklst(nsh*(nsh+1)/2,2)
      Double Precision Thre, Tol2E
      double precision time_iter_start
      integer time_prev_iter
      logical oconverged        ! Return convergence status
      logical occd
      logical, optional, intent(in) :: use_ccsd_omp
c
       double precision cmo(nbf,nbf),eorb(nbf),t1(*),
     &                  ht1(*),scra(*),scrb(*),
     &                  hiu(*),hia(*),habe(*),giu(*),gabe(*),
     &                  bbkp(*),bb(*),eccsd, refener, eccsdos
c
c$$$      double precision rtc
c$$$      external rtc
c
#include "ccsd_len.fh"
*rak-s
#include "ccsd_time.fh"
*rak-e
c
      integer iter,idiis,ndiis,bev,jhi,jlo,ihi,ilo,offt2,nvp1,
     &        i,av,lhoo,lhoov
      integer max2e,mem2
      double precision eold,dabs,tx(4),rms,deccsd,
     &                 tstart,tend

      Integer IAm
      logical orestart, oprint
      character*255 t2file, filename

      logical ccsd_read_restart, ccsd_write_restart
      external ccsd_read_restart, ccsd_write_restart


C
#include "msgids.fh"
#include "tcgmsg.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "ccsdps.fh"
#include "ccsd_data.fh"
#include "chemio.h"
#include "eaf.fh"
#ifdef USE_F90INTERFACE
      interface
      subroutine ccsd_pampt2(basis,idiis,nsh,ncor,nocc,nvir,nact,nbf,
     $                      iprt,cmo,eorb,t1,ht1,
     $                      scra,scrb,hiu,hia,habe,giu,gabe,
     $                      g_t2,g_ht2,g_nt2,g_nht2,
     $                      eccsd,g_ncoul,g_nexch,
     &        tklst, Tol2E, occd, eccsdos,use_ccsd_omp)
      integer basis,idiis,nsh,ncor,nocc,nvir,nact,nbf,iprt,
     &        g_t2,g_ht2,g_ncoul,g_nexch,g_nt2,g_nht2,
     &        tklst(*)
      double precision cmo(nbf,*),eorb(*),t1(*),
     &                 ht1(*),scra(*),scrb(*),
     &                 hiu(*),hia(*),habe(*),giu(*),gabe(*)
      double precision tol2e
      double precision eccsd,fac, eccsdos
      logical occd
      logical, optional, intent(in) :: use_ccsd_omp
      end subroutine
      end interface
#endif         
c
      double precision zip
      data zip/0.0d00/
c
      double precision util_scratch_dir_avail_for_me
      external util_scratch_dir_avail_for_me
      if (.not.present(use_ccsd_omp)) then
          call errquit('ccsd_iterdrv2: use_ccsd_omp not present!',0,0)
      endif
c
      IAm = GA_NodeID()
c
      oprint = util_print('information', print_low)
C
      call ga_sync()
      if (occsdps) then
         call pstat_on(ps_itrdrv)
      else
         call qenter('iterdrv',0)
      endif
c
      lhoo=nocc*(nocc+1)/2
      lhoov=nvir*nocc*(nocc+1)/2
c
      nvp1=nvir+1

c - integrals for ccsd_sxy (calculate as much as possible once and read after that)
c - input switch to avoid using disk
      if (.not. rtdb_get(rtdb, 'ccsd:usedisk', mt_log, 1, use_disk))
#ifdef NOIO
c     direct by default to avoid io
     $   use_disk=.false.
#else
     $   use_disk=.true.
#endif
      if (use_disk) then
         store_avail=util_scratch_dir_avail_for_me()
c - convert back to bytes
         store_avail=store_avail*1024
         call util_file_name('ccsdints',.true.,.true.,filename)
         if (eaf_open(filename,EAF_RW,sxy_hl).ne.0) call
     $      errquit('ccsd_iterdrv2: failed to open file',0,DISK_ERR)
      endif
      repeat=.false.

c - starting vectors

      call ga_distribution(g_t2,iam,jlo,jhi,ilo,ihi)

      call util_file_name('t2', .false.,.false.,t2file)

c - if converged, read vectors and return
c - else check if one can restart
      if (.not.rtdb_get(rtdb, 'ccsd:restart',MT_LOG, 1, orestart))
     $    orestart = .true.   ! By default always try to restart
      if (.not.rtdb_get(rtdb, 'ccsd:converged',MT_LOG, 1, oconverged))
     $    oconverged = .false.  ! Not converged yet
      if(.not.rtdb_get(rtdb,'ccsd:useinmemst2',mt_log, 1,use_inmem_st2))
     $   use_inmem_st2=.true. ! Old default kept ST2 in GA
      if(.not.rtdb_get(rtdb,'ccsd:rdonly_st2',mt_log, 1,use_rdonly_st2))
     $   use_rdonly_st2=.false. ! use ga_set_property read_only
      if(.not.rtdb_get(rtdb,'ccsd:mirr_st2',mt_log, 1,use_mirr_st2))
     $   use_mirr_st2=.false. ! use ga_set_property read_only

      if (orestart.or.oconverged) then
          orestart = ccsd_read_restart(t2file, g_t2, nvir, lhoov)
          if (oconverged) oconverged = orestart
      endif

c Note that if restart fails orestart will be set false
      if (.not. orestart) then
         if (iam .eq. 0 .and. oprint) write(luout,*) 
     M        ' Using MP2 initial guess vector '
        call dfill(lnov,zip,t1,1)
        if (nvp1.ge.ilo.and.nvp1.le.ihi)then
         call ga_put(g_t2,1,lnov,nvp1,nvp1,t1,lnov)
        endif
        call ccsd_mkt2(ncor,nocc,nvir,nact,eorb,scra,
     &                 scrb,g_nexch,g_t2,iprt)
      endif
      call ccsd_mktask(basis,nsh,tklst)
c
c ------------------------------------------------------------
c if the calculation is already converged get out of here
c ------------------------------------------------------------
c
      if(oconverged) then
         if (iam.eq.0 .and. oprint)write(luout,*) 
     $      'CCSD is already converged'
         offt2 = 0
         if(.not.rtdb_get(rtdb, 'ccsd:corlenerg',MT_DBL, 1, eccsd))
     $      eccsd=zip
         if(.not.rtdb_get(rtdb, 'ccsd:corlenergos',MT_DBL, 1, eccsdos))
     $      eccsdos=zip
         call ga_sync()
         if (occsdps) then
            call pstat_off(ps_itrdrv)
         else
            call qexit('iterdrv',0)
         endif
         return
      endif
c
c ------------------------------------------------------------
c create square GAs with proposed final ordering
c ------------------------------------------------------------
c
      if (.not.ga_create(MT_DBL,lnov,lnov,'nt2',
     &                   nvir,nvir,g_nt2))
     &     call errquit('ga_create g_nt2 failed',0, GA_ERR)
      if (.not.ga_create(MT_DBL,lnov,lnov,'nht2',
     &                   nvir,nvir,g_nht2))
     &     call errquit('ga_create g_nht2 failed',0, GA_ERR)
c
      call ga_sync()
c
      iter=0
      idiis=0
      ndiis=0
      eccsd=zip
      eold=zip
      rms=thre+thre
      deccsd=thre+thre
c
c ----------------------------------------------------------
c begin iterations
c ----------------------------------------------------------
      if (iam.eq.0 .and. oprint)write(luout,1234)
*rak 1234 format(//,71(1h-),/,' iter',8x,'correlation',11x,'delta',9x,
*rak     &       'rms',/,15x,'energy',14x,'energy',/,71(1h-))
 1234 format(//,73('-'),/,' iter',5x,'correlation',5x,'delta',7x,
     &       'rms',7x,'T2',5x,'Non-T2',6x,'Main',
     &    /,13x,'energy',6x,'energy',6x,'error',
     &    6x,'ampl',5x,'ampl',6x,'Block',
     &    /,48x,'time',5x,'time',6x,'time',
     &    /,73('-'))
   99 continue
      time_iter_start = util_wallsec()
      iter=iter+1
      idiis=idiis+1
      if (idiis.gt.mxvec)idiis=1
      ndiis=ndiis+1
      if (ndiis.gt.mxvec)ndiis=mxvec
      offt2=(idiis-1)*lhoov

c
c     Write out the current T1 and T2 for possible restart
c     (note that T1 is kept at end of T2)
c
      if (.not. ccsd_write_restart(t2file, g_t2, nvir, lhoov, offt2))
     $     call errquit('ccsd: failed saving vector for restart',0,
     &       UNKNOWN_ERR)
c
      oconverged = .false.
      if (dabs(rms).lt.thre.and.dabs(deccsd).lt.thre)then
         oconverged = .true.
         if(.not.rtdb_put(rtdb, 'ccsd:corlenerg',MT_DBL, 1, eccsd))
     $     call errquit('ccsd_iterdrv2: rtdb_put failed', 0, RTDB_ERR)
         if(.not.rtdb_put(rtdb, 'ccsd:corlenergos',MT_DBL, 1, eccsdos))
     $     call errquit('ccsd_iterdrv2: rtdb_put failed', 3, RTDB_ERR)
         if(.not.rtdb_put(rtdb, 'ccsd:converged',MT_LOG, 1, oconverged))
     $     call errquit('ccsd_iterdrv2: rtdb_put failed', 0, RTDB_ERR)
        if (iam.eq.0 .and. oprint)write(luout,1236)
 1236   format(18x,'*************converged*************',18x,
     &         /,73('-'),/)
        goto 999 
      endif
      if (iter.gt.maxit)then
        if (iam.eq.0 .and. oprint)write(luout,1237)
 1237   format(18x,'****maximum iterations exceeded****',18x,
     &         /,73('-'),/)
        goto 999
      endif
      if (iprt.gt.5.and.iam.eq.0)write(luout,134)iter,idiis,ndiis
 134  format(' iteration ',i3,' idiis ',i3,' ndiis ',i3)

*
*     RJH ... hack for CCD
      if (occd) then
      if (iam.eq.0)then
         call dfill(lnov, 0.0d0, t1, 1)
         call ga_put(g_t2,offt2+1,offt2+lnov,nvp1,nvp1,t1,lnov)
      endif
      endif
*
*

      if (iam.eq.0)then
       call ga_get(g_t2,offt2+1,offt2+lnov,nvp1,nvp1,t1,lnov)
      endif
c
      if (iam.eq.0.and.iprt.gt.5)then
       print *,'starting t1 vector'
       write(luout,139)(t1(i),i=1,lnov)
  139   format(4e14.4)
       do av=1,nvir
        call ga_get(g_t2,offt2+1,offt2+lhoov,av,av,scra,lhoov)
        print *,'starting t2 for av ',av
        write(luout,139)(scra(i),i=1,lhoov)
       enddo
      endif
c
      call ga_brdcst(msg_cc_t1a,t1,mdtob(lnov),0)
      tx(3)=tcgtime()
c
c ----------------------------------------------------------
c form ``sigma'' vector
c ----------------------------------------------------------
      call dfill(lnov,zip,ht1,1)
      call dfill(lhoov,zip,scra,1)
      if(g_ht2.eq.123456789) then
c     ga_created moved here from aoccsd
      if (.not.ga_create(MT_DBL,mxvec*lhoov,nvir+1,'g_ht2', 
     &                   mxvec*lhoov,0,g_ht2))
     &     call errquit('ga_create failed for g_ht2',0, GA_ERR)
      endif
      call ga_distribution(g_ht2,iam,jlo,jhi,ilo,ihi)
      call ga_sync
      do bev=1,nvir+1
         if (bev.ge.ilo.and.bev.le.ihi)then
            call ga_put(g_ht2,offt2+1,offt2+lhoov,bev,bev,scra,lhoov)
         endif
      enddo
      if (nvp1.ge.ilo.and.nvp1.le.ihi)then
         call ga_put(g_ht2,offt2+1,offt2+lnov,nvp1,nvp1,scra,lnov)
      endif
      call ga_sync()
      tx(1)=tcgtime()
      call ccsd_pampt2(basis,idiis,nsh,ncor,nocc,nvir,nact,nbf,
     &                iprt,cmo,eorb,t1,ht1,
     $                scra,scrb,hiu,hia,habe,giu,gabe,
     &                g_t2,g_ht2,g_nt2,g_nht2,eccsd,
     &                g_ncoul,g_nexch,tklst, Tol2E, occd, eccsdos,
     &                use_ccsd_omp)
      call ga_sync()
*
*     RJH more CCD hacks
*
      if (occd) then
      if (iam .eq. 0) then
         call dfill(lnov,zip,ht1,1)
         call ga_put(g_ht2,offt2+1,offt2+lnov,nvp1,nvp1,t1,lnov)
      endif
      call ga_sync()
      endif


      tx(2)=tcgtime()
      tstart=tx(2)
      call ccsd_pdiis(nocc,nvir,idiis,ndiis,t1,ht1,scra,scrb,mxvec,
     &                 bbkp,bb,iprt,iter,rms,g_t2,g_ht2,lhoov)
      call ga_sync()
      deccsd=eccsd-eold
      eold=eccsd
      tx(4)=tcgtime()
      tend=tx(4)
      if (iprt.gt.5.and.iam.eq.0)print *,' time pdiis ',tend-tstart
      if (iam.eq.0 .and. oprint)then
          write(luout,1235)iter,eccsd,deccsd,rms,
     &      tx(2)-tx(1),
     &      tx(4)-tx(3)-(tx(2)-tx(1)),
     &      main_block_time
      call util_flush(luout)
 1235 format(i4,f18.10,2(1p,d11.3,0p),3(f9.2))
c
      call ecce_print1('ccsd total energy', mt_dbl, eccsd+refener, 1)
      call ecce_print1('ccsd delta energy', mt_dbl, deccsd, 1)
      call ecce_print1('ccsd rms error', mt_dbl, rms, 1)
c
      endif
c
      time_prev_iter = nint((util_wallsec() - time_iter_start)*1.2d0)
c
      if (util_test_time_remaining(rtdb,time_prev_iter)) then
        goto 99
      else
        if (.not. ccsd_write_restart(t2file, g_t2, nvir, lhoov, offt2))
     $       call errquit('ccsd: failed saving vector for restart',1,
     &         UNKNOWN_ERR)
        goto 999
      endif
c
 999  continue
c
      if (.not.ga_destroy(g_nht2))
     &    call errquit('ga_dest g_nht2 fail',0, GA_ERR)
      if (.not.ga_destroy(g_nt2))
     &    call errquit('ga_dest g_nt2 fail',0, GA_ERR)
c
      if (use_disk) then
         if (eaf_close(sxy_hl).ne.0) call
     $   errquit('ccsd_iterdrv2: failed to close file',0,DISK_ERR)
      endif
c
      call ga_sync()
      if (occsdps) then
         call pstat_off(ps_itrdrv)
      else
         call qexit('iterdrv',0)
      endif
      return
      end
      logical function ccsd_read_restart(t2file, g_t2, nvir, lhoov)
      implicit none
#include "global.fh"
#include "inp.fh"
#include "util.fh"
#include "stdio.fh"
      character*(*) t2file
      integer nvir, lhoov
      integer g_t2
      logical oprint
c
      logical file_read_ga_patch
      external file_read_ga_patch
c
      call ga_sync
c
      oprint = util_print('information', print_low)
c
      ccsd_read_restart = 
     $     file_read_ga_patch(t2file,g_t2,1,lhoov,1,nvir+1)

      if (ga_nodeid().eq. 0 .and. oprint) then
         if (ccsd_read_restart) then
            write(luout,*) ' Restarting from old vector in ',
     $           t2file(1:inp_strlen(t2file))
         else
            write(luout,*) ' Failed reading restart vector from ',
     $           t2file(1:inp_strlen(t2file))
         endif
      endif
c
      end
      logical function ccsd_write_restart(t2file, g_t2, nvir, lhoov, 
     $     offt2)
      implicit none
#include "global.fh"
#include "util.fh"
#include "stdio.fh"
      character*(*) t2file
      integer nvir, lhoov, offt2
      integer g_t2
c
      logical file_write_ga_patch
      external file_write_ga_patch
c
c     file_write_ga_patch does a ga_sync immediately upon entry, 
c     so this call is unnecessary
c
c      call ga_sync
c
      ccsd_write_restart = 
     $     file_write_ga_patch(t2file,g_t2,offt2+1,offt2+lhoov,1,nvir+1)

      if ((.not. ccsd_write_restart) .and. ga_nodeid().eq. 0
     $     .and. util_print('information', print_low)) 
     $     write(luout,*) ' Write of T2 failed'
c
      end

         
