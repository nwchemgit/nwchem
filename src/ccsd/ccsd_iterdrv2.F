      subroutine ccsd_iterdrv2(rtdb,basis,nsh,ncor,nocc,nvir,nact,nbf,
     &                        maxit,Thre,iprt,cmo,eorb,t1,
     &                        ht1,scra,scrb,hiu,hia,habe,giu,gabe,
     &                        bbkp,bb,g_t2,g_ht2,
     &                        mxvec,eccsd,max2e,mem2,
     &                        g_ncoul,g_nexch,tklst,offt2,nvp1,
     $     refener, Tol2E, occd, oconverged)
C     $Id: ccsd_iterdrv2.F,v 2.18 2000-10-23 23:37:44 windus Exp $
      implicit none
#include "util.fh"
c
      integer rtdb,basis,nsh,ncor,nocc,nvir,nact,nbf,maxit,iprt,
     & g_t2,g_ht2,g_nt2,g_nht2,mxvec,g_ncoul,g_nexch,
     &        tklst(nsh*(nsh+1)/2,2)
      Double Precision Thre, Tol2E
      double precision time_iter_start
      integer time_prev_iter
      logical oconverged        ! Return convergence status
      logical occd
c
       double precision cmo(nbf,nbf),eorb(nbf),t1(*),
     &                  ht1(*),scra(*),scrb(*),
     &                  hiu(*),hia(*),habe(*),giu(*),gabe(*),
     &                  bbkp(*),bb(*),eccsd, refener
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
      logical orestart
      character*255 t2file

      logical ccsd_read_restart, ccsd_write_restart
      external ccsd_read_restart, ccsd_write_restart


C
#include "msgids.fh"
#include "tcgmsg.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
c
      double precision zip
      data zip/0.0d00/
c
      IAm = GA_NodeID()
C
      call ga_sync()
      call qenter('iterdrv',0)
c
      lhoo=nocc*(nocc+1)/2
      lhoov=nvir*nocc*(nocc+1)/2
c
      nvp1=nvir+1
c - starting vectors

      call ga_distribution(g_t2,iam,jlo,jhi,ilo,ihi)

      call util_file_name('t2', .false.,.false.,t2file)

      if (.not.rtdb_get(rtdb, 'ccsd:restart',MT_LOG, 1, orestart))
     $    orestart = .true.   ! By default always try to restart

      if (orestart) orestart = 
     $     ccsd_read_restart(t2file, g_t2, nvir, lhoov)

c Note that if restart fails orestart will be set false
      if (.not. orestart) then
         if (iam .eq. 0) write(6,*) ' Using MP2 initial guess vector '
        call dfill(lnov,zip,t1,1)
        if (nvp1.ge.ilo.and.nvp1.le.ihi)then
         call ga_put(g_t2,1,lnov,nvp1,nvp1,t1,1)
        endif
        call ccsd_mkt2(ncor,nocc,nvir,nact,eorb,scra,
     &                 scrb,g_nexch,g_t2,iprt)
      endif
      call ccsd_mktask(basis,nsh,tklst)
c
c ------------------------------------------------------------
c create square GAs with proposed final ordering
c ------------------------------------------------------------
      if (.not.ga_create(MT_DBL,lnov,lnov,'nt2',
     &                   nvir,nvir,g_nt2))
     &     call errquit('ga_create g_nt2 failed',0)
      if (.not.ga_create(MT_DBL,lnov,lnov,'nht2',
     &                   nvir,nvir,g_nht2))
     &     call errquit('ga_create g_nht2 failed',0)
c
      call ga_sync()
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
      if (iam.eq.0)write(6,1234)
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
      if (.not. ccsd_write_restart(t2file, g_t2, 
     $     nvir, lhoov, offt2))
     $     call errquit('ccsd: failed saving vector for restart',0)
c
      oconverged = .false.
      if (dabs(rms).lt.thre.and.dabs(deccsd).lt.thre)then
         oconverged = .true.
        if (iam.eq.0)write(6,1236)
 1236   format(18x,'*************converged*************',18x,
     &         /,73('-'),/)
        goto 999 
      endif
      if (iter.gt.maxit)then
        if (iam.eq.0)write(6,1237)
 1237   format(18x,'****maximum iterations exceeded****',18x,
     &         /,73('-'),/)
        goto 999
      endif
      if (iprt.gt.5.and.iam.eq.0)write(6,134)iter,idiis,ndiis
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
       write(6,139)(t1(i),i=1,lnov)
  139   format(4e14.4)
       do av=1,nvir
        call ga_get(g_t2,offt2+1,offt2+lhoov,av,av,scra,lhoov)
        print *,'starting t2 for av ',av
        write(6,139)(scra(i),i=1,lhoov)
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
      call ga_distribution(g_ht2,iam,jlo,jhi,ilo,ihi)
      call ga_sync
      do bev=1,nvir+1
         if (bev.ge.ilo.and.bev.le.ihi)then
            call ga_put(g_ht2,offt2+1,offt2+lhoov,bev,bev,scra,1)
         endif
      enddo
      if (nvp1.ge.ilo.and.nvp1.le.ihi)then
         call ga_put(g_ht2,offt2+1,offt2+lnov,nvp1,nvp1,scra,1)
      endif
      call ga_sync()
      tx(1)=tcgtime()
      call ccsd_pampt2(basis,idiis,nsh,ncor,nocc,nvir,nact,nbf,
     &                iprt,cmo,eorb,t1,ht1,
     $                scra,scrb,hiu,hia,habe,giu,gabe,
     &                g_t2,g_ht2,g_nt2,g_nht2,eccsd,
     &                g_ncoul,g_nexch,tklst, Tol2E, occd)
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
      if (iam.eq.0)then
*rak      write(6,1235)iter,eccsd,deccsd,rms,tx(2)-tx(1),tx(4)-tx(3)-
*rak     &                                       (tx(2)-tx(1))
*rak 1235 format(i4,d24.14,2d12.4,2d12.4)
          write(6,1235)iter,eccsd,deccsd,rms,
     &      tx(2)-tx(1),
     &      tx(4)-tx(3)-(tx(2)-tx(1)),
     &      main_block_time
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
         goto 999
      endif
c
 999  continue
c
      if (.not.ga_destroy(g_nht2))
     &    call errquit('ga_dest g_nht2 fail',0)
      if (.not.ga_destroy(g_nt2))
     &    call errquit('ga_dest g_nt2 fail',0)
c
      call ga_sync()
      call qexit('iterdrv',0)
      return
      end
      logical function ccsd_read_restart(t2file, g_t2, nvir, lhoov)
      implicit none
#include "global.fh"
#include "inp.fh"
      character*(*) t2file
      integer nvir, lhoov
      integer g_t2
c
      logical file_read_ga_patch
      external file_read_ga_patch
c
      call ga_sync
c
      ccsd_read_restart = 
     $     file_read_ga_patch(t2file,g_t2,1,lhoov,1,nvir+1)

      if (ga_nodeid().eq. 0) then
         if (ccsd_read_restart) then
            write(6,*) ' Restarting from old vector in ',
     $           t2file(1:inp_strlen(t2file))
         else
            write(6,*) ' Failed reading restart vector from ',
     $           t2file(1:inp_strlen(t2file))
         endif
      endif
c
      end
      logical function ccsd_write_restart(t2file, g_t2, nvir, lhoov, 
     $     offt2)
      implicit none
#include "global.fh"
      character*(*) t2file
      integer nvir, lhoov, offt2
      integer g_t2
c
      logical file_write_ga_patch
      external file_write_ga_patch
c
      call ga_sync
c
      ccsd_write_restart = 
     $     file_write_ga_patch(t2file,g_t2,offt2+1,offt2+lhoov,1,nvir+1)

      if ((.not. ccsd_write_restart) .and. ga_nodeid().eq. 0) 
     $     write(6,*) ' Write of T2 failed'
c
      end

         
