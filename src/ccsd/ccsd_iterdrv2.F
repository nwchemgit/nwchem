      subroutine ccsd_iterdrv2(rtdb,basis,nsh,ncor,nocc,nvir,nact,nbf,
     &                        maxit,Thre,iprt,cmo,eorb,t1,
     &                        ht1,scra,scrb,hiu,hia,habe,giu,gabe,
     &                        bbkp,bb,g_t2,g_ht2,
     &                        mxvec,eccsd,max2e,mem2,
     &                        g_ncoul,g_nexch,tklst,offt2,nvp1,
     $     refener, Tol2E)
C     $Id: ccsd_iterdrv2.F,v 2.14 1997-10-22 16:40:55 d3g681 Exp $
      implicit none
c
      integer rtdb,basis,nsh,ncor,nocc,nvir,nact,nbf,maxit,iprt,
     & g_t2,g_ht2,g_nt2,g_nht2,mxvec,g_ncoul,g_nexch,
     &        tklst(nsh*(nsh+1)/2,2)
      Double Precision Thre, Tol2E
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
      character*255 t1file, t2file

      logical ccsd_read_restart, ccsd_write_restart
      external ccsd_read_restart, ccsd_write_restart


C
#include "msgids.fh"
#include "tcgmsg.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
c
      double precision zip,ten
      data zip/0.0d00/,ten/10.0d00/
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

      call util_file_name('t1', .false.,.false.,t1file)
      call util_file_name('t2', .false.,.false.,t2file)

      if (.not.rtdb_get(rtdb, 'ccsd:restart',MT_LOG, 1, orestart))
     $    orestart = .false.

      if (orestart) orestart = 
     $     ccsd_read_restart(t1file, t2file, t1, g_t2, lnov)
c Note that if restart fails orestart will be set false
      if (.not. orestart) then
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
      iter=iter+1
      idiis=idiis+1
      if (idiis.gt.mxvec)idiis=1
      ndiis=ndiis+1
      if (ndiis.gt.mxvec)ndiis=mxvec
      offt2=(idiis-1)*lhoov

c
c     Write out the current T1 and T2 for possible restart
c
      if (.not. ccsd_write_restart(t1file, t2file, t1, g_t2, lnov))
     $     call errquit('ccsd: failed saving vector for restart',0)
c
      if (dabs(rms).lt.thre.and.dabs(deccsd).lt.thre)then
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
     &                g_ncoul,g_nexch,tklst, Tol2E)
      call ga_sync()
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
      goto 99
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

      logical function ccsd_read_restart(t1file, t2file, t1, g_t2, lnov)
      implicit none
#include "global.fh"
#include "eaf.fh"
      character*(*) t1file, t2file
      integer lnov
      double precision t1(lnov)
      integer g_t2
c
      integer fd, iam, ierr
      logical file_read_ga
      external file_read_ga
c
      call ga_sync
c
      iam = ga_nodeid()
c
      ccsd_read_restart = .false.
c
      if (.not. file_read_ga(t2file,g_t2)) then
         if (iam .eq. 0) write(6,*) ' Restart from T2 failed'
         return
      endif
c
      ierr = 0
c
      if (iam .eq. 0) then
         ierr = eaf_open(t1file, eaf_rw, fd)
         if (ierr .eq. 0) then
            ierr = eaf_read(fd, 0.0d0, t1, 8*lnov)
            if (ierr .ne. 0) write(6,*) ' Reading from T1 failed'
            if (eaf_close(fd) .ne. 0) call errquit('ccsd:closing t1',0)
         else
            write(6,*) ' Opening T1 failed'
         endif
      endif
c
      call ga_igop(321,ierr,1,'+')
      if (ierr .ne. 0) return
c
      call ga_brdcst(322,t1,8*lnov,0)
c
      call ga_sync
c
      ccsd_read_restart = .true.
c
      end

         
      logical function ccsd_write_restart(t1file, t2file, t1, g_t2,lnov)
      implicit none
#include "global.fh"
#include "eaf.fh"
      character*(*) t1file, t2file
      integer lnov
      double precision t1(lnov)
      integer g_t2
c
      integer fd, iam, ierr
      logical file_write_ga
      external file_write_ga
c
      call ga_sync
c
      iam = ga_nodeid()
c
      ccsd_write_restart = .false.
c
      if (.not. file_write_ga(t2file,g_t2)) then
         if (iam .eq. 0) write(6,*) ' Writing to T2 failed'
         return
      endif
c
      ierr = 0
c
      if (iam .eq. 0) then
         ierr = eaf_open(t1file, eaf_rw, fd)
         if (ierr .eq. 0) then
            ierr = eaf_write(fd, 0.0d0, t1, 8*lnov)
            if (ierr .ne. 0) write(6,*) ' Writing to T1 failed'
            if (eaf_close(fd) .ne. 0) call errquit('ccsd:closing t1',0)
         else
            write(6,*) ' Opening T1 failed'
         endif
      endif
c
      call ga_igop(321,ierr,1,'+')
      if (ierr .ne. 0) return
c
      call ga_sync
c
      ccsd_write_restart = .true.
c
      end

         
