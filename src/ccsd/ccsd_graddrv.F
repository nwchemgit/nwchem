      subroutine ccsd_graddrv(basis,nsh,ncor,nocc,nvir,nact,nbf,
     &                        maxit,Thre,iprt,cmo,eorb,t1,
     &                        hz1,scra,scrb,z1,hiu,hia,habe,
     &                        giu,gabe,bbkp,bb,g_nt2,g_z2,g_hz2,
     &                        mxvec,eccsd,max2e,mem2,
     &                        g_ncoul,g_nexch,tklst, Tol2E)
C     $Id$
      implicit none
#include "errquit.fh"
c
      integer msg_cc_z1a
c
      integer basis,nsh,ncor,nocc,nvir,nact,nbf,maxit,iprt,
     &        g_nt2,g_z2,g_hz2,g_nz2,g_nhz2,mxvec,g_ncoul,g_nexch,
     &        tklst(nsh*(nsh+1)/2,2)
      Double Precision Thre, Tol2E
c
       double precision cmo(nbf,nbf),eorb(nbf),t1(*),
     &                  hz1(*),scra(*),scrb(*),z1(*),
     &                  hiu(*),hia(*),habe(*),giu(*),gabe(*),
     &                  bbkp(*),bb(*),eccsd
c
c
c$$$      double precision rtc
c$$$      external rtc
c
c
      integer iter,idiis,ndiis,bev,jhi,jlo,ihi,ilo,offt2,nvp1,
     &        i,av,lhoo,lhoov
      integer max2e,mem2
      double precision eold,dabs,tx(4),rms,
     &                 tstart,tend
      double precision zip
      Integer IAm
      Logical oprint
C
#include "msgids.fh"
#include "tcgmsg.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "ccsd_len.fh"
#include "ccsdps.fh"
#include "util.fh"
c
* DATA MUST FOLLOW ALL DELCARATIONS 
      data zip/0.0d00/
c
      IAM = GA_NodeID()
C
      call ga_sync()
      if (occsdps) then
         call pstat_on(ps_graddrv)
      else
         call qenter('graddrv',0)
      endif
c
      oprint = util_print('information', print_low)
c
      lhoo=nocc*(nocc+1)/2
      lhoov=nvir*nocc*(nocc+1)/2
c
      nvp1=nvir+1
c - starting vectors
      call dfill(lnov,zip,z1,1)
c
      call ga_distribution(g_z2,iam,jlo,jhi,ilo,ihi)
      if (nvp1.ge.ilo.and.nvp1.le.ihi)then
       call ga_put(g_z2,1,lnov,nvp1,nvp1,z1,1)
      endif
      do av=1,nvir
        call dfill(lhoov,0.0d0,scra,1)
        call ga_put(g_z2,1,lhoov,av,av,scra,lhoov)
      enddo
c
      call ccsd_mktask(basis,nsh,tklst)
c
c ------------------------------------------------------------
c create square GAs with proposed final ordering
c ------------------------------------------------------------
      if (.not.ga_create(MT_DBL,lnov,lnov,'nz2',
     &                   nvir,nvir,g_nz2))
     &     call errquit('ga_create g_nz2 failed',0, GA_ERR)
      if (.not.ga_create(MT_DBL,lnov,lnov,'nhz2',
     &                   nvir,nvir,g_nhz2))
     &     call errquit('ga_create g_nhz2 failed',0, GA_ERR)
c
      call ga_sync()
      iter=0
      idiis=0
      ndiis=0
      eccsd=zip
      eold=zip
      rms=thre+thre
c
c ----------------------------------------------------------
c begin iterations
c ----------------------------------------------------------
      if (iam.eq.0.and.oprint)write(6,1234)
 1234 format(//,71(1h-),/,' iter',9x,'rms',/,71(1h-))
   99 continue
      iter=iter+1
      idiis=idiis+1
      if (idiis.gt.mxvec)idiis=1
      ndiis=ndiis+1
      if (ndiis.gt.mxvec)ndiis=mxvec
      offt2=(idiis-1)*lhoov
      if (dabs(rms).lt.thre)then
        if (iam.eq.0.and.oprint)write(6,1236)
 1236   format(18x,'*************converged*************',18x,
     &         /,71(1h-),/)
        goto 999 
      endif
      if (iter.gt.maxit)then
        if (iam.eq.0.and.oprint)write(6,1237)
 1237   format(18x,'****maximum iterations exceeded****',18x,
     &         /,71(1h-),/)
        goto 999
      endif
      if (iam.eq.0.and.iprt.gt.5.and.oprint)write(6,134)iter,idiis,ndiis
 134  format(' iteration ',i3,' idiis ',i3,' ndiis ',i3)
      if (iam.eq.0)then
       call ga_get(g_z2,offt2+1,offt2+lnov,nvp1,nvp1,z1,lnov)
      endif
c
      if (iam.eq.0.and.iprt.gt.5.and.oprint)then
       print *,'starting z1 vector'
       write(6,139)(z1(i),i=1,lnov)
       do av=1,nvir
        call ga_get(g_z2,offt2+1,offt2+lhoov,av,av,scra,lhoov)
        print *,'starting z2 for av ',av
        write(6,139)(scra(i),i=1,lhoov)
  139   format(4e14.4)
       enddo
      endif
c
      call ga_brdcst(msg_cc_z1a,z1,mdtob(lnov),0)
      tx(3)=tcgtime()
c
c ----------------------------------------------------------
c form ``sigma'' vector
c ----------------------------------------------------------
      call dfill(lnov,zip,hz1,1)
      call dfill(lhoov,zip,scra,1)
      call ga_distribution(g_hz2,iam,jlo,jhi,ilo,ihi)
      do bev=1,nvir+1
      if (bev.ge.ilo.and.bev.le.ihi)then
        call ga_put(g_hz2,offt2+1,offt2+lhoov,bev,bev,scra,1)
      endif
      enddo
      if (nvp1.ge.ilo.and.nvp1.le.ihi)then
       call ga_put(g_hz2,offt2+1,offt2+lnov,nvp1,nvp1,scra,1)
      endif
      call ga_sync()
      tx(1)=tcgtime()
      call ccsd_zvecs(basis,idiis,nsh,ncor,nocc,nvir,nact,nbf,
     &                iprt,cmo,eorb,t1,hz1,
     $                scra,scrb,z1,hiu,hia,habe,giu,gabe,
     &                g_nt2,g_z2,g_hz2,g_nz2,g_nhz2,
     &                g_ncoul,g_nexch,tklst, Tol2E)
      call ga_sync()
      tx(2)=tcgtime()
      tstart=tx(2)
      call ccsd_pdiis(nocc,nvir,idiis,ndiis,z1,hz1,scra,scrb,mxvec,
     &                 bbkp,bb,iprt,iter,rms,g_z2,g_hz2,lhoov)
      call ga_sync()
      eold=eccsd
      tx(4)=tcgtime()
      tend=tx(4)
      if (iprt.gt.5.and.iam.eq.0)print *,' time pdiis ',tend-tstart
      if (iam.eq.0.and.oprint)then
      write(6,1235)iter,rms,tx(2)-tx(1),tx(4)-tx(3)-
     &                                       (tx(2)-tx(1))
 1235 format(i4,d12.4,2d12.4)
      endif
      goto 99
 999  continue
c
      if (.not.ga_destroy(g_nhz2))
     &    call errquit('ga_dest g_nhz2 fail',0, GA_ERR)
      if (.not.ga_destroy(g_nz2))
     &    call errquit('ga_dest g_nz2 fail',0, GA_ERR)
c
      call ga_sync()
      if (occsdps) then
         call pstat_off(ps_graddrv)
      else
         call qexit('graddrv',0)
      endif
      return
      end

