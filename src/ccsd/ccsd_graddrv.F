      subroutine ccsd_graddrv(basis,nsh,ncor,nocc,nvir,nact,nbf,
     &                        maxit,convi,iprt,cmo,eorb,t1,
     &                        hz1,hz2,scra,scrb,z1,hiu,hia,habe,
     &                        giu,gabe,bbkp,bb,g_t2,g_z2,g_hz2,
     &                        g_moint,mxvec,eccsd,max2e,mem2,
     &                        g_ncoul,g_nexch,tklst)
      implicit none
c
      integer msg_cc_z1a
c
      integer basis,nsh,ncor,nocc,nvir,nact,nbf,maxit,convi,iprt,
     &        g_t2,g_z2,g_hz2,g_moint,mxvec,g_ncoul,g_nexch,
     &        tklst(nsh*(nsh+1)/2,2)
c
       double precision cmo(nbf,nbf),eorb(nbf),t1(*),
     &                  hz1(*),hz2(*),scra(*),scrb(*),z1(*),
     &                  hiu(*),hia(*),habe(*),giu(*),gabe(*),
     &                  bbkp(*),bb(*),eccsd
c
c
c$$$      double precision rtc
c$$$      external rtc
c
c
      integer iter,idiis,ndiis,bev,jhi,jlo,ihi,ilo,offt2,nvp1,
     &        i,av
      integer max2e,mem2,max1e,mem1
      double precision eold,thre,dabs,tx(4),rx(2),rms,deccsd,
     &                 tstart,tend
      double precision zip,one,half,two,ten
      data zip/0.0d00/,one/1.0d00/,ten/10.0d00/
      integer ii,jj
      Integer Nodes, IAm
C
#include "msgids.fh"
#include "tcgmsg.fh"
#include "global.fh"
#include "ccsd_len.fh"
c
c
      Nodes = GA_NNodes()
      IAM = GA_NodeID()
C
      call ga_sync()
      call qenter('graddrv',0)
      call ga_distribution(g_hz2,iam,jlo,jhi,ilo,ihi)
c
c      print *,' arrays ',g_z2,g_hz2,g_moint,g_aoint
c      print *,'basis set handle ',basis
      nvp1=nvir+1
c - starting vectors
      call dfill(lnov,zip,z1,1)
      if (nvp1.ge.ilo.and.nvp1.le.ihi)then
       call ga_put(g_z2,1,lnov,nvp1,nvp1,z1,1)
      endif
      do av=1,nvir
        call dfill(lnoov,0.0d0,scra,1)
        call ga_put(g_z2,1,lnoov,av,av,scra,lnoov)
      enddo
c
      call ccsd_mktask(basis,nsh,tklst)

c
      call ga_sync()
      iter=0
      idiis=0
      ndiis=0
      eccsd=zip
      eold=zip
      thre=ten**convi
      rms=thre+thre
      deccsd=thre+thre
c
c ----------------------------------------------------------
c begin iterations
c ----------------------------------------------------------
      if (iam.eq.0)write(6,1234)
 1234 format(//,71(1h-),/,' iter',9x,'rms',/,71(1h-))
   99 continue
      iter=iter+1
      idiis=idiis+1
      if (idiis.gt.mxvec)idiis=1
      ndiis=ndiis+1
      if (ndiis.gt.mxvec)ndiis=mxvec
      offt2=(idiis-1)*lnoov
      if (dabs(rms).lt.thre)then
        if (iam.eq.0)write(6,1236)
 1236   format(18x,'*************converged*************',18x,
     &         /,71(1h-),/)
        goto 999 
      endif
      if (iter.gt.maxit)then
        if (iam.eq.0)write(6,1237)
 1237   format(18x,'****maximum iterations exceeded****',18x,
     &         /,71(1h-),/)
        goto 999
      endif
      if (iprt.gt.5)write(6,134)iter,idiis,ndiis
 134  format(' iteration ',i3,' idiis ',i3,' ndiis ',i3)
      if (iam.eq.0)then
       call ga_get(g_z2,offt2+1,offt2+lnov,nvp1,nvp1,z1,lnov)
      endif
c
      if (iam.eq.0.and.iprt.gt.5)then
       print *,'starting z1 vector'
       write(6,139)(z1(i),i=1,lnov)
       do av=1,nvir
        call ga_get(g_z2,offt2+1,offt2+lnoov,av,av,scra,lnoov)
        print *,'starting z2 for av ',av
        write(6,139)(scra(i),i=1,lnoov)
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
      call dfill(lnoov,zip,hz2,1)
      do bev=1,nvir+1
      if (bev.ge.ilo.and.bev.le.ihi)then
        call ga_put(g_hz2,offt2+1,offt2+lnoov,bev,bev,hz2,1)
      endif
      enddo
      if (nvp1.ge.ilo.and.nvp1.le.ihi)then
       call ga_put(g_hz2,offt2+1,offt2+lnov,nvp1,nvp1,hz2,1)
      endif
      call ga_sync()
      tx(1)=tcgtime()
c$$$      rx(1)=rtc()
      call ccsd_zvecs(basis,idiis,nsh,ncor,nocc,nvir,nact,nbf,
     &                iprt,cmo,eorb,t1,hz1,hz2,
     $                scra,scrb,z1,hiu,hia,habe,giu,gabe,
     &                g_t2,g_z2,g_hz2,g_moint,eccsd,
     &                max2e,mem2,g_ncoul,g_nexch,tklst)
      call ga_sync()
      tx(2)=tcgtime()
c$$$      rx(2)=rtc()
      tstart=tx(2)
      call ccsd_pdiis(nocc,nvir,idiis,ndiis,z1,hz1,scra,scrb,mxvec,
     &                 bbkp,bb,iprt,iter,rms,g_z2,g_hz2)
      call ga_sync()
      deccsd=eccsd-eold
      eold=eccsd
      tx(4)=tcgtime()
      tend=tx(4)
      if (iprt.gt.5.and.iam.eq.0)print *,' time pdiis2 ',tend-tstart
      if (iam.eq.0)then
      write(6,1235)iter,rms,tx(2)-tx(1),tx(4)-tx(3)-
     &                                       (tx(2)-tx(1))
 1235 format(i4,d12.4,2d12.4)
c$$$      write(6,*)'Wall clock time',(rx(2)-rx(1))*3.33e-09
      endif
      goto 99
 999  continue
c
      call ga_sync()
      call qexit('graddrv',0)
      return
      end

