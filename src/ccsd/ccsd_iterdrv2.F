      subroutine ccsd_iterdrv2(basis,ncor,nocc,nvir,ndel,nbf,nact,nsh,
     &                        maxit,convi,iprt,cmo,eorb,t1,t2,
     &                        ht1,ht2,scra,scrb,hiu,hia,habe,giu,gabe,
     &                        bbkp,bb,g_t2,g_ht2,
     &                        g_moint,mxvec,eccsd,max2e,mem2,
     &                        g_ncoul,g_nexch,tklst)
      implicit none
c
      integer basis,ncor,nocc,nvir,ndel,nbf,nact,nsh,maxit,convi,iprt,
     &        g_t2,g_ht2,g_moint,mxvec,g_ncoul,g_nexch,
     &        tklst(nsh*(nsh+1)/2,2)
c
       double precision cmo(nbf,nbf),eorb(nbf),t1(*),t2(*),
     &                  ht1(*),ht2(*),scra(*),scrb(*),
     &                  hiu(*),hia(*),habe(*),giu(*),gabe(*),
     &                  bbkp(*),bb(*),eccsd
c
c$$$      double precision rtc
c$$$      external rtc
c
#include "ccsd_len.fh"
c
      integer iter,idiis,ndiis,bev,jhi,jlo,ihi,ilo,offt2,nvp1,
     &        i,av
      integer max2e,mem2,max1e,mem1
      double precision eold,thre,dabs,tx(4),rx(2),rms,deccsd,
     &                 tstart,tend

      integer ii,jj
      Integer Nodes, IAm

      double precision zip,one,half,two,ten
      data zip/0.0d00/,one/1.0d00/,ten/10.0d00/
C
#include "msgids.fh"
#include "tcgmsg.fh"
#include "global.fh"
c
c
      Nodes = GA_NNodes()
      IAm = GA_NodeID()
C
      call ga_sync()
      call qenter('iterdrv',0)
      call ga_distribution(g_ht2,iam,jlo,jhi,ilo,ihi)
c
c      print *,' arrays ',g_t2,g_ht2,g_moint,g_aoint
c      print *,'basis set handle ',basis
      nvp1=nvir+1
c - starting vectors
      call dfill(lnov,zip,t1,1)
capr
c      print *,'putting some rubbish in t1 to start with'
c      rms=one
c      do i=1,lnov
c       t1(i)=rms
c       rms=rms+one
c      enddo
capr
      if (nvp1.ge.ilo.and.nvp1.le.ihi)then
       call ga_put(g_t2,1,lnov,nvp1,nvp1,t1,1)
      endif
      call ccsd_mkt2(ncor,nocc,nvir,ndel,nbf,eorb,g_t2,g_moint,t2,iprt)
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
 1234 format(//,71(1h-),/,' iter',8x,'correlation',11x,'delta',9x,
     &       'rms',/,15x,'energy',14x,'energy',/,71(1h-))
   99 continue
      iter=iter+1
      idiis=idiis+1
      if (idiis.gt.mxvec)idiis=1
      ndiis=ndiis+1
      if (ndiis.gt.mxvec)ndiis=mxvec
      offt2=(idiis-1)*lnoov
      if (dabs(rms).lt.thre.and.dabs(deccsd).lt.thre)then
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
       call ga_get(g_t2,offt2+1,offt2+lnov,nvp1,nvp1,t1,lnov)
      endif
c
      if (iam.eq.0.and.iprt.gt.5)then
       print *,'starting t1 vector'
       write(6,139)(t1(i),i=1,lnov)
       do av=1,nvir
        call ga_get(g_t2,offt2+1,offt2+lnoov,av,av,scra,lnoov)
        print *,'starting t2 for av ',av
        write(6,139)(scra(i),i=1,lnoov)
  139   format(4e14.4)
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
      call dfill(lnoov,zip,ht2,1)
      do bev=1,nvir+1
      if (bev.ge.ilo.and.bev.le.ihi)then
        call ga_put(g_ht2,offt2+1,offt2+lnoov,bev,bev,ht2,1)
      endif
      enddo
      if (nvp1.ge.ilo.and.nvp1.le.ihi)then
       call ga_put(g_ht2,offt2+1,offt2+lnov,nvp1,nvp1,ht2,1)
      endif
      call ga_sync()
      tx(1)=tcgtime()
c$$$      rx(1)=rtc()
      call ccsd_pampt2(basis,idiis,ncor,nocc,nvir,ndel,nbf,nact,nsh,
     &                iprt,cmo,eorb,t1,t2,t2(lnoov+1),ht1,ht2,
     $                scra,scrb,hiu,hia,habe,giu,gabe,
     &                g_t2,g_ht2,g_moint,eccsd,
     &                max2e,mem2,g_ncoul,g_nexch,tklst)
      call ga_sync()
      tx(2)=tcgtime()
c$$$      rx(2)=rtc()
      tstart=tx(2)
      call ccsd_pdiis(nocc,nvir,idiis,ndiis,t1,ht1,scra,scrb,mxvec,
     &                 bbkp,bb,iprt,iter,rms,g_t2,g_ht2)
      call ga_sync()
      deccsd=eccsd-eold
      eold=eccsd
      tx(4)=tcgtime()
      tend=tx(4)
      if (iprt.gt.5.and.iam.eq.0)print *,' time pdiis2 ',tend-tstart
      if (iam.eq.0)then
      write(6,1235)iter,eccsd,deccsd,rms,tx(2)-tx(1),tx(4)-tx(3)-
     &                                       (tx(2)-tx(1))
 1235 format(i4,d24.14,2d12.4,2d12.4)
c$$$      write(6,*)'Wall clock time',(rx(2)-rx(1))*3.33e-09
      endif
      goto 99
 999  continue
c
c - put vectors in g_moint for use in triples, overwrite [vvoo] ints
      lnobj=lnoov+lnooo
      if (iam.eq.0)then
       call ga_get(g_t2,offt2+1,offt2+lnov,nvp1,nvp1,t1,lnov)
       call ga_put(g_moint,lnobj+1,lnobj+lnov,nvp1,nvp1,t1,lnov)
      endif
      do av=1,nvir
       if (av.ge.ilo.and.av.le.ihi)then
        call ga_get(g_t2,offt2+1,offt2+lnoov,av,av,scra,lnoov)
        call ga_put(g_moint,lnobj+1,lnobj+lnoov,av,av,scra,lnoov)
       endif
      enddo
      call ga_sync()
      call qexit('iterdrv',0)
      return
      end

