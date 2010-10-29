      subroutine ccsd_pdiis(nocc,nvir,idiis,ndiis,t1,ht1,scra,
     $                      scrb,mxvec,bbkp,bb,iprt,iter,rms,
     $                      g_newt2,g_newht2,lhoov)
c
c -- on input t1a contains t1: t1b contains error vector
      implicit none
C     $Id$
c
#include "ccsd_len.fh"
      integer lhoo,lhoov
      double precision t1(lnov),ht1(lnov),scra(lhoov),scrb(lhoov),
     &                 bbkp(*),bb(*),rms
      integer nocc,nvir,idiis,ndiis,mxvec,iprt,iter,g_newt2,g_newht2
c
      double precision sum(2),xmax,ddot,det
      integer i,j,av,bv,jlo,jhi,ilo,ihi,idp1,ndp1,nvp1,oft2,
     &        oft2p1,it,offbb,ad1
C
#include "msgids.fh"
#include "global.fh"
#include "ccsdps.fh"
c
      INTEGER INFO
      INTEGER ITMP(100)
      Integer IAm
c
      double precision zip,one
      data zip/0.0d00/,one/1.0d00/
c
      det = 0.0d00
      IAm = GA_NodeID()
      if (occsdps) then
         call pstat_on(ps_pdiis)
      else
         call qenter('pdiis',0)
      endif
c
      lhoo=nocc*(nocc+1)/2
      lhoov=nvir*nocc*(nocc+1)/2
      call ga_distribution(g_newt2,iam,jlo,jhi,ilo,ihi)
c
      idp1=idiis+1
      if (idp1.gt.mxvec)idp1=1
      if (iprt.gt.5.and.iam.eq.0)then
        print *,'iter,idiis,idp1,ndiis ',iter,idiis,idp1,ndiis
      endif
      ndp1=ndiis+1
      nvp1=nvir+1
      sum(1)=zip
      sum(2)=zip
      oft2=(idiis-1)*lhoov
      oft2p1=(idp1-1)*lhoov
c -- t1 error already in t1
      if (iam.eq.0)then
       do i=1,lnov
        sum(1)=sum(1)+t1(i)
        sum(2)=sum(2)+t1(i)*t1(i)
       enddo
      endif
      if (iam.eq.0.and.iprt.gt.5)print *,'sum1,sum2',sum(1),sum(2)
c
c -- t2 error in global array
      do av=1,nvir
       if (av.ge.ilo.and.av.le.ihi)then
        call ga_get(g_newht2,oft2+1,oft2+lhoov,av,av,scra,lhoov)
c        print *,' t2 error ',scra
        do bv=1,av
         do i=(bv-1)*lhoo+1,bv*lhoo
          sum(1)=sum(1)+scra(i)
          sum(2)=sum(2)+scra(i)*scra(i)
         enddo
        enddo
       endif
      enddo
      call ga_dgop(msg_cc_diis1,sum(1),2, '+')
      rms=dsqrt(dabs(sum(2)-(sum(1)*sum(1))))
      if (iam.eq.0.and.iprt.gt.5)
     &        print *,'sum1,sum2,rms ',sum(1),sum(2),rms

c -- compute new bb matrix elements
      if (iam.eq.0)then      
       do it=1,ndiis
        ad1=(it-1)*lhoov
        call ga_get(g_newht2,ad1+1,ad1+lnov,nvp1,nvp1,ht1,lnov)
        ad1=max(idiis,it)*(max(idiis,it)-1)/2+min(idiis,it)
        bbkp(ad1)=ddot(lnov,ht1,1,t1,1)
       enddo
      else
       call dfill(ndiis*(ndiis+1)/2,zip,bbkp,1)
      endif
c
      do av=1,nvir
       if (av.ge.ilo.and.av.le.ihi)then
        call ga_get(g_newht2,oft2+1,oft2+lhoov,av,av,scra,lhoov)
        do it=1,ndiis
         offbb=max(idiis,it)*(max(idiis,it)-1)/2+min(idiis,it)
         ad1=(it-1)*lhoov
         call ga_get(g_newht2,ad1+1,ad1+lhoov,av,av,scrb,lhoov)
         bbkp(offbb)=bbkp(offbb)+ddot(av*lhoo,scra,1,scrb,1)
        enddo
       endif
      enddo
      call ga_dgop(msg_cc_diis2,bbkp,ndiis*(ndiis+1)/2, '+')
c
      xmax=zip
      do i=1,ndiis*(ndiis+1)/2
       xmax=max(abs(bbkp(i)),xmax)
      enddo
      if (iprt.gt.5.and.iam.eq.0)print *,'xmax ',xmax
      xmax=one/xmax
      ad1=0
      offbb=ndp1*ndp1
      do i=1,ndiis
       do j=1,i
        ad1=ad1+1
        bb((i-1)*ndp1+j)=bbkp(ad1)*xmax
        bb((j-1)*ndp1+i)=bbkp(ad1)*xmax
       enddo
       bb((i-1)*ndp1+ndp1)=-one
       bb((ndp1-1)*ndp1+i)=-one
       bb(offbb+i)=zip
      enddo
      bb(ndp1*ndp1)=zip
      bb(offbb+ndp1)=-one
c
      if (iam.eq.0.and.iprt.gt.5)then
        print *,'error matrix for inversion'
        do i=1,ndp1
         write(*,201)(bb((i-1)*ndp1+j),j=1,ndp1)
        enddo
        print *,'right hand side'
        write(6,201)(bb(offbb+j),j=1,ndp1)
 201    format(5e15.5)
      endif
c
C      call flinq(bb,ndp1,ndp1,1,det)
c
      CALL dgesv( NDP1, 1, BB, NDP1, ITMP, BB(OFFBB+1), OFFBB+1, INFO )
c
      if (iam.eq.0.and.iprt.gt.5)then
        print *,'det = ',det
        print *,'solution vector'
        write(*,201)(bb(ndp1*ndp1+j),j=1,ndiis)
      endif
c
c -- compute new t1 vector
      call dfill(lnov,zip,t1,1)
      if (iam.eq.0)then
       offbb=ndp1*ndp1
       do it=1,ndiis
        ad1=(it-1)*lhoov
        call ga_get(g_newt2,ad1+1,ad1+lnov,nvp1,nvp1,ht1,lnov)
        call daxpy(lnov,bb(offbb+it),ht1,1,t1,1)
       enddo
       call ga_put(g_newt2,oft2p1+1,oft2p1+lnov,nvp1,nvp1,t1,lnov)
      endif
c
      offbb=ndp1*ndp1
      do av=1,nvir
       if (av.ge.ilo.and.av.le.ihi)then
        call dfill(lhoov,zip,scra,1)
        do it=1,ndiis
         ad1=(it-1)*lhoov
         call ga_get(g_newt2,ad1+1,ad1+lhoov,av,av,scrb,lhoov)
         call daxpy(lhoov,bb(offbb+it),scrb,1,scra,1)
        enddo
        call ga_put(g_newt2,oft2p1+1,oft2p1+lhoov,av,av,scra,lhoov)
       endif
      enddo
c
      if (occsdps) then
         call pstat_off(ps_pdiis)
      else
         call qexit('pdiis',0)
      endif
c
      return
      end
