C> \ingroup selci
C> @{
C>
      subroutine selci_eijil(e,i,j,l,ns,indbar,iocc,
     $     w1,w2,work,numf,numf2)
*
* $Id$
*
#include "implicit.fh"
#include "errquit.fh"
#include "ciinfo.fh"
      dimension e(*),indbar(*),iocc(*),work(*),
     $     w1(nfmax,nfmax,nsmax),w2(nfmax2,nfmax,nsmax-1)
c
c     e(u,v) = <Iu|Eij,il|Iv>, where u,v label the spin functions
c     and the orbital occupation of I is specified by ns, iocc
c     and indbar. Only for i.ne.j.ne.l. j and l are assumed
c     to be in order of increasing occupation in I
c
c     Note that this value is independent of i and can be simplified later
c
c     possible cases are:
c
c     iocc(j)  iocc(l)  case
c       0         0       1
c       0         1       2
c       1         1       3
c
      if (i.eq.j .or. i.eq.l .or. j.eq.l)
     $     call errquit('eijil: one of i=j,i=l,j=l',i, UNKNOWN_ERR)
      if (iocc(j).gt.iocc(l)) 
     $     call errquit('eijil: iocc ordering',iocc(j), UNKNOWN_ERR)
      numf = nf(ns)
      numf2 = nf(ns+2)
      ib = indbar(i) + (nsmax-ns-2)
*     call selci_mxma(w2(1,1,ib),1,nfmax2,w1(1,1,ib),1,nfmax,
*     $     e,1,numf,numf,numf2,numf2)
      call selci_axb(w2(1,1,ib),nfmax2,w1(1,1,ib),nfmax,
     $     e,numf,numf,numf2,numf2)
c
      icase = iocc(j) + iocc(l) + 1
      goto (10,20,30) icase
      call errquit('eijil: invalid case',icase, UNKNOWN_ERR)
c
c     nj=0, nl=0: e(u,v) = eiaeib(u,p)*w1(q,p,lb)*w1(v,q,jb)
c                                    ab     ->  la   ->     jl
c
 10   continue
      lb = indbar(l) + (nsmax-ns-2)
      jb = indbar(j) + (nsmax-ns-2)
      if (j.gt.l) jb = jb + 1
*      call selci_mxma(e,1,numf,w1(1,1,lb),nfmax,1,
*     $     work,1,numf,numf,numf2,numf2)
      call selci_axbt(e,numf,w1(1,1,lb),nfmax,
     $     work,numf,numf,numf2,numf2)
*      call selci_mxma(work,1,numf,w1(1,1,jb),nfmax,1,
*     $     e,1,numf,numf,numf2,numf2)
      call selci_axbt(work,numf,w1(1,1,jb),nfmax,
     $     e,numf,numf,numf2,numf2)
      return
c
c     nj=0, nl=1: e(u,v) = eiaeib(u,p)*w2(q,p,lb)*w1(v,q,jb)
c                                   lab     -> lla   ->    llj
c
 20   continue
      lb = indbar(l) + (nsmax-ns-2)
      jb = indbar(j) + (nsmax-ns)
      if (j.gt.l) jb = jb - 1
*      call selci_mxma(e,1,numf,w2(1,1,lb),nfmax2,1,
*     $     work,1,numf,numf,numf2,numf)
      call selci_axbt(e,numf,w2(1,1,lb),nfmax2,
     $     work,numf,numf,numf2,numf)
*      call selci_mxma(work,1,numf,w1(1,1,jb),nfmax,1,
*     $     e,1,numf,numf,numf,numf)
      call selci_axbt(work,numf,w1(1,1,jb),nfmax,
     $     e,numf,numf,numf,numf)
      numf2 = numf
      return
c
c     nj=1, nl=1: e(u,v) = eiaeib(u,p)*w2(q,p,lb)*w1(v,q,jb)
c                                  jlab    -> llja   ->   jjll
c
 30   continue
      lb = indbar(l) + (nsmax-ns-2)
      jb = indbar(j) + (nsmax-ns)
      if (j.gt.l) jb = jb - 1
      numf2m = nf(ns-2)
*      call selci_mxma(e,1,numf,w2(1,1,lb),nfmax2,1,
*     $     work,1,numf,numf,numf2,numf)
      call selci_axbt(e,numf,w2(1,1,lb),nfmax2,
     $     work,numf,numf,numf2,numf)
*      call selci_mxma(work,1,numf,w2(1,1,jb),nfmax2,1,
*     $     e,1,numf,numf,numf,numf2m)
      call selci_axbt(work,numf,w2(1,1,jb),nfmax2,
     $     e,numf,numf,numf,numf2m)
      numf2 = numf2m
c
      end
C>
C> @}
