C> \ingroup selci
C> @{
C>
      subroutine selci_eijkj(e,i,k,j,ns,indbar,iocc,
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
c     e(u,v) = <Iu|Eij,kj|Iv>, where u,v label the spin functions
c     and the orbital occupation of I is specified by ns, iocc
c     and indbar. Only for i.ne.j.ne.k. i and k are assumed
c     to be in order of decreasing occupaion in I
c
c     Note that this value is independent of j and can be simplified later
c
c     possible cases are:
c
c     iocc(i)  iocc(k)  case
c       1         1       1
c       3         1       2
c       3         3       3
c
      numf = nf(ns)
      icase = (iocc(i)+iocc(k))/2
      goto (10,20,30) icase
      call errquit('eijkj: invalid case',icase, INPUT_ERR)
c
c     ni=1, nk=1: e(u,v) = w1(u,q,ib)*w1(q,p,kb)*w1(s,p,jb)*w2(v,s,jb)
c
 10   continue
      numf2 = nf(ns-2)
      ib = indbar(i) + (nsmax-ns)
      kb = indbar(k) + (nsmax-ns)
      if (k.gt.i) kb = kb - 1
      jb = indbar(j) + (nsmax-ns)
      if (j.gt.i) jb = jb - 1
      if (j.gt.k) jb = jb - 1
*      call selci_mxma(w1(1,1,ib),1,nfmax,w1(1,1,kb),1,nfmax,
*     $     e,1,numf,numf,numf,numf)
      call selci_axb(w1(1,1,ib),nfmax,w1(1,1,kb),nfmax,
     $     e,numf,numf,numf,numf)
*      call selci_mxma(e,1,numf,w1(1,1,jb),nfmax,1,
*     $     work,1,numf,numf,numf,numf)
      call selci_axbt(e,numf,w1(1,1,jb),nfmax,
     $     work,numf,numf,numf,numf)
*      call selci_mxma(work,1,numf,w2(1,1,jb),nfmax2,1,
*     $     e,1,numf,numf,numf,numf2)
      call selci_axbt(work,numf,w2(1,1,jb),nfmax2,
     $     e,numf,numf,numf,numf2)
      return
c
c     ni=2, nk=1: e(u,v) = w2(u,q,ib)*w1(q,p,kb)*w1(s,p,jb)*w2(v,s,jb)
c
 20   continue
      numf2 = nf(ns+2)
      ib = indbar(i) + (nsmax-ns-2)
      kb = indbar(k) + (nsmax-ns-2)
      if (k.gt.i) kb = kb + 1
      jb = indbar(j) + (nsmax-ns-2)
      if (j.gt.i) jb = jb + 1
      if (j.gt.k) jb = jb - 1
*      call selci_mxma(w2(1,1,ib),1,nfmax2,w1(1,1,kb),1,nfmax,
*     $     e,1,numf,numf,numf2,numf2)
      call selci_axb(w2(1,1,ib),nfmax2,w1(1,1,kb),nfmax,
     $     e,numf,numf,numf2,numf2)
*      call selci_mxma(e,1,numf,w1(1,1,jb),nfmax,1,
*     $     work,1,numf,numf,numf2,numf2)
      call selci_axbt(e,numf,w1(1,1,jb),nfmax,
     $     work,numf,numf,numf2,numf2)
*      call selci_mxma(work,1,numf,w2(1,1,jb),nfmax2,1,
*     $     e,1,numf,numf,numf2,numf)
      call selci_axbt(work,numf,w2(1,1,jb),nfmax2,
     $     e,numf,numf,numf2,numf)
      numf2 = numf
      return
c
c     ni=2, nk=2: e(u,v) = w2(u,q,ib)*w2(q,p,kb)*w1(s,p,jb)*w2(v,s,jb)
c                       iikk   ->  kkia   ->  ikab    ->  ijka   ->  jjik
c
 30   continue
      numf2 = nf(ns+2)
      numf4 = nf(ns+4)
      ib = indbar(i) + (nsmax-ns-2)
      kb = indbar(k) + (nsmax-ns-4)
      if (k.gt.i) kb = kb + 1
      jb = indbar(j) + (nsmax-ns-4)
      if (j.gt.i) jb = jb + 1
      if (j.gt.k) jb = jb + 1
*      call selci_mxma(w2(1,1,ib),1,nfmax2,w2(1,1,kb),1,nfmax2,
*     $     e,1,numf,numf,numf2,numf4)
      call selci_axb(w2(1,1,ib),nfmax2,w2(1,1,kb),nfmax2,
     $     e,numf,numf,numf2,numf4)
*      call selci_mxma(e,1,numf,w1(1,1,jb),nfmax,1,
*     $     work,1,numf,numf,numf4,numf4)
      call selci_axbt(e,numf,w1(1,1,jb),nfmax,
     $     work,numf,numf,numf4,numf4)
*      call selci_mxma(work,1,numf,w2(1,1,jb),nfmax2,1,
*     $     e,1,numf,numf,numf4,numf2)
      call selci_axbt(work,numf,w2(1,1,jb),nfmax2,
     $     e,numf,numf,numf4,numf2)
      return
c
      end
C>
C> @}
