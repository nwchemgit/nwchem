C> \ingroup selci
C> @{
C>
      subroutine selci_eirerj(e,i,ir,j,ns,indbar,iocc,w1,w2,work)
*
* $Id$
*
#include "implicit.fh"
#include "errquit.fh"
#include "ciinfo.fh"
      dimension e(*),indbar(*),iocc(*),work(*),
     $     w1(nfmax,nfmax,nsmax),w2(nfmax2,nfmax,nsmax-1)
c
c     e(u,v) = <Iu|EirErj|Iv>, where u,v label the spin functions
c     and the orbital occupation of I is specified by ns, iocc
c     and indbar. Only for i.ne.j.ne.r, and r singly occupied.
c
c     e(u,v) = sum(q) <Iu|Eir|Sq><Sq|Erj|Iv>
c
c     <Iu|Eir|Sq> = sum(p) <Iu|Eia|Tp><Tp|Ear|Sq> 
c
c     possible cases are:
c
c     iocc(i)  iocc(j)  case
c       1         0       1
c       1         1       2
c       3         0       3
c       3         1       4
c
      if (iocc(ir).ne.1 .or. ir.eq.i .or. ir.eq.j .or. i.eq.j)
     $     call errquit('eirerj: wrong i,r,j',ir, INPUT_ERR)
c
      icase = iocc(i)+iocc(j)
      goto (10,20,30,40) icase
      call errquit('eirerj: fell thru goto',icase, INPUT_ERR)
c
 10   continue
c
c     ni=1, nj=0 : e(u,v) = w1(u,p,ib)*w2(q,p,rb)*w2(q,s,rb)*w1(v,s,jb)
c                         ir    ->    ra    ->   rr    ->   ra   ->   jr
c
      numf = nf(ns)
      numf2 = nf(ns-2)
      ib = indbar(i) + (nsmax-ns)
      irb = indbar(ir) + (nsmax-ns)
      if (ir.gt.i) irb = irb - 1
      jb = indbar(j) + (nsmax-ns)
      if (j.gt.i) jb = jb - 1
*      call selci_mxma(w1(1,1,ib),1,nfmax,w2(1,1,irb),nfmax2,1,
*     $     e,1,numf,numf,numf,numf2)
      call selci_axbt(w1(1,1,ib),nfmax,w2(1,1,irb),nfmax2,
     $     e,numf,numf,numf,numf2)
*      call selci_mxma(e,1,numf,w2(1,1,irb),1,nfmax2,
*     $     work,1,numf,numf,numf2,numf)
      call selci_axb(e,numf,w2(1,1,irb),nfmax2,
     $     work,numf,numf,numf2,numf)
*      call selci_mxma(work,1,numf,w1(1,1,jb),nfmax,1,
*     $     e,1,numf,numf,numf,numf)
      call selci_axbt(work,numf,w1(1,1,jb),nfmax,
     $     e,numf,numf,numf,numf)
      return
c
 20   continue
c
c     ni=1, nj=1 : e(u,v) = w1(u,p,ib)*w2(q,p,rb)*w2(q,s,rb)*w2(v,s,jb)
c                        ijr    ->   jra    ->  rrj   ->   jra  ->   jjr
c     same as 10 except for last mxma
      numf = nf(ns)
      numf2 = nf(ns-2)
      ib = indbar(i) + (nsmax-ns)
      irb = indbar(ir) + (nsmax-ns)
      if (ir.gt.i) irb = irb - 1
      jb = indbar(j) + (nsmax-ns)
      if (j.gt.i) jb = jb - 1
*      call selci_mxma(w1(1,1,ib),1,nfmax,w2(1,1,irb),nfmax2,1,
*     $     e,1,numf,numf,numf,numf2)
      call selci_axbt(w1(1,1,ib),nfmax,w2(1,1,irb),nfmax2,
     $     e,numf,numf,numf,numf2)
*      call selci_mxma(e,1,numf,w2(1,1,irb),1,nfmax2,
*     $     work,1,numf,numf,numf2,numf)
      call selci_axb(e,numf,w2(1,1,irb),nfmax2,
     $     work,numf,numf,numf2,numf)
*      call selci_mxma(work,1,numf,w2(1,1,jb),nfmax2,1,
*     $     e,1,numf,numf,numf,numf2)
      call selci_axbt(work,numf,w2(1,1,jb),nfmax2,
     $     e,numf,numf,numf,numf2)
      return
c
 30   continue
c
c     ni=2, nj=0 : e(u,v) = w2(u,p,ib)*w2(q,p,rb)*w2(q,s,rb)*w1(v,s,jb)
c                        iir    ->   ira    ->  rri   ->   ira  ->   ijr
c
      numf = nf(ns)
      numf2 = nf(ns+2)
      ib = indbar(i) + (nsmax-ns-2)
      irb = indbar(ir) + (nsmax-ns-2)
      if (ir.gt.i) irb = irb + 1
      jb = indbar(j) + (nsmax-ns-2)
      if (j.gt.i) jb = jb + 1
*      call selci_mxma(w2(1,1,ib),1,nfmax2,w2(1,1,irb),nfmax2,1,
*     $     e,1,numf,numf,numf2,numf)
      call selci_axbt(w2(1,1,ib),nfmax2,w2(1,1,irb),nfmax2,
     $     e,numf,numf,numf2,numf)
*      call selci_mxma(e,1,numf,w2(1,1,irb),1,nfmax2,
*     $     work,1,numf,numf,numf,numf2)
      call selci_axb(e,numf,w2(1,1,irb),nfmax2,
     $     work,numf,numf,numf,numf2)
*      call selci_mxma(work,1,numf,w1(1,1,jb),nfmax,1,
*     $     e,1,numf,numf,numf2,numf2)
      call selci_axbt(work,numf,w1(1,1,jb),nfmax,
     $     e,numf,numf,numf2,numf2)
      return
c
 40   continue
c
c     ni=2, nj=1 : e(u,v) = w2(u,p,ib)*w2(q,p,rb)*w2(q,s,rb)*w2(v,s,jb)
c                        iijr   ->   ijra   ->  rrij  ->   ijra  ->   jjir
c     same as 30 except for last mxma
      numf = nf(ns)
      numf2 = nf(ns+2)
      ib = indbar(i) + (nsmax-ns-2)
      irb = indbar(ir) + (nsmax-ns-2)
      if (ir.gt.i) irb = irb + 1
      jb = indbar(j) + (nsmax-ns-2)
      if (j.gt.i) jb = jb + 1
*      call selci_mxma(w2(1,1,ib),1,nfmax2,w2(1,1,irb),nfmax2,1,
*     $     e,1,numf,numf,numf2,numf)
      call selci_axbt(w2(1,1,ib),nfmax2,w2(1,1,irb),nfmax2,
     $     e,numf,numf,numf2,numf)
*      call selci_mxma(e,1,numf,w2(1,1,irb),1,nfmax2,
*     $     work,1,numf,numf,numf,numf2)
      call selci_axb(e,numf,w2(1,1,irb),nfmax2,
     $     work,numf,numf,numf,numf2)
*      call selci_mxma(work,1,numf,w2(1,1,jb),nfmax2,1,
*     $     e,1,numf,numf,numf2,numf)
      call selci_axbt(work,numf,w2(1,1,jb),nfmax2,
     $     e,numf,numf,numf2,numf)
      return
c
      end
C>
C> @}
