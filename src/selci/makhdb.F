C> \ingroup selci
C> @{
C>
      subroutine selci_makhdb(hij,ind,icase,ns,indbar,iocc,
     $     w1,w2,work1,work2,g,int12,int34,numf,numf2)
*
* $Id$
*
#include "implicit.fh"
#include "errquit.fh"
#include "ciinfo.fh"
      dimension hij(*),ind(4),indbar(*),iocc(*),
     $     w1(nfmax,nfmax,nsmax),w2(nfmax2,nfmax,nsmax-1),
     $     work1(*),work2(*),g(*),int12(*),int34(*)
#include "stfunc.fh"
c
c     |I> is an orbital double replacement from |J>.
c     Compute the hamiltonian matrix element in hij.
c
c     icase=1
c        ii -> jj , ind(1)=i, ind(2)= j
c        hij(u,v) = (ij|ij) delta(u,v) (u,v label spin functions)
c     icase=2
c        ik -> jj , ind(1)=i, ind(2)=k, ind(3)=j
c        hij = (ij|kj)<u|Eij,kj|v> ... can be simplified
c     icase=3
c        ii -> jl , ind(1)=i, ind(2)=j, ind(3)=l
c        hij = (ij|il)<u|Eij,il|v> ... can be simplified
c     icase=4
c        ik -> jl , ind(1)=i, ind(2)=k, ind(3)=j, ind(4)=l
c        hij = (ij|kl)<u|Eij,kl|v> + (il|kj)<u|Eil,kj|v>
c
      numf = nf(ns)
      goto (10,20,30,40) icase
      call errquit('makhdb: invalid icase',icase, INPUT_ERR)
c
c     icase=1
c        ii -> jj , ind(1)=i, ind(2)= j
c        hij(u,v) = (ij|ij) delta(u,v) (u,v label spin functions)
c
 10   continue
      numf2 = nf(ns)
      call dfill(numf*numf,0.0d0,hij,1)
      ij = itrian(ind(1),ind(2))
      call dfill(numf,g(int12(ij)+int34(ij)),hij,numf+1)
      return
c
c     icase=2
c        ik -> jj , ind(1)=i, ind(2)=k, ind(3)=j
c        hij = (ij|kj)<u|Eij,kj|v>
c
 20   continue
      call selci_eijkj(hij,ind(1),ind(2),ind(3),ns,indbar,iocc,w1,w2,
     $     work1,numf,numf2)
      gg = g(intadr(itrian(ind(1),ind(3)),itrian(ind(2),ind(3))))
      call dscal(numf*numf2,gg,hij,1)
      return
c
c     icase=3
c        ii -> jl , ind(1)=i, ind(2)=j, ind(3)=l
c        hij = (ij|il)<u|Eij,il|v> ... can be simplified
c
 30   continue
      call selci_eijil(hij,ind(1),ind(2),ind(3),ns,indbar,iocc,w1,w2,
     $     work1,numf,numf2)
      gg = g(intadr(itrian(ind(1),ind(2)),itrian(ind(1),ind(3))))
      call dscal(numf*numf2,gg,hij,1)
      return
c
c     icase=4
c        ik -> jl , ind(1)=i, ind(2)=k, ind(3)=j, ind(4)=l
c        hij = (ij|kl)<u|Eij,kl|v> + (il|kj)<u|Eil,kj|v>
c
 40   continue
      call selci_eijkl(hij,ind(1),ind(3),ind(2),ind(4),ns,indbar,
     $     iocc,w1,w2,work1,numf,numf2)
      gg = g(intadr(itrian(ind(1),ind(3)),itrian(ind(2),ind(4))))
      call dscal(numf*numf2,gg,hij,1)
      call selci_eijkl(work2,ind(1),ind(4),ind(2),ind(3),ns,indbar,iocc,
     $     w1,w2,work1,numf,numf2)
      gg = g(intadr(itrian(ind(1),ind(4)),itrian(ind(3),ind(2))))
      call daxpy(numf*numf2,gg,work2,1,hij,1)
c
      end
C>
C> @}
