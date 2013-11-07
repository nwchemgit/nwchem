C> \ingroup selci
C> @{
C>
      subroutine selci_makehs(hij,i,j,ns,indbar,iocc,lists,
     $     w1,w2,work1,work2,f,g,int12,int34,numf,numf2)
*
* $Id$
*
#include "implicit.fh"
#include "ciinfo.fh"
      dimension hij(*),indbar(*),iocc(*),lists(*),
     $     w1(nfmax,nfmax,nsmax),w2(nfmax2,nfmax,nsmax-1),
     $     work1(*),work2(*),f(*),g(*),int12(*),int34(*)
#include "stfunc.fh"
c
c     (i->j) |I> = |J>. Compute hij(u,v) = <Iu|H|Jv>.
c
c     hij(u,v) = { f(ij) + nj(jj|ij) }<Iu|Eij|Jv> +
c                sum(r socc .ne. i,j)[(ir|rj)<Iu|EirErj|Jv>]
c
c     do the <u|Eij|v> term
c
      call selci_eij(hij,i,j,ns,indbar,iocc,w1,w2)
      ii = iky(i) + i
      ij = itrian(i,j)
      jj = iky(j) + j
      fac = f(ij)
      nns = ns
      if (iocc(j).eq.1) then
         nns = nns - 2
         fac = fac + g(intadr(jj,ij))
      endif
      if (iocc(i).eq.3) nns = nns + 2
      numf = nf(ns)
      numf2 = nf(nns)
      nnf = numf*numf2
      call dscal(nnf,fac,hij,1)
c
c     do the two-electron term
c
      if (.not.(iocc(i).eq.1 .and. nf(ns-2).le.0)) then
         do 10 irrr = 1,ns
            ir = lists(irrr)
            if (ir.eq.i .or. ir.eq.j) goto 10
            call selci_eirerj(work2,i,ir,j,ns,indbar,iocc,w1,w2,work1)
            gg = g(intadr(itrian(ir,i),itrian(ir,j)))
            call daxpy(nnf,gg,work2,1,hij,1)
 10      continue
      endif
c
      end
C>
C> @}
