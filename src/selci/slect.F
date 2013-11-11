C> \ingroup selci
C> @{
      subroutine selci_slect(q,h,g,int12,int34,w1,w2,ioconf,indxci,
     $     roots,ci,nconmx,ncold,thresh,ept,enew,irange, ptnorm,
     &     ptnorm_mp,ept_mp,roots_mp,min2c)
*
* $Id$
*
#include "implicit.fh"
#include "ciinfo.fh"
#include "mptr.fh"
#include "global.fh"
#include "stdio.fh"
      dimension q(*)
      dimension h(nnorbs),g(numint),int12(nnorbs),int34(nnorbs),
     $     w1(nfmax,nfmax,nsmax),w2(nfmax2,nfmax,nsmax-1),
     $     ioconf(nintpo,nconmx),indxci(nconmx),roots(nroot),
     $     ci(nci,nroot),ept(nroot),enew(nroot),irange(21),
     $     ptnorm(nroot),
     &     ptnorm_mp(nroot), ept_mp(nroot), roots_mp(nroot)
c     
      dimension iocc(255), itemp(32)
#include "bitops.fh"
c     
c     loop thru single and double replacements from the current set
c     of orbital configurations. If this replacement interacts with
c     a previous configuration then we have already done this one.
c     Compute <*|H|I> where <*| label the spin functions in the new
c     configuration and I is the Ith root wavefunction. Then compute
c     the predicted energy lowering for each root. Add this to the
c     total PT corrected CI energy. If the energy lowering for any
c     root is greater than the input threshold then add this to the
c     list and add the incremental lowering to the prediction for the
c     new ci energy.
c     
      ivc = selci_mptr(nfmax*nroot)
      call dfill(nroot,0.0d0,ept,1)
      call dfill(nroot,0.0d0,ept_mp,1)
      call dfill(nroot,0.0d0,enew,1)
      call dfill(nroot,0.0d0,ptnorm,1)
      call dfill(nroot,0.0d0,ptnorm_mp,1)
      call ifill(21,0,irange,1)
c     
c     loop through orbital configurations already in ci
c     
      iconfm = 0
      nproc  = ga_nnodes()
      me     = ga_nodeid()
      master = 0
      icount = -1
c     
c     If first pass through selci then kill all restart files
c     
      if (ncold.eq.1.and.me.eq.0) call restk
      i10str = 1
      call selci_restin(i10str,ept,ept_mp,enew,ptnorm,ptnorm_mp,
     &     nroot,irange,21,iwpt,
     &     noconf,ioconf,nintpo,nconmx,ncold)
c     
c     Checkpoint after each 5% of configurations
c     
      iconfm = max(100,ncold/20)
c     
      do 10 iconf = i10str,ncold
c     
         if(min2c .gt. 0) then
            if(mod(iconf,iconfm).eq.0) then
               call selci_stool(iconf,ept,ept_mp,enew,ptnorm,ptnorm_mp,
     &              nroot,irange,21,iwpt,
     &              noconf,ioconf,nintpo,nconmx,ncold)
            endif
         endif
         if(ga_nodeid().eq.0) call util_flush(luout)
c
         call selci_upkcon(norbs, iocc, ioconf(1,iconf), nintpo, nbitpi)
c     
c     single replacements i -> j
c     
         do 20 i = 1,norbs
            if (iocc(i).eq.0) goto 20
c
            icount = icount + 1
            if (mod(icount,nproc).ne.me) goto 20
c
            iold = iocc(i)
            iocc(i) = iocc(i) - 1
            if (iocc(i).eq.2) iocc(i) = 1
            do 30 j = 1,norbs
               if (isym(i).ne.isym(j) .or.
     $              j.eq.i .or. iocc(j).eq.3) goto 30
               jold = iocc(j)
               iocc(j) = iocc(j) + 1
               if (iocc(j).eq.2) iocc(j) = 3
               call selci_pkcon(norbs, iocc, itemp, nintpo, nbitpi)
c     
               call selci_tester(q, h, g, int12, int34, w1, w2, ioconf,
     $              indxci, roots, ci, nconmx, thresh, ept, enew,
     $              iconf, ncold, q(ivc), iocc, itemp, irange, ptnorm,
     &              ptnorm_mp,ept_mp,roots_mp)
c     
               iocc(j) = jold
 30         continue
            iocc(i) = iold
 20      continue
c     
c     double replacements (i,j) -> (k,l) ... what makes you think
c     that this is going to be slooooow.
c     
         do 40 i = 1,norbs
            if (iocc(i).eq.0) goto 40
            iold = iocc(i)
            iocc(i) = iocc(i) - 1
            if (iocc(i).eq.2) iocc(i) = 1
            do 50 j = i,norbs
               if (iocc(j).eq.0) goto 50
               jold = iocc(j)
               iocc(j) = iocc(j) - 1
               if (iocc(j).eq.2) iocc(j) = 1
               do 60 k = 1,norbs
                  if (k.eq.i .or. k.eq.j .or. iocc(k).eq.3) goto 60
                  kold = iocc(k)
                  iocc(k) = iocc(k) + 1
                  if (iocc(k).eq.2) iocc(k) = 3
                  ijksym = ieor(ieor(isym(i),isym(j)),isym(k))
                  do 70 l = k,norbs
                     if (isym(l).ne.ijksym .or. l.eq.i .or. 
     $                    l.eq.j .or. iocc(l).eq.3) goto 70
c
                          icount = icount + 1
                          if (mod(icount,nproc).ne.me) goto 70
c
                          lold = iocc(l)
                          iocc(l) = iocc(l) + 1
                          if (iocc(l).eq.2) iocc(l) = 3
                          call selci_pkcon(norbs,iocc,itemp,
     $                         nintpo,nbitpi)
c     
                          call selci_tester(q, h, g, int12, int34,
     $                         w1, w2, 
     $                         ioconf, indxci, roots, ci, nconmx, 
     $                         thresh, ept, enew, iconf, ncold, 
     $                         q(ivc), iocc, itemp, irange, ptnorm,
     $                         ptnorm_mp,ept_mp,roots_mp)
c     
                          iocc(l) = lold
 70                    continue
                       iocc(k) = kold
 60                 continue
                    iocc(j) = jold
 50              continue
                 iocc(i) = iold
 40           continue
c     
 10        continue
c     
           junk = selci_mfree(ivc)
c     
           end
C> @}
