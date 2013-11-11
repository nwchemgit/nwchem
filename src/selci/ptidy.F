C> \ingroup selci
C> @{
      subroutine selci_ptidy(ioconf, nintpo, noconf, ncold, nconmx,
     $     nroot, ept, enew, irange, ptnorm, ptnorm_mp, ept_mp)
*
* $Id$
*
#include "implicit.fh"
#include "errquit.fh"
#include "global.fh"
#include "stdio.fh"
      dimension ioconf(nintpo,nconmx), ept(nroot), enew(nroot),
     &     ptnorm(nroot), irange(21)
      dimension ept_mp(nroot), ptnorm_mp(nroot)
c
      parameter (nnn = 20)
      dimension eptn(nnn), enewn(nnn), ptnormn(nnn), irangen(21)
      dimension eptn_mp(nnn), ptnormn_mp(nnn)
      integer pmap(0:1023)        ! 1024 max process count
      data isync/1/
c
c     collect info from the other processes at end of selection
c
c
      if (nroot.gt.nnn) call errquit('PTIDY: hard dimension fail',nroot,
     &       CALC_ERR)
c
      if (ga_nnodes() .gt. 1024) call errquit('ptidy: too many procs',0,
     &       CALC_ERR)
      call ga_list_nodeid(pmap, ga_nnodes()) ! Map to msg passing layer
*      write(6,*) ga_nodeid(),' processor map ', 
*     $     (pmap(i),i=0,ga_nnodes()-1)
      node0 = pmap(0)
c
      if (ga_nodeid().gt.0) then
c         write(6,*) nodeid(),noconf, ncold, nintpo, nroot
c         write(6,*) nodeid(),ept
c         write(6,*) nodeid(),ept_mp
c         write(6,*) nodeid(),enew
c         write(6,*) nodeid(),irange
c         write(6,*) nodeid(),ptnorm
c         write(6,*) nodeid(),ptnorm_mp
         call snd(21, noconf-ncold, mitob(1), node0, isync)
         call snd(22, ioconf(1,ncold+1),mitob(nintpo*(noconf-ncold)),
     $        node0, isync)
         call snd(23, ept, mdtob(nroot), node0, isync)
         call snd(24, enew, mdtob(nroot), node0, isync)
         call snd(25, ptnorm, mdtob(nroot), node0, isync)
         call snd(26, irange, mitob(21), node0, isync)
         call snd(27, ptnorm_mp, mdtob(nroot), node0, isync)
         call snd(28, ept_mp, mdtob(nroot), node0, isync)
      else
c        write(6,*) nodeid(),noconf, ncold, nintpo, nroot
c        write(6,*) nodeid(),ept
c        write(6,*) nodeid(),ept_mp
c        write(6,*) nodeid(),enew
c        write(6,*) nodeid(),irange
c        write(6,*) nodeid(),ptnorm
c        write(6,*) nodeid(),ptnorm_mp
c        call util_flush(6)
         ipt = noconf+1
         do 10 iislave=1,ga_nnodes()-1
            islave = pmap(iislave)
            call rcv(21, new, mitob(1), len, islave, node, isync)
            noconf = noconf + new
            if (noconf.le.nconmx) then
               call rcv(22, ioconf(1,ipt), mitob(nintpo*new), len,
     $              islave, node, isync)
               ipt = ipt + new
            else
               call errquit('ptidy: exceeded nconmx ',noconf,
     &       UNKNOWN_ERR)
            endif
            call rcv(23, eptn, mdtob(nroot), len, islave, node, isync)
            call rcv(24, enewn, mdtob(nroot), len, islave, node, isync)
            call rcv(25, ptnormn, mdtob(nroot), len, islave, node,
     &               isync)
            call rcv(26, irangen, mitob(21), len, islave, node, isync)
            call rcv(27, ptnormn_mp, mdtob(nroot), len, islave, node,
     &           isync)
            call rcv(28, eptn_mp, mdtob(nroot), len, islave, node,
     &           isync)
            do 20 i = 1,nroot
               ept(i) = ept(i) + eptn(i)
               enew(i) = enew(i) + enewn(i)
               ptnorm(i) = ptnorm(i) + ptnormn(i)
               ept_mp(i) = ept_mp(i) + eptn_mp(i)
               ptnorm_mp(i) = ptnorm_mp(i) + ptnormn_mp(i)
 20         continue
            do 30 i = 1,21
               irange(i) = irange(i) + irangen(i)
 30         continue
 10      continue
c        write(6,*) nodeid(),noconf, ncold, nintpo, nroot
c        write(6,*) nodeid(),ept
c        write(6,*) nodeid(),enew
c        write(6,*) nodeid(),irange
c        write(6,*) nodeid(),ptnorm
c        call util_flush(6)
      endif
c
      end
C> @}
