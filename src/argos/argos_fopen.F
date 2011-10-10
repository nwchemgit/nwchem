      subroutine argos_fopen(lclose)
c
      implicit none
c
#include "argos.fh"
c
      logical lclose
      character*5 gid
      integer iun,ibl
      logical lexist
c
c
c     Open debug file on every process if requested
c     The debug filename has the world process id
c     ---------------------------------------------
c
      if(idebug.gt.0) then
      lfndbg=18
      write(fildbg,'(a,i5.5,a)') 'nwchem_',mew,'.dbg'
      open(unit=lfndbg,file=fildbg,form='formatted',status='unknown')
      endif
c
      if(me.ne.0) return
c
c     All other files are only opened for process id 0 in each process group
c     ----------------------------------------------------------------------
c
      if(.not.lclose) then
      if(me.eq.0) then
      if(npg.gt.1) then
      write(gid,3300) meg
 3300 format(i5.5)
      ibl=index(filnam,' ')-1
      if(ibl.lt.1) then
      filnam='nwmd'
      ibl=4
      endif
      iun=index(filnam,'_')-1
      if(iun.lt.0) iun=ibl
      if(iun.eq.0) then
      ibl=index(filnam,' ')-1
      iun=5
      endif
      filtop=filnam(1:iun)//'.top'
      filrst=filnam(1:ibl)//'.rst'
      filout=filnam(1:ibl)//gid//'.out'
      filtrj=filnam(1:ibl)//gid//'.trj'
      filprp=filnam(1:ibl)//gid//'.prp'
      filmro=filnam(1:ibl)//gid//'.mro'
      filgib=filnam(1:ibl)//gid//'.gib'
      if(npg.gt.1) then
      filrst=filnam(1:ibl)//gid//'.rst'
      inquire(file=filrst(1:index(filrst,' ')-1),exist=lexist)
      if(.not.lexist) filrst=filnam(1:ibl)//'.rst'
      endif
      endif
      endif
      endif
c
      if(nfnewf.gt.0) then
      ibatch=ibatch+1
      write(filtrj,1000) root(1:index(root,' ')-1),
     + ibatch,'.trj'
      write(filhop,1000) root(1:index(root,' ')-1),
     + ibatch,'.hop'
      write(filprp,1000) root(1:index(root,' ')-1),
     + ibatch,'.prp'
      write(filpmf,1000) root(1:index(root,' ')-1),
     + ibatch,'.pmf'
      write(filtim,1000) root(1:index(root,' ')-1),
     + ibatch,'.tim'
      write(rfile,1000) root(1:index(root,' ')-1),
     + ibatch,'.rst'
 1000 format(a,i3.3,a)
      endif
c
      if(ntype.ne.3) then
      if(nfcoor.gt.0.or.nfscoo.gt.0.or.nfvelo.gt.0.or.nfsvel.gt.0) then
      if(lclose) close(lfntrj)
      open(unit=lfntrj,file=filtrj(1:index(filtrj,' ')-1),
     + form='formatted',status='unknown')
      call argos_cafe_trjhdr(lfntrj)
      endif
      if(nfhop.gt.0) then
      if(lclose) close(lfnhop)
      open(unit=lfnhop,file=filhop(1:index(filhop,' ')-1),
     + form='formatted',status='unknown')
      endif
      if(nfprop.gt.0) then
      if(lclose) close(lfnprp)
      open(unit=lfnprp,file=filprp(1:index(filprp,' ')-1),
     + form='formatted',status='unknown')
      call argos_prop_header()
      endif
      endif
      if(ntype.eq.2.and.iprpmf.ne.0) then
      if(lclose) close(lfnpmf)
      open(unit=iabs(lfnpmf),file=filpmf(1:index(filpmf,' ')-1),
     + form='formatted',status='unknown')
      endif
      if(nftime.gt.0) then
      if(lclose) close(lfntim)
      open(unit=lfntim,file=filtim(1:index(filtim,' ')-1),
     + form='formatted',status='unknown')
      call argos_hdrtim()
      endif
c
      if(.not.lclose) then
c
      if(ntype.eq.3) then
      open(unit=lfngib,file=filgib(1:index(filgib,' ')-1),
     + form='formatted',status='unknown')
      endif
      if(ntype.eq.0.and.nftri.gt.0) then
      open(unit=lfntri,file=filtri(1:index(filtri,' ')-1),
     + form='formatted',status='unknown')
      endif
c
      if(itest.gt.0) then
      open(unit=lfntst,file=filtst(1:index(filtst,' ')-1),
     + form='formatted',status='unknown')
      endif
c
      if(npg.gt.1.and.me.eq.0) then
      rfile=filnam(1:ibl)//gid//'.rst'
      endif
c
      endif
c
      return
      end
c $Id$
