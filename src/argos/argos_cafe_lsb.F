      subroutine argos_cafe_lsb(lself,iga,idt,istemp,
     + isfr,isto,jsfr,jsto,
     + idb,nsb,nbonds,mbonds,nconst,lbonds,itemp)
c
      implicit none
c
#include "argos_cafe_common.fh"
c
      integer iga(msa),idt(msa)
      integer nsb,mbonds
      integer idb(nsb,4)
      integer lbonds(*)
      integer itemp(nsatot),istemp(msa)
      integer isa,jsa,isfr,isto,jsfr,jsto,isafr,isato,jsafr,jsato
      integer isb,nbi,nbj,nbonds,nconst,ibloc,jbloc,isbloc,iq,jq
      logical lself
c
#include "bitops.fh"
c
      if(numb(2).eq.0) return
c
c     This subroutine evaluates the solute bonds list
c
      isafr=isfr
      isato=isto
      jsafr=jsfr
      jsato=jsto
c
      nbonds=0
      nconst=0
c
      do 1 isa=1,nsatot
      itemp(isa)=0
    1 continue
      do 2 jsa=jsafr,jsato
      itemp(iga(jsa))=3
c      if(iand(idt(jsa),mdynam).ne.ldynam) itemp(iga(jsa))=-3
    2 continue
      do 3 isa=isafr,isato
      itemp(iga(isa))=1
c      if(iand(idt(isa),mdynam).ne.ldynam) itemp(iga(isa))=-1
    3 continue
c
      do 4 isb=1,numb(2)
      if(idb(isb,3).ge.0) then
      nbi=itemp(idb(isb,1))
c      lbi=nbi.ge.0
c      nbi=iabs(nbi)
      if(nbi.gt.0) then
      nbj=itemp(idb(isb,2))
c      lbj=nbj.ge.0
      nbi=nbi+iabs(nbj)
      if(nbi.gt.1) then
      if((lself.and.nbi.eq.2).or.(.not.lself.and.nbi.eq.4)) then
c      if(lbi.or.lbj) then
      nbonds=nbonds+1
      if(idb(isb,3).eq.1) nconst=nconst+1
      if(nbonds.gt.mbonds)
     + call md_abort('Insufficient memory for bond list',0)
      lbonds(nbonds)=isb
c      endif
      endif
      endif
      endif
      endif
    4 continue
c
      return
      end
c $Id$
