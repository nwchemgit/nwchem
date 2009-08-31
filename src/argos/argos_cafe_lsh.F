      subroutine argos_cafe_lsh(lself,iga,idt,istemp,
     + isfr,isto,jsfr,jsto,
     + idh,nsh,nangls,mangls,nconst,langls,itemp)
c
      implicit none
c
#include "argos_cafe_common.fh"
c
      integer nsh,mangls
      integer idh(nsh,5)
      integer iga(msa),idt(msa)
      integer langls(1)
      integer itemp(nsatot),istemp(msa)
      integer isa,jsa,isfr,isto,jsfr,jsto,isafr,isato,jsafr,jsato
      integer ish,nhi,nangls,nconst
      integer ishloc,ihloc,jhloc,khloc,iq,jq,kq
      logical lself
c
#include "bitops.fh"
c
      if(numh(2).eq.0) return
c
c     This subroutine evaluates the solute angles list
c
      isafr=isfr
      isato=isto
      jsafr=jsfr
      jsato=jsto
c
      nangls=0
      nconst=0
c
      do 1 isa=1,nsatot
      itemp(isa)=0
    1 continue
      do 2 jsa=jsafr,jsato
      itemp(iga(jsa))=4
    2 continue
      do 3 isa=isafr,isato
      itemp(iga(isa))=1
    3 continue
c
      do 4 ish=1,numh(2)
c
c     skip angles with negative ishs
c
      if(idh(ish,4).ge.0) then
c
      nhi=itemp(idh(ish,1))
      if(nhi.gt.0) then
      nhi=nhi+itemp(idh(ish,2))
      if(nhi.gt.1) then
      nhi=nhi+itemp(idh(ish,3))
      if(nhi.gt.2) then
      if((lself.and.nhi.eq.3).or.
     + (.not.lself.and.(nhi.eq.6.or.nhi.eq.9))) then
      nangls=nangls+1
      if(idh(ish,4).eq.1) nconst=nconst+1
      if(nangls.gt.mangls)
     + call md_abort('Insufficient memory for angle list',0)
      langls(nangls)=ish
      endif
      endif
      endif
      endif
      endif
    4 continue
c
      if(iqmmm.ne.0) then
      do 5 isa=1,nsatot
      itemp(isa)=0
    5 continue
      do 6 jsa=jsafr,jsato
      itemp(iga(jsa))=jsa
    6 continue
      do 7 isa=isafr,isato
      itemp(iga(isa))=isa
    7 continue
      do 8 ishloc=1,nangls
      ish=langls(ishloc)
      ihloc=itemp(idh(ish,1))
      jhloc=itemp(idh(ish,2))
      khloc=itemp(idh(ish,3))
      iq=iand(idt(ihloc),lquant)
      jq=iand(idt(jhloc),lquant)
      kq=iand(idt(khloc),lquant)
      if(iq.eq.lquant.or.jq.eq.lquant.or.kq.eq.lquant) then
      if(iq.ne.lquant) istemp(ihloc)=1
      if(jq.ne.lquant) istemp(jhloc)=1
      if(kq.ne.lquant) istemp(khloc)=1
      endif
    8 continue
      endif
c
      return
      end
