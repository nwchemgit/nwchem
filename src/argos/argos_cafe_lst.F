      subroutine argos_cafe_lst(lself,iga,idt,istemp,
     + isfr,isto,jsfr,jsto,
     + idtn,nst,nthird,mthird,lthird,itemp)
c
      implicit none
c
#include "argos_cafe_common.fh"
c
      integer nst,mthird
      integer idtn(0:nst,2)
      integer iga(msa),idt(msa)
c      integer idst(mst),jdst(mst)
      integer lthird(1)
c      integer isl(msa,mis2)
      integer isa,jsa,isfr,isto,jsfr,jsto,isafr,isato,jsafr,jsato
      integer ist,nti,nthird,itloc,jtloc,istloc,iq,jq
      integer itemp(nsatot),istemp(msa)
      logical lself
c
      nthird=0
c
      if(numt(2).eq.0) return
c
c     This subroutine evaluates the solute third neighbor list
c
      isafr=isfr
      isato=isto
      jsafr=jsfr
      jsato=jsto
c
      do 1 isa=1,nsatot
      itemp(isa)=0
    1 continue
      do 2 jsa=jsafr,jsato
      itemp(iga(jsa))=3
    2 continue
      do 3 isa=isafr,isato
      itemp(iga(isa))=1
    3 continue
c
      do 4 isa=1,nsatot
      if(idtn(isa,1).gt.idtn(isa-1,1)) then
      do 5 ist=idtn(isa-1,1)+1,idtn(isa,1)
      nti=itemp(isa)
      if(nti.gt.0) then
      nti=nti+itemp(idtn(ist,2))
      if(nti.gt.1) then
      if((lself.and.nti.eq.2).or.(.not.lself.and.nti.eq.4)) then
      nthird=nthird+1
      if(nthird.gt.mthird)
     + call md_abort('Insufficient memory for third neighbor list',0)
      lthird(nthird)=ist
      endif
      endif
      endif
    5 continue
      endif
    4 continue
c
      return
      end
c $Id$
