      subroutine argos_cafe_lsd(lself,iga,idt,istemp,
     + isfr,isto,jsfr,jsto,
     + idd,nsd,ndihes,mdihes,nconst,ldihes,itemp)
c
      implicit none
c
#include "argos_cafe_common.fh"
c
      integer nsd,mdihes
      integer idd(nsd,6)
      integer iga(msa),idt(msa)
      integer ldihes(1)
      integer itemp(nsatot),istemp(msa)
      integer isa,jsa,isfr,isto,jsfr,jsto,isafr,isato,jsafr,jsato
      integer isd,ndi,ndihes,nconst
      integer isdloc,idloc,jdloc,kdloc,ldloc,iq,jq,kq,lq
      logical lself
c
#include "bitops.fh"
c
      if(numd(2).eq.0) return
c
c     This subroutine evaluates the solute dihedrals list
c
      isafr=isfr
      isato=isto
      jsafr=jsfr
      jsato=jsto
c
      ndihes=0
      nconst=0
c
      do 1 isa=1,nsatot
      itemp(isa)=0
    1 continue
      do 2 jsa=jsafr,jsato
      itemp(iga(jsa))=5
    2 continue
      do 3 isa=isafr,isato
      itemp(iga(isa))=1
    3 continue
c
      do 4 isd=1,numd(2)
c
c     skip dihedrals with negative isds
c
      if(idd(isd,5).ge.0) then
c
      ndi=itemp(idd(isd,1))
      if(ndi.gt.0) then
      ndi=ndi+itemp(idd(isd,2))
      if(ndi.gt.1) then
      ndi=ndi+itemp(idd(isd,3))
      if(ndi.gt.2) then
      ndi=ndi+itemp(idd(isd,4))
      if(ndi.gt.3) then
      if((lself.and.ndi.eq.4).or.
     + (.not.lself.and.(ndi.eq.8.or.ndi.eq.12.or.ndi.eq.16))) then
      ndihes=ndihes+1
      if(idd(isd,5).eq.1) nconst=nconst+1
      if(ndihes.gt.mdihes)
     + call md_abort('Insufficient memory for torsion list',0)
      ldihes(ndihes)=isd
      endif
      endif
      endif
      endif
      endif
      endif
    4 continue
c
      return
      end
c $Id$
