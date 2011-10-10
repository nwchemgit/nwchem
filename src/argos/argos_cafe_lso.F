      subroutine argos_cafe_lso(lself,iga,idt,istemp,
     + isfr,isto,jsfr,jsto,
     + ido,nso,nimprs,mimprs,nconst,limprs,itemp)
c
      implicit none
c
#include "argos_cafe_common.fh"
c
      integer nso,mimprs
      integer ido(nso,6)
      integer iga(msa),idt(msa)
      integer limprs(1)
      integer itemp(nsatot),istemp(msa)
      integer isa,jsa,isfr,isto,jsfr,jsto,isafr,isato,jsafr,jsato
      integer iso,noi,nimprs,nconst
      integer isoloc,ioloc,joloc,koloc,loloc,iq,jq,kq,lq
      logical lself
c
#include "bitops.fh"
c
      if(numi(2).eq.0) return
c
c     This subroutine evaluates the solute impropers list
c
      isafr=isfr
      isato=isto
      jsafr=jsfr
      jsato=jsto
c
      nimprs=0
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
      do 4 iso=1,numi(2)
c
c     skip dihedrals with negative isos
c
      if(ido(iso,5).ge.0) then
c
      noi=itemp(ido(iso,1))
      if(noi.gt.0) then
      noi=noi+itemp(ido(iso,2))
      if(noi.gt.1) then
      noi=noi+itemp(ido(iso,3))
      if(noi.gt.2) then
      noi=noi+itemp(ido(iso,4))
      if(noi.gt.3) then
      if((lself.and.noi.eq.4).or.
     + (.not.lself.and.(noi.eq.8.or.noi.eq.12.or.noi.eq.16))) then
      nimprs=nimprs+1
      if(ido(iso,5).eq.1) nconst=nconst+1
      if(nimprs.gt.mimprs)
     + call md_abort('Insufficient memory for improper list',0)
      limprs(nimprs)=iso
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
