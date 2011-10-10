      subroutine argos_cafe_lss(lself,lpbc,lpbcs,xs,iga,idt,igr,
     + isgm,iss,ish,isz,
     + isfr,isto,jsfr,jsto,nlocs,mpairs,npairs,lssjpt,lssin,lssj,list,
     + rwx,rw,ndxt,nt,ndxx,nx,lda,rda,lseq,lsthop,mprot)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "bitops.fh"
#include "util.fh"
#include "mafdecls.fh"
c
      real*8 rwx(mscr,3),rw(mscr)
      real*8 xs(msa,3)
      integer nlocs,nt,nx
      integer ndxt(0:nt,2),ndxx(0:nx,2)
c
      integer lssj(*),lssjpt(nlocs,2),lssin(nlocs,2)
      integer list(msa,3),iga(msa),igr(msa),idt(msa),iss(msa),isz(msa)
      integer ish(msa),isgm(msa)
      integer lseq(mseq),lsthop(2,*),mprot(mseq)
      logical lself,lpbc,lpbcs,lsame
c
      integer isafr,isfr,isato,isto,jsafr,jsfr,jsato,jsto
      integer mpairs,npairs,isa,jsa,nlist1,nlist2,idoff,iaoff
      integer ifrsx,itosx,ifrst,itost,jfrsx,jtosx,jfrst,jtost
      integer iskip1,iskip2,istsum,ist,jst,isxsum,isx,jsx,ilist
      integer nlist0,isa0,iisa,isglob,ipss,jsglob,nlist,leq,nssl0
      integer lda(16,*)
      real*8 rda(11,*)
c
      real*8 xacc(11)
      integer idon,iacc,i,idown,idfr,idto,iaown,iafr,iato
c
c
c     this subroutine evaluates the solute-solute pairlist
c
      isafr=isfr
      isato=isto
      jsafr=jsfr
      jsato=jsto
      if(lself) isato=isato-1
c
c      write(*,'(i5,a,i5)') me,' argos_cafe_lss nseq=',nseq
c      write(*,'(20i5)') (lseq(i),i=1,nseq)
c
      npairs=0
c
c     for lstype=0 the pairlist is atom based
c     ---------------------------------------
c
      if(lstype.eq.0) then
      do 1 isa=isafr,isato
      if(lself) then
      jsafr=isa+1
      jsato=isato+1
      endif
      do 2 jsa=jsafr,jsato
      rwx(jsa,1)=xs(isa,1)-xs(jsa,1)
      rwx(jsa,2)=xs(isa,2)-xs(jsa,2)
      rwx(jsa,3)=xs(isa,3)-xs(jsa,3)
    2 continue
      if(lpbc.or.lpbcs)
     + call argos_cafe_pbc(1,rwx,mscr,rwx,mscr,0,jsafr,jsato)
      do 3 jsa=jsafr,jsato
      rw(jsa)=rwx(jsa,1)**2+rwx(jsa,2)**2+rwx(jsa,3)**2
    3 continue
      nlist1=0
      nlist2=0
      isglob=iga(isa)
      ifrsx=ndxx(isglob-1,1)+1
      itosx=ndxx(isglob,1)
      ifrst=ndxt(isglob-1,1)+1
      itost=ndxt(isglob,1)
      do 4 jsa=jsafr,jsato
      iskip1=1
      iskip2=1
      if(rw(jsa).lt.rshrt2.and.
     + (iand(idt(isa),mdynam).eq.ldynam.or.includ.eq.1.or.
     + iand(idt(jsa),mdynam).eq.ldynam)) iskip1=0
      if(rw(jsa).lt.rlong2.and.
     + (iand(idt(isa),mdynam).eq.ldynam.or.includ.eq.1.or.
     + iand(idt(jsa),mdynam).eq.ldynam)) iskip2=0
      if(rw(jsa).lt.rrest2.and.
     + (iand(idt(isa),mrestr).eq.lrestr.or.
     + iand(idt(jsa),mrestr).eq.lrestr) ) iskip1=0
      istsum=0
      isxsum=0
      if(iskip1.eq.0.or.iskip2.eq.0) then
      jsglob=iga(jsa)
      jfrsx=ndxx(jsglob-1,1)+1
      jtosx=ndxx(jsglob,1)
      jfrst=ndxt(jsglob-1,1)+1
      jtost=ndxt(jsglob,1)
      istsum=1
      do 5 ist=ifrst,itost
      if(jsglob.eq.ndxt(ist,2)) goto 7
    5 continue
      do 6 jst=jfrst,jtost
      if(isglob.eq.ndxt(jst,2)) goto 7
    6 continue
      istsum=0
    7 continue
      isxsum=0
      if(istsum.eq.0) then
      isxsum=1
      do 8 isx=ifrsx,itosx
      if(jsglob.eq.ndxx(isx,2)) goto 10
    8 continue
      do 9 jsx=jfrsx,jtosx
      if(isglob.eq.ndxx(jsx,2)) goto 10
    9 continue
      isxsum=0
   10 continue
      endif
      endif
      if(istsum.eq.0.and.isxsum.eq.0) then
      if(iskip1.eq.0) then
      nlist1=nlist1+1
      list(nlist1,1)=jsa
      else
      if(iskip2.eq.0) then
      nlist2=nlist2+1
      list(nlist2,2)=jsa
      endif
      endif
      endif
    4 continue
      if(npairs+nlist1+nlist2.gt.mpairs)
     + call md_abort('Insufficient memory for pairlist',0)
      do 11 ilist=1,nlist1
      lssj(npairs+ilist)=list(ilist,1)
   11 continue
      lssjpt(isa-isafr+1,1)=npairs+1
      npairs=npairs+nlist1
      lsss=lsss+nlist1
      lssin(isa-isafr+1,1)=nlist1
      if(npss.eq.2) then
      do 12 ilist=1,nlist2
      lssj(npairs+ilist)=list(ilist,2)
   12 continue
      lssjpt(isa-isafr+1,2)=npairs+1
      npairs=npairs+nlist2
      llss=llss+nlist2
      lssin(isa-isafr+1,2)=nlist2
      endif
    1 continue
      endif
c
c     for lstype=1 the pairlist is charge group based
c         lstype=2 the pairlist is segment based
c
      if(lstype.eq.1.or.lstype.eq.2) then
      isa=isafr
      isa0=isafr
   13 continue
      if(lself) then
      jsafr=isa+1
      jsato=isato+1
      endif
      do 14 jsa=jsafr,jsato
      list(jsa,1)=0
      list(jsa,2)=0
   14 continue
   15 continue
      do 16 jsa=jsafr,jsato
      rwx(jsa,1)=xs(isa,1)-xs(jsa,1)
      rwx(jsa,2)=xs(isa,2)-xs(jsa,2)
      rwx(jsa,3)=xs(isa,3)-xs(jsa,3)
   16 continue
      if(lpbc.or.lpbcs)
     + call argos_cafe_pbc(1,rwx,mscr,rwx,mscr,0,jsafr,jsato)
      do 17 jsa=jsafr,jsato
      rw(jsa)=rwx(jsa,1)**2+rwx(jsa,2)**2+rwx(jsa,3)**2
   17 continue
      if(iand(idt(isa),mdynam).eq.ldynam.or.includ.eq.1) then
      do 18 jsa=jsafr,jsato
      if(rw(jsa).lt.rshrt2) list(jsa,1)=1
   18 continue
      endif
      if(iand(idt(isa),mrestr).eq.lrestr) then
      do 20 jsa=jsafr,jsato
      if(iand(idt(jsa),mdynam).eq.ldynam.and.rw(jsa).lt.rshrt2)
     + list(jsa,1)=1
      if(iand(idt(jsa),mdynam).ne.ldynam.and.rw(jsa).lt.rrest2)
     + list(jsa,1)=1
   20 continue
      endif
      if(iand(idt(isa),mfixed).eq.lfixed) then
      do 22 jsa=jsafr,jsato
      if(iand(idt(jsa),mdynam).eq.ldynam.and.rw(jsa).lt.rshrt2)
     + list(jsa,1)=1
      if(iand(idt(jsa),mrestr).eq.lrestr.and.rw(jsa).lt.rrest2)
     + list(jsa,1)=1
   22 continue
      endif
c
      if(npss.eq.2) then
      if(iand(idt(isa),mdynam).eq.ldynam.or.includ.eq.1) then
      do 26 jsa=jsafr,jsato
      if(rw(jsa).lt.rlong2) list(jsa,2)=1
   26 continue
      endif
      if(iand(idt(isa),mrestr).eq.lrestr) then
      do 28 jsa=jsafr,jsato
      if(iand(idt(jsa),mdynam).eq.ldynam.and.rw(jsa).lt.rlong2)
     + list(jsa,2)=1
   28 continue
      endif
      if(iand(idt(isa),mfixed).eq.lfixed) then
      do 30 jsa=jsafr,jsato
      if(iand(idt(jsa),mdynam).eq.ldynam.and.rw(jsa).lt.rlong2)
     + list(jsa,2)=1
   30 continue
cx      endif
c      if(nhop.gt.0) then
c      if(ish(isa).gt.0.and.iand(ish(isa),1).ne.0) then
c      do 426 jsa=jsafr,jsato
c      list(jsa,1)=0
c  426 continue
c      endif
c      do 427 jsa=jsafr,jsato
c      if(ish(jsa).gt.0.and.iand(ish(jsa),1).ne.0)
c     + list(jsa,1)=0
cx  427 continue
      endif
      endif
      if(isa.lt.isato) then
      if(lstype.eq.1.and.igr(isa+1).eq.igr(isa)) then
      isa=isa+1
      goto 15
      endif
c      if(lstype.eq.2.and.isl(isa+1,lssgm).eq.isl(isa,lssgm)) then
c      isa=isa+1
c      goto 15
c      endif
      endif
c
      do 32 jsa=isa0,isa
      if(iand(idt(jsa),mrestr).ne.lrestr) goto 33
   32 continue
      goto 34
   33 continue
c
c     solute atoms in the same charge group (lstype=1) or
c     segment (lstype=2) always all interact in short range list
c     this can only happen when lself is true
c
      if(lself) then
      do 35 jsa=isa0+1,isa
      if(iand(idt(isa),mdynam).eq.ldynam.or.includ.eq.1.or.
     + iand(idt(jsa),mdynam).eq.ldynam) list(jsa,1)=1
   35 continue
      endif
   34 continue
c
      if(lstype.eq.1) then
c
c     following two loops make sure that all atoms in a charge group
c     are counted if at least one atom is within short range cutoff
c
      do 36 jsa=jsafr+1,jsato
      if(igr(jsa).eq.igr(jsa-1).and.list(jsa-1,1).eq.1)
     + list(jsa,1)=1
   36 continue
      do 37 jsa=jsato-1,jsafr,-1
      if(igr(jsa).eq.igr(jsa+1).and.list(jsa+1,1).eq.1)
     + list(jsa,1)=1
   37 continue
c
c
      if(npss.eq.2) then
c
c     remove all short range interaction from long range list
c
      do 42 jsa=jsafr,jsato
      if(list(jsa,1).eq.1) list(jsa,2)=0
   42 continue
c
c     following two loops make sure that all atoms in a charge group
c     are counted if at least one atom is within long range cutoff
c
      do 43 jsa=jsafr+1,jsato
      if(igr(jsa).eq.igr(jsa-1).and.list(jsa-1,2).eq.1)
     + list(jsa,2)=1
   43 continue
      do 44 jsa=jsato-1,jsafr,-1
      if(igr(jsa).eq.igr(jsa+1).and.list(jsa+1,2).eq.1)
     + list(jsa,2)=1
   44 continue
c
      endif
      endif
c      if(lstype.eq.2) then
cc
cc     following two loops make sure that all atoms in a segment
cc     are counted if at least one atom is within short range cutoff
cc
c      do 43 jsa=jsafr+1,jsato
c      if(isl(jsa,lssgm).eq.isl(jsa-1,lssgm).and.list(jsa-1,1).eq.1)
c     + list(jsa,1)=1
c   43 continue
c      do 44 jsa=jsato-1,jsafr,-1
c      if(isl(jsa,lssgm).eq.isl(jsa+1,lssgm).and.list(jsa+1,1).eq.1)
c     + list(jsa,1)=1
c   44 continue
cc
c      if(npss.eq.2) then
cc
cc     following two loops make sure that all atoms in a charge group
cc     are counted if at least one atom is within long range cutoff
cc
c      do 45 jsa=jsafr+1,jsato
c      if(isl(jsa,lssgm).eq.isl(jsa-1,lssgm).and.list(jsa-1,2).eq.1)
c     + list(jsa,2)=1
c   45 continue
c      do 46 jsa=jsato-1,jsafr,-1
c      if(isl(jsa,lssgm).eq.isl(jsa+1,lssgm).and.list(jsa+1,2).eq.1)
c     + list(jsa,2)=1
c   46 continue
cc
c      endif
c      endif
      nlist0=0
      do 47 iisa=isa0,isa
c
c     remove third neighbor pairs and excluded pairs from the lists
c
      isglob=iga(iisa)
c
      do 48 ipss=1,npss
c
      do 49 jsa=jsafr,jsato
      list(jsa,3)=list(jsa,ipss)
   49 continue
c
      do 50 jsa=jsafr,jsato
      if(iga(iisa).eq.iga(jsa)) list(jsa,3)=0
      if(lself.and.jsa.lt.iisa) list(jsa,3)=0
   50 continue
      ifrst=ndxt(isglob-1,1)+1
      itost=ndxt(isglob,1)
      do 51 ist=ifrst,itost
      jsglob=ndxt(ist,2)
      do 52 jsa=jsafr,jsato
      if(iga(jsa).eq.jsglob) list(jsa,3)=0
   52 continue
   51 continue
      do 53 jsa=jsafr,jsato
      if(list(jsa,3).eq.1) then
      jsglob=iga(jsa)
      jfrst=ndxt(jsglob-1,1)+1
      jtost=ndxt(jsglob,1)
      do 54 jst=jfrst,jtost
      if(isglob.eq.ndxt(jst,2)) list(jsa,3)=0
   54 continue
      endif
   53 continue
      ifrsx=ndxx(isglob-1,1)+1
      itosx=ndxx(isglob,1)
      do 55 isx=ifrsx,itosx
      jsglob=ndxx(isx,2)
      do 56 jsa=jsafr,jsato
      if(iga(jsa).eq.jsglob) list(jsa,3)=0
   56 continue
   55 continue
      do 57 jsa=jsafr,jsato
      if(list(jsa,3).eq.1) then
      jsglob=iga(jsa)
      jfrsx=ndxx(jsglob-1,1)+1
      jtosx=ndxx(jsglob,1)
      do 58 jsx=jfrsx,jtosx
      if(isglob.eq.ndxx(jsx,2)) list(jsa,3)=0
   58 continue
      endif
   57 continue
c
c     remove interactions between dummy atoms
c
c     bits set in iss : 0 set 1 s-w is dummy               (1)
c                       1 set 2 s-w is dummy               (2)
c                       2 set 3 s-w is dummy               (4)
c                       3 set 1 s-s non-self is dummy      (8)
c                       4 set 2 s-s non-self is dummy     (16)
c                       5 set 3 s-s non-self is dummy     (32)
c                       6 set 1 s-s self is dummy         (64)
c                       7 set 2 s-s self is dummy        (128)
c                       8 set 3 s-s self is dummy        (256)
c
cz      if(.not.ithint) then
cz      do 150 jsa=jsafr,jsato
cz      if(iand(iss(iisa),72).eq.72) list(jsa,3)=0
cz      if(iand(iss(jsa),72).eq.72) list(jsa,3)=0
cz  150 continue
cz      else
cz      do 151 jsa=jsafr,jsato
cz      if(iand(iss(iisa),144).eq.144.and.iand(iss(iisa),288).eq.288)
cz     + list(jsa,3)=0
cz      if(iand(iss(jsa),144).eq.144.and.iand(iss(jsa),288).eq.288)
cz     + list(jsa,3)=0
cz      if(iand(iss(iisa),144).eq.144.and.iand(iss(jsa),288).eq.288)
cz     + list(jsa,3)=0
cz      if(iand(iss(jsa),144).eq.144.and.iand(iss(iisa),288).eq.288)
cz     + list(jsa,3)=0
cz  151 continue
cz      endif
c
c     make the actual pairlist
c
      nlist=0
      do 59 jsa=jsafr,jsato
      if(list(jsa,3).ne.0) then
c      write(*,'(a,8i5,2(1x,9l1))') 'LSS ',iisa,jsa,iga(iisa),iga(jsa),
c     + idt(iisa),idt(jsa),iss(iisa),iss(jsa),
c     + iand(iss(iisa),1).ne.0,iand(iss(iisa),2).ne.0,
c     + iand(iss(iisa),4).ne.0,iand(iss(iisa),8).ne.0,
c     + iand(iss(iisa),16).ne.0,iand(iss(iisa),32).ne.0,
c     + iand(iss(iisa),64).ne.0,iand(iss(iisa),128).ne.0,
c     + iand(iss(iisa),256).ne.0,
c     + iand(iss(jsa),1).ne.0,iand(iss(jsa),2).ne.0,
c     + iand(iss(jsa),4).ne.0,iand(iss(jsa),8).ne.0,
c     + iand(iss(jsa),16).ne.0,iand(iss(jsa),32).ne.0,
c     + iand(iss(jsa),64).ne.0,iand(iss(jsa),128).ne.0,
c     + iand(iss(jsa),256).ne.0
      nlist=nlist+1
      list(nlist,3)=jsa
      endif
   59 continue
      if(npairs+nlist.gt.mpairs)
     + call md_abort('Insufficient memory for pairlist',0)
      leq=0
c
c     check if this pairlist section was generated before
c
      if(ipss.eq.2.and.nlist.eq.nlist0) then
      nssl0=lssjpt(iisa-1,ipss)-1
      do 60 ilist=1,nlist
      if(lssj(nssl0+ilist).ne.list(ilist,3)) leq=1
   60 continue
      endif
c
      if(leq.eq.0) then
      if(npairs+nlist.gt.mpairs)
     + call md_abort('Insufficient memory for pairlist',0)
      do 61 ilist=1,nlist
      lssj(npairs+ilist)=list(ilist,3)
   61 continue
      lssjpt(iisa-isfr+1,ipss)=npairs+1
      npairs=npairs+nlist
      if(ipss.eq.1) then
      lsss=lsss+nlist
      else
      llss=llss+nlist
      endif
      else
      lssjpt(iisa-isfr+1,ipss)=lssjpt(iisa-isfr,ipss)
      endif
c
      lssin(iisa-isfr+1,ipss)=nlist
   48 continue
   47 continue
      isa=isa+1
      isa0=isa
      if(isa.le.isato) goto 13
      endif
c
      if(lself) then
      lssin(isto-isfr+1,1)=0
      lssin(isto-isfr+1,2)=0
      endif
c
      return
      end
c $Id$
