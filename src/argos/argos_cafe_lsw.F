      subroutine argos_cafe_lsw(lpbc,lpbcs,xs,isdt,isgr,isfr,isto,
     + xwm,iwdt,iwz,
     + iwfr,iwto,ndim,mpairs,npairs,lswjpt,lswin,lswj,list,rwx,rw)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "bitops.fh"
c
      real*8 rwx(mscr,3),rw(mscr)
      real*8 xwm(mwm,3),xs(msa,3)
      integer ndim
      integer lswj(*),lswjpt(ndim,2),lswin(ndim,2)
      integer list(mwm,2),iwdt(mwm),isdt(msa),isgr(msa)
      integer iwz(mwm)
c
      integer iwmfr,iwfr,iwmto,iwto,isafr,isfr,isato,isto
      integer mpairs,npairs,isa,iwm,nlist,ilist,isa0,iisa
c
      integer nswi1,nswi2
      logical lpbc,lpbcs
c
c     this subroutine evaluates the solute-solvent pairlist
c
      iwmfr=iwfr
      iwmto=iwto
      isafr=isfr
      isato=isto
c
c      maxlst=2*mavail
c
      npairs=0
c
c     for lstype=0 the pairlist is atom based
c
      if(lstype.eq.0) then
      do 1 isa=isafr,isato
      do 3 iwm=iwmfr,iwmto
      rwx(iwm,1)=xs(isa,1)-xwm(iwm,1)
      rwx(iwm,2)=xs(isa,2)-xwm(iwm,2)
      rwx(iwm,3)=xs(isa,3)-xwm(iwm,3)
    3 continue
      if(lpbc.or.lpbcs)
     + call argos_cafe_pbc(1,rwx,mscr,rwx,mscr,0,iwmfr,iwmto)
      do 8 iwm=iwmfr,iwmto
      rw(iwm)=rwx(iwm,1)**2+rwx(iwm,2)**2+rwx(iwm,3)**2
    8 continue
      do 4 iwm=iwmfr,iwmto
      list(iwm,1)=0
      list(iwm,2)=0
      if(rw(iwm).lt.rshrt2.and.
     + (iand(iwdt(iwm),mdynam).eq.ldynam.or.
     + iand(isdt(isa),mdynam).eq.ldynam)) list(iwm,1)=1
      if(rw(iwm).lt.rrest2.and.
     + (iand(iwdt(iwm),mrestr).eq.lrestr.or.
     + iand(isdt(isa),mrestr).eq.lrestr) ) list(iwm,1)=1
    4 continue
      if(npsw.eq.2) then
      do 401 iwm=iwmfr,iwmto
      if(rw(iwm).lt.rlong2.and.
     + (iand(iwdt(iwm),mdynam).eq.ldynam.or.
     + iand(isdt(isa),mdynam).eq.ldynam)) list(iwm,2)=1
      if(list(iwm,1).eq.1) list(iwm,2)=0
  401 continue
      endif
c
c     short range pairlist
c
      nlist=0
      do 402 iwm=iwmfr,iwmto
      if(list(iwm,1).eq.1) then
      nlist=nlist+1
      list(nlist,1)=iwm
      endif
  402 continue
      if(npairs+nlist.gt.mpairs)
     + call md_abort('Insufficient memory for pairlist',0)
      do 5 ilist=1,nlist
      lswj(npairs+ilist)=list(ilist,1)
    5 continue
      lswjpt(isa-isafr+1,1)=npairs+1
      npairs=npairs+nlist
      lssw=lssw+nlist
      lswin(isa-isafr+1,1)=nlist
c
c     long range pairlist
c
      if(npsw.eq.2) then
      nlist=0
      do 8402 iwm=iwmfr,iwmto
      if(list(iwm,2).eq.1) then
      nlist=nlist+1
      list(nlist,2)=iwm
      endif
 8402 continue
      if(npairs+nlist.gt.mpairs)
     + call md_abort('Insufficient memory for pairlist',0)
      do 8005 ilist=1,nlist
      lswj(npairs+ilist)=list(ilist,2)
 8005 continue
      lswjpt(isa-isafr+1,2)=npairs+1
      npairs=npairs+nlist
      llsw=llsw+nlist
      lswin(isa-isafr+1,2)=nlist
      endif
    1 continue
      endif
c
c     for lstype=1 the pairlist is charge group based
c         lstype=1 the pairlist is segment based
c
      if(lstype.eq.1.or.lstype.eq.2) then
      isa=isafr
      isa0=isafr
   21 continue
      nlist=0
      do 15 iwm=iwmfr,iwmto
      list(iwm,1)=0
      list(iwm,2)=0
   15 continue
   10 continue
      do 11 iwm=iwmfr,iwmto
      rwx(iwm,1)=xs(isa,1)-xwm(iwm,1)
      rwx(iwm,2)=xs(isa,2)-xwm(iwm,2)
      rwx(iwm,3)=xs(isa,3)-xwm(iwm,3)
   11 continue
      if(lpbc.or.lpbcs)
     + call argos_cafe_pbc(1,rwx,mscr,rwx,mscr,0,iwmfr,iwmto)
      do 13 iwm=iwmfr,iwmto
      rw(iwm)=rwx(iwm,1)**2+rwx(iwm,2)**2+rwx(iwm,3)**2
   13 continue
      do 9012 iwm=iwmfr,iwmto
      if(iand(iwdt(iwm),mdynam).eq.ldynam.and.
     + rw(iwm).lt.rshrt2) list(iwm,1)=1
      if(iand(iwdt(iwm),mrestr).eq.lrestr.and.
     + rw(iwm).lt.rrest2) list(iwm,1)=1
 9012 continue
      if(iand(isdt(isa),mdynam).eq.ldynam) then
      do 9013 iwm=iwmfr,iwmto
      if(rw(iwm).lt.rshrt2) list(iwm,1)=1
 9013 continue
      endif
      if(iand(isdt(isa),mrestr).eq.lrestr) then
      do 9014 iwm=iwmfr,iwmto
      if(rw(iwm).lt.rrest2) list(iwm,1)=1
 9014 continue
      endif
      if(npsw.eq.2) then
      do 9812 iwm=iwmfr,iwmto
      if(iand(iwdt(iwm),mdynam).eq.ldynam.and.
     + rw(iwm).lt.rlong2) list(iwm,2)=1
 9812 continue
      if(iand(isdt(isa),mdynam).eq.ldynam) then
      do 9813 iwm=iwmfr,iwmto
      if(rw(iwm).lt.rlong2) list(iwm,2)=1
 9813 continue
      endif
      endif
      do 9815 iwm=iwmfr,iwmto
      if(list(iwm,1).eq.1) list(iwm,2)=0
 9815 continue
      if(isa.eq.isato) goto 16
      if(lstype.eq.1) then
      if(isgr(isa+1).eq.isgr(isa)) then
      isa=isa+1
      goto 10
      endif
      endif
c      if(lstype.eq.2) then
c      if(isl(isa+1,lssgm).eq.isl(isa,lssgm)) then
c      isa=isa+1
c      goto 10
c      endif
c      endif
   16 continue
c
c     at this point the partial list for a solute charge group or
c     segment with all solvent has been evaluated
c
c     short range pairlist
c
      nlist=0
      do 17 iwm=iwmfr,iwmto
      if(list(iwm,1).eq.1) then
      nlist=nlist+1
      list(nlist,1)=iwm
      endif
   17 continue
      do 18 iisa=isa0,isa
      lswjpt(iisa-isafr+1,1)=npairs+1
      lswin(iisa-isafr+1,1)=nlist
   18 continue
      if(npairs+nlist.gt.mpairs)
     + call md_abort('Insufficient memory for pairlist',0)
      do 19 ilist=1,nlist
      lswj(npairs+ilist)=list(ilist,1)
   19 continue
      npairs=npairs+nlist
      lssw=lssw+nlist
c
c     long range pairlist
c
      if(npsw.eq.2) then
      nlist=0
      do 8017 iwm=iwmfr,iwmto
      if(list(iwm,2).eq.1) then
      nlist=nlist+1
      list(nlist,2)=iwm
      endif
 8017 continue
      do 8018 iisa=isa0,isa
      lswjpt(iisa-isafr+1,2)=npairs+1
      lswin(iisa-isafr+1,2)=nlist
 8018 continue
      if(npairs+nlist.gt.mpairs)
     + call md_abort('Insufficient memory for pairlist',0)
      do 8019 ilist=1,nlist
      lswj(npairs+ilist)=list(ilist,2)
 8019 continue
      npairs=npairs+nlist
      llsw=llsw+nlist
      endif
c
      isa=isa+1
      isa0=isa
      if(isa.le.isato) goto 21
      endif
      nswi1=0
      nswi2=0
      do 9191 isa=isafr,isato
      nswi1=nswi1+lswin(isa-isafr+1,1)
      if(npsw.eq.2) nswi2=nswi2+lswin(isa-isafr+1,2)
 9191 continue
c
      return
      end
c $Id$
