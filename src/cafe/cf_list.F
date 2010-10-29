      subroutine cf_lww(lself,lpbc,xwm,idt,iwfr,iwto,jwfr,jwto,
     + nlocw,mpairs,npairs,lwwjpt,lwwin,lwwj,list,rwx,rw)
c $Id$
      implicit none
c
#include "cf_common.fh"
#include "bitops.fh"
#include "mafdecls.fh"
c
      integer nlocw
      real*8 rwx(mscr,3),rw(mscr)
      real*8 xwm(mwm,3)
      integer lwwj(*),lwwjpt(nlocw,2)
      integer lwwin(nlocw,2)
      integer list(mwm,2),idt(mwm)
      logical lself,lpbc
c
      integer iwm,jwm,iwmfr,iwfr,iwmto,iwto,jwmfr,jwfr,jwmto,jwto
      integer ilist,nlist,mpairs,npairs
c
      iwmfr=iwfr
      if(lself) then
      lwwin(nlocw,1)=0
      lwwin(nlocw,2)=0
      lwwjpt(nlocw,1)=0
      lwwjpt(nlocw,2)=0
      iwmto=iwto-1
      jwmfr=0
      jwmto=iwto
      else
      iwmto=iwto
      jwmfr=jwfr
      jwmto=jwto
      endif
c
      npairs=0
      do 1 iwm=iwmfr,iwmto
      if(lself) jwmfr=iwm+1
c
      do 3 jwm=jwmfr,jwmto
      rwx(jwm,1)=xwm(iwm,1)-xwm(jwm,1)
      rwx(jwm,2)=xwm(iwm,2)-xwm(jwm,2)
      rwx(jwm,3)=xwm(iwm,3)-xwm(jwm,3)
    3 continue
      if(lpbc) call cf_pbc(1,rwx,mscr,rwx,mscr,0,jwmfr,jwmto)
      do 8 jwm=jwmfr,jwmto
      rw(jwm)=rwx(jwm,1)**2+rwx(jwm,2)**2+rwx(jwm,3)**2
    8 continue
      do 4 jwm=jwmfr,jwmto
      list(jwm,1)=0
      list(jwm,2)=0
      if(rw(jwm).lt.rshrt2.and.
     + (iand(idt(iwm),mdynam).eq.ldynam.or.
     + iand(idt(jwm),mdynam).eq.ldynam)) list(jwm,1)=1
      if(rw(jwm).lt.rrest2 .and.
     + (iand(idt(iwm),mrestr).eq.lrestr.or.
     + iand(idt(jwm),mrestr).eq.lrestr)) list(jwm,1)=1
    4 continue
      if(npww.eq.2) then
      do 2 jwm=jwmfr,jwmto
      if(rw(jwm).lt.rlong2.and.
     + (iand(idt(iwm),mdynam).eq.ldynam.or.
     + iand(idt(jwm),mdynam).eq.ldynam)) list(jwm,2)=1
      if(list(jwm,1).eq.1) list(jwm,2)=0
    2 continue
      endif
c
c     short range pairlist
c
      nlist=0
      do 9 jwm=jwmfr,jwmto
      if(list(jwm,1).eq.1) then
      nlist=nlist+1
      list(nlist,1)=jwm
      endif
    9 continue
      if(npairs+nlist.gt.mpairs)
     + call md_abort('Insufficient memory for pairlist',0)
      do 5 ilist=1,nlist
      lwwj(npairs+ilist)=list(ilist,1)
    5 continue
      lwwjpt(iwm-iwmfr+1,1)=npairs+1
      npairs=npairs+nlist
      lsww=lsww+nlist
      lwwin(iwm-iwmfr+1,1)=nlist
c
c     long range pairs
c
      nlist=0
      if(npww.eq.2) then
      do 12 jwm=jwmfr,jwmto
      if(list(jwm,2).eq.1) then
      nlist=nlist+1
      list(nlist,2)=jwm
      endif
   12 continue
      if(npairs+nlist.gt.mpairs)
     + call md_abort('Insufficient memory for pairlist',0)
      do 13 ilist=1,nlist
      lwwj(npairs+ilist)=list(ilist,2)
   13 continue
      lwwjpt(iwm-iwmfr+1,2)=npairs+1
      npairs=npairs+nlist
      llww=llww+nlist
      lwwin(iwm-iwmfr+1,2)=nlist
      endif
    1 continue
c
      return
      end
      subroutine cf_lsw(lpbc,lpbcs,xs,isdt,isgr,isfr,isto,xwm,iwdt,iwz,
     + iwfr,iwto,ndim,mpairs,npairs,lswjpt,lswin,lswj,list,rwx,rw)
c
      implicit none
c
#include "cf_common.fh"
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
      if(lpbc.or.lpbcs) call cf_pbc(1,rwx,mscr,rwx,mscr,0,iwmfr,iwmto)
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
      if(lpbc.or.lpbcs) call cf_pbc(1,rwx,mscr,rwx,mscr,0,iwmfr,iwmto)
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
      if(iqmmm.ne.0) then
      if(iand(isdt(isa),lquant).eq.lquant) then
      do 9015 iwm=iwmfr,iwmto
      if(rw(iwm).lt.rqmmm2) iwz(iwm)=1
 9015 continue
      endif
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
      subroutine cf_lss(lself,lpbc,lpbcs,xs,iga,idt,igr,isgm,iss,ish,
     + isz,
     + isfr,isto,jsfr,jsto,nlocs,mpairs,npairs,lssjpt,lssin,lssj,list,
     + rwx,rw,ndxt,nt,ndxx,nx,lda,rda,lseq,lsthop,mprot)
c
      implicit none
c
#include "cf_common.fh"
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
c      write(*,'(i5,a,i5)') me,' cf_lss nseq=',nseq
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
      if(lpbc.or.lpbcs) call cf_pbc(1,rwx,mscr,rwx,mscr,0,jsafr,jsato)
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
      if(lpbc.or.lpbcs) call cf_pbc(1,rwx,mscr,rwx,mscr,0,jsafr,jsato)
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
      if(iqmmm.ne.0) then
      if(iand(idt(isa),lquant).eq.0) then
      do 24 jsa=jsafr,jsato
      if(iand(idt(jsa),lquant).ne.0.and.rw(jsa).lt.rqmmm2)
     + isz(isa)=1
   24 continue
      else
      do 25 jsa=jsafr,jsato
      if(iand(idt(jsa),lquant).eq.0.and.rw(jsa).lt.rqmmm2)
     + isz(jsa)=1
   25 continue
      endif
      endif
cx      if(nhop.gt.0) then
c      if(ish(isa).gt.0.and.iand(ish(isa),1).ne.0) then
c      do 326 jsa=jsafr,jsato
c      list(jsa,1)=0
c  326 continue
c      endif
c      do 327 jsa=jsafr,jsato
c      if(ish(jsa).gt.0.and.iand(ish(jsa),1).ne.0)
c     + list(jsa,1)=0
c  327 continue
cx      endif
c
c     donor acceptor pair list
c
      if(nhop.gt.0) then
      idon=0
      iacc=0
      idown=0
      iaown=0
      if(ish(isa).gt.0) then
      if(iand(ish(isa),1).eq.0) then
      idon=isa
      idown=me
      else
      iacc=isa
      iaown=me
      endif
      endif
      if(idon.gt.0.or.iacc.gt.0) then
      do 225 jsa=jsafr,jsato
      if(isgm(isa).ne.isgm(jsa)) then
      if(idon.eq.isa) iacc=0
      if(iacc.eq.isa) idon=0
      if(ish(jsa).gt.0) then
      if(iacc.gt.0.and.iand(ish(jsa),1).eq.0) then
      idon=jsa
      if(lself) idown=me
      elseif(idon.gt.0.and.iand(ish(jsa),1).eq.1) then
      iacc=jsa
      if(lself) iaown=me
      endif
      endif
      if(idon.gt.0.and.iacc.gt.0.and.rw(jsa).lt.rhop2) then
c
      if(lseq(isgm(idon))+ish(idon)/2.le.mprot(idon)
     + .or.lseq(isgm(iacc))-ish(iacc)/2.ge.1) then
c
      do 226 i=1,nhops
      if(lsthop(1,i).eq.isgm(iacc).and.lsthop(2,i).eq.isgm(idon)) then
      goto 227
      endif
  226 continue
c
c
c
      call qhop_asite(msa,xs,idon,iga(idon),isgm(idon),lseq(isgm(idon)),
     + iacc,iga(iacc),isgm(iacc),lseq(isgm(iacc)),xacc,idoff,iaoff,
     + lsame)
c
c     do not allow hops to the same protonation state if the type of
c     residue/molecule is identical (e.q. water dissociating into OH and H3O)
c
      if(.not.lsame.or.lseq(isgm(idon)).ne.lseq(isgm(iacc))) then
c
      nlda=nlda+1
c      write(*,2222) me,nlda,idon,iga(idon),isgm(idon),lseq(isgm(idon)),
c     + ish(idon),iacc,iga(iacc),isgm(iacc),lseq(isgm(iacc)),ish(iacc),
c     + sqrt(rw(jsa))
c 2222 format(12i5,f12.6)
      if(nlda.gt.3*nhop) call md_abort('Increase nhop',nhop)
      lda(1,nlda)=iga(idon)
      lda(2,nlda)=idon
      lda(3,nlda)=isgm(idon)
      lda(4,nlda)=lseq(isgm(idon))+ish(idon)/2
      lda(5,nlda)=idoff
      lda(6,nlda)=iga(iacc)
      lda(7,nlda)=iacc
      lda(8,nlda)=isgm(iacc)
      lda(9,nlda)=lseq(isgm(iacc))-ish(iacc)/2
      lda(10,nlda)=iaoff
      lda(11,nlda)=idown
      lda(12,nlda)=0
      lda(13,nlda)=0
      lda(14,nlda)=iaown
      lda(15,nlda)=0
      lda(16,nlda)=0
      rda(1,nlda)=xacc(1)
      rda(2,nlda)=xacc(2)
      rda(3,nlda)=xacc(3)
      rda(4,nlda)=xacc(4)
      rda(5,nlda)=xacc(5)
      rda(6,nlda)=xacc(6)
      rda(7,nlda)=xacc(7)
      rda(8,nlda)=xacc(8)
      rda(9,nlda)=xacc(9)
      rda(10,nlda)=xacc(10)
      rda(11,nlda)=xacc(11)
c
      endif
c
  227 continue
c
      endif
      endif
      endif
c
  225 continue
      endif
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
      if(iqmmm.ne.0) then
      do 38 jsa=isafr+1,isato
      if(igr(jsa).eq.igr(jsa-1)) then
      if(iand(idt(jsa),lquant).ne.lquant.and.isz(jsa-1).ne.0)
     + isz(jsa)=1
      if(isz(jsa).ne.0.and.iand(idt(jsa-1),lquant).ne.lquant)
     + isz(jsa-1)=1
      endif
   38 continue
      do 39 jsa=isato-1,isafr,-1
      if(igr(jsa).eq.igr(jsa+1)) then
      if(iand(idt(jsa),lquant).ne.lquant.and.isz(jsa+1).ne.0)
     + isz(jsa)=1
      if(isz(jsa).ne.0.and.iand(idt(jsa+1),lquant).ne.lquant)
     + isz(jsa+1)=1
      endif
   39 continue
      do 40 jsa=jsafr+1,jsato
      if(igr(jsa).eq.igr(jsa-1)) then
      if(iand(idt(jsa),lquant).ne.lquant.and.isz(jsa-1).ne.0)
     + isz(jsa)=1
      if(isz(jsa).ne.0.and.iand(idt(jsa-1),lquant).ne.lquant)
     + isz(jsa-1)=1
      endif
   40 continue
      do 41 jsa=jsato-1,jsafr,-1
      if(igr(jsa).eq.igr(jsa+1)) then
      if(iand(idt(jsa),lquant).ne.lquant.and.isz(jsa+1).ne.0)
     + isz(jsa)=1
      if(isz(jsa).ne.0.and.iand(idt(jsa+1),lquant).ne.lquant)
     + isz(jsa+1)=1
      endif
   41 continue
      endif
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
      subroutine cf_lsb(lself,iga,idt,istemp,isfr,isto,jsfr,jsto,
     + idb,nsb,nbonds,mbonds,nconst,lbonds,itemp)
c
      implicit none
c
#include "cf_common.fh"
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
      do 8 isbloc=1,nbonds
      isb=lbonds(isbloc)
      ibloc=itemp(idb(isb,1))
      jbloc=itemp(idb(isb,2))
      iq=iand(idt(ibloc),lquant)
      jq=iand(idt(jbloc),lquant)
      if(iq.ne.lquant.and.jq.eq.lquant) istemp(ibloc)=1
      if(iq.eq.lquant.and.jq.ne.lquant) istemp(jbloc)=1
    8 continue
      endif
c
      return
      end
      subroutine cf_lsh(lself,iga,idt,istemp,isfr,isto,jsfr,jsto,
     + idh,nsh,nangls,mangls,nconst,langls,itemp)
c
      implicit none
c
#include "cf_common.fh"
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
      subroutine cf_lsd(lself,iga,idt,istemp,isfr,isto,jsfr,jsto,
     + idd,nsd,ndihes,mdihes,nconst,ldihes,itemp)
c
      implicit none
c
#include "cf_common.fh"
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
      do 8 isdloc=1,ndihes
      isd=ldihes(isdloc)
      idloc=itemp(idd(isd,1))
      jdloc=itemp(idd(isd,2))
      kdloc=itemp(idd(isd,3))
      ldloc=itemp(idd(isd,4))
      iq=iand(idt(idloc),lquant)
      jq=iand(idt(jdloc),lquant)
      kq=iand(idt(kdloc),lquant)
      lq=iand(idt(ldloc),lquant)
      if(iq.eq.lquant.or.jq.eq.lquant.or.
     + kq.eq.lquant.or.lq.eq.lquant) then
      if(iq.ne.lquant) istemp(idloc)=1
      if(jq.ne.lquant) istemp(jdloc)=1
      if(kq.ne.lquant) istemp(kdloc)=1
      if(lq.ne.lquant) istemp(ldloc)=1
      endif
    8 continue
      endif
c
      return
      end
      subroutine cf_lso(lself,iga,idt,istemp,isfr,isto,jsfr,jsto,
     + ido,nso,nimprs,mimprs,nconst,limprs,itemp)
c
      implicit none
c
#include "cf_common.fh"
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
      do 8 isoloc=1,nimprs
      iso=limprs(isoloc)
      ioloc=itemp(ido(iso,1))
      joloc=itemp(ido(iso,2))
      koloc=itemp(ido(iso,3))
      loloc=itemp(ido(iso,4))
      iq=iand(idt(ioloc),lquant)
      jq=iand(idt(joloc),lquant)
      kq=iand(idt(koloc),lquant)
      lq=iand(idt(loloc),lquant)
      if(iq.eq.lquant.or.jq.eq.lquant.or.
     + kq.eq.lquant.or.lq.eq.lquant) then
      if(iq.ne.lquant) istemp(ioloc)=1
      if(jq.ne.lquant) istemp(joloc)=1
      if(kq.ne.lquant) istemp(koloc)=1
      if(lq.ne.lquant) istemp(loloc)=1
      endif
    8 continue
      endif
c
      return
      end
      subroutine cf_lst(lself,iga,idt,istemp,isfr,isto,jsfr,jsto,
     + idtn,nst,nthird,mthird,lthird,itemp)
c
      implicit none
c
#include "cf_common.fh"
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
      if(iqmmm.ne.0) then
      do 6 isa=1,nsatot
      itemp(isa)=0
    6 continue
      do 7 jsa=jsafr,jsato
      itemp(iga(jsa))=jsa
    7 continue
      do 8 isa=isafr,isato
      itemp(iga(isa))=isa
    8 continue
      do 9 istloc=1,nthird
      ist=lthird(istloc)
      itloc=-1
      do 10 isa=1,nsatot
      itloc=itemp(isa)
      if(idtn(isa-1,1).lt.ist.and.idtn(isa,1).ge.ist) goto 11
   10 continue
   11 continue
      if(itloc.lt.0) call md_abort('Error in cf_lst',me)
      jtloc=itemp(idtn(ist,2))
      iq=iand(idt(itloc),lquant)
      jq=iand(idt(jtloc),lquant)
      if(iq.ne.lquant.and.jq.eq.lquant) istemp(itloc)=1
      if(iq.eq.lquant.and.jq.ne.lquant) istemp(jtloc)=1
    9 continue
      endif
c
      return
      end
