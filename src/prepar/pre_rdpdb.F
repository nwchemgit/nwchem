      logical function pre_rdpdb(filpdb,lfnpdb,filtmp,lfntmp,lnams,
     + csa,isar,isgm,isfnd,xs,vs,msa,nsa,nmis,msm,nsm,iwmr,iwfnd,
     + xw,vw,mwm,mwa,nwm,nwa)
c
c $Id: pre_rdpdb.F,v 1.5 1998-04-29 18:35:22 d3j191 Exp $
c
      implicit none
c
      integer lfnpdb,lfntmp,msa,nsa,nmis,msm,nsm,mwm,mwa,nwm,nwa
      logical lnams,lnames
      character*255 filpdb,filtmp
      character*16 csa(msa)
      integer isar(msa),isgm(msa),iwmr(mwm),isfnd(msa),iwfnd(mwm)
      real*8 xs(3,msa),vs(3,msa)
      real*8 xw(3,mwa,mwm),vw(3,mwa,mwm)
c
      character*6 name1,name2
      character*80 card
      integer i,last,nprv,num,iw
      real*8 x,y,z
      character*3 namr
      character*4 nama
      character*6 hdr
      logical lsolut
c
      lnames=lnams
      iw=1
c
c     open pdb file
c     -------------
c
      open(unit=lfnpdb,file=filpdb(1:index(filpdb,' ')-1),
     + form='formatted',status='old',err=9999)
c
c     open the alternate names file
c     -----------------------------
c
      if(lnames) then
      lnames=.false.
      open(unit=lfntmp,file=filtmp(1:index(filtmp,' ')-1),
     + form='formatted',status='old',err=99)
      lnames=.true.
   99 continue
      endif
c
c     initializations
c     ---------------
c
      nmis=0
      last=0
      nprv=0
      do 1 i=1,msa
      isfnd(i)=0
      isar(i)=0
      xs(1,i)=0.0d0
      xs(2,i)=0.0d0
      xs(3,i)=0.0d0
      vs(1,i)=0.0d0
      vs(2,i)=0.0d0
      vs(3,i)=0.0d0
      if(isgm(i).gt.last) last=isgm(i)
    1 continue
      nsa=0
      lsolut=.true.
c
  100 continue
c
      read(lfnpdb,1000,end=9,err=9999) card
 1000 format(a)
c
      if(card(1:4).ne.'ATOM'.and.card(1:6).ne.'HETATM') goto 100
c
      read(card,1001,end=9999,err=9999) hdr,nama,namr,num,x,y,z
 1001 format(a6,6x,a4,1x,a3,2x,i4,4x,3f8.3)
c
      if(num.gt.last.or.num.lt.nprv) then
      if(lnames) then
      close(unit=lfntmp)
      lnames=.false.
      endif
      lsolut=.false.
      endif
c
      if(lnames) then
      read(lfntmp,1002,end=9999,err=9999) name1,name2
 1002 format(2a6)
      endif
c
      if(lsolut) then
      do 3 i=1,msa
      if(num.eq.isgm(i)) then
      if(lnames) then
      if(csa(i)(11:14).ne.name2(1:4)) goto 3
      if(nama.ne.name1(1:4)) goto 3
      else
      if(nama.ne.csa(i)(11:14)) goto 3
      endif
      isfnd(i)=1
      xs(1,i)=1.0d-01*x
      xs(2,i)=1.0d-01*y
      xs(3,i)=1.0d-01*z
      goto 100
      endif
    3 continue
      else
      iwfnd(iw)=iwfnd(iw+1)
      iwmr(iw)=0
      endif
c
      goto 100
c
    9 continue
c
      close(unit=lfnpdb)
c
      nmis=0
      do 4 i=1,msa
      if(isfnd(i).eq.0) nmis=nmis+1
    4 continue
      nsa=msa
c
      pre_rdpdb=.true.
      return
c
 9999 continue
      pre_rdpdb=.false.
      return
      end
