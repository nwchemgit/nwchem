      logical function pre_rdpdb(lfnout,filpdb,lfnpdb,filtmp,lfntmp,
     + lnams,csa,isar,isgm,isfnd,xs,vs,msa,nsa,nmis,msm,nsm,iwmr,iwfnd,
     + xw,vw,cwa,mwm,mwa,nwm,nwa,namslv,iwater,altloc,chain,model,box)
c
c $Id: pre_rdpdb.F,v 1.17 2000-02-20 16:05:14 d3j191 Exp $
c
      implicit none
c
#include "util.fh"
c
      integer lfnout
      integer lfnpdb,lfntmp,msa,nsa,nmis,msm,nsm,mwm,mwa,nwm,nwa
      logical lnams,lnames
      character*255 filpdb,filtmp
      character*16 cwa(mwa),csa(msa)
      integer isar(msa),isgm(msa),iwmr(mwm),isfnd(msa),iwfnd(mwm)
      real*8 xs(3,msa),vs(3,msa),box(3)
      real*8 xw(3,mwa,mwm),vw(3,mwa,mwm)
      character*3 namslv
      integer iwater,model
      character*1 altloc,chain
c
      character*6 name1,name2
      character*80 card
      integer i,last,nprv,num,iw,length,mdl,nump
      real*8 x,y,z
      character*3 namr
      character*4 nama
      character*6 hdr
      logical lsolut
      character*1 alocal,achain
c
      lnames=lnams
      iw=0
      iwater=0
      mdl=0
      nump=0
c
      length=index(filpdb,' ')-1
      if(util_print('where',print_debug)) then
      write(lfnout,2006) filpdb(1:length)
 2006 format(/,'pdb_read  ',a)
      endif
c
c     open pdb file
c     -------------
c
      open(unit=lfnpdb,file=filpdb(1:index(filpdb,' ')-1),
     + form='formatted',status='old',err=9999)
c
      if(util_print('where',print_debug)) then
      write(lfnout,2000) filpdb(1:length)
 2000 format(' PDB OPEN  ',a)
      endif
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
      nwm=0
      lsolut=.true.
c
  100 continue
c
      read(lfnpdb,1000,end=9,err=9999) card
 1000 format(a)
c
c     BOX
c
      if(card(1:6).eq.'BOX   ') then
      read(card(7:42),'(3f12.6)') box
      box(1)=0.1d0*box(1)
      box(2)=0.1d0*box(2)
      box(3)=0.1d0*box(3)
      goto 100
      endif
c
c     MODEL
c
      if(model.ne.mdl) then
      if(card(1:6).eq.'MODEL ') then
      read(card(7:14),1001) mdl
 1001 format(i8)
      if(model.eq.mdl.and.util_print('where',print_debug)) then
      write(lfnout,2005) model
 2005 format(10x,'MODEL ',i8)
      endif
      goto 100
      endif
      endif
c
c     END and ENDMDL
c
      if(card(1:6).eq.'ENDMDL') then
      if(model.eq.mdl) goto 9
      goto 100
      endif
      if(card(1:6).eq.'END   ') goto 9
c
      if(card(1:4).ne.'ATOM'.and.card(1:6).ne.'HETATM') goto 100
c
      if(mdl.ne.model) goto 100
c
      read(card,1002,end=9999,err=9999) hdr,nama,alocal,namr,achain,num,
     + x,y,z
 1002 format(a6,6x,a4,a1,a3,1x,a1,i4,4x,3f8.3)
      if(util_print('where',print_debug)) then
      write(lfnout,1002) hdr,nama,alocal,namr,achain,num,x,y,z
      endif
      if(alocal.ne.' '.and.alocal.ne.altloc.and.altloc.ne.'*') goto 100
      if(achain.ne.' '.and.achain.ne.chain.and.chain.ne.'*') goto 100
c
      if(chain.eq.'*') then
      if(num.lt.nump) then
      if(achain.eq.'B') num=num+1000
      if(achain.eq.'C') num=num+2000
      if(achain.eq.'D') num=num+3000
      if(achain.eq.'E') num=num+4000
      if(achain.eq.'F') num=num+5000
      if(achain.eq.'G') num=num+6000
      if(achain.eq.'H') num=num+7000
      if(achain.eq.'I') num=num+8000
      if(achain.eq.' ') num=num+9000
      endif
      endif
      if(card(18:20).eq.namslv) lsolut=.false.
      if(num.ne.nump) nump=num
c
      if(.not.lsolut.or.num.gt.last.or.num.lt.nprv) then
      if(lnames) then
      close(unit=lfntmp)
      lnames=.false.
      endif
      lsolut=.false.
      endif
c
      if(lnames) then
      read(lfntmp,1003,end=9999,err=9999) name1,name2
 1003 format(2a6)
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
      if(iw.ne.num) then
      nwm=nwm+1
      iw=num
      iwfnd(nwm)=0
      endif
      iwfnd(nwm)=iwfnd(nwm)+1
      iwmr(nwm)=0
      if(namslv.eq.'HOH') iwater=1
      do 5 i=1,nwa
      if(cwa(i)(11:14).eq.nama.or.
     + (namslv.eq.'HOH'.and.
     +  ((nama.eq.' O  '.and.cwa(i)(11:14).eq.' OW ').or.
     +   (nama.eq.' H  '.and.cwa(i)(11:14).eq.' HW ')))) then
      xw(1,i,nwm)=1.0d-01*x
      xw(2,i,nwm)=1.0d-01*y
      xw(3,i,nwm)=1.0d-01*z
      endif
    5 continue
      endif
c
      goto 100
c
    9 continue
c
      close(unit=lfnpdb)
c
      if(util_print('where',print_debug)) then
      write(lfnout,2003) filpdb(1:length)
 2003 format(' PDB CLOSE ',a)
      endif
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

