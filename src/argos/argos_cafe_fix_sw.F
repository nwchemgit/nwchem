      subroutine argos_cafe_fix_sw(lfnout,lfn,numfix,iwgmn,iwdt,numw,
     + isgan,isatt,isdt,issgm,nums,snam,xwm,xs,idb,msb,nsb,itemp)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "msgids.fh"
#include "bitops.fh"
c
      integer argos_cafe_atonum,argos_cafe_numcon
      external argos_cafe_atonum,argos_cafe_numcon
c
      integer lfnout,lfn,numfix,numw,nums,msb,nsb
      integer iwdt(mwm),iwgmn(mwm),isgan(msa),isatt(msa),isdt(msa)
      integer issgm(msa),idb(msb,4),itemp(nsatot)
      character*16 snam(nsatot)
      real*8 xwm(mwm,3),xs(msa,3)
c
      integer i,j,ndof(3),ifix,ifr,ito,l,len,number,itarget,ntotal(2)
      integer ib,jb,nsad
      real*8 fact(3),xtarget(3),rfix,dist2,rfix2
      character*80 string,target,residu
      character*16 tag,token
      logical lnumber
c
      integer icolon,icomma,idash,length,jfr,jto
c
      if(me.eq.0) then
      write(lfnout,1000) numfix
 1000 format(/,' POSITIONAL CONSTRAINTS',/,
     + ' Number of fix directives is ',i5,/)
      endif
c
      number=0
c
      do 1 l=1,numfix
c
      if(me.eq.0) then
      read(lfn,2000) string
 2000 format(a80)
c
      target=string(41:80)
      len=index(target,':')
      residu=' '
      if(len.gt.0) then
      residu=target(1:len-1)
      target=target(len+1:80)
      target(13-len:11)=residu
      endif
      do 33 i=1,4
      if(target(i:i).eq.'_') target(i:i)=' '
   33 continue
c
      itarget=0
      do 34 i=1,nsatot
      if(snam(i)(6:16).eq.target(1:11)) itarget=i
   34 continue
c
      endif
c
      call ga_brdcst(mrg_d47,itarget,ma_sizeof(mt_int,1,mt_byte),0)
      call util_char_ga_brdcst(mrg_d47,string,0)
c
      if(itarget.gt.0) then
      do 35 i=1,nums
      if(isgan(i).eq.itarget) then
      itarget=i
      xtarget(1)=xs(i,1)
      xtarget(2)=xs(i,2)
      xtarget(3)=xs(i,3)
      goto 36
      endif
   35 continue
      xtarget(1)=zero
      xtarget(2)=zero
      xtarget(3)=zero
      itarget=0
   36 continue
      call ga_dgop(mrg_d47,xtarget,3,'+')
      endif
c
      if(string(1:3).eq.'fix') then
c
      read(string(6:40),'(i3,2i10,f12.6)') ifix,ifr,ito,rfix
      target=string(41:80)
      len=index(target,':')
      residu=' '
      if(len.gt.0) then
      residu=target(1:len-1)
      target=target(len+1:80)
      endif
      len=4
      do 21 i=1,4
      if(target(i:i).eq.'_') target(i:i)=' '
   21 continue
c
c     solvent
c
      if(string(5:5).eq.'w') then
      if(numw.gt.0) then
      if(ifr.eq.0.or.ito.eq.0) then
      if(ifix.eq.0) then
      do 2 i=1,numw
      iwdt(i)=iand(iwdt(i),not(mfixed))
    2 continue
      endif
      if(ifix.eq.1) then
      do 3 i=1,numw
      iwdt(i)=ior(iwdt(i),lfixed)
    3 continue
      endif
      else
      if(ifr.gt.0) then
      do 4 i=1,numw
      if(iwgmn(i).ge.ifr.and.iwgmn(i).le.ito) then
      if(ifix.eq.0) iwdt(i)=iand(iwdt(i),not(mfixed))
      if(ifix.eq.1) iwdt(i)=ior(iwdt(i),lfixed)
      endif
    4 continue
      else
      rfix2=rfix*rfix
      do 41 i=1,numw
      dist2=(xwm(i,1)-xtarget(1))**2+(xwm(i,2)-xtarget(2))**2+
     + (xwm(i,3)-xtarget(3))**2
      if(ito.eq.-1.and.dist2.le.rfix2) then
      if(ifix.eq.0) iwdt(i)=iand(iwdt(i),not(mfixed))
      if(ifix.eq.1) iwdt(i)=ior(iwdt(i),lfixed)
      elseif(ito.eq.-2.and.dist2.ge.rfix2) then
      if(ifix.eq.0) iwdt(i)=iand(iwdt(i),not(mfixed))
      if(ifix.eq.1) iwdt(i)=ior(iwdt(i),lfixed)
      endif
   41 continue
      endif
      endif
      endif
      endif
c
c     solute
c
      if(string(5:5).eq.'s') then
      if(nums.gt.0) then
      if(ifr.eq.0.or.ito.eq.0) then
      if(ifix.eq.0) then
      do 5 i=1,nums
      tag=snam(isgan(i))
      do 51 j=1,16
      if(target(j:j).eq.'?') tag(j+5:j+5)='?'
   51 continue
      if((target(1:2).eq.'  '.or.target(1:len).eq.tag(6:len+5)).and.
     + (residu(1:1).eq.' '.or.residu(1:3).eq.tag(1:3))) then
      isdt(i)=iand(isdt(i),not(mfixed))
      number=number+1
      endif
    5 continue
      endif
      if(ifix.eq.1) then
      do 6 i=1,nums
      tag=snam(isgan(i))
      do 61 j=1,16
      if(target(j:j).eq.'?') tag(j+5:j+5)='?'
   61 continue
      if((target(1:2).eq.'  '.or.target(1:len).eq.tag(6:len+5)).and.
     + (residu(1:1).eq.' '.or.residu(1:3).eq.tag(1:3))) then
      isdt(i)=ior(isdt(i),lfixed)
      number=number+1
      endif
    6 continue
      endif
      else
      if(ifr.gt.0) then
      do 7 i=1,nums
      if(issgm(i).ge.ifr.and.issgm(i).le.ito) then
      tag=snam(isgan(i))
      do 71 j=1,16
      if(target(j:j).eq.'?') tag(j+5:j+5)='?'
   71 continue
      if((target(1:2).eq.'  '.or.target(1:len).eq.tag(6:len+5)).and.
     + (residu(1:1).eq.' '.or.residu(1:3).eq.tag(1:3))) then
      number=number+1
      if(ifix.eq.0) isdt(i)=ior(isdt(i),not(mfixed))
      if(ifix.eq.1) isdt(i)=ior(isdt(i),lfixed)
      endif
      endif
    7 continue
      else
      rfix2=rfix*rfix
      do 72 i=1,nums
      dist2=(xs(i,1)-xtarget(1))**2+(xs(i,2)-xtarget(2))**2+
     + (xs(i,3)-xtarget(3))**2
      if(ito.eq.-1.and.dist2.le.rfix2) then
      if(ifix.eq.0) isdt(i)=iand(isdt(i),not(mfixed))
      if(ifix.eq.1) isdt(i)=ior(isdt(i),lfixed)
      elseif(ito.eq.-2.and.dist2.ge.rfix2) then
      if(ifix.eq.0) isdt(i)=iand(isdt(i),not(mfixed))
      if(ifix.eq.1) isdt(i)=ior(isdt(i),lfixed)
      endif
   72 continue
      endif
      endif
      endif
      endif
c
c     heavy solute atoms only
c
      if(string(5:5).eq.'x') then
      if(nums.gt.0) then
      if(ifr.eq.0.or.ito.eq.0) then
      if(ifix.eq.0) then
      do 8 i=1,nums
      tag=snam(isgan(i))
      do 81 j=1,16
      if(target(j:j).eq.'?') tag(j+5:j+5)='?'
   81 continue
      if((target(1:2).eq.'  '.or.target(1:len).eq.tag(6:len+5)).and.
     + (residu(1:1).eq.' '.or.residu(1:3).eq.tag(1:3))) then
      if(argos_cafe_atonum(isatt(i),int_mb(i_num)).ne.1)
     + isdt(i)=iand(isdt(i),not(mfixed))
      number=number+1
      endif
    8 continue
      endif
      if(ifix.eq.1) then
      do 9 i=1,nums
      tag=snam(isgan(i))
      do 91 j=1,16
      if(target(j:j).eq.'?') tag(j+5:j+5)='?'
   91 continue
      if((target(1:2).eq.'  '.or.target(1:len).eq.tag(6:len+5)).and.
     + (residu(1:1).eq.' '.or.residu(1:3).eq.tag(1:3))) then
      if(argos_cafe_atonum(isatt(i),int_mb(i_num)).ne.1)
     + isdt(i)=ior(isdt(i),lfixed)
      number=number+1
      endif
    9 continue
      endif
      else
      do 10 i=1,nums
      if(issgm(i).ge.ifr.and.issgm(i).le.ito) then
      tag=snam(isgan(i))
      do 101 j=1,16
      if(target(j:j).eq.'?') tag(j+5:j+5)='?'
  101 continue
      if((target(1:2).eq.'  '.or.target(1:len).eq.tag(6:len+5)).and.
     + (residu(1:1).eq.' '.or.residu(1:3).eq.tag(1:3))) then
      if(argos_cafe_atonum(isatt(i),int_mb(i_num)).ne.1) then
      if(ifix.eq.0) isdt(i)=ior(isdt(i),not(mfixed))
      if(ifix.eq.1) isdt(i)=ior(isdt(i),lfixed)
      endif
      number=number+1
      endif
      endif
   10 continue
      endif
      endif
      endif
c
c     selection
c
      if(string(5:5).eq.'S') then
      target=string(41:80)
c
      do 55 i=1,nsatot
      itemp(i)=0
   55 continue
c
c     for the residues
c
  140 continue
      icolon=index(target,':')
      if(icolon.gt.0) then
      token='                '
      icomma=index(target(1:icolon-1),',')
      if(icomma.gt.0) then
      length=icomma-1
      else
      length=icolon-1
      endif
      token(1:length)=target(1:length)
      do 57 i=1,length
      lnumber=.false.
      if(token(i:i).eq.'1') lnumber=.true.
      if(token(i:i).eq.'2') lnumber=.true.
      if(token(i:i).eq.'3') lnumber=.true.
      if(token(i:i).eq.'4') lnumber=.true.
      if(token(i:i).eq.'5') lnumber=.true.
      if(token(i:i).eq.'6') lnumber=.true.
      if(token(i:i).eq.'7') lnumber=.true.
      if(token(i:i).eq.'8') lnumber=.true.
      if(token(i:i).eq.'9') lnumber=.true.
      if(token(i:i).eq.'0') lnumber=.true.
      if(token(i:i).eq.'-') lnumber=.true.
      if(.not.lnumber) goto 58
   57 continue
c     number
      idash=index(token(1:length),'-')
      if(idash.gt.0) then
c     range numbered residues
      tag='                '
      tag(18-idash:16)=token(1:idash-1)
      read(tag,'(i16)') jfr
      tag='                '
      tag(17+idash-length:16)=token(idash+1:length)
      read(tag,'(i16)') jto
      else
c     single numbered residue
      tag='                '
      tag(17-length:16)=token(1:length)
      read(tag,'(i16)') jfr
      jto=jfr
      endif
      do 60 i=1,nsatot
      read(snam(i)(11:16),'(i6)') j
      if(j.ge.jfr.and.j.le.jto) itemp(i)=1
   60 continue
      goto 59
   58 continue
c     single named residue
      do 141 i=1,nsatot
      if(snam(i)(1:length).eq.token(1:length)) itemp(i)=1
  141 continue
   59 continue
      if(icomma.gt.0) then
      target=target(icomma+1:80)
      goto 140
      else
      target=target(icolon+1:80)
      endif
      else
      do 56 i=1,nsatot
      itemp(i)=1
   56 continue
      endif
c
c     for the atoms
c
      if(index(target,' ')-1.eq.0) then
      do 145 i=1,nsatot
      itemp(i)=itemp(i)+1
  145 continue
      goto 144
      endif
c
  142 continue
      icomma=index(target,',')
      if(icomma.gt.0) then
      length=icomma-1
      else
      length=index(target,' ')-1
      endif
      do 147 i=1,length
      if(target(i:i).eq.'_') target(i:i)=' '
  147 continue
      do 143 i=1,nsatot
      token(1:length)=snam(i)(6:5+length)
      do 146 j=1,length
      if(target(j:j).eq.'?') token(j:j)='?'
  146 continue
      if(token(1:length).eq.target(1:length)) itemp(i)=itemp(i)+1
  143 continue
      if(icomma.gt.0) then
      target=target(icomma+1:80)
      goto 142
      endif
c
  144 continue
      do 66 i=1,nums
      if(itemp(isgan(i)).gt.1) then
      if(ifix.eq.0) isdt(i)=iand(isdt(i),not(mfixed))
      if(ifix.eq.1) isdt(i)=ior(isdt(i),lfixed)
      endif
   66 continue
c
      endif
c
c     solute atoms constraints
c
      do 76 i=1,nsatot
      itemp(i)=0
   76 continue
c
      do 62 i=1,nums
      itemp(isgan(i))=i
   62 continue
c
      do 64 j=1,4
      do 63 i=1,nsb
      ib=itemp(idb(i,1))
      jb=itemp(idb(i,2))
      if(ib.gt.0.and.jb.gt.0.and.idb(i,3).gt.0) then
      if(iand(isdt(ib),mfixed).eq.lfixed) isdt(jb)=ior(isdt(jb),lfixed)
      if(iand(isdt(jb),mfixed).eq.lfixed) isdt(ib)=ior(isdt(ib),lfixed)
      endif
   63 continue
   64 continue
c
c     make permanent
c
      if(string(5:5).eq.'p') then
      do 77 i=1,nums
      j=iand(isdt(i),mfixed)/4
      isdt(i)=lshift(rshift(isdt(i),2),2)
      isdt(i)=ior(isdt(i),j)
   77 continue
      do 78 i=1,numw
      j=iand(iwdt(i),mfixed)/4
      iwdt(i)=lshift(rshift(iwdt(i),2),2)
      iwdt(i)=ior(iwdt(i),j)
   78 continue
      endif
c
      endif
c
    1 continue
c
      ntotal(1)=0
      ntotal(2)=0
      do 201 i=1,nums
      if(iand(isdt(i),mfixed).eq.lfixed) ntotal(1)=ntotal(1)+1
  201 continue
      do 202 i=1,numw
      if(iand(iwdt(i),mfixed).eq.lfixed) ntotal(2)=ntotal(2)+1
  202 continue
c
      if(np.gt.1) then
      call ga_igop(mrg_d48,ntotal,2,'+')
      endif
c
      if(me.eq.0) then
      write(lfnout,'(a,i8,a,i8,a)')
     +  ' Solute fix applied to  ',ntotal(1),' out of ',nsa,' atoms'
      write(lfnout,'(a,i8,a,i8,a,/)')
     +  ' Solvent fix applied to ',ntotal(2),' out of ',nwm,' molecules'
      endif
c
c      if(numw.gt.0) then
c      if(iand(ifix,1).eq.1) then
c      do 1 i=1,numw
c      iwdt(i)=ior(iwdt(i),lfixed)
c    1 continue
c      endif
c      if(iand(iunfix,1).eq.1) then
c      do 2 i=1,numw
c      iwdt(i)=iand(iwdt(i),not(lfixed))
c    2 continue
c      endif
c      if(iand(ifix,2).eq.2) then
c      do 3 i=1,numw
c      iwdt(i)=ior(iwdt(i),lfixed)
c    3 continue
c      endif
c      if(iand(iunfix,2).eq.2) then
c      do 4 i=1,numw
c      iwdt(i)=iand(iwdt(i),not(lfixed))
c    4 continue
c      endif
c      endif
c
c      if(nums.gt.0) then
c      if(iand(ifix,1).eq.1) then
c      do 5 i=1,nums
c      isdt(i)=ior(isdt(i),lfixed)
c    5 continue
c      endif
c      if(iand(iunfix,1).eq.1) then
c      do 6 i=1,nums
c      isdt(i)=iand(isdt(i),not(lfixed))
c    6 continue
c      endif
c      if(iand(ifix,4).eq.4) then
c      do 7 i=1,nums
c      isdt(i)=ior(isdt(i),lfixed)
c    7 continue
c      endif
c      if(iand(iunfix,4).eq.4) then
c      do 8 i=1,nums
c      isdt(i)=iand(isdt(i),not(lfixed))
c    8 continue
c      endif
c      if(iand(ifix,8).eq.8) then
c      do 9 i=1,nums
c      if(argos_cafe_atonum(isatt(i),int_mb(i_num)).ne.1)
c     + isdt(i)=ior(isdt(i),lfixed)
c    9 continue
c      endif
c      if(iand(iunfix,8).eq.8) then
c      do 10 i=1,nums
c      if(argos_cafe_atonum(isatt(i),int_mb(i_num)).ne.1)
c     + isdt(i)=iand(isdt(i),not(lfixed))
c   10 continue
c      endif
c      endif
c
      ndof(1)=0
      ndof(2)=0
      ndof(3)=0
c
      do 11 i=1,numw
      if(iand(iwdt(i),lfixed).ne.lfixed) ndof(1)=ndof(1)+1
   11 continue
      do 12 i=1,nums
      if(iand(isdt(i),lfixed).ne.lfixed) then
      ndof(2)=ndof(2)+1
      ndof(3)=ndof(3)+argos_cafe_numcon(int_mb(i_ibnd(2)),
     + mbt(2),isgan(i))
      endif
   12 continue
c
      call ga_igop(mrg_d46,ndof,3,'+')
      nsc=ndof(3)/2
c
      factmw=zero
      factms=zero
      factmp=zero
      if(noshak.eq.2.or.noshak.eq.3) then
      nsad=3*ndof(2)-3*ndums
      else
      nsad=3*ndof(2)-2*ndums
      endif
      if(ndof(1)*(3*nwa-nwc)-3*islow.gt.0)
     + factmw=two/(rgas*dble(ndof(1)*(3*nwa-nwc)-3*islow))
      if(3*ndof(2)-nsc-3*islow.gt.0)
     + factms=two/(rgas*dble(nsad-nsc-3*islow))
      if(ndof(1)*(3*nwa-nwc)+nsad-nsc-3*islow.gt.0)
     + factmp=two/(rgas*dble(ndof(1)*(3*nwa-nwc)+nsad-nsc-3*islow))
c
      fact(1)=factmw
      fact(2)=factms
      fact(3)=factmp
      call ga_brdcst(mrg_d47,fact,3*ma_sizeof(mt_int,1,mt_byte),0)
      factmw=fact(1)
      factms=fact(2)
      factmp=fact(3)
c
      return
      end
c $Id$
