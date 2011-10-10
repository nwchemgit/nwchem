      subroutine argos_prep_atom_selection(atomname,isegment,
     + natoms,target,list)
c
      implicit none
c
      integer natoms
      character*16 atomname(natoms)
      integer isegment(natoms)
      character*255 target
      character*10 string,aname
      integer list(natoms)
c
      integer i,j,ip,ndxcol,length,isgmfr,isgmto
      logical lsgm,latm
c
c
      do 1 i=1,natoms
      list(i)=0
 1    continue
c
      length=index(target,' ')-1
      ndxcol=index(target,':')
      lsgm=ndxcol.gt.0
      latm=ndxcol.gt.0
      if(ndxcol.eq.0) then
      latm=.false.
      do 2 i=1,length
      lsgm=.false.
      if(target(i:i).eq.'1') lsgm=.true.
      if(target(i:i).eq.'2') lsgm=.true.
      if(target(i:i).eq.'3') lsgm=.true.
      if(target(i:i).eq.'4') lsgm=.true.
      if(target(i:i).eq.'5') lsgm=.true.
      if(target(i:i).eq.'6') lsgm=.true.
      if(target(i:i).eq.'7') lsgm=.true.
      if(target(i:i).eq.'8') lsgm=.true.
      if(target(i:i).eq.'9') lsgm=.true.
      if(target(i:i).eq.'0') lsgm=.true.
      if(target(i:i).eq.',') lsgm=.true.
      if(target(i:i).eq.'-') lsgm=.true.
      if(.not.lsgm) latm=.true.
      if(latm) goto 3
 2    continue
 3    continue
      if(lsgm) then
      ndxcol=length+1
      else
      ndxcol=0
      endif
      endif
c
      if(lsgm) then
      string='          '
      isgmto=0
      isgmfr=0
      ip=1
      do 4 i=1,ndxcol
      if(target(i:i).eq.','.or.i.eq.ndxcol) then
      string(10+ip-i:10)=target(ip:i-1)
      read(string,'(i10)') isgmto
      string='          '
      ip=i+1
      if(isgmfr.eq.0) isgmfr=isgmto
      elseif(target(i:i).eq.'-') then
      string(10+ip-i:10)=target(ip:i-1)
      read(string,'(i10)') isgmfr
      string='          '
      ip=i+1
      endif
      if(isgmto.gt.0) then
      do 5 j=1,natoms
      if(isegment(j).ge.isgmfr.and.isegment(j).le.isgmto) list(j)=-1
 5    continue
      isgmto=0
      isgmfr=0
      endif
 4    continue
      else
      if(latm) then
      do 6 j=1,natoms
      list(j)=-1
 6    continue
      else
      do 7 j=1,natoms
      list(j)=1
 7    continue
      endif
      endif
c
      if(latm) then
      string='          '
      ip=1
      do 8 i=ndxcol+1,length+1
      if(target(i:i).eq.','.or.i.eq.length+1) then
      do 9 j=1,natoms
      if(list(j).eq.-1) then
      aname(1:6)=atomname(j)(11:16)
      if(string(1:1).eq.'_') aname(1:1)='_'
      if(string(1:1).eq.'?') aname(1:1)='?'
      if(string(2:2).eq.'?') aname(2:2)='?'
      if(string(3:3).eq.'?') aname(3:3)='?'
      if(string(4:4).eq.'?') aname(4:4)='?'
      if(string(5:5).eq.'?') aname(5:5)='?'
      if(string(6:6).eq.'?') aname(6:6)='?'
      if(string(1:6).eq.aname(1:6)) list(j)=1
      endif
 9    continue
      string='          '
      ip=1
      else
      string(ip:ip)=target(i:i)
      ip=ip+1
      endif
 8    continue
      else
      do 10 i=1,natoms
      if(list(i).lt.0) list(i)=1
   10 continue
      endif
c
      return
      end
c $Id$
