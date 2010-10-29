      integer function argos_prep_shorts(xatm,latm,matm,ilo,ihi,jlo,jhi,
     + ishort,jshort,rshort)
c
c $Id$
c
c     function to find atom ishort in [ilo:ihi] and jshort in [jlo,jhi]
c              with the shortest distance rshort in nm
c
c     in  : xatm(3,matm) = atomic coordinates
c         : latm(2,matm) = atomic numbers
c         : matm         = dimension of coordinate list
c         : ilo          = atom index ilo
c         : ihi          = atom index ihi
c         : jlo          = atom index jlo
c         : jhi          = atom index jhi
c
c     out : ishort       = atom index i with shortest distance
c         : jshort       = atom index j with shortest distance
c         : rshort       = distance between atoms ishort and jshort
c
      implicit none
c
      real*8 argos_prep_atsiz
      external argos_prep_atsiz
c
      integer matm
      real*8 xatm(3,matm)
      integer latm(5,matm)
      integer ilo,ihi,jlo,jhi,ishort(10),jshort(10)
      real*8 rshort(10)
c
      integer i,j,k,l,number,ishrt,jshrt
      real*8 dist,short
      integer nsh,ish(100),jsh(100)
      real*8 rsh(100)
      logical new
c
      nsh=0
      do 1 i=ilo,ihi
      if(latm(2,i).gt.1) then
      do 2 j=jlo,jhi
      if(latm(2,j).gt.1) then
      dist=sqrt((xatm(1,i)-xatm(1,j))*(xatm(1,i)-xatm(1,j))+
     + (xatm(2,i)-xatm(2,j))*(xatm(2,i)-xatm(2,j))+
     + (xatm(3,i)-xatm(3,j))*(xatm(3,i)-xatm(3,j)))
      if(dist.lt.argos_prep_atsiz(latm(2,i))+
     + argos_prep_atsiz(latm(2,j))) then
      nsh=nsh+1
      if(nsh.gt.100)
     + call md_abort('argos_prep_shorts: increase dimension',1)
      ish(nsh)=i
      jsh(nsh)=j
      rsh(nsh)=dist
      endif
      endif
    2 continue
      endif
    1 continue
c
      number=0
c
      if(nsh.gt.0) then
      do 3 i=1,nsh-1
      do 4 j=i+1,nsh
      if(rsh(i).gt.rsh(j)) then
      dist=rsh(i)
      k=ish(i)
      l=jsh(i)
      rsh(i)=rsh(j)
      ish(i)=ish(j)
      jsh(i)=jsh(j)
      rsh(j)=dist
      ish(j)=k
      jsh(j)=l
      endif
    4 continue
    3 continue
      do 5 i=1,nsh
      do 6 j=1,number
      if(ish(i).eq.ishort(j).or.jsh(i).eq.jshort(j)) goto 5
    6 continue
      number=number+1
      if(number.gt.10)
     +  call md_abort('argos_prep_shorts:increase dimension',2)
      ishort(number)=ish(i)
      jshort(number)=jsh(i)
      rshort(number)=rsh(i)
    5 continue
      endif
c
cc
c
c
c      number=0
c
c      do 1 l=1,10
c      ishrt=0
c      jshrt=0
c      short=0.0d0
c      do 2 i=ilo,ihi
c      do 3 j=jlo,jhi
c      do 4 k=1,number
c      if(ishort(k).eq.i.and.jshort(k).eq.j) goto 3
c    4 continue
c      if(i.ne.j.and.latm(2,i).gt.1.and.latm(2,j).gt.1) then
c      dist=sqrt((xatm(1,i)-xatm(1,j))*(xatm(1,i)-xatm(1,j))+
c     + (xatm(2,i)-xatm(2,j))*(xatm(2,i)-xatm(2,j))+
c     + (xatm(3,i)-xatm(3,j))*(xatm(3,i)-xatm(3,j)))
c      if(dist.lt.argos_prep_atsiz(latm(2,i))+argos_prep_atsiz(latm(2,j))
c     + .and.(ishrt.eq.0.or.dist.lt.short)) then
c      ishrt=i
c      jshrt=j
c      short=dist
c      endif
c      endif
c    3 continue
c    2 continue
c      if(ishrt.gt.0) then
c      number=number+1
c      if(number.gt.10)
c     + call md_abort('Too many links between fragment pair',number)
c      ishort(number)=ishrt
c      jshort(number)=jshrt
c      rshort(number)=short
c      endif
c    1 continue
c
      argos_prep_shorts=number
c
      return
      end
