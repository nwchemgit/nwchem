      logical function pre_short(xatm,latm,matm,ilo,ihi,jlo,jhi,
     + ishort,jshort,rshort)
c
c $Id: pre_short.F,v 1.8 2001-06-11 22:01:06 d3j191 Exp $
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
      real*8 pre_atsiz
      external pre_atsiz
c
      integer matm
      real*8 xatm(3,matm)
      integer latm(5,matm)
      integer ilo,ihi,jlo,jhi,ishort,jshort
      real*8 rshort
c
      integer i,j
      real*8 dist,short
c
      ishort=0
      short=0.0d0
c
      do 1 i=ilo,ihi
      do 2 j=jlo,jhi
      if(i.ne.j.and.latm(2,i).gt.1.and.latm(2,j).gt.1) then
      dist=sqrt((xatm(1,i)-xatm(1,j))*(xatm(1,i)-xatm(1,j))+
     + (xatm(2,i)-xatm(2,j))*(xatm(2,i)-xatm(2,j))+
     + (xatm(3,i)-xatm(3,j))*(xatm(3,i)-xatm(3,j)))
      if(dist.lt.pre_atsiz(latm(2,i))+pre_atsiz(latm(2,j))
     + .and.(ishort.eq.0.or.dist.lt.short)) then
      ishort=i
      jshort=j
      short=dist
      endif
      endif
    2 continue
    1 continue
c
      rshort=short
c
      pre_short=.false.
      if(ishort.gt.0) pre_short=.true.
c
      return
      end
      integer function pre_shorts(xatm,latm,matm,ilo,ihi,jlo,jhi,
     + ishort,jshort,rshort)
c
c $Id: pre_short.F,v 1.8 2001-06-11 22:01:06 d3j191 Exp $
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
      real*8 pre_atsiz
      external pre_atsiz
c
      integer matm
      real*8 xatm(3,matm)
      integer latm(5,matm)
      integer ilo,ihi,jlo,jhi,ishort(10),jshort(10)
      real*8 rshort(10)
c
      integer i,j,k,l,number,ishrt,jshrt
      real*8 dist,short
c
      number=0
c
      do 1 l=1,10
      ishrt=0
      short=0.0d0
      do 2 i=ilo,ihi
      do 3 j=jlo,jhi
      do 4 k=1,number
      if(ishort(k).eq.i.and.jshort(k).eq.j) goto 3
    4 continue
      if(i.ne.j.and.latm(2,i).gt.1.and.latm(2,j).gt.1) then
      dist=sqrt((xatm(1,i)-xatm(1,j))*(xatm(1,i)-xatm(1,j))+
     + (xatm(2,i)-xatm(2,j))*(xatm(2,i)-xatm(2,j))+
     + (xatm(3,i)-xatm(3,j))*(xatm(3,i)-xatm(3,j)))
      if(dist.lt.pre_atsiz(latm(2,i))+pre_atsiz(latm(2,j))
     + .and.(ishrt.eq.0.or.dist.lt.short)) then
      ishrt=i
      jshrt=j
      short=dist
      endif
      endif
    3 continue
    2 continue
      if(ishrt.gt.0) then
      number=number+1
      if(number.gt.10)
     + call errquit('Too many links between fragment pair',number) 
      ishort(number)=ishrt
      jshort(number)=jshrt
      rshort(number)=short
      endif
    1 continue
c
      pre_shorts=number
c
      return
      end
