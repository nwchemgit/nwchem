      logical function argos_prepare_short(xatm,latm,matm,ilo,ihi,jlo,
     + jhi,
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
      real*8 argos_prepare_atsiz
      external argos_prepare_atsiz
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
      if(dist.lt.argos_prepare_atsiz(latm(2,
     + i))+argos_prepare_atsiz(latm(2,j))
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
      argos_prepare_short=.false.
      if(ishort.gt.0) argos_prepare_short=.true.
c
      return
      end
