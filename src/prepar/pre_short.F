      logical function pre_short(xatm,latm,matm,ilo,ihi,jlo,jhi,
     + ishort,jshort,rshort)
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
      real*8 pre_atrad
      external pre_atrad
c
      integer matm
      real*8 xatm(3,matm)
      integer latm(5,matm)
      integer ilo,ihi,jlo,jhi,ishort,jshort
      real*8 rshort
c
      integer i,j
      real*8 dist2,short2
c
      ishort=0
c
      do 1 i=ilo,ihi
      do 2 j=jlo,jhi
      if(i.ne.j) then
      dist2=(xatm(1,i)-xatm(1,j))*(xatm(1,i)-xatm(1,j))+
     + (xatm(2,i)-xatm(2,j))*(xatm(2,i)-xatm(2,j))+
     + (xatm(3,i)-xatm(3,j))*(xatm(3,i)-xatm(3,j))
      if(ishort.eq.0.or.dist2.lt.short2) then
      ishort=i
      jshort=j
      short2=dist2
      endif
      endif
    2 continue
    1 continue
c
      pre_short=.false.
c
      if(ishort.gt.0.and.
     + latm(2,ishort).gt.0.and.latm(2,jshort).gt.0) then
      rshort=sqrt(short2)
      if(rshort.lt.pre_atrad(latm(2,ishort))+pre_atrad(latm(2,jshort)))
     + pre_short=.true.
      endif
c
      return
      end
