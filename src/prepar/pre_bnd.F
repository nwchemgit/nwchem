      logical function pre_bnd(xatm,latm,matm,ilo,ihi,lbnd,mbnd,nbnd)
c
c     function to append the list of bonds from atomic coordinates 
c              for all atoms i in [ilo,ihi]
c     
c     in  : xatm(3,matm) = atomic coordinates
c         : latm(2,matm) = atomic numbers
c         : matm         = dimension of coordinate list
c         : natm         = number of atoms
c         : mbnd         = dimension of the bond index list
c     out : lbnd(2,mbnd) = list of bond indices
c         : nbnd         = number of bonds
c         : mang         = dimension of the angle index list
c
c     function returns .false. if mbnd < nbnd
c
      implicit none
c
      real*8 pre_atrad
      external pre_atrad
c
      integer matm,mbnd
      integer natm,nbnd
      real*8 xatm(3,matm)
      integer latm(5,matm),lbnd(2,mbnd)
      integer ilo,ihi
      integer i,j,k,ltemp
      real*8 dist
c
      pre_bnd=.false.
c
      do 1 i=ilo,ihi-1
      if(latm(2,i).gt.0) then
      do 2 j=i+1,ihi
      if(latm(2,j).gt.0) then
      dist=sqrt((xatm(1,j)-xatm(1,i))*(xatm(1,j)-xatm(1,i))+
     + (xatm(2,j)-xatm(2,i))*(xatm(2,j)-xatm(2,i))+
     + (xatm(3,j)-xatm(3,i))*(xatm(3,j)-xatm(3,i)))
      if(dist.le.pre_atrad(latm(2,j))+pre_atrad(latm(2,i))) then
      if(nbnd.ge.mbnd) return
      nbnd=nbnd+1
      lbnd(1,nbnd)=i
      lbnd(2,nbnd)=j
      endif
      endif
    2 continue
      endif
    1 continue 
c
      pre_bnd=.true.
      return
      end
