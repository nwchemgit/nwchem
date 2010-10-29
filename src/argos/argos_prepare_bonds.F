      logical function argos_prepare_bonds(xatm,latm,matm,ilo,ihi,lbnd,
     + mbnd,nbnd,
     + llnk,clnk,mlnk,nlnk)
c
c $Id$
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
      real*8 argos_prepare_atsiz
      external argos_prepare_atsiz
c
      integer matm,mbnd
      integer nbnd
      real*8 xatm(3,matm)
      integer latm(5,matm),lbnd(2,mbnd)
      integer ilo,ihi
      integer i,j
      real*8 dist
      integer mlnk,nlnk
      integer llnk(4,mlnk)
      character*4 clnk(2,mlnk)
      integer ilnk,jlnk
c
      argos_prepare_bonds=.false.
c
      do 1 i=ilo,ihi-1
      if(latm(2,i).gt.0) then
      do 2 j=i+1,ihi
      if(latm(2,j).gt.0) then
      dist=sqrt((xatm(1,j)-xatm(1,i))*(xatm(1,j)-xatm(1,i))+
     + (xatm(2,j)-xatm(2,i))*(xatm(2,j)-xatm(2,i))+
     + (xatm(3,j)-xatm(3,i))*(xatm(3,j)-xatm(3,i)))
      if(dist.gt.0.05d0.and.
     + dist.le.argos_prepare_atsiz(latm(2,
     + j))+argos_prepare_atsiz(latm(2,i))) then
      if(nbnd.ge.mbnd)
     + call md_abort('Dimension problem in argos_prepare_bonds',nbnd)
      nbnd=nbnd+1
      lbnd(1,nbnd)=i
      lbnd(2,nbnd)=j
      endif
      endif
    2 continue
      endif
    1 continue
c
      do 3 i=1,nlnk
      ilnk=llnk(3,i)
      jlnk=llnk(4,i)
      if(ilnk.gt.jlnk) then
      j=ilnk
      ilnk=jlnk
      jlnk=j
      endif
      if(ilnk.ge.ilo.and.jlnk.le.ihi) then
      do 4 j=1,nbnd
      if(ilnk.eq.lbnd(1,j).and.jlnk.eq.lbnd(2,j)) goto 3
    4 continue
      nbnd=nbnd+1
      lbnd(1,nbnd)=ilnk
      lbnd(2,nbnd)=jlnk
      endif
    3 continue
c
      do 5 i=1,nbnd-1
      do 6 j=i+1,nbnd
      if(lbnd(1,i).gt.lbnd(1,j)) then
      ilnk=lbnd(1,i)
      lbnd(1,i)=lbnd(1,j)
      lbnd(1,j)=ilnk
      ilnk=lbnd(2,i)
      lbnd(2,i)=lbnd(2,j)
      lbnd(2,j)=ilnk
      endif
    6 continue
    5 continue
c
      argos_prepare_bonds=.true.
      return
      end
