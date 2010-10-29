      logical function argos_prep_ang(lbnd,mbnd,nbnd,lang,mang,nang)
c
c $Id$
c
c     function to extend the list of angles from a list of bonds
c
c     in  : lbnd(2,mbnd) = list of bond indices
c         : mbnd         = dimension of the bond index list
c         : nbnd         = number of bonds
c         : mang         = dimension of the angle index list
c     out : lang(3,mang) = list of angle indices
c         : nang         = number of angles
c
c     function returns .false. if mang < nang
c
      implicit none
c
      integer mbnd,mang
      integer nbnd,nang
      integer lbnd(2,mbnd),lang(3,mang)
      integer i,j,k,ltemp
c
c     determine the list of angles
c     ----------------------------
c
      argos_prep_ang=.false.
      do 1 i=1,nbnd-1
      do 2 j=i+1,nbnd
      if(lbnd(1,i).eq.lbnd(1,j)) then
      if(nang.ge.mang) return
      nang=nang+1
      lang(1,nang)=lbnd(2,i)
      lang(2,nang)=lbnd(1,i)
      lang(3,nang)=lbnd(2,j)
      else
      if(lbnd(2,i).eq.lbnd(1,j)) then
      if(nang.ge.mang) return
      nang=nang+1
      lang(1,nang)=lbnd(1,i)
      lang(2,nang)=lbnd(2,i)
      lang(3,nang)=lbnd(2,j)
      else
      if(lbnd(1,i).eq.lbnd(2,j)) then
      if(nang.ge.mang) return
      nang=nang+1
      lang(1,nang)=lbnd(2,i)
      lang(2,nang)=lbnd(1,i)
      lang(3,nang)=lbnd(1,j)
      else
      if(lbnd(2,i).eq.lbnd(2,j)) then
      if(nang.ge.mang) return
      nang=nang+1
      lang(1,nang)=lbnd(1,i)
      lang(2,nang)=lbnd(2,i)
      lang(3,nang)=lbnd(1,j)
      endif
      endif
      endif
      endif
    2 continue
    1 continue
c
      do 3 i=1,nang
      if(lang(1,i).gt.lang(3,i)) then
      ltemp=lang(1,i)
      lang(1,i)=lang(3,i)
      lang(3,i)=ltemp
      endif
    3 continue
c
      do 4 i=1,nang-1
      do 5 j=i+1,nang
      if(lang(2,i).gt.lang(2,j).or.
     + (lang(2,i).eq.lang(2,j).and.lang(1,i).gt.lang(1,j)).or.
     + (lang(2,i).eq.lang(2,j).and.lang(1,i).eq.lang(1,j).and.
     + lang(3,i).gt.lang(3,j))) then
      do 6 k=1,3
      ltemp=lang(k,i)
      lang(k,i)=lang(k,j)
      lang(k,j)=ltemp
    6 continue
      endif
    5 continue
    4 continue
c
      argos_prep_ang=.true.
      return
      end
