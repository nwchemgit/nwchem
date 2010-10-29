      logical function pre_ang(lbnd,mbnd,nbmd,lang,mang,nang)
c
c $Id$
c
      implicit none
c
      integer mbnd,mang
      integer nbnd,nang
      integer lbnd(mbnd,2),lang(mang,3)
c
c     determine the list of angles
c     ----------------------------
c
      nang=0
      do 1 i=1,nbnd-1
      do 2 j=i+1,nbnd
      if(lbnd(i,1).eq.lbnd(j,1)) then
      nang=nang+1
      if(nang.gt.mang) goto 9999
      lang(nang,1)=lbnd(i,2)
      lang(nang,2)=lbnd(i,1)
      lang(nang,3)=lbnd(j,2)
      else
      if(lbnd(i,2).eq.lbnd(j,1)) then
      nang=nang+1
      if(nang.gt.mang) goto 9999
      lang(nang,1)=lbnd(i,1)
      lang(nang,2)=lbnd(i,2)
      lang(nang,3)=lbnd(j,2)
      else
      if(lbnd(i,1).eq.lbnd(j,2)) then
      nang=nang+1
      if(nang.gt.mang) goto 9999
      lang(nang,1)=lbnd(i,2)
      lang(nang,2)=lbnd(i,1)
      lang(nang,3)=lbnd(j,1)
      else
      if(lbnd(i,2).eq.lbnd(j,2)) then
      nang=nang+1
      if(nang.gt.mang) goto 9999
      lang(nang,1)=lbnd(i,1)
      lang(nang,2)=lbnd(i,2)
      lang(nang,3)=lbnd(j,1)
      endif
      endif
      endif
      endif
    2 continue
    1 continue
c
      pre_angl=.true.
      return
c
 9999 continue
      pre_angl=.false.
      return
      end
