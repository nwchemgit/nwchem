      logical function pre_third(lbnd,mbnd,nbnd,lang,mang,nang,
     + l3rd,m3rd,n3rd,ib,ia,i3)
c
c $Id: pre_third.F,v 1.3 1998-04-29 18:35:28 d3j191 Exp $
c
      implicit none
c
      integer mbnd,nbnd,mang,nang,m3rd,n3rd
      integer ia,ib,i3
      integer lbnd(4,mbnd),lang(5,mang),l3rd(2,m3rd)
c
      integer i,j,k,l,it,jt
c
      do 1 i=ia,nang-1
      do 2 j=i+1,nang
      if(lang(2,i).eq.lang(1,j).and.lang(3,i).eq.lang(2,j)) then
      it=lang(1,i)
      jt=lang(3,j)
      else
      if(lang(2,i).eq.lang(3,j).and.lang(3,i).eq.lang(2,j)) then
      it=lang(1,i)
      jt=lang(1,j)
      else
      if(lang(2,i).eq.lang(1,j).and.lang(1,i).eq.lang(2,j)) then
      it=lang(3,i)
      jt=lang(3,j)
      else
      if(lang(2,i).eq.lang(3,j).and.lang(1,i).eq.lang(2,j)) then
      it=lang(3,i)
      jt=lang(1,j)
      else
      goto 2
      endif
      endif
      endif
      endif
      if(it.gt.jt) then
      l=it
      it=jt
      jt=l
      endif
      do 3 k=i3,n3rd
      if(l3rd(1,k).eq.it.and.l3rd(2,k).eq.jt) goto 2
    3 continue
      do 4 k=ib,nbnd
      if(lbnd(1,k).eq.it.and.lbnd(2,k).eq.jt) goto 2
      if(lbnd(1,k).eq.jt.and.lbnd(2,k).eq.it) goto 2
    4 continue
      do 5 k=ia,nang
      if(lang(1,k).eq.it.and.lang(3,k).eq.jt) goto 2
      if(lang(1,k).eq.jt.and.lang(3,k).eq.it) goto 2
    5 continue
      n3rd=n3rd+1
      if(n3rd.gt.m3rd) call errquit('increase m3rd',9999)
      l3rd(1,n3rd)=it
      l3rd(2,n3rd)=jt
    2 continue
    1 continue
      do 6 i=i3,n3rd-1
      do 7 j=i+1,n3rd
      if(l3rd(1,i).gt.l3rd(1,j).or.
     + (l3rd(1,i).eq.l3rd(1,j).and.l3rd(2,i).gt.l3rd(2,j))) then
      l=l3rd(1,i)
      l3rd(1,i)=l3rd(1,j)
      l3rd(1,j)=l
      l=l3rd(2,i)
      l3rd(2,i)=l3rd(2,j)
      l3rd(2,j)=l
      endif
    7 continue
    6 continue
c
      pre_third=.true.
      return
      end
