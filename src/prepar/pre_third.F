      logical function pre_third(lbnd,mbnd,nbnd,lang,mang,nang,
     + l3rd,m3rd,n3rd,ib,ia,i3,nval,ival,ndx)
c
c $Id: pre_third.F,v 1.5 2000-06-15 07:16:25 d3j191 Exp $
c
      implicit none
c
      integer mbnd,nbnd,mang,nang,m3rd,n3rd
      integer ia,ib,i3
      integer lbnd(4,mbnd),lang(5,mang),l3rd(2,m3rd)
      integer nval
      integer ival(nval),ndx(nval)
c
      integer i,j,k,l,it,jt,in,jn,ln,ir,itmp
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
c
      if(n3rd.gt.0) then
      do 6 i=1,n3rd
      ndx(i)=i
    6 continue
      ln=n3rd/2+1
      ir=n3rd
    7 continue
      if(ln.gt.1) then
      ln=ln-1
      itmp=ndx(ln)
      else
      itmp=ndx(ir)
      ndx(ir)=ndx(1)
      ir=ir-1
      if(ir.eq.1) then
      ndx(1)=itmp
      goto 8
      endif
      endif
      in=ln
      jn=ln+ln
    9 continue
      if(jn.le.ir) then
      if(jn.lt.ir) then
      if(l3rd(1,ndx(jn)).gt.l3rd(1,ndx(jn+1)).or.
     + (l3rd(1,ndx(jn)).eq.l3rd(1,ndx(jn+1)).and.
     + l3rd(2,ndx(jn)).gt.l3rd(2,ndx(jn+1)))) then
      else
      jn=jn+1
      endif
      endif
      if(l3rd(1,itmp).gt.l3rd(1,ndx(jn)).or.
     + (l3rd(1,itmp).eq.l3rd(1,ndx(jn)).and.
     + l3rd(2,itmp).gt.l3rd(2,ndx(jn)))) then
      jn=ir+1
      else
      ndx(in)=ndx(jn)
      in=jn
      jn=jn+jn
      endif
      goto 9
      endif
      ndx(in)=itmp
      goto 7
    8 continue
c
      do 10 k=1,2
      do 11 i=1,n3rd
      ival(i)=l3rd(k,i)
   11 continue
      do 12 i=1,n3rd
      l3rd(k,i)=ival(ndx(i))
   12 continue
   10 continue
      endif
c
      pre_third=.true.
      return
      end
