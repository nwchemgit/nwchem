      logical function argos_prep_third(lbnd,mbnd,nbnd,lang,mang,nang,
     + l3rd,m3rd,n3rd,ib,ia,i3,nval,ival,ndx,natm)
c
c $Id$
c
      implicit none
c
      integer mbnd,nbnd,mang,nang,m3rd,n3rd
      integer ia,ib,i3
      integer lbnd(4,mbnd),lang(5,mang),l3rd(2,m3rd)
      integer nval,natm
      integer ival(nval),ndx(nval)
c
      integer i,j,k,l,it,jt,in,jn,ln,ir,itmp
      integer ii,ij,ik
c
      call timer_start(1)
      do 1 i=1,natm
      ival(i)=0
      ndx(i)=0
    1 continue
      do 2 i=ia,nang
      j=lang(2,i)
      if(ival(j).eq.0) ival(j)=i
      ndx(j)=i
    2 continue
      do 3 i=ia,nang
      ii=lang(1,i)
      ij=lang(2,i)
      ik=lang(3,i)
      if(ival(ii).gt.0) then
      do 4 j=ival(ii),ndx(ii)
      if(j.lt.i) then
      if(lang(1,j).eq.ij) then
      if(lang(3,j).ne.ik) then
      it=ik
      jt=lang(3,j)
      goto 5
      endif
      endif
      if(lang(3,j).eq.ij) then
      if(lang(1,j).ne.ik) then
      it=ik
      jt=lang(1,j)
      goto 5
      endif
      endif
      endif
      goto 4
    5 continue
      if(ival(it).gt.0) then
      do 6 k=ival(it),ndx(it)
      if(lang(1,k).eq.jt) goto 4
      if(lang(3,k).eq.jt) goto 4
    6 continue
      endif
      if(ival(jt).gt.0) then
      do 7 k=ival(jt),ndx(jt)
      if(lang(1,k).eq.it) goto 4
      if(lang(3,k).eq.it) goto 4
    7 continue
      if(ival(it).gt.0) then
      do 8 k=ival(it),ndx(it)
      do 9 l=ival(jt),ndx(jt)
      if(lang(1,k).eq.lang(1,l)) goto 4
      if(lang(1,k).eq.lang(3,l)) goto 4
      if(lang(3,k).eq.lang(1,l)) goto 4
      if(lang(3,k).eq.lang(3,l)) goto 4
    9 continue
    8 continue
      endif
      endif
      n3rd=n3rd+1
      l3rd(1,n3rd)=min(it,jt)
      l3rd(2,n3rd)=max(it,jt)
    4 continue
      endif
      if(ival(ik).gt.0) then
      do 10 j=ival(ik),ndx(ik)
      if(j.lt.i) then
      if(lang(1,j).eq.ij) then
      if(lang(3,j).ne.ii) then
      it=ii
      jt=lang(3,j)
      goto 11
      endif
      endif
      if(lang(3,j).eq.ij) then
      if(lang(1,j).ne.ii) then
      it=ii
      jt=lang(1,j)
      goto 11
      endif
      endif
      endif
      goto 10
   11 continue
      if(ival(it).gt.0) then
      do 12 k=ival(it),ndx(it)
      if(lang(1,k).eq.jt) goto 10
      if(lang(3,k).eq.jt) goto 10
   12 continue
      endif
      if(ival(jt).gt.0) then
      do 13 k=ival(jt),ndx(jt)
      if(lang(1,k).eq.it) goto 10
      if(lang(3,k).eq.it) goto 10
   13 continue
      if(ival(it).gt.0) then
      do 14 k=ival(it),ndx(it)
      do 15 l=ival(jt),ndx(jt)
      if(lang(1,k).eq.lang(1,l)) goto 10
      if(lang(1,k).eq.lang(3,l)) goto 10
      if(lang(3,k).eq.lang(1,l)) goto 10
      if(lang(3,k).eq.lang(3,l)) goto 10
   15 continue
   14 continue
      endif
      endif
      n3rd=n3rd+1
      l3rd(1,n3rd)=min(it,jt)
      l3rd(2,n3rd)=max(it,jt)
   10 continue
      endif
    3 continue
c
      if(n3rd.gt.nval)
     + call md_abort('argos_prep_third: increase dimension nval',n3rd)
      if(n3rd.gt.0) then
      do 16 i=1,n3rd
      ndx(i)=i
   16 continue
      ln=n3rd/2+1
      ir=n3rd
   17 continue
      if(ln.gt.1) then
      ln=ln-1
      itmp=ndx(ln)
      else
      itmp=ndx(ir)
      ndx(ir)=ndx(1)
      ir=ir-1
      if(ir.eq.1) then
      ndx(1)=itmp
      goto 18
      endif
      endif
      in=ln
      jn=ln+ln
   19 continue
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
      goto 19
      endif
      ndx(in)=itmp
      goto 17
   18 continue
c
      do 20 k=1,2
      do 21 i=1,n3rd
      ival(i)=l3rd(k,i)
   21 continue
      do 22 i=1,n3rd
      l3rd(k,i)=ival(ndx(i))
   22 continue
   20 continue
      endif
c
      if(n3rd.gt.1) then
      l=n3rd
      n3rd=1
      do 23 i=2,l
      if(l3rd(1,i).ne.l3rd(1,n3rd).or.l3rd(2,i).ne.l3rd(2,n3rd)) then
      n3rd=n3rd+1
      l3rd(1,n3rd)=l3rd(1,i)
      l3rd(2,n3rd)=l3rd(2,i)
      endif
   23 continue
      endif
c
      argos_prep_third=.true.
      return
      end
