      logical function pre_excl(l3rd,m3rd,n3rd,lbnd,mbnd,nbnd,
     + lang,mang,nang,lexc,mexc,nexc,nval,ival,ndx)
c
c $Id: pre_excl.F,v 1.5 2000-06-15 07:16:25 d3j191 Exp $
c
      implicit none
c
      integer m3rd,n3rd,mbnd,nbnd,mang,nang,mexc,nexc
      integer l3rd(2,m3rd),lbnd(4,mbnd),lang(5,mang),lexc(2,mexc)
      integer nval
      integer ival(nval),ndx(nval)
c
      integer i,j,k,l,in,jn,ln,itmp,ir
c
      if(n3rd.gt.mexc) call errquit('increase mexc',9999)
      do 1 i=1,n3rd
      lexc(1,i)=l3rd(1,i)
      lexc(2,i)=l3rd(2,i)
    1 continue
      nexc=n3rd
c
      do 2 i=1,nbnd
      do 3 j=1,nexc
      if(lexc(1,j).eq.lbnd(1,i).and.lexc(2,j).eq.lbnd(2,i)) goto 2
    3 continue
      nexc=nexc+1
      if(nexc.gt.mexc) call errquit('increase mexc',9999)
      lexc(1,nexc)=lbnd(1,i)
      lexc(2,nexc)=lbnd(2,i)
    2 continue
c
      do 4 i=1,nang
      do 5 j=1,nexc
      if(lexc(1,j).eq.lang(1,i).and.lexc(2,j).eq.lang(3,i)) goto 4
    5 continue
      nexc=nexc+1
      if(nexc.gt.mexc) call errquit('increase mexc',9999)
      lexc(1,nexc)=lang(1,i)
      lexc(2,nexc)=lang(3,i)
    4 continue
c
      if(nexc.gt.0) then
      do 6 i=1,nexc
      ndx(i)=i
    6 continue
      ln=nexc/2+1
      ir=nexc
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
      if(lexc(1,ndx(jn)).gt.lexc(1,ndx(jn+1)).or.
     + (lexc(1,ndx(jn)).eq.lexc(1,ndx(jn+1)).and.
     + lexc(2,ndx(jn)).gt.lexc(2,ndx(jn+1)))) then
      else
      jn=jn+1
      endif
      endif
      if(lexc(1,itmp).gt.lexc(1,ndx(jn)).or.
     + (lexc(1,itmp).eq.lexc(1,ndx(jn)).and.
     + lexc(2,itmp).gt.lexc(2,ndx(jn)))) then
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
      do 11 i=1,nexc
      ival(i)=lexc(k,i)
   11 continue
      do 12 i=1,nexc
      lexc(k,i)=ival(ndx(i))
   12 continue
   10 continue
      endif
c
      pre_excl=.true.
      return
      end
