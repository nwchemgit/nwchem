      logical function pre_excl(l3rd,m3rd,n3rd,lbnd,mbnd,nbnd,
     + lang,mang,nang,lexc,mexc,nexc,nval,ival,ndx)
c
c $Id: pre_excl.F,v 1.6 2000-06-21 03:48:15 d3j191 Exp $
c
      implicit none
c
      integer m3rd,n3rd,mbnd,nbnd,mang,nang,mexc,nexc
      integer l3rd(2,m3rd),lbnd(4,mbnd),lang(5,mang),lexc(2,mexc)
      integer nval
      integer ival(nval),ndx(nval)
c
      integer i,k,l,in,jn,ln,itmp,ir
c
      if(n3rd.gt.mexc) call errquit('increase mexc',9999)
      do 1 i=1,n3rd
      lexc(1,i)=l3rd(1,i)
      lexc(2,i)=l3rd(2,i)
    1 continue
      nexc=n3rd
c
      do 2 i=1,nbnd
      nexc=nexc+1
      if(nexc.gt.mexc) call errquit('increase mexc',9999)
      lexc(1,nexc)=lbnd(1,i)
      lexc(2,nexc)=lbnd(2,i)
    2 continue
      do 3 i=1,nang
      nexc=nexc+1
      if(nexc.gt.mexc) call errquit('increase mexc',9999)
      lexc(1,nexc)=min(lang(1,i),lang(3,i))
      lexc(2,nexc)=max(lang(1,i),lang(3,i))
    3 continue
c
      if(nexc.gt.0) then
      do 4 i=1,nexc
      ndx(i)=i
    4 continue
      ln=nexc/2+1
      ir=nexc
    5 continue
      if(ln.gt.1) then
      ln=ln-1
      itmp=ndx(ln)
      else
      itmp=ndx(ir)
      ndx(ir)=ndx(1)
      ir=ir-1
      if(ir.eq.1) then
      ndx(1)=itmp
      goto 6
      endif
      endif
      in=ln
      jn=ln+ln
    7 continue
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
      goto 7
      endif
      ndx(in)=itmp
      goto 5
    6 continue
c
      do 8 k=1,2
      do 9 i=1,nexc
      ival(i)=lexc(k,i)
    9 continue
      do 10 i=1,nexc
      lexc(k,i)=ival(ndx(i))
   10 continue
    8 continue
      endif
c
      if(nexc.gt.1) then
      l=nexc
      nexc=1
      do 11 i=2,l
      if(lexc(1,i).ne.lexc(1,nexc).or.lexc(2,i).ne.lexc(2,nexc)) then
      nexc=nexc+1
      lexc(1,nexc)=lexc(1,i)
      lexc(2,nexc)=lexc(2,i)
      endif
   11 continue
      endif
c
      pre_excl=.true.
      return
      end
