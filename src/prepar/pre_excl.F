      logical function pre_excl(l3rd,m3rd,n3rd,lbnd,mbnd,nbnd,
     + lang,mang,nang,lexc,mexc,nexc,nval,ival,ndx,catt,matt,
     + lats,nparms,mats,latm,matm,natm,lupdat)
c
c $Id$
c
      implicit none
c
      integer m3rd,n3rd,mbnd,nbnd,mang,nang,mexc,nexc
      integer matt,mats,matm,natm
      integer l3rd(2,m3rd),lbnd(4,mbnd),lang(5,mang),lexc(2,mexc)
      integer nval,nparms
      integer ival(nval),ndx(nval)
      character*6 catt(2,matt)
      integer lats(nparms,mats)
      integer latm(11,matm)
      logical lupdat
c
      integer i,j,k,l,in,jn,ln,itmp,ir
c
      if(n3rd.gt.mexc) call md_abort('increase mexc',9999)
      do 1 i=1,n3rd
      lexc(1,i)=l3rd(1,i)
      lexc(2,i)=l3rd(2,i)
    1 continue
      nexc=n3rd
c
      do 2 i=1,nbnd
      nexc=nexc+1
      if(nexc.gt.mexc) call md_abort('increase mexc',9999)
      lexc(1,nexc)=lbnd(1,i)
      lexc(2,nexc)=lbnd(2,i)
    2 continue
      do 3 i=1,nang
      nexc=nexc+1
      if(nexc.gt.mexc) call md_abort('increase mexc',9999)
      lexc(1,nexc)=min(lang(1,i),lang(3,i))
      lexc(2,nexc)=max(lang(1,i),lang(3,i))
    3 continue
c
      if(lupdat.and.nparms.eq.3) then
      do 4 i=1,natm-1
      if(catt(1,lats(2,latm(3,i)))(6:6).eq.'D') then
      do 5 j=i+1,natm
      if(catt(1,lats(3,latm(3,j)))(6:6).eq.'D') then
      nexc=nexc+1
      if(nexc.gt.mexc) call md_abort('increase mexc',9999)
      lexc(1,nexc)=i
      lexc(2,nexc)=j
      write(*,'(a,2i5)') 'Excluding ',i,j
      endif
    5 continue
      elseif(catt(1,lats(3,latm(3,i)))(6:6).eq.'D') then
      do 6 j=i+1,natm
      if(catt(1,lats(2,latm(3,j)))(6:6).eq.'D') then
      nexc=nexc+1
      if(nexc.gt.mexc) call md_abort('increase mexc',9999)
      lexc(1,nexc)=i
      lexc(2,nexc)=j
      write(*,'(a,2i5)') 'Excluding ',i,j
      endif
    6 continue
      endif
    4 continue
      endif
c
      if(nexc.gt.1) then
      do 7 i=1,nexc
      ndx(i)=i
    7 continue
      ln=nexc/2+1
      ir=nexc
    8 continue
      if(ln.gt.1) then
      ln=ln-1
      itmp=ndx(ln)
      else
      itmp=ndx(ir)
      ndx(ir)=ndx(1)
      ir=ir-1
      if(ir.eq.1) then
      ndx(1)=itmp
      goto 9
      endif
      endif
      in=ln
      jn=ln+ln
   10 continue
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
      goto 10
      endif
      ndx(in)=itmp
      goto 8
    9 continue
c
      do 11 k=1,2
      do 12 i=1,nexc
      ival(i)=lexc(k,i)
   12 continue
      do 13 i=1,nexc
      lexc(k,i)=ival(ndx(i))
   13 continue
   11 continue
      endif
c
      if(nexc.gt.1) then
      l=nexc
      nexc=1
      do 14 i=2,l
      if(lexc(1,i).ne.lexc(1,nexc).or.lexc(2,i).ne.lexc(2,nexc)) then
      nexc=nexc+1
      lexc(1,nexc)=lexc(1,i)
      lexc(2,nexc)=lexc(2,i)
      endif
   14 continue
      endif
c
      pre_excl=.true.
      return
      end
