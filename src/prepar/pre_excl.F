      logical function pre_excl(l3rd,m3rd,n3rd,lbnd,mbnd,nbnd,
     + lang,mang,nang,lexc,mexc,nexc)
c
      implicit none
c
      integer m3rd,n3rd,mbnd,nbnd,mang,nang,mexc,nexc
      integer l3rd(2,m3rd),lbnd(4,mbnd),lang(5,mang),lexc(2,mexc)
c
      integer i,j,l
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
      if(lexc(1,j).eq.lang(1,i).and.lexc(2,j).eq.lang(3,i)) goto 2
    5 continue
      nexc=nexc+1
      if(nexc.gt.mexc) call errquit('increase mexc',9999)
      lexc(1,nexc)=lang(1,i)
      lexc(2,nexc)=lang(3,i)
    4 continue
c
      do 6 i=1,nexc-1
      do 7 j=i+1,nexc
      if(lexc(1,i).gt.lexc(1,j).or.
     + (lexc(1,i).eq.lexc(1,j).and.lexc(2,i).gt.lexc(2,j))) then
      l=lexc(1,i)
      lexc(1,i)=lexc(1,j)
      lexc(1,j)=l
      l=lexc(2,i)
      lexc(2,i)=lexc(2,j)
      lexc(2,j)=l
      endif
    7 continue
    6 continue
c
      pre_excl=.true.
      return
      end
