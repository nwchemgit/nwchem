      logical function pre_transl(itran,jtran,xtran,is,xs,msa,nsa)
c
c $Id: pre_transl.F,v 1.2 2002-04-10 03:13:34 d3j191 Exp $
c
      implicit none
c
      integer msa,nsa,itran,jtran
      integer is(msa)
      real*8 xs(3,msa),xtran(3)
c
      integer i,j
c
      if(nsa.le.0.or.itran.le.0.or.jtran.gt.nsa) then
      pre_transl=.false.
      return
      endif
c
      do 1 i=1,nsa
      if(is(i).ge.itran.and.is(i).le.jtran) then
      do 2 j=1,3
      xs(j,i)=xs(j,i)+xtran(j)
    2 continue
      endif
    1 continue
c
      pre_transl=.true.
c
      return
      end
