      logical function argos_prep_transl(itran,jtran,
     + xtran,is,xs,msa,nsa)
c
c $Id$
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
      argos_prep_transl=.false.
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
      argos_prep_transl=.true.
c
      return
      end
