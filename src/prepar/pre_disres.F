      logical function pre_disres(card,isgm,csa,msa,nsa,inoe,dnoe,mnoe,
     + nnoe)
c
      implicit none
c
      integer mnoe,nnoe,msa,nsa
      integer isgm(msa),inoe(2,mnoe),is1,is2
      real*8 dnoe(5,mnoe),rval(5)
      character*16 csa(msa)
      character*4 nm1,nm2
      character*80 card
      integer i
c
      read(card,1000,end=9999,err=9999) is1,nm1,is2,nm2,rval
 1000 format(6x,i5,a4,i5,a4,3f12.6,2e12.5)
      if(is1.eq.0.or.is2.eq.0.or.nm1.eq.'    '.or.nm2.eq.'    ')
     + goto 9999
      nnoe=nnoe+1
      if(nnoe.gt.mnoe) call errquit('Increase mnoe',9999)
      inoe(1,nnoe)=0
      inoe(2,nnoe)=0
      do 1 i=1,nsa
      if(is1.eq.isgm(i).and.nm1.eq.csa(i)(11:14)) inoe(1,nnoe)=i
      if(is2.eq.isgm(i).and.nm2.eq.csa(i)(11:14)) inoe(2,nnoe)=i
    1 continue
      if(inoe(1,nnoe).eq.0.or.inoe(2,nnoe).eq.0)
     + call errquit('pre_disres error',9999)
      do 2 i=1,5
      dnoe(i,nnoe)=rval(i)
    2 continue
c
      pre_disres=.true.
      return
c
 9999 continue
      pre_disres=.false.
      return
      end
