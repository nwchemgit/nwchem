      logical function pre_charge(latm,qatm,matm,lbnd,mbnd,nbnd)
c
      implicit none
c
c     in  : latm(2,matm)      = atomic number
c           matm              = dimension atom list
c           lbnd(1:2,mbnd)    = bond indices
c           mbnd              = dimension bond list
c           nbnd              = length bond list
c
c     out : qatm(matm)        = guestimated atomic partial charges
c
      integer matm
      integer latm(5,matm)
      real*8 qatm(matm)
      integer mbnd,nbnd
      integer lbnd(2,mbnd)
c
      integer i,ia,ja,it,jt,itmp
c
      do 1 i=1,matm
      qatm(i)=0.0d0
    1 continue
c
      do 2 i=1,nbnd
      ia=lbnd(1,i)
      ja=lbnd(2,i)
      it=latm(2,ia)
      jt=latm(2,ja)
      if(it.ne.jt) then
c
c     smaller atom number in it/ia
c
      if(it.gt.jt) then
      itmp=ia
      ia=ja
      ja=itmp
      itmp=it
      it=jt
      jt=itmp
      endif
c
c     H - C : 0.05  : -0.05
c
      if(it.eq.1.and.jt.eq.6) then
      qatm(ia)=qatm(ia)+0.05d0
      qatm(ja)=qatm(ja)-0.05d0
      endif
c
c     H - N : 0.27  : -0.27
c
      if(it.eq.1.and.jt.eq.7) then
      qatm(ia)=qatm(ia)+0.27d0
      qatm(ja)=qatm(ja)-0.27d0
      endif
c
c     H - S : 0.19  : -0.19
c
      if(it.eq.1.and.jt.eq.16) then
      qatm(ia)=qatm(ia)+0.19d0
      qatm(ja)=qatm(ja)-0.19d0
      endif
c
c     C - N : 0.07  : -0.07
c
      if(it.eq.6.and.jt.eq.7) then
      qatm(ia)=qatm(ia)+0.07d0
      qatm(ja)=qatm(ja)-0.07d0
      endif
c
c     C - O : 0.30  : -0.30
c
      if(it.eq.6.and.jt.eq.8.and.latm(3,ja).eq.2) then
      qatm(ia)=qatm(ia)+0.30d0
      qatm(ja)=qatm(ja)-0.30d0
      endif
c
c     C = O : 0.57  : -0.57
c
      if(it.eq.6.and.jt.eq.8.and.latm(3,ja).eq.1) then
      qatm(ia)=qatm(ia)+0.57d0
      qatm(ja)=qatm(ja)-0.57d0
      endif
c
c     C - S : 0.11  : -0.11
c
      if(it.eq.6.and.jt.eq.16) then
      qatm(ia)=qatm(ia)+0.11d0
      qatm(ja)=qatm(ja)-0.11d0
      endif
c
      endif
    2 continue
c
      pre_charge=.true.
      return
c
 9999 pre_charge=.false.
      return
      end
