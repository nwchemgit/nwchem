      logical function pre_parcnv(icomb,latt,patt,matt,natt)
c
      implicit none
c
      integer icomb,matt,natt
      integer latt(2,matt)
      real*8 patt(4,2,matt,matt)
c
      integer i,j,k
      real*8 epsil,rstar
c
c     C6,C12 -> C6,C12
c
      if(icomb.eq.1) then
      endif
c
c     epsilon,sigma  -> C6,C12
c
      if(icomb.eq.2) then
      endif
c
c     epsilon,(1/2)R* -> C6,C12
c
      if(icomb.eq.3) then
      do 1 i=1,natt
      do 2 k=1,2
      epsil=patt(1,k,i,i)
      rstar=patt(2,k,i,i)
      patt(1,k,i,i)=2.0d0*epsil*epsil*((rstar+rstar)**6)
      patt(2,k,i,i)=epsil*epsil*((rstar+rstar)**12)
    2 continue
      do 3 j=i+1,natt
      if(i.ne.j) then
      do 4 k=1,2
      epsil=patt(1,k,i,j)
      rstar=patt(2,k,i,j)
      patt(1,k,i,j)=2.0d0*epsil*((rstar)**6)
      patt(1,k,j,i)=patt(1,k,i,j)
      patt(2,k,i,j)=epsil*((rstar)**12)
      patt(2,k,j,i)=patt(2,k,i,j)
    4 continue
      endif
    3 continue
    1 continue
      endif
c
      pre_parcnv=.true.
      return
      end
