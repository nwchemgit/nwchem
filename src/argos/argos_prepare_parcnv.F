      logical function argos_prepare_parcnv(icomb,latt,patt,catt,
     + matt,natt)
c
c $Id$
c
      implicit none
c
      integer icomb,matt,natt
      integer latt(2,matt)
      real*8 patt(4,2,matt,matt)
      character*6 catt(2,matt)
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
      patt(1,k,i,i)=2.0d0*epsil*((rstar+rstar)**6)
      patt(2,k,i,i)=epsil*((rstar+rstar)**12)
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
c     modify parameters for self and dummy
c     ------------------------------------
c
      do 5 i=1,natt
      do 6 j=i,natt
      if((catt(1,i)(6:6).eq.'S'.and.catt(1,j)(6:6).ne.'S').or.
     + (catt(1,i)(6:6).ne.'S'.and.catt(1,j)(6:6).eq.'S')) then
      patt(1,1,i,j)=0.0d0
      patt(1,2,i,j)=0.0d0
      patt(2,1,i,j)=0.0d0
      patt(2,2,i,j)=0.0d0
      patt(1,1,j,i)=0.0d0
      patt(1,2,j,i)=0.0d0
      patt(2,1,j,i)=0.0d0
      patt(2,2,j,i)=0.0d0
      endif
    6 continue
      if(catt(1,i)(6:6).eq.'D') then
      patt(1,1,i,i)=0.0d0
      patt(1,2,i,i)=0.0d0
      patt(2,1,i,i)=0.0d0
      patt(2,2,i,i)=0.0d0
      do 7 j=1,natt
      patt(1,1,i,j)=0.0d0
      patt(1,2,i,j)=0.0d0
      patt(2,1,i,j)=0.0d0
      patt(2,2,i,j)=0.0d0
      patt(1,1,j,i)=0.0d0
      patt(1,2,j,i)=0.0d0
      patt(2,1,j,i)=0.0d0
      patt(2,2,j,i)=0.0d0
    7 continue
      endif
      if(catt(1,i)(6:6).eq.'Q'.or.catt(1,i)(6:6).eq.'H') then
      patt(1,1,i,i)=0.0d0
      patt(1,2,i,i)=0.0d0
      patt(2,1,i,i)=0.0d0
      patt(2,2,i,i)=0.0d0
      do 8 j=1,natt
      if(catt(1,j)(6:6).eq.'Q'.or.catt(1,j)(6:6).eq.'H') then
      patt(1,1,i,j)=0.0d0
      patt(1,2,i,j)=0.0d0
      patt(2,1,i,j)=0.0d0
      patt(2,2,i,j)=0.0d0
      patt(1,1,j,i)=0.0d0
      patt(1,2,j,i)=0.0d0
      patt(2,1,j,i)=0.0d0
      patt(2,2,j,i)=0.0d0
      endif
    8 continue
      endif
    5 continue
c
      argos_prepare_parcnv=.true.
      return
      end
