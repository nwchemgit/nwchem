      subroutine argos_pme_fftest(grid1,grid2,nodep)
c
      implicit none
c
#include "argos_pme_common.fh"
#include "mafdecls.fh"
c
      real*8 grid1(2,mgx,mgy,mgz),grid2(2,mgx,mgy,mgz)
c
      integer nodep(nodpme)
      integer i,j,k
      intrinsic abs
c
c     initialize local portion of grid
c
      imfft = 1
      do 1 k=1,mgz
      do 2 j=1,ngy
      do 3 i=1,ngx
      grid1(1,i,j,k)=zero
      grid1(2,i,j,k)=zero
    3 continue
    2 continue
    1 continue
      if(me.eq.nodep(1)) then
      grid1(1,1,1,1)=one
      grid1(2,1,1,1)=two
      endif
      call argos_pme_fft(ngx,ngy,ngz,1,grid1,grid2,int_mb(i_map1),
     + int_mb(i_map2),int_mb(i_knod),imfft)
c
      if(ngzloc.gt.0) then
      do 4 k=1,ngzloc
      do 5 j=1,ngy
      do 6 i=1,ngx
      if(abs(grid2(1,i,j,k)-one).gt.small)
     + call md_abort('fftchk : forward fft failed',me)
      if(abs(grid2(2,i,j,k)-two).gt.small)
     + call md_abort('fftchk : forward fft failed',me)
    6 continue
    5 continue
    4 continue
      endif
c
      call argos_pme_fft(ngx,ngy,ngz,-1,grid2,grid1,int_mb(i_map1),
     + int_mb(i_map2),int_mb(i_knod),imfft)
c
      do 7 k=1,ngzloc
      do 8 j=1,ngy
      do 9 i=1,ngx
      if(me.eq.nodep(1).and.i+j+k.eq.3) then
      if(abs(grid1(1,i,j,k)-one).gt.small)
     + call md_abort('fftchk : reverse fft failed',me)
      if(abs(grid1(2,i,j,k)-two).gt.small)
     + call md_abort('fftchk : reverse fft failed',me)
      else
      if(abs(grid1(1,i,j,k)).gt.small)
     + call md_abort('fftchk : reverse fft failed',me)
      if(abs(grid1(2,i,j,k)).gt.small)
     + call md_abort('fftchk : reverse fft failed',me)
      endif
    9 continue
    8 continue
    7 continue
c
      return
      end
c $Id$
