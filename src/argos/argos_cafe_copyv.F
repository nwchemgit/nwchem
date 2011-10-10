      subroutine argos_cafe_copyv(ia,ja,vdwai,vdw)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "mafdecls.fh"
c
      integer ia,ja
      real*8 vdwai(24,4),vdw(mat,mat,map,mset)
c
      integer i,j
c
      if(.not.lfree) then
      do 1 i=1,4
      do 2 j=1,nparms
      vdw(ia,ja,i,j)=vdwai(j,i)
      if(ia.ne.ja) vdw(ja,ia,i,j)=vdw(ia,ja,i,j)
    2 continue
    1 continue
      else
      do 3 i=1,4
      do 4 j=1,3
      vdw(ia,ja,i,j)=vdwai(j,i)
    4 continue
    3 continue
      do 5 i=1,4
      vdw(ia,ja,i,4)=vdw(ia,ja,i,3)-vdw(ia,ja,i,2)
      vdw(ia,ja,i,5)=vdw(ia,ja,i,2)
      vdw(ia,ja,i,6)=vdw(ia,ja,i,3)
    5 continue
      if(ia.ne.ja) then
      do 6 i=1,4
      do 7 j=1,6
      vdw(ja,ia,i,j)=vdw(ia,ja,i,j)
    7 continue
    6 continue
      endif
      endif
c
c      do 55 i=1,4
c      write(*,9) ia,ja,i,(vdw(ia,ja,i,j),j=1,mset)
c   55 continue
c    9 format(3i5,10e12.5)
      return
      end
c $Id$
