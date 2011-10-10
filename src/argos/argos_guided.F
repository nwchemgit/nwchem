      subroutine argos_guided(fw,fs,gw,gs)
c
      implicit none
c
#include "argos.fh"
#include "mafdecls.fh"
#include "global.fh"
c
      real*8 fw(mwm,3,mwa,2),fs(msa,3,2)
      real*8 gw(mwm,3,mwa),gs(msa,3)
c
      integer i,j,k
c
      do 1 k=1,nwa
      do 2 j=1,3
      do 3 i=1,nwmloc
      gw(i,j,k)=factgg*gw(i,j,k)+factgf*(fw(i,j,k,1)+fguide*gw(i,j,k))
      fw(i,j,k,1)=fw(i,j,k,1)+fguide*gw(i,j,k)
    3 continue
    2 continue
    1 continue
c
      do 4 j=1,3
      do 5 i=1,nsaloc
      gs(i,j)=factgg*gs(i,j)+factgf*(fs(i,j,1)+fguide*gs(i,j))
      fs(i,j,1)=fs(i,j,1)+fguide*gs(i,j)
    5 continue
    4 continue
c
      call argos_space_copyg(dbl_mb(i_gw),dbl_mb(i_gs))
      call ga_sync()
c
      return
      end
c $Id$
