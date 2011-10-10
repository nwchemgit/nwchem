      subroutine argos_space_pack(numsa,isl,islp)
c
      implicit none
c
#include "argos_space_common.fh"
#include "global.fh"
#include "msgids.fh"
#include "bitops.fh"
c
      integer numsa,isl(msa,mis2),islp(msa,npack)
c
      integer i,j,k,ipck
c
      ipck=0
      do 1 j=1,npack
      do 2 i=1,numsa
      islp(i,j)=0
    2 continue
      do 3 k=1,ipack(j)
      ipck=ipck+1
      do 4 i=1,numsa
      islp(i,j)=lshift(islp(i,j),nbits(ipck))+isl(i,ipck)
    4 continue
    3 continue
    1 continue
c
      return
      end
c $Id$
