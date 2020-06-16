      subroutine argos_space_unpack(numsa,isl,islp)
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
      ipck=mis2+1
      do 1 j=npack,1,-1
      do 2 k=ipack(j),1,-1
      ipck=ipck-1
      do 3 i=1,numsa
      isl(i,ipck)=iand(islp(i,j),misl(ipck))
      islp(i,j)=ishft(islp(i,j)-isl(i,ipck),-nbits(ipck))
    3 continue
    2 continue
    1 continue
c
      return
      end
c $Id$
