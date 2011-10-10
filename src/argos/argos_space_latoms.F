      subroutine argos_space_latoms(isl,lst,nums)
c
      implicit none
c
#include "argos_space_common.fh"
#include "global.fh"
#include "msgids.fh"
c
      integer nums
      integer isl(msa,mis2),lst(nsa)
c
      integer i
c
      do 1 i=1,nsa
      lst(i)=0
    1 continue
c
      call argos_cafe_links(lquant,nums,isl(1,lsgan),isl(1,lsdyn),lst)
c
      call ga_igop(msp_21,lst,nsa,'+')
c
      do 2 i=1,nums
      if(iand(isl(i,lsdyn),lquant).eq.0.and.lst(isl(i,lsgan)).gt.0)
     + isl(i,lsdyn)=ior(isl(i,lsdyn),lqlink)
    2 continue
c
      do 3 i=1,nsa
      lst(i)=0
    3 continue
c
      call argos_cafe_links(lqlink,nums,isl(1,lsgan),isl(1,lsdyn),lst)
c
      call ga_igop(msp_26,lst,nsa,'+')
c
      do 4 i=1,nums
      if(iand(isl(i,lsdyn),lqlink).eq.0.and.lst(isl(i,lsgan)).gt.0)
     + isl(i,lsdyn)=ior(isl(i,lsdyn),lqsoft)
    4 continue
c
      return
      end
c $Id$
