      subroutine argos_space_gtda(jnode,jsfr,jsto,isl,islp,
     + isga,isq3,xs,ndxda)
c
      implicit none
c
#include "argos_space_common.fh"
#include "global.fh"
c
      integer jnode,jsfr,jsto,ndxda
      integer isl(msa,mis2),islp(msa,npack)
      integer isga(msa),isq3(msa)
      real*8 xs(msa,3)
c
      integer i,il,ih,jl,jh
c
      if(npack.eq.0) then
      call ga_distribution(ga_is,jnode,il,ih,jl,jh)
      call ga_get(ga_is,il+jsfr-1,il+jsto-1,jl,jl,
     + isga(ndxda+1),msa)
      call ga_get(ga_is,il+jsfr-1,il+jsto-1,jl+lsct3-1,jl+lsct3-1,
     + isq3(ndxda+1),msa)
      else
      call ga_distribution(ga_is,jnode,il,ih,jl,jh)
      call ga_get(ga_is,il+jsfr-1,il+jsto-1,jl,jl+npack-1,
     + islp(nsaloc+1,1),msa)
      call argos_space_unpack(jsto-jsfr+1,isl(nsaloc+1,1),
     + islp(nsaloc+1,1))
      do 1 i=1,jsto-jsfr
      isga(ndxda+i)=isl(nsaloc+i,lsgan)
      isq3(ndxda+i)=isl(nsaloc+i,lsct3)
    1 continue
      endif
      call ga_distribution(ga_s,jnode,il,ih,jl,jh)
      call ga_get(ga_s,il+jsfr-1,il+jsto-1,jl,jl+2,
     + xs(ndxda+1,1),msa)
c
      return
      end
c $Id$
