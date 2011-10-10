      subroutine argos_diana_pcovar(ndex,cov,isel)
c
      implicit none
c
#include "argos_diana_params.fh"
#include "argos_diana_common.fh"
#include "global.fh"
#include "msgids.fh"
c
      integer ndex(*)
      real*8 cov(*)
      integer isel(msa)
c
      integer i,j,k,ia,ja,ix,jx
      integer iproc
c
      if(me.ne.0) return
c
      write(*,'(a)') ' Covariance Matrix'
      j=0
      do 1 i=1,nsa
      if(isel(i).ne.0) then
      j=j+1
      ndex(j)=i
      endif
    1 continue
c
      return
      do 6 iproc=0,np-1
      call ga_distribution(ga_cov,iproc,il,ih,jl,jh)
      do 3 i=il,ih
      call ga_get(ga_cov,i,i,jl,jh,cov,1)
      ia=ndex((i+2)/3)
      ix=mod(i+2,3)+1
      do 4 j=jl,jh
      ja=ndex((j+2)/3)
      jx=mod(j+2,3)+1
      write(*,'(4i5,f12.6)') i,j,ia,ja,cov(j-jl+1)
    4 continue
    3 continue
    6 continue
c
      return
      end
c $Id$
