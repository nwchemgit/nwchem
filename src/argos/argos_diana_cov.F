      subroutine argos_diana_cov(ndex,cov,isel,xtrj,xa)
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
      real*8 xtrj(ntrj,3,nfrdim),xa(nsa,3)
      integer isel(msa)
c
      integer i,j,k,ia,ja,ix,jx
      integer ilc,ihc,jlc,jhc
      integer iproc
c
c      write(*,'(a)') 'ARGOS_DIANA_COVAR'
c
      j=0
      do 1 i=1,nsa
      if(isel(i).ne.0) then
      j=j+1
      ndex(j)=i
      endif
    1 continue
c
      do 6 iproc=0,np-1
c
      call ga_distribution(ga_cov,iproc,ilc,ihc,jlc,jhc)
c
      do 3 i=ilc,ihc
      ia=ndex((i+2)/3)
      ix=mod(i+2,3)+1
      do 5 j=jlc,jhc
      cov(j-jlc+1)=0.0d0
    5 continue
      do 2 k=1,nfrme
      do 4 j=jlc,jhc
      ja=ndex((j+2)/3)
      jx=mod(j+2,3)+1
      cov(j-jlc+1)=cov(j-jlc+1)+
     + (xtrj(ia,ix,k)-xa(ia,ix))*(xtrj(ja,jx,k)-xa(ja,jx))
    4 continue
    2 continue
      call ga_acc(ga_cov,i,i,jlc,jhc,cov,1,one)
    3 continue
    6 continue
c
      call ga_sync()
c
      return
      end
c $Id$
