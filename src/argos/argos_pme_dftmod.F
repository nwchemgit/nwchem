      subroutine argos_pme_dftmod(bmod,barr,ng)
c
      implicit none
c
#include "argos_pme_common.fh"
c
      integer ng
      real*8 bmod(ng),barr(ng)
      integer i,j
      real*8 sum1,sum2,arg
c
      do 1 i=1,ng
      sum1=zero
      sum2=zero
      do 2 j=1,ng
      arg=twopi*dble((i-1)*(j-1))/dble(ng)
      sum1=sum1+barr(j)*cos(arg)
      sum2=sum2+barr(j)*sin(arg)
    2 continue
      bmod(i)=sum1*sum1+sum2*sum2
    1 continue
      do 3 i=1,ng
      if(bmod(i).lt.small) bmod(i)=half*(bmod(i-1)+bmod(i+1))
    3 continue
c
      return
      end
c $Id$
