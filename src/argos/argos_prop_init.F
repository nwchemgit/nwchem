      subroutine argos_prop_init()
c
      implicit none
c
#include "argos_prop.fh"
c
      integer i
c
      nsum=0
      nsumt=0
      nsump=0
      do 1 i=1,maxpro
      psum(i)=zero
      p2sum(i)=zero
      pslop(i)=zero
      psumt(i)=zero
      p2sumt(i)=zero
      pslopt(i)=zero
      psump(i)=zero
    1 continue
      tsum=zero
      t2sum=zero
      tsumt=zero
      t2sumt=zero
c
      nderiv=0
      dfree=zero
      do 2 i=1,24
      deriv(i)=zero
    2 continue
c
      lhdr=.false.
      lfhdr=.false.
      lauto=.false.
      lappr=.false.
      samrat=one
c
      return
      end
c $Id$
