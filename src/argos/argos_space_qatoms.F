      subroutine argos_space_qatoms(isl,numsa)
c
      implicit none
c
#include "argos_space_common.fh"
c
      integer argos_cafe_quantuma
      external argos_cafe_quantuma
c
      integer numsa
      integer isl(msa,mis2)
c
      integer i,itype
c
      do 1 i=1,numsa
      itype=argos_cafe_quantuma(isl(i,lsatt))
      if(itype.eq.1.or.itype.eq.2) isl(i,lsdyn)=ior(isl(i,lsdyn),lquant)
      if(itype.eq.2) isl(i,lsdyn)=ior(isl(i,lsdyn),lqhigh)
    1 continue
c
      return
      end
c $Id$
