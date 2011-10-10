      subroutine argos_cafe_cntlnk(ltype,msb,ibnd,nums,isgan,isdt,lst)
c
      implicit none
c
#include "argos_cafe_common.fh"
c
      integer msb,nums,ltype
      integer ibnd(msb,3),isgan(msa),isdt(msa),lst(nsa)
c
      integer i,j
c
      do 1 j=1,nums
      if(iand(isdt(j),ltype).ne.0) then
      do 2 i=1,numb(2)
      if(ibnd(i,1).eq.isgan(j)) lst(ibnd(i,2))=lst(ibnd(i,2))+1
      if(ibnd(i,2).eq.isgan(j)) lst(ibnd(i,1))=lst(ibnd(i,1))+1
    2 continue
      endif
    1 continue
c
      return
      end
c $Id$
