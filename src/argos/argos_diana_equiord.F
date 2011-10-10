      subroutine argos_diana_equiord(iequi,imol,idordr,ndordr)
c
      implicit none
c
#include "argos_diana_common.fh"
#include "argos_diana_params.fh"
c
      integer iequi(msa),imol(msa)
      integer i,j,last,idordr,ndordr
c
      last=nord
c
      do 1 i=2,ndordr+1
      if(imol(iord(nord,i,1)).ne.imol(iord(nord,1,1))) return
    1 continue
c
    2 continue
      do 3 i=1,ndordr+1
      if(iequi(iord(nord,i,1)).eq.0) goto 9
    3 continue
c
      nord=nord+1
      iord(nord,1,3)=idordr
      iord(nord,3,3)=ndordr
      if(nord.gt.mord) call md_abort('Increase dimension mord',0)
c
      do 4 i=1,ndordr+1
      iord(nord,i,1)=iequi(iord(nord-1,i,1))
      iord(nord,i,2)=0
    4 continue
      goto 2
c
    9 continue
c
      return
      end
c $Id$
