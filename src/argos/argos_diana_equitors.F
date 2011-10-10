      subroutine argos_diana_equitors(iequi,imol,idtors)
c
      implicit none
c
#include "argos_diana_common.fh"
#include "argos_diana_params.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "msgids.fh"
#include "util.fh"
c
      integer iequi(msa),imol(msa)
      integer i,j,last,idtors
c
      last=ntors
c
      do 1 i=2,4
      if(imol(itors(ntors,i,1)).ne.imol(itors(ntors,1,1))) return
    1 continue
c
    2 continue
      do 3 i=1,4
      if(iequi(itors(ntors,i,1)).eq.0) goto 9
    3 continue
c
      ntors=ntors+1
      itors(ntors,1,3)=idtors
      if(ntors.gt.mxtors) call md_abort('Increase dimension mxtors',0)
c
      do 4 i=1,4
      itors(ntors,i,1)=iequi(itors(ntors-1,i,1))
      itors(ntors,i,2)=0
    4 continue
      goto 2
c
    9 continue
      return
      end
c $Id$
