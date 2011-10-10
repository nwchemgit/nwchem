      subroutine argos_diana_lochdr(sgmnam,imol,isel)
c
      implicit none
c
#include "argos_diana_common.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c
      character*16 sgmnam(nsa)
      integer imol(msa),isel(nsa)
c
      integer i,num,igrp,ia,j,ito
c
      num=0
      do 1 i=1,nsa
      if(isel(i).ne.0) num=num+1
    1 continue
c
      write(lfnloc,1000) ngroup
 1000 format(2i7)
      do 4 j=1,ngroup
      igrp=igroup(j,1)
      ito=ldef(igrp)
      write(lfnloc,1000) ito,nsa
      do 2 i=1,ito
      ia=idef(igrp,i)
      write(lfnloc,1001) ia,imol(ia),sgmnam(ia)(11:16),sgmnam(ia)(1:5),
     + sgmnam(ia)(6:10)
    2 continue
    4 continue
      write(lfnloc,1000) num,nsa
      do 3 i=1,nsa
      if(isel(i).ne.0) then
      write(lfnloc,1001) i,imol(i),sgmnam(i)(11:16),sgmnam(i)(1:5),
     + sgmnam(i)(6:10)
 1001 format(i7,i5,a6,' ',a5,':',a5)
      endif
    3 continue
c
      return
      end
c $Id$
