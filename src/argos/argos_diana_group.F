      subroutine argos_diana_group(sgmnam,imol,isel,wt,x,iwrk)
c
c $Id$
c
      implicit none
c
#include "argos_diana_common.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "rtdb.fh"
c
      character*16 sgmnam(nsa)
      real*8 wt(nsa)
      integer isel(nsa),imol(msa),iwrk(mxdef,mxnum,maxgrp)
      real*8 x(nsa,3)
c
      integer i,j
c
      ngroup=ngroup+1
      if(ngroup.gt.maxgrp) call md_abort('Increase maxgrp',maxgrp)
c
      read(card(8:46),1000) (igroup(ngroup,i),i=1,3),
     + (rgroup(ngroup,i),i=1,2)
 1000 format(i7,i5,i3,2f12.6)
c
      do 1 i=1,mxdef
      do 2 j=1,mxnum
      iwrk(i,j,ngroup)=0
    2 continue
    1 continue
c
      return
      end
