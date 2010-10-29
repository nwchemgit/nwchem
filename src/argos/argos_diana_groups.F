      subroutine argos_diana_groups(sgmnam,imol,isel,wt,x)
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
      integer isel(nsa),imol(msa)
      real*8 x(nsa,3)
c
      integer i
c
      ngroups=ngroups+1
      if(ngroups.gt.maxgrp) call md_abort('Increase maxgrp',maxgrp)
c
      if(card(1:6).eq.'groups') then
      read(card(8:67),1001) (igroups(ngroups,i),i=1,6),
     + (rgroups(ngroups,i),i=1,2)
 1001 format(4i7,i5,i3,2f12.6)
      endif
      print*,'++++++++++',igroups(ngroups,5)
c
      filgrp=card(68:80)
      if(filgrp(1:1).ne.' ') then
      open(unit=lfngrp,file=filgrp(1:index(filgrp,' ')-1),
     + form='formatted',status='unknown')
      if(igroups(ngroups,5).eq.1)
     + call argos_diana_grpcogdis(sgmnam,imol,isel,wt,x,ngroups)
      if(igroups(ngroups,5).eq.2)
     + call argos_diana_grpdis(sgmnam,imol,isel,wt,x,ngroups)
      if(igroups(ngroups,5).eq.4)
     + call argos_diana_grpcogang(sgmnam,imol,isel,wt,x,ngroups)
      if(igroups(ngroups,5).eq.5)
     + call argos_diana_grpvectors(sgmnam,imol,isel,wt,x,ngroups)
      ngroups=ngroups-1
      close(unit=lfngrp)
      endif
c
      return
      end
