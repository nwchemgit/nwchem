      subroutine argos_space_cascad(ltran,lbbl)
c
      implicit none
c
#include "argos_space_common.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c
      integer argos_space_btop
      external argos_space_btop
c
      integer ltran(np,3),lbbl(mbbl,mbb2)
c
      integer i,j,k,ltemp
      integer icell,jcell,iproc,jproc
      logical lfirst
c
c     load balancing based on cascade box-list redistribution
c     -------------------------------------------------------
c
      do 1 i=1,np
      ltran(i,1)=ltran(i,1)-1
      ltran(i,2)=ltran(i,2)-1
      ltran(i,3)=-1
      if(ltran(i,2).ge.0)
     + ltran(i,3)=argos_space_btop(ltran(i,2),int_mb(i_iown))
    1 continue
c
c     remove transfers involving double processors
c
      do 2 i=2,np
      if(ltran(i,1).ge.0) then
      iproc=argos_space_btop(ltran(i,1),int_mb(i_iown))
      jproc=argos_space_btop(ltran(i,2),int_mb(i_iown))
      do 3 j=1,i-1
      if(argos_space_btop(ltran(j,1),int_mb(i_iown)).eq.jproc.or.
     + argos_space_btop(ltran(j,2),int_mb(i_iown)).eq.jproc) then
      ltran(i,1)=-1
      ltran(i,2)=-1
      goto 2
      endif
    3 continue
      endif
    2 continue
c
c     process the transfer list
c
      lfirst=.true.
      do 4 i=1,np
c
      iproc=-1
      jproc=-1
      icell=ltran(i,1)
      jcell=ltran(i,2)
      if(icell.ge.0) then
      iproc=argos_space_btop(icell,int_mb(i_iown))
      jproc=argos_space_btop(jcell,int_mb(i_iown))
      if(lfirst) then
      lpipo=(iproc.eq.ipairt).and.(jproc.eq.ipairf)
      ipairf=iproc
      ipairt=jproc
      lfirst=.false.
      endif
c
c     remove pair
c
      if(iproc.eq.me) then
      k=0
      do 5 j=1,nbbl
      if(lbbl(j,2).ne.jcell.or.lbbl(j,3).ne.icell) then
      k=k+1
      lbbl(k,1)=lbbl(j,1)
      lbbl(k,2)=lbbl(j,2)
      lbbl(k,3)=lbbl(j,3)
      endif
    5 continue
      nbbl=k
      endif
c
c     add pair
c
      if(jproc.eq.me) then
      nbbl=nbbl+1
      lbbl(nbbl,1)=iproc
      lbbl(nbbl,2)=icell
      lbbl(nbbl,3)=jcell
      endif
c
      endif
    4 continue
c
c     order the new box-box list
c
      nbbloc=0
      do 6 i=1,nbbl-1
      do 7 j=i+1,nbbl
      if((lbbl(i,1).ne.me.and.lbbl(j,1).eq.me).or.
     + (lbbl(i,1).gt.lbbl(j,1).and.lbbl(i,1).ne.me).or.
     + (lbbl(i,1).eq.lbbl(j,1).and.lbbl(i,2).gt.lbbl(j,2)).or.
     + (lbbl(i,1).eq.lbbl(j,1).and.lbbl(i,2).eq.lbbl(j,2).and.
     + lbbl(i,3).gt.lbbl(j,3))) then
      ltemp=lbbl(i,1)
      lbbl(i,1)=lbbl(j,1)
      lbbl(j,1)=ltemp
      ltemp=lbbl(i,2)
      lbbl(i,2)=lbbl(j,2)
      lbbl(j,2)=ltemp
      ltemp=lbbl(i,3)
      lbbl(i,3)=lbbl(j,3)
      lbbl(j,3)=ltemp
      endif
    7 continue
      if(lbbl(i,1).eq.me) nbbloc=i
    6 continue
      if(lbbl(nbbl,1).eq.me) nbbloc=nbbl
c
      return
      end
c $Id$
