      subroutine argos_cafe_pmec(ipset,xw,nwmloc,xs,isq,isgm,nsaloc,
     + ecorr,iwq,chg,lseq)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "mafdecls.fh"
c
      integer ipset,nwmloc,nsaloc
      real*8 xw(mwm,3,mwa),xs(msa,3),ecorr
      integer isq(msa),iwq(mwa),isgm(msa),lseq(mseq)
      real*8 chg(mqt,mqp,mset)
c
      integer i,icset
c
      ecorr=zero
c
      if(nwmloc.gt.0) then
      do 1 i=1,nwa
      ecorr=ecorr-chg(iwq(i),1,ipset)*chg(iwq(i),1,ipset)
    1 continue
      ecorr=ecorr*dble(nwmloc)
      endif
c
      if(nfhop.eq.0) then
      if(nsaloc.gt.0) then
      do 2 i=1,nsaloc
      ecorr=ecorr-chg(isq(i),1,ipset)*chg(isq(i),1,ipset)
    2 continue
      endif
      else
      do 3 i=1,nsaloc
      icset=lseq(isgm(i))
      ecorr=ecorr-chg(isq(i),1,icset)*chg(isq(i),1,icset)
    3 continue
      endif
c
      ecorr=ealpha*ecorr/sqrt(pi)
c
      return
      end
c $Id$
