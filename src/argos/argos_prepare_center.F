      logical function argos_prepare_center(lfnout,xw,mwm,mwa,nwm,nwa,
     + xwc,mwmc,nwmc,xs,msa,nsa,iopt)
c
c $Id$
c
      implicit none
c
#include "util.fh"
c
      integer mwa,mwm,nwa,nwm,msa,nsa,mwmc,nwmc,iopt,lfnout
      real*8 xw(3,mwa,mwm),xwc(3,mwa,mwmc),xs(3,msa)
c
      integer i,j,k,imax,imin,kfr,kto
      real*8 s,smin,smax
c
      if(nsa.gt.0) then
c
      kfr=1
      kto=3
      if(iopt.ne.0) then
      kfr=iopt
      kto=iopt
      endif
c
      do 1 k=kfr,kto
      smin=xs(k,1)
      smax=xs(k,1)
      imin=1
      imax=1
      do 2 i=2,nsa
      if(xs(k,i).lt.smin) imin=i
      if(xs(k,i).gt.smax) imax=i
      smin=min(smin,xs(k,i))
      smax=max(smax,xs(k,i))
    2 continue
      s=(smax+smin)/2.0d0
      do 3 i=1,nsa
      xs(k,i)=xs(k,i)-s
    3 continue
      if(nwmc.gt.0) then
      do 4 i=1,nwmc
      do 5 j=1,nwa
      xwc(k,j,i)=xwc(k,j,i)-s
    5 continue
    4 continue
      endif
      if(nwm.gt.0) then
      do 6 i=1,nwm
      do 7 j=1,nwa
      xw(k,j,i)=xw(k,j,i)-s
    7 continue
    6 continue
      endif
      if(util_print('restart',print_high)) then
      write(lfnout,'(i2,2i5,3f12.6)') k,imin,imax,smin,smax,s
      endif
      if(util_print('restart',print_low)) then
      if(k.eq.0) write(lfnout,1000)
 1000 format(' Solute centered')
      if(k.eq.1) write(lfnout,1001)
 1001 format(' Solute centered in x-dimension')
      if(k.eq.2) write(lfnout,1002)
 1002 format(' Solute centered in y-dimension')
      if(k.eq.3) write(lfnout,1003)
 1003 format(' Solute centered in z-dimension')
      endif
    1 continue
c
      endif
c
      argos_prepare_center=.true.
      return
      end
