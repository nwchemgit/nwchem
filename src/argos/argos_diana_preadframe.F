      subroutine argos_diana_preadframe(x,w,t)
c
      implicit none
c
#include "argos_diana_common.fh"
#include "argos_diana_params.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "util.fh"
c
      integer ftell
      external ftell
      real*8 x(msa,3),w(mwm,mwa,3),t(6)
      integer i,j,k
      real*8 timp
      integer inwm,inwa
c
      read(lfntrj,1000,err=9999,end=9998) card
 1000 format(a)
      if(card(1:5).ne.'frame')
     + call md_abort('Could not find frame',0)
      read(lfntrj,1001,err=9999,end=9998) timp,temp,pres,datum,tijd
 1001 format(2f12.6,e12.5,2a10)
      t(1)=timp
      t(2)=temp
      t(3)=pres
      read(lfntrj,1000) card
      if(card(54:54).eq.'.') then
      read(card,1012) box
 1012 format(f12.6,36x,f12.6,36x,f12.6)
      else
      read(card,1002) box(1)
 1002 format(f12.6)
      read(lfntrj,1022) box(2),box(3)
 1022 format(12x,f12.6,/,24x,f12.6)
      endif
      t(4)=box(1)
      t(5)=box(2)
      t(6)=box(3)
      read(lfntrj,1000) card
      if(card(8:8).eq.'F'.or.card(8:8).eq.'T') then
      read(card,1006) lxw,lvw,lfw,lpw,lsx,lsv,lfs,lps,inwm,inwa,nsa
      elseif(card(6:6).eq.'F'.or.card(6:6).eq.'T') then
      read(card,1007) lxw,lvw,lfw,lsx,lsv,lfs,inwm,inwa,nsa
      lpw=.false.
      lps=.false.
      else
      read(card,1003) lxw,lvw,lsx,lsv,inwm,inwa,nsa
      lfw=.false.
      lfs=.false.
      endif
      if(inwm.gt.0) then
      nwm=inwm
      nwa=inwa
      endif
      if(nwm.gt.mwm)
     + call md_abort('Restart file not compatible with trajectory',nwm)
 1003 format(4l1,3i10,2l1)
 1006 format(8l1,3i10)
 1007 format(6l1,3i10)
      if(nsa.gt.msa) call md_abort('Error in trj 1',0)
      if((lxw.or.lvw).and.nwm.gt.0) then
      do 2 i=1,nwm
      do 3 j=1,nwa
      read(lfntrj,1004) (w(i,j,k),k=1,3)
 1004 format(3f8.3)
    3 continue
    2 continue
      endif
      if(lsx.and.nsa.gt.0) then
      do 4 i=1,nsa
      read(lfntrj,1005) (x(i,j),j=1,3)
 1005 format(3f8.3)
    4 continue
      endif
c
      return
 9998 continue
      call md_abort('EOF in argos_diana_preadfram',me)
 9999 continue
      call md_abort('Error in argos_diana_preadfram',me)
      return
      end
c $Id$
