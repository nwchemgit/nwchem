      subroutine argos_cafe_lww(lself,lpbc,xwm,idt,iwfr,iwto,jwfr,jwto,
     + nlocw,mpairs,npairs,lwwjpt,lwwin,lwwj,list,rwx,rw)
c $Id$
      implicit none
c
#include "argos_cafe_common.fh"
#include "bitops.fh"
#include "mafdecls.fh"
c
      integer nlocw
      real*8 rwx(mscr,3),rw(mscr)
      real*8 xwm(mwm,3)
      integer lwwj(*),lwwjpt(nlocw,2)
      integer lwwin(nlocw,2)
      integer list(mwm,2),idt(mwm)
      logical lself,lpbc
c
      integer iwm,jwm,iwmfr,iwfr,iwmto,iwto,jwmfr,jwfr,jwmto,jwto
      integer ilist,nlist,mpairs,npairs
c
      iwmfr=iwfr
      if(lself) then
      lwwin(nlocw,1)=0
      lwwin(nlocw,2)=0
      lwwjpt(nlocw,1)=0
      lwwjpt(nlocw,2)=0
      iwmto=iwto-1
      jwmfr=0
      jwmto=iwto
      else
      iwmto=iwto
      jwmfr=jwfr
      jwmto=jwto
      endif
c
      npairs=0
      do 1 iwm=iwmfr,iwmto
      if(lself) jwmfr=iwm+1
c
      do 3 jwm=jwmfr,jwmto
      rwx(jwm,1)=xwm(iwm,1)-xwm(jwm,1)
      rwx(jwm,2)=xwm(iwm,2)-xwm(jwm,2)
      rwx(jwm,3)=xwm(iwm,3)-xwm(jwm,3)
    3 continue
      if(lpbc) call argos_cafe_pbc(1,rwx,mscr,rwx,mscr,0,jwmfr,jwmto)
      do 8 jwm=jwmfr,jwmto
      rw(jwm)=rwx(jwm,1)**2+rwx(jwm,2)**2+rwx(jwm,3)**2
    8 continue
      do 4 jwm=jwmfr,jwmto
      list(jwm,1)=0
      list(jwm,2)=0
      if(rw(jwm).lt.rshrt2.and.
     + (iand(idt(iwm),mdynam).eq.ldynam.or.
     + iand(idt(jwm),mdynam).eq.ldynam)) list(jwm,1)=1
      if(rw(jwm).lt.rrest2 .and.
     + (iand(idt(iwm),mrestr).eq.lrestr.or.
     + iand(idt(jwm),mrestr).eq.lrestr)) list(jwm,1)=1
    4 continue
      if(npww.eq.2) then
      do 2 jwm=jwmfr,jwmto
      if(rw(jwm).lt.rlong2.and.
     + (iand(idt(iwm),mdynam).eq.ldynam.or.
     + iand(idt(jwm),mdynam).eq.ldynam)) list(jwm,2)=1
      if(list(jwm,1).eq.1) list(jwm,2)=0
    2 continue
      endif
c
c     short range pairlist
c
      nlist=0
      do 9 jwm=jwmfr,jwmto
      if(list(jwm,1).eq.1) then
      nlist=nlist+1
      list(nlist,1)=jwm
      endif
    9 continue
      if(npairs+nlist.gt.mpairs)
     + call md_abort('Insufficient memory for pairlist',0)
      do 5 ilist=1,nlist
      lwwj(npairs+ilist)=list(ilist,1)
    5 continue
      lwwjpt(iwm-iwmfr+1,1)=npairs+1
      npairs=npairs+nlist
      lsww=lsww+nlist
      lwwin(iwm-iwmfr+1,1)=nlist
c
c     long range pairs
c
      nlist=0
      if(npww.eq.2) then
      do 12 jwm=jwmfr,jwmto
      if(list(jwm,2).eq.1) then
      nlist=nlist+1
      list(nlist,2)=jwm
      endif
   12 continue
      if(npairs+nlist.gt.mpairs)
     + call md_abort('Insufficient memory for pairlist',0)
      do 13 ilist=1,nlist
      lwwj(npairs+ilist)=list(ilist,2)
   13 continue
      lwwjpt(iwm-iwmfr+1,2)=npairs+1
      npairs=npairs+nlist
      llww=llww+nlist
      lwwin(iwm-iwmfr+1,2)=nlist
      endif
    1 continue
c
      return
      end
