      logical function dia_rdfram(x,w,t)
c
c $Id$
c
      implicit none
c
#include "dia_common.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "util.fh"
c
      real*8 x(msa,3),w(mwm,mwa,3),t(6)
      integer ftell
      external ftell
c
      integer i,j,k,lq
      character*255 fname
c
      logical succes
      real*8 timp
      integer inwm,inwa
      integer flags(7)
c
      lxw=.false.
      lvw=.false.
      lfw=.false.
      lsx=.false.
      lsv=.false.
      lfs=.false.
c
      if(me.eq.0) then
  100 continue
      if(fmttrj.eq.'trj') then
    1 continue
      read(lfntrj,1000,err=9999,end=9999) card
 1000 format(a)
      if(card(1:5).ne.'frame') goto 1
      read(lfntrj,1001,err=9999,end=9999) timp,temp,pres,datum,tijd
 1001 format(2f12.6,e12.5,2a10)
      t(1)=timp
      t(2)=temp
      t(3)=pres
      if(lclean.and.timp.lt.timr) then
      write(*,2221) timp,temp,pres,datum,tijd
 2221 format(' skip  ',2f12.6,e12.5,2a10)
      goto 1
      endif
      if(util_print('frames_read',print_debug))
     + write(*,2201) timp,temp,pres,datum,tijd
 2201 format(' frame ',2f12.6,e12.5,2a10)
      if(timp.lt.timr) timoff=time
      time=timp+timoff
      timr=timp
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
      if(util_print('frames_read',print_debug)) then
      write(*,'(a,f12.6,3i10,1x,4l1)')
     + ' Frame read:',time,inwm,inwa,nsa,lxw,lvw,lsx,lsv
      endif
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
      succes=.true.
      goto 9000
      elseif(fmttrj.eq.'sco'.or.fmttrj.eq.'coo') then
      datum='          '
      tijd='          '
   11 continue
      read(lfntrj,1000,err=9999,end=9999) card
      if(card(1:5).ne.'time:') goto 11
      read(lfntrj,2001) timp
 2001 format(f12.6)
      if(timp.lt.timr) timoff=time
      time=timp+timoff
      timr=timp
   12 continue
      read(lfntrj,1000,err=9999,end=9999) card
      if(card(1:15).ne.'box_dimensions:') goto 12
      read(lfntrj,2002) box
 2002 format(3f12.6)
      lxw=fmttrj.eq.'coo'
      lvw=.false.
      lsx=.true.
      lsv=.false.
   13 continue
      read(lfntrj,1000,err=9999,end=9999) card
      if(card(1:10).ne.'atom_list:') goto 13
      do 14 i=1,nsa
      read(lfntrj,1005,err=9999,end=9999) (x(i,j),j=1,3)
      x(i,1)=x(i,1)*1.0d-1
      x(i,2)=x(i,2)*1.0d-1
      x(i,3)=x(i,3)*1.0d-1
   14 continue
      if(lxw) then
      do 15 i=1,nwm
      do 16 j=1,nwa
      read(lfntrj,1005,err=9999,end=9999) (w(i,j,k),k=1,3)
   16 continue
   15 continue
      endif
      succes=.true.
      goto 9000
      else
      call md_abort('Trajectory file format error',0)
      endif
c
 9999 continue 
c
      close(unit=lfntrj)
      write(*,'(a)') ' Closing trj file '
      if(lrdf) call dia_rdfwrt()
c
      fname=filtrj
c
      lq=index(filtrj,'?')
      if(lq.eq.0) goto 9998
c
      iscof=iscof+1
      if(iscof.gt.ilast) goto 9998
      write(cnum,'(i3.3)') iscof
c
      fname=filtrj(1:lq-1)//cnum//filtrj(lq+1:index(filtrj,' ')-1)
      open(unit=lfntrj,file=fname(1:index(fname,' ')-1),
     + status='old',err=9998)
      write(*,3333) fname(1:index(fname,' ')-1)
 3333 format(/,' Opening trj file ',a)
c
      if(lrdf) call dia_rdfhdr(int_mb(i_rdf))
      goto 100
c
 9998 continue
      succes=.false.
c
      endif
 9000 continue
c
c     center and fold and rotate
c     --------------------------
c
      if(me.eq.0) then
      if(icent.gt.0) then
      if(.not.ltop) call md_abort('No topology for center',0)
      call dia_center_fold(x,w,int_mb(i_imol),nwm)
      endif
      if(nrot.gt.0) then
      call dia_rotate(x,w,nwm)
      endif
      endif
c
c     broadcast to all nodes
c
      if(np.gt.1) then
      call ga_brdcst(mag_d03,succes,ma_sizeof(mt_log,1,mt_byte),0)
      call ga_brdcst(mag_d05,x,3*msa*ma_sizeof(mt_dbl,1,mt_byte),0)
      call ga_brdcst(mag_d06,w,3*mwm*mwa*ma_sizeof(mt_dbl,1,mt_byte),0)
      do 234 i=1,6
      flags(i)=0
  234 continue
      if(lxw) flags(1)=1
      if(lvw) flags(2)=1
      if(lfw) flags(3)=1
      if(lsx) flags(4)=1
      if(lsv) flags(5)=1
      if(lfs) flags(6)=1
      flags(7)=nwm
      call ga_brdcst(mag_d07,flags,7*ma_sizeof(mt_int,1,mt_byte),0)
      lxw=flags(1).ne.0
      lvw=flags(2).ne.0
      lfw=flags(3).ne.0
      lsx=flags(4).ne.0
      lsv=flags(5).ne.0
      lfs=flags(6).ne.0
      nwm=flags(7)
      endif
c
      dia_rdfram=succes
c
      return
      end
      subroutine dia_center_fold(xs,xw,imol,nw)
c
      implicit none
c
#include "dia_common.fh"
#include "dia_params.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "util.fh"
c
      real*8 xs(msa,3),xw(mwm,mwa,3)
      integer imol(msa),nw
c
      integer i,j,k,nmol,ncount
      real*8 xc(3),dx,factor
c
      real*8 boxh(3)
c
      boxh(1)=half*box(1)
      boxh(2)=half*box(2)
      boxh(3)=half*box(3)
c
      xc(1)=zero
      xc(2)=zero
      xc(3)=zero
      ncount=0
c
      nmol=0
      do 1 i=1,nsa
      nmol=max(nmol,imol(i))
      if(imol(i).ge.icent.and.imol(i).le.jcent) then
      ncount=ncount+1
      xc(1)=xc(1)+xs(i,1)
      xc(2)=xc(2)+xs(i,2)
      xc(3)=xc(3)+xs(i,3)
      endif
    1 continue
      if(ncount.eq.0) call md_abort('No atom in molecule',icent)
      xc(1)=xc(1)/dble(ncount)
      xc(2)=xc(2)/dble(ncount)
      xc(3)=xc(3)/dble(ncount)
c
      do 2 i=1,nsa
      xs(i,1)=xs(i,1)-xc(1)
      xs(i,2)=xs(i,2)-xc(2)
      xs(i,3)=xs(i,3)-xc(3)
    2 continue
      do 3 j=1,nwa
      do 4 i=1,nw
      xw(i,j,1)=xw(i,j,1)-xc(1)
      xw(i,j,2)=xw(i,j,2)-xc(2)
      xw(i,j,3)=xw(i,j,3)-xc(3)
    4 continue
    3 continue
c
      do 5 j=1,nmol
      xc(1)=zero
      xc(2)=zero
      xc(3)=zero
      ncount=0
      do 6 i=1,nsa
      nmol=max(nmol,imol(i))
      if(imol(i).eq.j) then
      ncount=ncount+1
      xc(1)=xc(1)+xs(i,1)
      xc(2)=xc(2)+xs(i,2)
      xc(3)=xc(3)+xs(i,3)
      endif
    6 continue
      if(ncount.eq.0) goto 5
      xc(1)=xc(1)/dble(ncount)
      xc(2)=xc(2)/dble(ncount)
      xc(3)=xc(3)/dble(ncount)
      do 7 k=1,3
      dx=zero
      if(xc(k).lt.-boxh(k)) dx=box(k)
      if(xc(k).gt.boxh(k)) dx=-box(k)
      do 8 i=1,nsa
      if(imol(i).eq.j) xs(i,k)=xs(i,k)+dx
    8 continue
    7 continue
    5 continue
      if(nwa.gt.0.and.nw.gt.0) then
      factor=one/dble(nwa)
      do 9 i=1,nw
      xc(1)=zero
      xc(2)=zero
      xc(3)=zero
      do 10 j=1,nwa
      xc(1)=xc(1)+xw(i,j,1)
      xc(2)=xc(2)+xw(i,j,2)
      xc(3)=xc(3)+xw(i,j,3)
   10 continue
      xc(1)=xc(1)*factor
      xc(2)=xc(2)*factor
      xc(3)=xc(3)*factor
      do 11 k=1,3
      dx=zero
      if(xc(k).lt.-boxh(k)) dx=box(k)
      if(xc(k).gt.boxh(k)) dx=-box(k)
      do 12 j=1,nwa
      xw(i,j,k)=xw(i,j,k)+dx
   12 continue
   11 continue
    9 continue
      endif
c
      return
      end
      subroutine dia_rotate(xs,xw,nw)
c
      implicit none
c
#include "dia_common.fh"
#include "dia_params.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "util.fh"
c
      real*8 xs(msa,3),xw(mwm,mwa,3)
      integer nw
c
      real*8 rv(3),rw(3),rotang,rx(3),ry(3),pifac
      integer i,j,k
c
      pifac=tan(1.0d0)/90.0d0
c
      do 1 k=1,nrot
      rv(1)=0.0d0
      rv(2)=0.0d0
      rv(3)=0.0d0
      rw(1)=0.0d0
      rw(2)=0.0d0
      rw(3)=0.0d0
      rw(irot(k))=1.0d0
      rotang=arot(k)*pifac
      do 2 i=1,nsa
      rx(1)=xs(i,1)
      rx(2)=xs(i,2)
      rx(3)=xs(i,3)
      ry(1)=xs(i,1)
      ry(2)=xs(i,2)
      ry(3)=xs(i,3)
      call rotate(rv,rw,rotang,ry,rx)
      xs(i,1)=rx(1)
      xs(i,2)=rx(2)
      xs(i,3)=rx(3)
    2 continue
      do 3 i=1,nw
      do 4 j=1,nwa
      rx(1)=xw(i,j,1)
      rx(2)=xw(i,j,2)
      rx(3)=xw(i,j,3)
      ry(1)=xw(i,j,1)
      ry(2)=xw(i,j,2)
      ry(3)=xw(i,j,3)
      call rotate(rv,rw,rotang,ry,rx)
      xw(i,j,1)=rx(1)
      xw(i,j,2)=rx(2)
      xw(i,j,3)=rx(3)
    4 continue
    3 continue
    1 continue
c
      return
      end
      subroutine dia_preadframe(x,w,t)
c
      implicit none
c
#include "dia_common.fh"
#include "dia_params.fh"
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
      call md_abort('EOF in dia_preadfram',me)
 9999 continue
      call md_abort('Error in dia_preadfram',me)
      return
      end
