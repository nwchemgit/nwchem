      logical function ana_rdfram(x,ix,w)
c
c $Id: ana_rdfram.F,v 1.11 2000-05-16 14:46:45 d3j191 Exp $
c
      implicit none
c
#include "ana_common.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c
      real*8 x(msa,3),w(mwm,mwa,3)
      integer ix(msa)
c
      character*80 card
      integer i,j,k,lq
      character*255 fname
      character*3 cnum
c
      logical succes
      real*8 timp
c
      lxw=.false.
      lvw=.false.
      lsx=.false.
      lsv=.false.
c
      if(me.eq.0) then
  100 continue
      if(fmttrj.eq.'trj') then
    1 continue
      read(lfntrj,1000,err=9999,end=9999) card
 1000 format(a)
      if(card(1:5).ne.'frame') goto 1
      read(lfntrj,1001,err=9999,end=9999) timp,temp,pres
 1001 format(2f12.6,e12.5)
      if(timp.lt.timr) timoff=time
      time=timp+timoff
      timr=timp
      read(lfntrj,1002) box
 1002 format(f12.6,36x,f12.6,36x,f12.6)
      read(lfntrj,1003) lxw,lvw,lsx,lsv,nwm,nwa,nsa
 1003 format(4l1,3i10)
      if(nsa.gt.msa) call errquit('Error in trj 1',0)
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
      ix(i)=i
    4 continue
      endif
      succes=.true.
      goto 9000
      elseif(fmttrj.eq.'sco') then
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
      lxw=.false.
      lvw=.false.
      lsx=.true.
      lsv=.false.
   13 continue
      read(lfntrj,1000,err=9999,end=9999) card
      if(card(1:10).ne.'atom_list:') goto 13
      do 14 i=1,nsa
      read(lfntrj,1005) (x(i,j),j=1,3)
      x(i,1)=x(i,1)*1.0d-1
      x(i,2)=x(i,2)*1.0d-1
      x(i,3)=x(i,3)*1.0d-1
      ix(i)=i
   14 continue
      succes=.true.
      goto 9000
      else
      call errquit('Trajectory file format error',0)
      endif
c
 9999 continue 
c
      close(unit=lfntrj)
      write(*,'(a)') ' Closing trj file '
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
      goto 100
c
 9998 continue
      succes=.false.
c
      endif
 9000 continue
c
c     broadcast to all nodes
c
      if(np.gt.1) then
      call ga_brdcst(mag_d03,succes,ma_sizeof(mt_log,1,mt_byte),0)
      call ga_brdcst(mag_d04,ix,msa*ma_sizeof(mt_int,1,mt_byte),0)
      call ga_brdcst(mag_d05,x,3*msa*ma_sizeof(mt_dbl,1,mt_byte),0)
      endif
c
      ana_rdfram=succes
c
      return
      end
