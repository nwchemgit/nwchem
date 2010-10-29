      logical function argos_diana_rdfram(x,w)
c
c $Id$
c
      implicit none
c
#include "argos_diana_common.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "util.fh"
c
      real*8 x(msa,3),w(mwm,mwa,3)
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
      if(lrdf) call argos_diana_rdfwrt()
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
      if(lrdf) call argos_diana_rdfhdr(int_mb(i_rdf))
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
      call argos_diana_center_fold(x,w,int_mb(i_imol),nwm)
      endif
      if(nrot.gt.0) then
      call argos_diana_rotate(x,w,nwm)
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
      argos_diana_rdfram=succes
c
      return
      end
