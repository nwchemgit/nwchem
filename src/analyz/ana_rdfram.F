      logical function ana_rdfram(x,ix)
c
c $Id: ana_rdfram.F,v 1.1 1999-01-29 18:10:41 d3j191 Exp $
c
      implicit none
c
#include "ana_common.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c
      real*8 x(msa,3)
      integer ix(msa)
c
      character*80 card
      integer i,j,na,lq
      character*255 fname
      character*3 cnum
c
      logical succes
c
      if(me.eq.0) then
    1 continue
      read(lfnsco,1000,err=9999,end=9999) card
 1000 format(a)
      if(card(1:5).ne.'time:') goto 1
      read(lfnsco,1001,err=9999,end=9999) time
 1001 format(f12.6)
      time=time+timoff
    2 continue
      read(lfnsco,1000,err=9999,end=9999) card
      if(card(1:15).ne.'box_dimensions:') goto 2
      read(lfnsco,1002,err=9999,end=9999) box
 1002 format(3f12.6)
    3 continue
      read(lfnsco,1000,err=9999,end=9999) card
      if(card(1:10).ne.'atom_list:') goto 3
      do 4 i=1,nsa
      read(lfnsco,1003,err=9999,end=9999) (x(i,j),j=1,3)
 1003 format(3f8.3)
      x(i,1)=0.1d0*x(i,1)
      x(i,2)=0.1d0*x(i,2)
      x(i,3)=0.1d0*x(i,3)
      ix(i)=i
    4 continue
      succes=.true.
      goto 9000
c
 9999 continue 
c
      close(unit=lfnsco)
c
      fname=filsco
c
      lq=index(filsco,'?')
      if(lq.eq.0) goto 9998
c
      iscof=iscof+1
      if(iscof.gt.ilast) goto 9998
      write(cnum,'(i3.3)') iscof
c
      fname=filsco(1:lq-1)//cnum//filsco(lq+1:index(filsco,' ')-1)
      open(unit=lfnsco,file=fname(1:index(fname,' ')-1),
     + status='old',err=9998)
      write(*,'(a,a)') 'Opening sco file ',fname(1:index(fname,' ')-1)
      goto 1
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
