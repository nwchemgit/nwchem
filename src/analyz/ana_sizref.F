      subroutine ana_sizref(nsai,msai)
c
c $Id: ana_sizref.F,v 1.1 1999-01-29 18:10:42 d3j191 Exp $
c
      implicit none
c
#include "ana_common.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c
      character*1 cdum
      integer i,nwm,nwa,nsai,msai
c
      integer ibroad(10)
c
      if(me.eq.0) then
      open(unit=lfnref,file=filref(1:index(filref,' ')-1),
     + status='old',err=9999)
      read(lfnref,1000) cdum
      read(lfnref,1000) cdum
      read(lfnref,1000) cdum
      read(lfnref,1000) cdum
 1000 format(a)
      read(lfnref,1001) nwm
 1001 format(48x,i5)
      read(lfnref,1002) nwa
 1002 format(40x,i5)
      read(lfnref,1000) cdum
      if(nwm.eq.1) then
      read(lfnref,1000) cdum
      read(lfnref,1000) cdum
      read(lfnref,1000) cdum
      endif
      if(nwa.ne.0) then
      read(lfnref,1000) cdum
      read(lfnref,1000) cdum
      read(lfnref,1000) cdum
      endif
      read(lfnref,1003) nwm,nwa,nsai
 1003 format(2i10,10x,i10)
      close(unit=lfnref)
      endif
c
      msai=nsai
c
c     broadcast options to all nodes
c
      if(np.gt.1) then
      ibroad(1)=nsai
      ibroad(2)=msai
      call ga_brdcst(mag_d01,ibroad,2*ma_sizeof(mt_int,1,mt_byte),0)
      nsai=ibroad(1)
      msai=ibroad(2)
      endif
c
      return
c
 9999 continue
      call errquit('Unable to open reference file',0)
      return
      end

