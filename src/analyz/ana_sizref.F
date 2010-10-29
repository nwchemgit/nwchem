      subroutine ana_sizref(nsai,msai,nwmi,mwmi,nwai,mwai)
c
c $Id$
c
      implicit none
c
#include "ana_common.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c
      character*1 cdum
      integer i,nsai,msai,nwmi,mwmi,nwai,mwai,nsm
c
      integer ibroad(10),nhist
      logical lforces
c
      if(me.eq.0) then
      open(unit=lfnref,file=filref(1:index(filref,' ')-1),
     + status='old',err=9999)
c
      read(lfnref,1000) cdum
      read(lfnref,1000) cdum
      read(lfnref,1000) cdum
 1000 format(a)
      read(lfnref,1001) nhist,lforces
 1001 format(32x,i5,4x,l1)
      if(nhist.gt.0) then
      do 1 i=1,nhist
      read(lfnref,1000) cdum
    1 continue
      endif
c
      read(lfnref,1000) cdum
c
      read(lfnref,1000) cdum
      read(lfnref,1000) cdum
      read(lfnref,1000) cdum
c
      read(lfnref,1000) cdum
      read(lfnref,1000) cdum
c
      read(lfnref,1002) nwm,nwa,nsm,nsai
 1002 format(4i10)
c
      close(unit=lfnref)
      endif
c
      msai=nsai
      nwmi=nwm
      mwmi=nwm
      nwai=nwa
      mwai=nwa
c
c     broadcast options to all nodes
c
      if(np.gt.1) then
      ibroad(1)=nsai
      ibroad(2)=msai
      ibroad(3)=nwmi
      ibroad(4)=mwmi
      ibroad(5)=nwai
      ibroad(6)=mwai
      call ga_brdcst(mag_d01,ibroad,6*ma_sizeof(mt_int,1,mt_byte),0)
      nsai=ibroad(1)
      msai=ibroad(2)
      nwmi=ibroad(3)
      mwmi=ibroad(4)
      nwai=ibroad(5)
      mwai=ibroad(6)
      endif
c
      return
c
 9999 continue
      call md_abort('Unable to open reference file',0)
      return
      end

