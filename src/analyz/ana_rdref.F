      subroutine ana_rdref(x)
c
c $Id$
c
      implicit none
c
#include "ana_common.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "rtdb.fh"
c
      real*8 x(nsa,3)
      character*1 cdum
      integer i,j
      real*8 rdum
      integer nbxt,npbt,nsm
      integer mdold,nhist
      logical lforces
c
      if(.not.ma_verify_allocator_stuff()) print*,'Oopsi 1'
c
      if(me.eq.0) then
c
      if(.not.rtdb_get(irtdb,'md:mdold',mt_int,1,mdold)) mdold=1
c
      open(unit=lfnref,file=filref(1:index(filref,' ')-1),
     + status='old',err=9999)
      rewind(unit=lfnref)
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
      read(lfnref,1002) nbxt,npbt
 1002 format(2i5)
c
      read(lfnref,1007) box(1)
      read(lfnref,1008) box(2)
      read(lfnref,1009) box(3)
 1007 format(f12.6)
 1008 format(12x,f12.6)
 1009 format(24x,f12.6)
c
      read(lfnref,1000) cdum
      read(lfnref,1000) cdum
c
      read(lfnref,1003) nwm,nwa,nsm,nsa
 1003 format(4i10)
c
      if(mdold.ne.0) then
c
      do 2 i=1,nwm
      read(lfnref,1004) (rdum,j=1,6*nwa)
 1004 format(2x,6f13.8)
      read(lfnref,1000) cdum
    2 continue
c
      do 3 i=1,nsa
      read(lfnref,1005) (x(i,j),j=1,3)
 1005 format(2x,3f13.8)
    3 continue
c
      else
c
      do 4 i=1,nwm*(nwa+1)
      read(lfnref,1000) cdum
    4 continue
c
      do 5 i=1,nsa
      read(lfnref,1006) (x(i,j),j=1,3)
 1006 format(2x,3f13.8)
    5 continue
c
      endif
c
      close(unit=lfnref)
      endif
c
c     broadcast coordinates to all nodes
c
      if(np.gt.1) then
      call ga_brdcst(mag_d02,x(1,1),3*nsa*ma_sizeof(mt_dbl,1,mt_byte),0)
      endif
c
      do 6 i=1,nsa
      if(i.eq.1) then
      xsmin=x(1,3)
      xsmax=x(1,3)
      else
      if(x(i,3).lt.xsmin) xsmin=x(i,3)
      if(x(i,3).gt.xsmax) xsmax=x(i,3)
      endif
    6 continue
c
      xsmin=xsmin-0.5d0
      xsmax=xsmax+0.5d0
c
      return
c
 9999 continue
      call md_abort('Unable to open reference file',0)
      return
      end
      subroutine ana_center_ref(x)
c
c $Id$
c
      implicit none
c
#include "ana_common.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "rtdb.fh"
c
      real*8 x(nsa,3)
      integer i
c
c     center and fold
c
      if(me.eq.0.and.icent.gt.0) then
      if(.not.ltop) call md_abort('No topology for center',0)
      call ana_center_fold(x,x,int_mb(i_imol),0)
      endif
c
c     broadcast coordinates to all nodes
c
      if(np.gt.1) then
      call ga_brdcst(mag_d02,x(1,1),3*nsa*ma_sizeof(mt_dbl,1,mt_byte),0)
      endif
c
      do 1 i=1,nsa
      if(i.eq.1) then
      xsmin=x(1,3)
      xsmax=x(1,3)
      else
      if(x(i,3).lt.xsmin) xsmin=x(i,3)
      if(x(i,3).gt.xsmax) xsmax=x(i,3)
      endif
    1 continue
c
      xsmin=xsmin-0.5d0
      xsmax=xsmax+0.5d0
c
      return
      end
