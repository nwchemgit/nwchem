      subroutine ana_rdref(x)
c
c $Id: ana_rdref.F,v 1.1 1999-01-29 19:09:34 d3j191 Exp $
c
      implicit none
c
#include "ana_common.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c
      real*8 x(nsa,3)
      character*1 cdum
      integer i,j,nwm,nwa,nsai,msai,npprev,nbxp,nbyp,nbzp
      real*8 rdum
c
      if(.not.ma_verify_allocator_stuff()) print*,'Oopsi 1'
c
      if(me.eq.0) then
      open(unit=lfnref,file=filref(1:index(filref,' ')-1),
     + status='old',err=9999)
      rewind(unit=lfnref)
c      print*,'opened ',filref(1:index(filref,' ')-1)
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
      read(lfnref,1003) nwm,nwa,npprev
 1003 format(2i10,50x,i5)
      if(npprev.gt.0) then
      read(lfnref,1004) nbxp,nbyp,nbzp
 1004 format(3i5)
      read(lfnref,1005) (rdum,i=1,nbxp)
      read(lfnref,1005) (rdum,i=1,nbyp)
      read(lfnref,1005) (rdum,i=1,nbzp)
 1005 format(4e20.12)
      endif
      read(lfnref,1006) (rdum,i=1,3*nwm*nwa)
 1006 format(6f13.8)
      read(lfnref,1006) ((x(i,j),j=1,3),i=1,nsa)
      close(unit=lfnref)
      endif
c
c     broadcast coordinates to all nodes
c
      if(np.gt.1) then
      call ga_brdcst(mag_d02,x(1,1),3*nsa*ma_sizeof(mt_dbl,1,mt_byte),0)
      endif
c
      return
c
 9999 continue
      call errquit('Unable to open reference file',0)
      return
      end
