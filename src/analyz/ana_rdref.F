      subroutine ana_rdref(x)
c
c $Id: ana_rdref.F,v 1.3 1999-09-20 20:36:09 d3j191 Exp $
c
      implicit none
c
#include "ana_common.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c
      real*8 x(nsa,3)
      character*1 cdum
      integer i,j,nsai,msai,npprev,nbxp,nbyp,nbzp
      real*8 rdum
      integer nbxt,npbt,nsm
c
      if(.not.ma_verify_allocator_stuff()) print*,'Oopsi 1'
c
      if(me.eq.0) then
      open(unit=lfnref,file=filref(1:index(filref,' ')-1),
     + status='old',err=9999)
      rewind(unit=lfnref)
c
      read(lfnref,1000) cdum
      read(lfnref,1000) cdum
      read(lfnref,1000) cdum
      read(lfnref,1000) cdum
 1000 format(a)
c
      read(lfnref,1001) nbxt,npbt
 1001 format(2i5)
c
      read(lfnref,1000) cdum
      read(lfnref,1000) cdum
      read(lfnref,1000) cdum
c
      read(lfnref,1000) cdum
      read(lfnref,1000) cdum
c
      read(lfnref,1002) nwm,nwa,nsm,nsa
 1002 format(4i10)
c
      do 1 i=1,nwm
      read(lfnref,1006) (rdum,j=1,6*nwa)
 1006 format(2x,6f13.8)
      read(lfnref,1000) cdum
    1 continue
c
      do 2 i=1,nsa
      read(lfnref,1007) (x(i,j),j=1,3)
 1007 format(2x,3f13.8)
    2 continue
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
      return
c
 9999 continue
      call errquit('Unable to open reference file',0)
      return
      end
