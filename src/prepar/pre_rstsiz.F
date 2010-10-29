      logical function pre_rstsz(lfnrst,filrst,nwm,nwa,nsa,nwmc)
c
c $Id$
c
      implicit none
c
      integer lfnrst,nwm,nwa,nwmc
      character*255 filrst
c
      character*1 cdummy
      integer i,nsa,npbtyp,kinsep
c
      open(unit=lfnrst,file=filrst(1:index(filrst,' ')-1),
     + status='old',err=9999)
      rewind(lfnrst)
c
      do 1 i=1,4
      read(lfnrst,1000) cdummy
 1000 format(a1)
    1 continue
      read(lfnrst,1001) npbtyp
 1001 format(i5)
      read(lfnrst,1002) kinsep
 1002 format(40x,i5)
      read(lfnrst,1000) cdummy
      if(kinsep.eq.1) then
      do 2 i=1,3
      read(lfnrst,1000) cdummy
    2 continue
      endif
c
      read(lfnrst,1003) nwm,nwa,nsa,nwmc
 1003 format(2i10,10x,2i10)
c
      close(unit=lfnrst)
c
c      nwm=nwm-nwmc
c
      pre_rstsz=.true.
      return
c
 9999 continue
      pre_rstsz=.false.
      return
      end
      logical function pre_rstsiz(lfnrst,filrst,nwm,nwa,nsa,nwmc)
c
c $Id$
c
      implicit none
c
      integer lfnrst,nwm,nwa,nwmc
      character*255 filrst
c
      character*1 cdummy
      integer i,j,nsa,npbtyp,nbxtyp,nsm,nhist
      real*8 rdummy
c
      open(unit=lfnrst,file=filrst(1:index(filrst,' ')-1),
     + status='old',err=9999)
      rewind(lfnrst)
c
      do 1 i=1,3
      read(lfnrst,1000) cdummy
 1000 format(a1)
    1 continue
      read(lfnrst,1006) nhist
 1006 format(32x,i5)
      if(nhist.gt.0) then
      do 2 i=1,nhist
      read(lfnrst,1000) cdummy
    2 continue
      endif
      read(lfnrst,1002) npbtyp,nbxtyp,(rdummy,j=1,9)
 1002 format(2i5,/,(3f12.6))
      read(lfnrst,1003) rdummy
 1003 format(e12.5)
      read(lfnrst,1004) rdummy,rdummy,rdummy
 1004 format(3f12.6)
      read(lfnrst,1005) nwm,nwa,nsm,nsa,nwmc
 1005 format(7i10,2i5)
c
      close(unit=lfnrst)
c
c      nwm=nwm-nwmc
c
      pre_rstsiz=.true.
      return
c
 9999 continue
      pre_rstsiz=.false.
      return
      end
