      logical function pre_slvsiz(lfnslv,filslv,nwm,nwa)
c
c $Id$
c
      implicit none
c
      integer lfnslv,nwm,nwa
      character*255 filslv
c
      character*1 cdummy
      integer i,nsa,npbtyp
c
      open(unit=lfnslv,file=filslv(1:index(filslv,' ')-1),
     + status='old',err=9999)
c
      do 1 i=1,4
      read(lfnslv,1000) cdummy
 1000 format(a1)
    1 continue
      read(lfnslv,1001) npbtyp
 1001 format(i5)
      do 2 i=1,5
      read(lfnslv,1000) cdummy
    2 continue
c
      read(lfnslv,1002) nwm,nwa,nsa
 1002 format(2i10,10x,i10)
c
      close(unit=lfnslv)
c
      if(nsa.ne.0) goto 9999
      if(npbtyp.eq.0) goto 9999
c
      pre_slvsiz=.true.
      return
c
 9999 continue
      pre_slvsiz=.false.
      return
      end
