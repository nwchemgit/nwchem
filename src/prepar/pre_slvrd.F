      logical function pre_slvrd(lfnslv,filslv,xw,vw,mwm,mwa,nwm,nwa,
     + box)
c
      implicit none
c
      integer lfnslv,mwm,mwa,nwm,nwa
      character*255 filslv
      real*8 xw(3,mwa,mwm),vw(3,mwa,mwm),box(3)
c
      character*1 cdummy
      integer i,j,k,nsa,npbtyp
c
      open(unit=lfnslv,file=filslv(1:index(filslv,' ')-1),
     + status='old',err=9999)
c
      do 1 i=1,4
      read(lfnslv,1000) cdummy
 1000 format(a1)
    1 continue
      read(lfnslv,1001) npbtyp,box
 1001 format(i5,3f12.6)
      if(npbtyp.eq.0) goto 9998
      do 2 i=1,2
      read(lfnslv,1000) cdummy
    2 continue
c
      read(lfnslv,1002) nwm,nwa,nsa
 1002 format(2i10,10x,i10)
      if(nwm.ne.mwm) goto 9998
      if(nsa.ne.0) goto 9998
c
      read(lfnslv,1003) (((xw(k,j,i),k=1,3),j=1,nwa),i=1,nwm)
      read(lfnslv,1003) (((vw(k,j,i),k=1,3),j=1,nwa),i=1,nwm)
 1003 format(6f13.8)
c
      close(unit=lfnslv)
c
      pre_slvrd=.true.
      return
c
c
 9998 continue
      close(unit=lfnslv)
 9999 continue
      pre_slvrd=.false.
      return
      end
