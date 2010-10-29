      logical function argos_prep_slvrd(lfnslv,filslv,
     + xw,vw,mwm,mwa,nwm,nwa,box)
c
c $Id$
c
      implicit none
c
      integer lfnslv,mwm,mwa,nwm,nwa
      character*255 filslv
      real*8 xw(3,mwa,mwm),vw(3,mwa,mwm),box(3),vlat(3,3)
c
      character*1 cdum
      integer i,j,k,idum,jdum,kdum,nsa,npbtyp
      real*8 rdum
c
      open(unit=lfnslv,file=filslv(1:index(filslv,' ')-1),
     + status='old',err=9999)
c
      do 1 i=1,4
      read(lfnslv,1000) cdum
 1000 format(a1)
    1 continue
      read(lfnslv,1001) npbtyp
 1001 format(i5)
      read(lfnslv,1002) ((vlat(i,j),j=1,3),i=1,3)
 1002 format(3f12.6)
      box(1)=vlat(1,1)
      box(2)=vlat(2,2)
      box(3)=vlat(3,3)
      read(lfnslv,1003) jdum
 1003 format(40x,i5)
      read(lfnslv,1000) cdum
      if(jdum.ne.0) then
      read(lfnslv,1000) cdum
      endif
      read(lfnslv,1004) nwm,nwa,jdum,nsa,idum
 1004 format(4i10,30x,i5)
      if(idum.gt.0) then
      read(lfnslv,1005) idum,jdum,kdum
 1005 format(3i5)
      read(lfnslv,1006) (rdum,i=1,idum)
      read(lfnslv,1006) (rdum,i=1,jdum)
      read(lfnslv,1006) (rdum,i=1,kdum)
 1006 format(4e20.12)
      endif
c
      if(nwm.ne.mwm) goto 9998
      if(nsa.ne.0) goto 9998
c
      do 3 k=1,nwm
      read(lfnslv,1007) ((xw(i,j,k),i=1,3),(vw(i,j,k),i=1,3),j=1,nwa)
 1007 format(2x,6f13.8)
      read(lfnslv,1000) cdum
    3 continue
c
      close(unit=lfnslv)
c
      argos_prep_slvrd=.true.
      return
c
c
 9998 continue
      close(unit=lfnslv)
 9999 continue
      argos_prep_slvrd=.false.
      return
      end
