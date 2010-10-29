      logical function pre_sgmsiz(lfnsgm,filsgm,
     + natm,nbnd,nang,ndih,nimp,nparms)
c
c $Id$
c
c     function to determine dimensions for segment
c
      implicit none
c
      integer lfnsgm,natm,nbnd,nang,ndih,nimp,nparms
      character*255 filsgm
c
      integer length,npar
      character*80 card
      real*8 version
c
      length=index(filsgm,' ')-1
      open(unit=lfnsgm,file=filsgm(1:length),form='formatted',
     + status='old',err=9999)
c      print*,filsgm(1:length)
c
    1 continue
      read(lfnsgm,1000,err=9999,end=9999) card
 1000 format(a)
      if(card(1:1).eq.'#'.or.card(1:1).eq.'$') goto 1
      read(card,1001) version
 1001 format(f12.6)
      read(lfnsgm,1002,err=9999) natm,nbnd,nang,ndih,nimp,npar
 1002 format(5i5,5x,i5)
c      WRITE(*,1002) natm,nbnd,nang,ndih,nimp,npar
      close(unit=lfnsgm)
      nparms=max(nparms,npar)
      pre_sgmsiz=.true.
      return
c
 9999 continue
      pre_sgmsiz=.false.
      return
      end
