      logical function pre_sgmsiz(lfnsgm,filsgm,
     + natm,nbnd,nang,ndih,nimp)
c
c $Id: pre_sgmsiz.F,v 1.3 1998-04-29 18:35:25 d3j191 Exp $
c
c     function to determine dimensions for segment
c
      implicit none
c
      integer lfnsgm,natm,nbnd,nang,ndih,nimp
      character*255 filsgm
c
      integer length
      character*80 card
c
      length=index(filsgm,' ')-1
      open(unit=lfnsgm,file=filsgm(1:length),form='formatted',
     + status='old',err=9999)
c
    1 continue
      read(lfnsgm,1000,err=9999,end=9999) card
 1000 format(a)
      if(card(1:1).eq.'#'.or.card(1:1).eq.'$') goto 1
      read(card,1001,err=9999) natm,nbnd,nang,ndih,nimp
 1001 format(5i5)
      close(unit=lfnsgm)
      pre_sgmsiz=.true.
      return
c
 9999 continue
      pre_sgmsiz=.false.
      return
      end
