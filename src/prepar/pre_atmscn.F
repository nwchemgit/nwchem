      logical function pre_atmscn(lfnsgm,filsgm,atmnam,num,link)
c
c $Id: pre_atmscn.F,v 1.3 1998-04-29 18:35:11 d3j191 Exp $
c
      implicit none
c
      integer lfnsgm,num,link
      character*255 filsgm
      character*4 atmnam
c
      character*80 card
      integer i,length,natm,ilnk
      character*6 name
c
      length=index(filsgm,' ')-1
      open(unit=lfnsgm,file=filsgm(1:length),form='formatted',
     + status='old',err=9999)
c
      num=0
      link=0
c
    1 continue
      read(lfnsgm,1000,err=9999,end=9999) card
 1000 format(a)
      if(card(1:1).eq.'#'.or.card(1:1).eq.'$') goto 1
      read(card,1001,err=9999) natm
 1001 format(i5)
      do 2 i=1,natm
      read(lfnsgm,1002) name,ilnk 
 1002 format(5x,a6,30x,i4)
      if(name(1:4).eq.atmnam) then
      link=ilnk
      num=i
      close(unit=lfnsgm)
      pre_atmscn=.true.
      return
      endif
      read(lfnsgm,1000,err=9999,end=9999) card
    2 continue
      close(unit=lfnsgm)
 9999 continue
      pre_atmscn=.false.
      return
      end
