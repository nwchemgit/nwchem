      logical function util_nwchemrc_put(key,value)
c
      implicit none
c
      character*(*) key,value
c
      character*255 home,nwchrc
      character*1 dummy
      integer lenkey
c
      lenkey=index(key,' ')-1
      if(lenkey.le.0) goto 9999
      call getenv('HOME',home)
      nwchrc=home(1:index(home,' ')-1)//'/.nwchemrc '
      open(unit=69,file=nwchrc(1:index(nwchrc,' ')-1),form='formatted',
     + status='old',err=9999)
    1 continue
      read(69,6900,end=9,err=9) dummy
 6900 format(a1)
      goto 1
    9 continue
      write(69,6901) key(1:index(key,' '))//value(1:index(value,' ')-1)
 6901 format(a)
      close(unit=69)
c
      util_nwchemrc_put=.true.
      return
c
 9999 continue
      util_nwchemrc_put=.false.
      return
      end
      logical function util_nwchemrc_get(key,value)
c
      implicit none
c
      character*(*) key,value
c
      character*255 home,nwchrc,card
      integer indexv,lenkey,lencrd
      logical found
c
      lenkey=index(key,' ')-1
      if(lenkey.le.0) goto 9999
      call getenv('HOME',home)
      nwchrc=home(1:index(home,' ')-1)//'/.nwchemrc '
      open(unit=69,file=nwchrc(1:index(nwchrc,' ')-1),form='formatted',
     + status='old',err=9999)
      rewind(unit=69)
      found=.false.
    1 continue
      read(69,6900,end=9,err=9) card
 6900 format(a)
      lencrd=index(card,' ')-1
      if(lencrd.eq.lenkey.and.card(1:lencrd).eq.key(1:lenkey)) then
      indexv=index(card,' ')+1
      value=card(indexv:indexv+index(card(indexv:255),' '))
      found=.true.
      endif
      goto 1
    9 continue
      close(unit=69)
c
      if(.not.found) goto 9999
      util_nwchemrc_get=.true.
      return
c
 9999 continue
      util_nwchemrc_get=.false.
      return
      end
