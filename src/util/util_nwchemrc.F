*
*
*
*
* these should use the inp modules like the rest of nwchem
* see bas_library_file in .../basis/bas_input.F
* right now it is key value but later on we might want key
* array of values etc.  The inp_a, inp_f, inp_i routines
* provide all the flexibility we need with the ability to have
* physical/logical line style input. 
*
*
      logical function util_nwchemrc_put(key,value)
c
      implicit none
c
* $Id: util_nwchemrc.F,v 1.3 1998-04-24 07:32:48 d3e129 Exp $
c
      character*(*) key,value
c
      character*255 home,nwchrc
      character*1020 mystring
      character*1 dummy
      integer lenkey
      intrinsic index
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
      mystring = key(1:index(key,' '))//value(1:index(value,' ')-1)
      write(69,6901) mystring
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
