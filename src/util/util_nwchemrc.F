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
c $Id: util_nwchemrc.F,v 1.6 2000-10-17 23:12:15 edo Exp $
c
c Eventually these routines should use the inp modules
c Currently these are (key,value) but later may have to be extended
c to array of values, etc.. 
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
      call util_getenv('HOME',home)
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
#include "inp.fh"

c
      character*(*) key,value
c
      character*255 home,nwchrc,card,rightcard
      integer indexv,lenkey,lencrd,istart,util_getblnk
      external util_getblnk
      logical found
c
      lenkey=index(key,' ')-1
      if(lenkey.le.0) goto 9999
      call util_getenv('HOME',home)
      nwchrc=home(1:index(home,' ')-1)//'/.nwchemrc '
      open(unit=69,file=nwchrc(1:index(nwchrc,' ')-1),form='formatted',
     + status='old',err=9999)
      rewind(unit=69)
      found=.false.
    1 continue
      read(69,6900,end=9,err=9) card
 6900 format(a)
c
c     get beginning of key after white spaces
c
      istart=util_getblnk(card,255)
      lencrd=index(card(istart:255),' ')-1
      if(lencrd.eq.lenkey.and.card(istart:lencrd+istart).eq.
     .     key(1:lenkey)) then
      indexv=lencrd+istart+1
      rightcard=card(indexv:255)
      istart=util_getblnk(rightcard,255)
      lencrd=index(rightcard(istart:255),' ')-1
      value=rightcard(istart:lencrd+istart)
      found=.true.
      endif
      goto 1
    9 continue
      close(unit=69)
c
      if(.not.found) goto 99999
      util_nwchemrc_get=.true.
c
c     clean up blank spaces
c
      return
c
 9999 write(6,*) ' .nwchemrc not found !!'
      
99999 continue
      util_nwchemrc_get=.false.
      return
      end
      integer function util_getblnk(card,mxlgth)
      implicit none 
      integer mxlgth
      character*(*) card
c
      integer i
      logical ois_ws            ! Inline funtion  
      character*1 xtest
c
      ois_ws(xtest) = (xtest.eq.' ')
c
      util_getblnk=0
      do i=1,mxlgth
         if(.not.ois_ws(card(i:i))) then
            util_getblnk=i
            return
         endif
      enddo
      return
      end
