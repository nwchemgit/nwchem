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
c     $Id$
c
c     Eventually these routines should use the inp modules
c     Currently these are (key,value) but later may have to be extended
c     to array of values, etc.. 
c
c     This routine has been made Windows proof by allowing spaces in file
c     names and values of the key-value pairs.
c
#include "inp.fh"
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
      nwchrc=home(1:inp_strlen(home))//'/.nwchemrc '
      open(unit=69,file=nwchrc(1:inp_strlen(nwchrc)),form='formatted',
     + status='old',err=9999)
    1 continue
      read(69,6900,end=9,err=9) dummy
 6900 format(a1)
      goto 1
    9 continue
      mystring = key(1:index(key,' '))//value(1:inp_strlen(value))
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
c
c-----------------------------------------------------------------------
c
      logical function util_nwchemrc_get(key,value)
c
c     This function looks for a particular key value pair in nwchemrc
c     files. The key is given in the variable with the corresponding
c     name. If this key is found its value will be returned in the
c     variable value and the function will return .true.,
c     otherwise the function will return .false.
c
c     Because on UNIX/LINUX it is common to install a system wide
c     configuration file in /etc, and perhaps a more specific one in
c     a user-s home directory the original logic has been extended.
c     Now we first look for $HOME/.nwchemrc if the key is found in there
c     we will return the corresponding value. If key was not found or
c     the file does not exist we try /etc/nwchemrc. If both sources
c     fail the function returns .false.
c
c     On the key-value lines all information following a # token will
c     be ignored. Furthermore, this routine has been made Windows proof
c     by allowing spaces in the value of the key-value pair.
c
      implicit none
#include "inp.fh"
#include "stdio.fh"
c
      character*(*) key,value
c
      character*255 home,nwchrc,card,rightcard
      integer indexv,lenkey,lencrd,istart,util_getblnk
      integer ipass
      external util_getblnk
      logical found
      intrinsic index
c
      util_nwchemrc_get=.false.
      lenkey=index(key,' ')-1
      if(lenkey.le.-1) lenkey=len(key) ! key does not contain space
      if(lenkey.eq.1)  return          ! key is empty
c
      ipass=0
      found=.false.
c
c     Loop over different nwchemrc location options as long as we have
c     not found key.
c
      do while (.not.found)
        ipass=ipass+1
        if (ipass.eq.1) then
          call util_getenv('HOME',home)
          nwchrc=home(1:inp_strlen(home))//'/.nwchemrc '
        elseif (ipass.eq.2) then
          nwchrc="/etc/nwchemrc"
        else
          return
        endif
        open(unit=69,file=nwchrc(1:inp_strlen(nwchrc)),
     +       form='formatted',status='old',action='read',err=9999)
        rewind(unit=69)
c
c       Loop over the lines in the file as long as we have not found
c       key.
c
        do while (.not.found)
          read(69,6900,end=999,err=999) card
 6900     format(a)
c
c         Find the start (after skipping any white spaces) and length
c         of the key.
c
          istart=util_getblnk(card,255)
          lencrd=index(card(istart:255),' ')-1
c
c         Check whether the length and text of the key on the card and
c         the requested key match. If so, we have found the key and
c         extract to corresponding value.
c
          if (lencrd.eq.lenkey.and.
     +        card(istart:lencrd+istart).eq.key(1:lenkey)) then
            indexv=lencrd+istart+1
            rightcard=card(indexv:255)
            istart=util_getblnk(rightcard,255)
            lencrd=index(rightcard(istart:255),'#')-1
            if (lencrd.le.-1) lencrd=255-istart+1
            lencrd=inp_strlen(rightcard(istart:lencrd+istart-1))-1
            value=rightcard(istart:lencrd+istart)
            found=.true.
          endif
        enddo ! while(.not.found)
 999    close(unit=69)
        goto 9996
 9999   continue
cdbg    write(luout,*)'*** Key: ',key(1:lenkey),' Not found: ',nwchrc
 9996   continue
c
      enddo ! while(ipass.lt.2.and..not.found)
cdbg  write(luout,*)'*** Key: ',key(1:lenkey),' Value: ',value
c
      util_nwchemrc_get=found
c
      return
      end
c
c-----------------------------------------------------------------------
c
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
