      logical function pre_admin(lfncmd,filcmd,lfnout,lfnfrg)
c
c $Id: pre_admin.F,v 1.1 2003-07-31 12:35:45 d3j191 Exp $
c
      implicit none
c
#include "util.fh"
#include "pre_common.fh"
c
      integer lfncmd,lfnout,lfnfrg
      character*255 filcmd,card,string
      character*255 dir,file1,file2
      integer i,length
      integer natm,nparms,ipardef
c
      logical ladmin
c
      ladmin=.false.
c
c     check if the command file exists
c     --------------------------------
c
      open(unit=lfncmd,file=filcmd(1:index(filcmd,' ')-1),
     + form='formatted',status='old',err=99)
c
    1 continue
c
      read(lfncmd,3000,end=9,err=999) card
 3000 format(a)
c
      if(card(1:6).eq.'master') then
      ladmin=.true.
      string=card(8:255)
      dir='./ '
      file1=string
      length=index(string,' ')-1
      do 2 i=length,1,-1
      if(string(i:i).eq.'/') then
      dir=string(1:i)//' '
      file1=string(i+1:length)//' '
      goto 3
      endif
    2 continue
    3 continue
c
      write(filmst,'(a,a)') dir(1:index(dir,' ')-1),
     + file1(1:index(file1,' '))
      length=index(filmst,' ')-1
c
      open(unit=lfnmst,file=filmst(1:length),form='formatted',
     + status='old',err=9999)
c
      if(util_print('sequence',print_high)) then
      write(lfnout,2001) filmst(1:length)
 2001 format(' Reading master fragment',t40,a,/)
      endif
c
    4 continue
c
      read(lfnmst,1000,end=9999,err=9999) card
 1000 format(a)
c
      if(card(1:1).eq.'#'.or.card(1:1).eq.'$') goto 4
c
      read(card,1001) natm,nparms,ipardef
 1001 format(3i5)
c
      close(lfnmst,status='keep')
c
      endif
c
      goto 1
c
    9 continue
      close(unit=lfncmd)
c
   99 continue
  999 continue
      pre_admin=ladmin
      return
 9999 continue
      call errquit('Error in pre_master',0)
      return
      end
