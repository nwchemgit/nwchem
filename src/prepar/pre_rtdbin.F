      logical function pre_rtdbin(irtdb,ffield,
     + dir_s,dir_x,dir_u,dir_t,source,sysnam,calc,slvnam,slvmdl,
     + newtop,newseq,newrst)
c
c     function to read input for prepare module from rtdb
c
      implicit none
c
#include "rtdb.fh"
#include "mafdecls.fh"
c
      integer irtdb
      character*3 slvnam
      character*80 source
      character*80 ffield,sysnam,calc
      character*10 slvmdl
      character*255 dir_s,dir_x,dir_u,dir_t
      integer newtop,newseq,newrst
      integer len,ndx
c
      character*255 home,filnam,card
c
      ffield=' '
      dir_s=' '
      dir_x=' '
      dir_u=' '
      dir_t=' '
      source=' '
c
c
c     read from rtdb
c     --------------
c
      if(.not.rtdb_cget(irtdb,'prep:ffield',1,ffield)) continue
c
      if(.not.rtdb_cget(irtdb,'prep:dir_s',1,dir_s)) continue
      if(.not.rtdb_cget(irtdb,'prep:dir_x',1,dir_x)) continue
      if(.not.rtdb_cget(irtdb,'prep:dir_u',1,dir_u)) continue
      if(.not.rtdb_cget(irtdb,'prep:dir_t',1,dir_t)) continue
c
      if(.not.rtdb_cget(irtdb,'prep:sysnam',1,sysnam)) then
      if(.not.rtdb_cget(irtdb,'arg:project',1,sysnam)) then
      if(.not.rtdb_cget(irtdb,'file_prefix',1,sysnam))
     + call errquit('rtdb_get failed on sysnam',9999)
      endif
      endif
c
      if(.not.rtdb_cget(irtdb,'prep:source',1,source)) source=' '
c
      if(.not.rtdb_cget(irtdb,'prep:slvnam',1,slvnam)) slvnam='HOH'
      if(.not.rtdb_cget(irtdb,'prep:slvmdl',1,slvmdl)) slvmdl=' '
c
      if(.not.rtdb_get(irtdb,'prep:newtop',mt_int,1,newtop)) newtop=0
      if(.not.rtdb_get(irtdb,'prep:newseq',mt_int,1,newseq)) newseq=0
      if(.not.rtdb_get(irtdb,'prep:newrst',mt_int,1,newrst)) newrst=0
c
      if(slvnam.eq.'HOH'.and.slvmdl(1:1).eq.' ') slvmdl='spce      '
c
      ndx=index(sysnam,'_')
      len=index(sysnam,' ')
      if(ndx.gt.1.and.ndx.lt.len-1) then
      calc=sysnam(ndx+1:len)
      sysnam=sysnam(1:ndx-1)//' '
      else
      calc='md '
      endif
c
      if(newseq.eq.1) newtop=1
c
c     get defaults for unspecified variables from $HOME/.nwchemrc
c     -----------------------------------------------------------
c 
      call getenv('HOME',home)
c
      len=index(home,' ')-1
      filnam=home(1:len)//'/.nwchemrc '
      len=index(filnam,' ')-1
      open(unit=69,file=filnam(1:len),form='formatted',status='old',
     + err=9)
   10 continue
      read(69,6900,end=99,err=99) card
 6900 format(a)
      if(ffield(1:1).eq.' '.and.card(1:6).eq.'ffield') then
      len=index(card(8:255),' ')
      ffield=card(8:8+len)
      endif
      goto 10
   99 continue
      rewind(unit=69)
      if(ffield(1:1).eq.' ') ffield='amber '
      len=index(ffield,' ')-1
   11 continue
      read(69,6900,end=199,err=199) card
      if(dir_s(1:1).eq.' '.and.card(1:len+2).eq.ffield(1:len)//'_s')
     + dir_s=card(len+4:index(card(len+4:255),' ')+len+3)
      if(dir_x(1:1).eq.' '.and.card(1:len+2).eq.ffield(1:len)//'_x')
     + dir_x=card(len+4:index(card(len+4:255),' ')+len+3)
      if(dir_u(1:1).eq.' '.and.card(1:len+2).eq.ffield(1:len)//'_u')
     + dir_u=card(len+4:index(card(len+4:255),' ')+len+3)
      if(dir_t(1:1).eq.' '.and.card(1:len+2).eq.ffield(1:len)//'_t')
     + dir_t=card(len+4:index(card(len+4:255),' ')+len+3)
      goto 11
  199 continue
      close(unit=69)
    9 continue
c
      pre_rtdbin=.true.
      return
      end
