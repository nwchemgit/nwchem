      subroutine pre_input(irtdb)
c
      implicit none
c
#include "rtdb.fh"
#include "inp.fh"
#include "global.fh"
#include "mafdecls.fh"
c
      integer irtdb
      character*255 item
c
      character*3 slvnam
      character*10 slvmdl
      character*80 ffield,project,sysnam
      character*255 calc,dir_s,dir_x,dir_u,dir_t
      character*255 filnam,home,card
      integer newtop,newseq
c
      integer len
c
      if(ga_nodeid().eq.0) then
      write(6,*)
      write(6,*)
      call util_print_centered(6,'PREPARE Input Module',40,.true.)
      write(6,*)
      write(6,*)
      endif
c
c     set the defaults
c     ----------------
c
      ffield=' '
      dir_s=' '
      dir_x=' '
      dir_u=' '
      dir_t=' '
      sysnam=' '
      slvnam='HOH'
      slvmdl='spce      '
      newtop=0
      newseq=0
c
c     read the prepare input block
c     ----------------------------
c
    1 continue
c
      if(.not.inp_read()) call errquit('pre_input: inp_read failed',0)
c
    2 continue
c
      if(.not.inp_a(item)) goto 1
c
    3 continue
c
c     print & noprint
c     ---------------
c
      if(inp_compare(.false.,'print',item).or.
     + inp_compare(.false.,'noprint',item)) then
      call util_print_input(irtdb,'prepare')
      goto 1
      endif
c
c
c     force field
c     -----------
c
      if(inp_compare(.false.,'amber',item)) then
      ffield='amber '
      goto 2
      endif
      if(inp_compare(.false.,'charmm',item)) then
      ffield='charmm '
      goto 2
      endif
      if(inp_compare(.false.,'gromos',item)) then
      ffield='gromos '
      goto 2
      endif
      if(inp_compare(.false.,'oplsa',item)) then
      ffield='oplsa '
      goto 2
      endif
      if(inp_compare(.false.,'oplsg',item)) then
      ffield='oplsg '
      goto 2
      endif
c
c     directories
c     -----------
c
      if(inp_compare(.false.,'standard',item)) then
      if(.not.inp_a(dir_s)) call errquit('pre_input: error dir_s',0)
      goto 2
      endif
      if(inp_compare(.false.,'extensions',item)) then
      if(.not.inp_a(dir_x)) call errquit('pre_input: error dir_x',0)
      goto 2
      endif
      if(inp_compare(.false.,'user',item)) then
      if(.not.inp_a(dir_u)) call errquit('pre_input: error dir_u',0)
      goto 2
      endif
      if(inp_compare(.false.,'temporary',item)) then
      if(.not.inp_a(dir_t)) call errquit('pre_input: error dir_t',0)
      goto 2
      endif
c
c     sysnam
c     ------
c
      if(inp_compare(.false.,'system',item)) then
      if(.not.inp_a(sysnam)) call errquit('pre_input: error sysnam',0)
      goto 2
      endif
c
c     slvnam & slvmdl
c     ---------------
c
      if(inp_compare(.false.,'solvent',item)) then
   20 continue
      if(.not.inp_a(item)) call errquit('pre_input: solvent',9999)
      if(inp_compare(.false.,'name',item)) then
      if(.not.inp_a(item)) call errquit('pre_input: solvent',9999)
      slvnam=item(1:3)
      goto 20
      else
      if(inp_compare(.false.,'model',item)) then
      if(.not.inp_a(slvmdl)) call errquit('pre_input: solvent',9999)
      goto 20
      endif
      endif
      endif
c
c     newtop & newseq
c     ---------------
c
      if(inp_compare(.false.,'newtop',item)) then
      newtop=1
      goto 2
      endif
      if(inp_compare(.false.,'newseq',item)) then
      newseq=1
      newtop=1
      goto 2
      endif
c
c     end
c     ---
c
      if(.not.inp_compare(.false.,'end',item)) then
      write(6,9999) item
 9999 format('pre_input: unknown field ',a)
      call errquit('pre_input: error reading input file',0)
      endif
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
      if(card(1:len+2).eq.ffield(1:len)//'_s')
     + dir_s=card(len+4:index(card(len+4:255),' ')+len+3)
      if(card(1:len+2).eq.ffield(1:len)//'_x')
     + dir_x=card(len+4:index(card(len+4:255),' ')+len+3)
      if(card(1:len+2).eq.ffield(1:len)//'_u')
     + dir_u=card(len+4:index(card(len+4:255),' ')+len+3)
      if(card(1:len+2).eq.ffield(1:len)//'_t')
     + dir_t=card(len+4:index(card(len+4:255),' ')+len+3)
      goto 11
  199 continue
      close(unit=69)
    9 continue
      if(sysnam(1:1).eq.' ')
     + call errquit('pre_input: no system name specified',0)
      if(ffield(1:1).eq.' ')
     + call errquit('pre_input: no force field specified',0)
c
      len=index(dir_s,' ')-1
      if(len.gt.0.and.dir_s(len:len).ne.'/') dir_s(len+1:len+2)='/ '
      len=index(dir_x,' ')-1
      if(len.gt.0.and.dir_x(len:len).ne.'/') dir_x(len+1:len+2)='/ '
      len=index(dir_u,' ')-1
      if(len.gt.0.and.dir_s(len:len).ne.'/') dir_u(len+1:len+2)='/ '
      len=index(dir_t,' ')-1
      if(len.gt.0.and.dir_t(len:len).ne.'/') dir_t(len+1:len+2)='/ '
c
c     put on rtdb
c     -----------
c
      if(.not.rtdb_cput(irtdb,'prep:ffield',1,ffield))
     + call errquit('pre_input: rtdb_put failed',0)
c
      if(.not.rtdb_cput(irtdb,'prep:dir_s',1,dir_s))
     + call errquit('pre_input: rtdb_put failed',0)
      if(.not.rtdb_cput(irtdb,'prep:dir_x',1,dir_x))
     + call errquit('pre_input: rtdb_put failed',0)
      if(.not.rtdb_cput(irtdb,'prep:dir_u',1,dir_u))
     + call errquit('pre_input: rtdb_put failed',0)
      if(.not.rtdb_cput(irtdb,'prep:dir_t',1,dir_t))
     + call errquit('pre_input: rtdb_put failed',0)
c
      if(.not.rtdb_cput(irtdb,'prep:sysnam',1,sysnam))
     + call errquit('pre_input: rtdb_put failed',0)
c
      if(.not.rtdb_cput(irtdb,'prep:slvnam',1,slvnam))
     + call errquit('pre_input: rtdb_put failed',0)
      if(.not.rtdb_cput(irtdb,'prep:slvmdl',1,slvmdl))
     + call errquit('pre_input: rtdb_put failed',0)
c
      if(.not.rtdb_put(irtdb,'prep:newtop',mt_int,1,newtop))
     + call errquit('pre_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'prep:newseq',mt_int,1,newseq))
     + call errquit('pre_input: rtdb_put failed',0)
c
      return
      end
