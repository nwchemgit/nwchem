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
      character*80 ffield,projct,calc,dir_s,dir_x,dir_u,dir_t
      character*80 sysnam
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
      ffield='amber '
      dir_s='/usr/people/d3j191/data/amber_s/ '
      dir_x='/usr/people/d3j191/data/amber_x/ '
      dir_u='/usr/people/d3j191/data/amber_u/ '
      dir_t='/usr/people/d3j191/test/ '
      sysnam=' '
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
c     end
c     ---
c
      if(.not.inp_compare(.false.,'end',item)) then
      write(6,9999) item
 9999 format('pre_input: unknown field ',a)
      call errquit('pre_input: error reading input file',0)
      endif
c
c     do some tests
c     -------------
c
      if(sysnam(1:1).eq.' ')
     + call errquit('pre_in put: no sysnam specified',0)
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
      return
      end

