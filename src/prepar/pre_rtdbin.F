      logical function pre_rtdbin(irtdb,ffield,
     + dir_s,dir_x,dir_u,dir_t,source,sysnam,calc,slvnam,slvmdl,
     + newtop,newseq,newrst)
c
c     function to read input for prepare module from rtdb
c
      implicit none
c
      logical util_nwchemrc_get
      external util_nwchemrc_get
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
      character*255 home,filnam,key,value
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
      if(ffield(1:1).eq.' ') then
      key='ffield '
      if(util_nwchemrc_get(key,value)) then
      ffield=value(1:index(value,' ')-1)
      endif
      endif
      if(dir_s(1:1).eq.' ') then
      key=ffield(1:index(ffield,' ')-1)//'_s '
      if(util_nwchemrc_get(key,value))
     + dir_s=value(1:index(value,' '))
      endif
      if(dir_x(1:1).eq.' ') then
      key=ffield(1:index(ffield,' ')-1)//'_x '
      if(util_nwchemrc_get(key,value))
     + dir_x=value(1:index(value,' '))
      endif
      if(dir_u(1:1).eq.' ') then
      key=ffield(1:index(ffield,' ')-1)//'_u '
      if(util_nwchemrc_get(key,value))
     + dir_u=value(1:index(value,' '))
      endif
      if(dir_t(1:1).eq.' ') then
      key=ffield(1:index(ffield,' ')-1)//'_t '
      if(util_nwchemrc_get(key,value))
     + dir_t=value(1:index(value,' '))
      endif
c
      write(*,'(a,a)') 'ffield=',ffield(1:index(ffield,' ')-1)
      write(*,'(a,a)') 'dir_s =',dir_s(1:index(dir_s,' ')-1)
      write(*,'(a,a)') 'dir_x =',dir_x(1:index(dir_x,' ')-1)
      write(*,'(a,a)') 'dir_u =',dir_u(1:index(dir_u,' ')-1)
      write(*,'(a,a)') 'dir_t =',dir_t(1:index(dir_t,' ')-1)
c
      pre_rtdbin=.true.
      return
      end
