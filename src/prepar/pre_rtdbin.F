      logical function pre_rtdbin(irtdb,ffield,
     + dir_s,dir_x,dir_u,dir_t,sysnam,slvnam,slvmdl,newtop,newseq)
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
      character*80 ffield,sysnam
      character*10 slvmdl
      character*255 dir_s,dir_x,dir_u,dir_t
      integer newtop,newseq
c
      if(.not.rtdb_cget(irtdb,'prep:ffield',1,ffield))
     + call errquit('rtdb_get failed on ffield',9999)
c
      if(.not.rtdb_cget(irtdb,'prep:dir_s',1,dir_s))
     + call errquit('rtdb_get failed on dir_s',9999)
      if(.not.rtdb_cget(irtdb,'prep:dir_x',1,dir_x))
     + call errquit('rtdb_get failed on dir_x',9999)
      if(.not.rtdb_cget(irtdb,'prep:dir_u',1,dir_u))
     + call errquit('rtdb_get failed on dir_u',9999)
      if(.not.rtdb_cget(irtdb,'prep:dir_t',1,dir_t))
     + call errquit('rtdb_get failed on dir_t',9999)
c
      if(.not.rtdb_cget(irtdb,'prep:sysnam',1,sysnam))
     + call errquit('rtdb_get failed on sysnam',9999)
c
      if(.not.rtdb_cget(irtdb,'prep:slvnam',1,slvnam))
     + call errquit('rtdb_get failed on slvnam',9999)
      if(.not.rtdb_cget(irtdb,'prep:slvmdl',1,slvmdl))
     + call errquit('rtdb_get failed on slvmdl',9999)
c
      if(.not.rtdb_get(irtdb,'prep:newtop',mt_int,1,newtop))
     + call errquit('rtdb_get failed on newtop',9999)
      if(.not.rtdb_get(irtdb,'prep:newseq',mt_int,1,newseq))
     + call errquit('rtdb_get failed on newseq',9999)
c
      if(slvnam.eq.'HOH'.and.slvmdl(1:1).eq.' ') slvmdl='spce      '
c
      pre_rtdbin=.true.
      return
      end
