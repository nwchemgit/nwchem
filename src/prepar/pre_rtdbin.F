      logical function pre_rtdbin(irtdb,ffield,
     + dir_s,dir_x,dir_u,dir_t,sysnam)
c
c     function to read input for prepare module from rtdb
c
      implicit none
c
#include "rtdb.fh"
c
      integer irtdb
      character*80 ffield,dir_s,dir_x,dir_u,dir_t,sysnam
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
      pre_rtdbin=.true.
      return
      end
