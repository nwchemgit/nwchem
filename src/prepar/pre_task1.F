      logical function pre_task1(irtdb)
c
      implicit none
c
#include "rtdb.fh"
c
      logical scf,esp
      external scf,esp
c
      integer irtdb
c
      pre_task1=.false.
c
      if(.not.rtdb_cput(irtdb,'scf:print',1,'none'))
     + call errquit('pre_task1: rtdb_cput failed',9999)
c
      if(scf(irtdb)) then
c
      if(.not.rtdb_cput(irtdb,'esp:print',1,'none'))
     + call errquit('pre_task1: rtdb_cput failed',9999)
c
      if(esp(irtdb)) then
      pre_task1=.true.
      endif
      endif
c
      return
      end
