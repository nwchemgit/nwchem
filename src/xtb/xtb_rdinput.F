      subroutine xtb_rdinput(rtdb)
      implicit none
#include "rtdb.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "xtb.fh"

      integer rtdb
      character(len=4) method

      if (.not.rtdb_get(rtdb,'xtb:acc',mt_dbl,1,acc)) acc = 1.0d0
      if (.not.rtdb_get(rtdb,'xtb:kt',mt_dbl,1,kt)) then
        kt = 300.0d0*3.166808578545117d-06
      endif
      if (.not.rtdb_get(rtdb,'xtb:uhf',mt_int,1,uhf)) uhf = 0
      if (.not.rtdb_cget(rtdb,'xtb:method',1,method)) method = 'gfn2'

      if (method.eq.'gfn2') then
        gfn2 = .true.
        gfn1 = .false.
        gff = .false.
      elseif (method.eq.'gfn1') then
        gfn2 = .false.
        gfn1 = .true.
        gff = .false.
      elseif (method.eq.'gff') then
        gfn2 = .false.
        gfn1 = .false.
        gff = .true.
      endif

      if (.not.rtdb_put(rtdb,'xtb:acc',mt_dbl,1,acc))
     $  call errquit('xtb_rdinput: cannot put in rtdb',1,RTDB_ERR) 
      if (.not.rtdb_put(rtdb,'xtb:kt',mt_dbl,1,kt))
     $  call errquit('xtb_rdinput: cannot put in rtdb',1,RTDB_ERR) 
      if (.not.rtdb_put(rtdb,'xtb:uhf',mt_int,1,uhf))
     $  call errquit('xtb_rdinput: cannot put in rtdb',2,RTDB_ERR) 
      if (.not.rtdb_cput(rtdb,'xtb:method',1,method))
     $  call errquit('xtb_rdinput: cannot put in rtdb',3,RTDB_ERR) 

      end subroutine
