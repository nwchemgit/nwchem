      subroutine xtb_rdinput(rtdb)
      implicit none
#include "rtdb.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "xtb.fh"
#include "util_params.fh"

      integer rtdb
      character(len=4) method

      if (.not.rtdb_get(rtdb,'xtb:acc',mt_dbl,1,acc)) acc=1.0d0
      if (.not.rtdb_get(rtdb,'xtb:kt',mt_dbl,1,kt)) then
        kt = 300.0d0*kb_au
      endif
      if (.not.rtdb_get(rtdb,'xtb:uhf',mt_int,1,uhf)) uhf=0
      if (.not.rtdb_cget(rtdb,'xtb:method',1,method)) method='gfn2'
      if (.not.rtdb_get(rtdb,'xtb:verb',mt_int,1,verb)) verb=1
      if (.not.rtdb_get(rtdb,'xtb:nspin',mt_int,1,nspin)) nspin=1
      if (.not.rtdb_get(rtdb,'xtb:guess',mt_log,1,guess)) guess=.false.
      if (.not.rtdb_get(rtdb,'xtb:broyden',mt_dbl,1,broyden)) 
     $   broyden=0.4d0

      ! Figure restart file name
      if(.not.rtdb_cget(rtdb,'xtb:restart_out',1,restart_out))
     $  restart_out = ' '

      if (restart_out.eq.' ') then
        call util_file_name('xtbrestart',.false.,.false.,restart_out)
      endif  

      !Resolve restart file names to full paths
      call util_file_name_resolve(restart_out,.false.)

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
     $  call errquit('xtb_rdinput: cannot put in rtdb',2,RTDB_ERR) 
      if (.not.rtdb_put(rtdb,'xtb:uhf',mt_int,1,uhf))
     $  call errquit('xtb_rdinput: cannot put in rtdb',3,RTDB_ERR) 
      if (.not.rtdb_put(rtdb,'xtb:verb',mt_int,1,verb))
     $  call errquit('xtb_rdinput: cannot put in rtdb',4,RTDB_ERR) 
      if (.not.rtdb_put(rtdb,'xtb:guess',mt_log,1,guess))
     $  call errquit('xtb_rdinput: cannot put in rtdb',5,RTDB_ERR) 
      if (.not.rtdb_cput(rtdb,'xtb:method',1,method))
     $  call errquit('xtb_rdinput: cannot put in rtdb',6,RTDB_ERR) 
      if (.not.rtdb_put(rtdb,'xtb:nspin',mt_int,1,nspin))
     $  call errquit('xtb_rdinput: cannot put in rtdb',7,RTDB_ERR) 
      if (.not.rtdb_put(rtdb,'xtb:broyden',mt_dbl,1,broyden)) 
     $  call errquit('xtb_rdinput: cannot put in rtdb',8,RTDB_ERR) 
      if(.not.rtdb_cput(rtdb,'xtb:restart_out',1,restart_out))
     $  call errquit('xtb_rdinput: cannot put in rtdb',9,RTDB_ERR) 

      end subroutine
