      function leps_energy(rtdb)
      implicit none
      integer rtdb

#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "util.fh"
#include "global.fh"
      logical leps_energy
c
      write(*,*) "doing leps energy"
      if (.not. rtdb_put(rtdb,'leps:energy',mt_dbl,1,0.0))
     & call errquit('leps_energy:  no energy ', 0,
     & RTDB_ERR)
      return
      END

