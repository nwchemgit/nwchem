       logical function tce_mrcc_energy(rtdb)
c
c Parallel version of multireference CC interface
c Jiri Brabec 2011
c
       implicit none
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "sym.fh"
#include "util.fh"
#include "msgids.fh"
#include "stdio.fh"
#include "sf.fh"
#include "inp.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
#include "tce_hetio.fh"
#include "tce_diis.fh"
#include "tce_prop.fh"
#include "tce_restart.fh"

#ifdef MRCC_METHODS
#include "tce_mrcc.fh"
#endif

      integer rtdb
      logical nodezero

      double precision eone,etwo,enrep,energy

#ifdef MRCC_METHODS
      tce_mrcc_energy = .false.
      nodezero=(ga_nodeid().eq.0)
      write(*,*) "I am here1"
      call ga_sync()
      tce_mrcc_energy = .true.

      if (.not.rtdb_put(rtdb,'tce:energy',mt_dbl,1,epsilon+enrep))
     + call errquit('tce_energy: RTDB problem',0,MA_ERR)

       if (nodezero) call util_flush(LuOut)

ckbn mrcc-r-1 CC      call tce_tidy(rtdb)
ckbn mrcc-r-1 CC      call util_print_pop
      return
#else
      CALL ERRQUIT('MRCC Methods not compiled (tce_mrcc_energy.F) 
     + Compile  with "export MRCC_METHODS=y"',0,CAPMIS_ERR)
#endif

      end
