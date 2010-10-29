      subroutine tce_reference_check(rtdb,ref)
c
c $Id$
c
c Written by Jeff Hammond, January 2008.
c
      implicit none
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "util.fh"
#include "msgids.fh"
#include "stdio.fh"
#include "errquit.fh"
c#include "tce.fh"
#include "tce_main.fh"
c#include "tce_prop.fh"
c#include "tce_ints.fh"
c#include "tce_amps.fh"
c#include "tce_diis.fh"
      integer rtdb
      double precision ref
      logical scf
      external scf
      logical dft_energy
      external dft_energy
c
c     =========================
c     Ground-state HF/DFT first
c     =========================
c
      if (.not.rtdb_get(rtdb,'tce:reference',mt_int,1,reference)) then
        reference = 1
        if (.not.rtdb_put(rtdb,'tce:reference',mt_int,1,reference))
     1    call errquit('tce_property: failed to write reference',0,
     2    RTDB_ERR)
      endif
      if (reference.eq.0) then
        if (.not.dft_energy(rtdb)) return
        if (.not.rtdb_get(rtdb,'dft:energy',mt_dbl,1,ref))
     1    call errquit('tce_property: failed to get dft energy',0,
     2    RTDB_ERR)
      else if (reference.eq.1) then
        if (.not.scf(rtdb)) return
        if (.not.rtdb_get(rtdb,'scf:energy',mt_dbl,1,ref))
     1    call errquit('tce_property: failed to get scf energy',0,
     2    RTDB_ERR)
      else
        call errquit('tce_property: illegal reference',reference,
     1    UNKNOWN_ERR)
      endif
      return
      end
