      subroutine tddft_dftinit(rtdb,oskel,oadapt)
c
c $Id$
c
c Set various parameters for DFT calculations for the 
c ground state in the context of TDDFT calculations.
c Specifically, oskel (use of symmetry in density matrix
c integral contraction) will be suppressed.
c Written by So Hirata, Feb 2002.
c
      implicit none
#include "errquit.fh"
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
c
      integer rtdb             ! Runtime database handle
      logical oskel            ! Use of symmetry in Fock build
      logical oadapt           ! Symmetry adaptation of orbitals
      logical dummy
c
c --------------------------------------------
c Unset the converged flag of ground-state DFT
c --------------------------------------------
c
      dummy=rtdb_delete(rtdb,'dft:converged')
c
c ---------------------
c Geometry and symmetry
c ---------------------
c
      oskel=.false.
      if (.not.rtdb_put(rtdb,'dft:skeleton',mt_log,1,.false.))
     1  call errquit('tddft_dftinit: failed to set rtdb param',0,
     &       RTDB_ERR)
c
c We need to symmetry adapt the orbitals to characterize excited state
c wave function symmetries
c
c     oadapt=.false.
c     if (.not.rtdb_put(rtdb,'dft:adapt',mt_log,1,oadapt))
c    1  call errquit('tddft_dftinit: failed to set rtdb param',0)
c
c Return
c
      return
      end
