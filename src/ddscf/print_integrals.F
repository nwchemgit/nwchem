c
c     wrapper utility function to print ao integrals
      subroutine print_integrals(basis, oskel)
c     
C$Id$
c
      implicit none
c     
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "stdio.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "util.fh"
c
      integer basis
      logical oskel
c
      if (util_print('kinetic',print_never)) then
          write(6,*) "===================================="
          write(6,*) "Begin kinetic 1-e integrals"
          write(6,*) "===================================="
          call print_1eints(basis, basis, 'kinetic', oskel)
          write(6,*) "===================================="
          write(6,*) "End kinetic 1-e integrals"
          write(6,*) "===================================="
      endif
c
      if (util_print('potential',print_never)) then
          write(6,*) "===================================="
          write(6,*) "Begin potential 1-e integrals"
          write(6,*) "===================================="
          call print_1eints(basis, basis, 'potential', oskel)
          write(6,*) "===================================="
          write(6,*) "End potential 1-e integrals"
          write(6,*) "===================================="
      endif
c
      if (util_print('overlap',print_never)) then
          write(6,*) "===================================="
          write(6,*) "Begin overlap 1-e integrals"
          write(6,*) "===================================="
          call print_1eints(basis, basis, 'overlap', oskel)
          write(6,*) "===================================="
          write(6,*) "End overlap 1-e integrals"
          write(6,*) "===================================="
      endif
c
      if (util_print('ao2eints',print_never)) then
          write(6,*) "===================================="
          write(6,*) "Begin all 2-e integrals"
          write(6,*) "===================================="
          call print_aoints(basis)
          write(6,*) "===================================="
          write(6,*) "End all 2-e integrals"
          write(6,*) "===================================="
      endif
c
      if (util_print('allaoints',print_never)) ! Debug
     $     call print_aoints(basis)
c
      return
      end
