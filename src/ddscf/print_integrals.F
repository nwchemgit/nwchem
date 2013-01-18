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
       if (ga_nodeid().eq.0) then
          write(luout,*) "===================================="
          write(luout,*) "Begin kinetic 1-e integrals"
          write(luout,*) "===================================="
       end if
       call print_1eints(basis, basis, 'kinetic', oskel)
       if (ga_nodeid().eq.0) then
          write(luout,*) "===================================="
          write(luout,*) "End kinetic 1-e integrals"
          write(luout,*) "===================================="
       end if
      endif
c
      if (util_print('potential',print_never)) then
       if (ga_nodeid().eq.0) then
          write(luout,*) "===================================="
          write(luout,*) "Begin potential 1-e integrals"
          write(luout,*) "===================================="
       end if
       call print_1eints(basis, basis, 'potential', oskel)
       if (ga_nodeid().eq.0) then
          write(luout,*) "===================================="
          write(luout,*) "End potential 1-e integrals"
          write(luout,*) "===================================="
       end if
      endif
c
      if (util_print('overlap',print_never)) then
       if (ga_nodeid().eq.0) then
          write(luout,*) "===================================="
          write(luout,*) "Begin overlap 1-e integrals"
          write(luout,*) "===================================="
       end if
       call print_1eints(basis, basis, 'overlap', oskel)
       if (ga_nodeid().eq.0) then
          write(luout,*) "===================================="
          write(luout,*) "End overlap 1-e integrals"
          write(luout,*) "===================================="
       end if
      endif
c
      if (util_print('ao2eints',print_never)) then
       if (ga_nodeid().eq.0) then
          write(luout,*) "===================================="
          write(luout,*) "Begin all 2-e integrals"
          write(luout,*) "===================================="
       end if
       call print_aoints(basis)
       if (ga_nodeid().eq.0) then
          write(luout,*) "===================================="
          write(luout,*) "End all 2-e integrals"
          write(luout,*) "===================================="
       end if
      endif
c
      if (util_print('allaoints',print_never)) ! Debug
     $     call print_aoints(basis)
c
      return
      end
