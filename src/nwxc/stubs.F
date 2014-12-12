C> \brief Stub routines for nwxc_fnc_print
C> 
C> These stubs are only to be used for the nwxc_fnc_print program.
C> They resolve symbols that are not needed to print the functional
C> expression. In now circumstances should this stuff end up in 
C> the nwxc library, otherwise NWChem will fail in many places.
C>
      subroutine stubs
#include "errquit.fh"
      entry fmemreq
      entry pdspev
      entry pdspgv
      call errquit("nwxc_fnc_print stub",0,UERR)
      end
