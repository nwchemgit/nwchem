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
#ifdef SCALAPACK
      entry numroc
      entry blacs_gridexit
      entry blacs_gridinfo
      entry blacs_gridinit
      entry descinit
      entry iceil
      entry indxg2p
      entry pdgetrf
      entry pdgetrs
      entry pdlamch
      entry pdpotrf
      entry pdpotri
      entry pdpotrs
      entry pdsyev
      entry pdsyevd
      entry pdsyevr
      entry pdsyevx
      entry pdsygvx
      entry pzheevd
      entry pzheevr
#endif
      call errquit("nwxc_fnc_print: stub",0,0)
      end

