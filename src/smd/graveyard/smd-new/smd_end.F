c $Id$
      subroutine smd_end()
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      integer srtdb
c
      character*32 pname
c
      pname = "smd_end"
      call smd_atom_end()
c
      return
      end

