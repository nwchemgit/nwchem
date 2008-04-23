c $Id: smd_end.F,v 1.2 2008-04-23 16:48:29 marat Exp $
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

