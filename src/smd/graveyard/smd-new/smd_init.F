      subroutine smd_init(operation,rtdb)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      character*(*) operation
      integer rtdb
c
      character*80 pname
      character*80 tag

      pname = "smd_init"
      tag = operation
c
c     init core components
c     -------------------
      call smd_srtdb_init()
      call smd_rdinp()
      call smd_atom_init()
      call smd_param_init()
c      call smd_message_flush(6)
      return
      end

c $Id$
