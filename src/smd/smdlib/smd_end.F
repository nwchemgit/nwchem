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
      character*32 pname
c
      pname = "smd_end"
c
c      write(*,*) "in "//pname
c
      call smd_data_destroy_all()
c
      return
      end

