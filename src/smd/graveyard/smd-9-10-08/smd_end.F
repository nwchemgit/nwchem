c $Id: smd_end.F,v 1.1 2008-10-01 22:32:29 marat Exp $
      subroutine smd_end(rtdb)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
c     
      integer rtdb
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

