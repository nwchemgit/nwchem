      logical function int_chk_init(msg)
c $Id: int_chk_init.F,v 1.1 1994-04-04 20:27:34 d3e129 Exp $
      implicit none
#include "apiP.fh"
c:: passed
      character*(*) msg
c
c check initialization
c
      int_chk_init = init_int.ne.0
      if (.not.int_chk_init) then
        write(6,*)' ERROR:: <',msg,'>'
        write(6,*)' int_init was not called '
      endif
      end
