      subroutine do_riscf(doit)
C$Id: do_riscf.F,v 1.1 1995-12-29 02:04:10 vg038 Exp $
      implicit none
#include "cfock.fh"
      logical doit
      
      rifock = doit

      return
      end
