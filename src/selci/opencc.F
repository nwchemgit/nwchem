      subroutine opencc
#include "errquit.fh"
*
* $Id: opencc.F,v 1.1 2003-10-17 22:29:07 carlfahl Exp $
*
      entry putcc
      entry getcc
      entry getccn
      entry srchcc
      entry closecc
c
      call errquit(' C I/O routines not implemented ',-1, CAPMIS_ERR)
c
      end
