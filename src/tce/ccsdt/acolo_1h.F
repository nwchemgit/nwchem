c
      LOGICAL FUNCTION ACOLO_1H(h1b)
      IMPLICIT NONE
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
      INTEGER h1b
      INTEGER iocc_o
      acolo_1h = .false.
      iocc_o=int_mb(k_active+h1b-1)
      if(iocc_o.ge.numact-2) acolo_1h=.true.
      RETURN
      END
c $Id$
