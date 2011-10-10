c
      LOGICAL FUNCTION ACOLO_2H(h1b,h2b)
      IMPLICIT NONE
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
      INTEGER h1b
      INTEGER h2b
      INTEGER iocc_o
      acolo_2h = .false.
      iocc_o=int_mb(k_active+h1b-1)+int_mb(k_active+h2b-1)
      if(iocc_o.ge.numact-1) acolo_2h=.true.
      RETURN
      END
c $Id$
