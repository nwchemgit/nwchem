c
      LOGICAL FUNCTION ACOLO_1P_1H(p1b,h1b)
      IMPLICIT NONE
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
      INTEGER h1b
      INTEGER p1b
      INTEGER iocc_o
      INTEGER iocc_v
      acolo_1p_1h = .false.
      iocc_o=int_mb(k_active+h1b-1)
      iocc_v=int_mb(k_active+p1b-1)
      if(iocc_o.ge.numact-2.and.iocc_v.ge.numact-2) acolo_1p_1h=.true.
      RETURN
      END
c $Id$
