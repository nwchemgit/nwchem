c
      LOGICAL FUNCTION ACOLO_2P_1H(p1b,p2b,h1b)
      IMPLICIT NONE
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
      INTEGER h1b
      INTEGER p1b
      INTEGER p2b
      INTEGER iocc_o
      INTEGER iocc_v
      acolo_2p_1h = .false.
      iocc_o=int_mb(k_active+h1b-1)
      iocc_v=int_mb(k_active+p1b-1)+int_mb(k_active+p2b-1)
      if(iocc_o.ge.numact-2.and.iocc_v.ge.numact-1) acolo_2p_1h=.true.
      RETURN
      END
c $Id$
