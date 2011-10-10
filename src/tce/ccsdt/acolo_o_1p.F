c
c
c
      LOGICAL FUNCTION ACOLO_O_1P(p1b,h1b,h2b,h3b)
      IMPLICIT NONE
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
      INTEGER h1b
      INTEGER h2b
      INTEGER h3b
      INTEGER p1b
      INTEGER iocc_o
      INTEGER iocc_v
      acolo_o_1p = .false.
      iocc_o=int_mb(k_active+h1b-1)+int_mb(k_active+h2b-1)+
     &     int_mb(k_active+h3b-1)
      iocc_v=int_mb(k_active+p1b-1)
      if(iocc_o.ge.numact.and.iocc_v.ge.numact-2) acolo_o_1p=.true.
      RETURN
      END
c $Id$
