c
c
c
      LOGICAL FUNCTION ACOLO_O(h1b,h2b,h3b)
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
      INTEGER iocc
      acolo_o = .false.
      iocc=int_mb(k_active+h1b-1)+int_mb(k_active+h2b-1)+
     &     int_mb(k_active+h3b-1)
      if(iocc.ge.numact) acolo_o=.true.
      RETURN
      END
c $Id$
