c
      LOGICAL FUNCTION ACOLO_1A_2H(h1b,h2b)
      IMPLICIT NONE
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
      INTEGER h1b
      INTEGER h2b
      INTEGER io1
      acolo_1a_2h = .false.
      io1=int_mb(k_active+h1b-1)+int_mb(k_active+h2b-1) 
      if(io1.ge.numact-2)  acolo_1a_2h = .true.
      RETURN
      END
c $Id$
