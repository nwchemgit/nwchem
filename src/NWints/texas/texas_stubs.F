C stub routines for texas integral code.
C these routines are only called by various api routines
C 
C RA Kendall 11/95 PNNL
C
* $Id$
#define USE_TEXAS
#if defined(USE_TEXAS)
      subroutine texas_using_it4pnl()
      end
#else
      subroutine texas_terminate()
      end
      subroutine texas_init()
      end
      subroutine texas_set_accy()
      end
      subroutine texas_hf2_m()
#include "errquit.fh"
      call errquit(' use_texas is NOPE. should never get here ',911,
     &       INT_ERR)
      end
#endif
