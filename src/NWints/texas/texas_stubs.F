C stub routines for texas integral code.
C these routines are only called by various api routines
C 
C RA Kendall 11/95 PNNL
C
* $Id: texas_stubs.F,v 1.2 1995-11-01 17:34:49 d3e129 Exp $
      subroutine texas_terminate()
      end
      subroutine texas_init()
      end
      subroutine texas_hf2_m()
      call errquit(' use_texas is NOPE. should never get here ',911)
      end
