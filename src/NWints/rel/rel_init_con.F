C $Id: rel_init_con.F,v 1.1 1998-09-22 23:07:52 mg201 Exp $
************************************************************************
*                                                                      *
      subroutine rel_init_con
*                                                                      *
*   Routine to initialize constants required by relativistic integrals *
*                                                                      *
*   Written by Kenneth G. Dyall                                        *
*                                                                      *
************************************************************************
      implicit none
#include "rel_consts.fh"
*
      cau = 137.0359895d00
      alpha = 1.0d00/cau
      halp = 0.5d00*alpha
      halsq = halp*alpha
      qalsq = half*halsq
      salf = qalsq*qalsq
*
      return
      end
