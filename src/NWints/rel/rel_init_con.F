C $Id: rel_init_con.F,v 1.3 1998-10-20 20:16:31 mg201 Exp $
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
      csq = cau*cau
      alpha = 1.0d00/cau
      halp = 0.5d00*alpha
      halsq = halp*alpha
      qalsq = halp*halp
      salf = qalsq*qalsq
*
      return
      end
