C $Id: ecp_init_con.F,v 1.3 2000-04-13 20:37:27 mg201 Exp $
************************************************************************
*                                                                      *
      subroutine ecp_init_con
*                                                                      *
*   Routine to initialize constants required by ECP code               *
*                                                                      *
*   Written by Kenneth G. Dyall                                        *
*                                                                      *
************************************************************************
      implicit none
#include "ecp_consts.fh"
*
      zero = 0.0d00
      half = 0.5d00
      one = 1.0d00
      two = 2.0d00
      four = 4.0d00
      pi = atan(one)*four
      qpirt = sqrt(pi)/four
*
      ln_thr_ecp = -50.0d0
      thr_ecp = exp(ln_thr_ecp)
*
      return
      end
