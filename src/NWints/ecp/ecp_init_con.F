C $Id: ecp_init_con.F,v 1.2 1996-10-16 18:26:07 mg201 Exp $
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
      return
      end
