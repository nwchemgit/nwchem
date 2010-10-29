C $Id$
************************************************************************
*                                                                      *
      subroutine rel_init_con (rtdb)
*                                                                      *
*   Routine to initialize constants required by relativistic integrals *
*                                                                      *
*   Written by Kenneth G. Dyall                                        *
*                                                                      *
************************************************************************
      implicit none
c     
      integer rtdb ! [input] handle to database
c
#include "rtdb.fh"
#include "mafdecls.fh"
#include "rel_consts.fh"
c
      double precision clight
c
      cau = 137.0359895d00
      if (rtdb_get(rtdb,'relativistic:clight',mt_dbl,1,clight))
     &    cau = clight
      csq = cau*cau
      alpha = 1.0d00/cau
      halp = 0.5d00*alpha
      halsq = halp*alpha
      qalsq = halp*halp
      salf = qalsq*qalsq
c
      return
      end
