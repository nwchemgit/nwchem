czora...compute the matrix elements for the zora corrections
c
      subroutine zora_getv(rtdb, g_dens, g_zora, nexc)
c
C$Id: zora_getv.F,v 1.2 2007-07-11 01:01:43 niri Exp $
c
      implicit none
c
#include "mafdecls.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "global.fh"
c
      integer g_dens(2), g_zora(2)
      integer rtdb
      integer nexc
      double precision rho_n   !integrated electron number
      integer itmat
      double precision dummy

      call grid_quadv0(rtdb, g_dens, g_zora, nexc,rho_n, dummy, 
     &  dbl_mb(itmat))

      return
      end
