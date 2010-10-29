czora...compute the matrix elements for the zora corrections
c
      subroutine zora_getv(rtdb, g_dens, g_zora, nexc)
c
C$Id$
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
      double precision tmat
      double precision dummy(2)
c
      dummy(1) = 0.d0
      dummy(2) = 0.d0
c
      call grid_quadv0(rtdb, g_dens, g_zora, nexc,rho_n, dummy, tmat)
c
      return
      end
