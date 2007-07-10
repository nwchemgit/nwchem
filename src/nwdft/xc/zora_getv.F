czora...compute the matrix elements for the zora corrections
c
      subroutine zora_getv(rtdb, nbf_ao, ipol, 
     &      g_dens, g_zora, nexc)
c
C$Id: zora_getv.F,v 1.1 2007-07-10 21:43:02 niri Exp $
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
      integer nbf_ao, ipol
      double precision rho_n
      integer itmat
      double precision dummy
      integer me

      me=ga_nodeid()
 
      call grid_quadv0(rtdb, g_dens, g_zora, nexc,rho_n, dummy, 
     &  dbl_mb(itmat))

      return
      end
