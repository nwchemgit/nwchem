c
c     mnsol_interface
c
      subroutine mnsol_interface(rtdb,geom,nat,ixmem,
     &        sola,solb,solc,solg,solh,soln,icds)
c
      implicit none
c
#include "errquit.fh"
#include "inp.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "stdio.fh"
#include "nwc_const.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util_params.fh"
c
      integer rtdb,geom,nat,ixmem,icds
      double precision sola,solb,solc,solg,solh,soln
c
      character*16 atsym(nat)
      double precision c(3,nat)
      double precision q(nat)
      integer atnum(nat)
c
      double precision gcds,areacds
      double precision dcds(3,nat)
      double precision x(ixmem)
c
      integer iat,ixyz
      integer l_c, k_c, l_atsym, k_atsym, l_z, k_z
c
c     preliminaries
      do iat = 1,nat
         do ixyz = 1,3
           c(ixyz,iat) = 0.d0
         end do
         atnum(iat) = 0
         q(iat) = 0
      end do
c
c     get the coordinates
      if(.not.geom_cart_get2(geom,nat,atsym,c,q,atnum))
     &  call errquit('mnsol_interface: geom_cart_get2 failed.',911, 
     &    GEOM_ERR)
c
C     call cavity dispersion solvent (CDS) structure driver
      call cdsset(icds,gcds,areacds,nat,c,atnum,dcds,x,
     &          sola,solb,solc,solg,solh,soln)
c
c     store non-electrostatic energy contribution from SMD model
      gcds = gcds / cau2kcal
      if (.not. rtdb_put(rtdb,'cosmo:gcds',mt_dbl,1,gcds))
     $ call errquit('mnsol_interface: rtdb put failed',911,RTDB_ERR)
c
c     store area from SMD model 
      if (.not. rtdb_put(rtdb,'cosmo:areacds',mt_dbl,1,areacds))
     $ call errquit('mnsol_interface: rtdb put failed',911,RTDB_ERR)
c
c     store gradients from SMD model 
      if (.not. rtdb_put(rtdb, 'cosmo:dcds', mt_dbl,3*nat,dcds))
     $ call errquit('mnsol_interface: rtdb put failed',911,RTDB_ERR)
c
c     print gradients
c     if (ga_nodeid().eq.0) then
c       write (6,*) 'GRADIENTSGRADIENTSGRADIENTSGRADIENTS'
c       do iat = 1,nat
c         write (6,*) iat,dcds(1,iat),dcds(2,iat),dcds(3,iat)
c       enddo
c       write (6,*) 'GRADIENTSGRADIENTSGRADIENTSGRADIENTS'
c     endif
c
      return
      end
c $Id$
