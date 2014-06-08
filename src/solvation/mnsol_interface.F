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
c
      integer rtdb,geom,nat,ixmem,icds
      double precision sola,solb,solc,solg,solh,soln
c
      double precision c(3,nat)
      integer ian(nat)
      double precision gcds,areacds
      double precision dcds(3,nat)
      double precision x(ixmem)
c
      integer iat
      integer l_c, k_c, l_atsym, k_atsym, l_z, k_z
      logical status
c
c     preliminaries
      do iat = 1,nat
         c(1,iat) = 0.d0
         c(2,iat) = 0.d0
         c(3,iat) = 0.d0
         ian(iat) = 0
      end do
c
c     allocate c,z,atsym
      if(.not.ma_push_get(mt_dbl,nat*3,'c',l_c,k_c))
     & call errquit('mnsol_interface: k_c failed',911,MA_ERR)
      if(.not.ma_push_get(mt_dbl,nat,'z',l_z,k_z)) 
     & call errquit('mnsol_interface: k_z failed',911,MA_ERR)
      if(.not.ma_push_get(mt_byte,nat*16,'atsym',l_atsym,k_atsym)) 
     & call errquit('mnsol_interface: k_atsym failed',911,MA_ERR)
c
c     get the coordinates
      if(.not.geom_cart_get(geom,nat,byte_mb(k_atsym),
     &   dbl_mb(k_c),dbl_mb(k_z))) call errquit
     $    (' cosmo_initialize: geom_cart_get failed.',911, GEOM_ERR)
c
      do iat=1,nat
       c(1,iat) = dbl_mb(k_c  +3*(iat-1))
       c(2,iat) = dbl_mb(k_c+1+3*(iat-1))
       c(3,iat) = dbl_mb(k_c+2+3*(iat-1))
       ian(iat) = dbl_mb(k_z+iat-1)
c       if (ga_nodeid().eq.0)
c     &     write(luout,*) iat,ian(iat),c(1,iat),c(2,iat),c(3,iat)
      end do
C
C     call cavity dispersion solvent (CDS) structure driver
c     if (ga_nodeid().eq.0) write(luout,*) "calling cdsset"
      call cdsset(icds,gcds,areacds,nat,c,ian,dcds,x,
     &          sola,solb,solc,solg,solh,soln)
C
C     print out the output from CDSSET
C
c      if (ga_nodeid().eq.0) write (luout,*) "GCDS: ", GCDS
c      if (ga_nodeid().eq.0) write (luout,*) "AREACDS: ", AREACDS
c
c     store non-electrostatic energy contribution from SMD model
      gcds = gcds / 627.509451d+00
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
c     release local memory
      if(.not.ma_pop_stack(l_atsym)) 
     & call errquit('mnsol_interface: k_atsym failed',911,MA_ERR)
      if(.not.ma_pop_stack(l_z)) 
     & call errquit('mnsol_interface: k_z failed',911,MA_ERR)
      if(.not.ma_pop_stack(l_c)) 
     & call errquit('mnsol_interface: k_c failed',911,MA_ERR)
c
      return
      end
c $Id$
