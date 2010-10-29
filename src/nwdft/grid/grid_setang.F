      Subroutine grid_setang(rtdb)
c
c$Id$
c
      implicit none
#include "errquit.fh"
c      
#include "rtdb.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "cdft.fh"
c
      integer rtdb
      logical oprint
      integer l_ang, iang, n, mxnang
      double precision phi_inc, wphi, phi
      double precision pi
c
      pi=acos(-1d0)
c
c     Set print options.
c
      oprint = util_print('quadrature', print_high)
      if (.not. rtdb_get(rtdb, 'dft:mxnang', mt_int, 1, mxnang))
     &   call errquit('grid_setang: rtdb_put failed', 119, RTDB_ERR)
c
c     Define the Lebedev or Gauss-Legendre angular quadrature points and weights.
c

      if (ga_nodeid().ne.0) oprint = .false.
      if (leb) then ! lebedev
         call xc_lebedev(l_ang, oprint)
         nang_leb_pts = l_ang
      else ! gauss-legendre
c
c        Define the angular quadrature points and weights for mxnang possibilities.
c
         do 230 iang = 1, mxnang
c
c           Define the number of quadrature points in each dimension.
c
            ntheta(iang) = iang
c
            nphi(iang) = 2*ntheta(iang)
c
c           Define the Gauss-Legendre points and weights for the theta integration.
c
            call grid_gaussleg(ntheta(iang), costheta(1,iang), 
     &                    awghts(1,iang))
c
            do 200 n = 1, ntheta(iang)
               sintheta(n,iang) = sqrt(1.D0 -
     .            costheta(n,iang)*costheta(n,iang))
  200       continue
c
c           Define the Gauss-Chebyshev points and weights for the phi integration.
c
            phi_inc = (2.D0*PI)/nphi(iang)
            wphi = phi_inc
c
            do 210 n = 1, nphi(iang)
               phi = (dble(n) - 0.5D0)*phi_inc
               cosphi(n,iang) = cos(phi)
               sinphi(n,iang) = sin(phi)
  210       continue
c
c           Combine the quadrature weights for the theta and phi integrations into
c           one angular weight.
c
            do 220 n = 1,  ntheta(iang)
               awghts(n,iang) = awghts(n,iang)*wphi
  220       continue
  230    continue
      endif
      return
      end
