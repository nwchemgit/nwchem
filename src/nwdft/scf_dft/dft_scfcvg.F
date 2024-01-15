      subroutine dft_scfcvg(rms, derr, etold, etnew, e_conv, 
     &                      d_conv, g_conv, ipol, iter, iterations, 
     &     idone, rtdb, converged, diising,
     B     back_to_atguess)
c
c     $Id$
c
      implicit none
#include "errquit.fh"
c
      double precision rms(2)   ! [input]
      double precision derr(2)  ! [input]
      double precision etold    ! [input]
      double precision etnew    ! [input]
      double precision e_conv   ! [input]
      double precision d_conv   ! [input]
      double precision g_conv   ! [input]
      integer ipol              ! [input]
      integer iter              ! [input]
      integer iterations        ! [input]
      integer idone             ! [output]
      integer rtdb              ! [input]
      logical converged         ! [output]
      logical diising           ! [input]
      integer back_to_atguess   ! [in/out]
c
#include "mafdecls.fh"
#include "rtdb.fh"
c     
      logical e_conv_logical, d_conv_logical, g_conv_logical
      logical ENERGY, DENSITY, GRADIENT
      double precision de, abde
c
      converged = .false.
c
      e_conv_logical = .false.
      d_conv_logical = .false.
      g_conv_logical = .false.
c
      ENERGY = e_conv.gt.0
      DENSITY = d_conv.gt.0
      GRADIENT = g_conv.gt.0
c
      idone = 0
c
c     Evaluate change in energy.
c     
      de = etnew - etold
      etold=etnew
      abde = dabs(de)
c     
c     Check to see if energy is converged.
c
      if (ENERGY)then     
         if (abde.lt.e_conv)e_conv_logical = .true.
      else
         e_conv_logical = .true.
      endif
c     
c     Check for density matrix convergence.
c     
      if (DENSITY)then
         if (dsqrt(rms(1)).le.d_conv)d_conv_logical = .true.
         if (ipol.eq.2)then
           if (dsqrt(rms(2)).le.d_conv) then
             d_conv_logical = d_conv_logical.and..true.
          else
             d_conv_logical = d_conv_logical.and..false.
          endif
         endif
      else
         d_conv_logical = .true.
      endif
c     
c     Check for gradient convergence.
c     
      if (GRADIENT.and.diising)then
         if (derr(1).le.g_conv)g_conv_logical = .true.
         if (ipol.eq.2)then
           if (derr(2).le.g_conv) then
              g_conv_logical = g_conv_logical.and..true.
             else
                g_conv_logical = g_conv_logical.and..false.
             endif
         endif
      else
         g_conv_logical = .true.
      endif
c
c     Check over-all convergence.
c
      converged = e_conv_logical.and.d_conv_logical.and.g_conv_logical
      if (converged)idone = 1
c
c     try at guess if d_e > 1d-1 after back_to_atguess iterations
c
      if(.not.converged.and.back_to_atguess.gt.0) then
         if(abde.gt.1d-1.and.iter.gt.back_to_atguess)
     B        back_to_atguess=-1
      endif
c     
c     Check iteration value.
c     
      if (iter.lt.1)then
         return
      elseif (iter.eq.iterations)then
         idone = 1
      endif
c
c     If all convergence criterion met or number of iterations has been
c     exceeded, write "converged" to RTDB.
c
      if (idone.eq.1)then
         if (.not.rtdb_put(rtdb, 'dft:converged', MT_LOG, 1, converged))
     &      call errquit('dft_scfcvg: rtdb_put failed', 1, RTDB_ERR)
      endif
      return
      end
