      Subroutine grid_qpts(xyz, ictr, irsh, iang, aqbuf, nq, 
     *     rad)
c
c$Id: grid_qpts.F,v 1.5 2000-03-07 02:02:32 edo Exp $
c
      implicit none
c      
      integer iang ! [input]
      integer ictr ! [input]
      integer irsh ! [input]
      integer nq   ! [input]
c      
#include "cdft.fh"
#include "util.fh"
#include "global.fh"
c
c---> Cartesian Coordinates of Integration Center
c
      double precision xyz(3) ! [input]
c
c---> Cartesian Coordinates of Quadrature Points & Weights
c
      logical  oprint_grid
      double precision aqbuf(4,nq), qwght_nq,rad
c
c     Compute the quadrature points & weights for a given 
c     set of radial shells.
c     compute rad
c
      integer k, l, ia_ictr
      double precision r, tx, ty, tz, xor, yor, zor
c
      oprint_grid = util_print('griddebug', print_debug)
c
      nq = 0
c
      ia_ictr = iatype(ictr)
c
      r = rpts(irsh,ia_ictr)
      rad=r
      
c
c     debug writes
c
      if (ga_nodeid().eq.0.and.oprint_grid)then
         write(6,*)' xyz: ', xyz
         write(6,*)' ictr, irsh, iang, ia_ictr, rad: ',
     &               ictr, irsh, iang, ia_ictr, rad
      endif
c
      if (leb.or.daul) then
c         
         do k = 1, ntheta(iang)
c           
            tx=za_leb(1,k,iang)
            ty=za_leb(2,k,iang)
            tz=za_leb(3,k,iang)
c
c           Rotate the directional cosines.
c
            xor = qrot(1,1)*tx + qrot(1,2)*ty + qrot(1,3)*tz
            yor = qrot(2,1)*tx + qrot(2,2)*ty + qrot(2,3)*tz
            zor = qrot(3,1)*tx + qrot(3,2)*ty + qrot(3,3)*tz
c
c           Define the quadrature point and weight.
c            
            nq = nq + 1
            aqbuf(1,nq) = xyz(1) + r*xor
            aqbuf(2,nq) = xyz(2) + r*yor
            aqbuf(3,nq) = xyz(3) + r*zor
c             
            aqbuf(4,nq) = wa_leb(k,iang)*rwghts(irsh,ia_ictr)
         enddo
      else
c
c        Loop over the angular quadrature points.
c
         do 20 k = 1, ntheta(iang)
c
            qwght_nq = awghts(k,iang)*rwghts(irsh,ia_ictr)
c
            do 10 l = 1, nphi(iang)
c
c              Define the directional cosines.
c
               tx = sintheta(k,iang)*cosphi(l,iang)
               ty = sintheta(k,iang)*sinphi(l,iang)
               tz = costheta(k,iang)
c
c              Rotate the directional cosines.
c
               xor = qrot(1,1)*tx + qrot(1,2)*ty + qrot(1,3)*tz
               yor = qrot(2,1)*tx + qrot(2,2)*ty + qrot(2,3)*tz
               zor = qrot(3,1)*tx + qrot(3,2)*ty + qrot(3,3)*tz
c
c              Define the quadrature point and weight.
c             
               nq = nq + 1
               aqbuf(1,nq) = xyz(1) + r*xor
               aqbuf(2,nq) = xyz(2) + r*yor
               aqbuf(3,nq) = xyz(3) + r*zor
c             
               aqbuf(4,nq) = qwght_nq
c
   10       continue
c
   20    continue
c
      endif
c
      return
      end
