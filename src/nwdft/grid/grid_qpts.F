      Subroutine grid_qpts(r,xyz, ictr, irsh, iang, aqbuf, nq, 
     *     fnxyz,nxyz)
c
c$Id$
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
#include "stdio.fh"
      double precision xyz(3) ! [input] quadr center
      logical  oprint_grid
      double precision aqbuf(4,*), qwght_nq ! quadr pts and wght
      integer fnxyz      ! [in] cube def
      integer nx,ny,nz      !cube def
      double precision r ! [in]
      integer nxyz
      double precision nxyz2
      logical grid_incube
      external grid_incube
c
c     Compute the quadrature points & weights for a given 
c     set of radial shells.
c     Keep only point belonging to cube
c
      integer k, l, ia_ictr
      double precision  tx, ty, tz, xor, yor, zor
c
      oprint_grid = util_print('griddebug', print_debug)
c
      nq = 0
c
      call grid_n2nxyz(nxyz,fnxyz,nx,ny,nz)
      nxyz2=nxyz*0.5d0
      ia_ictr = iatype(ictr)
c
c     debug writes
c
      if (ga_nodeid().eq.0.and.oprint_grid)then
         write(LuOut,*)' xyz: ', xyz
         write(LuOut,*)' ictr, irsh, iang, ia_ictr: ',
     &               ictr, irsh, iang, ia_ictr
      endif
c
      if (leb) then
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
            xor=r*xor
            yor=r*yor
            zor=r*zor
c
c           Define the quadrature point and weight.
c            
            if(fnxyz.eq.-99.or.
     .           grid_incube(xor,yor,zor,nx,ny,nz,fnxyz,nxyz2,
     ,              lcube(ia_ictr))) then
               nq = nq + 1
               
               aqbuf(1,nq) = xyz(1) + xor
               aqbuf(2,nq) = xyz(2) + yor
               aqbuf(3,nq) = xyz(3) + zor
c             
               aqbuf(4,nq) = wa_leb(k,iang)*rwghts(irsh,ia_ictr)
            endif
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
               xor=r*xor
               yor=r*yor
               zor=r*zor
c
c              Define the quadrature point and weight.
c             
            if(fnxyz.eq.-99.or.
     .              grid_incube(xor,yor,zor,nx,ny,nz,fnxyz,nxyz2,
     ,              lcube(ia_ictr))) then
               nq = nq + 1
               aqbuf(1,nq) = xyz(1) + xor
               aqbuf(2,nq) = xyz(2) + yor
               aqbuf(3,nq) = xyz(3) + zor
c             
               aqbuf(4,nq) = qwght_nq
            endif
c
   10       continue
c
   20    continue
c
      endif
c
      return
      end
