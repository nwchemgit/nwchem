      subroutine cosmo_thomson(n,xyz)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
      integer n
      integer, parameter :: m = 100
      double precision,dimension(3*n)   :: xyz, cgrad, cgrad2, xyz2
      double precision,dimension(2*n)   :: scoor, sgrad, sgrad2, scoor2

      double precision,dimension(2*n) :: p, s, ys, ay
      double precision,dimension(2*n) :: hess

      double precision dxyz(3)
      double precision energy, r, ri
      double precision energy0,grad0,energy1,grad1,energy2
      double precision dum,alpha,sy,sinp,cosp,sint,cost
      double precision mindist, step

      integer i,j,k,l,tries

      integer ixyz,jxyz
      integer itheta, iphi
      integer x, y ,z
      integer ix, iy, iz
      integer jx, jy, jz
      double precision pi

      pi = acos(-1d0)

      x = 1
      y = 2
      z = 3

      hess(:) = 0.01d0

      energy = 0d0
      cgrad(:) = 0d0
      mindist = 2d0
      do i=1+ga_nodeid(),n,ga_nnodes()
        ix = 3*(i-1) + 1
        iy = 3*(i-1) + 2
        iz = 3*(i-1) + 3
        do j=1,i-1
          if (i.eq.j) cycle
          jx = 3*(j-1) + 1
          jy = 3*(j-1) + 2
          jz = 3*(j-1) + 3
          dxyz(x:z) = xyz(ix:iz) - xyz(jx:jz)
          r = sqrt(sum(dxyz**2))
          mindist = min(mindist,r)
          ri = 1d0/r
          energy = energy + ri
          cgrad(ix:iz) = cgrad(ix:iz) - 0.5d0*dxyz(x:z)*ri**3
          cgrad(jx:jz) = cgrad(jx:jz) + 0.5d0*dxyz(x:z)*ri**3
        enddo
      enddo
      call dgop(mt_dbl,energy,1,'+')
      call dgop(mt_dbl,cgrad,3*n,'+')

      do i=1,n
        ix = 3*(i-1) + 1
        iy = 3*(i-1) + 2
        iz = 3*(i-1) + 3
        itheta = 2*(i-1) + 1
        iphi   = 2*(i-1) + 2

        scoor(itheta)  = acos(xyz(iz))
        scoor(iphi)    = atan2(xyz(iy),xyz(ix))

        dxyz(1) = xyz(iz)*cos(scoor(iphi))
        dxyz(2) = xyz(iz)*sin(scoor(iphi))
        dxyz(3) = -sqrt(1d0-xyz(iz)**2)

        sgrad(itheta) = dot_product(dxyz,cgrad(ix:iz))
        sgrad(iphi) =  -cgrad(ix)*sin(scoor(iphi))
     $                 +cgrad(iy)*cos(scoor(iphi))
      enddo

      l = 0
      do

      dum = dot_product(sgrad,sgrad)
      if (l.eq.0) then
        energy0 = energy
        grad0 = sqrt(dum)
      endif
      energy1 = energy
      grad1 = sqrt(dum)

      if (sqrt(dum).lt.1d-3) exit

      alpha = 1d0
      p(:) = -hess(:)*sgrad(:)

      tries = 0

  100 continue
      tries = tries + 1
      do i=1,n
        sint = sin(scoor(2*(i-1)+1)+alpha*p(2*(i-1)+1))
        cost = cos(scoor(2*(i-1)+1)+alpha*p(2*(i-1)+1))
        sinp = sin(scoor(2*(i-1)+2)+alpha*p(2*(i-1)+2))
        cosp = cos(scoor(2*(i-1)+2)+alpha*p(2*(i-1)+2))
        xyz2(3*(i-1)+1) = sint*cosp
        xyz2(3*(i-1)+2) = sint*sinp
        xyz2(3*(i-1)+3) = cost
      enddo

      energy2 = 0d0
      cgrad2(:) = 0d0
      do i=1+ga_nodeid(),n,ga_nnodes()
        ix = 3*(i-1) + 1
        iy = 3*(i-1) + 2
        iz = 3*(i-1) + 3
        do j=1,i-1
          if (i.eq.j) cycle
          jx = 3*(j-1) + 1
          jy = 3*(j-1) + 2
          jz = 3*(j-1) + 3
          dxyz(x:z) = xyz2(ix:iz) - xyz2(jx:jz)
          r = sqrt(sum(dxyz**2))
          ri = 1d0/r
          energy2 = energy2 + ri
          cgrad2(ix:iz) = cgrad2(ix:iz) - 0.5d0*dxyz(x:z)*ri**3
          cgrad2(jx:jz) = cgrad2(jx:jz) + 0.5d0*dxyz(x:z)*ri**3
        enddo
      enddo
      call dgop(mt_dbl,energy2,1,'+')
      call dgop(mt_dbl,cgrad2,3*n,'+')
      do i=1,n
        ix = 3*(i-1) + 1
        iy = 3*(i-1) + 2
        iz = 3*(i-1) + 3
        itheta = 2*(i-1) + 1
        iphi   = 2*(i-1) + 2

        sinp = sin(atan2(xyz2(iy),xyz2(ix)))
        cosp = cos(atan2(xyz2(iy),xyz2(ix)))

        dxyz(1) = xyz2(iz)*cosp
        dxyz(2) = xyz2(iz)*sinp
        dxyz(3) = -sqrt(1d0-xyz2(iz)**2)

        sgrad2(itheta) = dot_product(dxyz,cgrad2(ix:iz))
        sgrad2(iphi) =  -cgrad2(ix)*sinp + cgrad2(iy)*cosp
      enddo

      if (tries.lt.4) then
      if ((energy2.le.energy+0.0001d0*alpha*dot_product(p,sgrad)).and.
     $ (abs(dot_product(p,sgrad2)).le.abs(0.9d0*dot_product(p,sgrad))))
     $then
        alpha = alpha
      else
        alpha = alpha*0.3d0
        goto 100
      endif
      endif

      if (tries.eq.4 .and. energy2.gt.energy) then
        hess = 0.002d0
      endif

      s(:) = alpha*p(:)
      ys(:) = sgrad2(:) - sgrad(:)
      sy = dot_product(s,ys)
      ay = hess(:)*ys(:)

      if (sy.gt.1d-6) then
        alpha = (sy+dot_product(ys,ay))/sy**2
        hess = hess + alpha*s(:)*s(:) - 2d0/sy * ay(:)*s(:)
      endif

      scoor(:) = scoor(:) + s(:)
      sgrad(:) = sgrad2(:)
      xyz(:) = xyz2(:)
      energy = energy2
      
      l = l + 1
      if (l.eq.1000) exit
      enddo

      if (ga_nodeid().eq.0) then
        write(*,*)
        write(*,*) "  Solving Thomson problem locally "
        write(*,'("    Initial Energy: ",F20.12)') energy0
        write(*,'("  Initial Gradient: ",F20.12)') grad0
        write(*,'("      Final Energy: ",F20.12)') energy1
        write(*,'("    Final Gradient: ",F20.12)') grad1
      endif

      call ga_brdcst(0,xyz,ma_sizeof(mt_dbl,3*n,mt_byte),0)

      end


      
