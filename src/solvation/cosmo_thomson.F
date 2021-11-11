      subroutine cosmo_thomson(n,xyz)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
      integer n
      integer, parameter :: m = 100
      double precision,dimension(3*n)   :: xyz
      double precision,dimension(2*n)   :: scoor, sgrad, sgrad2, scoor2
      double precision,dimension(2*n) :: p, s, ys, ay
      double precision,dimension(2*n) :: hess

      double precision energy, r, ri
      double precision energy0,grad0,energy1,grad1,energy2
      double precision dum,alpha,sy,sinp,cosp,sint,cost
      double precision mindist, step

      integer i,j,k,l,tries

      integer ixyz,jxyz
      integer itheta, iphi, jtheta, jphi
      integer x, y ,z
      integer ix, iy, iz
      integer jx, jy, jz
      double precision pi

      double precision sin1,sin2,sin12,cos1,cos2,cos12
      double precision theta2,phi2,theta1,phi1

      pi = acos(-1d0)

      do i=1,n
        ix = 3*(i-1) + 1
        iy = 3*(i-1) + 2
        iz = 3*(i-1) + 3
        itheta = 2*(i-1) + 1
        iphi   = 2*(i-1) + 2
        scoor(itheta)  = acos(xyz(iz))
        scoor(iphi)    = atan2(xyz(iy),xyz(ix))
      enddo

      hess(:) = 0.01d0
      call thomson_grad(sgrad,scoor,energy,n)

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
      scoor2(:) = scoor(:) + alpha*p(:)
      call thomson_grad(sgrad2,scoor2,energy2,n)

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

      scoor(:) = scoor2(:)
      sgrad(:) = sgrad2(:)
      energy = energy2
      
      l = l + 1
      if (l.eq.1000) exit
      enddo

      do i=1,n
        sint = sin(scoor(2*(i-1)+1))
        cost = cos(scoor(2*(i-1)+1))
        sinp = sin(scoor(2*(i-1)+2))
        cosp = cos(scoor(2*(i-1)+2))
        xyz(3*(i-1)+1) = sint*cosp
        xyz(3*(i-1)+2) = sint*sinp
        xyz(3*(i-1)+3) = cost
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

      subroutine thomson_grad(grad,coor,energy,n)
      implicit none
#include "global.fh"
      integer n
      double precision energy
      double precision,dimension(2*n) :: grad, coor

      integer i, j
      integer ix, iy, iz
      integer jx, jy, jz
      integer itheta, iphi
      integer jtheta, jphi

      double precision r, ri
      double precision theta1,theta2,phi1,phi2
      double precision sin1,cos1,sin2,cos2,sin12,cos12

      grad(:) = 0d0
      energy = 0d0
      do i=1+ga_nodeid(),n,ga_nnodes()
        ix = 3*(i-1) + 1
        iy = 3*(i-1) + 2
        iz = 3*(i-1) + 3
        itheta = 2*(i-1) + 1
        iphi   = 2*(i-1) + 2

        theta1 = coor(itheta)
        phi1   = coor(iphi)

        sin1 = sin(theta1)
        cos1 = cos(theta1)
        do j=1,i-1
          if (j.eq.i) cycle
          jx = 3*(j-1) + 1
          jy = 3*(j-1) + 2
          jz = 3*(j-1) + 3
          jtheta = 2*(j-1) + 1
          jphi   = 2*(j-1) + 2

          theta2 = coor(jtheta)
          phi2   = coor(jphi)

          sin2 = sin(theta2)
          cos2 = cos(theta2)
          cos12 = cos(phi1-phi2)
          sin12 = sin(phi1-phi2)

          r = sqrt(2 - 2*(sin1*sin2*cos12+cos1*cos2))
          ri = 1d0/r
          energy = energy + ri

          grad(itheta) = grad(itheta) - 
     $               0.5d0*(sin1*cos2 - cos1*sin2*cos12)*ri**3
          grad(jtheta) = grad(jtheta) -
     $               0.5d0*(cos1*sin2 - sin1*cos2*cos12)*ri**3     
          grad(iphi) = grad(iphi) - 0.5d0*sin2*sin12*ri**3
          grad(jphi) = grad(jphi) + 0.5d0*sin1*sin12*ri**3
        enddo
      enddo

      call dgop(0,energy,1,'+')
      call dgop(0,grad,2*n,'+')

      end subroutine
      
