c     Minimize the electrostatic energy of a collection
c     of point charges on the surface of a unit sphere. 
c
c     The implementation uses a Quasi-Newton method based
c     on the inverse BFGS update. Currently, only the
c     diagonal entries of the inverse hessian are updated.
c     The values of the constants c1 and c2 were taken
c     from the L-BFGS-B code of Zhu, Byrd, Nocedal and Morales.
c
c     Arfken convention (r,\theta,\phi) -> (radial,polar,azimuthal)
c     is used. The radial coordinate is constrained to r=1.

      subroutine cosmo_thomson(n,xyz)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
      integer n
      integer, parameter :: m = 100
      double precision,dimension(3*n)   :: xyz

      double precision energy, r, ri
      double precision energy0,grad0,energy1,grad1,energy2
      double precision dum, alpha, sy, ps

      integer l_scoor,l_scoor2,l_sgrad,l_sgrad2,l_p,l_s,l_ys,l_ay,l_hess
      integer iscoor,iscoor2,isgrad,isgrad2,ip,is,iys,iay,ihess

      integer i,iter,tries,maxiter,maxtries

      integer itheta, iphi
      integer ix, iy, iz

      double precision c1, c2
      double precision,external :: ddot

      parameter (maxiter = 1000, maxtries = 4)
      parameter (c1 = 0.0001d0, c2 = 0.9d0)

      if (.not.ma_push_get(mt_dbl,2*n,'scoor',l_scoor,iscoor))
     $  call errquit('cosmo_thomson: allocation failed',1,MA_ERR) 
      if (.not.ma_push_get(mt_dbl,2*n,'scoor2',l_scoor2,iscoor2))
     $  call errquit('cosmo_thomson: allocation failed',1,MA_ERR) 
      if (.not.ma_push_get(mt_dbl,2*n,'sgrad',l_sgrad,isgrad))
     $  call errquit('cosmo_thomson: allocation failed',1,MA_ERR) 
      if (.not.ma_push_get(mt_dbl,2*n,'sgrad2',l_sgrad2,isgrad2))
     $  call errquit('cosmo_thomson: allocation failed',1,MA_ERR) 
      if (.not.ma_push_get(mt_dbl,2*n,'s',l_s,is))
     $  call errquit('cosmo_thomson: allocation failed',1,MA_ERR) 
      if (.not.ma_push_get(mt_dbl,2*n,'p',l_p,ip))
     $  call errquit('cosmo_thomson: allocation failed',1,MA_ERR) 
      if (.not.ma_push_get(mt_dbl,2*n,'ys',l_ys,iys))
     $  call errquit('cosmo_thomson: allocation failed',1,MA_ERR) 
      if (.not.ma_push_get(mt_dbl,2*n,'ay',l_ay,iay))
     $  call errquit('cosmo_thomson: allocation failed',1,MA_ERR) 
      if (.not.ma_push_get(mt_dbl,2*n,'hess',l_hess,ihess))
     $  call errquit('cosmo_thomson: allocation failed',1,MA_ERR) 

      !Cartesian -> Spherical transformation
      do i=1,n
        ix = 3*(i-1) + 1
        iy = 3*(i-1) + 2
        iz = 3*(i-1) + 3
        itheta = 2*(i-1) + 0
        iphi   = 2*(i-1) + 1
        dbl_mb(iscoor+itheta) = acos(xyz(iz))
        dbl_mb(iscoor+iphi)   = atan2(xyz(iy),xyz(ix))
      enddo

      ! Initialize the elements of the inverse Hessian
      ! (Value was obtained by trial and error)
      if (n.gt.1000) then
        call dfill(2*n,0.001d0,dbl_mb(ihess),1)
      else
        call dfill(2*n,0.01d0,dbl_mb(ihess),1)
      endif

      ! Compute the initial repulsion energy and gradient 
      call thomson_grad(dbl_mb(isgrad),dbl_mb(iscoor),energy,n)

      do iter=0,maxiter

        ! Save energy and norm2 of the gradient
        dum = ddot(2*n,dbl_mb(isgrad),1,dbl_mb(isgrad),1)
        if (iter.eq.0) then
          energy0 = energy
          grad0 = sqrt(dum)
        endif
        energy1 = energy
        grad1 = sqrt(dum)

        ! Exit if norm2 is small
        if (sqrt(dum).lt.1d-3) exit

        ! Quasi-Newton step
        do i=0,2*n-1
          dbl_mb(ip+i) = -dbl_mb(ihess+i)*dbl_mb(isgrad+i)
        enddo
        ps = ddot(2*n,dbl_mb(ip),1,dbl_mb(isgrad),1)

        ! Currently, we do not perform an exact line search.
        ! Instead, we guess the value of alpha and check
        ! if the strong Wolfe conditions are met. If the 
        ! conditions are not met, alpha is divided by 3 and
        ! a new check is performed.
        alpha = 1d0

        do tries=1,maxtries
          ! Update coordinates according to current step size
          call dcopy(2*n,dbl_mb(iscoor),1,dbl_mb(iscoor2),1)
          call daxpy(2*n,alpha,dbl_mb(ip),1,dbl_mb(iscoor2),1)

          ! Compute new energy and gradient
          call thomson_grad(dbl_mb(isgrad2),dbl_mb(iscoor2),energy2,n)

          ! Check strong Wolfe conditions
          if ( (energy2.le.energy+c1*alpha*ps) .and.
     $         (abs(ddot(2*n,dbl_mb(ip),1,dbl_mb(isgrad2),1))
     $         .le.abs(c2*ps) )) then
            exit

          ! Scale alpha and try again
          elseif (tries.lt.maxtries) then
            alpha = alpha*0.3d0

          ! If energy is growing, restart inverse Hessian diagonal
          elseif (energy2.gt.energy) then
            if (n.gt.1000) then
              call dfill(2*n,0.002d0,dbl_mb(ihess),1)
            else
              call dfill(2*n,0.0002d0,dbl_mb(ihess),1)
            endif
          endif
        enddo

        ! BFGS update on the diagonal
        call dcopy(2*n,dbl_mb(ip),1,dbl_mb(is),1)
        call dscal(2*n,alpha,dbl_mb(is),1)

        alpha = 0d0
        do i=0,2*n-1
          dbl_mb(iys+i) = dbl_mb(isgrad2+i) - dbl_mb(isgrad+i)
          dbl_mb(iay+i) = dbl_mb(ihess+i)*dbl_mb(iys+i)
          alpha = alpha + dbl_mb(iay+i)*dbl_mb(iys+i)
        enddo
        sy = ddot(2*n,dbl_mb(is),1,dbl_mb(iys),1)

        ! Don't update if denominator gets too small
        if (sy.gt.1d-6) then
          alpha = (sy+alpha)/sy**2
          do i=0,2*n-1
            dbl_mb(ihess+i) = dbl_mb(ihess+i) + 
     $      alpha*dbl_mb(is+i)**2 - 2d0/sy*dbl_mb(iay+i)*dbl_mb(is+i)
          enddo
        endif

        ! Update coordinates, gradients and energy
        call dcopy(2*n,dbl_mb(iscoor2),1,dbl_mb(iscoor),1)
        call dcopy(2*n,dbl_mb(isgrad2),1,dbl_mb(isgrad),1)
        energy = energy2
      
      enddo

      ! Transform Spherical -> Cartesian coordinates
      do i=1,n
        ix = 3*(i-1) + 1
        iy = 3*(i-1) + 2
        iz = 3*(i-1) + 3
        itheta = 2*(i-1) + 0
        iphi   = 2*(i-1) + 1
        xyz(ix) = sin(dbl_mb(iscoor+itheta))*cos(dbl_mb(iscoor+iphi))
        xyz(iy) = sin(dbl_mb(iscoor+itheta))*sin(dbl_mb(iscoor+iphi))
        xyz(iz) = cos(dbl_mb(iscoor+itheta))
      enddo

      ! Print Minimization Summary
      if (ga_nodeid().eq.0) then
        write(*,*)
        write(*,*) "  Solving Thomson problem locally "
        write(*,'("    Initial Energy: ",F22.8)') energy0
        write(*,'("  Initial Gradient: ",F22.8)') grad0
        write(*,'("      Final Energy: ",F22.8)') energy1
        write(*,'("    Final Gradient: ",F22.8)') grad1
      endif

      ! Make sure all the processes have the same set of
      ! Cartesian coordiantes
      call ga_brdcst(0,xyz,ma_sizeof(mt_dbl,3*n,mt_byte),0)

      if (.not.ma_chop_stack(l_scoor))
     $  call errquit('cosmo_thomson: could not chop stack',0,MA_ERR) 

      end


c     Computes the electrostatic repulsion energy and
c     gradient for a collection a charges on the surface
c     of a unit sphere. Spherical coordinates are used.

      subroutine thomson_grad(grad,coor,energy,n)
      implicit none
#include "global.fh"
      integer n
      double precision energy
      double precision,dimension(2*n) :: grad, coor

      integer i, j
      integer itheta, iphi
      integer jtheta, jphi

      double precision r, ri, ri3
      double precision theta1,theta2,phi1,phi2
      double precision sin1,cos1,sin2,cos2,sin12,cos12

      grad(:) = 0d0
      energy = 0d0
      do i=1+ga_nodeid(),n,ga_nnodes()
        itheta = 2*(i-1) + 1
        iphi   = 2*(i-1) + 2

        theta1 = coor(itheta)
        phi1   = coor(iphi)

        sin1 = sin(theta1)
        cos1 = cos(theta1)
        do j=1,i-1
          jtheta = 2*(j-1) + 1
          jphi   = 2*(j-1) + 2

          theta2 = coor(jtheta)
          phi2   = coor(jphi)

          sin2 = sin(theta2)
          cos2 = cos(theta2)
          cos12 = cos(phi1-phi2)
          sin12 = sin(phi1-phi2)

          ! Distance between two points on the
          ! surface of the unit sphere
          r = sqrt(2d0 - 2*(sin1*sin2*cos12+cos1*cos2))
          ri = 1d0/r
          ri3 = ri**3

          energy = energy + ri

          grad(itheta) = grad(itheta) - 
     $               (sin1*cos2 - cos1*sin2*cos12)*ri3
          grad(jtheta) = grad(jtheta) -
     $               (cos1*sin2 - sin1*cos2*cos12)*ri3     

          ! Azimuthal gradient has a 1/sin(theta) factor
          grad(iphi) = grad(iphi) - sin2*sin12*ri3
          grad(jphi) = grad(jphi) + sin1*sin12*ri3
        enddo
      enddo

      call ga_dgop(0,energy,1,'+')
      call ga_dgop(0,grad(1),2*n,'+')
      call dscal(2*n,0.5d0,grad,1)

      end subroutine
      
