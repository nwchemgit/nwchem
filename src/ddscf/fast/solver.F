      subroutine solver(g, nx, ny, nz, poles, depth, lmax, 
     $     xlo, xhi, ylo, yhi, zlo, zhi)
      implicit none
#include "mafdecls.fh"
      integer nx, ny, nz        ! [input] Grid dimension
      integer depth, lmax       ! [input] FMM parameters
      double precision g(nx,ny,nz)
      double precision 
     $     poles(1:(lmax+1)**2,0:2**depth-1,0:2**depth-1,0:2**depth-1)
      double precision xlo, xhi, ylo, yhi, zlo, zhi
c
c     Given as input a function that is zero at the boundary
c     tabulated on a grid (zero boundary not included) and
c     a set of multipole moments in the volume solve the
c     free space Poisson equation
c     
      integer i, j, k
      double precision pi, kx, ky, kz, scale
      integer lenwork, lenworkfmm, l_work, k_work
      double precision poles0(-25:25,0:25)
      real*4 etime, tt(2), start
c
      pi = 4.0d0*atan(1.0d0)

      start = etime(tt)
      call vsint3d(nx,ny,nz,g)
      write(6,*) ' First FFT ', etime(tt)-start
      start = etime(tt)
      if (.not. ma_verify_allocator_stuff()) call errquit('fred',3)
c     
      scale = 4.0d0*pi*0.125d0/(dble(nx+1)*dble(ny+1)*dble(nz+1))
c
      do k = 1, nz
         kz = k*pi/(zhi-zlo)
         do j = 1, ny
            ky = j*pi/(yhi-ylo)
            do i = 1, nx
               kx = i*pi/(xhi-xlo)
               g(i,j,k) = scale*g(i,j,k)/(kx*kx+ky*ky+kz*kz)
            end do
         end do
      end do
      if (.not. ma_verify_allocator_stuff()) call errquit('fred',4)
c
      write(6,*) ' Scaling ', etime(tt)-start
      start = etime(tt)
c
      call vsint3d(nx,ny,nz,g)
      write(6,*) ' Second FFT ', etime(tt) - start
      start = etime(tt)
      if (.not. ma_verify_allocator_stuff()) call errquit('fred',5)
c
      lenwork = (max(nx,ny,nz)+2)**2
      lenworkfmm = (lmax+1)**2*2*(8**(depth+1)-1)/7
      lenwork = max(lenwork,lenworkfmm)
      if (.not. ma_push_get(mt_dbl, lenwork, 'solver: work',
     $     l_work, k_work)) call errquit('solver: work ', lenwork)
c
      call boundary_layer_multipoles(g,nx,ny,nz,
     $     xlo,xhi,ylo,yhi,zlo,zhi,poles,lmax,depth,dbl_mb(k_work),
     $     'poles')
      write(6,*) ' Boundary poles ', etime(tt) - start
      start = etime(tt)
      if (.not. ma_verify_allocator_stuff()) call errquit('fred',6)
c
      call fmm(lmax, depth, xhi-xlo, yhi-ylo, zhi-zlo, poles,
     $     poles0, dbl_mb(k_work), lenwork)
      write(6,*) ' FMM ', etime(tt) - start
      start = etime(tt)
      if (.not. ma_verify_allocator_stuff()) call errquit('fred',7)
c
      call boundary_layer_multipoles(g,nx,ny,nz,
     $     xlo,xhi,ylo,yhi,zlo,zhi,poles,lmax,depth,dbl_mb(k_work),
     $     'near')
      write(6,*) ' Boundary potentials ', etime(tt) - start
      start = etime(tt)
      if (.not. ma_verify_allocator_stuff()) call errquit('fred',8)
c
      if (.not. ma_pop_stack(l_work)) call errquit
     $     ('solver: ma corrupt?',0)
c
      end
      subroutine boundary_layer_multipoles(g,nx,ny,nz,
     $    xlo, xhi, ylo, yhi, zlo, zhi, poles, lmax, depth, work, task)
      implicit none
      character*(*) task
      integer nx, ny, nz
      double precision g(nx, ny, nz)
      double precision xlo, xhi, ylo, yhi, zlo, zhi
      integer lmax, depth
      double precision work(0:*),
     $     poles(1:(lmax+1)**2,0:2**depth-1,0:2**depth-1,0:2**depth-1)
c
c     if (task .eq. 'poles') then
c     .   In grid is a solution to the Poisson equation for the 
c     .   "smooth" density.  The solution is zero on the boundary, odd
c     .   and periodic (computed using a sine transform).
c     .
c     .   Add into poles the multipole moments arising from the
c     .   boundary layer potential (solution to Laplaces eqn in the volume)
c     .   that when added to the period potential yeilds the free-space 
c     .   potential.  These poles arise from the surface charge but
c     .   the multipolar expansion is about the center of the surface volumes.
c     .
c     .   Note that the smooth grid does NOT include the boundary (since
c     .   the potential is known to be zero there).
c     else if (task .eq. 'near') then
c     .   The FMM package has been called to integrate the far field
c     .   potential arising from the poles computed with task='poles'
c     .   and the result of this is stored in poles.
c     .   However, the well-separated criterion of 2 means that the 
c     .   outermost three layers of the volume have no potential
c     .   in them.  So here we again compute the surface MPOLES
c     .   (this time centered at the surface) and add into poles
c     .   the corresponding surface potential in all boxes that
c     .   are from 1 to well_separated from the surface.  The potential
c     .   is then 
c     .       not computed in the outermost layer       (0)
c     .       accurate to (1.5)**lmax in the next layer (1)
c     .       accurate to (2.5)**lmax in the next layer (2)
c     .   Thus to be very safe the potential should be not be used in the
c     .   outermost two layers, however, assuming that that any densities
c     .   to be contracted with the potential are small near the 
c     .   boundary, it might be safe to only discard the outermost layer.
c     end if
c
c     Work should be dimensioned max((nx+2)*(ny+2),...)
c
      integer i, j, k, l, m, nx1, ny1, nz1, nx2, ny2, ind
      integer kx, ky, kz        ! FMM box indices
      double precision x, y, z, hx, hy, hz, hx_fmm, hy_fmm, hz_fmm,
     $     fourpi, xi, yj, zk, cartijk, dx, dy, dz
c
      logical opoles, onear
      integer maxl
      parameter (maxl = 25)
      double precision q2d(0:maxl,0:maxl), q(-maxl:maxl,0:maxl)
      double precision d(((maxl+1)*(maxl+2))/2, -maxl:maxl, 0:maxl)
c
      double precision diff6,diff10,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10
      diff6(f1,f2,f3,f4,f5,f6) = (
     $     (12.0d0/7.0d0)*f1  - (15.0d0/28.0d0)*f2 +
     $     (10.0d0/63.0d0)*f3 - (1.0d0/28.0d0)*f4  +
     $     (2.0d0/385.0d0)*f5 - (1.0d0/2772.0d0)*f6) ! err = h^13
      diff10(f1,f2,f3,f4,f5,f6,f7,f8,f9,f10) = (
     $     (20.0d0/11.0d0)*f1     - (15.0d0/22.0d0)*f2   +
     $     (40.0d0/143.0d0)*f3    - (15.0d0/143.0d0)*f4  +
     $     (24.0d0/715.0d0)*f5    - (5.0d0/572.0d0)*f6   + 
     $     (30.0d0/17017.0d0)*f7  - (5.0d0/19448.0d0)*f8 + 
     $     (10.0d0/415701.0d0)*f9 - (1.0d0/923780.0d0)*f10) ! err = h^21
c
      if (lmax .gt. maxl) call errquit('blm: lmax>maxl', lmax)
      onear = .false.
      opoles= .false.
      if (task .eq. 'poles') then
         opoles = .true.
      else if (task .eq. 'near') then
         onear = .true.
      else
         call errquit('boundary_layer_multipoles: unknown task ', 0)
      end if
c
      fourpi = 16.0d0*atan(1.0d0)
c
      nx1 = nx + 1
      ny1 = ny + 1
      nz1 = nz + 1
      nx2 = nx + 2
      ny2 = ny + 2
c
      hx  = (xhi-xlo)/dble(nx+1) ! Grid spacing
      hy  = (yhi-ylo)/dble(ny+1)
      hz  = (zhi-zlo)/dble(nz+1)
c
      hx_fmm = (xhi-xlo)/dble(2**depth) ! FMM box size
      hy_fmm = (yhi-ylo)/dble(2**depth)
      hz_fmm = (zhi-zlo)/dble(2**depth)
c
      call xlm_init
      call xlm_coeff(maxl,d)
c
c     For each face the outward normal derivative is computed in work.
c     Although the input 3-D grid does not include the zero boundary
c     values it is convenient to put these into work to simplify the
c     interpolation and quadrature to generate the moments.
c
c     x = xlo surface
c
      call fill_zero_boundary(work,ny1,nz1)
      do k = 1, nz
         do j = 1, ny
*            work(j + k*ny2) = diff6(g(1,j,k),g(2,j,k),g(3,j,k),g(4,j,k),
*     $           g(5,j,k),g(6,j,k))/(hx*fourpi)
            work(j + k*ny2) = diff10(
     $           g(1,j,k),g(2,j,k),g(3,j,k),g(4,j,k),g(5,j,k),
     $           g(6,j,k),g(7,j,k),g(8,j,k),g(9,j,k),g(10,j,k)) / 
     $           (hx*fourpi)
         end do
      end do
c
      kx = 0
      x = xlo
      if (opoles) then
         dx = -0.5d0*hx_fmm
      else
         dx = 0.0d0             ! Relies upon 0**0 being 1
      endif
      do kz = 0,2**depth-1
         z = (kz+0.5d0)*hz_fmm + zlo
         do ky = 0,2**depth-1
            y = (ky+0.5d0)*hy_fmm + ylo
            call surface_2d_cart_moments(work,ny1,nz1,hy,hz,
     $           lmax,ky,kz,hy_fmm,hz_fmm,q2d,maxl)
            do l = 0, lmax
               do m = -l,l
                  q(m,l) = 0.0d0
               end do
               ind = 1
               i = 0
               do i = l,0,-1
                  xi = dx**i
                  do j = l-i,0,-1
                     k = l-i-j
                     cartijk = xi*q2d(j,k)
                     do m = -l,l
                        q(m,l) = q(m,l) +  d(ind,m,l)*cartijk
                     end do
                     ind = ind + 1
                  end do
               end do
            end do
c
            if (opoles) then
               call xlm_accumulate_to_packed(lmax, poles(1,kx,ky,kz),
     $              q, maxl)
            else
               call boundary_layer_near(poles,lmax,depth,q,maxl,
     $              xlo,xhi,ylo,yhi,zlo,zhi,kx,ky,kz,x,y,z)
            end if
         end do
      end do
c     
c     x = xhi surface
c
      call fill_zero_boundary(work,ny1,nz1)
      do k = 1, nz
         do j = 1, ny
*            work(j + k*ny2) = diff6(g(nx,j,k),g(nx-1,j,k),g(nx-2,j,k),
*     $           g(nx-3,j,k),g(nx-4,j,k),g(nx-5,j,k))/(hx*fourpi)
            work(j + k*ny2) = diff10(
     $           g(nx,j,k),  g(nx-1,j,k),g(nx-2,j,k),g(nx-3,j,k),
     $           g(nx-4,j,k),g(nx-5,j,k),g(nx-6,j,k),g(nx-7,j,k),
     $           g(nx-8,j,k),g(nx-9,j,k)) / (hx*fourpi)
         end do
      end do
c
      kx = 2**depth-1
      x = xhi
      if (opoles) then
         dx = 0.5d0*hx_fmm
      else
         dx = 0.0d0
      endif
      do kz = 0,2**depth-1
         z = (kz+0.5d0)*hz_fmm + zlo
         do ky = 0,2**depth-1
            y = (ky+0.5d0)*hy_fmm + ylo
            call surface_2d_cart_moments(work,ny1,nz1,hy,hz,
     $           lmax,ky,kz,hy_fmm,hz_fmm,q2d,maxl)
            do l = 0, lmax
               do m = -l,l
                  q(m,l) = 0.0d0
               end do
               ind = 1
               i = 0
               do i = l,0,-1
                  xi = dx**i
                  do j = l-i,0,-1
                     k = l-i-j
                     cartijk = xi*q2d(j,k)
                     do m = -l,l
                        q(m,l) = q(m,l) +  d(ind,m,l)*cartijk
                     end do
                     ind = ind + 1
                  end do
               end do
            end do
c
            if (opoles) then
               call xlm_accumulate_to_packed(lmax, poles(1,kx,ky,kz),
     $           q, maxl)
            else
               call boundary_layer_near(poles,lmax,depth,q,maxl,
     $              xlo,xhi,ylo,yhi,zlo,zhi,kx,ky,kz,x,y,z)
            end if
         end do
      end do
c     
c     y = ylo surface
c
      call fill_zero_boundary(work,nx1,nz1)
      do k = 1, nz
         do i = 1, nx
*            work(i + k*nx2) = diff6(g(i,1,k),g(i,2,k),g(i,3,k),
*     $           g(i,4,k),g(i,5,k),g(i,6,k))/(hy*fourpi)
            work(i + k*nx2) = diff10(
     $           g(i,1,k),g(i,2,k),g(i,3,k),g(i,4,k),g(i,5,k),
     $           g(i,6,k),g(i,7,k),g(i,8,k),g(i,9,k),g(i,10,k)) / 
     $           (hy*fourpi)
         end do
      end do
c
      ky = 0
      y = ylo
      if (opoles) then
         dy = -0.5d0*hy_fmm
      else
         dy = 0.0d0
      endif
      do kx = 0,2**depth-1
         x = (kx+0.5d0)*hx_fmm + xlo
         do kz = 0,2**depth-1
            z = (kz+0.5d0)*hz_fmm + zlo
            call surface_2d_cart_moments(work,nx1,nz1,hx,hz,
     $           lmax,kx,kz,hx_fmm,hz_fmm,q2d,maxl)
            do l = 0, lmax
               do m = -l,l
                  q(m,l) = 0.0d0
               end do
               ind = 1
               i = 0
               do i = l,0,-1
                  do j = l-i,0,-1
                     yj = dy**j
                     k = l-i-j
                     cartijk = yj*q2d(i,k)
                     do m = -l,l
                        q(m,l) = q(m,l) +  d(ind,m,l)*cartijk
                     end do
                     ind = ind + 1
                  end do
               end do
            end do
c
            if (opoles) then
               call xlm_accumulate_to_packed(lmax, poles(1,kx,ky,kz),
     $              q, maxl)
            else
               call boundary_layer_near(poles,lmax,depth,q,maxl,
     $              xlo,xhi,ylo,yhi,zlo,zhi,kx,ky,kz,x,y,z)
            end if
         end do
      end do
c     
c     y = yhi surface
c
      call fill_zero_boundary(work,nx1,nz1)
      do k = 1, nz
         do i = 1, nx
*            work(i + k*nx2) = diff6(g(i,ny,k),g(i,ny-1,k),g(i,ny-2,k),
*     $           g(i,ny-3,k),g(i,ny-4,k),g(i,ny-5,k))/(hy*fourpi)
            work(i + k*nx2) = diff10(
     $           g(i,ny  ,k),g(i,ny-1,k),g(i,ny-2,k),g(i,ny-3,k),
     $           g(i,ny-4,k),g(i,ny-5,k),g(i,ny-6,k),g(i,ny-7,k),
     $           g(i,ny-8,k),g(i,ny-9,k))/(hy*fourpi)
         end do
      end do
c
      ky = 2**depth-1
      y = yhi
      if (opoles) then
         dy = 0.5d0*hy_fmm
      else
         dy = 0.0d0
      endif
      do kx = 0,2**depth-1
         x = (kx+0.5d0)*hx_fmm + xlo
         do kz = 0,2**depth-1
            z = (kz+0.5d0)*hz_fmm + zlo
            call surface_2d_cart_moments(work,nx1,nz1,hx,hz,
     $           lmax,kx,kz,hx_fmm,hz_fmm,q2d,maxl)
            do l = 0, lmax
               do m = -l,l
                  q(m,l) = 0.0d0
               end do
               ind = 1
               i = 0
               do i = l,0,-1
                  do j = l-i,0,-1
                     yj = dy**j
                     k = l-i-j
                     cartijk = yj*q2d(i,k)
                     do m = -l,l
                        q(m,l) = q(m,l) +  d(ind,m,l)*cartijk
                     end do
                     ind = ind + 1
                  end do
               end do
            end do
c
            if (opoles) then
               call xlm_accumulate_to_packed(lmax, poles(1,kx,ky,kz),
     $              q, maxl)
            else
               call boundary_layer_near(poles,lmax,depth,q,maxl,
     $              xlo,xhi,ylo,yhi,zlo,zhi,kx,ky,kz,x,y,z)
            end if
         end do
      end do
c     
c     z = zlo surface
c
      call fill_zero_boundary(work,nx1,ny1)
      do j = 1, ny
         do i = 1, nx
*            work(i + j*nx2) = diff6(g(i,j,1),g(i,j,2),g(i,j,3),
*     $           g(i,j,4),g(i,j,5),g(i,j,6))/(hz*fourpi)
            work(i + j*nx2) = diff10(
     $           g(i,j,1),g(i,j,2),g(i,j,3),g(i,j,4),g(i,j,5),
     $           g(i,j,6),g(i,j,7),g(i,j,8),g(i,j,9),g(i,j,10)) / 
     $           (hz*fourpi)
         end do
      end do
c
      kz = 0
      z = zlo
      if (opoles) then
         dz = -0.5d0*hz_fmm
      else
         dz = 0.0d0
      endif
      do ky = 0,2**depth-1
         y = (ky+0.5d0)*hy_fmm + ylo
         do kx = 0,2**depth-1
            x = (kx+0.5d0)*hx_fmm + xlo
            call surface_2d_cart_moments(work,nx1,ny1,hx,hy,
     $           lmax,kx,ky,hx_fmm,hy_fmm,q2d,maxl)
            do l = 0, lmax
               do m = -l,l
                  q(m,l) = 0.0d0
               end do
               ind = 1
               i = 0
               do i = l,0,-1
                  do j = l-i,0,-1
                     k = l-i-j
                     zk = dz**k
                     cartijk = q2d(i,j)*zk
                     do m = -l,l
                        q(m,l) = q(m,l) +  d(ind,m,l)*cartijk
                     end do
                     ind = ind + 1
                  end do
               end do
            end do
c
            if (opoles) then
               call xlm_accumulate_to_packed(lmax, poles(1,kx,ky,kz),
     $              q, maxl)
            else
               call boundary_layer_near(poles,lmax,depth,q,maxl,
     $              xlo,xhi,ylo,yhi,zlo,zhi,kx,ky,kz,x,y,z)
            end if
         end do
      end do
c     
c     z = zhi surface
c
      call fill_zero_boundary(work,nx1,ny1)
      do j = 1, ny
         do i = 1, nx
*            work(i + j*nx2) = diff6(g(i,j,nx),g(i,j,nx-1),g(i,j,nx-2),
*     $           g(i,j,nx-3),g(i,j,nx-4),g(i,j,nx-5))/(hz*fourpi)
            work(i + j*nx2) = diff10(
     $           g(i,j,nx  ),g(i,j,nx-1),g(i,j,nx-2),g(i,j,nx-3),
     $           g(i,j,nx-4),g(i,j,nx-5),g(i,j,nx-6),g(i,j,nx-7),
     $           g(i,j,nx-8),g(i,j,nx-9))/(hz*fourpi)
         end do
      end do
c
      kz = 2**depth-1
      z = zhi
      if (opoles) then
         dz = 0.5d0*hz_fmm
      else
         dz = 0.0d0
      endif
      do ky = 0,2**depth-1
         y = (ky+0.5d0)*hy_fmm + ylo
         do kx = 0,2**depth-1
            x = (kx+0.5d0)*hx_fmm + xlo
            call surface_2d_cart_moments(work,nx1,ny1,hx,hy,
     $           lmax,kx,ky,hx_fmm,hy_fmm,q2d,maxl)
            do l = 0, lmax
               do m = -l,l
                  q(m,l) = 0.0d0
               end do
               ind = 1
               i = 0
               do i = l,0,-1
                  do j = l-i,0,-1
                     k = l-i-j
                     zk = dz**k
                     cartijk = q2d(i,j)*zk
                     do m = -l,l
                        q(m,l) = q(m,l) +  d(ind,m,l)*cartijk
                     end do
                     ind = ind + 1
                  end do
               end do
            end do
c
            if (opoles) then
               call xlm_accumulate_to_packed(lmax, poles(1,kx,ky,kz),
     $              q, maxl)
            else
               call boundary_layer_near(poles,lmax,depth,q,maxl,
     $              xlo,xhi,ylo,yhi,zlo,zhi,kx,ky,kz,x,y,z)
            end if
         end do
      end do
c     
      end
c$$$      subroutine oldsint3d(nx,ny,nz,g)
c$$$      implicit none
c$$$      integer nx, ny, nz
c$$$      double precision g(nx,ny,nz)
c$$$c
c$$$c     Destructive 3D sin transform - no blocking for the transpose.
c$$$c
c$$$      integer maxn
c$$$      parameter (maxn = 1024)
c$$$      double precision wsave(10*maxn), a(maxn)
c$$$      integer i, j, k
c$$$c
c$$$      if (ny .gt. maxn) call errquit('sint3d: ny>maxn', ny)
c$$$      if (nz .gt. maxn) call errquit('sint3d: nz>maxn', nz)
c$$$c
c$$$      call sinti(nx,wsave)
c$$$      do k = 1, nz
c$$$         do j = 1, ny
c$$$            call sint(nx,g(1,j,k),wsave)
c$$$         end do
c$$$      end do
c$$$c
c$$$      call sinti(ny,wsave)
c$$$      do k = 1, nz
c$$$         do i = 1, nx
c$$$            do j = 1, ny
c$$$               a(j) = g(i,j,k)
c$$$            end do
c$$$            call sint(ny,a,wsave)
c$$$            do j = 1, ny
c$$$               g(i,j,k) = a(j)
c$$$            end do
c$$$         end do
c$$$      end do
c$$$c
c$$$      call sinti(nz,wsave)
c$$$      do j = 1, ny
c$$$         do i = 1, nx
c$$$            do k = 1, nz
c$$$               a(k) = g(i,j,k)
c$$$            end do
c$$$            call sint(nz,a,wsave)
c$$$            do k = 1, nz
c$$$               g(i,j,k) = a(k)
c$$$            end do
c$$$         end do
c$$$      end do
c$$$c
c$$$      end
c$$$      subroutine sint3d(nx,ny,nz,g)
c$$$      implicit none
c$$$      integer nx, ny, nz
c$$$      double precision g(nx,ny,nz)
c$$$c
c$$$c     Destructive 3D sin transform - with blocking
c$$$c
c$$$      integer maxn, ncol
c$$$      parameter (maxn = 1025, ncol = 31)
c$$$      double precision wsave(10*maxn), a(maxn,ncol)
c$$$      integer i, j, k, ilo, ihi
c$$$c
c$$$      if (ny .gt. maxn) call errquit('sint3d: ny>maxn', ny)
c$$$      if (nz .gt. maxn) call errquit('sint3d: nz>maxn', nz)
c$$$c
c$$$      call sinti(nx,wsave)
c$$$      do k = 1, nz
c$$$         do j = 1, ny
c$$$            call sint(nx,g(1,j,k),wsave)
c$$$         end do
c$$$      end do
c$$$c
c$$$      call sinti(ny,wsave)
c$$$      do k = 1, nz
c$$$         do ilo = 1, nx, ncol
c$$$            ihi = min(nx,ilo+ncol-1)
c$$$            do j = 1, ny
c$$$               do i = ilo, ihi
c$$$                  a(j,i-ilo+1) = g(i,j,k)
c$$$               enddo
c$$$            enddo
c$$$            do i = ilo, ihi
c$$$               call sint(ny,a(1,i-ilo+1),wsave)
c$$$            end do
c$$$            do j = 1, ny
c$$$               do i = ilo, ihi
c$$$                  g(i,j,k) = a(j,i-ilo+1)
c$$$               enddo
c$$$            enddo
c$$$         end do
c$$$      end do
c$$$c
c$$$      call sinti(nz,wsave)
c$$$      do j = 1, ny
c$$$         do ilo = 1, nx, ncol
c$$$            ihi = min(nx,ilo+ncol-1)
c$$$            do k = 1, nz
c$$$               do i = ilo, ihi
c$$$                  a(k,i-ilo+1) = g(i,j,k)
c$$$               enddo
c$$$            end do
c$$$            do i = ilo, ihi
c$$$               call sint(nz,a(1,i-ilo+1),wsave)
c$$$            enddo
c$$$            do k = 1, nz
c$$$               do i = ilo, ihi
c$$$                  g(i,j,k) = a(k,i-ilo+1)
c$$$               enddo
c$$$            end do
c$$$         end do
c$$$      end do
c$$$c
c$$$      end
      subroutine vsint3d(nx,ny,nz,g)
      implicit none
      integer nx, ny, nz
      double precision g(nx,ny,nz)
c
c     Destructive 3D sin transform - using V-FFTPACK
c
      integer maxn, ncol
      parameter (maxn = 1025, ncol = 63)
      double precision wsave(10*maxn), scale
      double precision x(ncol,maxn+1), xt(ncol,maxn+1)
      integer i, j, k, ilo, ihi, jlo, jhi
c
      if (ny .gt. maxn) call errquit('sint3d: ny>maxn', ny)
      if (nz .gt. maxn) call errquit('sint3d: nz>maxn', nz)
c
      scale = sqrt(2.0d0*dble(nx+1))
      call vsinti(nx,wsave)
      do k = 1, nz
         do jlo = 1, ny, ncol
            jhi = min(ny,jlo+ncol-1)
            do j = jlo, jhi
               do i = 1, nx
                  x(j-jlo+1,i) = g(i,j,k)
               enddo
            enddo
            call vsint(jhi-jlo+1,nx,x,xt,ncol,wsave)
            do j = jlo, jhi
               do i = 1, nx
                  g(i,j,k) = x(j-jlo+1,i) * scale
               enddo
            enddo
         end do
      end do
c
      scale = sqrt(2.0d0*dble(ny+1))
      call vsinti(ny,wsave)
      do k = 1, nz
         do ilo = 1, nx, ncol
            ihi = min(nx,ilo+ncol-1)
            do j = 1, ny
               do i = ilo, ihi
                  x(i-ilo+1,j) = g(i,j,k)
               enddo
            enddo
            call vsint(ihi-ilo+1,ny,x,xt,ncol,wsave)
            do j = 1, ny
               do i = ilo, ihi
                  g(i,j,k) = x(i-ilo+1,j)*scale
               enddo
            enddo
         end do
      end do
c
      scale = sqrt(2.0d0*dble(nz+1))
      call vsinti(nz,wsave)
      do j = 1, ny
         do ilo = 1, nx, ncol
            ihi = min(nx,ilo+ncol-1)
            do k = 1, nz
               do i = ilo, ihi
                  x(i-ilo+1,k) = g(i,j,k)
               enddo
            end do
            call vsint(ihi-ilo+1,nz,x,xt,ncol,wsave)
            do k = 1, nz
               do i = ilo, ihi
                  g(i,j,k) = x(i-ilo+1,k)*scale
               enddo
            end do
         end do
      end do
c$$$      scale = 2.0*dble(nz+1)
c$$$      call vsinti(nz,wsave)
c$$$      do j = 1, ny
c$$$         do ilo = 1, nx, ncol
c$$$            ihi = min(nx,ilo+ncol-1)
c$$$            do k = 1, nz
c$$$               do i = ilo, ihi
c$$$                  x(i-ilo+1,k) = g(i,j,k)
c$$$               enddo
c$$$            end do
c$$$            call vsint(ihi-ilo+1,nz,x,xt,ncol,wsave)
c$$$            do k = 1, nz
c$$$               do i = ilo, ihi
c$$$                  g(i,j,k) = x(i-ilo+1,k)*scale
c$$$               enddo
c$$$            end do
c$$$         end do
c$$$      end do
c
      end
      subroutine fill_zero_boundary(work,ny1,nz1)
      implicit none
      integer ny1, nz1
      double precision work(0:ny1,0:nz1)
c
c     Fill the outside edges of work in with zeroes
c
      integer i
c
      do i = 0, ny1
         work(i,0) = 0.0d0
         work(i,nz1) = 0.0d0
      end do
      do i = 0,nz1
         work(0,i) = 0.0d0
         work(ny1,i) = 0.0d0
      end do
c
      end
      subroutine surface_2d_cart_moments(work,ny1,nz1,hy,hz,
     $     lmax,ky,kz,hy_fmm,hz_fmm,q2d,lq2d)
      implicit none
      integer ny1,nz1
      double precision work(0:ny1,0:nz1)
      double precision hy, hz, hy_fmm, hz_fmm
      integer lmax, ky, kz, lq2d
      double precision q2d(0:lq2d,0:lq2d)
c
c     Work() defines a uniform grid over (0:ny1*hy,0:nz1*hz).
c     Return in q2d(i,j) the cartesian multipoles of the
c     sub-region (ky*hy_fmm:(ky+1)*hy_fmm,kz*hz_fmm:(kz+1)*hz_fmm)
c     about the center of that region.
c
      double precision box_zlo, box_zhi, box_zo, vol_zlo, vol_zhi
      integer pts_zlo, pts_zhi, nptz
      double precision box_ylo, box_yhi, box_yo, vol_ylo, vol_yhi
      integer pts_ylo, pts_yhi, npty
      integer i, j, k, l, m
      double precision yy, zz
c
      integer order, ng, maxnpt, maxng, maxl, npt
      parameter (maxnpt=65,maxng=16,maxl=25)
      double precision tmp(maxnpt,maxnpt), f(maxng,maxng)
      double precision y(maxng),z(maxng),yw(maxng),zw(maxng)
      double precision yp(maxng,0:maxl),zp(maxng,0:maxl),yf(0:maxl)
c
      if (lmax .gt. maxl) call errquit('s2dcm: lmax>maxl',lmax)
c
      order = lmax + 2          ! Order of interpolation
      ng    = (order-1)/2 + 2   ! Order of quadrature
      npt   = 1.5*order         ! Min no. of points for interpolation
      if (ng .gt. maxng) call errquit
     $     ('surf_2d_cart_mpole: maxng too small',ng)
c
c     Determine the box and the points contained in it making
c     sure that we use sufficient points to enclose the entire box.
c
      box_zlo = kz*hz_fmm       ! FMM box info
      box_zhi = box_zlo + hz_fmm
      box_zo  = box_zlo + hz_fmm*0.5d0
      pts_zlo = box_zlo/hz
      pts_zhi = (box_zhi-1d-12*hz)/hz + 1
      if ((pts_zhi-pts_zlo) .lt. npt) then
         pts_zlo = max(0,pts_zlo-(npt-(pts_zhi-pts_zlo))/2)
         pts_zhi = min(nz1,pts_zlo+npt)
         pts_zlo = max(0,pts_zlo-(npt-(pts_zhi-pts_zlo)))
      end if
      nptz = pts_zhi - pts_zlo + 1
      call gauleg(box_zlo, box_zhi, z, zw, ng)
      vol_zlo = pts_zlo*hz
      vol_zhi = pts_zhi*hz
c     
      box_ylo = ky*hy_fmm
      box_yhi = box_ylo + hy_fmm
      box_yo  = box_ylo + hy_fmm*0.5d0
      pts_ylo = box_ylo/hy
      pts_yhi = (box_yhi-1d-12*hy)/hy + 1
      if ((pts_yhi-pts_ylo) .lt. npt) then
         pts_ylo = max(0,pts_ylo-(npt-(pts_yhi-pts_ylo))/2)
         pts_yhi = min(ny1,pts_ylo+npt)
         pts_ylo = max(0,pts_ylo-(npt-(pts_yhi-pts_ylo)))
      end if
      npty = pts_yhi - pts_ylo + 1
      call gauleg(box_ylo, box_yhi, y, yw, ng)
      vol_ylo = pts_ylo*hy
      vol_yhi = pts_yhi*hy
c
      if (npty.gt.maxnpt .or. nptz.gt.maxnpt) call errquit
     $     ('surf_2d_cart_mpole: maxnpt too small',max(npty,nptz))
c
c     Copy points in box to workspace
c
      do k = pts_zlo,pts_zhi
         do j = pts_ylo,pts_yhi
            tmp(j-pts_ylo+1,k-pts_zlo+1) = work(j,k)
         end do
      end do
c
c     Interpolate to the quadrature points
c
      call tn_interp_2d(npty,nptz,vol_ylo, vol_yhi,
     $     vol_zlo, vol_zhi, ng, ng, y, z, tmp, maxnpt,
     $     f, maxng, order)
*      call lagrange_interp_2d(npty,nptz,vol_ylo, vol_yhi,
*     $     vol_zlo, vol_zhi, ng, ng, y, z, tmp, maxnpt, maxnpt,
*     $     f, maxng, maxng)
c
c     Generate powers of cartesian coords at quadrature points
c     including the quadrature weights
c
      do i = 1,ng
         yy = yw(i)
         zz = zw(i)
         do j = 0,lmax
            yp(j,i) = yy
            zp(j,i) = zz
            yy = yy*(y(i)-box_yo)
            zz = zz*(z(i)-box_zo)
         end do
      end do
c
c     Do the quadrature to generate the cartesian moments
c
c     q2d(j,k) = sum(l,m) (yp(j,l)*zp(k,m)*f(l,m))
c     .        = sum(m) zp(k,m)* [ sum(l) yp(j,l)*f(l,m) ]
c     .        = sum(m) zp(k,m)*yf(j,m)
c
      do k = 0, lmax
         do j = 0, lmax
            q2d(j,k) = 0.0d0
         end do
      end do
c
      do m = 1, ng
         do j = 0,lmax
            yf(j) = 0.0d0
         end do
         do l = 1, ng
            do j = 0, lmax
               yf(j) = yf(j) + yp(j,l)*f(l,m)
            end do
         end do
         do k = 0,lmax
            do j = 0,lmax
               q2d(j,k) = q2d(j,k) + yf(j)*zp(k,m)
            end do
         end do
      end do
c
      end
      subroutine boundary_layer_near(poles, lmax, depth, q, lq,
     $     xlo, xhi, ylo, yhi, zlo, zhi, kx, ky, kz, x, y, z)
      implicit none
#include "fmmP.fh"
      integer lmax, depth, kx, ky, kz, lq
      double precision xlo, xhi, ylo, yhi, zlo, zhi, q(-lq:lq,0:lq),
     $     poles(1:(lmax+1)**2,0:2**depth-1,0:2**depth-1,0:2**depth-1)
c
      integer dkx, dky, dkz, kxn, kyn, kzn, twod
      double precision hx_fmm, hy_fmm, hz_fmm, x, y, z, xn, yn, zn
c
      integer maxl
      parameter (maxl = 25)
      double precision pot(-maxl:maxl,0:maxl)
c
      if (lmax .gt. maxl) call errquit('bln: lmax>maxl',lmax)
c
      twod = 2**depth
      hx_fmm = (xhi-xlo)/dble(twod) ! FMM box size
      hy_fmm = (yhi-ylo)/dble(twod)
      hz_fmm = (zhi-zlo)/dble(twod)
c
      do dkz = -well_separated, well_separated
         kzn = kz + dkz
         zn = zlo + (kzn+0.5d0)*hz_fmm
         if (kzn.ge.1 .and. kzn.lt.twod-1) then
            do dky = -well_separated, well_separated
               kyn = ky + dky
               yn = ylo + (kyn+0.5d0)*hy_fmm
               if (kyn.ge.1 .and. kyn.le.(twod-2)) then
                  do dkx = -well_separated, well_separated
                     kxn = kx + dkx
                     xn = xlo + (kxn+0.5d0)*hx_fmm
                     if (kxn.ge.1 .and. kxn.lt.twod-1) then
                        call xlm_new_multipole_to_local(lmax,
     $                       xn, yn, zn, pot, maxl, 
     $                       x , y , z , q  , lq)
c
*                        write(6,*) ' from ', kx, ky, kz, ' to ',
*     $                       kxn, kyn, kzn
*                        write(6,*) x, y, z, xn, yn, zn,
*     $                       sqrt((x-xn)**2+(y-yn)**2+(z-zn)**2)/hx_fmm
                        call xlm_accumulate_to_packed(lmax, 
     $                       poles(1,kxn,kyn,kzn),pot, maxl)
                     end if
                  end do
               end if
            end do
         end if
      end do
c
      end
      subroutine fft_interp_3d(g, nx, ny, nz, gg, nnx, nny, nnz)
      implicit none
c     
      integer nx, ny, nz, nnx, nny, nnz
      double precision g(nx,ny,nz), gg(nnx,nny,nnz)
c     
c     Given two uniform grids over the same volume interpolate
c     from the coarser (g) to the finer (gg) using FFT.
c
      double precision scale
      integer i, j, k
c
      call sint3d(nx,ny,nz,g)
      scale = 0.125d0/(dble(nx+1)*dble(ny+1)*dble(nz+1))
      call dfill(nnx*nny*nnz,0.0d0,gg,1)
      do k = 1, nz
         do j = 1, ny
            do i = 1, nx
               gg(i,j,k) = scale*g(i,j,k)
            enddo
         enddo
      enddo
c
      call sint3d(nnx,nny,nnz,gg)
c
      end
