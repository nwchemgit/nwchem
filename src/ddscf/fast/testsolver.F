      implicit none
#include "mafdecls.fh"
      integer nx, ny, nz, lmax, depth, lenxlm, nnx, nny, nnz
      parameter (nx=242, ny=242, nz=242, lmax=16, depth=3)
      parameter (nnx=75, nny=75, nnz=75)
c
      double precision g(nx,ny,nz), gg(nnx, nny, nnz), value,
     $     poles(1:(lmax+1)**2,0:2**depth-1,0:2**depth-1,0:2**depth-1)
c
c     solve Poisson equation subject to potential (and density)
c     being zero on the boundary (also periodic) ... sine transform.
c
c     Use a difference between two gaussians to produce a potential
c     that decays exponentially so that boundary conditions are exact.
c
      double precision xlo, xhi, ylo, yhi, zlo, zhi, x, y, z, pi
      double precision hx, hy, hz, hx_fmm, hy_fmm, hz_fmm
      double precision hhx, hhy, hhz
      double precision errmax, errrms, err, pot, dx, dy, dz
      integer i, j, k, ii, jj, kk
c
      integer n, l, m           ! Test function
      double precision alpha, a, b, c, func
c
      double precision fnlm, vnlm, fmm_potential
      external fnlm, vnlm, fmm_potential
c
      if (.not. ma_init(mt_dbl, 64000000, 64000000)) call errquit
     $     (' failed to init ma ', 128000000)
c
**      call initfm
      call anl_fit_init
      call xlm_init
c
      lenxlm = (lmax+1)**2
      call dfill(lenxlm*8**depth, 0.0d0, poles, 1)
c
      n = 3
      l = 1
      m = -1
      alpha = 0.95189d0
      a = 0.8d0
      b = -1.11d0
      c = 0.77777d0
c
      pi = 4.0d0*atan(1.0d0)
      xlo = -8.0d0
      xhi =  8.0d0
      ylo = -8.0d0
      yhi =  8.0d0
      zlo = -8.0d0
      zhi =  8.0d0
c
      hx  = (xhi-xlo)/dble(nx+1)
      hy  = (yhi-ylo)/dble(ny+1)
      hz  = (zhi-zlo)/dble(nz+1)
c
      hhx  = (xhi-xlo)/dble(nnx+1)
      hhy  = (yhi-ylo)/dble(nny+1)
      hhz  = (zhi-zlo)/dble(nnz+1)
c     
      hx_fmm = (xhi-xlo)/dble(2**depth) ! FMM box size
      hy_fmm = (yhi-ylo)/dble(2**depth)
      hz_fmm = (zhi-zlo)/dble(2**depth)
c
      do k = 1, nz
         z = dble(k)*hz + zlo
         do j = 1, ny
            y = dble(j)*hy + ylo
            do i = 1, nx
               x = dble(i)*hx + xlo
*              g(i,j,k) = f(x,y,z)
               g(i,j,k) = fnlm(n,l,m,alpha,a,b,c,x,y,z)
            end do
         end do
      end do
c
c$$$      call fft_interp_3d(g, nx, ny, nz, gg, nnx, nny, nnz)
c$$$      do k = 1, nnz
c$$$         z = dble(k)*hhz + zlo
c$$$         do j = 1, nny
c$$$            y = dble(j)*hhy + ylo
c$$$            do i = 1, nnx
c$$$               x = dble(i)*hhx + xlo
c$$$               value = fnlm(n,l,m,alpha,a,b,c,x,y,z)
c$$$               if (abs(value) .gt. 0.0001d0) then
c$$$                  if (abs(value-gg(i,j,k)).gt.1d-12) then
c$$$                     write(6,*) i, j, k, value, gg(i,j,k),
c$$$     $                    gg(i,j,k)/value
c$$$                     stop 44
c$$$                  endif
c$$$               endif
c$$$            end do
c$$$         end do
c$$$      end do
c$$$      stop 31
c
      call solver(g, nx, ny, nz, poles, depth, lmax, 
     $     xlo, xhi, ylo, yhi, zlo, zhi)
c
c     Test points must be located 2 FMM boxes away from the
c     boundary since we are not computing the direct contribution 
c     of the boundary layer charges.
c
      ii = 0
      jj = 0
      kk = 0
      errmax = 0.0d0
      errrms = 0.0d0
      do k = 1, nz
         z = dble(k)*hz + zlo
         do j = 1, ny
            y = dble(j)*hy + ylo
            do i = 1, nx
               x = dble(i)*hx + xlo
c
               dx = min(x-xlo, xhi-x)/2.0d0 ! Half distance from closest edge
               dy = min(y-ylo, yhi-y)/2.0d0
               dz = min(z-zlo, zhi-z)/2.0d0
               if (dx.gt.hx_fmm.and.dy.gt.hy_fmm.and.dz.gt.hz_fmm) then
*                  pot = v(x,y,z) 

                  pot = vnlm(n,l,m,alpha,a,b,c,x,y,z)

                  err = g(i,j,k)+ fmm_potential(lmax,depth,
     $                 xhi-xlo,yhi-ylo,zhi-zlo,poles,x-xlo,y-ylo,z-zlo)
     $                 -pot
*     write(6,78) x, g(i,j,k), pot
*     78            format(1x,f12.6,2f20.10)
                  errrms = errrms + err*err

                  if (k.eq.(nz/2) .and. j.eq.(ny/2)) then
                     func= fnlm(n,l,m,alpha,a,b,c,x,y,z)
                     write(6,891) i,pot,err+pot,(err+pot)/pot, func
 891                 format(i5,5f14.10)
                  endif

                  if (abs(err) .gt. errmax) then
                     ii = i
                     jj = j
                     kk = k
                     errmax = abs(err)
                  end if
               end if
            end do
         end do
      end do
c
      x = xlo + hx*nx*.5d0
      y = ylo + hy*ny*.5d0
      z = zlo + hz*nz*.5d0
      pot = vnlm(n,l,m,alpha,a,b,c,x,y,z)
      write(6,*) ' actual boundary layer potential ', 
     $     pot - g(nx/2,ny/2,nz/2) ! Assumes even no. points
      write(6,*) '    FMM boundary layer potential ',
     $     fmm_potential(lmax,depth,
     $     xhi-xlo,yhi-ylo,zhi-zlo,poles,x-xlo,y-ylo,z-zlo)
      errrms = sqrt(errrms/(nx*ny*nz))
      write(6,1) nx,ny,nz, xlo, xhi, ylo,yhi, zlo, zhi, 
     $     errmax, errrms
 1    format(' nx, ny, nz = ', 3i5/
     $       ' (xlo, xhi) = ', 2f8.2/
     $       ' (ylo, yhi) = ', 2f8.2/
     $       ' (zlo, zhi) = ', 2f8.2/
     $       '    max err = ',1p,d9.2,0p/
     $       '    rms err = ',1p,d9.2,0p)
      write(6,*) ii, jj, kk
c
      end
