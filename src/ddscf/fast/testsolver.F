*
* $Id$
*
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
      integer nx, ny, nz, lmax, depth, lenxlm, nnx, nny, nnz
      parameter (nx=134, ny=134, nz=134, lmax=16, depth=3)
      parameter (nnx=75, nny=75, nnz=75)
c     
      double precision g(nx,ny,nz) ! gg(nnx, nny, nnz)
      integer tree(300)
c     
c     solve Poisson equation subject to potential (and density)
c     being zero on the boundary (also periodic) ... sine transform.
c     
c     Use a difference between two gaussians to produce a potential
c     that decays exponentially so that boundary conditions are exact.
c     
      double precision q(-lmax:lmax,0:lmax), qn(-lmax:lmax,0:lmax)
      double precision xlo, xhi, ylo, yhi, zlo, zhi, x, y, z, pi
      double precision hx, hy, hz, hx_fmm, hy_fmm, hz_fmm
      double precision hhx, hhy, hhz
      double precision errmax, errrms, err, pot
      integer i, j, k, ii, jj, kk, lll, numsample
      parameter (numsample = 2000)
c     
      integer n, l, m           ! Test function
      integer maxorder
      parameter (maxorder = 17)
      double precision alpha, a, b, c, poti, potfmm,
     $     qfit((maxorder+1)*(maxorder+3))
      integer order, ninterp, ilo, jlo, klo
c     
      double precision fnlm, vnlm, fmm_potential,
     $     tn_interp_3d_point, util_random
      external fnlm, vnlm, fmm_potential
c     
      if (.not. ma_init(mt_dbl, 64000000, 64000000)) call errquit
     $     (' failed to init ma ', 128000000, MA_ERR)
c     
*     *      call initfm
      call anl_fit_init
      call util_print_enable('fastjtime')
c     
      lenxlm = (lmax+1)**2
      call xlm_init
c
      call xlm(lmax,1.0d0,2.0d0,3.0d0, q, lmax)
      call xlm_slower(lmax,1.0d0,2.0d0,3.0d0, qn, lmax)
      do l = 0, lmax
         do m = -l, l
            if (abs((qn(m,l)-q(m,l))/qn(m,l)) .gt. 1d-12) then
               write(6,*) m, l, qn(m,l), q(m,l)
               stop 66
            endif
         enddo
      enddo
      write(6,*) ' faster XLM is OK ?'
c     
      n =  3
      l =  1
      m =  1
      alpha = 1.35189d0
      a = 0.0d0 ! 0.8d0
      b = 0.0d0 ! -1.11d0
      c = 0.0d0 ! 0.77777d0
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
               g(i,j,k) = fnlm(n,l,m,alpha,a,b,c,x,y,z)
            end do
         end do
      end do
c
*      write(6,*) ' slice of g'
*      do i = 1, nx
*         write(6,*) i, g(i,ny/2,nz/2)
*      enddo
c     
      call fmm_tree_create(depth, lmax, tree)
      call solver(g, nx, ny, nz, tree, depth, lmax, 
     $     xlo, xhi, ylo, yhi, zlo, zhi, .false., .false.)
c     
c     Test points must be located 2 FMM boxes away from the
c     boundary since we are not computing the direct contribution 
c     of the boundary layer charges.
c     
      ii = 0
      jj = 0
      kk = 0
      do order = 3,maxorder,2
         call tn_fitting_matrix(order,order+3,qfit)
         errmax = util_random(12345)  ! So that always look at the same points
         errmax = 0.0d0
         errrms = 0.0d0
         do lll = 1, numsample
            x = xlo + 2.0d0*hx_fmm + util_random(0)*
     X           (xhi-xlo-4.0d0*hx_fmm)
            y = ylo + 2.0d0*hy_fmm + util_random(0)*
     Y           (yhi-ylo-4.0d0*hy_fmm)
            z = zlo + 2.0d0*hz_fmm + util_random(0)*
     Z           (zhi-zlo-4.0d0*hz_fmm)
c     
            ninterp = order+1+2 ! Works best of all options tried
            ilo = max(1, int((x-xlo)/hx)-ninterp/2)
            ilo = min(ilo,nx-ninterp+1)
            jlo = max(1, int((y-ylo)/hy)-ninterp/2)
            jlo = min(jlo,ny-ninterp+1)
            klo = max(1, int((z-zlo)/hz)-ninterp/2)
            klo = min(klo,nz-ninterp+1)
            poti = tn_interp_3d_point(g,nx,ny,nz,hx,hy,hz,x-xlo,y-ylo,
     $           z-zlo,ninterp,order,qfit)
            potfmm = fmm_potential(depth,lmax,tree,
     $           xhi-xlo,yhi-ylo,zhi-zlo,x-xlo,y-ylo,z-zlo)
c     
            pot = vnlm(n,l,m,alpha,a,b,c,x,y,z)
c
            err = poti + potfmm - pot
            errrms = errrms + err*err
            
            if (abs(err) .gt. errmax) then
               ii = nint((x-xlo)/hx)
               jj = nint((y-ylo)/hy)
               kk = nint((z-zlo)/hz)
               errmax = abs(err)
            end if
         end do
         write(6,*) ' order ', order
         errrms = sqrt(errrms/dble(numsample))
         write(6,1) order, errmax, errrms
 1       format(i5, 1p, 2d9.2)
      enddo 
c     
      x = xlo + hx*nx*.5d0
      y = ylo + hy*ny*.5d0
      z = zlo + hz*nz*.5d0
      pot = vnlm(n,l,m,alpha,a,b,c,x,y,z)
      write(6,*) n,l,m,alpha,a,b,c,x,y,z
      write(6,*) ' POT ', pot
      write(6,*) ' G ', g(nx/2,ny/2,nz/2)
      write(6,*) ' actual boundary layer potential ', 
     $     pot - g(nx/2,ny/2,nz/2) ! Assumes even no. points
      write(6,*) '    FMM boundary layer potential ',
     $     fmm_potential(depth,lmax,tree,
     $     xhi-xlo,yhi-ylo,zhi-zlo,x-xlo,y-ylo,z-zlo)
      write(6,*) ii, jj, kk
c
*      write(6,*) ' slice of g'
*      do i = 1, nx
*         write(6,*) i, g(i,ny/2,nz/2)
*      enddo
c$$$c
c$$$c     Test phase of xlm_potential versus vnlm
c$$$c
c$$$      call xlm_norm(lmax, qn, lmax)
c$$$      x = xhi
c$$$      y = yhi + 3.1519
c$$$      z = zhi + 2.71
c$$$      do n = 0, 5
c$$$         do l = n, 0, -2
c$$$            do m = -l, l
c$$$               call dfill((2*lmax+1)*(lmax+1),0.0d0,q,1)
c$$$               q(m,l) = sqrt(pi)*double_factorial(n+l+1) /
c$$$     $              ((2.0d0**((n+l+4)/2))*(alpha**((n+l+3.0d0)/2.0d0))*
c$$$     $              qn(m,l)**2)
c$$$               poti = xlm_potential_test(x,y,z,a,b,c,q,lmax,lmax)
c$$$               pot  = vnlm(n,l,m,alpha,a,b,c,x,y,z)
c$$$               write(6,771) n, l, m, poti, pot, abs(pot-poti)
c$$$ 771           format(1x,3i3,2x,2f12.8,2x,1p,d9.2)
c$$$            enddo
c$$$         enddo
c$$$      enddo
c     
      end
