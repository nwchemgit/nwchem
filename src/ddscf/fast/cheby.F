      subroutine tn_fitting_coeff(n,m,q,f,c,r)
*
* $Id$
*
      implicit none
#include "errquit.fh"
      integer n, m
      double precision q(0:n,1:m), f(1:m), c(0:n), r
c
c     Given the matrix from tn_fitting_matrix return
c     the fitting coeffs for a given function f() and
c     the residual at the fitting points
c
      integer i, j
      integer maxn
      parameter (maxn=64)
      double precision t(0:maxn), x, sum
c
      if (n .gt. maxn) call errquit('tn_f_c: n>maxn ', n, UNKNOWN_ERR)
c
      do i = 0,n
         c(i) = 0.0d0
      enddo
c
      do j = 1,m
         do i = 0,n
            c(i) = c(i) + q(i,j)*f(j)
         enddo
      enddo
c
      r = 0.0d0
      do j = 1, m
         x = -1.0d0 + 2.0d0*(j-1)/dble(m-1)
         call tn_evaluate(n,x,t)
         sum = 0.0d0
         do i = 0,n
            sum = sum + c(i)*t(i)
         enddo
         r = r + (f(j)-sum)**2
      enddo
c
      r = sqrt(r)
c
      end
      subroutine tn_collate_matrix(nx,xlo,xhi,nptx,x,order,q,lq)
      implicit none
#include "errquit.fh"
      integer nx, nptx, lq, order
      double precision xlo, xhi, x(nptx), q(lq,nx)
c
c     Return in q(nx,nptx) the matrix that uses Chebyshev interpolation
c     of the requested order to interpolate from nx points uniformly 
c     distributed in [xlo,xhi] to the list of points x(i), i=1,nptx
c
      integer maxorder, maxn
      parameter (maxorder=35,maxn=64)
      double precision qq((maxorder+1)*maxn),t(0:maxorder)
      double precision xx
      integer i,j,k,ind
c
      if (nx .gt. maxn) call errquit('tn_collate: nx>maxn',nx,
     &       UNKNOWN_ERR)
      if (order.gt.maxorder) call errquit('tn_collate: order',order,
     &       UNKNOWN_ERR)
c
      call tn_fitting_matrix(order,nx,qq)
      do j = 1, nx
         do i = 1, nptx
            q(i,j) = 0.0d0
         enddo
      enddo
      do i = 1, nptx
         xx = 2.0d0*(x(i)-xlo)/(xhi-xlo) - 1.0d0
         call tn_evaluate(order,xx,t)
         ind = 1
         do j = 1, nx
            do k = 0,order
               q(i,j) = q(i,j) + t(k)*qq(ind+k)
            enddo
            ind = ind + order + 1
         enddo
      enddo
c
      end
      subroutine tn_fitting_matrix(n,m,q)
      implicit none
#include "errquit.fh"
      integer n,m
      double precision q(0:n,1:m)
c
c     Return in q the matrix that forms a least-squares fit of
c     m points uniformly distributed on [-1,1] to a Chebyshev
c     expansion of order n.
c
c     Q = (At.A)^-1.A where Aij = Tj(xi), xi = -1 + 2*(i-1)/(m-1)
c
c     Then, given a function f() evaluated at the m uniform points 
c     a fit is computed as
c
c     ci = sum(j) Qij.f(xj) , or c = Qf
c
c     f(x) ~= sum(i) Ti(x)*ci
c
c     This calls the LAPACK routine dgels to use QR or a related
c     method ... DON'T use the normal equations ... I tried it
c     and it is too unstable for the present purpose.
c
      integer maxn, maxm, lwork
      parameter (maxn=64,maxm=100,lwork=10000)
      double precision at((maxn+1)*maxm),
     $     b((maxm+1)*(maxm+1)),h, x, work(lwork)
      integer i, j, info
c
      if (n.gt.maxn .or. m.gt.maxm) call errquit('tn_f_m: n or m ',
     &       UNKNOWN_ERR,
     $     n*10000+m)
c
      h = 2.0d0/(m-1)
      do i = 1, m
         x = -1.0d0 + h*(i-1)
         call tn_evaluate(n,x,at(1+(i-1)*(n+1)))
      enddo
c
      call dfill(m*m,0.0d0,b,1)
      call dfill(m,1.0d0,b,m+1)
      call dgels('t',n+1,m,m,at,n+1,b,m,work,lwork,info)
      if (info .ne. 0) call errquit('tn_f_m: dgels failed ', info,
     &       UNKNOWN_ERR)
c
      do i = 1, m
         do j = 0,n
            q(j,i) = b(1 + j + (i-1)*m)
         enddo
      enddo
c
      end
      subroutine tn_evaluate(n,x,t)
      implicit none
      double precision x
      integer n
      double precision t(0:n)
c
c     Compute the Chebyshev polynomials defined on -1,1 at the
c     point x.  n must be at least 1.
c
      integer j
c
      t(0)=1.0d0
      t(1)=x
      if (n.gt.1) then
         do j=1,n-1
            t(j+1) = 2.0d0*x*t(j) - t(j-1)
         enddo
      endif
c
      end
      subroutine tn_interp_2d(nx, ny, xlo, xhi, ylo, yhi, nptx, npty, 
     $     x, y, f, ldx, ff, ldxff, order)
      implicit none
#include "errquit.fh"
c     
c     Given a discretization of the area (xlo:xhi,ylo:yhi)
c     in f(1:nx,1:ny) return in ff(1:nptx,1:npty)
c     the values of f() interpolated onto the coordinates 
c     (x(i),y(j)) i=1,nptx, j=1,npty using a least-squares
c     Chebyshev approximation of the given order in each dimension.
c     
c     The input array f() is DESTROYED and must be dimensioned
c     ldx>=max(nx,ptx)
c     
      integer nx, ny, nptx, npty, ldx, ldxff, order
      double precision xlo,xhi,ylo,yhi
      double precision f(ldx, *), ff(ldxff,*)
      double precision x(nptx), y(npty)
c     
      integer i, j, l
      integer maxnpt, maxn
      parameter (maxnpt=100, maxn=100)
      double precision c(maxnpt,maxn), fij, tmp(maxnpt)
c     
      if (nx.gt.ldx) call errquit('tn_i_2d: nx>ldx',nx*10000+ldx,
     &       UNKNOWN_ERR)
      if (nptx.gt.ldx) call errquit('tn_i_2d: nptx>ldx',nptx*10000+ldx,
     &       UNKNOWN_ERR)
      if (nx.gt.maxn .or. nptx.gt.maxnpt) call errquit
     $     ('tn_i_2d: nx or nptx', 10000*nx+nptx,
     &       UNKNOWN_ERR)
      if (ny.gt.maxn .or. npty.gt.maxnpt) call errquit
     $     ('tn_i_2d: ny or npty', 10000*ny+npty,
     &       UNKNOWN_ERR)
c     
      do l = 1, nptx
         if (x(l).lt.xlo .or. x(l).gt.xhi) call errquit
     $        ('tn_interp_2d: x is out of range',l,
     &       UNKNOWN_ERR)
      enddo
      do l = 1, npty
         if (y(l).lt.ylo .or. y(l).gt.yhi) call errquit
     $        ('tn_interp_2d: y is out of range',l,
     &       UNKNOWN_ERR)
      enddo
c     
      call tn_collate_matrix(nx, xlo, xhi, nptx, x, order, c, maxnpt)
      do j = 1, ny
         do l = 1, nptx
            tmp(l) = 0.0d0
         end do
         do i = 1, nx
            fij = f(i,j)
            if (abs(fij) .gt. 0.0d0) then
               do l = 1, nptx
                  tmp(l) = tmp(l)+c(l,i)*fij
               end do
            end if
         end do
         do l = 1, nptx
            f(l,j) = tmp(l)
         end do
      end do
c     
      call tn_collate_matrix(ny, ylo, yhi, npty, y, order, c, maxnpt)
      do i = 1, nptx
         do l = 1, npty
            tmp(l) = 0.0d0
         end do
         do j = 1, ny
            fij = f(i,j)
            if (abs(fij) .gt. 0.0d0) then
               do l = 1, npty
                  tmp(l) = tmp(l) + c(l,j)*fij
               end do
            end if
         end do
         do l = 1, npty
            ff(i,l) = tmp(l)
         end do
      end do
c     
      end
      double precision function tn_interp_3d_point(
     $     g, nx, ny, nz, hx, hy, hz,
     $     x, y, z, 
     $     n, order, qq)
      implicit none
#include "errquit.fh"
c
c     Given a function tablulated on a grid with given spacing return
c     the value interpolated at (x,y,z) within the cube 
c     (0:hx*(nx-1),0:hy*(ny-1),0:hz*(nz-1))
c     a Chebyshev fit of the given order to the given number of points
c     in all dimensions. 
c
c     Note that the exterior points on the grid are ommited so
c     the grid is numbered rather oddly.  See solver.F.
c
c     The input array of points is PRESERVED.
c     
      integer nx, ny, nz, n, order
      double precision x, y, z, g(nx,ny,nz), hx, hy, hz, qq(*)
c     
      integer i, j, k, ind, ilo, jlo, klo, nhalf, jtop
      integer maxn, maxorder
      parameter (maxn=63, maxorder=53)
      double precision cx(maxn), cy(maxn), cz(maxn), sumi, sumj, sum,
     $     sumi0, sumi1, sumi2, sumi3
      double precision xx, yy, zz, twornm1,
     $     tx(0:maxorder), ty(0:maxorder), tz(0:maxorder)
c     
      if (n.gt.maxn) call errquit('tn_i_3d_pt: n', n, UNKNOWN_ERR)
      if (order.gt.maxorder) call errquit('tn_i_3d_pt:order',order,
     &       UNKNOWN_ERR)
c
      nhalf = n/2
      twornm1 = 2.0d0/dble(n-1)
      xx = x/hx                 ! Coords within whole volume
      yy = y/hy
      zz = z/hz
      ilo = max(1, int(xx)-nhalf) ! Corner of interpolation volume
      ilo = min(ilo,nx-n+1)
      jlo = max(1, int(yy)-nhalf)
      jlo = min(jlo,ny-n+1)
      klo = max(1, int(zz)-nhalf)
      klo = min(klo,nz-n+1)
      xx = (xx-dble(ilo))*twornm1 - 1.0d0 ! Rescale to [-1,1]
      yy = (yy-dble(jlo))*twornm1 - 1.0d0
      zz = (zz-dble(klo))*twornm1 - 1.0d0
c     
      tx(0)=1.0d0
      ty(0)=1.0d0
      tz(0)=1.0d0
      tx(1)=xx
      ty(1)=yy
      tz(1)=zz
      do j=1,order-1
         tx(j+1) = 2.0d0*xx*tx(j) - tx(j-1)
         ty(j+1) = 2.0d0*yy*ty(j) - ty(j-1)
         tz(j+1) = 2.0d0*zz*tz(j) - tz(j-1)
      enddo
c
      ind = 1
      do j = 1, n
         cx(j) = 0.0d0
         cy(j) = 0.0d0
         cz(j) = 0.0d0
         do k = 0,order
            cx(j) = cx(j) + tx(k)*qq(ind+k)
            cy(j) = cy(j) + ty(k)*qq(ind+k)
            cz(j) = cz(j) + tz(k)*qq(ind+k)
         enddo
         ind = ind + order + 1
      enddo
c
      ilo = ilo - 1
      jlo = jlo - 1
      klo = klo - 1
      sum = 0.0d0
c
c     Manually unroll the j loop into the i loop
c
      jtop = (n/4)*4
      do k = 1, n
         sumj = 0.0d0
         do j = 1, jtop, 4
            sumi0 = 0.0d0
            sumi1 = 0.0d0
            sumi2 = 0.0d0
            sumi3 = 0.0d0
            do i = 1, n
               sumi0 = sumi0 + cx(i)*g(i+ilo,j+jlo  ,k+klo)
               sumi1 = sumi1 + cx(i)*g(i+ilo,j+jlo+1,k+klo)
               sumi2 = sumi2 + cx(i)*g(i+ilo,j+jlo+2,k+klo)
               sumi3 = sumi3 + cx(i)*g(i+ilo,j+jlo+3,k+klo)
            end do
            sumj = sumj + cy(j)*sumi0 + cy(j+1)*sumi1 + cy(j+2)*sumi2 +
     $           cy(j+3)*sumi3
         end do
         do j = jtop+1,n
            sumi = 0.0d0
            do i = 1, n
               sumi = sumi + cx(i)*g(i+ilo,j+jlo,k+klo)
            end do
            sumj = sumj + cy(j)*sumi
         end do
         sum = sum + cz(k)*sumj
      end do
c$$$      do k = 1, n
c$$$         sumj = 0.0d0
c$$$         do j = 1, n
c$$$            sumi = 0.0d0
c$$$            do i = 1, n
c$$$               sumi = sumi + cx(i)*g(i+ilo,j+jlo,k+klo)
c$$$            end do
c$$$            sumj = sumj + cy(j)*sumi
c$$$         end do
c$$$         sum = sum + cz(k)*sumj
c$$$      end do
c     
      tn_interp_3d_point = sum
c
      end
      subroutine tn_lsq_fit(m,nfunc,order,f)
      implicit none
#include "errquit.fh"
      integer m, nfunc, order
      double precision f(m,nfunc)
c
c     Input is a set of nfunc functions on a uniform grid (1:m) 
c     tablulated in f(1:m,1:nfunc).
c
c     Return in f(1:order+1,1:nfunc) a least-squares fit to
c     Chebyshev polynomials up to the specified order.
c
c     order+1 <= m
c
      integer maxorder, maxm, lwork
      parameter (maxorder=35,maxm=63,lwork=10000)
      double precision at((maxorder+1)*maxm), h, x, work(lwork)
      integer i, info
c
      if (order+1.gt.m) call errquit('order+1>m',order,
     &       UNKNOWN_ERR)
      if (order.gt.maxorder) call errquit('order ', order,
     &       UNKNOWN_ERR)
      if (m .gt. maxm) call errquit('m', m,
     &       UNKNOWN_ERR)
c
      h = 2.0d0/(m-1)
      do i = 1, m
         x = -1.0d0 + h*dble(i-1)
         call tn_evaluate(order,x,at(1+(i-1)*(order+1)))
      enddo
c
      call dgels('t',order+1,m,nfunc,at,order+1,f,m,work,lwork,info)
      if (info .ne. 0) call errquit('info ', info,
     &       UNKNOWN_ERR)
c
*      write(6,*) ' solution from dgels'
*      call doutput(f,1,order+1,1,nfunc,order+1,nfunc,1)
c
      end
      subroutine tn_lsq_fit_cube(m,order,f)
      implicit none
      integer m, order
      double precision f(m,m,m)
c
c     Input in f(1:m,1:m,1:m) is a function tabulated in a cube.  
c
c     Return in f(1:order+1,1:order+1,1:order+1) a least squares
c     fit to Chebychev polynomials of the given order.
c     Also return in r the residual.
c
      integer maxorder, maxm
      parameter (maxorder=33, maxm=63)
      double precision g(maxm*maxm*maxm)
      integer i, j, k, op1, ind
c
      op1 = order+1
c
c     Transform the first dimension
c
      call tn_lsq_fit(m,m*m,order,f)
c
c     Rotate the dimensions and do the next
c
      call rotate_dims(f,m,m,g,m,op1,op1,m,m)
      call tn_lsq_fit(m,(op1)*m,order,g)
c
c     And the last one
c
      call rotate_dims(g,m,op1,f,m,op1,op1,op1,m)
      call tn_lsq_fit(m,(op1)**2,order,f)
c
c     Rotate again and copy back
c
      call rotate_dims(f,m,op1,g,op1,op1,op1,op1,op1)
      ind = 1
      do k = 1, op1
         do j = 1, op1
            do i = 1, op1
               f(i,j,k) = g(ind)
               ind = ind + 1
            enddo
         enddo
      enddo
c
      end
      subroutine rotate_dims(f, ldf1, ldf2, g, ldg1, ldg2,
     $     ni, nj, nk)
      implicit none
      integer ldf1, ldf2, ldg1, ldg2, ni, nj, nk
      double precision f(ldf1,ldf2,*), g(ldg1,ldg2,*)
c
c     Copy f into g permuting the dimensions so that g(k,i,j) = f(i,j,k)
c
      integer i, j, k
c
      do k = 1, nk
         do j = 1, nj
            do i = 1, ni
               g(k,i,j) = f(i,j,k)
            enddo
         enddo
      enddo
c
      end
      double precision function tn_cube_eval(g,ldg1,ldg2,order,x,y,z)
      implicit none
      integer ldg1, ldg2, order
      double precision g(ldg1,ldg2,*), x, y, z
c
c     Given in g a Chebychev fit of given order of a function defined
c     in a cube [-1,1] compute the value at the interior point x,y,z
c
      integer maxorder
      parameter (maxorder = 33)
      double precision tx(maxorder+1), ty(maxorder+1), tz(maxorder+1)
      double precision sum, sumi, sumj
      integer i, j, k, op1
c
      tx(1)=1.0d0
      ty(1)=1.0d0
      tz(1)=1.0d0
      tx(2)=x
      ty(2)=y
      tz(2)=z
      do j=2,order
         tx(j+1) = 2.0d0*x*tx(j) - tx(j-1)
         ty(j+1) = 2.0d0*y*ty(j) - ty(j-1)
         tz(j+1) = 2.0d0*z*tz(j) - tz(j-1)
      enddo
*      call tn_evaluate(order,x,tx)
*      call tn_evaluate(order,y,ty)
*      call tn_evaluate(order,z,tz)
c
      op1 = order+1
      sum = 0.0d0
      do k = 1, op1
         sumj = 0.0d0
         do j = 1, op1
            sumi = 0.0d0
            do i = 1, op1
               sumi = sumi + tx(i)*g(i,j,k)
            end do
            sumj = sumj + ty(j)*sumi
         end do
         sum = sum + tz(k)*sumj
      end do
c
      tn_cube_eval = sum
c
      end
