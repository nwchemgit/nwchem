      subroutine fastj(rtdb, basis, g_dens, eps, data)
      implicit none
#include "rtdb.fh"
#include "bas.fh"
#include "geom.fh"
#include "mafdecls.fh"
      integer rtdb
      integer basis
      integer g_dens
      double precision eps
      integer data
c
c     Given a normalized basis set and a density matrix
c     solve Poisson's equation so that subsequent calls to
c     fastj_evaluate() can be used to evaluate the potential
c     at arbitrary interior points.  
c
c     eps controls the precision in an as yet unspecified fashion
c
c     datais returned and must be passed unchanged to fastj_evaluate()
c     
c     call fastj_tidy(info) ... used to free allocated resources
c     (MA heap and GA).
c     
      integer iat, ish, i, nsh, nbf, geom, nat, ishlo, ishhi, type,
     $     nprim, ngen, n, nx, ny, nz, l_grid, k_grid,
     $     depth, lmax
      integer maxprim, maxgen, maxgrid
      parameter (maxprim = 40, maxgen=1, maxgrid=16)
      double precision alpha(maxprim)
      logical osph
      character*16 tag
      integer gridinfo(5,maxgrid)
      double precision gridh(3,maxgrid)
      integer ngrid, nnx, nny, nnz, l_poles, k_poles, lenpoles
      double precision c(3), charge, x, y, z, xlo, xhi, ylo, yhi,
     $     zlo, zhi, range, side, resolution, pi, exponent,
     $     hx, hy, hz, hx_fmm, hy_fmm, hz_fmm
      double precision fastj_gaussian_range
      external fastj_gaussian_range
c
c     Initialize radial fits of potentials due to gaussians
c     and also the multipoles
c
      call anl_fit_init
      call xlm_init
c
c     Get info
c
      pi = 4.0d0*atan(1.0d0)
      if (.not. bas_numcont(basis, nsh)) call errquit('fastJ: basis?',0)
      if (.not. bas_numbf(basis, nbf)) call errquit('fastJ: basis?',0)
      if (.not. bas_geom(basis, geom)) call errquit('fastJ: basis?',0)
      if (.not. geom_ncent(geom, nat)) call errquit('fastJ: geom?',0)
c
      write(6,76) eps
 76   format(' eps = ',1p,d9.1)
      do i = -7,7               ! alpha = 2**i
         range = fastj_gaussian_range(0, 2.0d0**i, eps)
         exponent = 1.0d0/ (4.0d0 * 2.0d0**i)
         resolution = fastj_gaussian_range(0, exponent, eps)/pi
         write(6,77) i, 2.0d0**i, range, resolution
 77      format(1x,i2,' alpha=',f8.4,' range=',f8.4,' resolution=',f8.4)
      end do
c
c     Determine the solution volume.  It must enlose the space of
c     all basis function products.  Look at the range of the square
c     of functions on each center.
c
      do iat = 1, nat
         if (.not. bas_ce2cnr(basis, iat, ishlo, ishhi))
     $        call errquit('fastJ: bas_ce2cnr', 0)
         if (.not. geom_cent_get(geom, iat, tag, c, charge))
     $        call errquit('fastJ: geom_cart_get ', iat)
         x = c(1)
         y = c(2)
         z = c(3)
         if (iat .eq. 1) then
            xlo = x
            xhi = x
            ylo = y
            yhi = y
            zlo = z
            zhi = z
         end if
         do ish = ishlo, ishhi
            if (.not.bas_continfo(basis, ish, type, nprim, ngen, osph))
     $           call errquit('fastJ: bas_cont_info', ish)
            type = abs(type)
            if (nprim.gt.maxprim .or. ngen.gt.maxgen)
     $           call errquit('fastJ: nprim or ngen?',nprim*100+ngen)
            if (.not.bas_get_exponent(basis, ish, alpha))
     $           call errquit('fastJ: bas_get_exponent?',ish)
c
            do i = 1, nprim
               range = fastj_gaussian_range(2*type, 2.0d0*alpha(i), eps)
               xlo = min(x-range,xlo)
               xhi = max(x+range,xhi)
               ylo = min(y-range,ylo)
               yhi = max(y+range,yhi)
               zlo = min(z-range,zlo)
               zhi = max(z+range,zhi)
            end do
         end do
      end do
c
c     Eventually the FMM box size will be determined by the
c     range of the functions that it must contain.  For now
c     set the depth to be fixed.  Similarly, lmax will be
c     set from the required precision.  For now set to 16.
c
      depth = 3
      lmax = 16
c
c     The FMM handles non-cubic volumes except for the caching
c     of rotation operators.  Until this is fixed use a cubic volume
c
c     Also, adjust the domain size so that we have one FMM box on
c     either side so that we have accurate potentials everywhere
c     the basis functions are non-zero
c
      side = max(xhi-xlo,yhi-ylo,zhi-zlo) * 0.5d0
c
      side = side * 2.0d0**depth / (2.0d0**depth - 2.0d0)
      x = (xhi+xlo)/2.0d0
      y = (yhi+ylo)/2.0d0
      z = (zhi+zlo)/2.0d0
      xlo = x - side
      xhi = x + side
      ylo = y - side
      yhi = y + side
      zlo = z - side
      zhi = z + side
c      
      write(6,1) xlo, xhi, ylo, yhi, zlo, zhi
 1    format(' Solution volume: ',3(2f7.2,'  '))
c
c     In this initial version all functions will be evaluated on
c     the grid.  Determine the necessary resolution ... must be
c     able to represent the square of each function.  Again,
c     resolution could be decreased by scaling with the density etc.
c
      resolution = 0.0d0        ! No. of grid points per unit
      do iat = 1, nat
         if (.not. bas_ce2cnr(basis, iat, ishlo, ishhi))
     $        call errquit('fastJ: bas_ce2cnr', 0)
         if (.not. geom_cent_get(geom, iat, tag, c, charge))
     $        call errquit('fastJ: geom_cart_get ', iat)
         do ish = ishlo, ishhi
            if (.not.bas_continfo(basis, ish, type, nprim, ngen, osph))
     $           call errquit('fastJ: bas_cont_info', ish)
            type = abs(type)
            if (nprim.gt.maxprim .or. ngen.gt.maxgen)
     $           call errquit('fastJ: nprim or ngen?',nprim*100+ngen)
            if (.not.bas_get_exponent(basis, ish, alpha))
     $           call errquit('fastJ: bas_get_exponent?',ish)
            do i = 1, nprim
               n = 2 * type
               exponent = 1.0d0 / (8.0d0 * alpha(i))
               resolution = max(resolution,
     $              fastj_gaussian_range(n, exponent, eps)/pi)
            end do
         end do
      end do
c
c     Now determine how many grid points we will have in each direction.
c     Assume that the FFT can handle radix 2,3,5 so round to the nearest
c     power of these.  Also, since we are using a sine transform the
c     zeroes at the boundary are not represented on the grid.
c     Also, using sint from FFTPACK we must have n+1 = power of radix.
c
c     Since we are temporaily restricted to cubic only do one dim
c
      n = int(2.0d0*side*resolution-1.0d0) + 1
      write(6,*) ' raw n ', n
      call fastj_optimize_n(n, 'upperbound')
      nx = n
      ny = n
      nz = n
      hx  = (xhi-xlo)/dble(nx+1)
      hy  = (yhi-ylo)/dble(ny+1)
      hz  = (zhi-zlo)/dble(nz+1)
c     
      hx_fmm = (xhi-xlo)/dble(2**depth) ! FMM box size
      hy_fmm = (yhi-ylo)/dble(2**depth)
      hz_fmm = (zhi-zlo)/dble(2**depth)
c
      write(6,2) resolution, nx, ny, nz
 2    format(' FFT grid points per unit ', f7.2/
     $       '       No. of grid points ', 3i5)
c
c     Allocate the grid and tabulate the smooth component
c     of the density ... eventually this must be on a
c     nest of grids and we must compute the multipoles
c     and build the direct lists at the same time.
c
c     For efficient evaluation we allocate a nest of grids
c     so that we don't oversample smooth functions.
c
      nnx = nx
      nny = ny
      nnz = nz
      do ngrid = 1, maxgrid
         if (.not. ma_push_get(mt_dbl, nnx*nny*nnz, 'fastj: grid',
     $        l_grid, k_grid)) 
     $        call errquit('fastj: ma grid',nnx*nny*nnz)
         call dfill(nnx*nny*nnz, 0.0d0, dbl_mb(k_grid), 1)
c
         gridinfo(1,ngrid) = nnx
         gridinfo(2,ngrid) = nny
         gridinfo(3,ngrid) = nnz
         gridinfo(4,ngrid) = k_grid
         gridinfo(5,ngrid) = l_grid
         gridh(1,ngrid) = (xhi-xlo)/dble(nnx+1)
         gridh(2,ngrid) = (xhi-xlo)/dble(nny+1)
         gridh(3,ngrid) = (xhi-xlo)/dble(nnz+1)
c
         write(6,81) ngrid,
     $        (gridinfo(i,ngrid),i=1,3),(gridh(i,ngrid),i=1,3)
 81      format(1x,i2,' Allocated grid(',i4,',',i4,',',i4,')  with h(',
     $        f6.2,',',f6.2,',',f6.2,')')
c         
         nnx = (nnx+1)/1.26 - 1
         nny = (nny+1)/1.26 - 1
         nnz = (nnz+1)/1.26 - 1
         call fastj_optimize_n(nnx, 'nearby')
         call fastj_optimize_n(nny, 'nearby')
         call fastj_optimize_n(nnz, 'nearby')
c
         if (nnx.lt.32 .or. nny.lt.32 .or. nnz.lt.32) goto 10
      enddo
 10   if (ngrid .gt. maxgrid) ngrid = maxgrid
c
      call fastj_density(ngrid, gridinfo, gridh,
     $     xlo, xhi, ylo, yhi, zlo, zhi,
     $     basis, geom, nat, nsh, nbf, g_dens, eps)
c
c     Free up all of the memory we are no longer using
c
      if (.not. ma_chop_stack(gridinfo(5,2))) call errquit
     $     ('fastJ: ma chop failed ', 1)
c
c     Only finest grid (no. 1) left now
c
      nx = gridinfo(1,1)
      ny = gridinfo(2,1)
      nz = gridinfo(3,1)
      k_grid = gridinfo(4,1)
      hx = gridh(1,1)
      hy = gridh(2,1)
      hz = gridh(3,1)
c
c     Allocate the FMM boxes ... just zero for now
c
      lenpoles = (lmax+1)**2 * (2**depth)**3
      if (.not. ma_push_get(mt_dbl, lenpoles, 'fastj: grid',
     $     l_poles, k_poles)) 
     $     call errquit('fastj: ma grid',nnx*nny*nnz)
      call dfill(lenpoles, 0.0d0, dbl_mb(k_poles), 1)
      if (.not. ma_verify_allocator_stuff()) call errquit('fred',0)
c
c     Solve the sucker
c
      call solver(
     $     dbl_mb(k_grid), nx, ny, nz, 
     $     dbl_mb(k_poles), depth, lmax,
     $     xlo, xhi, ylo, yhi, zlo, zhi)
c     
      if (.not. ma_chop_stack(gridinfo(5,1))) call errquit
     $     ('fastJ: ma chop failed ', 0)
c
      end
      subroutine fastj_density(ngrid, gridinfo, gridh,
     $     xlo, xhi, ylo, yhi, zlo, zhi,
     $     basis, geom, nat, nsh, nbf, g_dens, eps)
      implicit none
#include "mafdecls.fh"
#include "bas.fh"
#include "geom.fh"
      integer ngrid
      integer gridinfo(5,ngrid)
      double precision gridh(3,ngrid)
      double precision xlo, xhi, ylo, yhi, zlo, zhi
      integer basis, geom, nat, nsh, nbf
      integer g_dens
      double precision eps
c
c     Parameters
c
      integer maxprim           ! Max. no. primitives in a contraction
      parameter (maxprim=25)
      integer lmax, lmax2, lencart, lenprod, maxn
      parameter (lmax = 3)      ! Max. angular momentum of shell
      parameter (lmax2 = 2*lmax)
      parameter (lencart = ((lmax2+1)*(lmax2+2))/2)
      parameter (lenprod = ((lmax2+1)*(lmax2+2)*(lmax2+3))/6)
      parameter (maxn = 1024)   ! Max. grid dimension
c
      double precision d(lencart,-lmax2:lmax2,0:lmax2),
     $     dinv(lenprod, -lmax2:lmax2,0:lmax2)
      double precision work(lenprod**2), ndens(lenprod**2), 
     $     dens(lencart**2)
      double precision xi, yi, zi, xj, yj, zj, s
      double precision expi(maxprim), coeffi(maxprim)
      double precision expj(maxprim), coeffj(maxprim)
      double precision factor, charge, coords(3)
      double precision a, b, c
      double precision alpha, prefac, tol, dmax, gmin, pi
      double precision h, hx, hy, hz, resolution, scale
c
      character*16 tag
c
      integer ish,iat,ishlo,ishhi,itype,inprim,ingen,numi,ibflo,ibfhi
      integer jsh,jat,jshlo,jshhi,jtype,jnprim,jngen,numj,jbflo,jbfhi
      integer i, iprim 
      integer j, jprim 
      integer k, nx, ny, nz, nnx, nny, nnz, jk, jjkk
      integer ishtop, ijtype, igrid, k_grid, kk_grid
c
      logical iosph, josph
c
      double precision fastj_gaussian_range
      external fastj_gaussian_range
c     
      call xlm_coeff_inv(lmax2,d,dinv)
      pi = 4.0d0*atan(1.0d0)
c
*      write(6,*) ' Input density '
*      call ga_print(g_dens)
c
      tol = eps * 1d-3          ! Screen very tight on prefactor
c
      if (nx.gt.maxn .or. ny.gt.maxn .or. nz.gt.maxn) call errquit
     $     ('fastj: hard dimension of maxn exceeded ', maxn)
c
      do jat = 1, nat
         if (.not. bas_ce2cnr(basis, jat, jshlo, jshhi))
     $        call errquit('fastJ dens: bas_ce2cnr', jat)
         if (.not. geom_cent_get(geom, jat, tag, coords, charge))
     $        call errquit('fastJ: geom_cart_get ', jat)
         xj = coords(1)
         yj = coords(2)
         zj = coords(3)
*         write(6,*) ' JAT: tag,jat,xj,yj,zj,jshlo,jshhi ', 
*     $        tag, '  ', jat,xj, yj, zj, jshlo, jshhi
         do iat = 1, jat
            if (.not. bas_ce2cnr(basis, iat, ishlo, ishhi))
     $           call errquit('fastJ dens: bas_ce2cnr', iat)
            if (.not. geom_cent_get(geom, iat, tag, coords, charge))
     $           call errquit('fastJ: geom_cart_get ', iat)
            xi = coords(1)
            yi = coords(2)
            zi = coords(3)
*            write(6,*) ' IAT: tag,iat,xi,yi,zi,ishlo,ishhi ', 
*     $           tag, '  ', iat,xi, yi, zi, ishlo, ishhi
            do jsh = jshlo, jshhi
               if (.not.bas_continfo(basis, jsh, jtype, jnprim,jngen, 
     $              josph)) call errquit('fastJ: bas_cont_info', jsh)
               if (jtype .lt. 0 .or. jngen.ne.1 .or. josph) call errquit
     $              ('fastJ: no spherical/general/SP functions yet',0)
               if (jtype .gt. lmax) call errquit
     $              ('fastJ: hard dimension lmax exceeded ', lmax)
               if (.not. bas_cn2bfr(basis, jsh, jbflo, jbfhi))
     $              call errquit('fastJ dens: bas_cn2bfr', jsh)
               if (.not.bas_get_exponent(basis, jsh, expj))
     $              call errquit('fastJ: bas_get_exponent?',jsh)
               if (.not.bas_get_coeff(basis, jsh, coeffj))
     $              call errquit('fastJ: bas_get_coeff?',jsh)
               numj = jbfhi - jbflo + 1
c
*               write(6,*) ' JSH: jsh,jtype,jnprim,jbflo,jbfhi,numj ', 
*     $              jsh, jtype, jnprim, jbflo, jbfhi, numj
*               write(6,*) ' Exponents    ',(expj(j),j=1,jnprim)
*               write(6,*) ' Coefficients ',(coeffj(j),j=1,jnprim)
c
               ishtop = ishhi
               if (iat .eq. jat) ishtop = jsh
               do ish = ishlo, ishtop
                  if (.not.bas_continfo(basis,ish,itype,inprim,ingen, 
     $                 iosph)) call errquit('fastJ: bas_cont_info',ish)
                  if (itype .gt. lmax) call errquit
     $                 ('fastJ: hard dimension lmax exceeded ', lmax)
                  if (.not. bas_cn2bfr(basis, ish, ibflo, ibfhi))
     $                 call errquit('fastJ dens: bas_cn2bfr', ish)
                  if (.not.bas_get_exponent(basis, ish, expi))
     $                 call errquit('fastJ: bas_get_exponent?',ish)
                  if (.not.bas_get_coeff(basis, ish, coeffi))
     $                 call errquit('fastJ: bas_get_coeff?',ish)
                  numi = ibfhi - ibflo + 1
c
                  ijtype = itype + jtype
c
*               write(6,*) ' ISH: ish,itype,inprim,ibflo,ibfhi,numi ', 
*     $              ish, itype, inprim, ibflo, ibfhi, numi
*               write(6,*) ' Exponents    ',(expi(i),i=1,inprim)
*               write(6,*) ' Coefficients ',(coeffi(i),i=1,inprim)
c
c
c     Copy the block of density into local workspace
c
                  factor = 2.0d0 ! This assumes symmetric density!
                  if (ish .eq. jsh) factor = 1.0d0
                  call ga_get(g_dens,ibflo,ibfhi,jbflo,jbfhi,dens,numi)
                  dmax = 0.0d0
                  do i = 1, numi*numj
                     dmax = max(dmax,abs(dens(i)))
                  enddo
c
*                  write(6,*) ' IJ density ', dmax
*                  call output(dens,1,numi,1,numj,numi,numj,1)
c
                  do jprim = 1, jnprim
                     do iprim = 1, inprim
                        call gaussian_product(
     $                       expi(iprim), xi, yi, zi,
     $                       expj(jprim), xj, yj, zj,
     $                       alpha, a, b, c, prefac)
c
*                        write(6,*) ' PRIM: jp,ip,ei,ej,alpha,a,b,c,',
*     $                       'pre,fac,ci,cj ', jprim,iprim,
*     $                       expi(iprim), expj(jprim), alpha, a,b,c,
*     $                       prefac, factor, coeffi(iprim),coeffj(jprim)
c
                        prefac = prefac*factor*
     $                       coeffi(iprim)*coeffj(jprim) 
c
                        if (abs(prefac*dmax) .gt. tol) then
                           call cart_dens_translate(
     $                          itype, xi, yi, zi, 
     $                          jtype, xj, yj, zj,
     $                          a, b, c, dens, work, ndens)
                           call cart_dens_product(
     $                          itype, jtype, work, ndens)
c
c     Pick the most appropriate grid
c
                           resolution = fastj_gaussian_range(ijtype,
     $                          0.25d0/alpha, eps)/pi
                           h = 1.0d0/resolution
                           do igrid = ngrid,1,-1
                              if (h .ge. gridh(1,igrid)) goto 333
                           enddo
                           call errquit('NO GRID?',0)
c     
 333                       nx = gridinfo(1,igrid)
                           ny = gridinfo(2,igrid)
                           nz = gridinfo(3,igrid)
                           k_grid = gridinfo(4,igrid)
                           hx = gridh(1,igrid)
                           hy = gridh(2,igrid)
                           hz = gridh(3,igrid)
**                           write(6,*) ' evaluating on grid ', igrid
c
                           call fastj_density_eval(
     $                          ndens,
     $                          ijtype, alpha, prefac, a, b, c, 
     $                          dbl_mb(k_grid), nx, ny, nz, hx, hy, hz,
     $                          xlo, xhi, ylo, yhi, zlo, zhi, tol)
                        end if
                     end do      ! do iprim
                  end do         ! do jprim
               end do            ! do ish
            end do               ! do jsh
         end do                  ! do iat
      end do                     ! do jat
c
c     Interpolate all of the grids down to the finest using the FFT
c
      do igrid = ngrid,1,-1
         nx = gridinfo(1,igrid)
         ny = gridinfo(2,igrid)
         nz = gridinfo(3,igrid)
         k_grid = gridinfo(4,igrid)
         call vsint3d(nx,ny,nz,dbl_mb(k_grid))
         scale = 0.125d0/(dble(nx+1)*dble(ny+1)*dble(nz+1))
         call dscal(nx*ny*nz,scale,dbl_mb(k_grid),1)
      enddo
      nnx = gridinfo(1,1)
      nny = gridinfo(2,1)
      nnz = gridinfo(3,1)
      kk_grid = gridinfo(4,1)
      do igrid = 2, ngrid       ! Grid 1 is the finest
         nx = gridinfo(1,igrid)
         ny = gridinfo(2,igrid)
         nz = gridinfo(3,igrid)
         k_grid = gridinfo(4,igrid)
         do k = 1, nz
            do j = 1, ny
               jk   =  k_grid +  nx*(j-1 +  ny*(k-1)) - 1
               jjkk = kk_grid + nnx*(j-1 + nny*(k-1)) - 1
               do i = 1, nx
                  dbl_mb(i+jjkk) = dbl_mb(i+jjkk) + dbl_mb(i+jk)
               enddo
            enddo
         enddo
      enddo
      call vsint3d(nnx,nny,nnz,dbl_mb(kk_grid)) ! can be optimized away
c
c     Debug
c
      s = 0.0d0
      gmin = 0.0d0
      do k = 1, nnz
         do j = 1, nny
            do i = 1, nnx
               gmin = min(gmin,dbl_mb(kk_grid))
*               if (g(i,j,k).lt.-eps) then
*                  write(6,*) ' negative ', i, j, k, g(i,j,k)
*                  call errquit('33',33)
*               end if
               s = s + dbl_mb(kk_grid)
               kk_grid = kk_grid + 1
            end do
         end do
      end do
      s = s * gridh(1,1)**3.0d0
      write(6,*) ' Largest negative element ', gmin
      write(6,*) ' Trapezoidal integration of the density ', s
c
      end
      subroutine fastj_density_eval(
     $     ndens,
     $     ijtype, alpha, prefac, a, b, c, 
     $     g, nx, ny, nz, hx, hy, hz,
     $     xlo, xhi, ylo, yhi, zlo, zhi, tol)
      implicit none
      double precision ndens(*)
      integer ijtype
      double precision alpha, prefac, a, b, c, tol
      integer nx, ny, nz
      double precision g(nx, ny, nz), hx, hy, hz
      double precision xlo, xhi, ylo, yhi, zlo, zhi
c     
      integer maxn, lmax2
      parameter (maxn=1024, lmax2=6)
      double precision xx(maxn,0:lmax2), yy(maxn,0:lmax2),
     $     zz(maxn,0:lmax2)
      double precision ndmax, range
      double precision x, y, z, xxmax, s, sk
      integer i, j, k, l, p, ii, jj, kk, ind
      integer ilo, ihi, jlo, jhi, klo, khi
      double precision fastj_gaussian_range
      external fastj_gaussian_range
c
      ndmax = 0.0d0
      do i = 1,((ijtype+1)*(ijtype+2)*(ijtype+3))/6
         ndmax = max(ndmax,abs(ndens(i)))
      enddo
c     
      range = fastj_gaussian_range(ijtype, alpha, 
     $     min(tol/abs(ndmax*prefac),1d-4))
c     
*     write(6,*) ' primitive range ', range
c     
c     Precompute powers of cartesians and gaussians on the grid
c     
      ilo = max(int((a-range-xlo)/hx),1)
      ihi = min(int((a+range-xlo)/hx)+1,nx)
      jlo = max(int((b-range-ylo)/hy),1)
      jhi = min(int((b+range-ylo)/hy)+1,ny)
      klo = max(int((c-range-zlo)/hz),1)
      khi = min(int((c+range-zlo)/hz)+1,nz)
c     
*     write(6,*) ' restricted grid ',
*     $                          ilo, ihi, jlo, jhi, klo, khi
c     
      xxmax = 0.0d0
      x = xlo + ilo*hx - a
      do i = ilo, ihi
         xx(i,0)  = exp(-alpha*x*x)
         xxmax = max(xxmax,abs(xx(i,0)))
         x = x + hx
      end do
      y = ylo + jlo*hy - b
      do j = jlo, jhi
         yy(j,0) = exp(-alpha*y*y)
         y = y + hy
      end do
      z = zlo + klo*hz - c
      do k = klo, khi
         zz(k,0) = exp(-alpha*z*z)
         z = z + hz
      end do
c     
      do p = 1, ijtype
         x = xlo + ilo*hx - a
         do i = ilo, ihi
            xx(i,p) = xx(i,p-1)*x
            xxmax = max(xxmax,abs(xx(i,p)))
            x = x + hx
         end do
         y = ylo + jlo*hy - b
         do j = jlo, jhi
            yy(j,p) = yy(j,p-1)*y
            y = y + hy
         end do
         z = zlo + klo*hz - c
         do k = klo, khi
            zz(k,p) = zz(k,p-1)*z
            z = z + hz
         end do
      end do
c     
c     Try to get good pipelining and cache utilization by vectorizing
c     down the x direction of the grid in the inner loop.  By
c     evaluating over the cube rather than the sphere we do about
c     2x as many evaluations, but they are much more efficient.
c     
      if (ijtype .eq. 0) then
         do k = klo, khi
            sk = prefac*ndens(1)*zz(k,0)
            do j = jlo, jhi
               s = sk*yy(j,0)
               if (abs(s*xxmax) .gt. tol) then
                  do i = ilo, ihi
                     g(i,j,k) = g(i,j,k) + s*xx(i,0)
                  end do
               end if
            end do
         end do
      else
         do k = klo, khi
            ind = 0 
            do l = 0, ijtype
               do ii = l,0,-1
                  do jj = l-ii,0,-1
                     kk = l-ii-jj
                     ind = ind + 1
                     sk = prefac*ndens(ind)*zz(k,kk)
                     do j = jlo, jhi
                        s = sk*yy(j,jj)
                        if (abs(s*xxmax) .gt. tol) then
                           do i = ilo, ihi
                              g(i,j,k) = g(i,j,k) + s*xx(i,ii)
                           end do
                        end if
                     end do
                  end do
               end do
            end do
         end do
      end if                    ! if (ijtype)
c     
      end
      double precision function fastj_gaussian_range(n, alpha, eps)
      implicit none
c
      integer n
      double precision alpha, eps
c
c     Return an approximation to the outer solution of 
c     .     r^n*exp(-ar^2) = eps
c     .     r = (n*ln(-ln(eps)) - n*ln(a) - 4*ln(eps)) /
c     .         4*sqrt(-alpha*ln(eps))
c
c     Accuracy improves with smaller eps.
c
      double precision logeps
c
      logeps = log(eps)
c      
      fastj_gaussian_range = 
     $     (n*log(-logeps) - n*log(alpha) - 4.0d0*logeps) /
     $     sqrt(-16.0d0*alpha*logeps)
c
      end

      double precision function fnlm(n,l,m,alpha,a,b,c,x,y,z)
      implicit none
c
      integer n, l, m
      double precision alpha, a, b, c, x, y, z
c
c     r^(n-l) Xlm exp(-a*r*r) centered about (a,b,c)
c     
      double precision rsq, r, radial, q(-10:10,0:10)
c
      rsq = (x-a)**2 + (y-b)**2 + (z-c)**2
      r   = sqrt(rsq)
c
      radial = exp(-alpha*rsq)
      if (n .ne. l) radial = radial * r**(n-l)
c
      call xlm(l, x-a, y-b, z-c, q, 10)
c
      fnlm = radial * q(m,l)
c
*      if (r .lt. 0.5d0) then
*         write(6,*) n, l, m, alpha, a, b, c, x, y, z, r
*         call xlm_print(l, q, 10)
*         write(6,*) radial, q(m,l), radial
*      end if
c
**      fnlm = exp(-alpha*rsq) * (z-c)
c      
      end
      double precision function vnlm(n,l,m,alpha,a,b,c,x,y,z)
      implicit none
c
      integer n, l, m
      double precision alpha, a, b, c, x, y, z
c
c     Return the potential due to 
c         r^(n-l) Xlm exp(-a*r*r) centered about (a,b,c)
c
c     A(sqrt(a)*r) * 4*pi/(a**((n+2)/2)) * Xlm / r^l
c
c     A(x) = (0.5*x**(n-l) (n+l+1)F((n+l)/2,x**2) + (n-l)**I(n-l-1)/2)
c     .         / (2l+1)
c
c     
      double precision rsq, r, radial, q(-10:10,0:10)
      double precision anl_fit
      external anl
c
      rsq = (x-a)**2 + (y-b)**2 + (z-c)**2
      r   = sqrt(rsq)
c
      radial = anl_fit(
     $     n,l,sqrt(alpha)*r) * 4.0d0 * 3.1415926535897932d0 /
     $     (alpha**(dble(n-l+2)/2.0d0))
c
      call xlm(l, x-a, y-b, z-c, q, 10)
      vnlm = radial * q(m,l)

**      vnlm = radial * (z-c)
c
      end
      subroutine fastj_optimize_n(n, mode)
      implicit none
      integer n
      character*(*) mode
c
c   Find a number of the form 2^i 3^j 5^k - 1 that is close to n.
c   Optionally force it to be an upperbound.
c
      integer n2, n3, n5, ii2, ii3, ii5, i2, i3, i5
      double precision test, diff
c
      n = n + 1                 ! Temporarily to force to power of radix
c
      n2 = log(dble(n))/log(2.0d0) + 1
      n3 = log(dble(n))/log(3.0d0) + 1
      n5 = log(dble(n))/log(5.0d0) + 1
      ii2 = 0
      ii3 = 0
      ii5 = 0
      diff = n
      do i2 = 0, n2
         do i3 = 0, n3
            do i5 = 0, n5
               test = (2.0d0**i2)*(3.0d0**i3)*(5.0d0**i5) + 1d-6
               if (test.ge.n .or. mode.ne.'upperbound') then
                  if (abs(test-n) .lt. diff) then
                     diff = abs(test - n)
                     ii2 = i2
                     ii3 = i3
                     ii5 = i5
                  endif
               endif
            enddo
         enddo
      enddo
c
      n = (2.0d0**ii2)*(3.0d0**ii3)*(5.0d0**ii5)  - 1
c
      end
