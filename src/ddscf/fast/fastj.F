      subroutine fastj(rtdb, basis, g_dens, eps, data)
*
* $Id$
*
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "geom.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "global.fh"
      integer rtdb
      integer basis
      integer g_dens
      double precision eps
      integer data(2)
c
c     Given a normalized basis set and a density matrix
c     solve Poisson's equation so that subsequent calls to
c     fastj_evaluate() can be used to evaluate the potential
c     at arbitrary interior points.  
c
c     eps controls the precision in an as yet unspecified fashion
c
c     data is returned and must be passed unchanged to fastj_evaluate()
c     
c     call fastj_tidy(data) ... used to free allocated resources
c     (MA heap and GA).
c     
      integer iat, ish, i, nsh, nbf, geom, nat, ishlo, ishhi, type,
     $     nprim, ngen, n, nx, ny, nz, l_grid, k_grid,
     $     depth, lmax, levelbc
      integer maxprim, maxgen, maxgrid, ntest, maxl
      parameter (maxprim = 31, maxgen=1, maxgrid=16,ntest=100000,
     $     maxl = FMM_LMAX)
      double precision alpha(maxprim), start, used
      double precision rtest(3,ntest), ptest(ntest), ptest2(ntest)
      double precision q(-maxl:maxl,0:maxl)
      logical osph, status
      character*16 tag
      integer gridinfo(5,maxgrid)
      double precision gridh(3,maxgrid)
      integer ngrid, nnx, nny, nnz, lenpoles
      double precision c(3), charge, x, y, z, xlo, xhi, ylo, yhi,
     $     zlo, zhi, range, side, resolution, pi, exponent, xx, yy, zz,
     $     hx, hy, hz, hx_fmm, hy_fmm, hz_fmm, 
*     $     test1, 
     $     test2, test3,
     $     errmax, errrelmax, errrelrms, errrms, err, rneglect, resn
      integer nbig, nprimtot
c
      integer lenpairs, lenxlm, l_list, k_list, l_pairs, k_pairs, 
     $     l_pairs2, k_pairs2, npairs, tree(300)
      double precision oversample, total_time
c
      double precision fastj_gaussian_range, potential
      double precision  fastj_r_neglected
      external fastj_gaussian_range, potential
      external fastj_r_neglected
c
c     Initialize radial fits of potentials due to gaussians
c     and also the multipoles
c
      call anl_fit_init
      call xlm_init
      total_time = util_cpusec()
c
      write(6,*) ' before free g_dens ', g_dens
      call ga_summarize(0)
      if (.not. ga_destroy(-1000)) call errquit('-1000',0, GA_ERR)
      if (.not. ga_destroy(-999)) call errquit('-999',0, GA_ERR)
      if (.not. ga_destroy(-998)) call errquit('-998',0, GA_ERR)
      if (.not. ga_destroy(-996)) call errquit('-996',0, GA_ERR)
      if (.not. ga_destroy(-995)) call errquit('-995',0, GA_ERR)
      call ga_summarize(0)
c
c     Get info
c
      pi = 4.0d0*atan(1.0d0)
      if (.not. bas_numcont(basis, nsh)) call errquit('fastJ: basis?',0,
     &       BASIS_ERR)
      if (.not. bas_numbf(basis, nbf)) call errquit('fastJ: basis?',0,
     &       BASIS_ERR)
      if (.not. bas_geom(basis, geom)) call errquit('fastJ: basis?',0,
     &       BASIS_ERR)
      if (.not. geom_ncent(geom, nat)) call errquit('fastJ: geom?',0,
     &       GEOM_ERR)
c
      if (util_print('fastjranges',print_debug)) then
         write(6,76) eps
 76      format(' Ranges/resolutions for accuracy = ',1p,d9.1)
         do i = -7,7            ! alpha = 2**i
            range = fastj_gaussian_range(0, 2.0d0**i, eps)
            rneglect= fastj_r_neglected(0, 2.0d0**i, eps)
            exponent = 1.0d0/ (4.0d0 * 2.0d0**i)
            resolution = fastj_gaussian_range(0, exponent, eps)/pi
            resn = fastj_r_neglected(0, exponent, eps)/pi
            write(6,77) i, 2.0d0**i, range, rneglect, resolution, resn
 77         format(1x,i2,' alpha=',f8.4,' range=',f8.4, ' rneg=',f8.4,
     $           ' res=',f8.4, ' resneg=',f8.4)
         end do
      end if
c
c     Determine the solution volume.  It must enlose the space of
c     all basis function products.  Look at the range of the square
c     of functions on each center.
c
      nprimtot = 0
      do iat = 1, nat
         if (.not. bas_ce2cnr(basis, iat, ishlo, ishhi))
     $        call errquit('fastJ: bas_ce2cnr', 0, BASIS_ERR)
         if (.not. geom_cent_get(geom, iat, tag, c, charge))
     $        call errquit('fastJ: geom_cart_get ', iat, GEOM_ERR)
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
     $           call errquit('fastJ: bas_cont_info', ish, BASIS_ERR)
            type = abs(type)
            if (nprim.gt.maxprim .or. ngen.gt.maxgen)
     $           call errquit('fastJ: nprim or ngen?',nprim*100+ngen,
     &       UNKNOWN_ERR)
            if (.not.bas_get_exponent(basis, ish, alpha))
     $           call errquit('fastJ: bas_get_exponent?',ish, BASIS_ERR)
c
            nprimtot = nprimtot + nprim
            do i = 1, nprim
               range = fastj_r_neglected(2*type,2.0d0*alpha(i),eps)
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
c     Sidemin is the minimum dimension of a cube that contains everything.
c     The actual volume will include a border of 1 FMM box
c
c     The FMM handles non-cubic volumes except for the caching
c     of rotation operators.  Until this is fixed use a cubic volume
c
c     Pick initial guesses for the FMM tree depth and FFT grid points
c     but allow full override by the user.
c
      if (.not. rtdb_get(rtdb, 'fastj:depth', mt_int, 1, depth))
     $     depth = 4
      if (.not. rtdb_get(rtdb, 'fastj:levelbc', mt_int, 1, levelbc))
     $     levelbc = depth
c
      if (.not. rtdb_get(rtdb, 'fastj:lmax', mt_int, 1, lmax)) then
         lmax = nint(-log(eps/3.26e-3)/log(1.7783d0))
         lmax = max(7,lmax)     ! Required for eps=1e-4
      endif
      lmax = min(maxl,lmax)
c
      if (.not. rtdb_get(rtdb, 'fastj:oversample', mt_dbl, 1, 
     $     oversample)) oversample = 1.8d0
c
c     Adjust the domain size so that we have one FMM box on
c     either side so that we have accurate potentials everywhere
c     the products of basis functions are non-zero
c
      side = max(xhi-xlo,yhi-ylo,zhi-zlo) * 0.5d0
      side = side * 2.0d0**levelbc / (2.0d0**levelbc - 2.0d0)
      x = (xhi+xlo)/2.0d0       ! Volume center
      y = (yhi+ylo)/2.0d0
      z = (zhi+zlo)/2.0d0
      xlo = x - side
      xhi = x + side
      ylo = y - side
      yhi = y + side
      zlo = z - side
      zhi = z + side
c
c     Determine the resolution or grid spacing which will be used to
c     determine if functions can be represented on the grid or not
c
      if (.not. rtdb_get(rtdb, 'fastj:resolution', mt_dbl, 1, 
     $     resolution)) resolution = 6.66667d0
c
      if (.not. rtdb_get(rtdb, 'fastj:ngrid', mt_int, 1, n))
     $     n = (xhi-xlo)*resolution*oversample
c
      n = min(n,539)
c
      call fastj_optimize_n(n, 'nearby')
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
      if (util_print('fastjinfo',print_debug)) 
     $     write(6,1) xlo, xhi, ylo, yhi, zlo, zhi
 1    format(' fastj: Solution volume: ',3(2f7.2,'  '))
c
c     Since we are temporaily restricted to cubic only do one dim
c
      resolution = 1d0/hz
      if (util_print('fastjinfo',print_debug)) 
     $     write(6,2) oversample, resolution/oversample, resolution, 
     $     nx, ny, nz,  depth, levelbc, lmax
 2    format(
     $       ' fastj:            Oversample ', f7.2/
     $       ' fastj:  FFT sampling density ', f7.2/
     $       ' fastj: FFT fine grid density ', f7.2/
     $       ' fastj:    No. of grid points ', 3i5/
     $       ' fastj:        FMM tree depth ', i5/
     $       ' fastj:  FMM tree level of bc ', i5/
     $       ' fastj:                  lmax ', i5)
c
c     Allocate the grid and tabulate the smooth component
c     of the density ... eventually this must be on a
c     nest of grids and we must compute the multipoles
c     and build the direct lists at the same time.
c
c     For efficient evaluation we allocate a nest of grids
c     so that we don't oversample smooth functions.
c
c     Since evaluating the density is so fast compared to the
c     FFT interpolation we don't evaluate on the very finest grid
c
      nnx = nx
      nny = ny
      nnz = nz
      do ngrid = 1, maxgrid
         if (ngrid .eq. 1) then
            status = ma_alloc_get(mt_dbl, nnx*nny*nnz, 'fastj: grid',
     $           l_grid, k_grid)
         else
            status = ma_push_get(mt_dbl, nnx*nny*nnz, 'fastj: grid',
     $           l_grid, k_grid) 
         end if
         if (.not. status) call errquit('fastj: ma grid',nnx*nny*nnz,
     &       UNKNOWN_ERR)
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
         if (util_print('fastjnest',print_debug)) 
     $        write(6,81) ngrid,
     $        (gridinfo(i,ngrid),i=1,3),(gridh(i,ngrid),i=1,3)
 81      format(' fastj: ',i2,' Allocated grid(',
     $        i4,',',i4,',',i4,')  with h(',f6.2,',',f6.2,',',f6.2,')')
c         
         if (ngrid .eq. 1) then
            nnx = (nnx+1)/oversample - 1
            nny = (nny+1)/oversample - 1
            nnz = (nnz+1)/oversample - 1
         else
            nnx = (nnx+1)/1.26 - 1
            nny = (nny+1)/1.26 - 1
            nnz = (nnz+1)/1.26 - 1
         endif
         call fastj_optimize_n(nnx, 'nearby')
         call fastj_optimize_n(nny, 'nearby')
         call fastj_optimize_n(nnz, 'nearby')
c
         if (nnx.lt.32 .or. nny.lt.32 .or. nnz.lt.32) goto 10
      end do
 10   if (ngrid .gt. maxgrid) ngrid = maxgrid
c
c     Also allocate space for the multipoles (entire octree) and
c     the linked list of shell-pairs assigned to each level.
c     The data linked by the list will be reallocated on the
c     heap when we know how big it is.
c
      lenpairs = nprimtot*(nprimtot+1)/2
*      lenpairs = 12*(nsh*(nsh+1)) ! 24*#pairs for primtives
      lenpoles = (8**(depth+1)-1)/7
      lenxlm   = (lmax+1)**2
      if (.not. ma_alloc_get(mt_int, lenpoles, 'fastj: list',
     $     l_list, k_list)) call errquit
     $     ('fastj: cannot allocate list ', lenpoles, MA_ERR)
      if (.not. ma_push_get(mt_int, lenpairs*2, 'fastj: pairs',
     $     l_pairs, k_pairs)) call errquit
     $     ('fastj: cannot allocate pairs for plan b ', lenpairs*2,
     &       MA_ERR)
c
      call fmm_tree_create(depth, lmax, tree)
c
      call fastj_density(ngrid, gridinfo, gridh,
     $     xlo, xhi, ylo, yhi, zlo, zhi,
     $     basis, geom, nat, nsh, nbf, g_dens, eps, oversample,
     $     q, lmax, depth, tree, int_mb(k_list), 
     $     int_mb(k_pairs), lenpairs, npairs)
c
      write(6,*) ' FMM tree after creating poles '
      call fmm_tree_summarize(tree)
      write(6,*)
c
c     Allocate minimal space for linked list on heap and copy
c
      if (.not. ma_alloc_get(mt_int, npairs*2, 'fastj: pairs2',
     $     l_pairs2, k_pairs2)) call errquit
     $     ('fastj: cannot allocate pairs ', npairs*2, MA_ERR)
      call icopy(npairs*2, int_mb(k_pairs), 1, int_mb(k_pairs2), 1)
c
c     Free up all of the memory we are no longer using
c
      if (ngrid .gt. 1) then
         if (.not. ma_chop_stack(gridinfo(5,2))) call errquit
     $        ('fastJ: ma chop failed ', 1, MA_ERR)
      end if
c
c     Only finest grid (no. 1) left now
c
      nx = gridinfo(1,1)
      ny = gridinfo(2,1)
      nz = gridinfo(3,1)
      k_grid = gridinfo(4,1)
      l_grid = gridinfo(5,1)
      hx = gridh(1,1)
      hy = gridh(2,1)
      hz = gridh(3,1)
c
c     Solve the sucker
c
      call solver(
     $     dbl_mb(k_grid), nx, ny, nz, 
     $     tree, depth, lmax, levelbc,
     $     xlo, xhi, ylo, yhi, zlo, zhi, .true., .true.)
c
      write(6,*) ' FMM tree after solving '
      call fmm_tree_summarize(tree)
      write(6,*)
c
      call fastj_pack_data(data,xlo,xhi,ylo,yhi,zlo,zhi,
     $     hx,hy,hz,k_grid,l_grid,tree,
     $     depth,lmax,nx,ny,nz,eps,l_list,k_list,l_pairs2,k_pairs2,
     $     npairs, g_dens, basis)
c
      total_time = util_cpusec() - total_time
      if (util_print('fastjinfo',print_debug)) write(6,1111) total_time
 1111 format(' fastj: Total solution time ', f8.2)
c
c     Now check the errors
c
      if (util_print('fastjerror',print_debug)) then
         errmax = 0.0d0
         errrms = 0.0d0
         errrelmax = 0.0d0
         errrelrms = 0.0d0
         nbig = 0
         do i = 1, ntest
*            x = xlo + 2.0d0*hx_fmm + util_random(0)*(xhi-xlo-4.0d0*hx_fmm)
*            y = ylo + 2.0d0*hy_fmm + util_random(0)*(yhi-ylo-4.0d0*hy_fmm)
*            z = zlo + 2.0d0*hz_fmm + util_random(0)*(zhi-zlo-4.0d0*hz_fmm)
 333        x = xlo + hx_fmm + util_random(0)*(xhi-xlo-2.0d0*hx_fmm)
            y = ylo + hy_fmm + util_random(0)*(yhi-ylo-2.0d0*hy_fmm)
            z = zlo + hz_fmm + util_random(0)*(zhi-zlo-2.0d0*hz_fmm)
            if ((x-xlo).lt.2*hx_fmm .or. (y-ylo).lt.2*hy_fmm .or.
     $           (z-zlo).lt.2*hz_fmm .or. (xhi-x).lt.2*hx_fmm .or.
     $           (yhi-y).lt.2*hy_fmm .or. (zhi-z).lt.2*hz_fmm) then
               rtest(1,i) = x
               rtest(2,i) = y
               rtest(3,i) = z
            else
               goto 333
            endif
         end do
c
         start = util_cpusec()
         call potential_list(basis, g_dens, ntest, rtest, ptest2, 
     $        eps*1d-2)
         used = (util_cpusec()-start)/ntest
         write(6,8321) used
 8321    format(' Time per point for potential_list  ', 1p, d9.2)

         start = util_cpusec()
         call fastj_eval_list_from_data(data, ntest, rtest, ptest)
         used = (util_cpusec()-start)/ntest
         write(6,8322) used
 8322    format(' Time per point for fastj_eval_list ', 1p, d9.2)
c
         do i = 1, ntest
            x = rtest(1,i)
            y = rtest(2,i)
            z = rtest(3,i)
*            test1 = fastj_eval_point(data, x, y, z)
            test2 = ptest2(i)   ! Old but faster
*            test3 = potential(basis, g_dens, x, y, z) ! Old and v. slow
            test3 = ptest(i)    ! New
            err = abs(test2 - test3)
            errrms = errrms + err*err
            if (abs(test2).gt.eps) then
               errrelmax = max(errrelmax, abs(err/test2))
               errrelrms = errrelrms + (err/test2)**2
            endif
            if (abs(err) .gt. errmax) then
               errmax = abs(err)
               xx = x
               yy = y
               zz = z
            end if
            if (abs(test2-test3) .gt. eps) then
               nbig = nbig + 1
*               write(6,781) i, x, y, z, test2, test3, test2-test3
* 781           format(1x,i5,2x,3f8.2,2x,2f16.10,1p,d9.2)
            end if
         end do
         errrms = sqrt(errrms / ntest)
         errrelrms = sqrt(errrelrms / ntest)
         write(6,82) errrms, errmax, errrelrms, errrelmax,
     $        xx, yy, zz, nbig
 82      format(' fastj: errrms', 1p, d9.2, 
     $        ', errmax=',d9.2,', errrelrms',d9.2, ', errrelmax',d9.2,
     $        0p/
     $         ' fastj: max error at ', 3f10.4, ' nbig', i5)
c     
         stop
c$$$         start = util_cpusec()
c$$$         do i = 1, 100
c$$$            x = xlo + 2.0d0*hx_fmm + util_random(0)*(xhi-xlo-4.0d0*hx_fmm)
c$$$            y = ylo + 2.0d0*hy_fmm + util_random(0)*(yhi-ylo-4.0d0*hy_fmm)
c$$$            z = zlo + 2.0d0*hz_fmm + util_random(0)*(zhi-zlo-4.0d0*hz_fmm)
c$$$            test1 = fastj_eval_point(data, x, y, z)
c$$$         end do
c$$$         used = (util_cpusec()-start)/100
c$$$         write(6,83) used
c$$$ 83      format(' fastj: Time per function evaluation ', 1p,d9.2,' sec')
      end if
c
      end
      subroutine fastj_density(ngrid, gridinfo, gridh,
     $     xlo, xhi, ylo, yhi, zlo, zhi,
     $     basis, geom, nat, nsh, nbf, g_dens, eps, oversample, q,
     $     lmax, depth, tree, list, pairs, lenpairs, npairs)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "bas.fh"
#include "geom.fh"
#include "util.fh"
      integer ngrid
      integer gridinfo(5,ngrid)
      double precision gridh(3,ngrid)
      double precision xlo, xhi, ylo, yhi, zlo, zhi
      integer basis, geom, nat, nsh, nbf
      integer g_dens
      double precision eps, oversample
      integer lmax, depth, lenpairs
      double precision q(-lmax:lmax,0:lmax)
      integer tree(300)
      integer list((8**(depth+1)-1)/7)
      integer pairs(2, lenpairs) ! Plan B
      integer npairs            ! [output]
c
c     Parameters
c
      integer maxprim           ! Max. no. primitives in a contraction
      parameter (maxprim=31)
      integer ijtmax, ijtmax2, lencart, lenprod, maxn
      parameter (ijtmax = 5)      ! Max. angular momentum of shell
      parameter (ijtmax2 = 2*ijtmax)
      parameter (lencart = ((ijtmax2+1)*(ijtmax2+2))/2)
      parameter (lenprod = ((ijtmax2+1)*(ijtmax2+2)*(ijtmax2+3))/6)
      parameter (maxn = 1024)   ! Max. grid dimension
c
      double precision d(lencart,-ijtmax2:ijtmax2,0:ijtmax2),
     $     dinv(lenprod, -ijtmax2:ijtmax2,0:ijtmax2)
      double precision work(lenprod**2), ndens(lenprod**2), 
     $     dens(lencart**2),qn(-ijtmax2:ijtmax2,0:ijtmax2)
      double precision xi, yi, zi, xj, yj, zj, ri(3), rj(3), s
      double precision expi(maxprim), coeffi(maxprim)
      double precision expj(maxprim), coeffj(maxprim)
      double precision factor, charge, test
      double precision a, b, c, range
      double precision alpha, prefac, tol, dmax, pi
      double precision h, hx, hy, hz, resolution, scale, start, used
c
      character*16 tag
c
      integer ish,iat,ishlo,ishhi,itype,inprim,ingen,numi,ibflo,ibfhi
      integer jsh,jat,jshlo,jshhi,jtype,jnprim,jngen,numj,jbflo,jbfhi
      integer i, iprim, numi_cart
      integer j, jprim, numj_cart
      integer k, nx, ny, nz, nnx, nny, nnz, jk, jjkk, ijk, ind, ptr
      integer ishtop, ijtype, igrid, k_grid, kk_grid, l, n
c
      integer level, numlevel(0:8), numgrid(1:16)
      double precision hfmm, dist, x, y, z, totalq, rootpi, radial,
     $     tolres, tolgrid
c
      double precision fuzzy
      parameter (fuzzy = .5)    ! Fraction of neighbour for leaks .. was 0.5
c
      logical iosph, josph

      double precision fastj_gaussian_range, double_factorial,
     $     fastj_r_neglected, fastj_resolution
      external fastj_gaussian_range, double_factorial,
     $     fastj_r_neglected, fastj_resolution
c
      integer packshprim, packpair, ishprim, jshprim
#include "bitops.fh"
      packshprim(ish,iprim) = ior(ishft(ish,5),iprim)
      packpair(ishprim,jshprim) = ior(ishft(ishprim,16),jshprim)
c     
      call xlm_coeff_inv(ijtmax2,d,dinv)
      call xlm_norm(ijtmax2,qn,ijtmax2)
      pi = 4.0d0*atan(1.0d0)
      totalq = 0.0d0
      call ifill(9,0,numlevel,1)
      call ifill(16,0,numgrid,1)
c
c     Initialize the linked list assigning shell pairs to the FMM tree
c
      call ifill((8**(depth+1)-1)/7, -1, list, 1)
      npairs = 0
c
      rootpi = sqrt(4.0d0*atan(1.0d0))
      tol = eps * 1d-3          ! Screen very tightly on radial prefactor
c     .                         ! 1d-2 not quite enuf for adz 5water
c
      if (nx.gt.maxn .or. ny.gt.maxn .or. nz.gt.maxn) call errquit
     $     ('fastj: hard dimension of maxn exceeded ', maxn,
     &       UNKNOWN_ERR)
c
      start = util_cpusec()
      do jat = 1, nat
         if (.not. bas_ce2cnr(basis, jat, jshlo, jshhi))
     $        call errquit('fastJ dens: bas_ce2cnr', jat, BASIS_ERR)
         if (.not. geom_cent_get(geom, jat, tag, rj, charge))
     $        call errquit('fastJ: geom_cart_get ', jat, GEOM_ERR)
         xj = rj(1)
         yj = rj(2)
         zj = rj(3)
         do iat = 1, jat
            if (.not. bas_ce2cnr(basis, iat, ishlo, ishhi))
     $           call errquit('fastJ dens: bas_ce2cnr', iat, BASIS_ERR)
            if (.not. geom_cent_get(geom, iat, tag, ri, charge))
     $           call errquit('fastJ: geom_cart_get ', iat, GEOM_ERR)
            xi = ri(1)
            yi = ri(2)
            zi = ri(3)
            do jsh = jshlo, jshhi
               if (.not.bas_continfo(basis, jsh, jtype, jnprim,jngen, 
     $              josph)) call errquit('fastJ: bas_cont_info', jsh,
     &       BASIS_ERR)
               if (jtype .lt. 0 .or. jngen.ne.1) call errquit
     $              ('fastJ: no general/SP functions yet',0, BASIS_ERR)
               if (jtype .gt. ijtmax) call errquit
     $              ('fastJ: hard dimension ijtmax exceeded ', ijtmax,
     &       BASIS_ERR)
               if (.not. bas_cn2bfr(basis, jsh, jbflo, jbfhi))
     $              call errquit('fastJ dens: bas_cn2bfr', jsh,
     &       BASIS_ERR)
               if (.not.bas_get_exponent(basis, jsh, expj))
     $              call errquit('fastJ: bas_get_exponent?',jsh,
     &       BASIS_ERR)
               if (.not.bas_get_coeff(basis, jsh, coeffj))
     $              call errquit('fastJ: bas_get_coeff?',jsh,
     &       BASIS_ERR)
               numj = jbfhi - jbflo + 1
c
               ishtop = ishhi
               if (iat .eq. jat) ishtop = jsh
               do ish = ishlo, ishtop
                  if (.not.bas_continfo(basis,ish,itype,inprim,ingen, 
     $                 iosph)) call errquit('fastJ: bas_cont_info',ish,
     &       BASIS_ERR)
                  if (itype .gt. ijtmax) call errquit
     $                 ('fastJ: hard dimension ijtmax exceeded ',ijtmax,
     &       BASIS_ERR)
                  if (.not. bas_cn2bfr(basis, ish, ibflo, ibfhi))
     $                 call errquit('fastJ dens: bas_cn2bfr', ish,
     &       BASIS_ERR)
                  if (.not.bas_get_exponent(basis, ish, expi))
     $                 call errquit('fastJ: bas_get_exponent?',ish,
     &       BASIS_ERR)
                  if (.not.bas_get_coeff(basis, ish, coeffi))
     $                 call errquit('fastJ: bas_get_coeff?',ish,
     &       BASIS_ERR)
                  numi = ibfhi - ibflo + 1
c
                  ijtype = itype + jtype
c
c     Copy the block of density into local workspace and transform
c     to cartesians if necessary
c
                  factor = 2.0d0 ! This assumes symmetric density!
                  if (ish .eq. jsh) factor = 1.0d0
                  factor = factor
c
                  numi_cart = (itype+1)*(itype+2)/2
                  numj_cart = (jtype+1)*(jtype+2)/2
                  if (iosph) then
                     call ga_get(g_dens,ibflo,ibfhi,jbflo,jbfhi,
     $                    ndens,numi)
                     call dens_sph_to_cart(itype, jtype, 
     $                    ndens, numi, dens, work, numi_cart, numj_cart)
                     dmax = 0.0d0
                     do i = 1, numi_cart*numj_cart
                        dmax = max(dmax,abs(dens(i)))
                     end do
                  else
                     call ga_get(g_dens,ibflo,ibfhi,jbflo,jbfhi,
     $                    dens,numi)
                     dmax = 0.0d0
                     do i = 1, numi*numj
                        dmax = max(dmax,abs(dens(i)))
                     end do
                  end if
c
                  do jprim = 1, jnprim
                     do iprim = 1, inprim
                        call gaussian_product(
     $                       expi(iprim), xi, yi, zi,
     $                       expj(jprim), xj, yj, zj,
     $                       alpha, a, b, c, prefac)
c
                        prefac = prefac*factor*
     $                       coeffi(iprim)*coeffj(jprim) 
c
c     Test estimates the magnitude of any component of the product function
c
                        test = 0.0d0
                        scale = abs(prefac*dmax*rootpi/
     $                       (4.0d0*alpha*sqrt(alpha)))
                        do n = 0, ijtype
                           radial = scale
                           do l = n, 0, -2
                              test = max(test,
     $                             radial*double_factorial(n+l+1)/
     $                             (qn(l,l)*qn(l,l)))
                              radial = radial * (2.0d0*alpha)
                           end do
                           scale = scale / (2.0d0*alpha)
                        end do
c
                        if (test.gt.tol) then
*                           range = fastj_r_neglected(
*     $                          ijtype,alpha, 
*     $                          min(tol/abs(dmax*prefac),1d-4))

                           range = fastj_r_neglected(
     $                          ijtype,alpha, 
     $                          min(1d-2,eps/abs(dmax*prefac)))


*                           write(6,*) ish, jsh, iprim, jprim,
*     $                          ' RANGE ', range, dmax, prefac
c
c     Pick the most appropriate grid, reserving the finest grid for
c     oversampling.
c
                           tolres = min(1d-4,0.1d0*eps/abs(dmax*prefac))
                           resolution = fastj_r_neglected(ijtype,
     $                          0.25d0/alpha, tolres)/pi

*     Seems to work but no great advantage on systems tested
*                           tolres = min(1d-4,0.1d0*eps/abs(dmax*prefac))
*                           resolution = 1.1d0*
*     $                          fastj_resolution(alpha,ijtype,tolres)

                           h = 1.0d0/resolution
                           do igrid = ngrid,2,-1 ! Don't eval on grid 1
                              if (h .ge. gridh(1,igrid)) goto 333
                           end do
 333                       if (igrid .le. 1) then ! Don't eval on 1
c
c     The function cannot be represented on any available grid.
c     Find the smallest FMM box that will hold at least a large
c     fraction ... the distance of leakage into neighbouring
c     cells must be less than fuzzy*width
c     
                              ind = 1
c
                              do level = depth, 0, -1
                                 hfmm = (xhi-xlo)/2.0d0**level
                                 i = int((a-xlo)/hfmm)
                                 j = int((b-ylo)/hfmm)
                                 k = int((c-zlo)/hfmm)
                                 x = xlo + (dble(i)+0.5d0)*hfmm ! Box center
                                 y = ylo + (dble(j)+0.5d0)*hfmm
                                 z = zlo + (dble(k)+0.5d0)*hfmm
c     
c     To compute the distance to the nearest edge, compute 
c     max(dist. from center along each axis) and then is just h/2-dist.
c
                                 dist = max(abs(a-x),abs(b-y),abs(c-z))
                                 dist = hfmm*0.5d0 - dist ! Distance to edge
                                 dist = max(0d0,range-dist) ! Distance of spill
                                 if (dist .lt. fuzzy*hfmm) goto 222
                                 ind = ind + 8**level
                              end do
 222                          if (level .lt. 0) call errquit
     $                             ('fastj: solution volume?',0,
     &       UNKNOWN_ERR)
* 31                           format(1x,a,4i5,a,4i5,2f6.2)
*                              write(6,31) ' Assigning ', ish, iprim,
*     $                             jsh, jprim, ' to fmm box ', 
*     $                             level, i, j, k, hfmm, range
                              numlevel(level) = numlevel(level) + 1
c
c     Compute the multipole moments about the GPT center
c     then translate to the box center (expensive!)
c     
                              call pot_prim_shell_pair_multipoles(
     $                             dens,
     $                             itype, ri, numi_cart,
     $                             jtype, rj, numj_cart,
     $                             dinv, ijtmax2,
     $                             a, b, c, alpha, prefac,
     $                             q, lmax)
                              totalq = totalq + q(0,0)
*                              write(6,*) ' Q before translate '
*                              call xlm_print(lmax,q,lmax)
                              call xlm_translate2(lmax, 
     $                             x-a, y-b, z-c, q, lmax)
*                              write(6,*) ' Q after translate '
*                              call xlm_print(lmax,q,lmax)
*                              call xlm_translate2(lmax, 
*     $                             a-x, b-y, c-z, q, lmax)
*                              write(6,*) ' Q after translate back '
*                              call xlm_print(lmax,q,lmax)
*                              call xlm_translate2(lmax, 
*     $                             x-a, y-b, z-c, q, lmax)

c     
                              ind = ind+i+(2**level)*(j+2**level*k)
c
         if (ind .le. 0 .or. ind.gt.(8**(depth+1)-1)/7) call errquit
     $        ('ind ? ', ind, UNKNOWN_ERR)
*                              call xlm_accumulate_to_packed(lmax,
*     $                             poles(1,ind), q, lmax)
         call fmm_tree_add_cell(tree, level, i, j, k, q, lmax, lmax)
c     
c     Add the shell pair to the linked list of that box.  Rather
c     than add the primitive shell pair, add the contracted shell pair
c     if not already there.  The near field code must reproduce
c     the logic immediately above to figure out if a given primitive
c     pair is in that box.
c
c     Well, the above plan seems too prone to error ... instead we
c     save the shell and primitive data packed to save memory
c     
c$$$     This implements plan A ... note memory alloc must be changed
c$$$                              ptr = list(ind)
c$$$                              if (ptr .gt. 0) then
c$$$                                 odoit = ( pairs(1,ptr).ne.ish .or.
c$$$     $                                pairs(2,ptr).ne.jsh) 
c$$$                              else
c$$$                                 odoit = .true.
c$$$                              endif
c$$$                              if (odoit) then
c$$$                                 npairs = npairs + 1
c$$$                                 if (npairs.gt.lenpairs)call errquit
c$$$     $                                ('fastj:too many pairs',npairs)
c$$$                                 pairs(1,npairs) = ish
c$$$                                 pairs(2,npairs) = jsh
c$$$                                 pairs(3,npairs) = list(ind)
c$$$                                 list(ind) = npairs
c$$$                                 write(6,*) ' Adding ', ish, jsh,
c$$$     $                                level, i, j, k, ind, npairs
c$$$                              end if
c
c     Plan B
c
                              ptr = list(ind)
                              npairs = npairs + 1
                              if (npairs.gt.lenpairs)call errquit
     $                             ('fastj:too many pairs',npairs,
     &       UNKNOWN_ERR)
                              pairs(1,npairs) = packpair(
     $                             packshprim(ish,iprim),
     $                             packshprim(jsh,jprim))

c$$$
c$$$                              iijj = pairs(1,npairs)
c$$$                              iish  = iand(rshift(iijj,21),2047)
c$$$                              iiprim= iand(rshift(iijj,16),31)
c$$$                              ijsh  = iand(rshift(iijj,5),2047)
c$$$                              ijprim= iand(iijj,31)



                              pairs(2,npairs) = list(ind)
                              list(ind) = npairs
*                              write(6,*) ' Adding ', ish, iprim,
*     $                             jsh, jprim,
*     $                             level, i, j, k, ind, npairs
c$$$                              write(6,*) '*Adding ', iish, iiprim,
c$$$     $                             ijsh, ijprim
c
                           else
c     
c     The function is representable on a grid ... tabulate it
c     
                              numgrid(igrid) = numgrid(igrid) + 1
*                              write(6,*) ' Assigning ', ish, iprim,
*     $                             jsh, jprim, ' to grid ', igrid
                              call cart_dens_translate(
     $                             itype, xi, yi, zi, 
     $                             jtype, xj, yj, zj,
     $                             a, b, c, dens, work, ndens)
                              call cart_dens_product(
     $                             itype, jtype, work, ndens)
                              nx = gridinfo(1,igrid)
                              ny = gridinfo(2,igrid)
                              nz = gridinfo(3,igrid)
                              k_grid = gridinfo(4,igrid)
                              hx = gridh(1,igrid)
                              hy = gridh(2,igrid)
                              hz = gridh(3,igrid)
c     
                              range = fastj_r_neglected(
     $                             ijtype,alpha, 
     $                             min(1d-4,0.1d0*eps/abs(dmax*prefac)))
*                              tolgrid = 1d-2*eps/
*     $                             (4d0*3.14d0*range*range*hx)
                              tolgrid = 1d-4*eps
                              call fastj_density_eval(
     $                             ndens,
     $                             ijtype, alpha, prefac, range, 
     $                             a, b, c, 
     $                             dbl_mb(k_grid), nx, ny, nz, 
     $                             hx, hy, hz,
     $                             xlo, xhi, ylo, yhi, zlo, zhi, 
     $                             tolgrid)

                           end if !if (igrid)
                        end if   ! if (screen)
                     end do      ! do iprim
                  end do         ! do jprim
               end do            ! do ish
            end do               ! do jsh
         end do                  ! do iat
      end do                     ! do jat
c
      if (util_print('fastjinfo', print_debug)) then
         write(6,301) 'fastj: Numgrid  ', (numgrid(i),i=1,ngrid)
 301     format(1x,a,20i6)
         write(6,301) 'fastj: Numlevel ', (numlevel(i),i=0,depth)
      endif
c
      used = util_cpusec() - start
      if (util_print('fastjtime',print_debug))
     $     write(6,12) used
 12   format(' fastj: Evaluation of the density took ', f6.2)
c
c     Interpolate all of the grids down to the finest using the FFT
c
c     Note that we did not evaluate on grid 1 to save the cost of the FFT
c
      start = util_cpusec()
      do igrid = ngrid,2,-1     ! No need to do grid 1
         nx = gridinfo(1,igrid)
         ny = gridinfo(2,igrid)
         nz = gridinfo(3,igrid)
         k_grid = gridinfo(4,igrid)
         call vsint3d(nx,ny,nz,dbl_mb(k_grid))
         scale = 0.125d0/(dble(nx+1)*dble(ny+1)*dble(nz+1))
         call dscal(nx*ny*nz,scale,dbl_mb(k_grid),1)
      end do
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
               end do
            end do
         end do
      end do
c
c     Here can either transform back and let the solver do a redundant
c     forward transformation, OR, stay in Fourier space.  However, if
c     stay in Fourier space need to rescale after the interpolation.
c     Note that we cannot mess with the density without doing the
c     FFT.
c
**      call vsint3d(nnx,nny,nnz,dbl_mb(kk_grid)) ! can be optimized away
c
      scale = 8.0d0*dble(nnx+1)*dble(nny+1)*dble(nnz+1)
      call dscal(nnx*nny*nnz,scale,dbl_mb(kk_grid),1)
c      
      used = util_cpusec() - start
      if (util_print('fastjtime',print_debug))
     $     write(6,13) used
 13   format(' fastj: Fourier interpolation to fine grid ',f6.2)
c
      if (util_print('fastjinfo',print_debug)) then
         s = 0.0d0
         do k = 1, nnz, 2
            do j = 1, nny, 2
               do i = 1, nnx, 2
                  ijk = kk_grid + i-1 + nnx*(j-1 + nny*(k-1))
                  s = s + dbl_mb(ijk)/(dble(i)*dble(j)*dble(k))
               end do
            end do
         end do
         s = s*8.0d0*gridh(1,1)**3 / pi**3
         write(6,33) s, totalq, s+totalq
 33      format(
     $        ' fastj: Norm of density on grid = ',f20.12/
     $        ' fastj: Norm of density in fmm  = ',f20.12/
     $        ' fastj: Norm of density total   = ',f20.12)
      end if
c
c     Screen the density for -ve values and check norm ... only
c     if we did the FFT!
c
c$$$      s = 0.0d0
c$$$      gmin = 0.0d0
c$$$      do k = 1, nnz
c$$$         do j = 1, nny
c$$$            do i = 1, nnx
c$$$               gmin = min(gmin,dbl_mb(kk_grid))
c$$$               if (dbl_mb(kk_grid) .lt. 0.0d0) dbl_mb(kk_grid) = 0.0d0
c$$$               s = s + dbl_mb(kk_grid)
c$$$               kk_grid = kk_grid + 1
c$$$            end do
c$$$         end do
c$$$      end do
c$$$      if (util_print('fastjdnorm',print_debug)) then
c$$$         s = s * gridh(1,1)**3.0d0
c$$$         write(6,31)  gmin
c$$$ 31      format(' fastj: Largest negative element ',1p,d9.2)
c$$$         write(6,32) s, s-nint(s)
c$$$ 32      format(' fastj: Trapezoidal integration of the density ',
c$$$     $        f15.10,2x,1p,d9.2)
c$$$      end if
c     
      end
      subroutine fastj_density_eval(
     $     ndens,
     $     ijtype, alpha, prefac, range, a, b, c, 
     $     g, nx, ny, nz, hx, hy, hz,
     $     xlo, xhi, ylo, yhi, zlo, zhi, tol)
      implicit none
      double precision ndens(*)
      integer ijtype
      double precision alpha, prefac, range, a, b, c, tol
      integer nx, ny, nz
      double precision g(nx, ny, nz), hx, hy, hz
      double precision xlo, ylo, zlo, xhi, yhi, zhi
c     
      integer maxn, ijtmax2
      parameter (maxn=1024, ijtmax2=10)
      double precision xx(maxn,0:ijtmax2), yy(maxn,0:ijtmax2),
     $     zz(maxn,0:ijtmax2)
      double precision ndmax
      double precision x, y, z, xxmax, s, sk
      integer i, j, k, l, p, ii, jj, kk, ind
      integer ilo, ihi, jlo, jhi, klo, khi
*      double precision fastj_gaussian_range
*      external fastj_gaussian_range
c
      ndmax = 0.0d0
      do i = 1,((ijtype+1)*(ijtype+2)*(ijtype+3))/6
         ndmax = max(ndmax,abs(ndens(i)))
      end do
c     
*      range = 1.1d0*fastj_gaussian_range(ijtype, alpha, 
*     $     min(tol/abs(ndmax*prefac),1d-4))
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
      subroutine fastj_tidy(data)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
      integer data(2)
c
      double precision xlo,xhi,ylo,yhi,zlo,zhi,eps,hx,hy,hz
      integer k_grid, depth,lmax,nx,ny,nz, l_grid, tree(300)
      integer l_list,k_list,l_pairs2,k_pairs2,
     $     npairs,g_dens,basis
c
      call fastj_unpack_data(data,xlo,xhi,ylo,yhi,zlo,zhi,
     $     hx,hy,hz,k_grid,l_grid,tree,
     $     depth,lmax,nx,ny,nz,eps,l_list,k_list,l_pairs2,k_pairs2,
     $     npairs,g_dens,basis)
c
      if (.not. ma_free_heap(l_grid)) call errquit
     $     ('fastj_tidy: ma?',-1, MA_ERR)
      if (.not. ma_free_heap(l_list)) call errquit
     $     ('fastj_tidy: ma?',-1, MA_ERR)
      if (.not. ma_free_heap(l_pairs2)) call errquit
     $     ('fastj_tidy: ma?',-1, MA_ERR)
      if (.not. ma_free_heap(data(1))) call errquit
     $     ('fastj_tidy: ma?',1, MA_ERR)
      if (.not. ma_free_heap(data(2))) call errquit
     $     ('fastj_tidy: ma?',2, MA_ERR)
c
      call fmm_tree_destroy(tree)
c
      end
      subroutine fastj_optimize_n(n, mode)
      implicit none
      integer n
      character*(*) mode
c
c     Find a number of the form 2^i 3^j 5^k 7^l - 1 that is close to n.
c     Optionally force it to be an upperbound.
c
c     For ESSL need to impose additional constraints.
c
      integer n2, n3, n5, n7, ii2, ii3, ii5, ii7, i2, i3, i5, i7
      double precision test, diff
c
      n = n + 1                 ! Temporarily to force to power of radix
c
      n2 = log(dble(n))/log(2.0d0) + 1
      n3 = log(dble(n))/log(3.0d0) + 1
      n5 = log(dble(n))/log(5.0d0) + 1
      n7 = log(dble(n))/log(7.0d0) + 1
#ifdef ESSL
      n3 = min(2,n3)
      n5 = min(1,n5)
      n7 = min(1,n7)
#else
      n7 = 0
#endif
      ii2 = 0
      ii3 = 0
      ii5 = 0
      ii7 = 0
      diff = n
      do i2 = 0, n2
         do i3 = 0, n3
            do i5 = 0, n5
               do i7 = 0, n7
                  test = (2d0**i2)*(3d0**i3)*(5d0**i5)*(7d0**i7) + 1d-6
                  if (test.ge.n .or. mode.ne.'upperbound') then
                     if (abs(test-n) .lt. diff) then
                        diff = abs(test - n)
                        ii2 = i2
                        ii3 = i3
                        ii5 = i5
                        ii7 = i7
                     end if
                  end if
               end do
            end do
         end do
      end do
c
      n = (2d0**ii2)*(3d0**ii3)*(5d0**ii5)*(7d0**ii7)  - 1
c
      end
c$$$      double precision function fastj_eval_point(data,x,y,z)
c$$$      implicit none
c$$$#include "mafdecls.fh"
c$$$c
c$$$      double precision x, y, z
c$$$      integer data(2)
c$$$c
c$$$c     Evaluate the potential at a point ... not as efficient
c$$$c     as doing multiple points
c$$$c
c$$$      integer k_grid, lmax, depth, nx, ny, nz, l_grid
c$$$      double precision xlo, xhi, ylo, yhi, zlo, zhi, hx, hy, hz, eps
c$$$      double precision potgrid, potfmm, hfmm
c$$$      integer order, ninterp
c$$$      double precision tn_interp_3d_point, fmm_potential
c$$$      integer tree(300)
c$$$c
c$$$c     Get the data
c$$$c
c$$$      call fastj_unpack_data(
c$$$     $     data,
c$$$     $     xlo,xhi,ylo,yhi,zlo,zhi,
c$$$     $     hx, hy, hz,
c$$$     $     k_grid, l_grid, tree,
c$$$     $     depth,lmax,nx,ny,nz,eps)
c$$$c
c$$$      hfmm = (xhi - xlo)/2.0d0**depth
c$$$      if ( x.lt.(xlo+hfmm) .or. x.gt.(xhi-hfmm) .or.
c$$$     $     y.lt.(ylo+hfmm) .or. y.gt.(yhi-hfmm) .or.
c$$$     $     z.lt.(zlo+hfmm) .or. z.gt.(zhi-hfmm) ) then
c$$$         fastj_eval_point = 0.0d0
c$$$         return
c$$$      end if
c$$$c
c$$$*      write(6,*) ' xlo,xhi,ylo,yhi,zlo,zhi ',xlo,xhi,ylo,yhi,zlo,zhi
c$$$*      write(6,*) ' hx, hy, hz ', hx, hy, hz
c$$$*      write(6,*) ' k_grid, k_poles,depth,lmax,nx,ny,nz,eps ',
c$$$*     $     k_grid, k_poles,depth,lmax,nx,ny,nz,eps
c$$$c
c$$$      if (eps.ge.1d-6) then
c$$$         order = 17
c$$$      else if (eps.ge.1d-8) then
c$$$         order = 25
c$$$      else
c$$$         order = 33
c$$$      end if
c$$$      ninterp = order + 1 + 2
c$$$c
c$$$      potgrid = tn_interp_3d_point(dbl_mb(k_grid),nx,ny,nz,hx,hy,hz,
c$$$     $     x-xlo,y-ylo, z-zlo,ninterp,order)
c$$$c
c$$$      potfmm = fmm_potential(depth,lmax,tree,
c$$$     $     xhi-xlo,yhi-ylo,zhi-zlo,x-xlo,y-ylo,z-zlo)
c$$$c
c$$$**      write(6,*) ' potgrid, potfmm ', potgrid, potfmm
c$$$c
c$$$      fastj_eval_point = -(potgrid + potfmm)
c$$$c      
c$$$      end
      subroutine fastj_eval_list_from_data(data, n, coords, v)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "util.fh"
c
      integer n
      double precision coords(3,n), v(n)
      integer data(2)
c
c     Evaluate the potential at a list of points within the volume
c
      integer k_grid, lmax, depth, nx, ny, nz, l_grid
      double precision xlo, xhi, ylo, yhi, zlo, zhi, hx, hy, hz, eps
      integer order, ninterp, twol
      integer tree(300)
      integer l_list, k_list, l_pairs2, k_pairs2, npairs, g_dens, basis
      double precision tn_interp_3d_point
*      double precision fmm_potential
c
      integer k_map, l_map, k_list_info, l_list_info, k_next, l_next,
     $     k_coords_r, kx, ky, kz, ptr, numr, k_v_r, l_v_r,
     $     ish, jsh, l_coords_r, ind, level, k_near, l_near,
     $     k_v_near, l_v_near, jj, k_xyz, iijj, iprim, jprim
      integer kxp, kyp, kzp, luse, levelp
      integer lq
      parameter (lq = FMM_LMAX)
      double precision q(-lq:lq,0:lq)
      double precision hxp, hyp, hzp, xo, yo, zo, rtwol
      double precision x, y, z, potfmm, potgrid, fact
      logical present(31*31)    ! Max. prim. sq
c
      integer sh_info(4,2047)   ! 1=type; 2=nprim; 3=ngen; 4=osph
      integer sh_bfr(2,2047)    ! 1=lo; 2=hi
      integer ishprev, jshprev, nprimj, nprimi, numi, numj, i, j
      
      integer level_do_grid
c
      integer l_fit, k_fit
c
      integer lmax2, lencart, lenprod
      parameter (lmax2 = 2*5)   ! Twice max angular momentum - h functions
      parameter (lencart = ((lmax2+1)*(lmax2+2))/2)
      parameter (lenprod = ((lmax2+1)*(lmax2+2)*(lmax2+3))/6)
c
      double precision start, used, tfar, tgrid
      integer ilo, ihi, jlo, jhi, klo, khi
      double precision vol_xlo, vol_ylo, vol_zlo ! Volume for interp
      integer k_gtmp, l_gtmp, ngrid_per_box, ngx, ngy, ngz,
     $     ngrid_per_vol, nborder, maxorder, minorder,
     $     max_ngrid_per_vol, stride, ind_save
      double precision hfmm, potgrid1, xlm_local_potential
c
      double precision 
     $     d(lencart,-lmax2:lmax2,0:lmax2),
     $     dinv(lenprod, -lmax2:lmax2,0:lmax2),
     $     dens(lencart**2)

      double precision obin(0:34), tbin(0:34)

      integer numfit
#include "bitops.fh"
      numfit(i) = (2*i+9)*(i+2)*(i+1)/6 ! Sum of (n+1)*(n+3) n=0..i

      call xlm_init()
      call anl_fit_init()
      call xlm_coeff_inv(lmax2,d,dinv)
      call dfill(35,0.0d0,obin,1)
      call dfill(35,0.0d0,tbin,1)
c
      ngx = 1         ! get rid of compiler warnings
      ngy = 1
      ngz = 1
c
c     Get the data
c
      call fastj_unpack_data(
     $     data,
     $     xlo,xhi,ylo,yhi,zlo,zhi,
     $     hx, hy, hz,
     $     k_grid, l_grid, tree,
     $     depth,lmax,nx,ny,nz,eps,l_list,k_list,l_pairs2,k_pairs2,
     $     npairs, g_dens, basis)
c
      write(6,*) ' Unpacking'
      write(6,*) ' xlo,xhi,ylo,yhi,zlo,zhi ',xlo,xhi,ylo,yhi,zlo,zhi
      write(6,*) ' hx, hy, hz ', hx, hy, hz
      write(6,*) ' k_grid, depth,lmax,nx,ny,nz,eps ',
     $     k_grid, depth,lmax,nx,ny,nz,eps
      write(6,*) ' l_list, k_list, l_pairs2, k_pairs2, npairs ',
     $     l_list, k_list, l_pairs2, k_pairs2, npairs
      write(6,*) ' g_dens ', g_dens
      call fmm_tree_summarize(tree)
c

      level_do_grid = min(4,depth)
      if (nx .eq. 0) level_do_grid = 999

      if (eps .ge. 1d-4) then
         minorder = 3           ! Actually will use 5.
         maxorder = 7
         nborder = 3
         stride = 2
      else if (eps .ge. 1d-6) then
         minorder = 3           ! Actually will use 5.
         maxorder = 15
         nborder = 3
         stride = 2
      else if (eps .ge. 1d-8) then
         minorder = 3           ! Actually will use 7
         maxorder = 23          ! was 23
         nborder = 5
         stride = 4
      else
         minorder = 11           ! Actually will use 11
         maxorder = 27          
         nborder = 5             ! ?
         stride = 4
      endif
c
      call dfill(n, 0.0d0, v, 1)
c
      call bas_all_cn_info(basis, sh_info, sh_bfr)
c
      if (.not. ma_push_get(mt_int, n, 'map', l_map, k_map))
     $     call errquit('fj: map ', n, MA_ERR)
      if (.not. ma_push_get(mt_int,2*8**depth,'info',l
     $     _list_info,k_list_info)) call errquit('fj: info ',2*8**depth,
     &       GA_ERR)
      if (.not. ma_push_get(mt_int, n, 'next', l_next, k_next))
     $     call errquit('fj: next ', n, MA_ERR)
      if (.not. ma_push_get(mt_dbl, 3*n, 'near', l_near, k_near))
     $     call errquit('fj: near ', 3*n, MA_ERR)
      if (.not. ma_push_get(mt_dbl, n, 'v_near', l_v_near, k_v_near))
     $     call errquit('fj: v_near', n, MA_ERR)
      if (.not. ma_push_get(mt_dbl, 3*n, 'cr', l_coords_r, k_coords_r))
     $     call errquit('fj: coords_r ', 3*n, MA_ERR)
      if (.not. ma_push_get(mt_dbl, n, 'vr', l_v_r, k_v_r))
     $     call errquit('fj: vr ', n, MA_ERR)
      if (.not. ma_push_get(mt_dbl, numfit(maxorder), 
     $     'fit', l_fit, k_fit))
     $     call errquit('fj: fit ',numfit(maxorder), MA_ERR)
c
      do i = minorder, maxorder
         if (i .le. 9) then
            call tn_fitting_matrix(i,i+1,dbl_mb(k_fit+numfit(i-1)))
         else
            call tn_fitting_matrix(i,i+3,dbl_mb(k_fit+numfit(i-1)))
         endif
      enddo
c
      hfmm = (xhi-xlo)/2.0d0**level_do_grid ! Assumes cubic box
      ngrid_per_box = hfmm/hx  + 1
      ninterp = maxorder + 1 + 2 ! Dynamically reset below
      if (ngrid_per_box .lt. ninterp) 
     $     nborder = nborder + (ninterp-ngrid_per_box-1)/2 + 1
      max_ngrid_per_vol = ngrid_per_box + 2*nborder
      if (mod(max_ngrid_per_vol,2).eq.0) 
     $     max_ngrid_per_vol = max_ngrid_per_vol + 1
      write(6,*) ' max_ngrid_per_vol ', max_ngrid_per_vol
c      
      if (.not. ma_push_get(mt_dbl, (max_ngrid_per_vol)**3, 
     $     'gtmp', l_gtmp, k_gtmp)) call errquit('fj: gtmp ', 
     $     max_ngrid_per_vol**3, MA_ERR)
c
c     Pairs may have been assigned to the top few levels of the
c     octree if they span all/most of the domain.  Also, ind
c     is being accumulated so must start from the finest
c
      ind = 0
      do level = depth, 0, -1
         ind_save = ind
         start = util_cpusec()

         twol = 2**level

         call fmm_points_to_boxes(level, n, coords,
     $        xlo, xhi, ylo, yhi, zlo, zhi,
     $        int_mb(k_next), int_mb(k_list_info), 
     $        dbl_mb(k_coords_r), int_mb(k_map))
         call dfill(n, 0.0d0, dbl_mb(k_v_r), 1)

         do kz = 0, twol-1
            do ky = 0, twol-1
               do kx = 0, twol-1
c
                  ptr = int_mb(k_list+ind) ! Head of shell pair linked list
c
                  if (ptr .gt. 0) then 
                     call fmm_make_list_of_near_points(
     $                    level, kx, ky, kz,
     $                    n, dbl_mb(k_coords_r), int_mb(k_list_info),
     $                    dbl_mb(k_near), numr)
                  else
                     numr = 0
                  endif
c
                  if (numr .gt. 0) then
*                     write(6,*) level, kx, ky, kz, ind+1, numr
                     call dfill(numr, 0.0d0, dbl_mb(k_v_near), 1)
c
c     Traverse the linked list of primtive shell pairs for this box/level
c     and evaluate the near field potentials.
c     Also bundle back into contracted shell as an optimization.
c     Note that some pairs of primtives may be missing so must
c     form a mask to see if pairs are present or not.
c
c
                     ishprev = -1
                     jshprev = -1
 10                  if (ptr .gt. 0) then
                        ptr = (ptr - 1)*2
                        iijj = int_mb(k_pairs2+ptr  )
                        ptr  = int_mb(k_pairs2+ptr+1)
c
                        ish  = iand(ishft(iijj,21),-2047)
                        iprim= iand(ishft(iijj,16),-31)
                        jsh  = iand(ishft(iijj,5),-2047)
                        jprim= iand(iijj,31)
                        if (ish.ne.ishprev .or. jsh.ne.jshprev) then
                           if (ishprev .ne. -1) then
c     EVALUATE HERE using ishprev and jshprev
*                              write(6,*) '   Evaluating1 ',
*     $                             ishprev,jshprev, nprimi, nprimj

                              call ga_get(g_dens, 
     $                             sh_bfr(1,ishprev), sh_bfr(2,ishprev),
     $                             sh_bfr(1,jshprev), sh_bfr(2,jshprev),
     $                             dens, numi)

                              fact = 1.0d0
                              if (ishprev.ne.jshprev) fact = 2.0d0

                              call potential_from_shell_pair(
     $                             basis, ishprev, jshprev,
     $                             dens, numi, eps*1d-1, d, dinv, lmax2,
     $                             numr, 
     $                             dbl_mb(k_near), dbl_mb(k_v_near),
     $                             fact, .true., present, nprimi)
     $                             
                           endif
                           ishprev = ish
                           jshprev = jsh
                           nprimi = sh_info(2,ish)
                           nprimj = sh_info(2,jsh)
                           numi = sh_bfr(2,ish) - sh_bfr(1,ish) + 1
                           numj = sh_bfr(2,jsh) - sh_bfr(1,jsh) + 1
                           do i = 1, nprimi*nprimj
                              present(i) = .false.
                           enddo
                        endif
c
*                        write(6,*)  '     ', ish, iprim, jsh, jprim
                        present(iprim + (jprim-1)*nprimi) = .true.
                        goto 10
                     endif
                     if (ishprev .ne. -1) then
c     EVALUATE HERE using ishprev and jshprev
*                        write(6,*) '   Evaluating2 ',
*     $                       ishprev,jshprev, nprimi, nprimj
                        call ga_get(g_dens, 
     $                       sh_bfr(1,ishprev), sh_bfr(2,ishprev),
     $                       sh_bfr(1,jshprev), sh_bfr(2,jshprev),
     $                       dens, numi)
                        
                        fact = 1.0d0
                        if (ishprev.ne.jshprev) fact = 2.0d0
                        
                        call potential_from_shell_pair(
     $                       basis, ishprev, jshprev,
     $                       dens, numi, eps*1d-1, d, dinv, lmax2,
     $                       numr, 
     $                       dbl_mb(k_near), dbl_mb(k_v_near),
     $                       fact, .true., present, nprimi)
     $                       
                     endif
                     call fmm_add_values_from_near_points(
     $                    level, kx, ky, kz,
     $                    n, dbl_mb(k_v_r), int_mb(k_list_info),
     $                    dbl_mb(k_v_near), numr)
c     add the potential from the near field into reordered points
                  endif
                  ind = ind + 1
               enddo
            enddo
         enddo
c
c     Do the grid
c
         tgrid = 0.0d0
         if (level.eq.level_do_grid .and. nx.gt.0) then
            tgrid = util_cpusec()
            ind = ind_save
            do kz = 0, twol-1
               do ky = 0, twol-1
                  do kx = 0, twol-1
                     k_xyz = k_list_info + 2*(kx + twol*(ky + twol*kz))
c
c     Gather a cube of points with nborder points on either side
c     in order to improve memory locality.  
c
c     Tried examining necessary precision at box center but it did
c     not work ... now examine each point.
c
                     if (int_mb(k_xyz) .le. int_mb(k_xyz+1)) then
c
                        ngrid_per_vol = max_ngrid_per_vol
c     
                        x = xlo + dble(kx)*hfmm ! Coords of box corner
                        y = ylo + dble(ky)*hfmm
                        z = zlo + dble(kz)*hfmm
                        ilo = max(1, int((x-xlo)/hx)-nborder)
                        ilo = min(ilo,nx-ngrid_per_vol+1)
                        ilo = max(1,ilo)
                        jlo = max(1, int((y-ylo)/hy)-nborder)
                        jlo = min(jlo,ny-ngrid_per_vol+1)
                        jlo = max(1,jlo)
                        klo = max(1, int((z-zlo)/hz)-nborder)
                        klo = min(klo,nz-ngrid_per_vol+1)
                        klo = max(1,klo)
                        ihi = min(nx, ilo+ngrid_per_vol-1)
                        jhi = min(ny, jlo+ngrid_per_vol-1)
                        khi = min(nz, klo+ngrid_per_vol-1)
                        
                        vol_xlo = xlo + (ilo-1)*hx
                        vol_ylo = ylo + (jlo-1)*hy
                        vol_zlo = zlo + (klo-1)*hz
                        
                        ngx = (ihi - ilo + 1)
                        ngy = (jhi - jlo + 1)
                        ngz = (khi - klo + 1)
                        
c$$$                        call fastj_copy3d(nx,ny,nz,
c$$$     $                       ilo,ihi,jlo,jhi,klo,khi,
c$$$     $                       dbl_mb(k_grid), dbl_mb(k_gtmp))

c$$$                        ijk = k_gtmp - ilo
c$$$                        do k = klo, khi
c$$$                           do j = jlo, jhi
c$$$                              jkgrid = k_grid-1 + nx*(j-1 + ny*(k-1))
c$$$                              do i = ilo, ihi
c$$$                                 dbl_mb(i + ijk) = dbl_mb(i + jkgrid)
c$$$                              enddo
c$$$                              ijk = ijk + (ihi - ilo + 1)
c$$$                           enddo
c$$$                        enddo
                     endif
c
                     do j = int_mb(k_xyz), int_mb(k_xyz+1)
                        jj = k_coords_r + (j-1)*3
c
                        x = dbl_mb(jj  )
                        y = dbl_mb(jj+1)
                        z = dbl_mb(jj+2)

                        order = minorder
                        if (order .le. 9) then
                           ninterp = order + 1
                        else
                           ninterp = order + 3
                        endif
                        if (ninterp.gt.ngx .or.
     $                       ninterp.gt.ngy .or. ninterp.gt.ngz)
     $                       call errquit('aajsflk ',
     $                       ninterp + 100*(ngx+100*(ngy+100*ngz)),
     &       UNKNOWN_ERR)
c$$$                        potgrid1 = tn_interp_3d_point(
c$$$     $                       dbl_mb(k_gtmp), ngx, ngy, ngz, hx, hy, hz,
c$$$     $                       x-vol_xlo, y-vol_ylo, z-vol_zlo, 
c$$$     $                       ninterp, order, 
c$$$     $                          dbl_mb(k_fit+numfit(order-1)))
                        potgrid1 = tn_interp_3d_point(
     $                       dbl_mb(k_grid), nx, ny, nz, hx, hy, hz,
     $                       x-xlo, y-ylo, z-zlo, 
     $                       ninterp, order, 
     $                          dbl_mb(k_fit+numfit(order-1)))
                        do order = minorder+stride,maxorder, stride
                           if (order .le. 9) then
                              ninterp = order + 1
                           else
                              ninterp = order + 3
                           endif
                           if (ninterp.gt.ngx .or.
     $                          ninterp.gt.ngy .or. ninterp.gt.ngz)
     $                          call errquit('aajsflk ',
     $                          ninterp + 100*(ngx+100*(ngy+100*ngz)),
     &       UNKNOWN_ERR)
c$$$                           potgrid = tn_interp_3d_point(
c$$$     $                          dbl_mb(k_gtmp), ngx, ngy, ngz,hx, hy,hz,
c$$$     $                          x-vol_xlo, y-vol_ylo, z-vol_zlo, 
c$$$     $                          ninterp, order,
c$$$     $                          dbl_mb(k_fit+numfit(order-1)))
                        potgrid = tn_interp_3d_point(
     $                       dbl_mb(k_grid), nx, ny, nz, hx, hy, hz,
     $                       x-xlo, y-ylo, z-zlo, 
     $                       ninterp, order, 
     $                          dbl_mb(k_fit+numfit(order-1)))

                           if (abs(potgrid-potgrid1).lt.0.1d0*eps) 
     $                          goto 31
                           potgrid1 = potgrid
                        enddo
*                        write(6,*) ' Used maxorder ', x, y, z, potgrid
                        order = order - stride
 31                     continue
                        obin(order) = obin(order) + 1
                        tbin(order) = tbin(order) + 2*(ninterp)**3 +
     $                       3*order + 22 + ninterp*(order+1)*2
*                        potgrid = tn_interp_3d_point(
*     $                       dbl_mb(k_grid), nx, ny, nz, hx, hy, hz,
*     $                       x-xlo, y-ylo, z-zlo, ninterp, order)
c
                        dbl_mb(k_v_r+j-1) = dbl_mb(k_v_r+j-1) - potgrid
                     enddo
                  ind = ind + 1
                  enddo
               enddo
            enddo
            tgrid = util_cpusec() - tgrid
         endif
c
         tfar = 0.0d0
         if (level .eq. depth) then
c     
c     Evaluate the potential from the far field
c
            tfar = util_cpusec()
            ind = 0
            do kz = 0, twol-1
               do ky = 0, twol-1
                  do kx = 0, twol-1
                     k_xyz = k_list_info + 2*(kx + twol*(ky + twol*kz))
                     if (int_mb(k_xyz) .le. int_mb(k_xyz+1)) then
c     
c     Walk up tree until find cell with local potential
c     
                        kxp = kx
                        kyp = ky
                        kzp = kz
                        rtwol = 1.0d0 / twol
                        do levelp = level, 2, -1
                           call fmm_tree_get_cell(tree, levelp, 
     $                          kxp, kyp, kzp, q, lq, luse)
                           if (luse .ge. 0) goto 666
                           kxp = kxp / 2
                           kyp = kyp / 2
                           kzp = kzp / 2
                           rtwol = rtwol + rtwol
                        enddo
 666                    continue
                        if (luse .ge. 0) then
                           hxp = (xhi-xlo)*rtwol
                           hyp = (yhi-ylo)*rtwol
                           hzp = (zhi-zlo)*rtwol
                           xo  = (dble(kxp)+0.5d0)*hxp + xlo ! Coords of box center
                           yo  = (dble(kyp)+0.5d0)*hyp + ylo
                           zo  = (dble(kzp)+0.5d0)*hzp + zlo
c     
                           do j = int_mb(k_xyz), int_mb(k_xyz+1)
                              jj = k_coords_r + (j-1)*3
                              x = dbl_mb(jj  )
                              y = dbl_mb(jj+1)
                              z = dbl_mb(jj+2)
*                              potfmm = fmm_potential(
*     $                             depth,lmax,tree,
*     $                             xhi-xlo, yhi-ylo, zhi-zlo,
*     $                             x-xlo, y-ylo, z-zlo)
                              
                              potfmm = xlm_local_potential(
     $                             x-xo,y-yo,z-zo,
     $                             q,lq,luse)

                              dbl_mb(k_v_r+j-1) = dbl_mb(k_v_r+j-1) - 
     $                             potfmm
                           enddo
                        endif
                     endif
                     ind = ind + 1
                  enddo
               enddo
            enddo
            tfar = util_cpusec() - tfar
         endif
c
c     Add the potential from this level back into the output order
c
         do j = 0,n-1
*            write(6,*) ' mapping ', j, int_mb(k_map+j), 
*     $           v(int_mb(k_map+j)), dbl_mb(k_v_r+j)
            v(int_mb(k_map+j)) = v(int_mb(k_map+j)) + dbl_mb(k_v_r+j)
         enddo
c
         used = util_cpusec() - start
         write(6,781) level, used, tfar, tgrid
 781     format(' Level ', i2,': used=', f8.2, '  far=',f8.2,
     $        '  grid=',f8.2)
      enddo                     ! level
c
      write(6,*) ' obin '
      call doutput(obin,1,34,1,1,34,1,1)
      write(6,*) ' tbin '
      call doutput(tbin,1,34,1,1,34,1,1)
c
      if (.not. ma_chop_stack(l_map)) call errquit('sigh',0, MA_ERR)
c
      end
      subroutine fastj_unpack_data(data,xlo,xhi,ylo,yhi,zlo,zhi,
     $     hx,hy,hz,k_grid,l_grid,tree,
     $     depth,lmax,nx,ny,nz,eps,l_list,k_list,l_pairs2,k_pairs2,
     $     npairs,g_dens,basis)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
c
      integer data(2)
      double precision xlo,xhi,ylo,yhi,zlo,zhi,eps,hx,hy,hz
      integer k_grid,depth,lmax,nx,ny,nz,k_r,k_i,l_grid
      integer l_list,k_list,l_pairs2,k_pairs2, npairs, g_dens, basis
      integer tree(300)
c
***      write(6,*) ' UNPACKED handles ', data(1), data(2)
      if (.not. ma_get_index(data(1),k_r))
     $     call errquit('fastj_unpack: invalid data',1, MA_ERR)
      if (.not. ma_get_index(data(2),k_i))
     $     call errquit('fastj_pack: invalid data',2, MA_ERR)
c
      xlo     =  dbl_mb(k_r   )
      xhi     =  dbl_mb(k_r+ 1)
      ylo     =  dbl_mb(k_r+ 2)
      yhi     =  dbl_mb(k_r+ 3)
      zlo     =  dbl_mb(k_r+ 4)
      zhi     =  dbl_mb(k_r+ 5)
      eps     =  dbl_mb(k_r+ 6)
      hx      =  dbl_mb(k_r+ 7)
      hy      =  dbl_mb(k_r+ 8)
      hz      =  dbl_mb(k_r+ 9)
                               
      k_grid  =  int_mb(k_i   )
      l_grid  =  int_mb(k_i+ 1)
*      k_poles =  int_mb(k_i+ 2)
*      l_poles =  int_mb(k_i+ 3)
      depth   =  int_mb(k_i+ 4)
      lmax    =  int_mb(k_i+ 5)
      nx      =  int_mb(k_i+ 6)
      ny      =  int_mb(k_i+ 7)
      nz      =  int_mb(k_i+ 8)
      l_list  =  int_mb(k_i+ 9)
      k_list  =  int_mb(k_i+10)
      l_pairs2=  int_mb(k_i+11)
      k_pairs2=  int_mb(k_i+12)
      npairs  =  int_mb(k_i+13)
      g_dens  =  int_mb(k_i+14)
      basis   =  int_mb(k_i+15)
      call icopy(300, int_mb(k_i+16), 1, tree, 1)
c
      end
      subroutine fastj_pack_data(data,xlo,xhi,ylo,yhi,zlo,zhi,
     $     hx, hy, hz, k_grid, l_grid, tree,
     $     depth,lmax,nx,ny,nz,eps,l_list,k_list,l_pairs2,k_pairs2,
     $     npairs, g_dens,basis)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
c
      integer data(2)
      double precision xlo,xhi,ylo,yhi,zlo,zhi,eps,hx,hy,hz
      integer k_grid, depth,lmax,nx,ny,nz, l_grid
      integer l_list,k_list,l_pairs2,k_pairs2, npairs, g_dens, basis
      integer tree(300)
c
      integer l_r, k_r, l_i, k_i
c
      if (.not. ma_alloc_get(mt_dbl,10, 'fastjr', l_r, k_r))
     $     call errquit('fastj_pack: unable to alloc data',1, MA_ERR)
      if (.not. ma_alloc_get(mt_int,300+16, 'fastji', l_i, k_i))
     $     call errquit('fastj_pack: unable to alloc data',2, MA_ERR)
c
      data(1) = l_r
      data(2) = l_i
c
      dbl_mb(k_r   ) = xlo
      dbl_mb(k_r+ 1) = xhi
      dbl_mb(k_r+ 2) = ylo
      dbl_mb(k_r+ 3) = yhi
      dbl_mb(k_r+ 4) = zlo
      dbl_mb(k_r+ 5) = zhi
      dbl_mb(k_r+ 6) = eps
      dbl_mb(k_r+ 7) = hx
      dbl_mb(k_r+ 8) = hy
      dbl_mb(k_r+ 9) = hz

      int_mb(k_i   ) = k_grid
      int_mb(k_i+ 1) = l_grid
*      int_mb(k_i+ 2) = k_poles
*      int_mb(k_i+ 3) = l_poles
      int_mb(k_i+ 4) = depth
      int_mb(k_i+ 5) = lmax
      int_mb(k_i+ 6) = nx
      int_mb(k_i+ 7) = ny
      int_mb(k_i+ 8) = nz
      int_mb(k_i+ 9) = l_list
      int_mb(k_i+10) = k_list
      int_mb(k_i+11) = l_pairs2
      int_mb(k_i+12) = k_pairs2
      int_mb(k_i+13) = npairs
      int_mb(k_i+14) = g_dens
      int_mb(k_i+15) = basis
      call icopy(300, tree, 1, int_mb(k_i+16), 1)
c
      end
      subroutine pot_prim_shell_pair_multipoles(
     $     dens,
     $     itype, ri, numi,
     $     jtype, rj, numj,
     $     dinv, lmaxd,
     $     a, b, c, alpha, prefac,
     $     mpoles, maxl)
      implicit none
#include "errquit.fh"
c
c     Given a block of the density matrix over a pair of CARTESIAN
c     PRIMITIVE shells and a detailed description of those shells
c     compute the multipole moments about the given coords (a,b,c)
c     which MUST be the GPT center.
c
c     .   mpoles() = multipole moments about center WITHOUT field terms
c
c     If radius is returned as 0.0d0 then the is no significant overlap
c     between the shells and the entire pair may be neglected.
c     
c     If spherical basis functions are used, the density matrix block
c     must first be transformed using dens_cart_to_sph().
c
c     Input shell parameters (numi/j = dimension of cartesian shells)
c
      integer itype, numi
      integer jtype, numj
      double precision ri(3)
      double precision rj(3)
c
      double precision dens(numi,numj) ! Cartesian density block
      integer lmaxd
      double precision dinv(*)  ! [input] dinv from xlm_coeff_inv
      integer maxl
      double precision mpoles(-maxl:maxl,0:maxl) ! [output] mpoles about center
      double precision a, b, c, alpha, prefac ! [input]
c
c     Local variables
c
      integer lmax, lmax2, lenprod
      parameter (lmax = 5)      ! Max. angular momentum of shell
      parameter (lmax2 = 2*lmax)
      parameter (lenprod = ((lmax+1)*(lmax+2)*(lmax+3))/6)
      double precision df(lmax2+1)
      double precision work(lenprod**2), ndens(lenprod**2) 
      double precision qn(-lmax2:lmax2,0:lmax2) ! [scratch]
c
      integer ijtype
      integer i, l, m, n, nlm
      double precision factor, pi, rootpi, scale
      double precision radial
c
c     Externals
c
      double precision double_factorial, fastj_gaussian_range
      external double_factorial, fastj_gaussian_range
c
c     Sanity checks
c
      if (itype.gt.lmax .or. jtype.gt.lmax) call errquit
     $     ('pot_shell_pair_multipole: lmax?',itype*1000+jtype,
     &       UNKNOWN_ERR)
c     
      ijtype = itype + jtype
      pi = 4.0d0*atan(1.0d0)
      rootpi = sqrt(pi)
      do i = 1, 2*ijtype + 1
         df(i) = double_factorial(i)
      end do
c
c     Compute the multipoles about the natrual center
c
      call xlm_norm(ijtype, qn, lmax2)
      do l = 0, ijtype
         do m = -l,l
            qn(m,l) = 1.0d0/(qn(m,l)*qn(m,l))
         end do
      end do
c
c     Translate density to product center and reform as spherical
c
      if (ijtype .eq. 0) then
         ndens(1) = dens(1,1)
      else
         call cart_dens_trans_prod_sph(
     $        itype, ri(1), ri(2), ri(3), 
     $        jtype, rj(1), rj(2), rj(3), 
     $        a, b, c, work, dinv, lmaxd, 
     $        dens, ndens)
      end if
c
c     Form multipoles at the GPT center (a,b,c) for this shell
c     
      do l = 0, maxl
         do m = -l, l
            mpoles(m,l) = 0.0d0
         end do
      end do
      nlm = 0
      scale = prefac*rootpi/(4.0d0*alpha*sqrt(alpha))
      do n = 0, ijtype
         radial = scale
         do l = n, 0, -2
            factor = radial*df(n+l+1)
            do m = -l, l
               nlm = nlm + 1
               mpoles(m,l) = mpoles(m,l) + ndens(nlm)*factor
            end do
            radial = radial * (2.0d0*alpha)
         end do
         scale = scale / (2.0d0*alpha)
      end do
      do l = 0, ijtype
         do m = -l, l
            mpoles(m,l) = mpoles(m,l)*qn(m,l)
         end do
      end do
c
      end

      subroutine dummy(a)
      end
      
      subroutine fastj_copy3d(nx,ny,nz,ilo,ihi,jlo,jhi,klo,khi,grid,tmp)
      implicit none
      integer nx,ny,nz,ilo,ihi,jlo,jhi,klo,khi
      double precision grid(nx,ny,nz),tmp(ilo:ihi,jlo:jhi,klo:khi)
c
c     Copy a sub-volume from the grid into a dense array
c
      integer i, j, k
      do k = klo, khi
         do j = jlo, jhi
            do i = ilo, ihi
               tmp(i,j,k) = grid(i,j,k)
            enddo
         enddo
      enddo
      end
      double precision function fastj_resolution(a,n,eps)
      implicit none
      double precision a
      double precision n
      double precision eps
c
c     Estimate the required resolution by examining the
c     decay of the fourier transform of a gaussian
c
      double precision t0
      double precision t1
      double precision pi
      parameter (pi = 3.141593d0)
      double precision double_factorial
      external double_factorial
c
      t0 = sqrt(log(2d0*a/(pi*eps**4)))
      t1 = sqrt(t0**2-2d0*n*log(a)-
     $     2d0*log(double_factorial(n))+
     $     2d0*(2d0*n-1d0)*log(sqrt(a)*t0))

c$$$      t0 = sqrt(log(2d0*a/pi/eps**4)-2d0*n*log(a)-
c$$$     $     2d0*log(double_factorial(n)))
c$$$      t1 = sqrt(t0**2+2d0*(2d0*n-1d0)*log(sqrt(a)*t0))

      fastj_resolution = sqrt(a)*t1 / pi

      return
      end
 
