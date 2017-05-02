      program test
*
* $Id$
*
      implicit none
c
      double precision total_flops
      common / totalf / total_flops
c
      integer lmax, depth, lwork
      parameter (lmax=32, depth=3)
      parameter (lwork=(lmax+1)**2*2*(8**(depth+1)-1)/7)
      double precision poles((lmax+1)**2,
     $     0:2**depth-1,0:2**depth-1,0:2**depth-1)
      double precision work(lwork), poles0(-lmax:lmax,0:lmax)
      double precision fmm_potential, direct_sum, util_random
      double precision fmm_direct_potential
      external fmm_potential, fmm_direct_potential, direct_sum, 
     $     util_random
      double precision x, y, z, dimx, dimy, dimz, time
      real*4 etime, stuff(2)
c
      integer ncharge, i, ntest, j
      parameter (ncharge = 10, ntest=100)
      double precision charges(ncharge), coords(3,ncharge)
      double precision chargesr(ncharge), coordsr(3,ncharge)
      double precision pot1(ncharge), pot2(ncharge), pot3(ncharge)
      integer list_next(ncharge), map(ncharge),
     $     list_info(2,0:2**depth-1,0:2**depth-1,0:2**depth-1)
c
      double precision far, near, err, relerr, l2err, l2norm, sum
c
      call xlm_init
      total_flops = 0.0d0
      write(6,*) ' lwork = ', lwork
      write(6,*) ' poles = ', (lmax+1)**2 * 8**depth
c
      dimx = 1.0d0
      dimy = 1.0d0
      dimz = 1.0d0
      do i = 1, ncharge
         charges(i) = util_random(0)-0.5d0
         coords(1,i)= util_random(0)
         coords(2,i)= util_random(0)
         coords(3,i)= util_random(0)
*         write(6,77) charges(i), (coords(j,i),j=1,3)
* 77      format(f12.6,2x,3f12.6)
      enddo
c
      write(6,*) ' Depth ', depth-1, ' lmax ', lmax
      time = etime(stuff)
      call fmm_charges_to_poles(lmax, depth-1, dimx, dimy, dimz,
     $     poles, coords, charges, ncharge, list_info, list_next,
     $     map)
      time = etime(stuff) - time
      write(6,61) ' Time to form multipoles ', time
 61   format(1x,a,f12.6)
      time = etime(stuff)
      call fmm(lmax, depth-1, dimx, dimy, dimz, poles, .false., poles0, 
     $     work, lwork)
      time = etime(stuff) - time
      write(6,61) ' Time to perform far field summation ', time
c
      write(6,*)
      write(6,*)
      write(6,*) ' Depth ', depth, ' lmax ', lmax
c
      time = etime(stuff)
      call fmm_charges_to_poles(lmax, depth, dimx, dimy, dimz,
     $     poles, coords, charges, ncharge, list_info, list_next,
     $     map)
c
      time = etime(stuff) - time
      write(6,61) ' Time to form multipoles ', time
c     
      time = etime(stuff)
      call fmm(lmax, depth, dimx, dimy, dimz, poles, .false., poles0, 
     $     work, lwork)
      time = etime(stuff) - time
      write(6,61) ' Time to perform far field summation ', time
c
      l2err = 0.0d0
      l2norm = 0.0d0
      err = 0.0d0
      relerr = 0.0d0
      do i = 1, ntest
         x = util_random(0)
         y = util_random(0)
         z = util_random(0)
         far = fmm_potential(lmax,depth,dimx,dimy,dimz,poles,x,y,z)
         near = fmm_direct_potential(
     $        lmax, depth, dimx, dimy, dimz, x, y, z,
     $        coords, charges, ncharge, list_info, map)
         sum = direct_sum(coords, charges, ncharge, x, y, z)
         err = max(err,abs(sum-far-near))
         relerr = max(relerr, abs(sum-far-near)/sum)
         l2err = l2err + (sum-far-near)**2
         l2norm = l2norm + sum**2
      enddo
      write(6,62) ' MAX abs err = ', err
 62   format(1x,a,1p,d9.2)
      write(6,62) ' MAX rel err = ', relerr
      write(6,62) ' L2 err norm = ', sqrt(l2err/l2norm)
c
**      write(6,*) ' total flops ', total_flops
c
      time = etime(stuff)
      do j = 1, 1000
         far = fmm_potential(lmax,depth,dimx,dimy,dimz,poles,x,y,z)
         near = fmm_direct_potential(lmax, depth,
     $        dimx, dimy, dimz, x, y, z,
     $        coords, charges, ncharge, list_info, map)
      enddo
      time = etime(stuff)-time
      write(6,61) ' time per FMM potential eval ', time/1000
      time = etime(stuff)
      do j = 1, 1000
         sum = direct_sum(coords, charges, ncharge, x, y, z)
      enddo
      time = etime(stuff)-time
      write(6,61) ' time per direct potential eval ', time/1000
c
      time = etime(stuff)
      do j = 1, ncharge
         chargesr(j) = charges(map(j))
         coordsr(1,j) = coords(1,map(j))
         coordsr(2,j) = coords(2,map(j))
         coordsr(3,j) = coords(3,map(j))
      enddo
      call fmm_potential_at_charges2(lmax, depth,dimx,dimy,dimz,
     $     poles, coordsr, chargesr, ncharge, list_info, pot3)
      do j = 1, ncharge
         pot1(map(j)) = pot3(j)
      enddo
      time = etime(stuff)-time
      write(6,61) ' Time to evaluate all FMM potentials ', time
c
      time = etime(stuff)
      call direct_sum_at_charges(coords, charges, ncharge, pot2)
      time = etime(stuff)-time
      write(6,61) ' Time to evaluate all direct potentials ', time
c
      l2err = 0.0d0
      l2norm = 0.0d0
      err = 0.0d0
      do i = 1, ncharge
         err = max(err,abs(pot1(i)-pot2(i)))
         l2err = l2err + (pot1(i)-pot2(i))**2
         l2norm = l2norm + pot1(i)**2
      enddo
      write(6,62) ' MAX abs err at charges = ', err
      write(6,62) ' L2 err norm at charges = ', sqrt(l2err/l2norm)
c 
*      write(6,*) ' TOP LEVEL MULTIPOLES '
*      call xlm_print(lmax,poles0,lmax)
c
      end
      double precision function fmm_potential(lmax, depth,
     $     dimx, dimy, dimz, potential, x, y, z)
      implicit none
#include "errquit.fh"
      integer lmax, depth
      double precision dimx, dimy, dimz, x, y, z
      double precision potential((lmax+1)**2,
     $     0:2**depth-1,0:2**depth-1,0:2**depth-1)
c
c     Given a potential previously evaluated by the routine fmm
c     and augmented appropriatly for the near field terms,
c     evaluate the potential at the point (x,y,z).  The volume
c     is assumed to be defined on (0:dimx,0:dimy,0:dimz)
c
      integer kx, ky, kz
      double precision hx, hy, hz, xo, yo, zo
      integer lq
      parameter (lq = FMM_LMAX)
      double precision q(-lq:lq,0:lq)
      double precision xlm_local_potential
      external xlm_local_potential
c
      hx = dimx/dble(2**depth)  ! Size of boxes
      hy = dimy/dble(2**depth)
      hz = dimz/dble(2**depth)
      kx = int(x/hx) ! Index of box containing point
      ky = int(y/hy) ! Index of box containing point
      kz = int(z/hz) ! Index of box containing point
      xo = (dble(kx)+0.5d0)*hx   ! Coords of box center
      yo = (dble(ky)+0.5d0)*hy
      zo = (dble(kz)+0.5d0)*hz
c
      if (kx.lt.0 .or. kx.ge.2**depth) call errquit('x?',0, FMM_ERR)
      if (ky.lt.0 .or. ky.ge.2**depth) call errquit('y?',0, FMM_ERR)
      if (kz.lt.0 .or. kz.ge.2**depth) call errquit('z?',0, FMM_ERR)

*      if (kx.lt.2 .or. kx.ge.2**depth-2) call errquit('x2?',0)
*      if (ky.lt.2 .or. ky.ge.2**depth-2) call errquit('y2?',0)
*      if (kz.lt.2 .or. kz.ge.2**depth-2) call errquit('z2?',0)

c
      call xlm_unpack(lmax, potential(1,kx,ky,kz), q, lq)
      fmm_potential = xlm_local_potential(x-xo,y-yo,z-zo,
     $     q,lq,lmax)
c
      end
      double precision function direct_sum(coords, charges, ncharge,
     $     x, y, z)
      implicit none
      integer ncharge
      double precision coords(3,*), charges(*), x, y, z
c
c     Return the Coulomb potential at the specified point due to the
c     given charges
c
      integer i
      double precision sum
      sum = 0.0d0
      do i = 1, ncharge
         sum = sum + charges(i)/sqrt(
     $        (coords(1,i)-x)**2 +
     $        (coords(2,i)-y)**2 + 
     $        (coords(3,i)-z)**2)
      end do
c
      direct_sum = sum
c
      end
      subroutine direct_sum_at_charges(coords, charges, ncharge, pot)
      implicit none
      integer ncharge
      double precision coords(3,*), charges(*), x, y, z, pot(ncharge)
c
c     Evaluate the potential at each charge due to the other charges
c
      integer i, j
      double precision sum, rdist, xx, yy, zz
c
      do j = 1, ncharge
         pot(j) = 0.0d0
      enddo
c
      do j = 1, ncharge
         sum = 0.0d0
         x = coords(1,j)
         y = coords(2,j)
         z = coords(3,j)
         do i = 1, j-1
            xx = coords(1,i)-x
            yy = coords(2,i)-y
            zz = coords(3,i)-z
            rdist = 1.0d0/sqrt(xx*xx+yy*yy+zz*zz)
            sum = sum + charges(i)*rdist
            pot(i) = pot(i) + charges(j)*rdist
         enddo
         pot(j) = pot(j) + sum
      end do
c
      end
      subroutine fmm_charges_to_poles(lmax, depth, dimx, dimy, dimz, 
     $     poles, coords, charges, ncharge, list_info, list_next,
     $     map)
      implicit none
#include "errquit.fh"
      integer lmax, depth, ncharge
      double precision dimx, dimy, dimz
      double precision poles((lmax+1)**2,
     $     0:2**depth-1,0:2**depth-1,0:2**depth-1)
      double precision charges(ncharge), coords(3,ncharge)
      integer list_info(2,0:2**depth-1,0:2**depth-1,0:2**depth-1),
     $     list_next(ncharge), map(ncharge)
c
c     Given a set of charges at the specified coordinates, 
c     build the input for FMM by migrating the sources
c     to the nearest center on the grid.  The coordinates are
c     assumed to be relative to the frame of the box which 
c     is defined on (0:dimx,0:dimy,0:dimz)
c
c     Initially build a linked list associating charges with each box.
c
c     list_info(1,box)  = first charge in box (or 0 if none)
c     list_next(charge) = next charge in box (or 0 if none)
c
c     Then construct a map that reorders the charges so that
c
c     list_info(1,box) = first reordered charge in box
c     list_info(2,box) = last reordered charge in box
c     map(reordered charge) = actual charge
c
      integer total
      integer kx, ky, kz, i, l, m
      double precision hx, hy, hz, xo, yo, zo, x, y, z
      integer lq
      parameter (lq = FMM_LMAX)
      double precision q(-lq:lq,0:lq)
      integer max_per_box, min_per_box, ind, n, twod
      double precision avgsq_per_box
c
      if (lmax .gt. lq) call errquit('fmm_ctop: lmax>lq', lmax, FMM_ERR)
      twod = 2**depth
      call xlm_init
c
      hx = dimx/dble(twod)  ! Size of boxes
      hy = dimy/dble(twod)
      hz = dimz/dble(twod)
      call dfill((lmax+1)**2*8**depth,0.0d0,poles,1)
      call ifill(8**depth * 2, 0, list_info, 1)
c
      do i = 1, ncharge
         x = coords(1,i)
         y = coords(2,i)
         z = coords(3,i)
         kx = int(x/hx)         ! Index of box containing point
         ky = int(y/hy)
         kz = int(z/hz)
         xo = (dble(kx)+0.5d0)*hx ! Coords of box center
         yo = (dble(ky)+0.5d0)*hy
         zo = (dble(kz)+0.5d0)*hz
         write(6,*) x, y, z, kx, ky, kz, depth, twod
         if (kx.lt.0 .or. kx.ge.twod) call errquit('x!',0, FMM_ERR)
         if (ky.lt.0 .or. ky.ge.twod) call errquit('y!',0, FMM_ERR)
         if (kz.lt.0 .or. kz.ge.twod) call errquit('z!',0, FMM_ERR)
c
         call xlm(lmax, x-xo, y-yo, z-zo, q, lq)
         do l = 0, lmax
            do m = -l, l
               q(m,l) = q(m,l)*charges(i)
            end do
         end do
*         call xlm_translate(lmax, xo-x, yo-y, zo-z, q, lq)
         call xlm_accumulate_to_packed
     $        (lmax,poles(1,kx,ky,kz),q,lq)
c
         list_next(i) = list_info(1,kx,ky,kz)
         list_info(1,kx,ky,kz) = i
      end do
c
      max_per_box = 0
      min_per_box = ncharge
      avgsq_per_box = 0
      total = 0
c      
      do kz = 0, twod-1
         do ky = 0, twod-1
            do kx = 0, twod-1
               ind = list_info(1,kx,ky,kz)
               list_info(1,kx,ky,kz) = total+1
               n = 0
 10            if (ind .gt. 0) then
                  n = n + 1
                  total = total + 1
                  map(total) = ind
                  ind = list_next(ind)
                  goto 10
               endif
               list_info(2,kx,ky,kz) = total
               if (list_info(2,kx,ky,kz) .lt. list_info(1,kx,ky,kz))then
                  list_info(1,kx,ky,kz) = 0
                  list_info(2,kx,ky,kz) = -1
               endif
               max_per_box = max(n,max_per_box)
               min_per_box = min(n,min_per_box)
               avgsq_per_box = avgsq_per_box + n*n
            enddo
         enddo
      enddo
c
      if (total .ne. ncharge) then
         write(6,*) ' total, ncharge ', total, ncharge
         stop 23
      endif
      write(6,1) max_per_box, min_per_box, ncharge/8**depth, 
     $     sqrt(avgsq_per_box/8.0d0**depth - (ncharge/8.0d0**depth)**2)
 1    format(' Particles/box:',
     $     ' max =', i4, '; min =', i4, '; mean =', i4,
     $     '; rms =', f6.1)
c
      end
      subroutine fmm_points_to_boxes(depth, n, coords,
     $     xlo, xhi, ylo, yhi, zlo, zhi,
     $     list_next, list_info, coords_r, map)
      implicit none
#include "errquit.fh"
c
      integer depth, n
      double precision coords(3,n), coords_r(3,n)
      double precision xlo, xhi, ylo, yhi, zlo, zhi
      integer list_next(n)      ! [scratch]
      integer list_info(2,0:2**depth-1,0:2**depth-1,0:2**depth-1)
      integer map(n)
c
c     Given a list of points in the given solution volume assign
c     the points to the FMM octree of given depth. 
c
c     Initially build a linked list associating points with each box.
c
c     list_info(1,box)  = first point in box (or 0 if none)
c     list_next(point) = next point in box (or 0 if none)
c
c     Then construct a map that reorders the points so that we return
c
c     list_info(1,box) = first reordered point in box, or 0 if none
c     list_info(2,box) = last reordered point in box, or -1 if none
c     map(reordered point) = index of actual point in coords()
c     coords_r(3,n)    = reordered point
c
c     Values computed at the reordered points can be put
c     back into the original order using
c
c     do j = 1, n
c     .    value(map(j)) = value_r(j)
c     enddo
c
      integer total
      integer kx, ky, kz, i
      double precision hx, hy, hz, xo, yo, zo, x, y, z
      integer max_per_box, min_per_box, ind, twod, ninbox
      double precision avgsq_per_box
c
      twod = 2**depth
c
      hx = (xhi-xlo)/dble(twod)  ! Size of boxes
      hy = (yhi-ylo)/dble(twod)
      hz = (yhi-ylo)/dble(twod)
      call ifill(8**depth * 2, 0, list_info, 1)
c
      do i = 1, n
         x = coords(1,i)
         y = coords(2,i)
         z = coords(3,i)
         kx = int((x-xlo)/hx)         ! Index of box containing point
         ky = int((y-ylo)/hy)
         kz = int((z-zlo)/hz)
         xo = xlo + (dble(kx)+0.5d0)*hx ! Coords of box center
         yo = ylo + (dble(ky)+0.5d0)*hy
         zo = zlo + (dble(kz)+0.5d0)*hz
         if (kx.lt.0 .or. kx.ge.twod) call errquit('fmm:p2b: x!',0,
     &       FMM_ERR)
         if (ky.lt.0 .or. ky.ge.twod) call errquit('fmm:p2b: y!',0,
     &       FMM_ERR)
         if (kz.lt.0 .or. kz.ge.twod) call errquit('fmm:p2b: z!',0,
     &       FMM_ERR)
c
         list_next(i) = list_info(1,kx,ky,kz)
         list_info(1,kx,ky,kz) = i
      end do
c
      max_per_box = 0
      min_per_box = n
      avgsq_per_box = 0
      total = 0
c      
      do kz = 0, twod-1
         do ky = 0, twod-1
            do kx = 0, twod-1
               ind = list_info(1,kx,ky,kz)
               list_info(1,kx,ky,kz) = total+1
               ninbox = 0
 10            if (ind .gt. 0) then
                  ninbox = ninbox + 1
                  total = total + 1
                  map(total) = ind
                  ind = list_next(ind)
                  goto 10
               endif
               list_info(2,kx,ky,kz) = total
               if (list_info(2,kx,ky,kz) .lt. list_info(1,kx,ky,kz))then
                  list_info(1,kx,ky,kz) = 0
                  list_info(2,kx,ky,kz) = -1
               endif
               max_per_box = max(ninbox,max_per_box)
               min_per_box = min(ninbox,min_per_box)
               avgsq_per_box = avgsq_per_box + dble(ninbox)**2
            enddo
         enddo
      enddo
c
      do i = 1, n
         coords_r(1,i) = coords(1,map(i))
         coords_r(2,i) = coords(2,map(i))
         coords_r(3,i) = coords(3,map(i))
      enddo
c
      write(6,1) depth, max_per_box, min_per_box, n/8**depth, 
     $     sqrt(abs(avgsq_per_box/8.0d0**depth - (n/8.0d0**depth)**2))
 1    format(' Level ', i2,' particles/box:',
     $     ' max =', i6, '; min =', i6, '; mean =', i6,
     $     '; rms =', 1p,d9.1)
c
      end
      subroutine fmm_add_values_from_near_points(
     $     depth, kx, ky, kz,
     $     n, v_r, list_info,
     $     v_near, oldnumr)
      implicit none
#include "errquit.fh"
#include "fmmP.fh"
      integer kx, ky, kz, depth, n, numr, oldnumr
      double precision v_r(n), v_near(*)
      integer list_info(2,0:2**depth-1,0:2**depth-1,0:2**depth-1)
c
c     Take values computed at the points formed by calling
c     fmm_make_list_of_near_points() and add them into the
c     values at coords_r(), the reordered points from 
c     fmm_points_to_boxes().
c
      integer kxn, kyn, kzn, dkx, dky, dkz, twod, ninbox, i, ii
c
      twod = 2**depth
      if (kx.lt.0 .or. kx.ge.twod) call errquit('x?',0, FMM_ERR)
      if (ky.lt.0 .or. ky.ge.twod) call errquit('y?',0, FMM_ERR)
      if (kz.lt.0 .or. kz.ge.twod) call errquit('z?',0, FMM_ERR)
c
      numr = 1
      do dkz = -well_separated, well_separated
         kzn = kz + dkz
         if (kzn.ge.0 .and. kzn.lt.twod) then
            do dky = -well_separated, well_separated
               kyn = ky + dky
               if (kyn.ge.0 .and. kyn.lt.twod) then
                  do dkx = -well_separated, well_separated
                     kxn = kx + dkx
                     if (kxn.ge.0 .and. kxn.lt.twod) then
                        ninbox = list_info(2,kxn,kyn,kzn) - 
     $                       list_info(1,kxn,kyn,kzn) + 1
                        if (ninbox .gt. 0) then
                           ii = list_info(1,kxn,kyn,kzn) 
                           do i = 0, ninbox-1
                              v_r(ii+i) = v_r(ii+i) + v_near(i+numr)
                           enddo
                           numr = numr + ninbox
                        endif
                     endif
                  enddo
               endif
            enddo
         endif
      enddo
      numr = numr - 1
      if (numr .ne. oldnumr) call errquit('fmm:avfnp confused',numr,
     &       FMM_ERR)
c
      end
      subroutine fmm_make_list_of_near_points(depth, kx, ky, kz, 
     $     n, coords_r, list_info, r, numr)
      implicit none
#include "errquit.fh"
#include "fmmP.fh"
      integer kx, ky, kz, depth, n, numr
      double precision coords_r(3,n), r(3,*)
      integer list_info(2,0:2**depth-1,0:2**depth-1,0:2**depth-1)
c
c     From the points that have been ordered by fmm_points_to_boxes
c     return a list of points at which the box (kx, ky, kz) must
c     compute its near field potential.
c
c     r(3,1:numr) and numr are returned.  r() should probably be
c     the same dimension as coords to be safe.
c
      integer kxn, kyn, kzn, dkx, dky, dkz, twod, ninbox
c
      twod = 2**depth
      if (kx.lt.0 .or. kx.ge.twod) call errquit('x?',0, FMM_ERR)
      if (ky.lt.0 .or. ky.ge.twod) call errquit('y?',0, FMM_ERR)
      if (kz.lt.0 .or. kz.ge.twod) call errquit('z?',0,
     &       FMM_ERR)
c
      numr = 1
      do dkz = -well_separated, well_separated
         kzn = kz + dkz
         if (kzn.ge.0 .and. kzn.lt.twod) then
            do dky = -well_separated, well_separated
               kyn = ky + dky
               if (kyn.ge.0 .and. kyn.lt.twod) then
                  do dkx = -well_separated, well_separated
                     kxn = kx + dkx
                     if (kxn.ge.0 .and. kxn.lt.twod) then
                        ninbox = list_info(2,kxn,kyn,kzn) - 
     $                       list_info(1,kxn,kyn,kzn) + 1
                        if (ninbox .gt. 0) then
                           call dcopy(3*ninbox,
     $                          coords_r(1,list_info(1,kxn,kyn,kzn)), 1,
     $                          r(1,numr), 1)
                           numr = numr + ninbox
                        endif
                     endif
                  enddo
               endif
            enddo
         endif
      enddo
      numr = numr - 1
c
      end
      double precision function fmm_direct_potential(
     $     lmax, depth,
     $     dimx, dimy, dimz, x, y, z,
     $     coords, charges, ncharge, list_info, map)
      implicit none
#include "errquit.fh"
#include "fmmP.fh"
      integer lmax, depth
      double precision dimx, dimy, dimz, x, y, z
      integer ncharge
      double precision charges(ncharge), coords(3,ncharge)
      integer list_info(2,0:2**depth-1,0:2**depth-1,0:2**depth-1)
      integer map(ncharge)
c
c     Evaluate the direct local potential at the given point
c     ... i.e., directly sum interactions within the box
c     containing the point and immediate neighbours
c
      integer kx, ky, kz, kxn, kyn, kzn, dkx, dky, dkz, twod
      integer i, ii
      double precision hx, hy, hz, xx, yy, zz, potential
c
      twod = 2**depth
      hx = dimx/dble(twod)  ! Size of boxes
      hy = dimy/dble(twod)
      hz = dimz/dble(twod)
      kx = int(x/hx) ! Index of box containing point
      ky = int(y/hy) ! Index of box containing point
      kz = int(z/hz) ! Index of box containing point
c
      if (kx.lt.0 .or. kx.ge.twod) call errquit('x?',0, FMM_ERR)
      if (ky.lt.0 .or. ky.ge.twod) call errquit('y?',0,
     &       FMM_ERR)
      if (kz.lt.0 .or. kz.ge.twod) call errquit('z?',0, FMM_ERR)
c
      potential = 0.0d0
c
      do dkz = -well_separated, well_separated
         kzn = kz + dkz
         if (kzn.ge.0 .and. kzn.lt.twod) then
            do dky = -well_separated, well_separated
               kyn = ky + dky
               if (kyn.ge.0 .and. kyn.lt.twod) then
                  do dkx = -well_separated, well_separated
                     kxn = kx + dkx
                     if (kxn.ge.0 .and. kxn.lt.twod) then
                        do ii = list_info(1,kxn,kyn,kzn), 
     $                       list_info(2,kxn,kyn,kzn)
                           i = map(ii)
                           xx = coords(1,i) - x
                           yy = coords(2,i) - y
                           zz = coords(3,i) - z
                           potential = potential + charges(i)/
     $                          sqrt(xx*xx + yy*yy + zz*zz)
                        enddo
                     endif
                  enddo
               endif
            enddo
         endif
      enddo
c
      fmm_direct_potential = potential
c
      end
      subroutine fmm_potential_at_charges2(lmax, depth,dimx,dimy,dimz,
     $     potential, coords, charges, ncharge, list_info, pot)
      implicit none
#include "errquit.fh"
#include "fmmP.fh"
      integer lmax, depth
      double precision dimx, dimy, dimz
      double precision potential((lmax+1)**2,
     $     0:2**depth-1,0:2**depth-1,0:2**depth-1)
      integer ncharge
      double precision charges(ncharge), coords(3,ncharge), pot(ncharge)
      integer list_info(2,0:2**depth-1,0:2**depth-1,0:2**depth-1)
c     
c     Evaluate the near and far field potentials at all charges
c     returning the result in pot
c     
      integer kx, ky, kz, i, l, m, kxn, kyn, kzn, ind
      integer twod, j, dkyhi, dkxhi
      integer kxnlo, kynlo, kznlo
      double precision hx, hy, hz, xo, yo, zo, sum, 
     $     x, y, z, xx, yy, zz, rdist

      integer lq
      parameter (lq = FMM_LMAX)
      double precision w(-lq:lq,0:lq)
      double precision xlm_local_potential
      external xlm_local_potential
c     
      if (lmax .gt. lq) call errquit('fmm_pac2: lmax>lq', lmax, FMM_ERR)
      twod = 2**depth
      hx = dimx/dble(twod)  ! Size of boxes
      hy = dimy/dble(twod)
      hz = dimz/dble(twod)
c
      do i = 1, ncharge
         pot(i) = 0.0d0
      enddo
c
      do kz = 0, twod-1
         do ky = 0, twod-1
            do kx = 0, twod-1
               xo = (dble(kx)+0.5d0)*hx ! Coords of box center
               yo = (dble(ky)+0.5d0)*hy
               zo = (dble(kz)+0.5d0)*hz
c     
c     Compute the far field term
c     
               do i = list_info(1,kx,ky,kz),list_info(2,kx,ky,kz)
                  x = coords(1,i)
                  y = coords(2,i)
                  z = coords(3,i)
                  sum = 0.0d0
                  call xlm(lmax, x-xo, y-yo, z-zo, w, lq)
                  ind = 1
                  do l = 0, lmax
                     do m = -l,l
                        sum = sum + w(m,l)*potential(ind,kx,ky,kz)
                        ind = ind + 1
                     end do
                  end do
                  pot(i) = pot(i) + sum
               enddo
c     
c     Near field terms
c
               do i = list_info(1,kx,ky,kz),list_info(2,kx,ky,kz)
                  x = coords(1,i)
                  y = coords(2,i)
                  z = coords(3,i)
                  sum = 0.0d0
                  do j = list_info(1,kx,ky,kz), i-1
                     xx = coords(1,j) - x
                     yy = coords(2,j) - y
                     zz = coords(3,j) - z
                     rdist = 1.0d0/sqrt(xx*xx+yy*yy+zz*zz)
                     sum = sum + charges(j)*rdist
                     pot(j) = pot(j) + charges(i)*rdist
                  enddo
                  pot(i) = pot(i) + sum
               enddo
c
               kxnlo = max(0,kx-well_separated)
               kynlo = max(0,ky-well_separated)
               kznlo = max(0,kz-well_separated)
               do kzn = kznlo, kz
                  dkyhi = well_separated
                  if (kz .eq. kzn) dkyhi = 0
                  do kyn = kynlo, min(twod-1,ky+dkyhi)
                     dkxhi = well_separated
                     if (kz.eq.kzn .and. ky.eq.kyn) dkxhi = -1
                     do kxn = kxnlo, min(twod-1,kx+dkxhi)
                        do j = list_info(1,kxn,kyn,kzn),
     $                       list_info(2,kxn,kyn,kzn)
                           x = coords(1,j)
                           y = coords(2,j)
                           z = coords(3,j)
                           sum = 0.0d0
                           do i = list_info(1,kx,ky,kz), 
     $                          list_info(2,kx,ky,kz)
                              xx = coords(1,i) - x
                              yy = coords(2,i) - y
                              zz = coords(3,i) - z
                              rdist = 1.0d0/sqrt(xx*xx+yy*yy+zz*zz)
                              sum = sum + charges(i)*rdist
                              pot(i) = pot(i) + charges(j)*rdist
                           enddo
                           pot(j) = pot(j) + sum
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
      end
      subroutine fmm_potential_at_charges(lmax, depth,dimx,dimy,dimz,
     $     potential, coords, charges, ncharge, list_info, pot)
      implicit none
#include "errquit.fh"
#include "fmmP.fh"
      integer lmax, depth
      double precision dimx, dimy, dimz
      double precision potential((lmax+1)**2,
     $     0:2**depth-1,0:2**depth-1,0:2**depth-1)
      integer ncharge
      double precision charges(ncharge), coords(3,ncharge), pot(ncharge)
      integer list_info(2,0:2**depth-1,0:2**depth-1,0:2**depth-1)
c     
c     Evaluate the near and far field potentials at all charges
c     returning the result in pot
c     
      integer kx, ky, kz, i, l, m, kxn, kyn, kzn, ind
      integer twod, j, dkyhi, dkxhi
      integer kxnlo, kynlo, kznlo
      double precision hx, hy, hz, xo, yo, zo, sum, 
     $     x, y, z, xx, yy, zz, rdist

      integer lq
      parameter (lq = FMM_LMAX)
      double precision w(-lq:lq,0:lq)
      double precision xlm_local_potential
      external xlm_local_potential
c     
      if (lmax .gt. lq) call errquit('fmm_pac: lmax>lq', lmax, FMM_ERR)
      twod = 2**depth
      hx = dimx/dble(twod)  ! Size of boxes
      hy = dimy/dble(twod)
      hz = dimz/dble(twod)
c
      do i = 1, ncharge
         pot(i) = 0.0d0
      enddo
c
      do i = 1, ncharge
         x = coords(1,i)
         y = coords(2,i)
         z = coords(3,i)
         kx = int(x/hx)         ! Index of box containing point
         ky = int(y/hy)         ! Index of box containing point
         kz = int(z/hz)         ! Index of box containing point
         xo = (dble(kx)+0.5d0)*hx ! Coords of box center
         yo = (dble(ky)+0.5d0)*hy
         zo = (dble(kz)+0.5d0)*hz
c
         sum = 0.0d0
c
c     Compute the far field term
c
         call xlm(lmax, x-xo, y-yo, z-zo, w, lq)
         ind = 1
         do l = 0, lmax
            do m = -l,l
               sum = sum + w(m,l)*potential(ind,kx,ky,kz)
               ind = ind + 1
            end do
         end do
c     
c     Near field terms ... avoiding current charge
c
         do j = list_info(1,kx,ky,kz), list_info(2,kx,ky,kz)
            if (j .gt. i) then
               xx = coords(1,j) - x
               yy = coords(2,j) - y
               zz = coords(3,j) - z
               rdist = 1.0d0/sqrt(xx*xx+yy*yy+zz*zz)
               sum = sum + charges(j)*rdist
               pot(j) = pot(j) + charges(i)*rdist
            endif
         enddo
c
         kxnlo = max(0,kx-well_separated)
         kynlo = max(0,ky-well_separated)
         kznlo = max(0,kz-well_separated)
         do kzn = kznlo, kz
            dkyhi = well_separated
            if (kz .eq. kzn) dkyhi = 0
            do kyn = kynlo, min(twod-1,ky+dkyhi)
               dkxhi = well_separated
               if (kz.eq.kzn .and. ky.eq.kyn) dkxhi = -1
               do kxn = kxnlo, min(twod-1,kx+dkxhi)
                  do j = list_info(1,kxn,kyn,kzn), 
     $                 list_info(2,kxn,kyn,kzn)
                     xx = coords(1,j) - x
                     yy = coords(2,j) - y
                     zz = coords(3,j) - z
                     rdist = 1.0d0/sqrt(xx*xx+yy*yy+zz*zz)
                     sum = sum + charges(j)*rdist
                     pot(j) = pot(j) + charges(i)*rdist
                  enddo
               enddo
            enddo
         enddo
c
         pot(i) = pot(i) + sum
c
      enddo
      end
      subroutine fmm(lmax, depth, dimx, dimy, dimz, 
     $     mpoles, otree, mpoles0, work, lwork)
      implicit none
#include "errquit.fh"
      double precision flops, flopsv, flopst, calls
      common /xlmflops/flops, flopsv, flopst, calls
c
c     Given a regular power-of-two grid of unnormalized, real
c     multipoles, sum the far-field interactions using FMM.
c     Optionally allow multipoles to be defined on coarser grids
c     going up the octree (otree=.true.).
c
c     The FMM results in a local spherical Taylor series associated
c     with each cell at the finest level from which the far-field 
c     potential may be computed by calling xlm_local_potenial.
c     
c     The direct interaction of each cell with its eight immediate 
c     neighbours and the interactions within a cell are NOT computed 
c     here and must be formed by other means.
c
c     Arguments
c
c     lmax  = Maximum order of multipole to use.  The error is
c     .       proportional to 2**(-lmax).
c
c     depth = level of the tree at which multipoles are provided
c     .       There will be 8**depth boxes.  Depth can be anything
c     .       greater than zero, but for efficiency should be at
c     .       least 3 and probably no more than 7.
c     .       (depth=0 means one box covers the whole domain)
c
c     dimx/y/z= The volume of the summation.  Boxes at
c     .       a given level are labelled k=0,...,2**level-1 for
c     .       each of x,y,z.  The width of each box is 
c     .       hx=dimx/2**level.  The center of each box, at which
c     .       both multipoles and local taylor series are centered,
c     .       is given by xlo+(k+0.5)*hx where xlo is the coordinate
c     .       of the bottom-left-front corner of the box in user
c     .       space.
c     
c     mpoles(lenxlm,2**depth,2**depth,2**depth)
c     .       Inputs the multipoles within each box about the box center.
c     .       The definition of the multipoles corresponds to the unnormalized
c     .       real solid spherical harmonics as computed by xlm().  The
c     .       routine xlm_coeff() can also be use to transform from simple
c     .       cartesian moments to this basis.  See xlm_pack/unpack
c     .       for the compressed storage format.
c     .       lenxlm = sum(l=0,lmax)(2*l+1) = (lmax+1)**2
c     .
c     .       Returns the spherical Taylor series expansions.  The potential
c     .       at a point within a box may be computed from the corresponding
c     .       series by
c     .           double precision q(-lq:lq,0:lmax) ! lq>= lmax
c     .           call xlm_unpack(lmax, mpoles(1,kx,ky,kz), q, lq)
c     .           u = xlm_local_potential(x,y,z,q,lq,lmax)
c     .       where x,y,z are the coordinates relative to the center of
c     .       the box, and kx,ky,kz label the box as described above.
c
c     otree = If false, then mpoles is as described above, which specifies
c     .       the multipole moments associated with just the finest level.
c     .
c     .       If true, then multipoles may be specified at all levels of
c     .       the octree and mpoles() consists of the concatenated grids
c     .       starting with the finest and will be of dimension 
c     .       lenxlm**8^(depth+1)-1)/7.  The effect of otree=false may
c     .       be emulated with otree=true by just specifying mpoles
c     .       at the finest level and providing zeroes at the other levels.
c
c     mpoles0(-lmax:lmax,0:lmax)
c     .       Returns the multipole moments at level=0 (the whole domain).
c
c     work(lwork)
c     .       Scratch workspace.  
c     .       lwork >= lenxlm*2*[sum(k=0,depth)8^k]
c     .             >= lenxlm*2*(8^(depth+1)-1)/7
c
      integer lmax, depth, lwork, mptr
      double precision dimx, dimy, dimz
      double precision mpoles((lmax+1)**2,*)
      double precision work((lmax+1)**2,*)
      double precision mpoles0(-lmax:lmax,0:lmax)
      logical otree
c
      integer lenxlm, poles, parent_poles, level,
     $     potential, parent_potential
      
      integer luse(32768)

      flops = 0
      flopsv= 0
      flopst= 0
      calls = 0

      call xlm_init

      lenxlm = (lmax+1)**2

      if (lwork/lenxlm .lt. 2*((8**(depth+1)-1)/7)) 
     $     call errquit('fmm: lwork must be ',0, 2*((8**(depth+1)-1)/7),
     &       FMM_ERR)

c     In the initial debug implementation copy all data into
c     workspace and copy results back.  Will do a memory
c     saving version later.  So for now, work is organized
c     as follows
c     level=depth   potential+mpoles
c     level=depth-1 potential+mpoles
c     ...
c     level=0       potential+mpoles
c
*****      call dfill(lwork, 0.0d0, work, 1)
      call dcopy(lenxlm*(8**depth), mpoles, 1, work, 1)
c
c     Upward pass
c
c
      parent_poles = 1
      do level = depth, 1, -1

         poles        = parent_poles
         potential    = poles + (8**level)
         parent_poles = poles + (8**level)*2
c
         if ((parent_poles+2*8**(level-1)-1)*lenxlm .gt. lwork) 
     $        call errquit('xxxxxxx?',0, FMM_ERR)
c
c     Pass multipoles up to the parent and pass short-range potential
c     to neighbours at the current level.  If we have source terms
c     on the parent level (otree=true) then copy them, otherwise
c     zero the parent poles.
c
         if (otree) then
            mptr = (parent_poles-1)/2 + 1
            call dcopy(lenxlm*(8**(level-1)), mpoles(1,mptr), 1,
     $           work(1,parent_poles), 1)
         else
            call dfill(lenxlm*(8**(level-1)), 0.0d0,
     $           work(1,parent_poles), 1)
         endif
         call dfill(lenxlm*(8**level), 0.0d0, work(1,potential), 1)
c
         call fmm_pass1(lmax, level, dimx, dimy, dimz,
     $        work(1,poles), work(1,potential), work(1,parent_poles),
     $        luse)
c
         if (level .eq. 1) then
            call xlm_unpack(lmax, work(1,parent_poles),
     $           mpoles0, lmax)
         end if
c
      end do
c
c     Downward pass ... note value of pointer potential is used from above
c     (the last level at which a local potential was computed was 1
c     ... this is actually zero since no suitable neighbours exist
c     until level 2).
c
      do level = 2, depth
         parent_potential = potential
         potential = parent_potential - 8**(level-1) - 8**level

         call fmm_pass2(lmax, level, 
     $        dimx, dimy, dimz, 
     $        work(1,parent_potential), work(1,potential))
      end do
c
      call dcopy(lenxlm*(8**depth), work(1,potential), 1, mpoles, 1)
c
*      write(6,*) ' flops calls ', flops,flopsv,flopst,calls
c
      end
      subroutine fmm_pass1(lmax, level, dimx, dimy, dimz,
     $     poles, potential, parent_poles, luse)
      implicit none
#include "errquit.fh"
#include "fmmP.fh"
      integer lmax, level
      double precision dimx, dimy, dimz
      double precision poles((lmax+1)**2,
     $     0:2**level-1,0:2**level-1,0:2**level-1)
      double precision potential((lmax+1)**2,
     $     0:2**level-1,0:2**level-1,0:2**level-1)
      double precision parent_poles((lmax+1)**2,
     $     0:2**(level-1)-1,0:2**(level-1)-1,0:2**(level-1)-1)
      integer luse(0:2**level-1,0:2**level-1,0:2**level-1)
c     
      double precision hx, hy, hz, x, y, z, xn, yn, zn, xp, yp, zp
      double precision xx, yy, zz
      integer dkx, dky, dkz, l, ln, m
      integer kx, ky, kz, kxn, kyn, kzn, kxp, kyp, kzp
      integer lq
      parameter (lq = FMM_LMAX)
      double precision q(-lq:lq,0:lq), qn(-lq:lq,0:lq)
      double precision p(-lq:lq,0:lq), pn(-lq:lq,0:lq)
      double precision scale, scale0
      integer lmaxeff, lmaxeffn, hi, lo, twol, ind, indn
      integer lusemin
      parameter (lusemin = 6)
c
c     If in_all_boxes is true, then the potential is transferred even
c     to boxes that do not contain charges.  If you only want to compute
c     the potential at the location of charges, set in_all_boxes false.
c
      logical in_all_boxes
      parameter (in_all_boxes = .true.)
c
c     Note that to localize memory references even though
c     q(), etc., are dimensioned lq>=lmax, we are telling the
c     subroutines below they are dimensioned lmax, thus the
c     arrays are NOT easily used in this routine.
c
      if (lq .lt. lmax) call errquit('fmm_p_m_u: lq ', lmax, FMM_ERR)
c
      twol = 2**level
c     
      hx = dimx/dble(twol)
      hy = dimy/dble(twol)
      hz = dimz/dble(twol)
c     
      scale0 = (well_separated+1d0)/0.87d0 - 1.0d0

      do kz = 0, twol-1
         do ky = 0, twol-1
            do kx = 0, twol-1
               call xlm_unpack(lmax,poles(1,kx,ky,kz),q,lq)
               do l = lmax, 0, -1
                  do m = -l,l
                     if (abs(q(m,l)).ne.0.0d0) goto 33
                  enddo
               enddo
 33            luse(kx,ky,kz) = l
            enddo
         enddo
      enddo
c     
      do kz = 0, twol-1
         do ky = 0, twol-1
            do kx = 0, twol-1
               x = (dble(kx)+0.5d0)*hx ! Coords of current box center
               y = (dble(ky)+0.5d0)*hy
               z = (dble(kz)+0.5d0)*hz
c     
               l = luse(kx,ky,kz)
               if (l .ge. 0) then
                  call xlm_unpack(lmax,poles(1,kx,ky,kz),q,lmax)
                  kxp = kx/2    ! Parent box indices
                  kyp = ky/2
                  kzp = kz/2
                  xp = (dble(kxp)+0.5d0)*hx*2.0d0 ! Coords of parent box
                  yp = (dble(kyp)+0.5d0)*hy*2.0d0
                  zp = (dble(kzp)+0.5d0)*hz*2.0d0
c     
c     Pass multipoles up to parent (note that this destroys q so
c     must unpack again afterwards)
c     
                  call xlm_translate(lmax, xp-x, yp-y, zp-z, q, lmax)
                  call xlm_accumulate_to_packed
     $                 (lmax,parent_poles(1,kxp,kyp,kzp),q,lmax)
               endif
c     
c     Propagate local taylor series expansion
c     
               ind = kx+twol*(ky+twol*kz)
               lo = 2*well_separated
               hi = 2*well_separated + 1
               if (l .ge. 0) call xlm_unpack
     $              (lmax,poles(1,kx,ky,kz),q,lmax)
               do kzn = max(0,int(kz/2)*2-lo),
     $              min(twol-1,int(kz/2)*2+hi)
                  dkz = abs(kzn - kz)
                  do kyn = max(0,int(ky/2)*2-lo),
     $                 min(twol-1,int(ky/2)*2+hi)
                     dky = abs(kyn - ky)
                     do kxn = max(0,int(kx/2)*2-lo),
     $                    min(twol-1,int(kx/2)*2+hi)
                        dkx = abs(kxn - kx)
                        indn = kxn+twol*(kyn+twol*kzn)
                        if ((dkx.gt.well_separated .or.
     $                       dky.gt.well_separated .or.
     $                       dkz.gt.well_separated) .and.
     $                       indn.gt.ind) then
c     
c     Limit precision if possible for distant boxes.  h*sqrt(3)/2 is the
c     distance of the corner of a cube from its center.
c     
                           ln = luse(kxn,kyn,kzn)
                           lmaxeff  = l
                           lmaxeffn = ln
                           if (ln.ge.0 .or. l.ge.0) then
                              xn = (dble(kxn)+0.5d0)*hx
                              yn = (dble(kyn)+0.5d0)*hy
                              zn = (dble(kzn)+0.5d0)*hz
                              xx = xn - x
                              yy = yn - y
                              zz = zn - z
                              scale = sqrt(xx*xx+yy*yy+zz*zz)
                              scale = scale / (0.87d0*hx) - 1.0d0
                              scale = log(scale0)/log(scale)
*                              scale = 1.0d0
                           endif
                           if (ln .ge. 0) then 
                              call xlm_unpack(lmax,
     $                             poles(1,kxn,kyn,kzn), qn, lmax)
                              lmaxeffn = min(ln,max(6,nint(lmax*scale)))
                           endif
                           if (l .ge. 0) then 
                              lmaxeff = min(l,max(6,nint(lmax*scale)))
                           endif
c
                           if (lmaxeff.ge.0 .and. lmaxeffn.ge.0) then
                              call xlm_new_multipole_to_local(lmaxeff, 
     $                             xn, yn, zn, pn, lmax,x,y,z,q,lmax)
                              call xlm_accumulate_to_packed(lmaxeff,
     $                             potential(1,kxn,kyn,kzn),pn,lmax)
                              call xlm_new_multipole_to_local(lmaxeffn, 
     $                             x, y, z, p, lmax, xn, yn, zn,qn,lmax)
                              call xlm_accumulate_to_packed(lmaxeffn,
     $                             potential(1,kx,ky,kz),p,lmax)
                           else if (lmaxeff.ge.0.and.in_all_boxes) then
                              call xlm_new_multipole_to_local(lmaxeff, 
     $                             xn, yn, zn, pn, lmax, x,y,z,q,lmax)
                              call xlm_accumulate_to_packed(lmaxeff,
     $                             potential(1,kxn,kyn,kzn),pn,lmax)
                           else if (lmaxeffn.ge.0.and.in_all_boxes) then
                              call xlm_new_multipole_to_local(lmaxeffn, 
     $                             x, y, z, p, lmax, xn, yn,zn,qn,lmax)
                              call xlm_accumulate_to_packed(lmaxeffn,
     $                             potential(1,kx,ky,kz),p,lmax)
                           endif
                        end if
                     end do
                  end do
               end do
            end do
         end do
      end do
c
      end
      subroutine fmm_pass2(lmax, level, 
     $     dimx, dimy, dimz, parent_potential, potential)
      implicit none
#include "errquit.fh"
      integer lmax, level
      double precision dimx, dimy, dimz
      double precision potential((lmax+1)**2,
     $     0:2**level-1,0:2**level-1,0:2**level-1)
      double precision parent_potential((lmax+1)**2,
     $     0:2**(level-1)-1,0:2**(level-1)-1,0:2**(level-1)-1)
c
      double precision hx, hy, hz, x, y, z, xp, yp, zp
      integer kx, ky, kz, kxp, kyp, kzp
      integer lq
      parameter (lq = FMM_LMAX)
      double precision q((2*lq+1)*(lq+1))
c
      if (lq .lt. lmax) call errquit('fmm_p_f_a: lq ', lmax, FMM_ERR)
c
      hx = dimx/dble(2**level)
      hy = dimy/dble(2**level)
      hz = dimz/dble(2**level)
c
      do kz = 0, 2**level-1
         do ky = 0, 2**level-1
            do kx = 0, 2**level-1
               x = (dble(kx)+0.5d0)*hx ! Coords of current box center
               y = (dble(ky)+0.5d0)*hy
               z = (dble(kz)+0.5d0)*hz
               kxp = kx/2       ! Parent box indices
               kyp = ky/2
               kzp = kz/2
               xp = (dble(kxp)+0.5d0)*hx*2.0d0 ! Coords of parent box
               yp = (dble(kyp)+0.5d0)*hy*2.0d0
               zp = (dble(kzp)+0.5d0)*hz*2.0d0
               call xlm_unpack(lmax, parent_potential(1,kxp,kyp,kzp), 
     $              q, lmax)
               call xlm_local_translate(lmax,x-xp,y-yp,z-zp,q,lmax)
               call xlm_accumulate_to_packed(lmax,
     $              potential(1,kx,ky,kz),q,lmax)
            end do
         end do
      end do
c
      end
c$$$      subroutine fmm_multipole_to_local_rot(lmax, 
c$$$     $     xn, yn, zn, pn, lp, x, y, z, q, lq)
c$$$      implicit none
c$$$#include "xlm.fh"
c$$$c     
c$$$      integer lmax, lq, lp
c$$$      double precision xn, yn, zn, x, y, z
c$$$      double precision q(-lq:lq,0:lmax), pn(-lp:lp,0:lp)
c$$$c
c$$$      double precision theta, phi, r, pi
c$$$      double precision qq(-maxl:maxl,0:maxl)
c$$$      integer l, m
c$$$c
c$$$c     Determine the rotation angle
c$$$c
c$$$      pi = 4.0d0*atan(1.0d0)
c$$$      call fmm_xyz_to_phi_theta(xn-x,yn-y,zn-z,phi,theta)
c$$$*      write(6,*) ' doing it ', xn, yn, zn, x, y, z, 
c$$$*     $     theta*180.0d0/pi, phi*180.0d0/pi
c$$$c      
c$$$c     Rotate the source multipoles
c$$$c
c$$$      do l = 0, lmax
c$$$         do m = -l, l
c$$$            qq(m,l) = 2.0d0*q(m,l) * (b(m,l)*c(m,l))
c$$$         enddo
c$$$         qq(0,l) = qq(0,l) * 0.5d0
c$$$      enddo
c$$$c
c$$$*      write(6,*) ' before phi ', phi
c$$$*      call xlm_print(lmax,q,lq)
c$$$      call xlm_apply_z_rotation(lmax,phi,qq,maxl)
c$$$*      write(6,*) ' after phi ', phi
c$$$*      call xlm_print(lmax,q,lq)
c$$$      call xlm_apply_y_rotation(lmax,theta,qq,maxl)
c$$$c
c$$$      do l = 0, lmax
c$$$         do m = -l, l
c$$$            qq(m,l) = 0.5d0*qq(m,l) / (b(m,l)*c(m,l))
c$$$         enddo
c$$$         qq(0,l) = qq(0,l) * 2.0d0
c$$$      enddo
c$$$c
c$$$c     Translate 0,0,r
c$$$c
c$$$      r = sqrt((xn-x)**2 + (yn-y)**2 + (zn-z)**2)
c$$$      call xlm_multipole_to_local_z(lmax, 
c$$$     $     r, pn, lp, qq, maxl)
c$$$*      call xlm_multipole_to_local(lmax, 
c$$$*     $     0.0d0, 0.0d0, r, pn, lp, 0.0d0, 0.0d0, 0.0d0, qq, maxl)
c$$$c
c$$$c     Rotate the potentials back
c$$$c     
c$$$*      write(6,*) ' before -phi ', -phi
c$$$*      call xlm_print(lmax,pn,lp)
c$$$      call xlm_apply_y_rotation(lmax,-theta,pn,lp)
c$$$*      write(6,*) ' after -phi ', -phi
c$$$*      call xlm_print(lmax,pn,lp)
c$$$      call xlm_apply_z_rotation(lmax,-phi,pn,lp)
c$$$c
c$$$      end



      subroutine fmm_xyz_to_phi_theta(xx, yy, zz,
     $     phi, theta)
      implicit none
      double precision xx, yy, zz
      double precision theta, phi
c
c     dx, dy, dz are integer coordinates of a box relative to
c     origin ... determine the spherical coords
c     
      double precision r, x, y, z
c
      x = xx
      y = yy
      z = zz
c
      if (x.eq.0.0d0 .and. y.eq.0.0d0) then
         phi = 0.0d0
         theta = 0.0d0
      else 
         r = sqrt(x*x + y*y + z*z)
         z = z/r
         if (abs(z) .gt. 1.0d0) z = sign(1.0d0,z)
         theta = acos(z)
         phi = atan2(y,x)
      endif
c
      end
