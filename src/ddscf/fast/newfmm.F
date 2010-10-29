      double precision function fmm_potential(depth,lmax,tree, 
     $     dimx, dimy, dimz, x, y, z)
*
* $Id$
*
      implicit none
#include "errquit.fh"
      integer tree(300)
      integer lmax, depth
      double precision dimx, dimy, dimz, x, y, z
c
c     Given a potential previously evaluated by the routine fmm
c     and augmented appropriatly for the near field terms,
c     evaluate the potential at the point (x,y,z).  The volume
c     is assumed to be defined on (0:dimx,0:dimy,0:dimz)
c
      integer kx, ky, kz, level
      double precision hx, hy, hz, xo, yo, zo
      integer lq, luse
      parameter (lq = FMM_LMAX)
      double precision q(-lq:lq,0:lq)
      double precision xlm_local_potential
      external xlm_local_potential
c
      double precision rtwol
c
      rtwol = 0.5d0**depth
      fmm_potential = 0.0d0
      do level = depth, 2, -1
         hx = dimx*rtwol
         hy = dimy*rtwol
         hz = dimz*rtwol
         kx = int(x/hx)         ! Index of box containing point
         ky = int(y/hy)         ! Index of box containing point
         kz = int(z/hz)         ! Index of box containing point
         xo = (dble(kx)+0.5d0)*hx ! Coords of box center
         yo = (dble(ky)+0.5d0)*hy
         zo = (dble(kz)+0.5d0)*hz
c     
*         if (kx.lt.0 .or. kx.ge.2**level) call errquit('newfmm:x?',0)
*         if (ky.lt.0 .or. ky.ge.2**level) call errquit('newfmm:y?',0)
*         if (kz.lt.0 .or. kz.ge.2**level) call errquit('newfmm:z?',0)
c     
         call fmm_tree_get_cell(tree, level, kx, ky, kz, q, lq, luse)
         if (luse .ge. 0) then
            fmm_potential = xlm_local_potential(x-xo,y-yo,z-zo,
     $           q,lq,luse)
            return
         endif
         rtwol = rtwol + rtwol
      enddo
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
      subroutine fmm_charges_to_poles(depth, lmax, tree, 
     $     dimx, dimy, dimz, 
     $     coords, charges, ncharge, list_info, list_next,
     $     map)
      implicit none
#include "errquit.fh"
      integer tree(300)
      integer lmax, depth, ncharge
      double precision dimx, dimy, dimz
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
         if (kx.lt.0 .or. kx.ge.twod) call errquit('newfmm:x!',0,
     &         FMM_ERR)
         if (ky.lt.0 .or. ky.ge.twod) call errquit('newfmm:y!',0, 
     &        FMM_ERR)
         if (kz.lt.0 .or. kz.ge.twod) call errquit('newfmm:z!',0,
     &         FMM_ERR)
c
         call xlm(lmax, x-xo, y-yo, z-zo, q, lq)
         do l = 0, lmax
            do m = -l, l
               q(m,l) = q(m,l)*charges(i)
            end do
         end do
*         write(6,*) ' Before add cell '
*         call xlm_print(lmax,q,lq)
         call fmm_tree_add_cell(tree,depth,kx,ky,kz,q,lq,lmax)
*         call dfill((2*lq+1)*(lq+1),0.0d0,q,1)
*         call fmm_tree_get_cell(tree,depth,kx,ky,kz,q,lq,luse)
*         write(6,*) ' After get cell'
*         call xlm_print(luse,q,lq)
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
         call errquit('newfmm:fmm_charges_to_poles: lost something!',0,
     &       FMM_ERR)
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
      if (kx.lt.0 .or. kx.ge.twod) call errquit('newfmm:x?',0, FMM_ERR)
      if (ky.lt.0 .or. ky.ge.twod) call errquit('newfmm:y?',0, FMM_ERR)
      if (kz.lt.0 .or. kz.ge.twod) call errquit('newfmm:z?',0, FMM_ERR)
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
      if (kx.lt.0 .or. kx.ge.twod) call errquit('newfmm:x?',0, FMM_ERR)
      if (ky.lt.0 .or. ky.ge.twod) call errquit('newfmm:y?',0, FMM_ERR)
      if (kz.lt.0 .or. kz.ge.twod) call errquit('newfmm:z?',0, FMM_ERR)
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
      if (kx.lt.0 .or. kx.ge.twod) call errquit('newfmm:x?',0, FMM_ERR)
      if (ky.lt.0 .or. ky.ge.twod) call errquit('newfmm:y?',0, FMM_ERR)
      if (kz.lt.0 .or. kz.ge.twod) call errquit('newfmm:z?',0, FMM_ERR)
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
      subroutine fmm_potential_at_charges2(depth,lmax, tree,
     $     dimx,dimy,dimz,
     $     coords, charges, ncharge, list_info, pot)
      implicit none
#include "errquit.fh"
#include "fmmP.fh"
      integer lmax, depth
      double precision dimx, dimy, dimz
      integer ncharge, tree(300)
      double precision charges(ncharge), coords(3,ncharge), pot(ncharge)
      integer list_info(2,0:2**depth-1,0:2**depth-1,0:2**depth-1)
c     
c     Evaluate the near and far field potentials at all charges
c     returning the result in pot
c     
      integer kx, ky, kz, i, l, m, kxn, kyn, kzn
      integer twod, j, dkyhi, dkxhi, luse
      integer kxnlo, kynlo, kznlo
      double precision hx, hy, hz, xo, yo, zo, sum, 
     $     x, y, z, xx, yy, zz, rdist

      integer lq, kxp, kyp, kzp, level
      parameter (lq = FMM_LMAX)
      double precision w(-lq:lq,0:lq),q(-lq:lq,0:lq)
      double precision xlm_local_potential
      external xlm_local_potential
c     
      if (lmax .gt. lq) call errquit('fmm_pac2: lmax>lq', lmax, FMM_ERR)
      twod = 2**depth
c     
      do i = 1, ncharge
         pot(i) = 0.0d0
      enddo
c     
      do kz = 0, twod-1
         do ky = 0, twod-1
            do kx = 0, twod-1
c     
c     Compute the far field term if there is one
c     
               if (list_info(1,kx,ky,kz).le.list_info(2,kx,ky,kz)) then
                  kxp = kx
                  kyp = ky
                  kzp = kz
                  hx = dimx/dble(twod) ! Size of boxes
                  hy = dimy/dble(twod)
                  hz = dimz/dble(twod)
                  do level = depth, 2, -1
                     xo = (dble(kxp)+0.5d0)*hx ! Coords of box center
                     yo = (dble(kyp)+0.5d0)*hy
                     zo = (dble(kzp)+0.5d0)*hz
                     call fmm_tree_get_cell(tree, level, kxp, kyp, kzp, 
     $                    q, lq, luse)
                     if (luse .ge. 0) goto 11
                     kxp = kxp /2
                     kyp = kyp / 2
                     kzp = kzp / 2
                     hx = hx + hx
                     hy = hy + hy
                     hz = hz + hz
                  enddo
 11               continue
                  if (luse .ge. 0) then
                     do i = list_info(1,kx,ky,kz),list_info(2,kx,ky,kz)
                        x = coords(1,i)
                        y = coords(2,i)
                        z = coords(3,i)
                        sum = 0.0d0
                        call xlm(lmax, x-xo, y-yo, z-zo, w, lq)
                        do l = 0, lmax
                           do m = -l,l
                              sum = sum + w(m,l)*q(m,l)
                           end do
                        end do
                        pot(i) = pot(i) + sum
                     enddo
                  endif
               endif
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
      subroutine fmm_potential_at_charges(depth,lmax,tree,
     $     dimx,dimy,dimz,
     $     coords, charges, ncharge, list_info, pot)
      implicit none
#include "errquit.fh"
#include "fmmP.fh"
      integer lmax, depth, tree(300)
      double precision dimx, dimy, dimz
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

      integer lq, luse
      parameter (lq = FMM_LMAX)
      double precision w(-lq:lq,0:lq), q(-lq:lq,0:lq)
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
         call fmm_tree_get_cell(tree, depth, kx, ky, kz, 
     $        q, lq, luse)
         if (luse .ge. 0) then
            call xlm(lmax, x-xo, y-yo, z-zo, w, lq)
            ind = 1
            do l = 0, lmax
               do m = -l,l
                  sum = sum + w(m,l)*q(m,l)
               end do
            end do
         endif
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
      subroutine fmm(depth, lmax, tree, dimx, dimy, dimz)
      implicit none
#include "fmmP.fh"
      double precision flops, flopsv, flopst, calls
      common /xlmflops/flops, flopsv, flopst, calls
c
c     Given an oct-tree of unnormalized, real multipoles, sum the 
c     far-field interactions using the FMM.
c
c     The FMM results in a local spherical Taylor series associated
c     with each cell at the finest level from which the far-field 
c     potential may be computed by calling xlm_local_potenial.
c     
c     The direct interaction of each cell with its immediate 
c     neighbours and the interactions within a cell are NOT computed 
c     here and must be formed by other means.
c
c     Arguments
c
c     lmax  = Maximum order of multipole to use.
c
c     depth = level of the tree at which multipoles are provided
c     .       (depth=0 is the box that covers the whole domain)
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
c     tree = the data structure returned by fmm_tree_create()
c     .      and containing multipole sources
c
      integer lmax, depth
      double precision dimx, dimy, dimz
      integer tree(300)
c
      integer level, wtree(300)
c
      do level = 0, 9
         leveloff(level) = (8**level-1)/7
      enddo
c
      flops = 0
      flopsv= 0
      flopst= 0
      calls = 0
c
      call xlm_init
c
      call fmm_tree_create(depth, lmax, wtree)
c
c     Upward pass
c
      call fmm_pass1(depth, lmax, dimx, dimy, dimz, tree, wtree)
c
c     Downward pass
      call fmm_pass2(depth, lmax, dimx, dimy, dimz, tree, wtree)
c
      call fmm_tree_destroy(wtree)
c
*      write(6,*) ' flops calls ', flops,flopsv,flopst,calls
c
      end
      subroutine fmm_pass1(depth, lmax, dimx, dimy, dimz, tree, wtree)
      implicit none
#include "errquit.fh"
#include "fmmP.fh"
#include "mafdecls.fh"
      integer lmax, depth, level
      integer tree(300)
      double precision dimx, dimy, dimz
c     
      double precision hx, hy, hz, x, y, z, xn, yn, zn, xp, yp, zp
      double precision xx, yy, zz
      integer dkx, dky, dkz, l, ln
      integer kx, ky, kz, kxn, kyn, kzn, kxp, kyp, kzp
      integer lq
      parameter (lq = FMM_LMAX)
      double precision q(-lq:lq,0:lq), qn(-lq:lq,0:lq)
      double precision p(-lq:lq,0:lq), pn(-lq:lq,0:lq)
      double precision scale, scale0
      integer lmaxeff, lmaxeffn, hi, lo, twol, ind, indn
      integer wtree(300)
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
c     
      do level = depth, 3, -1
         twol = 2**level
         hx = dimx/dble(twol)
         hy = dimy/dble(twol)
         hz = dimz/dble(twol)
         do kz = 0, twol-1
            do ky = 0, twol-1
               do kx = 0, twol-1
                  x = (dble(kx)+0.5d0)*hx ! Coords of current box center
                  y = (dble(ky)+0.5d0)*hy
                  z = (dble(kz)+0.5d0)*hz
c     
                  call fmm_tree_get_cell(tree,level,kx,ky,kz,q,lmax,l)
                  if (l .ge. 0) then
                     kxp = kx/2 ! Parent box indices
                     kyp = ky/2
                     kzp = kz/2
                     xp = (dble(kxp)+0.5d0)*hx*2.0d0 ! Coords of parent box
                     yp = (dble(kyp)+0.5d0)*hy*2.0d0
                     zp = (dble(kzp)+0.5d0)*hz*2.0d0
c     
c     Pass multipoles up to parent (note that this destroys q so
c     must unpack again afterwards)
c     
                     call xlm_translate2(lmax,xp-x,yp-y,zp-z,q,lmax)
                     call fmm_tree_add_cell(tree,level-1,
     $                    kxp,kyp,kzp,q,lmax,lmax)

*                     if (.not. ma_verify_allocator_stuff()) 
*     $                    call errquit('a',0)
                  endif
               enddo
            enddo
         enddo
      enddo
c
      scale0 = (well_separated+1d0)/0.87d0 - 1.0d0
c     
      do level = depth, 2, -1
         twol = 2**level
         hx = dimx/dble(twol)
         hy = dimy/dble(twol)
         hz = dimz/dble(twol)
         do kz = 0, twol-1
            do ky = 0, twol-1
               do kx = 0, twol-1
                  x = (dble(kx)+0.5d0)*hx ! Coords of current box center
                  y = (dble(ky)+0.5d0)*hy
                  z = (dble(kz)+0.5d0)*hz
c     
                  call fmm_tree_get_cell(tree,level,kx,ky,kz,q,lmax,l)
c     
c     Propagate local taylor series expansion
c     
                  ind = kx+twol*(ky+twol*kz)
                  lo = 2*well_separated
                  hi = 2*well_separated + 1
                  do kzn = max(0,int(kz/2)*2-lo),
     $                 min(twol-1,int(kz/2)*2+hi)
                     dkz = abs(kzn - kz)
                     do kyn = max(0,int(ky/2)*2-lo),
     $                    min(twol-1,int(ky/2)*2+hi)
                        dky = abs(kyn - ky)
                        do kxn = max(0,int(kx/2)*2-lo),
     $                       min(twol-1,int(kx/2)*2+hi)
                           dkx = abs(kxn - kx)
                           indn = kxn+twol*(kyn+twol*kzn)
                           if ((dkx.gt.well_separated .or.
     $                          dky.gt.well_separated .or.
     $                          dkz.gt.well_separated) .and.
     $                          indn.gt.ind) then
c     
c     Limit precision if possible for distant boxes.  h*sqrt(3)/2 is the
c     distance of the corner of a cube from its center.
c     
                              call fmm_tree_get_cell(tree,level,
     $                             kxn,kyn,kzn,qn,lmax,ln)
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
                              endif
                              if (ln .ge. 0) then 
                                 lmaxeffn = 
     $                             min(ln,max(6,nint(lmax*scale)))
                                 call xlm_new_multipole_to_local(
     $                                lmaxeffn, x, y, z, p, lmax,
     $                                xn, yn, zn, qn, lmax)
                                 call fmm_tree_add_cell(wtree, level,
     $                                kx, ky, kz, p, lmax, lmaxeffn)
*                     if (.not. ma_verify_allocator_stuff()) 
*     $                    call errquit('b',0)
                              endif
                              if (l .ge. 0) then
                                 lmaxeff  = 
     $                                min(l, max(6,nint(lmax*scale)))
                                 call xlm_new_multipole_to_local(
     $                                lmaxeff, xn, yn, zn, pn, lmax,
     $                                x, y, z, q, lmax)
                                 call fmm_tree_add_cell(wtree, level,
     $                                kxn, kyn, kzn, pn, lmax, lmaxeff)
*                     if (.not. ma_verify_allocator_stuff()) 
*     $                    call errquit('c',0)
                              endif
                           end if
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do
c
      end
      subroutine fmm_pass2(depth, lmax, dimx, dimy, dimz, tree, wtree)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "fmmP.fh"
      integer lmax, depth, level
      double precision dimx, dimy, dimz
      integer tree(300), wtree(300)
c
      double precision hx, hy, hz, x, y, z, xp, yp, zp
      integer kx, ky, kz, kxp, kyp, kzp, ind, indp
      integer lq, luse, k_probe, l_probe, ncells
      logical odoit
      parameter (lq = FMM_LMAX)
      double precision q((2*lq+1)*(lq+1))
      logical fmm_tree_probe_cell
      external fmm_tree_probe_cell
      integer cellindex
      cellindex(level,kx,ky,kz) = k_probe + leveloff(level) + 
     $     kx + (2**level)*(ky + (2**level)*kz)
c
      if (lq .lt. lmax) call errquit('fmm_pass2: lq ', lmax, FMM_ERR)
c
      ncells = (8**(depth)-1)/7 ! No. of cells excluding lowest level
      if (.not. ma_alloc_get(mt_log,ncells,'fmm:probe',
     $     l_probe,k_probe)) call errquit
     $     ('fmm_pass2: failed to allocate probe tree',ncells, MA_ERR)
      call ifill(ncells,.false.,log_mb(k_probe),1)
c
c     Determine which cells have descendents that have potentials
c     by passing a boolean up the tree.
c     
      do level = depth,3,-1
         do kz = 0, 2**level-1
            do ky = 0, 2**level-1
               do kx = 0, 2**level-1
                  kxp = kx/2    ! Parent box indices
                  kyp = ky/2
                  kzp = kz/2
                  indp = cellindex(level-1,kxp,kyp,kzp)
                  if (level.ne.depth) then
                     ind  = cellindex(level,kx,ky,kz)
                     log_mb(indp) = log_mb(indp) .or. log_mb(ind)
                  endif
                  log_mb(indp) = log_mb(indp) .or.
     $                 fmm_tree_probe_cell(wtree,level,kx,ky,kz)
               enddo
            enddo
         enddo
      enddo
c
      do level = 3, depth
         hx = dimx/dble(2**level)
         hy = dimy/dble(2**level)
         hz = dimz/dble(2**level)
         do kz = 0, 2**level-1
            do ky = 0, 2**level-1
               do kx = 0, 2**level-1
                  x = (dble(kx)+0.5d0)*hx ! Coords of current box center
                  y = (dble(ky)+0.5d0)*hy
                  z = (dble(kz)+0.5d0)*hz
                  kxp = kx/2    ! Parent box indices
                  kyp = ky/2
                  kzp = kz/2
                  xp = (dble(kxp)+0.5d0)*hx*2.0d0 ! Coords of parent box
                  yp = (dble(kyp)+0.5d0)*hy*2.0d0
                  zp = (dble(kzp)+0.5d0)*hz*2.0d0
c
c     Only move the parent's potential down if this cell or any of
c     its decendents already has a non-zero one.
c     
                  odoit = fmm_tree_probe_cell(wtree,level,kx,ky,kz)
                  if (level.ne.depth) odoit = odoit .or. 
     $                 log_mb(cellindex(level,kx,ky,kz))
                  if (odoit) then
                     call fmm_tree_get_cell(wtree,level-1,
     $                    kxp,kyp,kzp,q,lmax,luse)
                     if (luse .ge. 0) then
                        call xlm_local_translate(lmax,x-xp,y-yp,z-zp,
     $                       q,lmax)
                        call fmm_tree_add_cell(wtree,level,kx,ky,kz,
     $                       q,lmax,lmax)
*     if (.not. ma_verify_allocator_stuff()) 
*     $                    call errquit('d',0)
                     endif
                  endif
               end do
            end do
         end do
      enddo
c
      do level = 0, depth
         do kz = 0, 2**level-1
            do ky = 0, 2**level-1
               do kx = 0, 2**level-1
                  call fmm_tree_get_cell(wtree,level,kx,ky,kz,q,lmax,
     $                 luse)
                  if (luse .ge. 0) then
                     call fmm_tree_put_cell(tree,level,kx,ky,kz,q,lmax,
     $                    lmax)
*                     if (.not. ma_verify_allocator_stuff()) 
*     $                    call errquit('e',0)
                  else
                     call fmm_tree_zero_cell(tree,level,kx,ky,kz)
*                     if (.not. ma_verify_allocator_stuff()) 
*     $                    call errquit('f',0)
                  end if
               end do
            end do
         end do
      enddo
c
      if (.not. ma_free_heap(l_probe)) call errquit
     $     ('fmm_pass2: free heap failed ', 0, MA_ERR)
c     
      end
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


#define I_DEPTH 1
#define I_LMAX  2
#define I_NCPA  3
#define I_CURA  4
#define I_NCIA  5
#define I_KTREE 6
#define I_LTREE 7
#define I_KBUFS 100
#define I_LBUFS 200


      subroutine fmm_tree_create(depth, lmax, tree)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "fmmP.fh"
      integer depth, lmax, tree(300)
c     
      integer n_cells, n_cells_per_alloc, l_tree, k_tree, level
c
      n_cells = (8**(depth+1)-1)/7
      n_cells_per_alloc = max(512,n_cells / 100 + 1)
      if (.not. ma_alloc_get(mt_int, n_cells, 'fmm tree ptrs',
     $     l_tree, k_tree)) call errquit
     $     ('fmm_tree_create: unable to allocate ptrs ', n_cells,
     &       MA_ERR)
      call ifill(n_cells, -1, int_mb(k_tree), 1)
c      
      tree(I_DEPTH) = depth
      tree(I_LMAX)  = lmax
      tree(I_NCPA)  = n_cells_per_alloc
      tree(I_CURA)  = -1
      tree(I_NCIA)  = n_cells_per_alloc ! Start full to force alloc
      tree(I_KTREE) = k_tree
      tree(I_LTREE) = l_tree
c
      do level = 0, 9
         leveloff(level) = (8**level-1)/7
      enddo
c
      end
      subroutine fmm_tree_destroy(tree)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
      integer tree(300)
c
      integer i
c
      do i = 0, tree(I_CURA)
         if (.not. ma_free_heap(tree(I_LBUFS+I))) call errquit
     $        ('fmm_tree_destroy: failed freeing cells',i, MA_ERR)
      enddo
      if (.not. ma_free_heap(tree(I_LTREE))) call errquit
     $     ('fmm_tree_destroy: failed freeing tree',0, MA_ERR)
c
      end
      subroutine fmm_tree_put_cell(tree, level, kx, ky, kz, 
     $     q, lq, luse)
      implicit none
#include "errquit.fh"
#include "fmmP.fh"
#include "mafdecls.fh"
      integer tree(300), level, kx, ky, kz, lq, luse
      double precision q(-lq:lq,0:*)
c
      integer lmax, ind, ptr, n_cells_per_alloc, curalloc, 
     $     n_cells_in_alloc, l, m
c
      integer cellindex
      cellindex(level,kx,ky,kz) = tree(I_KTREE) + leveloff(level) + 
     $     kx + (2**level)*(ky + (2**level)*kz)
c
      lmax = tree(2)
      ind = cellindex(level, kx, ky, kz) ! Index of desired cell
      ptr = int_mb(ind)         ! Pointer to cell poles
c
*      write(6,*) ' Putting    ', level, kx, ky, kz, 
*     $     cellindex(level, kx, ky, kz)-tree(I_KTREE), ptr
c
      if (ptr .eq. -1) then
         n_cells_per_alloc = tree(I_NCPA)
         curalloc = tree(I_CURA)
         n_cells_in_alloc = tree(I_NCIA)
         if (n_cells_in_alloc .ge. n_cells_per_alloc) then
            curalloc = curalloc + 1
            n_cells_in_alloc = 0
            if (.not. ma_alloc_get(mt_dbl, 
     $           (lmax+1)**2*n_cells_per_alloc, 'fmm tree buf',
     $           tree(I_LBUFS+curalloc), tree(I_KBUFS+curalloc)))
     $           call errquit('fmm_tree_put_cell: ma failed ',
     $           (lmax+1)**2*n_cells_per_alloc, MA_ERR)
            call dfill((lmax+1)**2*n_cells_per_alloc, 0.0d0,
     $           dbl_mb(tree(I_KBUFS+curalloc)), 1)
         endif
         ptr = tree(I_KBUFS+curalloc) + n_cells_in_alloc*(lmax+1)**2

*         write(6,*) ' Allocating ', level, kx, ky, kz, 
*     $        cellindex(level, kx, ky, kz)-tree(I_KTREE), ptr

         n_cells_in_alloc = n_cells_in_alloc + 1
c
         tree(I_CURA) = curalloc
         tree(I_NCIA) = n_cells_in_alloc
         int_mb(ind)   = ptr
      endif
c
      do l = 0, luse
         do m = -l, l
            dbl_mb(ptr) = q(m,l)
            ptr = ptr + 1
         enddo
      enddo
      do l = luse+1,lmax
         do m = -l, l
            dbl_mb(ptr) = 0.0d0
            ptr = ptr + 1
         enddo
      enddo
c
      end
      subroutine fmm_tree_get_cell(tree, level, kx, ky, kz, 
     $     q, lq, luse)
      implicit none
#include "fmmP.fh"
#include "mafdecls.fh"
      integer tree(300), level, kx, ky, kz, lq, luse
      double precision q(-lq:lq,0:*)
c
      integer lmax, l, m, ptr
c
      integer cellindex
      cellindex(level,kx,ky,kz) = tree(I_KTREE) + leveloff(level) + 
     $     kx + (2**level)*(ky + (2**level)*kz)
c
      lmax = tree(I_LMAX)
      ptr = int_mb(cellindex(level,kx,ky,kz))
c
      if (ptr .ne. -1) then
         do l = 0, lmax
            do m = -l, l
               q(m,l) = dbl_mb(ptr)
               ptr = ptr + 1
            enddo
         enddo
         do l = lmax, 0, -1
            do m = -l, l
               if (abs(q(m,l)).ne.0.0d0) goto 10
            enddo
         enddo
 10      luse = l
      else
         luse = -1
      endif
c
      end
      subroutine fmm_tree_add_cell(tree, level, kx, ky, kz, 
     $     q, lq, luse)
      implicit none
#include "mafdecls.fh"
#include "fmmP.fh"
      integer tree(300), level, kx, ky, kz, lq, luse
      double precision q(-lq:lq,0:*)
c
      integer lmax, ptr, l, m
c
      integer cellindex
      cellindex(level,kx,ky,kz) = tree(I_KTREE) + leveloff(level) + 
     $     kx + (2**level)*(ky + (2**level)*kz)
c
      lmax = tree(2)
      ptr = int_mb(cellindex(level, kx, ky, kz))
c
*      write(6,*) ' Adding ', level, kx, ky, kz, 
*     $     cellindex(level, kx, ky, kz)-tree(I_KTREE), ptr
      if (ptr .eq. -1) then
         call fmm_tree_put_cell(tree, level, kx, ky, kz, 
     $        q, lq, luse)
      else
         do l = 0, luse
            do m = -l, l
               dbl_mb(ptr) = dbl_mb(ptr) + q(m,l)
               ptr = ptr + 1
            enddo
         enddo
      endif
c
      end
      subroutine fmm_tree_zero_cell(tree, level, kx, ky, kz)
      implicit none
#include "mafdecls.fh"
#include "fmmP.fh"
      integer tree(300), level, kx, ky, kz
c
      integer lmax, ptr, l, m
c
      integer cellindex
      cellindex(level,kx,ky,kz) = tree(I_KTREE) + leveloff(level) + 
     $     kx + (2**level)*(ky + (2**level)*kz)
c
      lmax = tree(2)
      ptr = int_mb(cellindex(level, kx, ky, kz))
      if (ptr .ne. -1) then
         do l = 0, lmax
            do m = -l, l
               dbl_mb(ptr) = 0.0d0
               ptr = ptr + 1
            enddo
         enddo
      endif
c
      end
      logical function fmm_tree_probe_cell(tree, level, kx, ky, kz)
      implicit none
#include "mafdecls.fh"
#include "fmmP.fh"
      integer tree(300), level, kx, ky, kz
c
      integer lmax, ptr
c
      integer cellindex
      cellindex(level,kx,ky,kz) = tree(I_KTREE) + leveloff(level) + 
     $     kx + (2**level)*(ky + (2**level)*kz)
c
      lmax = tree(2)
      ptr = int_mb(cellindex(level, kx, ky, kz))
      fmm_tree_probe_cell = ptr .ne. -1
c
      end
      subroutine fmm_tree_summarize(tree)
      implicit none
#include "mafdecls.fh"
#include "fmmP.fh"
      integer tree(300), level, kx, ky, kz, depth, lmax
      integer ncell, curalloc, ncell_per_alloc
      logical fmm_tree_probe_cell
      external fmm_tree_probe_cell
c
      depth = tree(I_DEPTH)
      lmax  = tree(I_LMAX)
      curalloc = tree(I_CURA)
      ncell_per_alloc = tree(I_NCPA)
      write(6,1) depth, lmax
 1    format(/' Summary of FMM tree:  depth=',i1,'  lmax=',i2)
      do level = 0, depth
         ncell = 0
         do kz = 0, 2**level-1
            do ky = 0, 2**level-1
               do kx = 0, 2**level-1
                  if (fmm_tree_probe_cell(tree,level,kx,ky,kz))
     $                 ncell = ncell + 1
               enddo
            enddo
         enddo
         write(6,2) level, ncell, (100.0d0*ncell)/8d0**level
 2       format(' Level',i2,': ncells=',i7,' (',f6.2,'%)')
      enddo
      write(6,3) dble((curalloc+1)*ncell_per_alloc)*((lmax+1)**2)*8d-6 +
     $     ((8d0**(depth+1)-1)/7d0)*4d-6
 3    format(' Allocated memory ', f6.2,' Mbytes')
c
      end
        
