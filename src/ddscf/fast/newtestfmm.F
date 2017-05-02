      program test
*
* $Id$
*
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
c
      double precision total_flops
      common / totalf / total_flops
c
      integer lmax, depth, tree(300)
      parameter (lmax=16, depth=4)
      double precision fmm_potential, direct_sum, util_random
      double precision fmm_direct_potential
      external fmm_potential, fmm_direct_potential, direct_sum, 
     $     util_random
      double precision x, y, z, dimx, dimy, dimz, time
      real*4 etime, stuff(2)
c
      integer ncharge, i, ntest, j
      parameter (ncharge = 5000, ntest=1000)
      double precision charges(ncharge), coords(3,ncharge)
      double precision chargesr(ncharge), coordsr(3,ncharge)
      double precision pot1(ncharge), pot2(ncharge), pot3(ncharge)
      integer list_next(ncharge), map(ncharge),
     $     list_info(2,0:2**depth-1,0:2**depth-1,0:2**depth-1)
c
      double precision far, near, err, relerr, l2err, l2norm, sum
c
      if (.not. ma_init(mt_dbl, 60000000, 60000000))
     $     call errquit('nwchem: ma_init failed', -1, MA_ERR)
c
      call xlm_init
      total_flops = 0.0d0
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
      write(6,*)
      write(6,*)
      write(6,*) ' Depth ', depth, ' lmax ', lmax
c
      call fmm_tree_create(depth, lmax, tree)
      write(6,*) ' RIght after create '
      call fmm_tree_summarize(tree)
c
      time = etime(stuff)
      call fmm_charges_to_poles(depth, lmax, tree, dimx, dimy, dimz,
     $     coords, charges, ncharge, list_info, list_next,
     $     map)
      write(6,*) ' RIght after charges to polse '
      call fmm_tree_summarize(tree)
c
      time = etime(stuff) - time
      write(6,61) ' Time to form multipoles ', time
 61   format(1x,a,1x,f8.2)
c     
      time = etime(stuff)
      call fmm(depth, lmax, tree, dimx, dimy, dimz)
      time = etime(stuff) - time
      write(6,61) ' Time to perform far field summation ', time
      write(6,*) ' RIght after FMM '
      call fmm_tree_summarize(tree)
c
      l2err = 0.0d0
      l2norm = 0.0d0
      err = 0.0d0
      relerr = 0.0d0
      do i = 1, ntest
         x = util_random(0)
         y = util_random(0)
         z = util_random(0)
         far = fmm_potential(depth,lmax,tree,dimx,dimy,dimz,
     $        x,y,z)
         near = fmm_direct_potential(
     $        lmax, depth, dimx, dimy, dimz, x, y, z,
     $        coords, charges, ncharge, list_info, map)
         sum = direct_sum(coords, charges, ncharge, x, y, z)
         err = max(err,abs(sum-far-near))
*         write(6,*) ' far, near, sum, err ', far, near, sum, err
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
c$$$      time = etime(stuff)
c$$$      do j = 1, 1000
c$$$         far = fmm_potential(depth,lmax,tree,dimx,dimy,dimz,x,y,z)
c$$$         near = fmm_direct_potential(lmax, depth,
c$$$     $        dimx, dimy, dimz, x, y, z,
c$$$     $        coords, charges, ncharge, list_info, map)
c$$$      enddo
c$$$      time = etime(stuff)-time
c$$$      write(6,61) ' time per FMM potential eval ', time/1000
c$$$      time = etime(stuff)
c$$$      do j = 1, 1000
c$$$         sum = direct_sum(coords, charges, ncharge, x, y, z)
c$$$      enddo
c$$$      time = etime(stuff)-time
c$$$      write(6,61) ' time per direct potential eval ', time/1000
c
      time = etime(stuff)
      do j = 1, ncharge
         chargesr(j) = charges(map(j))
         coordsr(1,j) = coords(1,map(j))
         coordsr(2,j) = coords(2,map(j))
         coordsr(3,j) = coords(3,map(j))
      enddo
      call fmm_potential_at_charges2(depth, lmax, tree, 
     $     dimx,dimy,dimz,
     $     coordsr, chargesr, ncharge, list_info, pot3)
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
      end
