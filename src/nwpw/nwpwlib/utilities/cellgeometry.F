*
* $Id$
*

      subroutine load_cellgeometry(geom,coords,amatrix)
      implicit none
      integer geom
      double precision coords(*)
      double precision amatrix(*)

#include "errquit.fh"
#include "geom.fh"

      if (.not. geom_cart_coords_get(geom, coords))
     >   call errquit('load_cellgeometry: geom cart?',
     >                0,GEOM_ERR)
      if (.not. geom_cart_to_frac(geom, coords))
     >    call errquit('load_cell_geometry: cart_to_frac?',
     >                 0,GEOM_ERR)
      if (.not. geom_amatrix_get(geom, amatrix))
     >   call errquit('load_cellgeometry: failed to get amatrix',
     >                 0,GEOM_ERR)
      return
      end

      subroutine store_cellgeometry(geom,coords,amatrix)
      implicit none
      integer geom
      double precision coords(*)
      double precision amatrix(3,3)

#include "errquit.fh"
#include "geom.fh"

*     **** external functions ****
      logical  geom_amatrix_set
      external geom_amatrix_set

      if (.not. geom_amatrix_set(geom,amatrix))
     >   call errquit('store_cellgeometry: failed to set amatrix',
     >                0,GEOM_ERR)
      if (.not. geom_frac_to_cart(geom,coords))
     >   call errquit('store_cell_geometry: frac_to_cart?',
     >                0,GEOM_ERR)
      if (.not. geom_cart_coords_set(geom,coords))
     >   call errquit('store_cellgeometry: geom cart?',
     >                0,GEOM_ERR)
      if (.not. geom_cart_to_frac(geom, coords))
     >    call errquit('store_cell_geometry: cart_to_frac?',
     >                 0,GEOM_ERR)

      !**** testing ****
      if (.not. geom_amatrix_get(geom, amatrix))
     >   call errquit('load_cellgeometry: failed to get amatrix',
     >                 0,GEOM_ERR)

      return
      end


 
      subroutine adjust_cellgeometry(rtdb,unita,unita2)
      implicit none
      integer rtdb
      real*8 unita(3,3),unita2(3,3)

#include "bafdecls.fh"
#include "btdb.fh"
#include "geom.fh"

      logical value
      integer i,nion,r1(2),geom
      real*8 rxyz(3),q
      character*16     t

      
      value = geom_create(geom,'geometry')
      value = geom_rtdb_load(rtdb,geom,'geometry')

*     **** get the number ions ****
      value = geom_ncent(geom,nion)

*     ***** allocate ion data structure *****
      value = BA_alloc_get(mt_dbl,(3*nion),'r1',r1(2),r1(1))

      do i=1,nion
         value = geom_cent_get(geom,i,t,dbl_mb(r1(1)+(i-1)*3),q)
      end do

      call adjust_ions(nion,dbl_mb(r1(1)),unita,unita2)

      do i=1,nion
         value = geom_cent_get(geom,i,t,rxyz,q)
         value = geom_cent_set(geom,i,t,dbl_mb(r1(1)+(i-1)*3),q)
      end do
      value = geom_rtdb_delete(rtdb,'geometry')
      value = geom_rtdb_store(rtdb,geom,'geometry')
      value = geom_destroy(geom)

      value = BA_free_heap(r1(2))
      return
      end


      subroutine adjust_ions(nion,rion,unita,unita2)
      implicit none
      integer nion
      real*8  rion(3,*)
      real*8  unita(3,3)
      real*8  unita2(3,3)

      integer i,j
      real*8 a(3,3),b(3,3),tion(3),volume

*     **** Convert from real coordinates to fractional coordinates. *****
      do j=1,3
      do i=1,3
        a(i,j) = unita(i,j)
      end do
      end do

      b(1,1) = a(2,2)*a(3,3) - a(3,2)*a(2,3)
      b(2,1) = a(3,2)*a(1,3) - a(1,2)*a(3,3)
      b(3,1) = a(1,2)*a(2,3) - a(2,2)*a(1,3)

      b(1,2) = a(2,3)*a(3,1) - a(3,3)*a(2,1)
      b(2,2) = a(3,3)*a(1,1) - a(1,3)*a(3,1)
      b(3,2) = a(1,3)*a(2,1) - a(2,3)*a(1,1)

      b(1,3) = a(2,1)*a(3,2) - a(3,1)*a(2,2)
      b(2,3) = a(3,1)*a(1,2) - a(1,1)*a(3,2)
      b(3,3) = a(1,1)*a(2,2) - a(2,1)*a(1,2)
      volume = a(1,1)*b(1,1)
     >       + a(2,1)*b(2,1)
     >       + a(3,1)*b(3,1)

      volume = 1.0d0/volume
      call dscal(9,volume,b,1)

      do i=1,nion
         tion(1) = rion(1,i)
         tion(2) = rion(2,i)
         tion(3) = rion(3,i)
         rion(1,i) = b(1,1)*tion(1)
     >             + b(2,1)*tion(2)
     >             + b(3,1)*tion(3)
         rion(2,i) = b(1,2)*tion(1)
     >             + b(2,2)*tion(2)
     >             + b(3,2)*tion(3)
         rion(3,i) = b(1,3)*tion(1)
     >             + b(2,3)*tion(2)
     >             + b(3,3)*tion(3)

      end do

*     **** Convert from fractional coordinates to real coordinates. *****
      do j=1,3
      do i=1,3
        a(i,j) = unita2(i,j)
      end do
      end do


      do i=1,nion
         tion(1) = rion(1,i)
         tion(2) = rion(2,i)
         tion(3) = rion(3,i)

         rion(1,i) = a(1,1)*tion(1)
     >             + a(1,2)*tion(2)
     >             + a(1,3)*tion(3)
         rion(2,i) = a(2,1)*tion(1)
     >             + a(2,2)*tion(2)
     >             + a(2,3)*tion(3)
         rion(3,i) = a(3,1)*tion(1)
     >             + a(3,2)*tion(2)
     >             + a(3,3)*tion(3)
      end do

      return
      end
