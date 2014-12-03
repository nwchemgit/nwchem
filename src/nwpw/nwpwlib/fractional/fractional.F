c     $Id$

      subroutine fractional_to_cartesian(rtdb)
      implicit none
#include "errquit.fh"
      integer rtdb

#include "global.fh"
#include "inp.fh"
#include "bafdecls.fh"
#include "btdb.fh"
#include "geom.fh"


*     **** local variables ****
      integer geom
      integer i,j,l
      logical value
      double precision rxyz(3),q
      character*16     t

      character*50 cell_name
      character*50 rtdb_unita

      integer nion,rion(2)      
      real*8  tion(3)
      real*8 unita(3,3)
      real*8 a(3,3)

      if (.not. inp_a(cell_name))
     >     call errquit(
     >          'fractional_to_cartesian: failed to read keyword', 0,
     &       INPUT_ERR)
      
*     *************************
*     **** read in lattice ****
*     *************************

*     ***** cell: stuff ****
      l = index(cell_name,' ') - 1
      rtdb_unita = cell_name(1:l)//':unita'
           
      value = btdb_get(rtdb,rtdb_unita,mt_dbl,9,unita)


*     ********************************
*     **** read in ion structure *****
*     ********************************
      value = geom_create(geom,'geometry')
      value = geom_rtdb_load(rtdb,geom,'geometry')

*     **** get the number ions ****
      value = geom_ncent(geom,nion)


*     ***** allocate ion data structure *****
      value = BA_push_get(mt_dbl,(3*nion),'rion',rion(2),rion(1))
      call dcopy(3*nion,0.0d0,0,dbl_mb(rion(1)),1)
 
*     **** get ion positions in fractional coordinates ****
      do i=1,nion
         value = geom_cent_get(geom,i,t,dbl_mb(rion(1)+(i-1)*3),q)
      end do


      write(*,*)
      write(*,*) "**** fractional --> cartesian ****"
      write(*,*)
      write(*,*) "Lattice Vectors"
      write(*,*) "---------------"
      write(*,114) "  unita  :",unita(1,1),unita(2,1),unita(3,1)
      write(*,114) "          ",unita(1,2),unita(2,2),unita(3,2)
      write(*,114) "          ",unita(1,3),unita(2,3),unita(3,3)
  114 format(a,3F10.3) 



*     ********************************************
*     ***** convert to cartesian coordinates *****
*     ********************************************
      do j=1,3
      do i=1,3
        a(i,j) = unita(i,j)
      end do
      end do

      do i=1,nion
         tion(1) = dbl_mb(rion(1)+0+(i-1)*3)
         tion(2) = dbl_mb(rion(1)+1+(i-1)*3)
         tion(3) = dbl_mb(rion(1)+2+(i-1)*3)

         dbl_mb(rion(1)+0+(i-1)*3) = a(1,1)*tion(1)
     >                             + a(1,2)*tion(2)
     >                             + a(1,3)*tion(3)
         dbl_mb(rion(1)+1+(i-1)*3) = a(2,1)*tion(1)
     >                             + a(2,2)*tion(2)
     >                             + a(2,3)*tion(3)
         dbl_mb(rion(1)+2+(i-1)*3) = a(3,1)*tion(1)
     >                             + a(3,2)*tion(2)
     >                             + a(3,3)*tion(3)
      end do



*     ************************************************
*     **** write out converted ion data structure ****
*     ************************************************
      do i=1,nion
         value = geom_cent_get(geom,i,t,rxyz,q)
         value = geom_cent_set(geom,i,t,dbl_mb(rion(1)+(i-1)*3),q)
      end do
      value = geom_rtdb_delete(rtdb,'geometry')
      value = geom_rtdb_store(rtdb,geom,'geometry')
      value = geom_destroy(geom)


*     ***************************************************
*     **** rewrite geometry data structure to stdout ****
*     ***************************************************
      if (.not. geom_create(geom,'geometry'))
     $        call errquit('geom_create:  failed ', 0, GEOM_ERR)
      if (.not. geom_rtdb_load(rtdb,geom,'geometry'))
     $        call errquit('geom_rtdb_load: ', 0, RTDB_ERR)
      if (.not. geom_print(geom))
     $        call errquit('geom_input: print failed ', 0, GEOM_ERR)
      if (.not.geom_print_distances(geom))
     &        call errquit('geom_input: print_distances failed ',911,
     &       GEOM_ERR)
      if (.not.geom_print_angles(geom))
     &        call errquit('geom_input: print_angles failed ',911,
     &       GEOM_ERR)
      if(.not. geom_verify_coords(geom))
     &     call errquit
     &     ('geom_input: geom_verify_coords failed',911, GEOM_ERR)
      if (.not. geom_destroy(geom))
     $     call errquit('geom_input: geom_destroy failed', 0, GEOM_ERR)

      value = BA_pop_stack(rion(2))
      return
      END




      subroutine cartesian_to_fractional(rtdb)
      implicit none
#include "errquit.fh"
      integer rtdb

#include "global.fh"
#include "inp.fh"
#include "bafdecls.fh"
#include "btdb.fh"
#include "geom.fh"


*     **** local variables ****
      integer geom
      integer i,l
      logical value
      double precision rxyz(3),q
      character*16     t

      character*50 cell_name
      character*50 rtdb_unita

      integer nion,rion(2)      
      real*8  tion(3)
      real*8 unita(3,3)
      real*8 b(3,3),volume

      if (.not. inp_a(cell_name))
     >     call errquit(
     >          'fractional_to_cartesian: failed to read keyword', 0,
     &       INPUT_ERR)
      
*     *************************
*     **** read in lattice ****
*     *************************

*     ***** cell: stuff ****
      l = index(cell_name,' ') - 1
      rtdb_unita = cell_name(1:l)//':unita'
           
      value = btdb_get(rtdb,rtdb_unita,mt_dbl,9,unita)
      

*     ********************************
*     **** read in ion structure *****
*     ********************************
      value = geom_create(geom,'geometry')
      value = geom_rtdb_load(rtdb,geom,'geometry')

*     **** get the number ions ****
      value = geom_ncent(geom,nion)


*     ***** allocate ion data structure *****
      value = BA_push_get(mt_dbl,(3*nion),'rion',rion(2),rion(1))
      call dcopy(3*nion,0.0d0,0,dbl_mb(rion(1)),1)
 
*     **** get ion positions in fractional coordinates ****
      do i=1,nion
         value = geom_cent_get(geom,i,t,dbl_mb(rion(1)+(i-1)*3),q)
      end do
 
      write(*,*)
      write(*,*) "**** cartesian --> fractional ****"
      write(*,*)
      write(*,*) "Lattice Vectors"
      write(*,*) "---------------"
      write(*,114) "  unita  :",unita(1,1),unita(2,1),unita(3,1)
      write(*,114) "          ",unita(1,2),unita(2,2),unita(3,2)
      write(*,114) "          ",unita(1,3),unita(2,3),unita(3,3)
  114 format(a,3F10.3) 

 
*     ********************************************
*     ***** convert to fractional coordinates *****
*     ********************************************
      
      b(1,1) = unita(2,2)*unita(3,3) - unita(3,2)*unita(2,3)
      b(2,1) = unita(3,2)*unita(1,3) - unita(1,2)*unita(3,3)
      b(3,1) = unita(1,2)*unita(2,3) - unita(2,2)*unita(1,3)
      b(1,2) = unita(2,3)*unita(3,1) - unita(3,3)*unita(2,1)
      b(2,2) = unita(3,3)*unita(1,1) - unita(1,3)*unita(3,1)
      b(3,2) = unita(1,3)*unita(2,1) - unita(2,3)*unita(1,1)
      b(1,3) = unita(2,1)*unita(3,2) - unita(3,1)*unita(2,2)
      b(2,3) = unita(3,1)*unita(1,2) - unita(1,1)*unita(3,2)
      b(3,3) = unita(1,1)*unita(2,2) - unita(2,1)*unita(1,2)
      volume = unita(1,1)*b(1,1)
     >       + unita(2,1)*b(2,1)
     >       + unita(3,1)*b(3,1)
      
      volume = 1.0d0/volume
      call dscal(9,volume,b,1)


      do i=1,nion
         tion(1) = dbl_mb(rion(1)+0+(i-1)*3)
         tion(2) = dbl_mb(rion(1)+1+(i-1)*3)
         tion(3) = dbl_mb(rion(1)+2+(i-1)*3)

         dbl_mb(rion(1)+0+(i-1)*3) = b(1,1)*tion(1)
     >                             + b(1,2)*tion(2)
     >                             + b(1,3)*tion(3)
         dbl_mb(rion(1)+1+(i-1)*3) = b(2,1)*tion(1)
     >                             + b(2,2)*tion(2)
     >                             + b(2,3)*tion(3)
         dbl_mb(rion(1)+2+(i-1)*3) = b(3,1)*tion(1)
     >                             + b(3,2)*tion(2)
     >                             + b(3,3)*tion(3)
      end do




*     ************************************************
*     **** write out converted ion data structure ****
*     ************************************************
      do i=1,nion
         value = geom_cent_get(geom,i,t,rxyz,q)
         value = geom_cent_set(geom,i,t,dbl_mb(rion(1)+(i-1)*3),q)
      end do
      value = geom_rtdb_delete(rtdb,'geometry')
      value = geom_rtdb_store(rtdb,geom,'geometry')
      value = geom_destroy(geom)



*     ***************************************************
*     **** rewrite geometry data structure to stdout ****
*     ***************************************************
      if (.not. geom_create(geom,'geometry'))
     $        call errquit('geom_create:  failed ', 0, GEOM_ERR)
      if (.not. geom_rtdb_load(rtdb,geom,'geometry'))
     $        call errquit('geom_rtdb_load: ', 0, RTDB_ERR)
      if (.not. geom_print(geom))
     $        call errquit('geom_input: print failed ', 0, GEOM_ERR)
      if (.not.geom_print_distances(geom))
     &        call errquit('geom_input: print_distances failed ',911,
     &       GEOM_ERR)
      if (.not.geom_print_angles(geom))
     &        call errquit('geom_input: print_angles failed ',911,
     &       GEOM_ERR)
      if(.not. geom_verify_coords(geom))
     &     call errquit
     &     ('geom_input: geom_verify_coords failed',911, GEOM_ERR)
      if (.not. geom_destroy(geom))
     $     call errquit('geom_input: geom_destroy failed', 0, GEOM_ERR)
      
      value = BA_pop_stack(rion(2))

      return
      END




