      logical function update_unitcell(rtdb)
*
* $Id: update_unitcell.F,v 1.6 2002-05-04 19:28:29 bylaska Exp $
*
      implicit none
#include "mafdecls.fh"
#include "rtdb.fh"
      integer rtdb
      logical value
 

      integer l

      character*30 cell_name,rtdb_unita

      double precision alpha
      double precision unita(3,3)
      double precision unita2(3,3)
      double precision stress(3,3)

*     **** cell: stuff ****
      value = rtdb_get(rtdb,'pspw:stress',mt_dbl,9,stress)
      value = value.and.
     >        rtdb_cget(rtdb,'cgsd:cell_name',1,cell_name)
      l = index(cell_name,' ') - 1
      rtdb_unita = cell_name(1:l)//':unita'
      value = value.and.
     >        rtdb_get(rtdb,rtdb_unita,mt_dbl,9,unita)

*     **** check unita for defaults ****
      call check_unita_for_default(rtdb,unita)



      write(*,*) " -------------------- "
      write(*,*) "| Updating unit cell |"
      write(*,*) " -------------------- "

*     *** simple steepest descent update ****
      alpha = -5.00d0
      call dcopy(9,unita,1,unita2,1)
      call daxpy(9,alpha,stress,1,unita2,1)



      write(*,*) "strain tensor:",stress(1,1),stress(2,1),stress(3,1)
      write(*,*) "              ",stress(1,2),stress(2,2),stress(3,2)
      write(*,*) "              ",stress(1,3),stress(2,3),stress(3,3)
      write(*,*) 
      write(*,*) "old unita:",unita(1,1),unita(2,1),unita(3,1)
      write(*,*) "          ",unita(1,2),unita(2,2),unita(3,2)
      write(*,*) "          ",unita(1,3),unita(2,3),unita(3,3)
      write(*,*) 
      write(*,*) "new unita:",unita2(1,1),unita2(2,1),unita2(3,1)
      write(*,*) "          ",unita2(1,2),unita2(2,2),unita2(3,2)
      write(*,*) "          ",unita2(1,3),unita2(2,3),unita2(3,3)
      write(*,*) 

*     **** update geometry object ****
      write(*,*) "Ion positions adjusted......"
      call adjust_cellgeometry(rtdb,unita,unita2)

      value = value.and.
     >        rtdb_put(rtdb,rtdb_unita,mt_dbl,9,unita2)

 
      update_unitcell = value
      return 
      end


 
      subroutine adjust_cellgeometry(rtdb,unita,unita2)
      implicit none
      integer rtdb
      real*8 unita(3,3),unita2(3,3)

#include "mafdecls.fh"
#include "rtdb.fh"
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
      value = MA_alloc_get(mt_dbl,(3*nion),'r1',r1(2),r1(1))

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

      value = MA_free_heap(r1(2))
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
