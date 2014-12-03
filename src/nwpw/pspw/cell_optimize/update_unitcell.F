*
* $Id$
*
      logical function update_unitcell(rtdb)
      implicit none
      integer rtdb
 
#include "bafdecls.fh"
#include "btdb.fh"
#include "geom.fh"

*     **** local variables ****
      logical value
      integer l,geom
      integer r1(2),nion,nionplus3
      double precision alpha
      double precision unita(3,3)
      double precision unita2(3,3)
      double precision stress(3,3)
      character*50 cell_name,rtdb_unita

      value = .true.

*     **** initialize and read in geometry object ****
      value     = value.and.geom_create(geom,'geometry')
      value     = value.and.geom_rtdb_load(rtdb,geom,'geometry')
      value     = value.and.geom_ncent(geom,nion)
      nionplus3 = nion+3
      value     = value.and.
     >            BA_alloc_get(mt_dbl,(3*nion),'r1',r1(2),r1(1))
      call load_cellgeometry(geom,dbl_mb(r1(1)),unita)

*     **** read in stress ****
      value = value.and.btdb_get(rtdb,'pspw:stress',mt_dbl,9,stress)


      write(*,*) " -------------------- "
      write(*,*) "| Updating unit cell |"
      write(*,*) " -------------------- "

*     *** simple steepest descent update ****
      alpha = -0.50d0
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


*     **** write out updated geometry object ****
      call store_cellgeometry(geom,dbl_mb(r1(1)),unita2)
      value = value.and.geom_rtdb_delete(rtdb,'geometry')
      value = value.and.geom_rtdb_store(rtdb,geom,'geometry')
      value = value.and.geom_destroy(geom)
      value = value.and.BA_free_heap(r1(2))
      write(*,*) "geometry adjusted......"

 
      update_unitcell = value
      return 
      end

