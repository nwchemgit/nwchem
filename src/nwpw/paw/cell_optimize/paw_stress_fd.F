*
* $Id$
*
*     *****************************************
*     *                                       *
*     *            paw_stress_fd              *
*     *                                       *
*     *****************************************
      logical function paw_stress_fd(rtdb)
      implicit none
      integer rtdb
      logical value

#include "bafdecls.fh"
#include "btdb.fh"
#include "inp.fh"
#include "util.fh"
#include "geom.fh"
#include "stdio.fh"

      integer MASTER,taskid
      parameter (MASTER=0)
      logical oprint,lstress_numerical,noit_energy

      integer l,i,j
      integer r1(2),geom,nion


      double precision alpha,E
      double precision unita(3,3)
      double precision unita2(3,3)
      double precision unita3(3,3)

      double precision stress(3,3)
      double precision ke_stress(3,3)
      double precision xc_stress(3,3)
      double precision ewald_stress(3,3)
      double precision hartree_stress(3,3)
      double precision local_stress(3,3)
      double precision alpha_stress(3,3)
      double precision RIJ_stress(3,3)
      double precision ewaldG_stress(3,3)

      double precision lattice(6)
      double precision lattice2(6)
      double precision lattice3(6)

      double precision lstress(6)
      double precision ke_lstress(6)
      double precision xc_lstress(6)
      double precision ewald_lstress(6)
      double precision hartree_lstress(6)
      double precision local_lstress(6)
      double precision alpha_lstress(6)
      double precision RIJ_lstress(6)
      double precision ewaldG_lstress(6)

      integer date_len
      character*40 date_string

*     *** external functions ****
      logical  paw_noit_energy,paw_energy
      external paw_noit_energy,paw_energy

      call Parallel_taskid(taskid)
      oprint = (taskid.eq.MASTER)

*     **** initialize and read in geometry object ****
      value     = .true.
      value     = value.and.geom_create(geom,'geometry')
      value     = value.and.geom_rtdb_load(rtdb,geom,'geometry')
      value     = value.and.geom_ncent(geom,nion)
      value     = value.and.
     >            BA_alloc_get(mt_dbl,(3*nion),'r1',r1(2),r1(1))
      call load_cellgeometry(geom,dbl_mb(r1(1)),unita)

      if (.not.
     >    btdb_get(rtdb,'nwpw:lstress_numerical',
     >             mt_log,1,lstress_numerical))
     >  lstress_numerical = .false.

      if (.not.
     >    btdb_get(rtdb,'nwpw:noit_stress_numerical',
     >             mt_log,1,noit_energy))
     >  noit_energy = .true.



*     **** get finite difference step and initialize ****
      if (.not.
     >    btdb_get(rtdb,'nwpw:stress_alpha',mt_dbl,1,alpha))
     >  alpha = 0.00100d0


*     **** print out finite difference stress options ****
      if (oprint) then
         write(luout,*)
         if (lstress_numerical) then
            write(luout,*) "numerical finite difference paw lstress"
         else
            write(luout,*) "numerical finite difference paw stress"
         end if
         if (noit_energy) then
            write(luout,*) "  - noit_energy"
         else
            write(luout,*) "  - energy mimimized"
         end if
         write(luout,*) "  - alpha=",alpha
         write(luout,*)
      end if


*     ******************************
*     ******************************
*     ** finite difference stress **
*     ******************************
*     ******************************
      if (.not.lstress_numerical) then

      call dcopy(9,0.0d0,0,ke_stress,      1)
      call dcopy(9,0.0d0,0,xc_stress,      1)
      call dcopy(9,0.0d0,0,hartree_stress, 1)
      call dcopy(9,0.0d0,0,ewald_stress,   1)
      call dcopy(9,0.0d0,0,local_stress,   1)


      do j=1,3
      do i=1,3


*     ***********************
*     * + finite difference *
*     ***********************
      call dcopy(9,unita,1,unita2,1)
      unita2(i,j) = unita2(i,j) + alpha

*     **** update geometry object ****
      call util_date(date_string)
      date_len = inp_strlen(date_string)
      if (oprint) then
      write(luout,*)
      write(luout,
     &   '(1x,a12,i2,i2,1x,f8.3,1x,a4,1x,a10,f10.1,6x,a5,a27)')
     &   'amatrix i,j:',i,j,unita2(i,j),'(+)',
     &   'wall time:',util_wallsec(),
     &   'date:',date_string(1:date_len)
      call util_flush(6)
      end if
      call store_cellgeometry(geom,dbl_mb(r1(1)),unita2)
      value = value.and.geom_rtdb_store(rtdb,geom,'geometry')
      if (noit_energy) then
         value = value.and.paw_noit_energy(rtdb)
      else
         value = value.and.paw_energy(rtdb)
      end if


      value = value.and.
     >        btdb_get(rtdb,'paw:E_ke',mt_dbl,1,E)
      ke_stress(i,j) = E

      value = value.and.
     >        btdb_get(rtdb,'paw:E_xc',mt_dbl,1,E)
      xc_stress(i,j) = E

      value = value.and.
     >        btdb_get(rtdb,'paw:E_hartree',mt_dbl,1,E)
      hartree_stress(i,j) = E

      value = value.and.
     >        btdb_get(rtdb,'paw:E_mult',mt_dbl,1,E)
      ewald_stress(i,j) = E

      value = value.and.
     >        btdb_get(rtdb,'paw:E_local',mt_dbl,1,E)
      local_stress(i,j) = E



*     ***********************
*     * - finite difference *
*     ***********************
      call dcopy(9,unita,1,unita3,1)
      unita3(i,j) = unita3(i,j) - alpha

*     **** update geometry object ****
      call util_date(date_string)
      date_len = inp_strlen(date_string)
      if (oprint) then
      write(luout,*)
      write(luout,
     &   '(1x,a12,i2,i2,1x,f8.3,1x,a4,1x,a10,f10.1,6x,a5,a27)')
     &   'amatrix i,j:',i,j,unita3(i,j),'(-)',
     &   'wall time:',util_wallsec(),
     &   'date:',date_string(1:date_len)
      call util_flush(6)
      end if
      call store_cellgeometry(geom,dbl_mb(r1(1)),unita3)
      value = value.and.geom_rtdb_store(rtdb,geom,'geometry')
      if (noit_energy) then
         value = value.and.paw_noit_energy(rtdb)
      else
         value = value.and.paw_energy(rtdb)
      end if


      value = value.and.btdb_get(rtdb,'paw:E_ke',mt_dbl,1,E)
      ke_stress(i,j) = ke_stress(i,j) - E

      value = value.and.btdb_get(rtdb,'paw:E_xc',mt_dbl,1,E)
      xc_stress(i,j) = xc_stress(i,j) - E

      value = value.and.btdb_get(rtdb,'paw:E_hartree',mt_dbl,1,E)
      hartree_stress(i,j) = hartree_stress(i,j) - E

      value = value.and.btdb_get(rtdb,'paw:E_mult',mt_dbl,1,E)
      ewald_stress(i,j) = ewald_stress(i,j) - E

      value = value.and.btdb_get(rtdb,'paw:E_local',mt_dbl,1,E)
      local_stress(i,j) = local_stress(i,j) - E


      end do
      end do
      call dscal(9,(0.5d0/alpha),ke_stress,1)
      call dscal(9,(0.5d0/alpha),xc_stress,1)
      call dscal(9,(0.5d0/alpha),hartree_stress,1)
      call dscal(9,(0.5d0/alpha),ewald_stress,1)
      call dscal(9,(0.5d0/alpha),local_stress,1)

      call write_stress(ke_stress,"ke FD")
      call write_stress(hartree_stress,"coulomb FD")
      call write_stress(local_stress,"local FD")
      call write_stress(xc_stress,"xc FD")
      call write_stress(ewald_stress,"ewald FD")

      call dcopy(9,ke_stress,1,stress,1)
      call daxpy(9,1.0d0,xc_stress,1,stress,1)
      call daxpy(9,1.0d0,hartree_stress,1,stress,1)
      call daxpy(9,1.0d0,ewald_stress,1,stress,1)
      call daxpy(9,1.0d0,local_stress,1,stress,1)
      call write_stress(stress,"TOTAL FD")


*     ********************************
*     * output stress for use driver *
*     ********************************
      value = value.and.btdb_put(rtdb,'paw:stress',mt_dbl,9,stress)



*     *******************************
*     *******************************
*     ** finite difference lstress **
*     *******************************
*     *******************************
      else
         call cell_unita_abc_abg(unita,lattice)
         do i=1,6
*           ***********************
*           * + finite difference *
*           ***********************
            call dcopy(6,lattice,1,lattice2,1)
            lattice2(i) = lattice2(i) + alpha
            call cell_abc_abg_unita(lattice2,unita2)

*           **** update geometry object ****
            call util_date(date_string)
            date_len = inp_strlen(date_string)
            if (oprint) write(6,*)
            if (oprint) write(6,
     &      '(1x,a12,i2,1x,f8.3,1x,a4,1x,a10,f10.1,6x,a5,a27)')
     &      'lattice i:',i,lattice2(i),'(+)',
     &      'wall time:',util_wallsec(),
     &      'date:',date_string(1:date_len)
            if (oprint) call util_flush(6)
            call store_cellgeometry(geom,dbl_mb(r1(1)),unita2)
            value = value.and.geom_rtdb_store(rtdb,geom,'geometry')
            if (noit_energy) then
               value = value.and.paw_noit_energy(rtdb)
            else
               value = value.and.paw_energy(rtdb)
            end if

            value = value.and.
     >               btdb_get(rtdb,'paw:E_ke',mt_dbl,1,E)
            ke_lstress(i) = E

            value = value.and.
     >               btdb_get(rtdb,'paw:E_xc',mt_dbl,1,E)
            xc_lstress(i) = E

            value = value.and.
     >                  btdb_get(rtdb,'paw:E_hartree',mt_dbl,1,E)
            hartree_lstress(i) = E

            value = value.and.
     >                  btdb_get(rtdb,'paw:E_mult',mt_dbl,1,E)
            ewald_lstress(i) = E

            value = value.and.
     >               btdb_get(rtdb,'paw:E_local',mt_dbl,1,E)
            local_lstress(i) = E

*           ***********************
*           * - finite difference *
*           ***********************
            call dcopy(6,lattice,1,lattice3,1)
            lattice3(i) = lattice3(i) - alpha
            call cell_abc_abg_unita(lattice3,unita3)

*           **** update geometry object ****
            call util_date(date_string)
            date_len = inp_strlen(date_string)
            if (oprint) write(luout,*)
            if (oprint) write(luout,
     &         '(1x,a12,i2,1x,f8.3,1x,a4,1x,a10,f10.1,6x,a5,a27)')
     &         'lattice i:',i,lattice3(i),'(-)',
     &         'wall time:',util_wallsec(),
     &         'date:',date_string(1:date_len)
            if (oprint) call util_flush(6)
            call store_cellgeometry(geom,dbl_mb(r1(1)),unita3)
            value = value.and.geom_rtdb_store(rtdb,geom,'geometry')
            if (noit_energy) then
               value = value.and.paw_noit_energy(rtdb)
            else
               value = value.and.paw_energy(rtdb)
            end if

            value = value.and.
     >               btdb_get(rtdb,'paw:E_ke',mt_dbl,1,E)
            ke_lstress(i) = ke_lstress(i) - E

            value = value.and.
     >               btdb_get(rtdb,'paw:E_xc',mt_dbl,1,E)
            xc_lstress(i) = xc_lstress(i) - E

            value = value.and.
     >                  btdb_get(rtdb,'paw:E_hartree',mt_dbl,1,E)
            hartree_lstress(i) = hartree_lstress(i) - E

            value = value.and.
     >                  btdb_get(rtdb,'paw:E_mult',mt_dbl,1,E)
            ewald_lstress(i) = ewald_lstress(i) - E

            value = value.and.
     >               btdb_get(rtdb,'paw:E_local',mt_dbl,1,E)
            local_lstress(i) = local_lstress(i) - E
         end do
         call dscal(6,(0.5d0/alpha),ke_lstress,1)
         call dscal(6,(0.5d0/alpha),xc_lstress,1)
         call dscal(6,(0.5d0/alpha),hartree_lstress,1)
         call dscal(6,(0.5d0/alpha),ewald_lstress,1)
         call dscal(6,(0.5d0/alpha),local_lstress,1)

         call write_lstress(ke_lstress,"ke FD")
         call write_lstress(hartree_lstress,"coulomb FD")
         call write_lstress(local_lstress,"local FD")
         call write_lstress(xc_lstress,"xc FD")
         call write_lstress(ewald_lstress,"ewald FD")

         call dcopy(6,ke_lstress,1,lstress,1)
         call daxpy(6,1.0d0,xc_lstress,1,lstress,1)
         call daxpy(6,1.0d0,hartree_lstress,1,lstress,1)
         call daxpy(6,1.0d0,ewald_lstress,1,lstress,1)
         call daxpy(6,1.0d0,local_lstress,1,lstress,1)
         call write_lstress(lstress,"TOTAL FD")

*        *********************************
*        * output lstress for use driver *
*        *********************************
         value = value.and.
     >           btdb_put(rtdb,'paw:lstress',mt_dbl,6,lstress)
      end if


      call store_cellgeometry(geom,dbl_mb(r1(1)),unita)
      !value = value.and.geom_rtdb_delete(rtdb,'geometry')
      value = value.and.geom_rtdb_store(rtdb,geom,'geometry')
      value = value.and.geom_destroy(geom)
      value = value.and.BA_free_heap(r1(2))
 
      paw_stress_fd = value
      return 
      end
