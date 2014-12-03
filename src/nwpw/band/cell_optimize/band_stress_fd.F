*
* $Id$
*
      logical function band_stress_fd(rtdb)
      implicit none
      integer rtdb
      logical value

#include "bafdecls.fh"
#include "btdb.fh"
#include "inp.fh"
#include "util.fh"
#include "geom.fh"

      integer taskid,MASTER
      parameter (MASTER=0)
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
      double precision nonlocal_stress(3,3)
      double precision alpha_stress(3,3)
      double precision RIJ_stress(3,3)
      double precision ewaldG_stress(3,3)
      double precision disp_stress(3,3)

      integer date_len
      character*40 date_string

*     *** external functions ****
      logical  band_noit_energy
      external band_noit_energy

      call Parallel_taskid(taskid)

*     **** initialize and read in geometry object ****
      value     = geom_create(geom,'geometry')
      value     = value.and.geom_rtdb_load(rtdb,geom,'geometry')
      value     = value.and.geom_ncent(geom,nion)
      value     = value.and.
     >            BA_alloc_get(mt_dbl,(3*nion),'r1',r1(2),r1(1))
      call load_cellgeometry(geom,dbl_mb(r1(1)),unita)


*     **** get finite difference step and initialize ****
      if (.not.
     >    btdb_get(rtdb,'pspw:stress_alpha',mt_dbl,1,alpha))
     >  alpha = 0.00100d0
      call dscal(9,0.0d0,ke_stress,      1)
      call dscal(9,0.0d0,xc_stress,      1)
      call dscal(9,0.0d0,hartree_stress, 1)
      call dscal(9,0.0d0,ewald_stress,   1)
      call dscal(9,0.0d0,local_stress,   1)
      call dscal(9,0.0d0,nonlocal_stress,1)
      call dscal(9,0.0d0,disp_stress,1)


      do j=1,3
      do i=1,3


*     *****************************
*     **** + finite difference ****
*     *****************************
      call dcopy(9,unita,1,unita2,1)
      unita2(i,j) = unita2(i,j) + alpha

*     **** update geometry object ****
      call util_date(date_string)
      date_len = inp_strlen(date_string)
      if (taskid.eq.MASTER) then
         write(6,*)
         write(6,
     &   '(1x,a12,i2,i2,1x,f8.3,1x,a4,1x,a10,f10.1,6x,a5,a27)')
     &   'amatrix i,j:',i,j,unita2(i,j),'(+)',
     &   'wall time:',util_wallsec(),
     &   'date:',date_string(1:date_len)
         call util_flush(6)
      end if
      call store_cellgeometry(geom,dbl_mb(r1(1)),unita2)
      value = value.and.geom_rtdb_store(rtdb,geom,'geometry')
      value = value.and.band_noit_energy(rtdb)


      value = value.and.
     >        btdb_get(rtdb,'band:E_ke',mt_dbl,1,E)
      ke_stress(i,j) = E

      value = value.and.
     >        btdb_get(rtdb,'band:E_xc',mt_dbl,1,E)
      xc_stress(i,j) = E

      value = value.and.
     >        btdb_get(rtdb,'band:E_hartree',mt_dbl,1,E)
      hartree_stress(i,j) = E

      value = value.and.
     >        btdb_get(rtdb,'band:E_ewald',mt_dbl,1,E)
      ewald_stress(i,j) = E

      value = value.and.
     >        btdb_get(rtdb,'band:E_local',mt_dbl,1,E)
      local_stress(i,j) = E

      value = value.and.
     >        btdb_get(rtdb,'band:E_nonlocal',mt_dbl,1,E)
      nonlocal_stress(i,j) = E

      value = value.and.
     >        btdb_get(rtdb,'band:E_disp',mt_dbl,1,E)
      disp_stress(i,j) = E

*     *****************************
*     **** - finite difference ****
*     *****************************
      call dcopy(9,unita,1,unita3,1)
      unita3(i,j) = unita3(i,j) - alpha

*     **** update geometry object ****
      call util_date(date_string)
      date_len = inp_strlen(date_string)
      if (taskid.eq.MASTER) then
         write(6,*)
         write(6,
     &   '(1x,a12,i2,i2,1x,f8.3,1x,a4,1x,a10,f10.1,6x,a5,a27)')
     &   'amatrix i,j:',i,j,unita3(i,j),'(-)',
     &   'wall time:',util_wallsec(),
     &   'date:',date_string(1:date_len)
         call util_flush(6)
      end if
      call store_cellgeometry(geom,dbl_mb(r1(1)),unita3)
      value = value.and.geom_rtdb_store(rtdb,geom,'geometry')
    
      value = band_noit_energy(rtdb)


      value = value.and.btdb_get(rtdb,'band:E_ke',mt_dbl,1,E)
      ke_stress(i,j) = ke_stress(i,j) - E

      value = value.and.btdb_get(rtdb,'band:E_xc',mt_dbl,1,E)
      xc_stress(i,j) = xc_stress(i,j) - E

      value = value.and.btdb_get(rtdb,'band:E_hartree',mt_dbl,1,E)
      hartree_stress(i,j) = hartree_stress(i,j) - E

      value = value.and.btdb_get(rtdb,'band:E_ewald',mt_dbl,1,E)
      ewald_stress(i,j) = ewald_stress(i,j) - E

      value = value.and.btdb_get(rtdb,'band:E_local',mt_dbl,1,E)
      local_stress(i,j) = local_stress(i,j) - E

      value = value.and.btdb_get(rtdb,'band:E_nonlocal',mt_dbl,1,E)
      nonlocal_stress(i,j) = nonlocal_stress(i,j) - E

      value = value.and.btdb_get(rtdb,'band:E_disp',mt_dbl,1,E)
      disp_stress(i,j) = disp_stress(i,j) - E

      end do
      end do
      call dscal(9,(0.5d0/alpha),ke_stress,1)
      call dscal(9,(0.5d0/alpha),xc_stress,1)
      call dscal(9,(0.5d0/alpha),hartree_stress,1)
      call dscal(9,(0.5d0/alpha),ewald_stress,1)
      call dscal(9,(0.5d0/alpha),local_stress,1)
      call dscal(9,(0.5d0/alpha),nonlocal_stress,1)
      call dscal(9,(0.5d0/alpha),disp_stress,1)

      call write_stress(ke_stress,"ke FD")
      call write_stress(hartree_stress,"coulomb FD")
      call write_stress(local_stress,"local FD")
      call write_stress(nonlocal_stress,"nonlocal FD")
      call write_stress(xc_stress,"xc FD")
      call write_stress(ewald_stress,"ewald FD")
      call write_stress(disp_stress,"disp FD")

      call dcopy(9,ke_stress,1,stress,1)
      call daxpy(9,1.0d0,xc_stress,1,stress,1)
      call daxpy(9,1.0d0,hartree_stress,1,stress,1)
      call daxpy(9,1.0d0,ewald_stress,1,stress,1)
      call daxpy(9,1.0d0,local_stress,1,stress,1)
      call daxpy(9,1.0d0,nonlocal_stress,1,stress,1)
      call daxpy(9,1.0d0,disp_stress,1,stress,1)
      call write_stress(stress,"TOTAL FD")


*     **************************************
*     **** output stress for use driver ****
*     **************************************
      value = value.and.btdb_put(rtdb,'band:stress',mt_dbl,9,stress)


      call store_cellgeometry(geom,dbl_mb(r1(1)),unita)
      !value = value.and.geom_rtdb_delete(rtdb,'geometry')
      value = value.and.geom_rtdb_store(rtdb,geom,'geometry')
      value = value.and.geom_destroy(geom)
      value = value.and.BA_free_heap(r1(2))
 
      band_stress_fd = value
      return 
      end
