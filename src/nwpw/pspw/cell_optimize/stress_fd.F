*
* $Id: stress_fd.F,v 1.3 2002-04-26 20:37:41 bylaska Exp $
*
      logical function stress_fd(rtdb)
      implicit none
#include "mafdecls.fh"
#include "rtdb.fh"
      integer rtdb
      logical value
 

      integer l,i,j

      character*30 cell_name,rtdb_unita

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

*     *** external functions ****
      logical  pspw_noit_energy
      external pspw_noit_energy

*     **** cell: stuff ****
      value = rtdb_cget(rtdb,'cgsd:cell_name',1,cell_name)
      l = index(cell_name,' ') - 1
      rtdb_unita = cell_name(1:l)//':unita'
      value = value.and.
     >        rtdb_get(rtdb,rtdb_unita,mt_dbl,9,unita)

*     **** check unita for defaults ****
      call check_unita_for_default(rtdb,unita)

*     **** get fractional coordinates ****

      if (.not.
     >    rtdb_get(rtdb,'pspw:stress_alpha',mt_dbl,1,alpha))
     >  alpha = 0.00100d0
      call dcopy(9,0.0d0,0,ke_stress,      1)
      call dcopy(9,0.0d0,0,xc_stress,      1)
      call dcopy(9,0.0d0,0,hartree_stress, 1)
      call dcopy(9,0.0d0,0,ewald_stress,   1)
      call dcopy(9,0.0d0,0,local_stress,   1)
      call dcopy(9,0.0d0,0,nonlocal_stress,1)


      do j=1,3
      do i=1,3


*     *****************************
*     **** + finite difference ****
*     *****************************
      call dcopy(9,unita,1,unita2,1)
      unita2(i,j) = unita2(i,j) + alpha

*     **** update geometry object ****
      write(*,*) "Ion positions adjusted......"
      call adjust_cellgeometry(rtdb,unita,unita2)
      value = value.and.
     >        rtdb_put(rtdb,rtdb_unita,mt_dbl,9,unita2)

      value = value.and.pspw_noit_energy(rtdb)

      call adjust_cellgeometry(rtdb,unita2,unita)
      value = value.and.
     >        rtdb_put(rtdb,rtdb_unita,mt_dbl,9,unita)


      value = value.and.
     >        rtdb_get(rtdb,'pspw:E_ke',mt_dbl,1,E)
      ke_stress(i,j) = E

      value = value.and.
     >        rtdb_get(rtdb,'pspw:E_xc',mt_dbl,1,E)
      xc_stress(i,j) = E

      value = value.and.
     >        rtdb_get(rtdb,'pspw:E_hartree',mt_dbl,1,E)
      hartree_stress(i,j) = E

      value = value.and.
     >        rtdb_get(rtdb,'pspw:E_ewald',mt_dbl,1,E)
      ewald_stress(i,j) = E

      value = value.and.
     >        rtdb_get(rtdb,'pspw:E_local',mt_dbl,1,E)
      local_stress(i,j) = E

      value = value.and.
     >        rtdb_get(rtdb,'pspw:E_nonlocal',mt_dbl,1,E)
      nonlocal_stress(i,j) = E


*     *****************************
*     **** - finite difference ****
*     *****************************
      call dcopy(9,unita,1,unita3,1)
      unita3(i,j) = unita3(i,j) - alpha

*     **** update geometry object ****
      write(*,*) "Ion positions adjusted......"
      call adjust_cellgeometry(rtdb,unita,unita3)
      value = value.and.
     >        rtdb_put(rtdb,rtdb_unita,mt_dbl,9,unita3)

      value = pspw_noit_energy(rtdb)

      call adjust_cellgeometry(rtdb,unita3,unita)
      value = value.and.
     >        rtdb_put(rtdb,rtdb_unita,mt_dbl,9,unita)

      value = value.and.
     >        rtdb_get(rtdb,'pspw:E_ke',mt_dbl,1,E)
      ke_stress(i,j) = ke_stress(i,j) - E

      value = value.and.
     >        rtdb_get(rtdb,'pspw:E_xc',mt_dbl,1,E)
      xc_stress(i,j) = xc_stress(i,j) - E

      value = value.and.
     >        rtdb_get(rtdb,'pspw:E_hartree',mt_dbl,1,E)
      hartree_stress(i,j) = hartree_stress(i,j) - E

      value = value.and.
     >        rtdb_get(rtdb,'pspw:E_ewald',mt_dbl,1,E)
      ewald_stress(i,j) = ewald_stress(i,j) - E

      value = value.and.
     >        rtdb_get(rtdb,'pspw:E_local',mt_dbl,1,E)
      local_stress(i,j) = local_stress(i,j) - E

      value = value.and.
     >        rtdb_get(rtdb,'pspw:E_nonlocal',mt_dbl,1,E)
      nonlocal_stress(i,j) = nonlocal_stress(i,j) - E

      end do
      end do
      call dscal(9,(0.5d0/alpha),ke_stress,1)
      call dscal(9,(0.5d0/alpha),xc_stress,1)
      call dscal(9,(0.5d0/alpha),hartree_stress,1)
      call dscal(9,(0.5d0/alpha),ewald_stress,1)
      call dscal(9,(0.5d0/alpha),local_stress,1)
      call dscal(9,(0.5d0/alpha),nonlocal_stress,1)

      call write_stress(ke_stress,"ke FD")
      call write_stress(xc_stress,"xc FD")
      call write_stress(hartree_stress,"hartree FD")
      call write_stress(ewald_stress,"ewald FD")
      call write_stress(local_stress,"local FD")
      call write_stress(nonlocal_stress,"nonlocal FD")

      call dcopy(9,ke_stress,1,stress,1)
      call daxpy(9,1.0d0,xc_stress,1,stress,1)
      call daxpy(9,1.0d0,hartree_stress,1,stress,1)
      call daxpy(9,1.0d0,ewald_stress,1,stress,1)
      call daxpy(9,1.0d0,local_stress,1,stress,1)
      call daxpy(9,1.0d0,nonlocal_stress,1,stress,1)
      call write_stress(stress,"TOTAL FD")
 
      stress_fd = value
      return 
      end
