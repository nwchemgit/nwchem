
*     ***************************
*     *				*
*     *	  c_cgsd_energy_stress	*
*     *				*
*     ***************************

      subroutine c_cgsd_energy_stress(stress,lstress)
      implicit none
      real*8 stress(3,3)
      real*8 lstress(6)

#include "stdio.fh"
#include "util.fh"

*     **** local variables ****
      integer taskid,MASTER
      parameter (MASTER=0)

      logical oprint,hprint
      integer u,v,s
      real*8  tstress(3,3),ht(3,3),scal
      real*8  sigma(3,3),xcstress(3,3)

*     **** external functions ****
      logical  cpsp_semicore,control_print
      real*8   lattice_omega,lattice_unita
      external cpsp_semicore,control_print
      external lattice_omega,lattice_unita
      logical  ion_disp_on
      external ion_disp_on

      call Parallel_taskid(taskid)
      oprint = ((taskid.eq.MASTER).and.control_print(print_medium))
      hprint = ((taskid.eq.MASTER).and.control_print(print_high))


      call dcopy(9,0.0d0,0,stress,1)

*     **** Kinetic energy component : dE_kin/dhuv ****
      call cpsi_1ke_stress(tstress)
      call daxpy(9,1.0d0,tstress,1,stress,1)
      if (hprint) call write_stress(tstress,"ke")

*     **** Coulomb energy component : dE_Coul/dhuv ****
      call cpsi_1coulomb_stress(tstress)
      call daxpy(9,1.0d0,tstress,1,stress,1)
      if (hprint) call write_stress(tstress,"coulomb")

*     **** Local pseudo energy component : dE_local/dhuv ****
      call c_dng_1vlocal_stress(tstress)
      call daxpy(9,1.0d0,tstress,1,stress,1)
      if (hprint) call write_stress(tstress,"local")

*     **** Nonlocal pseudo energy component : dE_nolocal/dhuv ****
      call dcopy(9,0.0d0,0,tstress,1)
      call cpsi_1vnonlocal_stress(tstress)
      call daxpy(9,1.0d0,tstress,1,stress,1)
      if (hprint) call write_stress(tstress,"nonlocal")

*     **** xc energy component : dE_xc/dhuv ****
      call dcopy(9,0.0d0,0,tstress,1)
      call c_rho_1exc_stress(tstress)
      call dcopy(9,tstress,1,xcstress,1)
      if (hprint) call write_stress(tstress,"xc1")

*     **** Core-correction Coulomb energy component : dE_core/dhuv ****
      if (cpsp_semicore(0)) then
        call c_rho_1semicore_stress(tstress)
        if (hprint) call write_stress(tstress,"xc2")
        call daxpy(9,1.0d0,tstress,1,xcstress,1)
      end if

      call daxpy(9,1.0d0,xcstress,1,stress,1)
      if (hprint) call write_stress(xcstress,"xc")
 
*     **** Ewald energy component : dE_ewald/dhuv ****
      call ewald_stress(tstress)
      call daxpy(9,1.0d0,tstress,1,stress,1)
      if (hprint) call write_stress(tstress,"ewald")

*     **** Dispersion energy component : dE_disp/dhuv ****
      if (ion_disp_on()) then
         call ion_disp_stress(tstress)
        call daxpy(9,1.0d0,tstress,1,stress,1)
        if (hprint) call write_stress(tstress,"Dispersion")
      end if



**     **** define ht ****
       do v=1,3
       do u=1,3
          ht(u,v)=lattice_unita(v,u)
       end do
       end do

       call cell_lattice_gradient(stress,lstress)

*      **** apply lattice_vector constraints ****
       call zero_cell_constraint(stress,lstress)

 
**     **** define stress tensor ****
       call dcopy(9,0.0d0,0,sigma,1)
       scal = -1.0d0/lattice_omega()
       do v=1,3
       do u=1,3
          do s=1,3
             sigma(u,v) = sigma(u,v) + scal*stress(u,s)*ht(s,v)
          end do
       end do
       end do

c     *** call dscal(9,-1.0d0,stress,1)
      if (oprint) call write_stress(stress,"total gradient")

      if (oprint) then
        write(luout,*)
        write(luout,1811) "dE/da     = ",lstress(1)
        write(luout,1811) "dE/db     = ",lstress(2)
        write(luout,1811) "dE/dc     = ",lstress(3)
        write(luout,1811) "dE/dalpha = ",lstress(4)
        write(luout,1811) "dE/dbeta  = ",lstress(5)
        write(luout,1811) "dE/dgamma = ",lstress(6)
        write(luout,*)
      end if

      if (hprint) call write_stress(sigma,"Internal Stress Tensor")
      return
 1800 FORMAT('Lattice Constraint: gradient(',i1,',',i1,')-->0')
 1801 FORMAT('Lattice Constraint: gradient(',i1,')-->0')
 1811 FORMAT(5X,A14,F11.5)
      end

c $Id$
