
*     ***************************
*     *				*
*     *	  c_cgsd_energy_stress	*
*     *				*
*     ***************************

      subroutine c_cgsd_energy_stress(stress)
      implicit none
      real*8 stress(3,3)

      integer taskid,MASTER
      parameter (MASTER=0)

*     **** local variables ****
      integer u,v,s
      real*8  tstress(3,3),ht(3,3),scal
      real*8  sigma(3,3),xcstress(3,3)

*     **** external functions ****
      logical  psp_semicore
      real*8   lattice_omega,lattice_unita
      external psp_semicore
      external lattice_omega,lattice_unita


      call dcopy(9,0.0d0,0,stress,1)

*     **** Kinetic energy component : dE_kin/dhuv ****
      call cpsi_1ke_stress(tstress)
      call daxpy(9,1.0d0,tstress,1,stress,1)
      call write_stress(tstress,"ke")

*     **** Coulomb energy component : dE_Coul/dhuv ****
      call cpsi_1coulomb_stress(tstress)
      call daxpy(9,1.0d0,tstress,1,stress,1)
      call write_stress(tstress,"coulomb")

c*     **** Local pseudo energy component : dE_local/dhuv ****
c      call dng_1vlocal_stress(tstress)
c      call daxpy(9,1.0d0,tstress,1,stress,1)
c      call write_stress(tstress,"local")
c
c*     **** Nonlocal pseudo energy component : dE_nolocal/dhuv ****
c      call dcopy(9,0.0d0,0,tstress,1)
c      call psi_1vnonlocal_stress(tstress)
c      call daxpy(9,1.0d0,tstress,1,stress,1)
c      call write_stress(tstress,"nonlocal")
c
c*     **** xc energy component : dE_xc/dhuv ****
c      call dcopy(9,0.0d0,0,xcstress,1)
c      call rho_1exc_stress(tstress)
c      call daxpy(9,1.0d0,tstress,1,xcstress,1)
c      call write_stress(tstress,"xc1")
c
c*     **** Core-correction Coulomb energy component : dE_core/dhuv ****
c      if (psp_semicore(0)) then
c        call rho_1semicore_stress(tstress)
c        call write_stress(tstress,"xc2")
c        call daxpy(9,1.0d0,tstress,1,xcstress,1)
c      end if
c
c      call daxpy(9,1.0d0,xcstress,1,stress,1)
c      call write_stress(xcstress,"xc")
 
*     **** Ewald energy component : dE_ewald/dhuv ****
      call c_ewald_stress(tstress)
      call daxpy(9,1.0d0,tstress,1,stress,1)
      call write_stress(tstress,"ewald")

**     **** define ht ****
       do v=1,3
       do u=1,3
          ht(u,v)=lattice_unita(v,u)
       end do
       end do
 
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
      call write_stress(stress,"total gradient")

      call write_stress(sigma,"Internal Stress Tensor")
      return
      end

