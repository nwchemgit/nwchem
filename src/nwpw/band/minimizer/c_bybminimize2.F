c
c     $Id$                       
c

*  ************************************************************
*  *                                                          *
*  *             Band by Band Kohn-Sham Minimizer             *
*  *                  (density mixing)                        *
*  *                                                          *
*  *                                                          *
*  ************************************************************
  
      subroutine c_bybminimize2(E,deltae,deltac,current_iteration,
     >                       set_iterations,iterations)
      implicit none
      real*8     E(*)
      real*8     deltae,deltac
      integer    current_iteration
      logical    set_iterations
      integer    iterations


#include "stdio.fh"
#include "bafdecls.fh"
#include "util.fh"

*     **** local variables ****
      integer MASTER,taskid
      parameter (MASTER=0)

      real*8  deltat_min
      parameter (deltat_min=1.0d-3)
       
      integer rho_in(2),rho_out(2),rho_junk(2)
      real*8  Ein,deltae_history(10),ehfx,phfx
      real*8  ks_deltae,deltav,dV,deltav_old,diis_error
      real*8  charge(2),tcharge
      integer nx,ny,nz,stalled_count,ms,nb,i


      real*8     tole,tolc
      real*8     ehartree,eorbit,exc,pxc,eion
      real*8     Enew,Eold,alpha
      common / c_cgsd_block / Enew,Eold,alpha

      integer it,it_in,ispin,sd_it,cg_it

      logical value,precondition,done,stalled,oprint
      integer nfft3d
      


*     **** external functions ****
      logical  band_HFX
      logical  control_print,control_precondition,control_fractional
      integer  control_ispin,control_scf_algorithm
      integer  control_it_in,control_it_out,control_version
      integer  cpsi_ispin,cpsi_ne,cpsi_nbrillq
      real*8   control_tole,control_tolc,control_ks_alpha
      real*8   c_rho_error,cpsi_1energy
      real*8   c_dng_1ehartree,lattice_omega
      real*8   cpsi_1ke
      real*8   cpsi_1vl
      real*8   cpsi_1vnl
      real*8   cpsi_1vnlso
      real*8   c_rho_1exc
      real*8   c_rho_1pxc
      real*8   ewald_e
      real*8   cpsi_1eorbit,cpsi_smearfermi,cpsi_smearcorrection
      real*8   cpsi_occupation,brillioun_weight
   
      external band_HFX
      external control_print,control_precondition,control_fractional
      external control_ispin,control_scf_algorithm
      external control_it_in,control_it_out,control_version
      external control_tole,control_tolc,control_ks_alpha
      external cpsi_ispin,cpsi_ne,cpsi_nbrillq
      external c_rho_error,cpsi_1energy
      external c_dng_1ehartree,lattice_omega
      external cpsi_1ke
      external cpsi_1vl,psi_1v_field
      external cpsi_1vnl
      external cpsi_1vnlso
      external c_rho_1exc
      external c_rho_1pxc
      external ewald_e
      external cpsi_1eorbit,cpsi_smearfermi,cpsi_smearcorrection
      external cpsi_occupation,brillioun_weight

      real*8   cpsi_1_noupdate_energy,cpsi_eigenvalue,ion_disp_energy
      external cpsi_1_noupdate_energy,cpsi_eigenvalue,ion_disp_energy
      logical  ion_disp_on
      external ion_disp_on

      logical  cpsi_try_read_density,control_use_fractional_rho
      external cpsi_try_read_density,control_use_fractional_rho
   

      Ein = E(1)
      call Parallel_taskid(taskid)
      oprint = (taskid.eq.MASTER).and.control_print(print_medium)

      call C3dB_nx(1,nx)
      call C3dB_ny(1,ny)
      call C3dB_nz(1,nz)
      dV = lattice_omega()/dble(nx*ny*nz)
      if (set_iterations) then
        it_in = iterations
        sd_it = 2 
        cg_it = 1 
      else
        it_in = control_it_in()*control_it_out()
        sd_it = 10
        cg_it = 10
      end if
      tole  = control_tole()
      tolc  = control_tolc()
      precondition = control_precondition()
      ispin = cpsi_ispin()
      deltav_old = 10.0d0
      deltav     = 0.0d0

      stalled       = .false.
      deltae_history(1) = 0.0d0
      deltae_history(2) = 0.0d0
      deltae_history(3) = 0.0d0
      deltae_history(4) = 0.0d0
      stalled_count     = 0

      call C3dB_nfft3d(1,nfft3d)

*     **** allocate rho_in and rho_out ****
      value = BA_push_get(mt_dbl,2*nfft3d,
     >                     'rho_in',rho_in(2),rho_in(1))
      value = value.and.
     >        BA_push_get(mt_dbl,2*nfft3d,
     >                     'rho_out',rho_out(2),rho_out(1))
      value = value.and.
     >        BA_push_get(mt_dbl,2*nfft3d,
     >                   'rho_junk',rho_junk(2),rho_junk(1))
      if (.not. value) 
     >   call errquit('c_bybminimize2:out of stack memory',0,0)

*     **** ion-ion energy ****
      eion = ewald_e()


*     **********************
*     **** bybminimizer ****
*     **********************


*     **** set the initial density ****
      if (current_iteration.eq.1) then
         Enew  = cpsi_1energy()
         alpha = control_ks_alpha()
         deltae = -9232323299.0d0
         ks_deltae = tole
         value = .false.
         if (control_fractional()) then
            if (control_use_fractional_rho()) then
               value = cpsi_try_read_density(1)
            end if
         end if
         call c_electron_gen_vall()
         call cpsi_get_density(1,dbl_mb(rho_in(1)))
         call cpsi_get_density(1,dbl_mb(rho_out(1)))

         call cpsi_1gen_hml()
         call cpsi_diagonalize_hml()
         if (.not.value) then
            !call cpsi_1gen_hml()
            !call cpsi_diagonalize_hml()
            call cpsi_1rotate2()
            call cpsi_2to1()
         end if
         if (control_fractional()) then
            call cpsi_1define_occupation(0.0d0,.false.)
            Enew = Enew + cpsi_smearcorrection()
         end if
      else
         call cpsi_get_density(1,dbl_mb(rho_in(1)))
         call cpsi_get_density(1,dbl_mb(rho_out(1)))
      end if

      !*** set the density ***
      call cpsi_set_density(1,dbl_mb(rho_out(1)))

*     **** set the charge from the occupations ****
      do ms=1,ispin
         charge(ms) = 0.0d0
         do nb=1,cpsi_nbrillq()
         do i=1,cpsi_ne(ms)
            charge(ms) = charge(ms)
     >                 + cpsi_occupation(nb,ms,i)*brillioun_weight(nb)
         end do
         end do
         call K1dB_SumAll(charge(ms))
      end do
  
*     **** iniitialize SCF Mixing ****    
      call nwpw_scf_mixing_init(control_scf_algorithm(),alpha,
     >                5,ispin,nfft3d,dbl_mb(rho_out(1)))


*     ***** scf loop ****
      it   = 0
 2    it   = it + 1
      
*     **** diaganolize KS matrix ****
      call cpsi_KS_update(1,precondition,ks_deltae)
      call c_rho_1to2()

*     **** rotate orbitals or define fractional occupation **** 
      call cpsi_1gen_hml()
      call cpsi_diagonalize_hml()
      call cpsi_1rotate2()
      call cpsi_2to1()

      if (control_fractional()) then
         call cpsi_1define_occupation(-1.0d0,.false.)
      end if

      Eold = Enew
      Enew = cpsi_1energy()
      if (control_fractional()) Enew = Enew + cpsi_smearcorrection()
      deltae = Enew-Eold


      call cpsi_get_density(1,dbl_mb(rho_in(1)))

*     **** compute delta rho ****
      call dcopy(ispin*nfft3d,
     >           dbl_mb(rho_in(1)),1,
     >           dbl_mb(rho_junk(1)),1)
      call daxpy(ispin*nfft3d,
     >              (-1.0d0),
     >              dbl_mb(rho_out(1)),1,
     >              dbl_mb(rho_junk(1)),1)
      deltav = ddot(ispin*nfft3d,
     >                 dbl_mb(rho_junk(1)),1,
     >                 dbl_mb(rho_junk(1)),1)
      call C3dB_SumAll(deltav)
      deltav = deltav*dV
      


*     **** update rho using density mixing ****
      if ((it.le.0) .or.
     >    ((dabs(deltae).lt.1.0d5) .and.
     >    (deltav       .lt.1.0d5) .and.
     >    (.not.stalled          ))) then

         call nwpw_scf_mixing(dbl_mb(rho_in(1)),dbl_mb(rho_out(1)),
     >                        deltae,diis_error)


*     **** bad convergence - try fixed step steepest descent ****
      else

  30     call c_sdminimize(sd_it)
         Eold   = Enew
         Enew   = cpsi_1energy()
         if (control_fractional()) Enew = Enew + cpsi_smearcorrection()
         !if (dabs(Enew-Eold).gt.1.0d2) go to 30

         call dcopy(ispin*2*nfft3d,
     >              dcpl_mb(rho_out(1)),1,
     >              dcpl_mb(rho_junk(1)),1)

         call cpsi_get_density(1,dbl_mb(rho_out(1)))
         call nwpw_scf_mixing_reset(dbl_mb(rho_out(1)))


         call daxpy(ispin*nfft3d,
     >              (-1.0d0),
     >              dbl_mb(rho_out(1)),1,
     >              dbl_mb(rho_junk(1)),1)
         deltav = ddot(ispin*nfft3d,
     >                 dbl_mb(rho_junk(1)),1,
     >                 dbl_mb(rho_junk(1)),1)
         call C3dB_SumAll(deltav)
         deltav = deltav*dV

         call cpsi_1gen_hml()
         call cpsi_diagonalize_hml()
         call cpsi_1rotate2()
         call cpsi_2to1()

c*        **** define fractional occupation ****
c         if (control_fractional()) 
c     >      call cpsi_1define_occupation(1.0d0,.false.)

         stalled       = .false.
         deltae_history(1) = 0.0d0
         deltae_history(2) = 0.0d0
         deltae_history(3) = 0.0d0
         deltae_history(4) = 0.0d0
         stalled_count     = 0

      end if


      !*** make sure charge is positive and normalized properly ***
      do ms=1,ispin
         call C3dB_r_abs1(1,dbl_mb(rho_out(1)+(ms-1)*nfft3d))
         call C3dB_r_dsum(1,dbl_mb(rho_out(1)+(ms-1)*nfft3d),tcharge)
         tcharge = tcharge*dV
         !write(*,*) "tcharge=",ms,tcharge,charge(ms)
         if (dabs(tcharge-charge(ms)).gt.1.0d-10) 
     >      call C3dB_r_SMul1(1,(charge(ms)/tcharge),
     >                        dbl_mb(rho_out(1)+(ms-1)*nfft3d))
      end do


      !*** set the density ***
      call cpsi_set_density(1,dbl_mb(rho_out(1)))
      
          
*     **** tolerance checks ****
      deltae = Enew-Eold   
      deltac = c_rho_error()
      E(1)   = Enew+eion

      if ((oprint).and.(.not.set_iterations)) then 
        write(luout,1310) it,E(1),deltae,deltac,deltav
        call util_flush(luout)
      end if
 1310 FORMAT(I8,E20.10,3E15.5)


      !**** set ks_deltae ****
      ks_deltae = 0.001d0*dabs(deltae)
      if (ks_deltae.lt.(0.1d0*tole)) ks_deltae = 0.1d0*tole
      if (ks_deltae.gt.1.0d-4) ks_deltae = 1.0d-4
      !ks_deltae = 0.1d0*tole



      deltav_old = deltav

      deltae_history(1)    = deltae_history(2)
      deltae_history(2)    = deltae_history(3)
      deltae_history(3)    = deltae_history(4)
      deltae_history(4)    = deltae        

      if (stalled_count .gt.4) then
        stalled = (deltae_history(4)
     >            +deltae_history(3)
     >            +deltae_history(2)
     >            +deltae_history(1)).gt.0.0d0
      else
         stalled = .false.
      end if
      stalled_count = stalled_count + 1
      stalled = .false.

      precondition = precondition.and.(dabs(deltae).gt.1*tole)

      done = ( (    (dabs(deltae).lt.tole)
     >         .and.(deltae.lt.0.0d0)
     >         .and.(deltac      .lt.tolc))
     >       .or. (it.ge.it_in))

      if (.not.done) go to 2


      if (control_fractional()) then
         if (control_use_fractional_rho()) then
            call cpsi_write_density(1)
         end if
      end if
      call c_rho_1to2()

*     **** free memory ****
      call nwpw_scf_mixing_end()
      value =           BA_pop_stack(rho_junk(2))
      value = value.and.BA_pop_stack(rho_out(2))
      value = value.and.BA_pop_stack(rho_in(2))
      if (.not. value) 
     >  call errquit('c_bybminimize: popping stack',1,0)

c      call cpsi_check()


*     **** set return energies **** - This is duplicate code
      !Enew     = cpsi_1energy()
      eorbit   = cpsi_1eorbit()
      ehartree = c_dng_1ehartree()
      exc      = c_rho_1exc()
      pxc      = c_rho_1pxc()

      E(1)  = Enew + eion
      E(2)  = eorbit
      E(3)  = ehartree
      E(4)  = exc
      E(5)  = eion
      E(6)  = cpsi_1ke()
      E(7)  = cpsi_1vl()
      E(8)  = cpsi_1vnl()
      E(9)  = 2.0d0*ehartree
      E(10) = pxc
      E(11) = cpsi_1vnlso()

      if (control_fractional()) then
        E(12) = cpsi_smearcorrection()
      end if

*     **** HFX terms ****
      if (band_HFX()) then
         call c_electron_HFX_energies(ehfx,phfx)
         E(26) = ehfx
         E(27) = phfx
      end if

*     **** Dispersion energy ****
      if (ion_disp_on()) then
         E(33) = ion_disp_energy()
         E(1)  = E(1) + E(33)
      end if


      return
      end
 

