c
c     $Id: bybminimize.F,v 1.14 2004-05-06 07:09:15 edo Exp $                       
c

*  ************************************************************
*  *                                                          *
*  *             Band by Band Kohn-Sham Minimizer             *
*  *                                                          *
*  *                                                          *
*  *                                                          *
*  ************************************************************
  
      subroutine bybminimize(E,deltae,deltac,current_iteration,
     >                       set_iterations,iterations)
      implicit none
      real*8     E(*)
      real*8     deltae,deltac
      integer    current_iteration
      logical    set_iterations
      integer    iterations


#include "stdio.fh"
#include "mafdecls.fh"
#include "util.fh"

*     **** local variables ****
      integer MASTER,taskid
      parameter (MASTER=0)

      real*8  deltat_min
      parameter (deltat_min=1.0d-3)
       
      integer vall_in(2),vall_out(2),vall_junk(2),rho_in(2)
      real*8  E0,dE0,deltae_old,Ein,deltae_history(10)
      real*8  ks_deltae,deltav,dV,deltav_old,diis_error
      integer nx,ny,nz,stalled_count


      real*8     sum0,sum1,sum3,scale,tole,tolc
      real*8     ehartree,eorbit,exc,pxc,eion
      real*8     Enew,Eold,alpha
      !common / cgsd_block / Enew,Eold,alpha

      integer it,it_in,i,ispin,bfgscount,icount,sd_it,cg_it

      logical value,precondition,done,stalled,deltav_bad(4),oprint
      integer n2ft3d
      real*8  e_ionmm,e_qmmm,e_mmmm,e_pol,e_vib,e_cav
      real*8  e_qmmm_e,e_qmmm_q,e_qmmm_lj,e_mmmm_q,e_mmmm_lj
      


*     **** external functions ****
      integer  control_ispin,control_scf_algorithm
      integer  control_it_in,control_it_out,psi_ne,control_version
      real*8   control_tole,control_tolc,control_ks_alpha
      real*8   rho_error,psi_1energy,psi_error
      real*8   dng_1ehartree,lattice_omega
      real*8   psi_1ke
      real*8   psi_1vl,psi_1v_field
      real*8   psi_1vnl
      real*8   rho_1exc
      real*8   rho_1pxc
      real*8   ewald_e,ion_ion_e
      real*8   psi_1eorbit
   
      external control_ispin,control_scf_algorithm,control_print
      logical control_print
      external control_it_in,control_it_out,psi_ne,control_version
      external control_tole,control_tolc,control_ks_alpha
      external rho_error,psi_1energy,psi_error
      external dng_1ehartree,lattice_omega
      external psi_1ke
      external psi_1vl,psi_1v_field
      external psi_1vnl
      external rho_1exc
      external rho_1pxc
      external ewald_e,ion_ion_e
      external psi_1eorbit

*     ***** QM/MM external functions ****
      logical  pspw_qmmm_found
      real*8   pspw_qmmm_Energy_QMMM_Q
      real*8   pspw_qmmm_Energy_QMMM_LJ
      real*8   pspw_qmmm_Energy_MMMM_Q
      real*8   pspw_qmmm_Energy_MMMM_LJ
      real*8   pspw_qmmm_Energy_pol
      real*8   pspw_qmmm_Energy_vib
      real*8   pspw_qmmm_Energy_cav
      external pspw_qmmm_found
      external pspw_qmmm_Energy_QMMM_Q
      external pspw_qmmm_Energy_QMMM_LJ
      external pspw_qmmm_Energy_MMMM_Q
      external pspw_qmmm_Energy_MMMM_LJ
      external pspw_qmmm_Energy_pol
      external pspw_qmmm_Energy_vib
      external pspw_qmmm_Energy_cav

*     ***** pspw_charge external functions ****
      logical  pspw_charge_found,control_precondition
      real*8   pspw_charge_Energy_ion,pspw_charge_Energy_charge
      external pspw_charge_found,control_precondition
      external pspw_charge_Energy_ion,pspw_charge_Energy_charge

      real*8   psi_1_noupdate_energy,psi_eigenvalue
      external psi_1_noupdate_energy,psi_eigenvalue
      double precision ddot
      external ddot
   
      Ein = E(1)
      call Parallel_taskid(taskid)
      oprint = (taskid.eq.MASTER).and.control_print(print_medium)

      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)
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
      ispin = control_ispin()
      deltav_old = 10.0d0
      deltav     = 0.0d0

      stalled       = .false.
      deltae_history(1) = 0.0d0
      deltae_history(2) = 0.0d0
      deltae_history(3) = 0.0d0
      deltae_history(4) = 0.0d0
      stalled_count     = 0

      call D3dB_n2ft3d(1,n2ft3d)

*     **** allocate rho_in and rho_out ****
      value = MA_push_get(mt_dbl,2*n2ft3d,
     >                     'vall_in',vall_in(2),vall_in(1))
      value = value.and.
     >        MA_push_get(mt_dbl,2*n2ft3d,
     >                     'vall_out',vall_out(2),vall_out(1))
      value = value.and.
     >        MA_push_get(mt_dbl,2*n2ft3d,
     >                   'vall_junk',vall_junk(2),vall_junk(1))
      value = value.and.
     >        MA_push_get(mt_dbl,2*n2ft3d,
     >                     'rho_in',rho_in(2),rho_in(1))
      if (.not. value) 
     >   call errquit('bybminimize:out of stack memory',0,0)

*     **** ion-ion energy ****
      eion = 0.0d0
      if (control_version().eq.3) eion = ewald_e()
      if (control_version().eq.4) eion = ion_ion_e()


*     **********************
*     **** bybminimizer ****
*     **********************


*     **** set the initial density ****
      if (current_iteration.eq.1) then
         Enew   = psi_1energy()
         alpha = control_ks_alpha()
         deltae = -9232323299.0d0
         ks_deltae = tole
         call electron_gen_vall()
         call electron_get_vall(dbl_mb(vall_out(1)))
         call electron_get_vall(dbl_mb(vall_in(1))) 

         call psi_1gen_hml()
         call psi_diagonalize_hml_assending()
         call psi_1rotate2()
         call psi_2to1()
      else
         call electron_get_vall(dbl_mb(vall_out(1)))
         call electron_get_vall(dbl_mb(vall_in(1))) 
         !call psi_get_density(1,dbl_mb(rho_in(1)))
      end if
  
*     **** iniitialize SCF Mixing ****    
      call nwpw_scf_mixing_init(control_scf_algorithm(),alpha,
     >                5,ispin,n2ft3d,dbl_mb(vall_out(1)))


*     ***** diis loop ****
      it   = 0
 2    it   = it + 1
      
*     **** diaganolize KS matrix ****
      call psi_KS_update(1,precondition,ks_deltae)

      call rho_1to2()
      Eold = Enew
      Enew = psi_1energy()

      deltae = Enew-Eold

      call electron_gen_vall()
      call electron_get_vall(dbl_mb(vall_in(1)))

*     **** compute deltaV ****
      call dcopy(ispin*n2ft3d,
     >           dbl_mb(vall_in(1)),1,
     >           dbl_mb(vall_junk(1)),1)
      call daxpy(ispin*n2ft3d,
     >              (-1.0d0),
     >              dbl_mb(vall_out(1)),1,
     >              dbl_mb(vall_junk(1)),1)
      deltav = ddot(ispin*n2ft3d,
     >                 dbl_mb(vall_junk(1)),1,
     >                 dbl_mb(vall_junk(1)),1)
      call D3dB_SumAll(deltav)
      deltav = deltav*dV
      


*     **** update vall using density mixing ****
      if ((it.le.0) .or.
     >    ((dabs(deltae).lt.1.0d5) .and.
     >    (deltav       .lt.1.0d5) .and.
     >    (.not.stalled          ))) then

         call nwpw_scf_mixing(dbl_mb(vall_in(1)),dbl_mb(vall_out(1)),
     >                        diis_error)

*     **** bad convergence - try fixed step steepest descent ****
      else

  30     call sdminimize(sd_it)
         Eold   = Enew
         Enew   = psi_1energy()
         !if (dabs(Enew-Eold).gt.1.0d2) go to 30

         call dcopy(ispin*n2ft3d,
     >              dbl_mb(vall_out(1)),1,
     >              dbl_mb(vall_junk(1)),1)

         call electron_gen_vall()
         call electron_get_vall(dbl_mb(vall_out(1))) 
         call nwpw_scf_mixing_reset(dbl_mb(vall_out(1)))


         call daxpy(ispin*n2ft3d,
     >              (-1.0d0),
     >              dbl_mb(vall_out(1)),1,
     >              dbl_mb(vall_junk(1)),1)
         deltav = ddot(ispin*n2ft3d,
     >                 dbl_mb(vall_junk(1)),1,
     >                 dbl_mb(vall_junk(1)),1)
         call D3dB_SumAll(deltav)
         deltav = deltav*dV

         call psi_1gen_hml()
         call psi_diagonalize_hml_assending()
         call psi_1rotate2()
         call psi_2to1()

         stalled       = .false.
         deltae_history(1) = 0.0d0
         deltae_history(2) = 0.0d0
         deltae_history(3) = 0.0d0
         deltae_history(4) = 0.0d0
         stalled_count     = 0

      end if
      call electron_set_vall(dbl_mb(vall_out(1))) 
      
          
*     **** tolerance checks ****
      deltae = Enew-Eold   
      deltac = rho_error()
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



*     **** free memory ****
      call nwpw_scf_mixing_end()
      value =           MA_pop_stack(rho_in(2))
      value = value.and.MA_pop_stack(vall_junk(2))
      value = value.and.MA_pop_stack(vall_out(2))
      value = value.and.MA_pop_stack(vall_in(2))
      if (.not. value) 
     >  call errquit('bybminimize: popping stack',1,0)

      call psi_check()


*     **** set return energies **** - This is duplicate code
      !Enew     = psi_1energy()
      eorbit   = psi_1eorbit()
      ehartree = dng_1ehartree()
      exc      = rho_1exc()
      pxc      = rho_1pxc()

      E(1)  = Enew + eion
      E(2)  = eorbit
      E(3)  = ehartree
      E(4)  = exc
      E(5)  = eion
      E(6)  = psi_1ke()
      E(7)  = psi_1vl()
      E(8)  = psi_1vnl()
      E(9)  = 2.0d0*ehartree
      E(10) = pxc
      if (pspw_qmmm_found()) then
*        **** get QM/MM  energies ****
         e_qmmm_e  = psi_1v_field()
         e_qmmm_q   = pspw_qmmm_Energy_QMMM_Q()
         e_qmmm_lj  = pspw_qmmm_Energy_QMMM_LJ()

         e_mmmm_q  = pspw_qmmm_Energy_MMMM_Q()
         e_mmmm_lj  = pspw_qmmm_Energy_MMMM_LJ()

         e_pol   = pspw_qmmm_Energy_pol()
         e_vib   = pspw_qmmm_Energy_vib()
         e_cav   = pspw_qmmm_Energy_cav()

         E(1) = E(1) + e_qmmm_q
     >               + e_qmmm_lj
     >               + e_mmmm_q
     >               + e_mmmm_lj
     >               + e_pol
     >               + e_vib
     >               + e_cav
         E(11) = e_qmmm_e
         E(12) = e_qmmm_q
         E(13) = e_qmmm_lj
         E(14) = e_mmmm_q
         E(15) = e_mmmm_lj
         E(16) = e_pol
         E(17) = e_vib
         E(18) = e_cav
      end if

*     **** get pspw_charge  energies ****
      if (pspw_charge_found()) then
         E(19)  = psi_1v_field()
         E(20)  = pspw_charge_Energy_ion()
         E(21)  = pspw_charge_Energy_charge()
         E(1)   = E(1) + E(20) + E(21)
      end if


      return
      end
 

