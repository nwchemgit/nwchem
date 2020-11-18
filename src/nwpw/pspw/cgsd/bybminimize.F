c
c     $Id$                       
c

*  ************************************************************
*  *                                                          *
*  *             Band by Band Kohn-Sham Minimizer             *
*  *                                                          *
*  *                                                          *
*  *                                                          *
*  ************************************************************
  
      subroutine bybminimize(E,deltae,deltac,current_iteration,
     >                       set_iterations,iterations,failed)
      implicit none
      real*8     E(*)
      real*8     deltae,deltac
      integer    current_iteration
      logical    set_iterations
      integer    iterations
      logical    failed


#include "stdio.fh"
#include "bafdecls.fh"
#include "util.fh"

*     **** local variables ****
      integer MAX_SD_COUNT
      parameter (MAX_SD_COUNT = 3)
      integer MASTER,taskid
      parameter (MASTER=0)

      real*8  deltat_min
      parameter (deltat_min=1.0d-3)
       
      integer vall_in(2),vall_out(2),vall_junk(2),rho_in(2)
      real*8  E0,dE0,deltae_old,Ein,deltae_history(10)
      real*8  ks_deltae,deltav,dV,deltav_old,diis_error,e00
      integer nx,ny,nz,stalled_count,sd_count


      real*8     tole,tolc
      real*8     ehartree,eorbit,exc,pxc,eion
      real*8     Enew,Eold,alpha
      !common / cgsd_block / Enew,Eold,alpha

      integer it,it_in,i,ispin,bfgscount,icount,sd_it,cg_it

      logical value,precondition,done,stalled,deltav_bad(4),oprint
      logical ks_block
      integer n2ft3d,n2ft3d_map
      !real*8  e_ionmm,e_qmmm,e_mmmm,e_pol,e_vib,e_cav
      !real*8  e_qmmm_e,e_qmmm_q,e_qmmm_lj,e_mmmm_q,e_mmmm_lj
      real*8 e_lj,e_q,e_spring
      real*8 ehfx,phfx
      
      logical cosmo_on,cosmo1_on,V_APC_on
      real*8  eapc,papc

*     **** external functions ****
      logical control_print
      integer  control_ispin,control_scf_algorithm,control_ks_algorithm
      integer  control_it_in,control_it_out,psi_ne,control_version
      real*8   control_tole,control_tolc,control_ks_alpha
      real*8   rho_error,psi_1energy,psi_error
      real*8   dng_1ehartree,lattice_omega
      real*8   psi_1ke
      real*8   psi_1vl,psi_1v_field,dng_1vl_mm
      real*8   psi_1vnl
      real*8   rho_1exc
      real*8   rho_1pxc
      real*8   ewald_e,ion_ion_e
      real*8   psi_1eorbit
   
      external control_print
      external control_ispin,control_scf_algorithm,control_ks_algorithm
      external control_it_in,control_it_out,psi_ne,control_version
      external control_tole,control_tolc,control_ks_alpha
      external rho_error,psi_1energy,psi_error
      external dng_1ehartree,lattice_omega
      external psi_1ke
      external psi_1vl,psi_1v_field,dng_1vl_mm
      external psi_1vnl
      external rho_1exc
      external rho_1pxc
      external ewald_e,ion_ion_e
      external psi_1eorbit

*     ***** QM/MM external functions ****
      logical  pspw_qmmm_found
      real*8   pspw_qmmm_LJ_E,pspw_qmmm_Q_E,pspw_qmmm_spring_E
      real*8   pspw_qmmm_LJ_Emix,pspw_qmmm_Q_Emix
      external pspw_qmmm_found
      external pspw_qmmm_LJ_E,pspw_qmmm_Q_E,pspw_qmmm_spring_E
      external pspw_qmmm_LJ_Emix,pspw_qmmm_Q_Emix

*     ***** pspw_charge external functions ****
      logical  pspw_charge_found,control_precondition,pspw_HFX
      real*8   pspw_charge_Energy_ion,pspw_charge_Energy_charge
      external pspw_charge_found,control_precondition,pspw_HFX
      external pspw_charge_Energy_ion,pspw_charge_Energy_charge

      real*8   psi_1_noupdate_energy,psi_eigenvalue
      external psi_1_noupdate_energy,psi_eigenvalue
      logical  psp_U_psputerm,meta_found
      external psp_U_psputerm,meta_found
      logical  nwpw_meta_gga_on,ion_disp_on
      external nwpw_meta_gga_on,ion_disp_on
      real*8   psi_1meta_gga_pxc,ion_disp_energy
      external psi_1meta_gga_pxc,ion_disp_energy

      logical  nwpw_cosmo_on,nwpw_cosmo1_on,pspw_V_APC_on
      external nwpw_cosmo_on,nwpw_cosmo1_on,pspw_V_APC_on
      real*8   psi_1vl_cosmo,nwpw_cosmo_EQionq,nwpw_cosmo_Eqq
      external psi_1vl_cosmo,nwpw_cosmo_EQionq,nwpw_cosmo_Eqq

      integer  control_ks_maxit_orb
      external control_ks_maxit_orb

   
      Ein = E(1)
      call Parallel_taskid(taskid)
      oprint = (taskid.eq.MASTER).and.control_print(print_medium)
      cosmo_on  = nwpw_cosmo_on()
      cosmo1_on = nwpw_cosmo1_on()
      V_APC_on  = pspw_V_APC_on()

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
      sd_count          = 0

      call D3dB_n2ft3d(1,n2ft3d)
      call D3dB_n2ft3d_map(1,n2ft3d_map)

*     **** allocate rho_in and rho_out ****
      value = BA_push_get(mt_dbl,2*n2ft3d,
     >                     'vall_in',vall_in(2),vall_in(1))
      value = value.and.
     >        BA_push_get(mt_dbl,2*n2ft3d,
     >                     'vall_out',vall_out(2),vall_out(1))
      value = value.and.
     >        BA_push_get(mt_dbl,2*n2ft3d,
     >                   'vall_junk',vall_junk(2),vall_junk(1))
      value = value.and.
     >        BA_push_get(mt_dbl,2*n2ft3d,
     >                     'rho_in',rho_in(2),rho_in(1))
      if (.not. value) 
     >   call errquit('bybminimize:out of stack memory',0,0)
c      call dcopy(2*n2ft3d,0.0d0,0,dbl_mb(vall_in(1)),1)
c      call dcopy(2*n2ft3d,0.0d0,0,dbl_mb(vall_out(1)),1)
c      call dcopy(2*n2ft3d,0.0d0,0,dbl_mb(vall_junk(1)),1)
c      call dcopy(2*n2ft3d,0.0d0,0,dbl_mb(rho_in(1)),1)
      call Parallel_shared_vector_zero(.false.,2*n2ft3d,
     >                                 dbl_mb(vall_in(1)))
      call Parallel_shared_vector_zero(.false.,2*n2ft3d,
     >                                 dbl_mb(vall_out(1)))
      call Parallel_shared_vector_zero(.false.,2*n2ft3d,
     >                                 dbl_mb(vall_junk(1)))
      call Parallel_shared_vector_zero(.true.,2*n2ft3d,
     >                                 dbl_mb(rho_in(1)))



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

*     **** iniitialize RMM-DIIS ****
      if (control_ks_algorithm().eq.1) call pspw_rmmdiis_init(5)

*     **** iniitialize blocked cg ****
      ks_block = .false.
      if (control_ks_algorithm().eq.-1)  then
         ks_block = .true.
         call linesearch_maxiter_set(control_ks_maxit_orb())
      end if



*     ***** diis loop ****
      it   = 0
 2    it   = it + 1
      
*     **** diaganolize KS matrix ****
      if (ks_block) then
         call psi_KS_block_update(e00,deltae,it,ks_deltae)
      else
         call psi_KS_update(1,
     >                      control_ks_algorithm(),
     >                      precondition,
     >                      ks_deltae)
      end if

c      call psi_KS_update(1,
c     >                   control_ks_algorithm(),
c     >                   precondition,
c     >                   ks_deltae)


      call rho_1to2()
      Eold = Enew
      Enew = psi_1energy()

      deltae = Enew-Eold

      call electron_gen_vall()
      call electron_get_vall(dbl_mb(vall_in(1)))

*     **** compute deltaV ****
c      call dcopy(ispin*n2ft3d_map,
c     >           dbl_mb(vall_in(1)),1,
c     >           dbl_mb(vall_junk(1)),1)
      call Parallel_shared_vector_copy(.true.,ispin*n2ft3d,
     >                                 dbl_mb(vall_in(1)),
     >                                 dbl_mb(vall_junk(1)))
      call DAXPY_OMP(ispin*n2ft3d_map,
     >              (-1.0d0),
     >              dbl_mb(vall_out(1)),1,
     >              dbl_mb(vall_junk(1)),1)
      call D3dB_rr_dot(1,dbl_mb(vall_junk(1)),
     >                   dbl_mb(vall_junk(1)),deltav)
      if (ispin.gt.1) then
         call D3dB_rr_dot(1,dbl_mb(vall_junk(1)+n2ft3d),
     >                      dbl_mb(vall_junk(1)+n2ft3d),e00)
         deltav = deltav + e00
      end if
c      deltav = ddot(ispin*n2ft3d_map,
c     >                 dbl_mb(vall_junk(1)),1,
c     >                 dbl_mb(vall_junk(1)),1)
c      call D3dB_SumAll(deltav)
      deltav = deltav*dV
      


*     **** update vall using density mixing ****
c      if ((it.le.0) .or.
c     >    ((dabs(deltae).lt.1.0d1) .and.
c     >    (deltav       .lt.1.0d1) .and.
c     >    (.not.stalled          ))) then
      if ((it.le.0) .or.
     >    ((deltae.lt.0.0d0) .and.
     >    (.not.stalled          ))) then

         call nwpw_scf_mixing(dbl_mb(vall_in(1)),dbl_mb(vall_out(1)),
     >                        deltae,diis_error)

*     **** bad convergence - try fixed step steepest descent ****
      else

  30     call sdminimize(sd_it)
         sd_count = sd_count + 1
         Eold   = Enew
         Enew   = psi_1energy()
        
c         if ((Enew.gt.Eold).or.(dabs(Enew-Eold).gt.1.0d-1)) go to 30
         if ((Enew.gt.Eold).and.(sd_count.lt.MAX_SD_COUNT)) go to 30

c         call dcopy(ispin*n2ft3d,
c     >              dbl_mb(vall_out(1)),1,
c     >              dbl_mb(vall_junk(1)),1)
         call Parallel_shared_vector_copy(.true.,ispin*n2ft3d,
     >                                    dbl_mb(vall_out(1)),
     >                                    dbl_mb(vall_junk(1)))
         call electron_gen_vall()
         call electron_get_vall(dbl_mb(vall_out(1))) 
         call nwpw_scf_mixing_reset(dbl_mb(vall_out(1)))


         call DAXPY_omp(ispin*n2ft3d,
     >                  (-1.0d0),
     >                  dbl_mb(vall_out(1)),1,
     >                  dbl_mb(vall_junk(1)),1)
         call D3dB_rr_dot(1,dbl_mb(vall_junk(1)),
     >                   dbl_mb(vall_junk(1)),deltav)
         if (ispin.gt.1) then
            call D3dB_rr_dot(1,dbl_mb(vall_junk(1)+n2ft3d),
     >                         dbl_mb(vall_junk(1)+n2ft3d),e00)
            deltav = deltav + e00
         end if
c         deltav = ddot(ispin*n2ft3d,
c     >                 dbl_mb(vall_junk(1)),1,
c     >                 dbl_mb(vall_junk(1)),1)
c         call D3dB_SumAll(deltav)
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
c      stalled = .false.
      if (deltae.gt.0.0d0) stalled = .true.

      precondition = precondition.and.(dabs(deltae).gt.1*tole)

      done = ( (    (dabs(deltae).lt.tole)
     >         .and.(deltae.lt.0.0d0)
     >         .and.(deltac      .lt.tolc))
     >       .or. (it.ge.it_in)
     >       .or. (sd_count.ge.MAX_SD_COUNT))

      if (.not.done) go to 2



*     **** free memory ****
      call nwpw_scf_mixing_end()
      if (control_ks_algorithm().eq.1) call pspw_rmmdiis_end()
      value =           BA_pop_stack(rho_in(2))
      value = value.and.BA_pop_stack(vall_junk(2))
      value = value.and.BA_pop_stack(vall_out(2))
      value = value.and.BA_pop_stack(vall_in(2))
      if (.not. value) 
     >  call errquit('bybminimize: popping stack',1,0)

c      call psi_check()


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
         e_lj     = pspw_qmmm_LJ_E()
         e_q      = pspw_qmmm_Q_E()
         e_spring = pspw_qmmm_spring_E()
         E(1)  = E(1) + e_lj + e_q + e_spring

         E(11) = e_lj
         E(12) = e_q
         E(13) = e_spring

*        **** Eqm-mm energy ***
         E(14) = pspw_qmmm_LJ_Emix()
         E(14) = E(14) + pspw_qmmm_Q_Emix()
         E(14) = E(14) + dng_1vl_mm()

      end if

*     **** COSMO terms ****
      if (cosmo_on) then

         !*** cosmo1 ****
         if (cosmo1_on) then
            E(46) = psi_1vl_cosmo()
            E(47) = nwpw_cosmo_EQionq()
            E(48) = nwpw_cosmo_Eqq()

         !*** cosmo2 ****
         else
            call electron_apc_energies(eapc,papc)
            E(22) = eapc
            E(23) = papc
            E(46) = eapc
            E(47) = nwpw_cosmo_EQionq() !** E(Qion|q)
            E(48) = nwpw_cosmo_Eqq()    !** E(q|q)

            !E(1)  = E(1) + E(22) - E(23) + E(47) + E(48)
            E(1)  = E(1) + E(22) - E(23)
         end if

      else if (V_APC_on) then
         call electron_apc_energies(eapc,papc)
         E(22) = eapc
         E(23) = papc
         E(1)  = E(1) + eapc - papc
      end if


*     **** get pspw_charge  energies ****
      if (pspw_charge_found()) then
         E(19)  = psi_1v_field()
         E(20)  = pspw_charge_Energy_ion()
         E(21)  = pspw_charge_Energy_charge()
         E(1)   = E(1) + E(20) + E(21)
      end if

*     **** HFX terms ****
      if (pspw_HFX()) then
         call electron_HFX_energies(ehfx,phfx)
         E(26) = ehfx
         E(27) = phfx
      end if

*     **** DFT+U terms ****
      if (psp_U_psputerm()) then
         call electron_U_energies(ehfx,phfx)
         E(29) =  ehfx
         E(30) =  phfx
      end if

*     **** Metadynamics potential terms ****
      if (meta_found()) then
         call electron_meta_energies(ehfx,phfx)
         E(31) =  ehfx
         E(32) =  phfx
      end if

*     **** Metadynamics GGA Tau potential term ****
      if (nwpw_meta_gga_on()) then
         E(10) = E(10) - psi_1meta_gga_pxc()
      end if

*     **** Dispersion energy ****
      if (ion_disp_on()) then
         E(33) = ion_disp_energy()
         E(1)  = E(1) + E(33)
      end if


      failed = (sd_count.ge.MAX_SD_COUNT)

      return
      end
 

