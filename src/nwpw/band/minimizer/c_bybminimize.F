c
c     $Id: c_bybminimize.F,v 1.1 2003-12-02 19:16:49 bylaska Exp $                       
c

*  ************************************************************
*  *                                                          *
*  *             Band by Band Kohn-Sham Minimizer             *
*  *                                                          *
*  *                                                          *
*  *                                                          *
*  ************************************************************
  
      subroutine c_bybminimize(E,deltae,deltac,
     >                       set_iterations,iterations)
      implicit none
      real*8     E(*)
      real*8     deltae,deltac
      logical    set_iterations
      integer    iterations


#include "stdio.fh"
#include "mafdecls.fh"

*     **** local variables ****
      integer MASTER,taskid
      parameter (MASTER=0)

      real*8  deltat_min
      parameter (deltat_min=1.0d-3)
       
      integer vall_in(2),vall_out(2),vall_junk(2)
      real*8  E0,dE0,deltae_old,Ein,deltae_history(10)
      real*8  ks_deltae,deltav,dV,deltav_old,diis_error
      integer nx,ny,nz,stalled_count


      real*8     sum0,sum1,sum3,scale,tole,tolc
      real*8     ehartree,eorbit,exc,pxc,eion
      real*8     Enew,Eold,alpha

      integer it,it_in,i,ispin,bfgscount,icount,sd_it,cg_it

      logical value,done,stalled,deltav_bad(4)
      integer nfft3d
      real*8  e_ionmm,e_qmmm,e_mmmm,e_pol,e_vib,e_cav
      real*8  e_qmmm_e,e_qmmm_q,e_qmmm_lj,e_mmmm_q,e_mmmm_lj
      


*     **** external functions ****
      integer  control_ispin
      integer  control_it_in,control_it_out,cpsi_ne,control_version
      real*8   control_tole,control_tolc,control_ks_alpha
      real*8   c_rho_error,cpsi_1energy
      real*8   c_dng_1ehartree,lattice_omega
      real*8   cpsi_1ke
      real*8   cpsi_1vl
      real*8   cpsi_1vnl
      real*8   c_rho_1exc
      real*8   c_rho_1pxc
      real*8   c_ewald_e
      real*8   cpsi_1eorbit
      real*8   ddot
   
      external control_ispin
      external control_it_in,control_it_out,cpsi_ne,control_version
      external control_tole,control_tolc,control_ks_alpha
      external c_rho_error,cpsi_1energy
      external c_dng_1ehartree,lattice_omega
      external cpsi_1ke
      external cpsi_1vl
      external cpsi_1vnl
      external c_rho_1exc
      external c_rho_1pxc
      external c_ewald_e
      external cpsi_1eorbit
      external ddot

   
      Ein = E(1)
      call Parallel_taskid(taskid)

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
      ispin = control_ispin()
      deltav_old = 10.0d0
      deltav     = 0.0d0

      stalled       = .false.
      deltae_history(1) = 0.0d0
      deltae_history(2) = 0.0d0
      deltae_history(3) = 0.0d0
      deltae_history(4) = 0.0d0
      stalled_count     = 0

      call C3dB_nfft3d(1,nfft3d)

*     **** allocate vall_in and vall_out ****
      value = MA_push_get(mt_dcpl,2*nfft3d,
     >                     'vall_in',vall_in(2),vall_in(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,2*nfft3d,
     >                     'vall_out',vall_out(2),vall_out(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,2*nfft3d,
     >                   'vall_junk',vall_junk(2),vall_junk(1))
      if (.not. value) 
     >   call errquit('c_bybminimize:out of stack memory',0,0)

*     **** ion-ion energy ****
      eion = c_ewald_e()


*     **********************
*     **** bybminimizer ****
*     **********************


*     **** set the initial density ****
      Enew   = cpsi_1energy()
      alpha = control_ks_alpha()
      deltae = -9232323299.0d0
      ks_deltae = tole
      call c_electron_gen_vall()
      call c_electron_get_vall(dbl_mb(vall_out(1)))
      call c_electron_get_vall(dbl_mb(vall_in(1))) 
      

      call cpsi_1gen_hml()
      call cpsi_diagonalize_hml()
      call cpsi_1rotate2()
      call cpsi_2to1()

  
*     **** iniitialize DIIS ****    
      call nwpw_diis_init(10,ispin*2*nfft3d,dcpl_mb(vall_out(1)))


*     ***** diis loop ****
      it   = 0
 2    it   = it + 1
      
*     **** diaganolize KS matrix ****
      call cpsi_KS_update(1,ks_deltae)

      call c_rho_1to2()
      Eold = Enew
      Enew = cpsi_1energy()

      deltae = Enew-Eold

      call dcopy(ispin*2*nfft3d,
     >           dbl_mb(vall_in(1)),1,
     >           dbl_mb(vall_junk(1)),1)
      call c_electron_gen_vall()
      call c_electron_get_vall(dbl_mb(vall_in(1)))
      call daxpy(ispin*2*nfft3d,
     >              (-1.0d0),
     >              dbl_mb(vall_in(1)),1,
     >              dbl_mb(vall_junk(1)),1)
      deltav = ddot(ispin*2*nfft3d,
     >                 dbl_mb(vall_junk(1)),1,
     >                 dbl_mb(vall_junk(1)),1)
      call C3dB_SumAll(deltav)
      deltav = deltav*dV


*     **** update vall ****
      if ((it.le.0) .or.
     >    ((dabs(deltae).lt.1.0d+1) .and.
     >    (deltav       .lt.1.0d+1) .and.
     >    (.not.stalled          ))) then
        call nwpw_diis(dbl_mb(vall_in(1)),
     >                 dbl_mb(vall_in(1)),
     >                 dbl_mb(vall_out(1)),
     >                 diis_error)

      else if ((dabs(deltae).lt.1.0d+2) .and.
     >         (deltav      .lt.1.0d+10) .and.
     >         (.not.stalled          )) then
         call nwpw_diis_reset(dbl_mb(vall_in(1)))
         call dscal(2*nfft3d,(1.0d0-alpha),dcpl_mb(vall_out(1)),1)
         call daxpy(2*nfft3d,alpha,dcpl_mb(vall_in(1)),1,
     >                             dcpl_mb(vall_out(1)),1)
        
      else
         call c_sdminimize(sd_it)
  30     Eold   = Enew
         Enew   = cpsi_1energy()

         call dcopy(ispin*2*nfft3d,
     >              dcpl_mb(vall_out(1)),1,
     >              dcpl_mb(vall_junk(1)),1)

         call c_electron_gen_vall()
         call c_electron_get_vall(dcpl_mb(vall_out(1))) 

         call nwpw_diis_reset(dcpl_mb(vall_out(1)))

         call daxpy(ispin*2*nfft3d,
     >              (-1.0d0),
     >              dcpl_mb(vall_out(1)),1,
     >              dcpl_mb(vall_junk(1)),1)
         deltav = ddot(ispin*2*nfft3d,
     >                 dcpl_mb(vall_junk(1)),1,
     >                 dcpl_mb(vall_junk(1)),1)
         call C3dB_SumAll(deltav)
         deltav = deltav*dV

         call cpsi_1gen_hml()
         call cpsi_diagonalize_hml_assending()
         call cpsi_1rotate2()
         call cpsi_2to1()

         stalled       = .false.
         deltae_history(1) = 0.0d0
         deltae_history(2) = 0.0d0
         deltae_history(3) = 0.0d0
         deltae_history(4) = 0.0d0
         stalled_count     = 0

      end if
      call c_electron_set_vall(dcpl_mb(vall_out(1))) 
          
*     **** tolerance checks ****
      deltae = Enew-Eold   
      deltac = c_rho_error()
      E(1)   = Enew+eion

      if ((taskid.eq.MASTER).and.(.not.set_iterations)) then 
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


      done = ( (    (dabs(deltae).lt.tole)
     >         .and.(deltae.lt.0.0d0)
     >         .and.(deltac      .lt.tolc))
     >       .or. (it.ge.it_in))

      if (.not.done) go to 2



*     **** free memory ****
      call nwpw_list_end()
      value =           MA_pop_stack(vall_junk(2))
      value = value.and.MA_pop_stack(vall_out(2))
      value = value.and.MA_pop_stack(vall_in(2))
      if (.not. value) 
     >  call errquit('c_bybminimize: popping stack',1,0)

      call cpsi_check()


*     **** set return energies **** - This is duplicate code
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




      return
      end
 

