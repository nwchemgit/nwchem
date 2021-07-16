*
*     $Id$                       
*
*  ************************************************************
*  *                MPI cgminimize2 routine                    *
*  *               (Fletcher-Reeves' steps)                   *
*  *                                                          *
*  *  This is a developing cgsdv3 parallel code wrtten for    *
*  *  NWChem                                                  *
*  *                                                          *
*  ************************************************************
  
      subroutine c_cgminimize2(E,deltae,deltac,current_iteration,it_in)
      implicit none
      real*8     E(*)
      real*8     deltae,deltac
      integer    current_iteration
      integer    it_in

#include "bafdecls.fh"
#include "errquit.fh"

*     **** local variables ****

      real*8  deltat_min
      parameter (deltat_min=1.0d-3)
       
c      complex*16 H0(nfft3d,nemax)
c      complex*16 G1(nfft3d,nemax)
      integer H0_tag,G0_tag,G1_tag
      real*8  E0,dE0

      real*8     sum0,sum1,scale,tole,tolc
      real*8     ehartree,eorbit,exc,pxc,eion
      real*8     Enew,Eold,Estart
      common / c_cgsd_block / Enew,Eold,Estart

      logical value,failed
      integer it
      real*8 tmin,deltat
      real*8 max_sigma,ehfx,phfx

      integer ispin,ne(2)
      integer neall,nfft3d,npack1,nbrill

*     **** external functions ****
      real*8   cpsi_geodesic2_energy,cpsi_geodesic2_denergy
      external cpsi_geodesic2_energy,cpsi_geodesic2_denergy

      logical  band_hfx
      integer  cpsi_neq,control_version,cpsi_ispinq,cpsi_nbrillq
      integer  cpsi_data_push_stack
      real*8   control_tole,control_tolc
      real*8   c_rho_error
      real*8   c_dng_1ehartree
      real*8   cpsi_1ke
      real*8   cpsi_1vl
      real*8   cpsi_1vnl
      real*8   cpsi_1vnlso
      real*8   c_rho_1exc
      real*8   c_rho_1pxc
      real*8   ewald_e
      real*8   cpsi_1eorbit
      real*8   linesearch
      external band_hfx
      external cpsi_neq,control_version,cpsi_ispinq,cpsi_nbrillq
      external cpsi_data_push_stack
      external control_tole,control_tolc
      external c_rho_error
      external c_dng_1ehartree
      external cpsi_1ke
      external cpsi_1vl
      external cpsi_1vnl
      external c_rho_1exc
      external c_rho_1pxc
      external ewald_e
      external cpsi_1eorbit
      external linesearch
      real*8   ion_disp_energy,cpsi_smearcorrection
      external ion_disp_energy,cpsi_smearcorrection
      logical  ion_disp_on,control_fractional
      external ion_disp_on,control_fractional
      real*8   cpsi_1energy
      external cpsi_1energy
       
  
      call C3dB_nfft3d(1,nfft3d)
      call Cram_max_npack(npack1)
      ispin  = cpsi_ispinq()
      ne(1)  = cpsi_neq(1)
      ne(2)  = cpsi_neq(2)
      neall  = ne(1)+ne(2)
      nbrill = cpsi_nbrillq()

      value = MA_set_error_print(.true.)

*     **** allocate H0, G0, and G1 ****
      H0_tag = cpsi_data_push_stack(nbrill,neall,2*npack1)
c      G0_tag = cpsi_data_push_stack(nbrill,neall,2*npack1)
      G1_tag = cpsi_data_push_stack(nbrill,neall,2*npack1)


      Estart = Enew
      Eold   = Enew

*     ***** get the initial gradient and direction ****
      call cpsi_1get_TSgradient(G1_tag,E0)

      call Pneb_ff_traceall_tag(0,0,G1_tag,G1_tag,sum1)
      call BGrsm_ff_Copy_tag(G1_tag,H0_tag)


*     ******************************************
*     ****                                  ****
*     **** Start of conjugate gradient loop ****
*     ****                                  ****
*     ******************************************

      tole  = control_tole()
      tolc  = control_tolc()
      tmin  = deltat_min

      do it=2,it_in

*        **** initialize the geodesic line data structure ****
         call cpsi_1geodesic2_start(H0_tag,max_sigma,dE0)

*        ******* line search *********
         if (tmin.gt.0.0d0) then
            deltat = tmin
         else
            deltat = deltat_min
         end if

         Enew = linesearch(0.0d0,E0,dE0,deltat,
     >                        cpsi_geodesic2_energy,
     >                        cpsi_geodesic2_denergy,
     >                        tole,tmin,deltae,1)
c         write(*,*)'e0= ',e0
c         write(*,*)'de0 = ',de0
c         write(*,*)'tmin = ',tmin
c         call flush(6)
         call cpsi_geodesic2_final(tmin)

         deltac = c_rho_error()

*        **** exit loop early ****
         if ((dabs(deltae).lt.tole).and.(deltac.lt.tolc)) 
     >      go to 30
     
   

*        **** transport the previous search directions ****
c         call BGrsm_ff_Copy_tag(G1_tag,G0_tag)

         call cpsi_1geodesic2_transport( tmin,H0_tag)
c         call cpsi_1geodesic2_Gtransport(tmin,G0_tag)

*        **** make psi1 <--- psi2(tmin) ****
         call cpsi_2to1()

*        **** get the new gradient - also updates densities****
         call cpsi_1get_TSgradient(G1_tag,E0)

c         call BGrsm_fff_Sub_tag(G1_tag,G0_tag,G0_tag)

         sum0  = sum1

         call Pneb_ff_traceall_tag(0,0,G1_tag,G1_tag,sum1)

*        **** the new direction using Fletcher-Reeves ****
         if ( (dabs(deltae).le.(1.0d-2)).and.
     >        (tmin.gt.deltat_min)) then

            if (sum0.gt.1.0d-15) then
               scale = sum1/sum0
            else
               scale = 0.0d0
            end if
            if (scale.gt.10.0d0) scale = 0.0d0

            call BGrsm_f_dscale_tag(scale,H0_tag)
            call BGrsm_fff_Sum_tag(G1_tag,H0_tag,H0_tag)

*        **** the new direction using steepest-descent ****
         else
            call BGrsm_ff_Copy_tag(G1_tag,H0_tag)
         end if

      end do

*     **** initialize the geodesic line data structure ****
      call cpsi_1geodesic2_start(H0_tag,max_sigma,dE0)

*     ******* line search *********
      if (tmin.gt.0.0d0) then
         deltat = tmin
      else
         deltat = deltat_min
      end if

      Enew = linesearch(0.0d0,E0,dE0,deltat,
     >                        cpsi_geodesic2_energy,
     >                        cpsi_geodesic2_denergy,
     >                        tole,tmin,deltae,1)

      call cpsi_geodesic2_final(tmin)

      deltac = c_rho_error()

*     **** free memory ****
  30  call cpsi_data_pop_stack(G1_tag)
c      call cpsi_data_pop_stack(G0_tag)
      call cpsi_data_pop_stack(H0_tag)

      call cpsi_2to1()
      call cpsi_1toelectron() 

      eion     = ewald_e()
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

      if (control_fractional()) then
        call cpsi_1gen_hml()
        call cpsi_diagonalize_hml()
        call cpsi_1rotate2()
        call cpsi_2to1()
        call cpsi_1define_occupation(-1.0d0,.false.)
        !call cpsi_1define_occupation(0.0d0,.true.)
        E(12) = cpsi_smearcorrection()
        E(1)  = E(1) + E(12)
      end if

      return
      end
 

