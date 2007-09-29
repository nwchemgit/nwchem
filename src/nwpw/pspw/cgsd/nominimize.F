c
c     $Id: nominimize.F,v 1.10 2007-09-29 00:33:31 bylaska Exp $ 
c
*     ******************************************
*     *                                        *
*     *             nominimize                 *
*     *                                        *
*     ******************************************
      subroutine nominimize(E)
      implicit none
      real*8     E(*)

#include "mafdecls.fh"
#include "errquit.fh"

*     **** local variables ****
      integer G1(2)
      real*8  E0

      real*8     ehartree,eorbit,exc,pxc,eion
      real*8     Enew


      logical value
      integer neall,npack1
      !real*8 e_pol,e_vib,e_cav
      !real*8  e_qmmm_e,e_qmmm_q,e_qmmm_lj,e_mmmm_q,e_mmmm_lj
      real*8  e_lj,e_q,e_spring
      real*8  ehsic,phsic,exsic,pxsic,ehfx,phfx
      

*     **** external functions ****
      integer  psi_neq,control_version
      real*8   dng_1ehartree
      real*8   psi_1ke
      real*8   psi_1vl,psi_1v_field
      real*8   psi_1vnl
      real*8   rho_1exc
      real*8   rho_1pxc
      real*8   ewald_e,ion_ion_e
      real*8   psi_1eorbit
   
      external psi_neq,control_version
      external dng_1ehartree
      external psi_1ke
      external psi_1vl,psi_1v_field
      external psi_1vnl
      external rho_1exc
      external rho_1pxc
      external ewald_e,ion_ion_e
      external psi_1eorbit
      external linesearch

*     ***** QM/MM external functions ****
      logical  pspw_qmmm_found
      real*8   pspw_qmmm_LJ_E,pspw_qmmm_Q_E,pspw_qmmm_spring_E
      external pspw_qmmm_found
      external pspw_qmmm_LJ_E,pspw_qmmm_Q_E,pspw_qmmm_spring_E

*     ***** pspw_charge external functions ****
      logical  pspw_charge_found
      real*8   pspw_charge_Energy_ion,pspw_charge_Energy_charge
      external pspw_charge_found
      external pspw_charge_Energy_ion,pspw_charge_Energy_charge

*     ***** pspw_SIC and pspw_HFX external functions ****
      logical  pspw_SIC,pspw_HFX
      external pspw_SIC,pspw_HFX

*     ***** PAW functions ****
      logical  psp_pawexist
      real*8   psp_kinetic_core,psi_1ke_atom
      external psp_pawexist
      external psp_kinetic_core,psi_1ke_atom

      call Pack_npack(1,npack1)
      neall = psi_neq(1)+psi_neq(2)


*     **** allocate G1 ****
      value = MA_alloc_get(mt_dcpl,npack1*neall,
     >                     'G1',G1(2),G1(1))
      if (.not. value) call errquit('nominimize:out of heap memory',0,
     &       MA_ERR)


      call psi_1get_Tgradient(dcpl_mb(G1(1)),E0)
      Enew = E0

c      call psi_1toelectron() 

      eion = 0.0d0
      if (control_version().eq.3) eion = ewald_e()
      if (control_version().eq.4) eion = ion_ion_e()

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
      end if

*     **** paw energies ****
      if (psp_pawexist()) then
         E(14) = psi_1ke_atom()      !*** kinetic atom ***
         E(15) = psp_kinetic_core()  !*** kinetic core ***
      end if

*     **** get pspw_charge  energies ****
      if (pspw_charge_found()) then
         E(19)  = psi_1v_field()
         E(20)  = pspw_charge_Energy_ion()
         E(21)  = pspw_charge_Energy_charge()
         E(1)   = E(1) + E(20) + E(21)
      end if


*     **** SIC corrections ****
      if (pspw_SIC()) then
         call electron_SIC_energies(ehsic,phsic,exsic,pxsic)
         E(22) = ehsic
         E(23) = exsic
         E(24) = phsic
         E(25) = pxsic
      end if

*     **** HFX terms ****
      if (pspw_HFX()) then
         call electron_HFX_energies(ehfx,phfx)
         E(26) = ehfx
         E(27) = phfx
      end if


      value = MA_free_heap(G1(2))
      if (.not. value) 
     >  call errquit('nominimize:error freeing heap memory',0, MA_ERR)

  
      return
      end
 

