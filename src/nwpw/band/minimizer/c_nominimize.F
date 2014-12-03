c
c     $Id$ 
c

      subroutine c_nominimize(E)
      implicit none
      real*8     E(*)

#include "bafdecls.fh"
#include "errquit.fh"

*     **** local variables ****
      integer neall,npack1,nbrillq,taskid
      integer G1_tag
      real*8  ehartree,eorbit,exc,pxc,eion,ehfx,phfx
      real*8  E0,Enew


*     **** external functions ****
      logical  band_HFX
      integer  cpsi_neq,Pneb_nbrillq,cpsi_data_alloc
      real*8   c_dng_1ehartree
      real*8   cpsi_1ke
      real*8   cpsi_1vl
      real*8   cpsi_1vnl
      real*8   c_rho_1exc
      real*8   c_rho_1pxc
      real*8   ewald_e
      real*8   cpsi_1eorbit,cpsi_smearcorrection
      real*8   cpsi_1vnlso
      external band_HFX
      external cpsi_neq,Pneb_nbrillq,cpsi_data_alloc
      external c_dng_1ehartree
      external cpsi_1ke
      external cpsi_1vl
      external cpsi_1vnl
      external c_rho_1exc
      external c_rho_1pxc
      external ewald_e
      external cpsi_1eorbit,cpsi_smearcorrection
      external cpsi_1vnlso
      logical  ion_disp_on
      external ion_disp_on
      real*8   ion_disp_energy
      external ion_disp_energy


      call Cram_max_npack(npack1)
      neall  = cpsi_neq(1)+cpsi_neq(2)
      nbrillq = Pneb_nbrillq()


*     **** allocate G1 ****
      G1_tag = cpsi_data_alloc(nbrillq,neall,2*npack1)

      call cpsi_1get_Tgradient(G1_tag,E0)
      Enew = E0

c      call cpsi_1toelectron() 

      eion = ewald_e()

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

c     **** HFX terms ****
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

      call cpsi_data_dealloc(G1_tag)

      return
      end
 

