c
c     $Id: c_nominimize.F,v 1.7 2008-11-17 17:25:40 bylaska Exp $ 
c

      subroutine c_nominimize(E)
      implicit none
      real*8     E(*)

#include "mafdecls.fh"
#include "errquit.fh"

*     **** local variables ****
      integer neall,npack1,nbrillq
      integer G1_tag
      real*8  ehartree,eorbit,exc,pxc,eion
      real*8  E0,Enew


*     **** external functions ****
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
      integer taskid

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

      call cpsi_data_dealloc(G1_tag)

      return
      end
 

