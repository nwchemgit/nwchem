c
c     $Id: c_nominimize.F,v 1.1 2004-11-13 02:25:29 bylaska Exp $ 
c

      subroutine c_nominimize(E)
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
      integer neall,npack1,nbrill
      

*     **** external functions ****
      integer  cpsi_ne,brillioun_nbrillioun
      real*8   c_dng_1ehartree
      real*8   cpsi_1ke
      real*8   cpsi_1vl
      real*8   cpsi_1vnl
      real*8   c_rho_1exc
      real*8   c_rho_1pxc
      real*8   c_ewald_e
      real*8   cpsi_1eorbit
   
      external cpsi_ne,brillioun_nbrillioun
      external c_dng_1ehartree
      external cpsi_1ke
      external cpsi_1vl
      external cpsi_1vnl
      external c_rho_1exc
      external c_rho_1pxc
      external c_ewald_e
      external cpsi_1eorbit


      call Cram_max_npack(npack1)
      neall  = cpsi_ne(1)+cpsi_ne(2)
      nbrill = brillioun_nbrillioun()


*     **** allocate G1 ****
      value = MA_alloc_get(mt_dcpl,npack1*neall*nbrill,
     >                     'G1',G1(2),G1(1))
      if (.not. value) call errquit('nominimize:out of heap memory',0,
     &       MA_ERR)


      call cpsi_1get_Tgradient(dcpl_mb(G1(1)),E0)
      Enew = E0

      call cpsi_1toelectron() 

      eion = c_ewald_e()

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


      value = MA_free_heap(G1(2))
      if (.not. value) 
     >  call errquit('c_nominimize:error freeing heap memory',0, MA_ERR)

  
      return
      end
 

