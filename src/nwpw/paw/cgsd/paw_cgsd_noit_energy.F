c     $Id: paw_cgsd_noit_energy.F,v 1.1 2004-11-29 16:05:19 bylaska Exp $                       


*     ***************************
*     *				*
*     *	   paw_cgsd_noit_energy	*
*     *				*
*     ***************************

      real*8 function paw_cgsd_noit_energy()
      implicit none

#include "mafdecls.fh"
#include "rtdb.fh"
#include "errquit.fh"

*     **** local variables ***
      logical value
      real*8  E(30)
      integer rtdb

*     **** external functions ****
      integer  control_rtdb
      external control_rtdb

*     **** generate phaze factors and local psp and core density ****
      call paw_electron_gen_vl_potential()


      call paw_nominimize(E)


*     **** debug - energies written to rtdb for numerical stress ****
      rtdb = control_rtdb()
      value =           rtdb_put(rtdb,'paw:E_ke',      mt_dbl,1,
     >                           E(2)+E(3))
      value = value.and.rtdb_put(rtdb,'paw:E_hartree', mt_dbl,1,
     >                           E(4)+E(5))
      value = value.and.rtdb_put(rtdb,'paw:E_mult',    mt_dbl,1,
     >                           E(6))
      value = value.and.rtdb_put(rtdb,'paw:E_xc',      mt_dbl,1,
     >                           E(7)+E(8))
      value = value.and.rtdb_put(rtdb,'paw:E_local',   mt_dbl,1,
     >                           E(9)+E(10))
      if (.not. value) call errquit(
     >   'paw_cgsd_noit_energy: numerical stress - error writing rtdb',
     >    0,RTDB_ERR)


      paw_cgsd_noit_energy = E(1)
      return
      end

