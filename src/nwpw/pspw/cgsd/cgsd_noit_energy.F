c     $Id: cgsd_noit_energy.F,v 1.1 2002-04-25 00:34:38 bylaska Exp $                       
*     ***************************
*     *				*
*     *	   cgsd_noit_energy	*
*     *				*
*     ***************************

      real*8 function cgsd_noit_energy()
      implicit none
#include "stdio.fh"

#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"

      integer MASTER
      parameter (MASTER=0)

      logical value
      real*8  E(30)
      integer rtdb

*     **** external functions ****
      logical  psp_semicore
      integer  control_rtdb
      external psp_semicore
      external control_rtdb

*     **** generate phaze factors and local psp and core density ****
      call phafac()
      call electron_gen_vl_potential()
      if (psp_semicore(0)) call semicore_density_update()


      call nominimize(E)

*     **** diagonalize hamiltonian and rotate psi ****
c     call psi_spin_density(en)
c     call psi_1gen_hml()
c     call psi_diagonalize_hml()
c
c     call psi_1rotate2()
c     call psi_2to1()


*     **** debug - energies written to rtdb for numerical stress ****
      rtdb = control_rtdb()
      value =           rtdb_put(rtdb,'pspw:E_ke',      mt_dbl,1,E(6))
      value = value.and.rtdb_put(rtdb,'pspw:E_hartree', mt_dbl,1,E(3))
      value = value.and.rtdb_put(rtdb,'pspw:E_xc',      mt_dbl,1,E(4))
      value = value.and.rtdb_put(rtdb,'pspw:E_ewald',   mt_dbl,1,E(5))
      value = value.and.rtdb_put(rtdb,'pspw:E_local',   mt_dbl,1,E(7))
      value = value.and.rtdb_put(rtdb,'pspw:E_nonlocal',mt_dbl,1,E(8))
      if (.not. value) call errquit(
     >   'cgsd_noit_energy: numerical stress - error writing rtdb',0)


      cgsd_noit_energy = E(1)
      return
      end

