*
* $Id: paw_comp_charge.F,v 1.2 2003-01-23 21:25:12 marat Exp $
*


!**************************************************
!
!       Name: paw_proj_init
!
!       Purpose: initializes the paw projectors
!
!       Created:        7/30/2002
!**************************************************
      subroutine paw_comp_charge_init()
      implicit none

#include "mafdecls.fh"
#include "paw_comp_charge_data.fh"
#include "paw_proj.fh"
#include "paw_utils.fh"

        integer  ion_nion,ion_natm,ion_nkatm,ion_katm
        external ion_nion,ion_natm,ion_nkatm,ion_katm

        integer nion
        integer nkatm
        
        !*** get total number of ions in the system ***
        nion = ion_nion()
        
        !*** get tital number of diff kinds of atoms ***
        nkatm = ion_nkatm()

        write(*,*) "compensation charge total number of basis functions", paw_proj_nbasis(1)
        write(*,*) "compensation charge total number of ions", nion, nkatm

      end subroutine
