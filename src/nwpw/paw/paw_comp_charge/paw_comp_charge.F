*
* $Id: paw_comp_charge.F,v 1.4 2003-01-24 19:58:22 marat Exp $
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
#include "paw_params.fh"
#include "paw_basis.fh"

      integer  ion_nion,ion_natm,ion_nkatm,ion_katm
      external ion_nion,ion_natm,ion_nkatm,ion_katm

      logical ok
      integer nion
      integer nkatm
      integer qlm_size
      integer ia
        
      double precision q00_tmp

      !*** get total number of ions in the system ***
      nion = ion_nion()
        
      !*** get total number of diff kinds of atoms ***
      nkatm = ion_nkatm()

      !*** calculate size of comp. charge array ***
      qlm_size = 0
      do ia=1,nkatm
         qlm_size = qlm_size
     >             + paw_proj_nbasis(ia)*ion_natm(ia)
      end do
      write(*,*) "comp charge size", qlm_size

      !*** allocate comp. charge arrays ***
      ok = .TRUE.
      ok = ok .AND. my_alloc(MT_DCPL,qlm_size,"qlm",qlm)
      ok = ok .AND. my_alloc(MT_DBL,nkatm,"q00",q00)

      !*** set q00 array ***
      do ia=1,nkatm
	q00_tmp = (get_core_charge(ia)-get_ion_charge(ia))/sqrt(fourpi)    
        dbl_mb(q00(1)+ia-1)=q00_tmp
        write(*,*) "comp_charge:4pi",sqrt(fourpi)
        write(*,*) "comp_charge:q00",q00_tmp,dbl_mb(q00(1)+ia-1)
        write(*,*) "comp_charge:ion",get_ion_charge(ia)
        write(*,*) "comp_charge:core",get_core_charge(ia)
      end do
      end subroutine
