*
* $Id: paw_comp_charge.F,v 1.6 2003-01-27 06:40:12 marat Exp $
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
#include "paw_ma.fh"
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
      end do

      end subroutine

      subroutine paw_comp_charge_update()
        implicit none

#include "mafdecls.fh"
#include "paw_comp_charge_data.fh"
#include "paw_proj.fh"
#include "paw_basis.fh"
  
        integer  ion_nion,ion_katm,psi_ne
        external ion_nion,ion_katm,psi_ne

        integer na
        integer ia
        integer i00
        integer i,i1,i2,nb
        integer l,l1,l2
        integer m
        integer tot_natom

        write(*,*) "paw_comp_charge_update"
        tot_natom = ion_nion()
    
c       !*** set values for l=0 m=0 elements ***
        i00 = 0
        do na=1,tot_natom
          ia = ion_katm(na)
          dcpl_mb(qlm(1)+i00) = dbl_mb(q00(1)+ia-1)
          i00 = i00 + paw_proj_nbasis(ia)
        end do
 
c        do ia=0,qlm(3)-1
c         write(*,*) dcpl_mb(qlm(1)+ia)
c        end do
        
        do na=1,tot_natom
          ia = ion_katm(na)
          do i=1,paw_proj_nbasis(ia)
            l = paw_proj_l(i,ia)
            m = paw_proj_m(i,ia)
            nb =paw_basis_nbasis(ia)
            do i1=1,nb
            do i2=1,nb
              l1 = paw_basis_orb_l(i1,ia)
              l2 = paw_basis_orb_l(i2,ia)
              write(*,*) "i1,i2,l1,l2",i1,i2,l1,l2
            end do 
            end do 
          end do
        end do
	end subroutine
	

