*
* $Id: paw_comp_charge.F,v 1.12 2003-02-11 01:41:23 edo Exp $
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
#include "paw_geom.fh"
#include "paw_ma.fh"
#include "paw_params.fh"
#include "paw_basis.fh"


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
     >             + (paw_basis_mult_l(ia) + 1)**2*ion_natm(ia)
      end do

      !*** allocate comp. charge arrays ***
      ok = .TRUE.
      ok = ok .AND. my_alloc(MT_DCPL,qlm_size,"qlm",qlm)
      ok = ok .AND. my_alloc(MT_DBL,nkatm,"q00",q00)

      !*** set q00 array ***
      do ia=1,nkatm
        q00_tmp = (paw_basis_core_charge(ia)-paw_basis_ion_charge(ia))
     >            /sqrt(fourpi)    
        dbl_mb(q00(1)+ia-1)=q00_tmp
      end do

      end

      subroutine paw_comp_charge_update()
        implicit none

#include "mafdecls.fh"
#include "paw_comp_charge_data.fh"
#include "paw_proj.fh"
#include "paw_basis.fh"
#include "paw_geom.fh"
#include "paw_comp_charge_matrix.fh"
#include "paw_ovlp_data.fh"

        integer  psi_ne
        external psi_ne

        integer na
        integer ia
        integer nlmi
        integer i1,i2
        integer l,l1,l2
        integer m
        integer tot_natom
        integer nproj
        integer i_matr,i_matr0
        integer i_w,i_w0,i_w1
        integer m1,m2
        integer nb2,nb
        integer mult_l
        complex qlm_tmp
        
        double precision gen_gaunt_coeff
        external gen_gaunt_coeff

        tot_natom = ion_nion()
    
c       !*** set values for l=0 m=0 elements ***
        nlmi = 0
        do na=1,tot_natom
          ia = ion_katm(na)
          dcpl_mb(qlm(1)+nlmi) = dbl_mb(q00(1)+ia-1)
          nlmi = nlmi + (paw_basis_mult_l(ia)+1)**2
        end do
 
        
        nlmi = 0
        do na=1,tot_natom
          ia = ion_katm(na)
          mult_l = paw_basis_mult_l(ia)
          i_matr0 = int_mb(i_comp_charge_matrix(1) + ia - 1)
          i_w0 = int_mb(i_paw_ovlp_w(1) + ia - 1)
          nb =paw_basis_nbasis(ia)
          nb2 = nb*nb
          nproj = paw_proj_nbasis(ia)
          do l=0,mult_l  
            do m=0,l
              qlm_tmp = cmplx(0.0,0.0)
              i_w1 = 0
              do i1=1,nb
                l1 = paw_basis_orb_l(i1,ia)
                do m1=-l1,l1
                  do i2=1,nb
                    l2 = paw_basis_orb_l(i2,ia)
                    do m2=-l2,l2
                      i_matr = i_matr0 +i1+(i2-1)*nb+l*nb2
                      i_w = i_w0 + i_w1
                      qlm_tmp = qlm_tmp + 
     >                     dbl_mb(comp_charge_matrix(1) + i_matr-1)* 
     >                     gen_gaunt_coeff(l,m,l2,m2,l1,m1)*
     >                     dcpl_mb(paw_ovlp_w(1) + i_w) 
                      i_w1 = i_w1 + 1
                     end do
                   end do
                 end do
               end do
               dcpl_mb(qlm(1) + nlmi) = dcpl_mb(qlm(1) + nlmi)+ qlm_tmp
               nlmi = nlmi+1
             end do
           end do
         end do 
 
        write(*,*) "paw_comp_charge_coeff"
        nlmi = 0
        do na=1,tot_natom
          ia = ion_katm(na)
          mult_l = paw_basis_mult_l(ia)
          do l=0,mult_l  
          do m=-l,l
            write(*,*) l,m
            write(*,*) dcpl_mb(qlm(1) + nlmi) 
            nlmi = nlmi+1
          end do
          end do
         end do 
      end !subroutine



      subroutine paw_comp_charge_end()
      implicit none

#include "paw_comp_charge_data.fh"
#include "paw_ma.fh"

      logical ok

      ok =        my_dealloc(qlm)
      ok = ok.and.my_dealloc(q00)
      if (.not.ok)
     > call errquit("paw_comp_charge_init: error freeing heap",0)

      return
      end !subroutine




