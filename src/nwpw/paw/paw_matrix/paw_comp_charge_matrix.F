*
* $Id: paw_comp_charge_matrix.F,v 1.2 2003-02-04 01:33:42 marat Exp $
*
!**************************************************
!
!	name:
!
!	purpose 
!
!	created:	7/9/2002
!**************************************************
      subroutine  init_paw_comp_charge_matrix()
      implicit none

#include "mafdecls.fh"
#include "paw_basis.fh"
#include "paw_comp_charge_matrix.fh"
#include "paw_utils.fh"
#include "paw_ma.fh"


      !*** local variables ***
      logical ok
      integer it,tot_ntype
      integer mult_l
      integer nbasis
      integer matr_size, matr_index

      tot_ntype = paw_basis_tot_ntype()
      ok = my_alloc(mt_int,tot_ntype,"i_comp_charge_matrix",
     >              i_comp_charge_matrix)
      ok = my_alloc(mt_int,tot_ntype,"s_comp_charge_matrix",
     >              i_comp_charge_matrix)
  
      !*** determine matrix size ***
      matr_size = 0
      do it=1,tot_ntype
        nbasis = paw_basis_nbasis(it)
        mult_l  = paw_basis_mult_l(it)
        int_mb(s_comp_charge_matrix(1) + it - 1) = nbasis**2*(mult_l+1)                  
        matr_size = matr_size + int_mb(s_comp_charge_matrix(1) + it - 1)                
      end do
  
      !*** set the indexing to access diff atom types ***
      int_mb(i_comp_charge_matrix(1) ) = 0
      do it=1,tot_ntype-1
        int_mb(i_comp_charge_matrix(1) +it) = 
     >         int_mb(i_comp_charge_matrix(1) + it - 1)  
     >       + int_mb(s_comp_charge_matrix(1) + it - 1)
      end do

      ok = my_alloc(mt_dbl,matr_size,"comp_charge_matrix",
     >              comp_charge_matrix)
  
      do it=1,tot_ntype
    
        matr_index = int_mb(i_comp_charge_matrix(1) +it-1)
        call find_comp_charge_matrix( 
     >            int_mb(paw_basis_i_i_cut(it)),
     >            int_mb(paw_basis_i_ngrid(it)),
     >            int_mb(paw_basis_i_nbasis(it)),
     >            int_mb(paw_basis_i_orb_l(it)),
     >            paw_basis_mult_l(it),
     >            dbl_mb(paw_basis_i_log_amesh(it)),
     >            dbl_mb(paw_basis_i_rgrid(it)),
     >            dbl_mb(paw_basis_i_phi_ae(it)),
     >            dbl_mb(paw_basis_i_phi_ps(it)),
     >            dbl_mb(comp_charge_matrix(1) + matr_index))
      end do

      return
      end !subroutine init_paw_comp_charge_matrix

!*************************************************
!
!   name :   find_comp_charge_matrix
!
!
!   purpose :
!
!
!   created :
!
!*************************************************
      subroutine find_comp_charge_matrix( ic,
     >                        ngrid,
     >                        nbasis,
     >                        orb_l,
     >                        mult_l,
     >                        log_amesh,
     >                        r,
     >                        phi_ae,
     >                        phi_ps,
     >                        matrix)
      implicit none
      integer    ic
      integer    ngrid
      integer    nbasis
      integer    orb_l(nbasis)
      integer    mult_l

      double precision log_amesh
      double precision  r(ngrid)
      double precision phi_ae(ngrid,nbasis)
      double precision phi_ps(ngrid,nbasis)
      double precision matrix(nbasis,nbasis,mult_l+1)

#include "mafdecls.fh"
#include "integrate.fh"

      !*** local variables ***
      logical  ok
      integer  i,j,ii
      integer  la
      integer  power_r
      integer  i_f,l_f
      double precision tmp


      power_r = 0

      !*** temporary array for calculatcomp_charges ***
      ok = ma_push_get(mt_dbl,ic,"tmp_array",l_f,i_f)

      do i=1,nbasis
        do j=1,i

          do ii=1,ic
          dbl_mb(i_f+ii-1) = phi_ae(ii,i)*phi_ae(ii,j) 
     >                         - phi_ps(ii,i)*phi_ps(ii,j)
          end do

          do la=0,mult_l

            tmp = def_integr(2*la+2,dbl_mb(i_f),la,r,log_amesh,ic)

            matrix(i,j,la+1) = tmp
            matrix(j,i,la+1) = tmp

          end do

        end do
      end do

      ok = ma_chop_stack(l_f)

      return
      end !subroutine find_comp_charge_matrix

