*
* $Id$
*

!**************************************************
!
!       name:
!
!       purpose 
!
!       created:        7/9/2002
!**************************************************
      subroutine  init_paw_comp_pot_matrix()
      implicit none

#include "bafdecls.fh"
#include "paw_basis.fh"
#include "paw_matrix_comp_pot.fh"
#include "paw_utils.fh"
#include "paw_ma.fh"


      !*** local variables ***
      logical ok
      integer it,tot_ntype
      integer mult_l
      integer nbasis
      integer mtr_size, matr_index

      tot_ntype = paw_basis_tot_ntype()
      ok = my_alloc(mt_int,tot_ntype,"i_paw_comp_pot_matrix",
     >              i_paw_comp_pot_matrix)
      ok = ok.and.
     >     my_alloc(mt_int,tot_ntype,"s_paw_comp_pot_matrix",
     >              s_paw_comp_pot_matrix)
      if (.not.ok)
     > call errquit(
     > "init_paw_paw_comp_pot_matrix: out of heap memory",0,1)
  

      !*** determine matrix size ***
      mtr_size = 0
      do it=1,tot_ntype
        nbasis = paw_basis_nbasis(it)
        mult_l  = paw_basis_mult_l(it)
        int_mb(s_paw_comp_pot_matrix(1)+it-1) = nbasis**2*(mult_l+1)
        mtr_size=mtr_size+int_mb(s_paw_comp_pot_matrix(1) + it - 1)
      end do
  
      !*** set the indexing to access diff atom types ***
      int_mb(i_paw_comp_pot_matrix(1) )  = 0
      do it=1,tot_ntype-1
        int_mb(i_paw_comp_pot_matrix(1) +it) = 
     >           int_mb(i_paw_comp_pot_matrix(1) + it - 1)
     >         + int_mb(s_paw_comp_pot_matrix(1) + it - 1)
      end do

      ok = my_alloc(mt_dbl,mtr_size,"paw_comp_pot_matrix",
     >                paw_comp_pot_matrix)
      if (.not.ok)
     > call errquit(
     > "init_paw_comp_pot_matrix: out of heap memory",0,1)
  
      do it=1,tot_ntype
    
        matr_index = int_mb(i_paw_comp_pot_matrix(1) +it-1)
        call find_comp_pot_matrix(
     >            int_mb(paw_basis_i_i_cut(it)),
     >            int_mb(paw_basis_i_ngrid(it)),
     >            int_mb(paw_basis_i_nbasis(it)),
     >            int_mb(paw_basis_i_orb_l(it)),
     >            paw_basis_mult_l(it),
     >            dbl_mb(paw_basis_i_log_amesh(it)),
     >            dbl_mb(paw_basis_i_rgrid(it)),
     >            dbl_mb(paw_basis_i_sigma(it)),
     >            dbl_mb(paw_basis_i_phi_ps(it)),
     >            dbl_mb(paw_comp_pot_matrix(1) + matr_index))
      end do

      return
      end !subroutine init_paw_comp_pot_matrix

      subroutine  end_paw_comp_pot_matrix()
      implicit none

#include "bafdecls.fh"
#include "paw_basis.fh"
#include "paw_matrix_comp_pot.fh"
#include "paw_utils.fh"
#include "paw_ma.fh"


      !*** local variables ***
      logical ok


      ok =        my_dealloc(i_paw_comp_pot_matrix)
      ok = ok.and.my_dealloc(s_paw_comp_pot_matrix)
      ok = ok.and.my_dealloc(paw_comp_pot_matrix)
      if (.not.ok)
     > call errquit(
     > "end_paw_comp_pot_matrix: error freeing heap",0,0)

      return
      end !subroutine end_paw_comp_pot_matrix

!*************************************************
!
!   name :   find_comp_pot_matrix
!
!
!   purpose :
!
!
!   created :
!
!*************************************************
      subroutine find_comp_pot_matrix( ic,
     >                        ngrid,
     >                        nbasis,
     >                        orb_l,
     >                        mult_l,
     >                        log_amesh,
     >                        r,
     >                        sigma,
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
      double precision  sigma
      double precision phi_ps(ngrid,nbasis)
      double precision matrix(nbasis,nbasis,mult_l+1)

#include "bafdecls.fh"
#include "integrate.fh"
#include "paw_utils.fh"
#include "paw_ma.fh"

      !*** local variables ***
      logical  ok
      integer  i,j,ii
      integer  la
      integer  power_r
      integer  i_f1,l_f1
      integer  i_f2,l_f2
      integer p1,p2
      double precision tmp


      power_r = 0

      !*** temporary array for calculatcomp_charges ***
      ok = ma_push_get(mt_dbl,ic,"tmp_array1",l_f1,i_f1)
      ok = ma_push_get(mt_dbl,ic,"tmp_array2",l_f2,i_f2)

      do la=0,mult_l
    
        call find_gaussian_rspace_array(la,sigma,ic,r,
     >                                  dbl_mb(i_f1))

        do ii=1,ic
        dbl_mb(i_f1+ii-1) =  dbl_mb(i_f1+ii-1)*r(ii)*r(ii)
        end do

          p1       = la +2
          p2       = 2*la + 2

        do i=1,nbasis
        do j=1,i

          do ii=1,ic
          dbl_mb(i_f2+ii-1) =  phi_ps(ii,i)*phi_ps(ii,j)
          end do

          tmp = find_multipole_energy(la,ic,r,
     >                                p1,dbl_mb(i_f1),
     >                                p2,dbl_mb(i_f2),log_amesh)

          matrix(i,j,la+1) = tmp
          matrix(j,i,la+1) = tmp

        end do
        end do

      end do

      ok = ma_chop_stack(l_f1)

      return
      end !subroutine find_comp_pot_matrix


      subroutine  find_gaussian_rspace_array(l,sigma_tmp,nr,r_tmp,gl)
      implicit none
      integer  l
      integer  nr
      double precision sigma_tmp
      double precision r_tmp(*)
      double precision gl(*)

#include "paw_params.fh"

      !*** local variables ***
      integer ii
      double precision  c
  
      double precision double_factorial
      external double_factorial


      c = 2**(l+2)/(sqrt(pi)*double_factorial(2*l+1)*sigma_tmp**(2*l+3))

      !*** this fixes possible underflow error ***
      do ii=1,nr
        gl(ii) = 0.0d0
      end do


      do ii=1,nr
        if (dabs(r_tmp(ii)) .lt. 8.0d0*sigma_tmp)
     >    gl(ii) = c*r_tmp(ii)**l*exp(-(r_tmp(ii)/sigma_tmp)**2)
      end do


      return
      end 
