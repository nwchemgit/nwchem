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
      subroutine  init_paw_comp_charge_matrix()
      implicit none

#include "bafdecls.fh"
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
      integer ic,istart

      tot_ntype = paw_basis_tot_ntype()
      ok = my_alloc(mt_int,tot_ntype,"i_comp_charge_matrix",
     >              i_comp_charge_matrix)
      ok = ok.and.
     >     my_alloc(mt_int,tot_ntype,"s_comp_charge_matrix",
     >              s_comp_charge_matrix)
      if (.not.ok)
     > call errquit(
     >  "init_paw_comp_charge_matrix: out of heap memory",0,0)

  
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
      if (.not.ok)
     > call errquit(
     >  "init_paw_comp_charge_matrix: out of heap memory",0,1)
  
      do it=1,tot_ntype
        ic    = paw_basis_npoints(it)
        istart = paw_basis_i_start(it)
    
        matr_index = int_mb(i_comp_charge_matrix(1) +it-1)
        call find_comp_charge_matrix(ic,istart,
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

      subroutine  end_paw_comp_charge_matrix()
      implicit none

#include "bafdecls.fh"
#include "paw_basis.fh"
#include "paw_comp_charge_matrix.fh"
#include "paw_utils.fh"
#include "paw_ma.fh"

      !*** local variables ***
      logical ok

      ok =        my_dealloc(i_comp_charge_matrix)
      ok = ok.and.my_dealloc(s_comp_charge_matrix)
      ok = ok.and.my_dealloc(comp_charge_matrix)
      if (.not.ok)
     > call errquit(
     >  "end_paw_comp_charge_matrix: error freeing heap",0,0)

      return
      end !subroutine end_paw_comp_charge_matrix

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
      subroutine find_comp_charge_matrix(ic,istart,
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
      integer    ic,istart
      integer    ngrid
      integer    nbasis
      integer    orb_l(nbasis)
      integer    mult_l

      double precision log_amesh
      double precision  r(ngrid)
      double precision phi_ae(ngrid,nbasis)
      double precision phi_ps(ngrid,nbasis)
      double precision matrix(nbasis,nbasis,mult_l+1)

#include "bafdecls.fh"
#include "integrate.fh"

      !*** local variables ***
      logical  ok
      integer  i,j,ig,igg,np
      integer  la
      integer  power_r
      integer  i_f,l_f
      double precision tmp


      call Parallel_np(np)
      power_r = 0

      !*** temporary array for calculatcomp_charges ***
      ok = ma_push_get(mt_dbl,ic,"tmp_array",l_f,i_f)
      if (.not.ok)
     > call errquit('find_comp_charge_matrix:out of stack memory',0,0)


      do i=1,nbasis
        do j=1,i

          do ig=1,ic
          igg=ig+istart-1
          dbl_mb(i_f+ig-1) = phi_ae(igg,i)*phi_ae(igg,j) 
     >                         - phi_ps(igg,i)*phi_ps(igg,j)
          end do

          do la=0,mult_l

c            tmp = def_integr(2*la+2,dbl_mb(i_f),la,r,log_amesh,ic)
            tmp = def_integr1(2*la+2,dbl_mb(i_f),la,r(istart),
     >                        log_amesh,ic,istart)

            
            matrix(i,j,la+1) = tmp
            matrix(j,i,la+1) = tmp

          end do

        end do
      end do
      if (np.gt.1) 
     > call D3dB_Vector_SumAll(nbasis*nbasis*(mult_l+1),matrix)

      ok = ma_pop_stack(l_f)
      if (.not.ok)
     > call errquit('find_comp_charge_matrix:error popping stack',0,0)

      return
      end !subroutine find_comp_charge_matrix

