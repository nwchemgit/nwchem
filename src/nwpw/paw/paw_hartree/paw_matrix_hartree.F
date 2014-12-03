*
* $Id$
*

!*************************************************
!
!   name :   find_paw_hartree_matrix
!
!
!   purpose :
!
!
!   created :
!
!*************************************************
      subroutine find_paw_hartree_matrix(ic,
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
      double precision r(ngrid)
      double precision phi_ae(ngrid,nbasis)
      double precision phi_ps(ngrid,nbasis)
      double precision matrix(nbasis,nbasis,nbasis,nbasis,mult_l+1)

#include "bafdecls.fh"
#include "integrate.fh"
#include "paw_utils.fh"   
#include "paw_ma.fh"   

      !*** local variables ***
      logical  ok
      integer  i,j,ip,jp,ii
      integer  la
      integer  power_r
      integer  i_f1,l_f1
      integer  i_f2,l_f2
      integer  i_f3,l_f3
      integer  i_f4,l_f4
      integer  sum_l,diff_l
      integer  sum_lp,diff_lp
      integer  p1,p2
      double precision  tmp


      power_r = 0
      la = nbasis*nbasis*nbasis*nbasis*(mult_l+1)
      call dcopy(la,0.0d0,0,matrix,1)

      !*** temporary array for calculations ***
      ok =        ma_push_get(mt_dbl,ic,"tmp_array1",l_f1,i_f1)
      ok = ok.and.ma_push_get(mt_dbl,ic,"tmp_array2",l_f2,i_f2)
      ok = ok.and.ma_push_get(mt_dbl,ic,"tmp_array3",l_f3,i_f3)
      ok = ok.and.ma_push_get(mt_dbl,ic,"tmp_array4",l_f4,i_f4)
      if (.not.ok)
     > call errquit('find_paw_hartree_matrix:out of stack memory',0,1)
  
      do i=1,nbasis
      do j=1,i

        sum_l  = orb_l(i) + orb_l(j)
        diff_l = ABS(orb_l(i) - orb_l(j))
        p1     = sum_l + 2

        do ii=1,ic
        dbl_mb(i_f1+ii-1) = phi_ae(ii,i)*phi_ae(ii,j)
        dbl_mb(i_f3+ii-1) = phi_ps(ii,i)*phi_ps(ii,j)
        end do

        do ip=1,nbasis
        do jp=1,ip       

          sum_lp  = orb_l(ip) + orb_l(jp)
          diff_lp = ABS(orb_l(ip) - orb_l(jp))
          p2    = sum_lp + 2
        
          do ii=1,ic
             dbl_mb(i_f2+ii-1) = phi_ae(ii,ip)*phi_ae(ii,jp)
             dbl_mb(i_f4+ii-1) = phi_ps(ii,ip)*phi_ps(ii,jp)
          end do

          do la=0,mult_l
            if(la .le. sum_l .and. la .ge. diff_l ) then
            if(la .le. sum_lp .and. la .ge. diff_lp ) then
  
            tmp = find_multipole_energy(la,ic,r,p1,
     >                                  dbl_mb(i_f1),p2,
     >                                  dbl_mb(i_f2),log_amesh)
     >          - find_multipole_energy(la,ic,r,p1,
     >                                  dbl_mb(i_f3),p2,
     >                                  dbl_mb(i_f4),log_amesh)

            matrix(i,j,ip,jp,la+1) = tmp
            matrix(j,i,ip,jp,la+1) = tmp
            matrix(i,j,jp,ip,la+1) = tmp
            matrix(j,i,jp,ip,la+1) = tmp

            end if
            end if
          end do
        end do
        end do
      end do
      end do


      ok =        ma_pop_stack(l_f4)
      ok = ok.and.ma_pop_stack(l_f3)
      ok = ok.and.ma_pop_stack(l_f2)
      ok = ok.and.ma_pop_stack(l_f1)
      if (.not.ok)
     > call errquit('find_paw_hartree_matrix:error popping stack',0,2)

!  write(*,*) "hartree matrix"
!  write(*,*) matrix


      return
      end

!**************************************************
!
!       name:
!
!       purpose 
!
!       created:        7/9/2002
!**************************************************
      subroutine  init_paw_hartree_matrix()
      implicit none

#include "bafdecls.fh"
#include "paw_basis.fh"
#include "paw_matrix_hartree_data.fh"
#include "paw_utils.fh"
#include "paw_ma.fh"


      !*** local variables ***
      logical ok
      integer it,tot_ntype
      integer mult_l
      integer nbasis
      integer matr_size, mtr_index

      tot_ntype = paw_basis_tot_ntype()
      ok = my_alloc(mt_int,tot_ntype,"i_paw_hartree_matrix",
     >              i_paw_hartree_matrix)
      ok = ok.and.
     >     my_alloc(mt_int,tot_ntype,"s_paw_hartree_matrix",
     >              s_paw_hartree_matrix)
      if (.not.ok) call errquit(
     > 'init_paw_paw_hartree_matrix: out of heap memory',0,0)
  
      !*** determine matrix size ***
      matr_size = 0
      do it=1,tot_ntype
        mult_l  = paw_basis_mult_l(it)
        nbasis = paw_basis_nbasis(it)

        int_mb(s_paw_hartree_matrix(1) + it - 1) = nbasis**4*(mult_l+1) 
        matr_size = matr_size + int_mb(s_paw_hartree_matrix(1) + it - 1)
      end do
  

      !*** set the indexing to access diff atom types ***
      int_mb(i_paw_hartree_matrix(1) )  = 0
      do it=1,tot_ntype-1
        int_mb(i_paw_hartree_matrix(1) +it) =
     >         int_mb(i_paw_hartree_matrix(1) + it - 1)
     >       + int_mb(s_paw_hartree_matrix(1) + it - 1)
      end do

      ok = my_alloc(mt_dbl,matr_size,"paw_hartree_matrix",
     >     paw_hartree_matrix)
      if (.not.ok)
     > call errquit("init_paw_hartree_matrix: out of heap memory",0,1)
  
      do it=1,tot_ntype
    
        mtr_index = int_mb(i_paw_hartree_matrix(1) +it-1)
        call find_paw_hartree_matrix( int_mb(paw_basis_i_i_cut(it)),
     >                        int_mb(paw_basis_i_ngrid(it)),
     >                        int_mb(paw_basis_i_nbasis(it)),
     >                        int_mb(paw_basis_i_orb_l(it)),
     >                        paw_basis_mult_l(it),
     >                        dbl_mb(paw_basis_i_log_amesh(it)),
     >                        dbl_mb(paw_basis_i_rgrid(it)),
     >                        dbl_mb(paw_basis_i_phi_ae(it)),
     >                        dbl_mb(paw_basis_i_phi_ps(it)),
     >                        dbl_mb(paw_hartree_matrix(1) + mtr_index))
      end do

      return
      end !subroutine init_paw_hartree_matrix


      subroutine  end_paw_hartree_matrix()
      implicit none

#include "bafdecls.fh"
#include "paw_basis.fh"
#include "paw_matrix_hartree_data.fh"
#include "paw_utils.fh"
#include "paw_ma.fh"

      !**** local variables ****
      logical ok

      ok = my_dealloc(i_paw_hartree_matrix)
      ok = ok.and.my_dealloc(s_paw_hartree_matrix)
      ok = ok.and.my_dealloc(paw_hartree_matrix)
      if (.not.ok)
     > call errquit("end_paw_hartree_matrix: error freeing heap",0,0)

      return
      end !subroutine end_paw_hartree_matrix



