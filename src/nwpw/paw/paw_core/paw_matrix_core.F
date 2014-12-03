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
      subroutine  init_paw_core_matrix()
      implicit none
  
#include "bafdecls.fh"
#include "paw_basis.fh"
#include "paw_matrix_core_data.fh"
#include "paw_utils.fh"
#include "paw_ma.fh"


      !*** local variables ***
      logical ok
      integer it,tot_ntype
      integer matr_size, matr_index

      tot_ntype = paw_basis_tot_ntype()
      ok = my_alloc(mt_int,tot_ntype,"i_paw_core_matrix",
     >                   i_paw_core_matrix)
      if (.not.ok) 
     > call errquit("init_paw_core_matrix: out of memory",0,0)
  
      !*** determine matrix size ***
      matr_size = 0
      do it=1,tot_ntype
        matr_size = matr_size + int_mb(paw_basis_i_nbasis(it))**2
      end do

      !*** set the indexing to access diff atom types ***
      int_mb(i_paw_core_matrix(1) )  = 0
      do it=1,tot_ntype-1
        int_mb(i_paw_core_matrix(1) +it) = 
     >         int_mb(i_paw_core_matrix(1) + it - 1) 
     >       + int_mb(paw_basis_i_nbasis(it))**2
      end do

      ok = my_alloc(mt_dbl,matr_size,"paw_core_matrix",paw_core_matrix)
      if (.not.ok) 
     > call errquit("init_paw_core_matrix: out of memory",0,1)
  
      do it=1,tot_ntype
    
        matr_index = int_mb(i_paw_core_matrix(1) +it-1)
        call find_core_matrix( int_mb(paw_basis_i_i_cut(it)),
     >                        int_mb(paw_basis_i_ngrid(it)),
     >                        int_mb(paw_basis_i_nbasis(it)),
     >                        int_mb(paw_basis_i_orb_l(it)),
     >                        dbl_mb(paw_basis_i_log_amesh(it)),
     >                        dbl_mb(paw_basis_i_rgrid(it)),
     >                        dbl_mb(paw_basis_i_core_density(it)),
     >                        dbl_mb(paw_basis_i_ps_core_density(it)),
     >                        dbl_mb(paw_basis_i_phi_ae(it)),
     >                        dbl_mb(paw_basis_i_phi_ps(it)),
     >                        dbl_mb(paw_core_matrix(1) + matr_index))
      end do

      return
      end !subroutine init_paw_core_matrix

      subroutine  end_paw_core_matrix()
      implicit none
  
#include "bafdecls.fh"
#include "paw_basis.fh"
#include "paw_matrix_core_data.fh"
#include "paw_utils.fh"
#include "paw_ma.fh"

      !*** local variables ***
      logical ok

      ok =        my_dealloc(i_paw_core_matrix)
      ok = ok.and.my_dealloc(paw_core_matrix)
      if (.not.ok) 
     > call errquit("end_paw_core_matrix: error freeing heap",0,0)

      return
      end !subroutine end_paw_core_matrix


!**************************************************
!
!       name:
!
!       purpose 
!
!       created:        7/9/2002
!**************************************************
      subroutine find_core_matrix( ic,
     >                        ngrid,
     >                        nbasis,
     >                        orb_l,
     >                        log_amesh,
     >                        r,
     >                        core_density,
     >                        ps_core_density,
     >                        phi_ae,
     >                        phi_ps,
     >                        matrix)
      implicit none
      integer    ic
      integer    ngrid
      integer    nbasis
      integer    orb_l(nbasis)

      double precision     log_amesh
      double precision     r(ngrid)
      double precision     core_density(ngrid)
      double precision     ps_core_density(ngrid)
      double precision     phi_ae(ngrid,nbasis)
      double precision     phi_ps(ngrid,nbasis)
      double precision     matrix(nbasis,nbasis)

#include "bafdecls.fh"
#include "integrate.fh"
#include "paw_utils.fh"   
#include "paw_ma.fh"   

      !*** local variables ***
      logical  ok
      integer  i,j,ii
      integer  l
      integer  power_r,power_f
      integer  i_f1,l_f1
      integer  i_f2,l_f2
      integer  i_f3,l_f3
      integer  i_f4,l_f4
      double precision  tmp


      power_r = 0

      !*** temporary array for calculations ***
      ok =        ma_push_get(mt_dbl,ic,"tmp_array1",l_f1,i_f1)
      ok = ok.and.ma_push_get(mt_dbl,ic,"tmp_array2",l_f2,i_f2)
      ok = ok.and.ma_push_get(mt_dbl,ic,"tmp_array3",l_f3,i_f3)
      ok = ok.and.ma_push_get(mt_dbl,ic,"tmp_array4",l_f4,i_f4)
      if (.not.ok)
     > call errquit('find_core_matrix: out of stack memory',0,0)


      do ii=1,ic
      dbl_mb(i_f3+ii-1) = core_density(ii)*r(ii)**2
      dbl_mb(i_f4+ii-1) = ps_core_density(ii)*r(ii)**2
      end do

      do i=1,nbasis
        l = orb_l(i)
        power_f = 2*orb_l(i) + 2
        do j=1,i
          if(orb_l(i).eq.orb_l(j)) then
    
            power_f = 2*l+2
            do ii=1,ic
            dbl_mb(i_f1+ii-1) =  phi_ae(ii,i)*phi_ae(ii,j)
            dbl_mb(i_f2+ii-1) =  phi_ps(ii,i)*phi_ps(ii,j)
            end do
   
            tmp = find_multipole_energy(0,ic,r,power_f,
     >                                  dbl_mb(i_f1),2,
     >                                  dbl_mb(i_f3),log_amesh)
     >          - find_multipole_energy(0,ic,r,power_f,
     >                                  dbl_mb(i_f2),2,
     >                                  dbl_mb(i_f4),log_amesh)

            matrix(i,j) = tmp
            matrix(j,i) = tmp
          else
            matrix(i,j) = 0.0d0
            matrix(j,i) = 0.0d0
          end if

        end do
      end do

      ok =        ma_pop_stack(l_f4)
      ok = ok.and.ma_pop_stack(l_f3)
      ok = ok.and.ma_pop_stack(l_f2)
      ok = ok.and.ma_pop_stack(l_f1)
      if (.not.ok)
     > call errquit('find_core_matrix:error popping stack',0,1)

      return
      end



      subroutine paw_core_matrix_ptr(ptr)
      implicit none
      integer ptr

#include "paw_matrix_core_data.fh"

      ptr = paw_core_matrix(1)

      return
      end


