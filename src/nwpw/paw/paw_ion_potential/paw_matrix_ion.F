*
* $Id: paw_matrix_ion.F,v 1.5 2003-03-06 01:46:40 bylaska Exp $
*

!**************************************************
!
!	name:
!
!	purpose 
!
!	created:	7/9/2002
!**************************************************
      subroutine  init_paw_ion_matrix()
      implicit none

#include "mafdecls.fh"
#include "paw_basis.fh"
#include "paw_matrix_ion_data.fh"
#include "paw_ma.fh"

      !*** local variables ***
      logical ok
      integer it,tot_ntype
      integer matr_size, matr_index

      tot_ntype = paw_basis_tot_ntype()
      ok = my_alloc(mt_int,tot_ntype,"i_paw_ion_matrix",
     >                                 i_paw_ion_matrix)
      if (.not.ok)
     > call errquit("init_paw_ion_matrix: out of heap memory",0)
  
      !*** determine matrix size ***
      matr_size = 0
      do it=1,tot_ntype
        matr_size = matr_size + int_mb(paw_basis_i_nbasis(it))**2
      end do

      !*** set the indexing to access diff atom types ***
      int_mb(i_paw_ion_matrix(1) )  = 0
      do it=1,tot_ntype-1
        int_mb(i_paw_ion_matrix(1) +it)  = 
     >                      int_mb(i_paw_ion_matrix(1) + it - 1) 
     >                    + int_mb(paw_basis_i_nbasis(it))**2
      end do

      ok = my_alloc(mt_dbl,matr_size,"paw_ion_matrix",paw_ion_matrix)
      if (.not.ok)
     > call errquit("init_paw_ion_matrix: out of heap memory",1)
  
      do it=1,tot_ntype
    
       matr_index = int_mb(i_paw_ion_matrix(1) +it-1)
       call find_ion_matrix( int_mb(paw_basis_i_i_cut(it)), 
     >                        int_mb(paw_basis_i_ngrid(it)),
     >                        int_mb(paw_basis_i_nbasis(it)),
     >                        int_mb(paw_basis_i_orb_l(it)),
     >                        dbl_mb(paw_basis_i_log_amesh(it)),
     >                        dbl_mb(paw_basis_i_rgrid(it)),
     >                        paw_basis_ion_charge(it),
     >                        dbl_mb(paw_basis_i_phi_ae(it)),
     >                        dbl_mb(paw_ion_matrix(1) + matr_index))
      end do

      return
      end !subroutine init_paw_ion_matrix

      subroutine  end_paw_ion_matrix()
      implicit none

#include "mafdecls.fh"
#include "paw_basis.fh"
#include "paw_matrix_ion_data.fh"
#include "paw_ma.fh"

      !*** local variables ***
      logical ok

      ok =        my_dealloc(i_paw_ion_matrix)
      ok = ok.and.my_dealloc(paw_ion_matrix)
      if (.not.ok)
     > call errquit("end_paw_ion_matrix: error freeing heap",0)

      return
      end !subroutine end_paw_ion_matrix
!*************************************************
!
!   name :   find_ion_matrix
!
!
!   purpose :
!
!
!   created :
!
!*************************************************
      subroutine find_ion_matrix( ic,
     >                        ngrid,
     >                        nbasis,
     >                        orb_l,
     >                        log_amesh,
     >                        r,
     >                        ion_charge,
     >                        phi_ae,
     >                        matrix)
      implicit none
      integer    ic
      integer    ngrid
      integer    nbasis
      integer    orb_l(nbasis)

      double precision     log_amesh
      double precision     r(ngrid)
      double precision     ion_charge
      double precision     phi_ae(ngrid,nbasis)
      double precision     matrix(nbasis,nbasis)

#include "mafdecls.fh"
#include "integrate.fh"

      !*** local variables ***
      logical  ok
      integer  i,j,ii
      integer  power_r,power_f
      integer  i_f,l_f
      double precision  tmp


      power_r = -1

      !*** temporary array for calculations ***
      ok = ma_push_get(mt_dbl,ic, "tmp_array",l_f,i_f)
      if (.not.ok)
     > call errquit('find_ion_matrix: out of stack memory',0)


      do i=1,nbasis
        power_f = 2*orb_l(i)+2
        do j=1,i
          if(orb_l(i).eq.orb_l(j)) then
            do ii=1,ic
              dbl_mb(i_f+ii-1) =  phi_ae(ii,i)*phi_ae(ii,j)
            end do

            tmp = def_integr(power_f,dbl_mb(i_f),power_r,r,
     >                       log_amesh,ic)
        
            tmp = -tmp*ion_charge
            matrix(i,j) = tmp
            matrix(j,i) = tmp
          else
            matrix(i,j) = 0.0d0
            matrix(j,i) = 0.0d0
          end if

        end do
      end do

      ok = ma_pop_stack(l_f)
      if (.not.ok)
     > call errquit('find_ion_matrix:error popping',1)

!  write(*,*) "ion matrix"
!  write(*,*) matrix

      return
      end !subroutine find_ion_matrix



      subroutine paw_ion_matrix_ptr(ptr)
      implicit none
      integer ptr

#include "paw_matrix_ion_data.fh"

      ptr = paw_ion_matrix(1)

      return
      end

