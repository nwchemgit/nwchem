*
* $Id$
*

!**************************************************
!
!       name: init_paw_vloc_matrix
!
!       purpose: initializes local pseudopotential 
!                matrix
!
!       created:        7/9/2002
!**************************************************
      subroutine  init_paw_vloc_matrix()
      implicit none

#include "bafdecls.fh"
#include "paw_basis.fh"
#include "paw_matrix_vloc_data.fh"
#include "paw_ma.fh"

      !*** local variables ***
      logical ok
      integer it,tot_ntype
      integer matr_size, matr_index
      integer ic,istart

      tot_ntype = paw_basis_tot_ntype()
      ok = my_alloc(mt_int,tot_ntype,"i_paw_vloc_matrix",
     >                                  i_paw_vloc_matrix)
      if (.not.ok) 
     > call errquit("init_paw_vloc_matrix: error allocating heap",0,0)
  
      !*** determine matrix size ***
      matr_size = 0
      do it=1,tot_ntype
        matr_size = matr_size + int_mb(paw_basis_i_nbasis(it))**2
      end do

      !*** set the indexing to access diff atom types ***
      int_mb(i_paw_vloc_matrix(1) )  = 0
      do it=1,tot_ntype-1
        int_mb(i_paw_vloc_matrix(1) +it)  = 
     >         int_mb(i_paw_vloc_matrix(1) + it - 1)
     >       + int_mb(paw_basis_i_nbasis(it))**2
      end do

      ok = my_alloc(mt_dbl,matr_size,"paw_vloc_matrix",paw_vloc_matrix)
      if (.not.ok) 
     > call errquit("init_paw_vloc_matrix: error allocating heap",0,1)
  
      do it=1,tot_ntype
       ic     = paw_basis_npoints(it)
       istart = paw_basis_i_start(it)
    
       matr_index = int_mb(i_paw_vloc_matrix(1) +it-1)
       call find_paw_vloc_matrix(ic,istart,
     >                        int_mb(paw_basis_i_ngrid(it)),
     >                        int_mb(paw_basis_i_nbasis(it)),
     >                        int_mb(paw_basis_i_orb_l(it)),
     >                        dbl_mb(paw_basis_i_log_amesh(it)),
     >                        dbl_mb(paw_basis_i_rgrid(it)),
     >                        dbl_mb(paw_basis_i_v_ps(it)),
     >                        dbl_mb(paw_basis_i_phi_ps(it)),
     >                        dbl_mb(paw_vloc_matrix(1) + matr_index))


      end do

      return
      end !subroutine init_paw_vloc_matrix


      subroutine  end_paw_vloc_matrix()
      implicit none

#include "bafdecls.fh"
#include "paw_basis.fh"
#include "paw_matrix_vloc_data.fh"
#include "paw_ma.fh"

      !*** local variables ***
      logical ok

      ok =        my_dealloc(i_paw_vloc_matrix)
      ok = ok.and.my_dealloc(paw_vloc_matrix)
      if (.not.ok) 
     > call errquit("end_paw_vloc_matrix: error freeing heap",0,0)


      return
      end !subroutine end_paw_vloc_matrix

!**************************************************
!
!       name:  find_paw_vloc_matrix
!
!       purpose: 
!
!       created:        7/9/2002
!**************************************************
      subroutine find_paw_vloc_matrix(ic,istart,
     >                        ngrid,
     >                        nbasis,
     >                        orb_l,
     >                        log_amesh,
     >                        r,
     >                        v_ps,
     >                        phi_ps,
     >                        matrix)
      implicit none
      integer    ic,istart
      integer    ngrid
      integer    nbasis
      integer    orb_l(nbasis)

      double precision     log_amesh
      double precision     r(ngrid)
      double precision     v_ps(ngrid)
      double precision     phi_ps(ngrid,nbasis)
      double precision     matrix(nbasis,nbasis)

#include "bafdecls.fh"
#include "integrate.fh"

      !*** local variables ***
      logical  ok
      integer  i,j,ig,igg,np
      integer  l
      integer  power_r,power_f
      integer  i_f,l_f
      double precision  tmp


      call Parallel_np(np)
      power_r = 0

      !*** temporary array for calculations ***
      ok = ma_push_get(mt_dbl,ic, "tmp_array",l_f,i_f)
      if (.not.ok)
     > call errquit('find_paw_vloc_matrix: out of stack memory',0,0)


      do i=1,nbasis
        l = orb_l(i)
        power_f = 2*orb_l(i) + 2
        do j=1,i
          if(orb_l(i).eq.orb_l(j)) then
            
            do ig=1,ic
            igg=ig+istart-1
            dbl_mb(i_f+ig-1) =  phi_ps(igg,i)
     >                             *phi_ps(igg,j)
     >                             *v_ps(igg)
            end do
c            tmp = def_integr(power_f,dbl_mb(i_f),power_r,
c     >                       r,log_amesh,ic)
            tmp = def_integr1(power_f,dbl_mb(i_f),power_r,
     >                       r(istart),log_amesh,ic,istart)
        
            matrix(i,j) = tmp
            matrix(j,i) = tmp
          else
            matrix(i,j) = 0.0d0
            matrix(j,i) = 0.0d0
          end if
          
        end do
      end do
      if (np.gt.1) call D3dB_Vector_SumAll(nbasis*nbasis,matrix)

      ok = ma_pop_stack(l_f)
      if (.not.ok)
     > call errquit('find_paw_vloc_matrix:error popping stack',0,0)

      return
      end !subroutine find_paw_vloc_matrix



      subroutine paw_vloc_matrix_ptr(ptr)
      implicit none
      integer ptr

#include "paw_matrix_vloc_data.fh"

      ptr = paw_vloc_matrix(1)

      return
      end
