c
c $Id: paw_overlap_matrix.F,v 1.2 2003-02-11 01:41:23 edo Exp $
c

!**************************************************
!
!	name: init_paw_overlap_matrix
!
!	purpose 
!
!	created:	7/9/2002
!**************************************************
      subroutine  init_paw_overlap_matrix()
      implicit none
  
#include "paw_basis.fh"
#include "mafdecls.fh"
#include "paw_overlap_matrix.fh"
#include "paw_ma.fh"


      !*** local variables ***
      logical ok
      integer it,tot_ntype
      integer nbasis
      integer matr_size, matr_index

      tot_ntype = paw_basis_tot_ntype()
      ok = my_alloc(mt_int,tot_ntype,"i_overlap_matrix",
     >              i_overlap_matrix)
      ok = ok.and.
     >     my_alloc(mt_int,tot_ntype,"s_overlap_matrix",
     >              s_overlap_matrix)
      if (.not.ok)
     > call errquit("init_paw_overlap_matrix: error getting heap",0)
  
      !*** determine matrix size ***
      matr_size = 0
      do it=1,tot_ntype
        nbasis = paw_basis_nbasis(it)
        int_mb(s_overlap_matrix(1) + it - 1) = nbasis**2                   
        matr_size = matr_size + int_mb(s_overlap_matrix(1) + it - 1)                
      end do
  
      print*,"matrix_size=",matr_size
      !*** set the indexing to access diff atom types ***
      int_mb(i_overlap_matrix(1) )  = 0
      do it=1,tot_ntype-1
        int_mb(i_overlap_matrix(1) +it)  = 
     >            int_mb(i_overlap_matrix(1) + it - 1)  
     >          + int_mb(s_overlap_matrix(1) + it - 1)
      end do

      ok = my_alloc(mt_dbl,matr_size,"overlap_matrix",overlap_matrix)
      if (.not.ok)
     > call errquit("init_paw_overlap_matrix: error getting heap",1)
  
      do it=1,tot_ntype
    
        matr_index = int_mb(i_overlap_matrix(1) +it-1)
        call find_overlap_matrix( int_mb(paw_basis_i_i_cut(it)),
     >                          int_mb(paw_basis_i_ngrid(it)),
     >                          int_mb(paw_basis_i_nbasis(it)),
     >                          int_mb(paw_basis_i_orb_l(it)),
     >                          dbl_mb(paw_basis_i_log_amesh(it)),
     >                          dbl_mb(paw_basis_i_rgrid(it)),
     >                          dbl_mb(paw_basis_i_phi_ae(it)),
     >                          dbl_mb(paw_basis_i_phi_ps(it)),
     >                          dbl_mb(overlap_matrix(1) + matr_index))
  

      end do

      return
      end !subroutine init_paw_overlap_matrix

!**************************************************
!
!	name: end_paw_overlap_matrix
!
!	purpose 
!
!	created:	2/10/2003
!**************************************************
      subroutine  end_paw_overlap_matrix()
      implicit none
  
#include "paw_basis.fh"
#include "mafdecls.fh"
#include "paw_overlap_matrix.fh"
#include "paw_ma.fh"


      !*** local variables ***
      logical ok


      ok =        my_dealloc(i_overlap_matrix)
      ok = ok.and.my_dealloc(s_overlap_matrix)
      ok = ok.and.my_dealloc(overlap_matrix)
      if (.not.ok)
     > call errquit("end_paw_overlap_matrix: error freeing heap",0)



      return
      end !subroutine end_paw_overlap_matrix

!*************************************************
!
!   name :   find_overlap_matrix
!
!
!   purpose :
!
!
!   created : 2/10/2003
!
!*************************************************
      subroutine find_overlap_matrix( ic,
     >                        ngrid,
     >                        nbasis,
     >                        orb_l,
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

      double precision log_amesh
      double precision  r(ngrid)
      double precision phi_ae(ngrid,nbasis)
      double precision phi_ps(ngrid,nbasis)
      double precision matrix(nbasis,nbasis)

#include "integrate.fh"
#include "mafdecls.fh"

      !*** local variables ***
      logical  ok
      integer  i,j,ii
      integer  l
      integer  power_r,power_f
      integer  i_f,l_f
      double precision tmp


       power_r = 0

       !*** temporary array for calculatoverlaps ***
       ok = ma_push_get(mt_dbl,ic, "tmp_array", l_f,i_f)
       if (.not.ok)
     >  call errquit("find_overlap_matrix: out of stack memory",0)

       do i=1,nbasis
         l = orb_l(i)
         power_f = 2*orb_l(i)+2
         do j=1,i
           if(orb_l(i).eq.orb_l(j)) then

             do ii=1,ic
             dbl_mb(i_f+ii-1) =  phi_ae(ii,i)*phi_ae(ii,j) 
     >                            - phi_ps(ii,i)*phi_ps(ii,j)

             end do

             tmp = def_integr(power_f,dbl_mb(i_f),power_r,r,
     >                        log_amesh,ic)
        
             matrix(i,j) = tmp
             matrix(j,i) = tmp
           else
             matrix(i,j) = 0.0d0
             matrix(j,i) = 0.0d0
           end if

         end do
       end do

       ok = ma_chop_stack(l_f)
       if (.not.ok)
     >  call errquit("find_overlap_matrix: error popping stack",1)


      end !subroutine find_overlap_matrix


!*************************************************
!
!   name :   paw_overlap_matrix_gen
!
!
!   purpose : computes S = overlap1 + overlap2
!
!   created : 2/10/2003
!
!*************************************************
c      subroutine paw_overlap_matrix_gen(ne,psi1,psi2,S)
c      implicit none
c      integer    ne
c      complex*16 psi1(*),psi2(*)
c      real*8     S(ne,ne)
c
c#include "mafdecls.fh"
c#include "paw_ma.fh"
c#include "paw_basis.fh"
c#include "paw_overlap_matrix.fh"
c
c      !***** local variables ****
c      integer npack1
c
c      !*** calculate overlap1 ****
c      call Pack_npack(1,npack1)
c      call Grsm_ggm_dot(npack1,ne,
c     >                  psi1,psi2,
c     >                  S)
c
c      !*** calculate overlap2 ****
c      call paw_overlap_calc_scramble_coeff(ne,psi1,dcpl_mb(ct(1)))
c      call paw_overlap_unscramble_coeff(ne,dcpl_mb(ct(1)),c1)
c      call paw_overlap_calc_scramble_coeff(ne,psi2,ct)
c      call paw_overlap_unscramble_coeff(ne,ct,c2)
c
c      return
c      end

