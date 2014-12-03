c
c $Id$
c

!**************************************************
!
!       name: init_paw_overlap_matrix
!
!       purpose 
!
!       created:        7/9/2002
!**************************************************
      subroutine  init_paw_overlap_matrix()
      implicit none
  
#include "paw_basis.fh"
#include "bafdecls.fh"
#include "paw_overlap_matrix.fh"
#include "paw_ma.fh"


      !*** local variables ***
      logical ok
      integer it,tot_ntype
      integer nbasis
      integer matr_size, matr_index
      integer ic,istart

      tot_ntype = paw_basis_tot_ntype()
      ok = my_alloc(mt_int,tot_ntype,"i_overlap_matrix",
     >              i_overlap_matrix)
      ok = ok.and.
     >     my_alloc(mt_int,tot_ntype,"s_overlap_matrix",
     >              s_overlap_matrix)
      if (.not.ok)
     > call errquit("init_paw_overlap_matrix: error getting heap",0,0)

  
      !*** determine matrix size ***
      matr_size = 0
      do it=1,tot_ntype
        nbasis = paw_basis_nbasis(it)
        int_mb(s_overlap_matrix(1) + it - 1) = nbasis**2
        matr_size = matr_size + int_mb(s_overlap_matrix(1) + it - 1)
      
      end do
  
      
      !*** set the indexing to access diff atom types ***
      int_mb(i_overlap_matrix(1) )  = 0
      do it=1,tot_ntype-1
        int_mb(i_overlap_matrix(1) +it)  = 
     >            int_mb(i_overlap_matrix(1) + it - 1)  
     >          + int_mb(s_overlap_matrix(1) + it - 1)
      end do

      ok = my_alloc(mt_dbl,matr_size,"overlap_matrix",overlap_matrix)
      if (.not.ok)
     > call errquit("init_paw_overlap_matrix: error getting heap",0,1)
  
      do it=1,tot_ntype
        ic     = paw_basis_npoints(it)
        istart = paw_basis_i_start(it)

        matr_index = int_mb(i_overlap_matrix(1) +it-1)
        call find_overlap_matrix(ic,istart,
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
!       name: end_paw_overlap_matrix
!
!       purpose 
!
!       created:        2/10/2003
!**************************************************
      subroutine  end_paw_overlap_matrix()
      implicit none
  
#include "paw_basis.fh"
#include "bafdecls.fh"
#include "paw_overlap_matrix.fh"
#include "paw_ma.fh"


      !*** local variables ***
      logical ok


      ok =        my_dealloc(i_overlap_matrix)
      ok = ok.and.my_dealloc(s_overlap_matrix)
      ok = ok.and.my_dealloc(overlap_matrix)
      if (.not.ok)
     > call errquit("end_paw_overlap_matrix: error freeing heap",0,0)



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
      subroutine find_overlap_matrix(ic,istart,
     >                        ngrid,
     >                        nbasis,
     >                        orb_l,
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

      double precision log_amesh
      double precision  r(ngrid)
      double precision phi_ae(ngrid,nbasis)
      double precision phi_ps(ngrid,nbasis)
      double precision matrix(nbasis,nbasis)

#include "integrate.fh"
#include "bafdecls.fh"

      !*** local variables ***
      logical  ok
      integer  i,j,ig,igg,np
      integer  l
      integer  power_r,power_f
      integer  i_f,l_f
      double precision tmp


       call Parallel_np(np)
       power_r = 0

       !*** temporary array for calculatoverlaps ***
       ok = ma_push_get(mt_dbl,ic, "tmp_array", l_f,i_f)
       if (.not.ok)
     >  call errquit("find_overlap_matrix: out of stack memory",0,0)

       do i=1,nbasis
         l = orb_l(i)
         power_f = 2*orb_l(i)+2
         do j=1,i
           if(orb_l(i).eq.orb_l(j)) then

             do ig=1,ic
             igg=ig+istart-1
             dbl_mb(i_f+ig-1) =  phi_ae(igg,i)*phi_ae(igg,j) 
     >                            - phi_ps(igg,i)*phi_ps(igg,j)

             end do

c             tmp = def_integr(power_f,dbl_mb(i_f),power_r,r,
c     >                        log_amesh,ic)
             tmp = def_integr1(power_f,dbl_mb(i_f),power_r,r(istart),
     >                        log_amesh,ic,istart)
        
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
     >  call errquit("find_overlap_matrix: error popping stack",0,1)


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
      subroutine paw_overlap_matrix_gen(n,ne,psi1,psi2,Sij)
      implicit none
      integer    n,ne
      complex*16 psi1(*),psi2(*)
      real*8     Sij(n,n)

#include "bafdecls.fh"
#include "paw_ma.fh"
#include "paw_basis.fh"
#include "paw_proj.fh"
#include "paw_geom.fh"
#include "paw_ovlp_data.fh"
#include "paw_overlap_matrix.fh"

      !*** local variables ***
      logical ok
      integer npack1,indxr,indxs,indx
      integer i,j,ii,ia,a,b,ma,mb,la,lb,r,s
      integer nbasis,basis_nbasis,nion,ii_shift
      integer matrx_ptr
      integer c1(2),c2(2),ct(2)


      !*** calculate overlap1 ****

      call Pack_npack(1,npack1)

*    **** neede to change Grsm_ggm_dot for unrestricted calculations ***
      call Grsm_ggm2_dot(npack1,ne,    
     >                  psi1,psi2,
     >                  Sij,n)

      !*** allocate space from stack ****
      ok = BA_push_get(mt_dcpl,ne*total_nbasis,'c1',c1(2),c1(1))
      ok = ok.and.
     >     BA_push_get(mt_dcpl,ne*total_nbasis,'c2',c2(2),c2(1))
      ok = ok.and.
     >     BA_push_get(mt_dcpl,ne*total_nbasis,'ct',ct(2),ct(1))
      if (.not.ok)
     > call errquit('paw_overlap_matrix_gen: out of stack memory',0,0)

      !*** calculate overlap2 ****

      !**** ct_r^a(j) = Sum(rs) (q_rs^a)*c_s^a(j) ****
      call paw_ovlp_calc_coef(ne,psi2,dcpl_mb(c2(1)))
      call  paw_overlap_pv_gen(ne,total_nbasis,
     >          dcpl_mb(c2(1)),
     >          dcpl_mb(ct(1)))


      !**** S(i,j) = Sum(i,j,r) dconjg( c1_r^a(i) ) * ct_r^a(j) ****
      call paw_ovlp_calc_coef(ne,psi1,dcpl_mb(c1(1)))
      ii_shift = 0
      nion = ion_nion()
      do ii=1,nion
        ia     = ion_katm(ii)
        nbasis = paw_proj_nbasis(ia)
        do r=1,nbasis
           do i=1,ne
           do j=1,ne
             indxr = (r-1+ii_shift)*ne + (i-1)
             indxs = (r-1+ii_shift)*ne + (j-1)
             
             Sij(i,j) = Sij(i,j) 
     >                + dconjg(dcpl_mb(c1(1)+indxr))
     >                 *dcpl_mb(ct(1)+indxs)
           end do 
           end do 
        end do
        ii_shift = ii_shift+nbasis
      end do


      ok =        BA_pop_stack(ct(2))
      ok = ok.and.BA_pop_stack(c2(2))
      ok = ok.and.BA_pop_stack(c1(2))
      if (.not.ok)
     > call errquit('paw_overlap_matrix_gen: error popping stack',0,0)
      return
      end


!*************************************************
!
!   name :   paw_overlap_sym_matrix_gen
!
!
!   purpose : computes S = overlap1 + overlap2
!
!   created : 2/10/2003
!
!*************************************************
      subroutine paw_overlap_sym_matrix_gen(n,ne,psi1,psi2,Sij)
      implicit none
      integer    n,ne
      complex*16 psi1(*),psi2(*)
      real*8     Sij(n,n)

#include "bafdecls.fh"
#include "paw_ma.fh"
#include "paw_basis.fh"
#include "paw_proj.fh"
#include "paw_geom.fh"
#include "paw_ovlp_data.fh"
#include "paw_overlap_matrix.fh"

      !*** local variables ***
      logical ok
      integer npack1,indxr,indxs,indx
      integer i,j,ii,ia,a,b,ma,mb,la,lb,r,s
      integer nbasis,basis_nbasis,nion,ii_shift
      integer matrx_ptr
      integer c1(2),c2(2),ct(2)


      !*** calculate overlap1 ****

      call Pack_npack(1,npack1)

*    **** neede to change Grsm_ggm_dot for unrestricted calculations ***
      call Grsm_ggm2_sym_dot(npack1,ne,
     >                  psi1,psi2,
     >                  Sij,n)

      !*** allocate space from stack ****
      ok = BA_push_get(mt_dcpl,ne*total_nbasis,'c1',c1(2),c1(1))
      ok = ok.and.
     >     BA_push_get(mt_dcpl,ne*total_nbasis,'c2',c2(2),c2(1))
      ok = ok.and.
     >     BA_push_get(mt_dcpl,ne*total_nbasis,'ct',ct(2),ct(1))
      if (.not.ok)
     > call errquit('paw_overlap_matrix_gen: out of stack memory',0,0)

      !*** calculate overlap2 ****

      !**** ct_r^a(j) = Sum(rs) (q_rs^a)*c_s^a(j) ****
      call paw_ovlp_calc_coef(ne,psi2,dcpl_mb(c2(1)))
      call  paw_overlap_pv_gen(ne,total_nbasis,
     >          dcpl_mb(c2(1)),
     >          dcpl_mb(ct(1)))


      !**** S(i,j) = Sum(i,j,r) dconjg( c1_r^a(i) ) * ct_r^a(j) ****
      call paw_ovlp_calc_coef(ne,psi1,dcpl_mb(c1(1)))
      ii_shift = 0
      nion = ion_nion()
      do ii=1,nion
        ia     = ion_katm(ii)
        nbasis = paw_proj_nbasis(ia)
        do r=1,nbasis
           do i=1,ne
           do j=1,ne
             indxr = (r-1+ii_shift)*ne + (i-1)
             indxs = (r-1+ii_shift)*ne + (j-1)

             Sij(i,j) = Sij(i,j)
     >                + dconjg(dcpl_mb(c1(1)+indxr))
     >                 *dcpl_mb(ct(1)+indxs)
           end do
           end do
        end do
        ii_shift = ii_shift+nbasis
      end do


      ok =        BA_pop_stack(ct(2))
      ok = ok.and.BA_pop_stack(c2(2))
      ok = ok.and.BA_pop_stack(c1(2))
      if (.not.ok)
     > call errquit('paw_overlap_matrix_gen: error popping stack',0,0)
      return
      end





!*************************************************
!
!   name :   paw_overlap_pv_gen
!
!
!   purpose : computes pv
!
!   created : 2/10/2003
!
!*************************************************

      subroutine paw_overlap_pv_gen(ne,total_nbasis,c1,pv)
      implicit none
      integer    ne,total_nbasis
      complex*16 c1(ne,total_nbasis)
      complex*16 pv(ne,total_nbasis) !* RESULT *

#include "bafdecls.fh"
#include "paw_basis.fh"
#include "paw_proj.fh"
#include "paw_geom.fh"
#include "paw_overlap_matrix.fh"

      !*** local variables ***
      integer indx
      integer j,ii,ia,a,b,ma,mb,la,lb,r,s
      integer nbasis,basis_nbasis,nion,ii_shift
      integer matrx_ptr


*     **** pv^a(j) = Sum(rs) (q_rs^a)*c_s^a(j) ****
      call dcopy(2*ne*total_nbasis,0.0d0,0,pv,1)
      nion = ion_nion()
      ii_shift=0
      do ii=1,nion
        ia           = ion_katm(ii)
        nbasis       = paw_proj_nbasis(ia)
        basis_nbasis = paw_basis_nbasis(ia)
        matrx_ptr    = int_mb(i_overlap_matrix(1)+ia-1)

        r = 0
        do a=1,basis_nbasis
          la=paw_basis_orb_l(a,ia)
          do ma=-la,la
            r = r + 1

            s = 0
            do b=1,basis_nbasis
              lb    = paw_basis_orb_l(b,ia)
              indx  = matrx_ptr + (b-1)*basis_nbasis + (a-1)
              do mb=-lb,lb
                 s = s + 1
                 if ((mb.eq.ma).and.(la.eq.lb)) then
                   do j=1,ne
                     pv(j,r+ii_shift)
     >                 = pv(j,r+ii_shift)
     >                 + dbl_mb(overlap_matrix(1)+indx)
     >                  *c1(j,s+ii_shift)
                   end do !*j*
                 end if
              end do !*mb*
            end do !*b*
          end do !*ma*
        end do   !*a*

        ii_shift = ii_shift+nbasis
      end do !*ii*

      return
      end


!**************************************************
!       Name: paw_overlap_matrix_ptr
!
!       Purpose:
!
!       Created:        2/10/2002
!**************************************************
      integer function paw_overlap_matrix_ptr(ia)
      implicit none
      integer ia

#include "bafdecls.fh"
#include "paw_overlap_matrix.fh"

      paw_overlap_matrix_ptr 
     >   = overlap_matrix(1) 
     >  + int_mb(i_overlap_matrix(1)+ia-1)

      return
      end

