c
c $Id$
c

!**************************************************
!
!       name:
!
!       purpose 
!
!       created:        7/9/2002
!**************************************************
      subroutine  init_paw_pot_comp()
      implicit none

#include "bafdecls.fh"
#include "paw_basis.fh"
#include "paw_proj.fh"
#include "paw_pot_comp_data.fh"
#include "paw_ma.fh"
#include "paw_geom.fh"

      !*** local variables ***
      logical ok
      integer in,ia
      integer paw_pot_comp_size

  
c     !*** calculate the size and indexing ***
      ok = my_alloc(mt_int,ion_nion(),"i_paw_pot_comp",
     >                                     i_paw_pot_comp)
      if (.not.ok) 
     > call errquit('init_paw_pot_comp:out of heap memory',0,0)
      paw_pot_comp_size = 0
      do in=1,ion_nion()
         int_mb(i_paw_pot_comp(1) + in - 1) = paw_pot_comp_size
         ia = ion_katm(in)
         paw_pot_comp_size = paw_pot_comp_size
     >                + paw_proj_nbasis(ia)**2
      end do


      ok = my_alloc(mt_dcpl,paw_pot_comp_size,
     >               "paw_pot_comp",paw_pot_comp)
      if (.not.ok)
     > call errquit("init_paw_pot_comp:out of heap memoryh",0,1)
  
      !write(*,*) "initialized paw_comp_pot"
      end !subroutine init_paw_pot_comp

!**************************************************
!
!       name:
!
!       purpose 
!
!       created:        7/9/2002
!**************************************************
      subroutine  end_paw_pot_comp()
      implicit none

#include "paw_basis.fh"
#include "bafdecls.fh"
#include "paw_pot_comp_data.fh"
#include "paw_ma.fh"

      !*** local variables ***
      logical ok

      ok =        my_dealloc(i_paw_pot_comp)
      ok = ok.and.my_dealloc(paw_pot_comp)
      if (.not.ok)
     > call errquit("end_paw_pot_comp: error freeing heap",0,0)

      return
      end !subroutine end_paw_pot_comp

!**************************************************
!
!       name:
!
!       purpose 
!
!       created:        7/9/2002
!************************************************
      subroutine paw_pot_comp_solve()
      
      implicit none

#include "bafdecls.fh"
#include "paw_geom.fh"
#include "paw_ma.fh"
#include "paw_matrix_comp_pot.fh"
#include "paw_basis.fh"
#include "paw_proj.fh"
#include "paw_pot_comp_data.fh"
#include "paw_comp_charge_data.fh"

      double precision gen_gaunt_coeff
      external gen_gaunt_coeff
      integer in,ia
      integer i_mtr,i_mtr0
      integer i,j,li,lj,mi,mj
      integer nb,nb2,nb3,nb4
      integer l,m,ilm
      integer mult_l
      integer i_cp,i_cp0
      integer nilm,njlm
      integer i_qlm,i_qlm0
      complex*16 tmp_comp_pot

*       *** initializing to zero ***
        do i=1,paw_pot_comp(3)
           dcpl_mb(paw_pot_comp(1)+i-1) = dcmplx(0.0d0,0.0d0)
        end do

      do in=1,ion_nion()
        ia = ion_katm(in)
        nb = paw_basis_nbasis(ia)
        nb2 = nb*nb
        nb3 = nb2*nb
        nb4 = nb3*nb
        mult_l = paw_basis_mult_l(ia)
        i_mtr0 = int_mb(i_paw_comp_pot_matrix(1) + ia - 1)
        i_qlm0 = int_mb(i_paw_qlm(1) + in - 1)
        i_cp0 = int_mb(i_paw_pot_comp(1) + in - 1)
        nilm = 0
        do i=1,nb
        li = paw_basis_orb_l(i,ia)
        njlm = 0
        do j=1,nb
        lj = paw_basis_orb_l(j,ia)
        ilm = 0
        do l=0,mult_l
        do m=-l,l
          do mi=-li,li
          do mj=-lj,lj
          if ( m.eq.(mi-mj)) then
             i_mtr = i_mtr0-1+i+(j-1)*nb+l*nb2
             i_cp  = i_cp0-1+(njlm+lj+mj+1)+
     >                (nilm+li+mi)*paw_proj_nbasis(ia)
             i_qlm = i_qlm0 + ilm

             tmp_comp_pot = dbl_mb(paw_comp_pot_matrix(1)+i_mtr)*
     >               gen_gaunt_coeff(l,m,li,mi,lj,mj)*
     >               dcpl_mb(paw_qlm(1) + i_qlm)

                  dcpl_mb(paw_pot_comp(1)+i_cp) = 
     >                  dcpl_mb(paw_pot_comp(1)+i_cp)-
     >                  tmp_comp_pot


          end if
          end do !mi
          end do !mj
          ilm=ilm+1
        end do !m
        end do !l

       njlm = njlm + 2*lj+1
        end do !j
        nilm = nilm + 2*li+1
        end do !i  
      end do !in

 
c*     *** print comp potential ***
c      write(*,*) "comp pot"
c      do in=1,ion_nion()
c        ia = ion_katm(in)
c        i_cp0 = int_mb(i_paw_pot_comp(1) + in - 1)
c        nb = paw_basis_nbasis(ia)
c        nilm = 0
c        do i=1,nb
c          li = paw_basis_orb_l(i,ia)
c          do mi=-li,li
c            njlm = 0
c            do j=1,nb
c              lj = paw_basis_orb_l(j,ia)
c                do mj=-lj,lj
c                 i_cp  = i_cp0-1+(njlm+lj+mj+1)+
c     >                   (nilm+li+mi)*paw_proj_nbasis(ia)
c
c                  write(*,*) i,mi,j,mj,in,dcpl_mb(paw_pot_comp(1)+i_cp),
c     >                     i_cp
c                end do
c                njlm=njlm + 2*lj+1
c             end do
c           end do   
c           nilm = nilm + 2*li+1
c         end do
c      end do

      return
      end 

!!**************************************************
!!
!!      name:
!!
!!      purpose 
!!
!!      created:        7/9/2002
!!************************************************
!      subroutine paw_pot_comp_print()
!      
!      implicit none
!
!#include "bafdecls.fh"
!#include "paw_geom.fh"
!#include "paw_ma.fh"
!#include "paw_matrix_comp_pot.fh"
!#include "paw_basis.fh"
!#include "paw_proj.fh"
!#include "paw_pot_comp_data.fh"
!#include "paw_comp_charge_data.fh"
!
!      double precision gen_gaunt_coeff
!      external gen_gaunt_coeff
!      integer in,ia
!      integer i_mtr,i_mtr0
!      integer i,j,li,lj,mi,mj
!!      integer nb,nb2,nb3,nb4
!      integer l,m,ilm
!      integer mult_l
!      integer i_cp,i_cp0
!      integer nilm,njlm
!      integer i_qlm,i_qlm0
!      complex*16 tmp_comp_pot
!
!
!*     *** print comp potential ***
!      write(38,*) paw_pot_comp(3)
!      do in=1,ion_nion()
!        ia = ion_katm(in)
!        i_cp0 = int_mb(i_paw_pot_comp(1) + in - 1)
!        nb = paw_basis_nbasis(ia)
!        nilm = 0
!        do i=1,nb
!          li = paw_basis_orb_l(i,ia)
!          do mi=-li,li
!            njlm = 0
!            do j=1,nb
!              lj = paw_basis_orb_l(j,ia)
!                do mj=-lj,lj
!                 i_cp  = i_cp0-1+(njlm+lj+mj+1)+
!     >                   (nilm+li+mi)*paw_proj_nbasis(ia)
!
!                  write(38,*) i,mi,j,mj,in,
!!     >                        dcpl_mb(paw_pot_comp(1)+i_cp)
!                end do
!                njlm=njlm + 2*lj+1
!             end do
!           end do   
!           nilm = nilm + 2*li+1
!         end do
!      end do
!
!      return
!      end 
!

      subroutine paw_comp_pot_ptr(ptr)
      implicit none
      integer ptr

#include "paw_pot_comp_data.fh"

      ptr = paw_pot_comp(1)
      return
      end

