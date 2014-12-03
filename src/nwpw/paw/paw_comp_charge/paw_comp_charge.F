*
* $Id$
*


!**************************************************
!
!       Name: paw_proj_init
!
!       Purpose: initializes the paw projectors
!
!       Created:        7/30/2002
!**************************************************
      subroutine paw_comp_charge_init()
      implicit none

#include "bafdecls.fh"
#include "paw_comp_charge_data.fh"
#include "paw_geom.fh"
#include "paw_ma.fh"
#include "paw_params.fh"
#include "paw_basis.fh"


      logical ok
      integer nion
      integer nkatm
      integer qlm_size
      integer ia,in
      integer i
      integer l,m
        
      double precision q00_tmp

      !*** get total number of ions in the system ***
      nion = ion_nion()
        
      !*** get total number of diff kinds of atoms ***
      nion = ion_nion()
      nkatm = ion_nkatm()

      !*** calculate size of comp. charge array ***
      ok = my_alloc(mt_int,nion,"i_paw_qlm",i_paw_qlm)
      if (.not.ok)
     > call errquit("paw_comp_charge_init:out of heap memory",0,0)

      qlm_size = 0
      do in=1,nion
       int_mb(i_paw_qlm(1) + in - 1) = qlm_size
        ia = ion_katm(in)
        qlm_size = qlm_size
     >             + (paw_basis_mult_l(ia) + 1)**2
      end do

      !*** allocate comp. charge arrays ***
      ok =        my_alloc(MT_DCPL,qlm_size,"paw_qlm",paw_qlm)
      ok = ok.and.my_alloc(MT_DCPL,qlm_size,"paw_qlm1",
     >                     paw_qlm_spin(1,1))
      ok = ok.and.my_alloc(MT_DCPL,qlm_size,"paw_qlm2",
     >                     paw_qlm_spin(1,2))
      ok = ok.and.my_alloc(MT_INT,qlm_size,"paw_qlm_l",paw_qlm_l)
      ok = ok.and.my_alloc(MT_INT,qlm_size,"paw_qlm_m",paw_qlm_m)
      ok = ok.and.my_alloc(MT_DBL,nkatm,"paw_q00",paw_q00)
      if (.not.ok)
     > call errquit("paw_comp_charge_init:out of heap memory",0,1)

      i = 0
      do in=1,nion
         ia = ion_katm(in)
         do l=0,paw_basis_mult_l(ia)
           do m=-l,l
             int_mb(paw_qlm_l(1)+i)=l
             int_mb(paw_qlm_m(1)+i)=m 
             i = i + 1
           end do
         end do
      end do

      !*** set paw_q00 array ***
      do ia=1,nkatm
        q00_tmp = (paw_basis_core_charge(ia)-paw_basis_ion_charge(ia))
     >            /sqrt(fourpi)    
        dbl_mb(paw_q00(1)+ia-1)=q00_tmp
      end do

      return
      end


      subroutine paw_comp_charge_update()
        implicit none

#include "bafdecls.fh"
#include "paw_comp_charge_data.fh"
#include "paw_proj.fh"
#include "paw_basis.fh"
#include "paw_geom.fh"
#include "paw_comp_charge_matrix.fh"
#include "paw_ovlp_data.fh"


        integer na
        integer ia
        integer nlmi
        integer tot_natom
        integer i_matr,i_matr0
        integer i_w,i_w0
        integer nb2,nb
        integer mult_l
        complex*16 qlm_tmp
        
        integer i,j,li,lj,mi,mj
        integer l,m
        integer i_qlm,i_qlm0
        integer nilm,njlm
        double precision gen_gaunt_coeff
        external gen_gaunt_coeff

        tot_natom = ion_nion()

*       *** zero the comp charges ***
        do i=1,paw_qlm(3)
          dcpl_mb(paw_qlm(1)+i-1) = dcmplx(0.0d0,0.0d0)
        end do
    
c       !*** set values for l=0 m=0 elements ***
        nlmi = 0
        do na=1,tot_natom
          ia = ion_katm(na)
          dcpl_mb(paw_qlm(1)+nlmi) = dbl_mb(paw_q00(1)+ia-1)
          nlmi = nlmi + (paw_basis_mult_l(ia)+1)**2
        end do
c        go to 999
 
        do na=1,tot_natom
          ia = ion_katm(na)
          mult_l = paw_basis_mult_l(ia)
          nb =paw_basis_nbasis(ia)
          nb2 = nb*nb
          i_matr0 = int_mb(i_comp_charge_matrix(1) + ia - 1)
          i_w0 = int_mb(i_paw_ovlp_w(1) + na - 1)
          i_qlm0 = int_mb(i_paw_qlm(1) + na - 1)
          nilm = 0
          do i=1,nb
          li = paw_basis_orb_l(i,ia)
          njlm = 0
          do j=1,nb
          lj = paw_basis_orb_l(j,ia)
          nlmi = 0
          do l=0,mult_l
              qlm_tmp = dcmplx(0.0d0,0.0d0)
              if(l.le.(li+lj) .and. l.ge.abs(li-lj)) then
                do m=-l,l
                do mi=-li,li
                 do mj=-lj,lj
                   if(m.eq.(mj-mi)) then
                   i_w   = i_w0-1+(njlm+lj+mj+1)+
     >                   (nilm+li+mi)*paw_proj_nbasis(ia)
                      
                   i_qlm = i_qlm0-1+nlmi+l+m+1
                   i_matr = i_matr0-1+j+(i-1)*nb+l*nb2
                   dcpl_mb(paw_qlm(1)+i_qlm) = 
     >                     dcpl_mb(paw_qlm(1)+i_qlm)+
     >                     dbl_mb(comp_charge_matrix(1) + i_matr)* 
     >                     gen_gaunt_coeff(l,m,lj,mj,li,mi)*
     >                     dcpl_mb(paw_ovlp_w(1) + i_w) 
                  end if
                 end do
                end do
                end do !m
              end if
             nlmi = nlmi + 2*l+1
          end do !l
          njlm = njlm + 2*lj+1
          end do !j
          nilm = nilm + 2*li+1
          end do !i
        end do !na

 999    continue
        
 
c       write(*,*) "paw_comp_charge_coeff"
c       do na=1,tot_natom
c         i_qlm0 = int_mb(i_paw_qlm(1) + na - 1)
c         ia = ion_katm(na)
c         mult_l = paw_basis_mult_l(ia)
c         nlmi = 0
c         do l=0,mult_l  
c         do m=-l,l
c           i_qlm = i_qlm0-1+nlmi+l+m+1
c           write(*,*) nlmi+l+m+1,na,dcpl_mb(paw_qlm(1) + i_qlm) 
c         end do
c         nlmi = nlmi + 2*l+1
c         end do
c        end do 

      return
      end !subroutine



      subroutine paw_comp_charge_end()
      implicit none

#include "paw_comp_charge_data.fh"
#include "paw_ma.fh"

      logical ok

      ok =        my_dealloc(i_paw_qlm)
      ok = ok.and.my_dealloc(paw_qlm)
      ok = ok.and.my_dealloc(paw_qlm_spin(1,1))
      ok = ok.and.my_dealloc(paw_qlm_spin(1,2))
      ok = ok.and.my_dealloc(paw_qlm_l)
      ok = ok.and.my_dealloc(paw_qlm_m)
      ok = ok.and.my_dealloc(paw_q00)
      if (.not.ok)
     > call errquit("paw_comp_charge_end: error freeing heap",0,0)

      return
      end !subroutine




      subroutine paw_comp_charge_spin_update()
      implicit none

#include "bafdecls.fh"
#include "paw_comp_charge_data.fh"
#include "paw_proj.fh"
#include "paw_basis.fh"
#include "paw_geom.fh"
#include "paw_comp_charge_matrix.fh"
#include "paw_ovlp_data.fh"


        integer na
        integer ia
        integer nlmi
        integer tot_natom
        integer i_matr,i_matr0
        integer i_w,i_w0
        integer nb2,nb
        integer mult_l
        complex*16 qlm_tmp
        
        integer i,j,li,lj,mi,mj,ms
        integer l,m
        integer i_qlm,i_qlm0
        integer nilm,njlm
        integer ptr_qlm,ptr_w
        double precision occ
        double precision gen_gaunt_coeff
        external gen_gaunt_coeff

        tot_natom = ion_nion()

c        occ = 1.0d0
c        if (ispin.eq.1) occ = 0.5d0  ! rescal ovlp_weight if a restricted calculation
c                                     ! since paw_ovlp_w_u is multiplied by 2 if restricted

        do ms=1,ispin
           ptr_qlm = paw_qlm_spin(1,ms)
           ptr_w   = paw_ovlp_w_spin(1,ms)

*       *** zero the comp charges ***
        do i=1,paw_qlm_spin(3,ms)
          dcpl_mb(ptr_qlm+i-1) = dcmplx(0.0d0,0.0d0)
        end do
    
 
        do na=1,tot_natom
          ia = ion_katm(na)
          mult_l = paw_basis_mult_l(ia)
          nb =paw_basis_nbasis(ia)
          nb2 = nb*nb
          i_matr0 = int_mb(i_comp_charge_matrix(1) + ia - 1)
          i_w0 = int_mb(i_paw_ovlp_w(1) + na - 1)
          i_qlm0 = int_mb(i_paw_qlm(1) + na - 1)
          nilm = 0
          do i=1,nb
          li = paw_basis_orb_l(i,ia)
          njlm = 0
          do j=1,nb
          lj = paw_basis_orb_l(j,ia)
          nlmi = 0
          do l=0,mult_l
              qlm_tmp = dcmplx(0.0d0,0.0d0)
              if(l.le.(li+lj) .and. l.ge.abs(li-lj)) then
                do m=-l,l
                do mi=-li,li
                 do mj=-lj,lj
                   if(m.eq.(mj-mi)) then
                   i_w   = i_w0-1+(njlm+lj+mj+1)+
     >                   (nilm+li+mi)*paw_proj_nbasis(ia)
                      
                   i_qlm = i_qlm0-1+nlmi+l+m+1
                   i_matr = i_matr0-1+j+(i-1)*nb+l*nb2
                   dcpl_mb(ptr_qlm+i_qlm) 
     >               =  dcpl_mb(ptr_qlm+i_qlm)
     >               +     dbl_mb(comp_charge_matrix(1) + i_matr) 
     >                    *gen_gaunt_coeff(l,m,lj,mj,li,mi)
     >                    *dcpl_mb(ptr_w + i_w) 
c     >                    *occ
                  end if
                 end do
                end do
                end do !m
              end if
             nlmi = nlmi + 2*l+1
          end do !l
          njlm = njlm + 2*lj+1
          end do !j
          nilm = nilm + 2*li+1
          end do !i
        end do !na
        
 
      end do !*ms*

      return
      end !subroutine



