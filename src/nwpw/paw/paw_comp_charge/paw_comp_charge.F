*
* $Id: paw_comp_charge.F,v 1.16 2003-02-15 22:15:46 marat Exp $
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

#include "mafdecls.fh"
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
      ok = my_alloc(mt_int,nion,"i_paw_qlm",
     >                                i_paw_qlm)

      qlm_size = 0
      do in=1,nion
       int_mb(i_paw_qlm(1) + in - 1) = qlm_size
        ia = ion_katm(in)
        qlm_size = qlm_size
     >             + (paw_basis_mult_l(ia) + 1)**2
      end do

      !*** allocate comp. charge arrays ***
      ok = .TRUE.
      ok = ok .AND. my_alloc(MT_DCPL,qlm_size,"paw_qlm",paw_qlm)
      ok = ok .AND. my_alloc(MT_INT,qlm_size,"paw_qlm_l",paw_qlm_l)
      ok = ok .AND. my_alloc(MT_INT,qlm_size,"paw_qlm_m",paw_qlm_m)
      ok = ok .AND. my_alloc(MT_DBL,nkatm,"paw_q00",paw_q00)

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

      end

      subroutine paw_comp_charge_update()
        implicit none

#include "mafdecls.fh"
#include "paw_comp_charge_data.fh"
#include "paw_proj.fh"
#include "paw_basis.fh"
#include "paw_geom.fh"
#include "paw_comp_charge_matrix.fh"
#include "paw_ovlp_data.fh"


        integer na
        integer ia
        integer nlmi
        integer i1,i2
        integer l1,l2
        integer tot_natom
        integer nproj
        integer i_matr,i_matr0
        integer i_w,i_w0,i_w1
        integer m1,m2
        integer nb2,nb
        integer mult_l
        complex qlm_tmp
        
        integer i,j,li,lj,mi,mj
        integer l,m
        integer i_qlm
        integer nilm,njlm
        double precision gen_gaunt_coeff
        external gen_gaunt_coeff

        tot_natom = ion_nion()
    
c       !*** set values for l=0 m=0 elements ***
        nlmi = 0
        do na=1,tot_natom
          ia = ion_katm(na)
          dcpl_mb(paw_qlm(1)+nlmi) = dbl_mb(paw_q00(1)+ia-1)
          nlmi = nlmi + (paw_basis_mult_l(ia)+1)**2
        end do
 
        do na=1,tot_natom
          ia = ion_katm(na)
          mult_l = paw_basis_mult_l(ia)
          nb =paw_basis_nbasis(ia)
          nb2 = nb*nb
          i_matr0 = int_mb(i_comp_charge_matrix(1) + ia - 1)
          i_w0 = int_mb(i_paw_ovlp_w(1) + na - 1)
          nilm = 0
          do i=1,nb
          li = paw_basis_orb_l(i,ia)
          njlm = 0
          do j=1,nb
          lj = paw_basis_orb_l(j,ia)
            do i_qlm = 1,paw_qlm(3)
              qlm_tmp = dcmplx(0.0,0.0)
              l = int_mb(paw_qlm_l(1) + i_qlm - 1)
              m = int_mb(paw_qlm_m(1) + i_qlm - 1)
              if(l.le.(li+lj) .and. l.ge.abs(li-lj)) then
                do mi=-li,li
                 do mj=-lj,lj
                   if(m.eq.(mj-mi)) then
                   i_w   = i_w0-1+(njlm+lj+mj+1)+
     >                   (nilm+li+mi)*paw_proj_nbasis(ia)
                      
                   i_matr = i_matr0-1+i+(j-1)*nb+l*nb2
                   qlm_tmp = qlm_tmp + 
     >                     dbl_mb(comp_charge_matrix(1) + i_matr)* 
     >                     gen_gaunt_coeff(l,m,lj,mj,li,mi)*
     >                     dcpl_mb(paw_ovlp_w(1) + i_w) 
        write(*,*) "i,j,li,lj,mi,mj,l,m",i,j,li,lj,mi,mj,l,m
        write(*,*) "charge-matrix",
     >              dbl_mb(comp_charge_matrix(1) + i_matr) 
        write(*,*) "gaunt",
     >              gen_gaunt_coeff(l,m,lj,mj,li,mi)
        write(*,*) "overlap,i_w",
     >               i_w,dcpl_mb(paw_ovlp_w(1) + i_w)
        write(*,*) "qlm_tmp",
     >               qlm_tmp
                   end if
                 end do
                end do
              end if
              dcpl_mb(paw_qlm(1)+i_qlm-1) = dcpl_mb(paw_qlm(1)+i_qlm-1)+
     >                                      qlm_tmp
            end do !i_qlm
          njlm = njlm + 2*lj+1
          end do !j
          nilm = nilm + 2*li+1
          end do !i
        end do !na
        
c        nlmi = 0
c        do na=1,tot_natom
c          ia = ion_katm(na)
c          mult_l = paw_basis_mult_l(ia)
c          i_matr0 = int_mb(i_comp_charge_matrix(1) + ia - 1)
cc          i_w0 = int_mb(i_paw_ovlp_w(1) + ia - 1)
c          i_w0 = int_mb(i_paw_ovlp_w(1) + na - 1)
c          nb =paw_basis_nbasis(ia)
c          nb2 = nb*nb
c          nproj = paw_proj_nbasis(ia)
c          do l=0,mult_l  
c            do m=0,l
c              qlm_tmp = cmplx(0.0,0.0)
c              i_w1 = 0
c              do i1=1,nb
c                l1 = paw_basis_orb_l(i1,ia)
c                do m1=-l1,l1
c                  do i2=1,nb
c                    l2 = paw_basis_orb_l(i2,ia)
c                    do m2=-l2,l2
c                      i_matr = i_matr0 +i1+(i2-1)*nb+l*nb2
c                      i_w = i_w0 + i_w1
c                      qlm_tmp = qlm_tmp + 
c     >                     dbl_mb(comp_charge_matrix(1) + i_matr-1)* 
c     >                     gen_gaunt_coeff(l,m,l2,m2,l1,m1)*
c     >                     dcpl_mb(paw_ovlp_w(1) + i_w) 
c                      i_w1 = i_w1 + 1
c                     end do
c                   end do
c                 end do
c               end do
c               dcpl_mb(paw_qlm(1) + nlmi) = dcpl_mb(paw_qlm(1) + nlmi)
c     >                                  + qlm_tmp
c               nlmi = nlmi+1
c             end do
c           end do
c         end do 
c 
        write(*,*) "paw_comp_charge_coeff"
        nlmi = 0
        do na=1,tot_natom
          ia = ion_katm(na)
          mult_l = paw_basis_mult_l(ia)
          do l=0,mult_l  
          do m=-l,l
            write(*,*) l,m
            write(*,*) dcpl_mb(paw_qlm(1) + nlmi) 
            nlmi = nlmi+1
          end do
          end do
         end do 
      end !subroutine



      subroutine paw_comp_charge_end()
      implicit none

#include "paw_comp_charge_data.fh"
#include "paw_ma.fh"

      logical ok

      ok =        my_dealloc(paw_qlm)
      ok = ok.and.my_dealloc(paw_q00)
      if (.not.ok)
     > call errquit("paw_comp_charge_init: error freeing heap",0)

      return
      end !subroutine




