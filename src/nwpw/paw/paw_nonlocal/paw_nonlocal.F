*
* $Id$
*

* $Log: not supported by cvs2svn $
* Revision 1.15  2005/01/31 02:09:38  bylaska
* Large nested loops in PAW  parallelized using an openMP-like strategy.
* ...EJB
*
* Revision 1.14  2004/08/12 18:39:41  bylaska
* A prototype of a Grassmann CG paw minimizer (i.e. nwpw:minimizer 1) has been added.
* The code is similar to the CG minimizer in pspw, but differences exist
* because the residual |R> = (1 - S|psi><psi|)|Hpsi> is not the same as the
* tangent vector |T> = (1 - |psi><psi|S)|R>.
*
* Forces still need to be implemented.
*
* ...EJB
*
* Revision 1.13  2003/10/21 02:05:17  marat
* switched to new errquit by running global replace operation
* see the script below (note it will not work on multiline errquit calls)
* *********************************************************
* #!/bin/sh
*
* e=`find . -name "*F" -print`
*
* for f in $e
* do
* cp $f $f.bak
* sed  's|\(^[ ].*call[ ]*errquit([^,]*\)\(,[^,]*\)\()\)|\1,0\2\3|' $f.bak > $f
* #rm $f.bak
* done
* **********************************************************
*
* Revision 1.12  2003/03/25 01:11:30  edo
* uninitialized variable
*
* Revision 1.11  2003/03/14 01:20:10  marat
* added routines for fetching pointers
* to nonlocal matrices
* MV
*
* Revision 1.10  2003/03/05 23:16:32  bylaska
* Commented out write statements and other minor fixes.....
* self-consistent loop looks like it is working.....
* ....EJB
*
* Revision 1.9  2003/03/05 20:36:17  bylaska
* multipole potential added...Eigenvalues agree between this code and
* F90 code....EJB
*
* Revision 1.8  2003/03/05 01:46:58  bylaska
* bug fixes....The eigenvalues now agree with the F90 code for all terms
* except the multipole potential....I have currently commented out
* the multipole potential term....EJB
*

!**************************************************
!
!       Name: paw_nonlocal_init
!
!       Purpose: initializes the paw nonlocal operator
!
!       Created:        2/25/2003
!**************************************************
      subroutine paw_nonlocal_init(ispin0,ne)
      implicit none
      integer ispin0
      integer ne(2)

#include "bafdecls.fh"
#include "paw_proj.fh"
#include "paw_nonlocal_data.fh" 
#include "paw_matrix_kin_data.fh"
#include "paw_geom.fh" 
#include "paw_ma.fh" 
#include "paw_basis.fh" 
  
*     **** local variables ****
      logical value
      integer ii,ia,nion
      integer i,li,mi,j,lj,mj,il,jl
      integer basis_nbasis,indx,indx_m,ptr
      integer ke_ptr,ion_ptr,core_ptr,vloc_ptr

      nion = ion_nion()

      !*** calculate total number of (n,l,m) projectors  - Yet Again ***
      total_nbasis = 0
      do ia=1,ion_nkatm()
         total_nbasis = total_nbasis 
     >                + paw_proj_nbasis(ia)*ion_natm(ia)
      end do

      !*** calculate the size and indexing for Gop  - Yet Again ***
      value = my_alloc(mt_int,nion,"i_paw_Gop",i_paw_Gop)
      if (.not.value) call errquit('paw_nonlocal_init: alloc heap',0,1)

      paw_Gop_size = 0
      do ii=1,nion
         int_mb(i_paw_Gop(1)+ii-1) = paw_Gop_size
         ia = ion_katm(ii)
         paw_Gop_size = paw_Gop_size
     >                + paw_proj_nbasis(ia)**2
      end do


      !*** allocate space for the non-local operator ***
      !-- total number of electron states
      ispin = ispin0
      neall = ne(1) + ne(2)
      ne_u=ne(1)
      ne_d=ne(2)

      !*** space for the non-local operator ***
      value = value.and.my_alloc(mt_dcpl,(total_nbasis*neall),
     >                          "paw_nl_coeff",paw_nl_coeff)
      value = value.and.my_alloc(mt_dcpl,paw_Gop_size,
     >                           "paw_Gop_u",paw_Gop_spin(1,1))
      if (ispin.eq.2)
     > value = value.and.my_alloc(mt_dcpl,paw_Gop_size,
     >                           "paw_Gop_d",paw_Gop_spin(1,2))
      value = value.and.my_alloc(mt_dcpl,paw_Gop_size,
     >                           "paw_Gop_static",paw_Gop_static)
      if (.not.value) call errquit('paw_nonlocal_init: alloc heap',0,1)

      do ii=1,paw_Gop_static(3)
         dcpl_mb(paw_Gop_static(1)+ii-1)=dcmplx(0.0,0.0)
      end do

      !*** set the static part of the nonlocal operator ****
      call paw_kin_matrix_ptr(ke_ptr)
      call paw_ion_matrix_ptr(ion_ptr)
      call paw_core_matrix_ptr(core_ptr)
      call paw_vloc_matrix_ptr(vloc_ptr)
      indx_m = 0
      do ii=1,nion
         ia           = ion_katm(ii)
         basis_nbasis = paw_basis_nbasis(ia)
         ptr = int_mb(i_paw_kin_matrix(1)+ia-1)

         jl = 0
         do j=1,basis_nbasis
           lj = paw_basis_orb_l(j,ia)
           do mj=-lj,lj

             il = 0
             do i=1,basis_nbasis
               li = paw_basis_orb_l(i,ia)
               do mi=-li,li
              
                if(mi.eq.mj) then
                 indx = ptr + (j-1)*basis_nbasis + (i-1)
                 dcpl_mb(paw_Gop_static(1)+indx_m) 
     >               = dbl_mb(ke_ptr  +indx)
     >               + dbl_mb(ion_ptr +indx)
     >               + dbl_mb(core_ptr+indx)
     >               - dbl_mb(vloc_ptr+indx)
                 end if

                 indx_m = indx_m + 1
               end do !*mi*
             end do !*i*

           end do !*mj*
         end do !*j*

      end do !*ii*

      return
      end

!**************************************************
!
!       Name: paw_nonlocal_end
!
!       Purpose: removes space used by the paw nonlocal
!
!       Created:        8/07/2002
!**************************************************
      subroutine paw_nonlocal_end()
      implicit none      

#include "bafdecls.fh"
#include "paw_ma.fh" 
#include "paw_nonlocal_data.fh"

      !*** local variables ***
      logical value
      value=.true.
      value = value.and.my_dealloc(i_paw_Gop)

      value = value.and.my_dealloc(paw_Gop_spin(1,1))
      if (ispin.eq.2) value = value.and.my_dealloc(paw_Gop_spin(1,2))

      value = value.and.my_dealloc(paw_Gop_static)
      value = value.and.my_dealloc(paw_nl_coeff)
      if (.not.value) call errquit('paw_nonlocal_end: dealloc heap',0,0)
      return
      end


!**************************************************
!
!       Name: paw_nonlocal_gen_nl_coeff
!
!       Purpose: 
!
!       Created:        8/07/2002
!**************************************************
      subroutine paw_nonlocal_gen_nl_coeff()
      implicit none      

#include "bafdecls.fh"
#include "paw_nonlocal_data.fh"
#include "paw_proj.fh"
#include "paw_basis.fh" 
#include "paw_geom.fh" 

      !*** local variables ***
      integer taskid,np,tcount
      integer k
      integer ii,ia,n,nion,nbasis,nlm,nlm1
      integer i_nbasis,indx,indx1,indx_m
      integer coeff_ptr
      integer hartree_ptr,mult_ptr,comp_ptr,xc_ptr
      integer n1(2),n2(2),ms
      real*8  scal1

      call Parallel_np(np)
      call Parallel_taskid(taskid)

      !**** add up Gop ****
      call paw_hartree_pot_ptr(hartree_ptr)
      call paw_mult_pot_ptr(mult_ptr)
      call paw_comp_pot_ptr(comp_ptr)


      call dcopy(2*paw_Gop_size,
     >           dcpl_mb(paw_Gop_static(1)),1,
     >           dcpl_mb(paw_Gop_spin(1,1)),1)
      call daxpy(2*paw_Gop_size,1.0d0,
     >           dcpl_mb(hartree_ptr),1,
     >           dcpl_mb(paw_Gop_spin(1,1)),1)

      call daxpy(2*paw_Gop_size,1.0d0,
     >           dcpl_mb(mult_ptr),1,
     >           dcpl_mb(paw_Gop_spin(1,1)),1)

      call daxpy(2*paw_Gop_size,1.0d0,
     >           dcpl_mb(comp_ptr),1,
     >           dcpl_mb(paw_Gop_spin(1,1)),1)

      if (ispin.eq.2) then
      call dcopy(2*paw_Gop_size,
     >           dcpl_mb(paw_Gop_spin(1,1)),1,
     >           dcpl_mb(paw_Gop_spin(1,2)),1)
      call paw_xc_pot_down_ptr(xc_ptr)
      call daxpy(2*paw_Gop_size,1.0d0,
     >           dcpl_mb(xc_ptr),1,
     >           dcpl_mb(paw_Gop_spin(1,2)),1)
      end if

      call paw_xc_pot_up_ptr(xc_ptr)
      call daxpy(2*paw_Gop_size,1.0d0,
     >           dcpl_mb(xc_ptr),1,
     >           dcpl_mb(paw_Gop_spin(1,1)),1)




      !**** paw_nl_coeff(alpha,ii) = Sum(ii,beta) Gop(alpha,beta,ii)*coeff(beta,ii) ****
      n1(1) = 1
      n2(1) = ne_u
      n1(2) = ne_u+1
      n2(2) = ne_u+ne_d

      nion = ion_nion()
      call paw_ovlp_coeff_ptr(coeff_ptr)
      call dcopy(2*paw_nl_coeff(3),0.0d0,0,
     >           dcpl_mb(paw_nl_coeff(1)),1)
       
      tcount = 0
      i_nbasis = 0
      indx_m = 0
      do ii=1,nion
         ia = ion_katm(ii)
         nbasis = paw_proj_nbasis(ia)

         do nlm =1,nbasis
         do nlm1=1,nbasis
         if (tcount.eq.taskid) then

           !*** up and down spin ***
           do ms=1,ispin
           do n=n1(ms),n2(ms)
             indx  = (i_nbasis+nlm-1) *neall + (n-1)
             indx1 = (i_nbasis+nlm1-1)*neall + (n-1)

             dcpl_mb(paw_nl_coeff(1)+indx) 
     >           = dcpl_mb(paw_nl_coeff(1)+indx) 
     >           + dcpl_mb(paw_Gop_spin(1,ms)+indx_m)
     >            *dcpl_mb(coeff_ptr+indx1)
           end do
           end do


         end if
         tcount = mod(tcount+1,np)
         indx_m = indx_m + 1
         end do
         end do
         i_nbasis = i_nbasis + nbasis
      end do
      if (np.gt.1) call D3dB_Vector_Sumall(2*paw_nl_coeff(3),
     >                              dcpl_mb(paw_nl_coeff(1)))

      scal1 = -1.0d0
      call dscal(2*paw_nl_coeff(3),scal1,
     >           dcpl_mb(paw_nl_coeff(1)),1)

      return
      end

      

!**************************************************
!
!       Name: paw_nonlocal
!
!       Purpose: 
!
!       Created: Wed Feb 12 15:18:24 PST 2003
!**************************************************
      subroutine paw_nonlocal(ispin0,ne0,psi,Hpsi,move,fion)
      implicit none
      integer ispin0,ne0(2)
      complex*16 psi(*)
      complex*16 Hpsi(*)
      logical    move
      real*8     fion(3,*)

#include "bafdecls.fh"
#include "paw_nonlocal_data.fh"
#include "paw_proj.fh"
#include "paw_basis.fh" 
#include "paw_geom.fh" 

      !**** local variables ***
      logical ok
      integer ii,ia,n,nbasis,nion,prj,npack1
      integer nlm,nlmii,indx_prj,indx_psi,indx_coeff
      integer exi(2),p(2)

      call paw_nonlocal_gen_nl_coeff()

      call Pack_npack(1,npack1)
      nion  = ion_nion()

      ok = BA_push_get(mt_dcpl,npack1,'exi', exi(2), exi(1))
      ok = ok.and.
     >     BA_push_get(mt_dcpl,npack1,'p',p(2),p(1))
      if (.not. ok) call errquit('paw_nonlocal:out of stack memory',0,0)

      nlmii = 0
      do ii=1,nion
        ia     = ion_katm(ii)
        nbasis = paw_proj_nbasis(ia)
        prj    = paw_proj_i_prj_atom(ia)
        call strfac_pack(1,ii,dcpl_mb(exi(1)))


        do nlm=1,nbasis
           nlmii = nlmii + 1
           
           indx_prj = (nlm-1)*npack1
           call Pack_cc_Mul(1,
     >                      dcpl_mb(exi(1)),
     >                      dcpl_mb(prj+indx_prj),
     >                      dcpl_mb(p(1)))


           do n=1,neall
             indx_psi   = (n-1)*npack1 + 1
             indx_coeff = paw_nl_coeff(1) 
     >                  + (nlmii-1)*neall 
     >                  + (n-1) 

             call Pack_cc_zaxpy(1,
     >                          dcpl_mb(indx_coeff),
     >                          dcpl_mb(p(1)),
     >                          Hpsi(indx_psi))
           end do

        end do

      end do
      
      ok =        BA_pop_stack(p(2))
      ok = ok.and.BA_pop_stack(exi(2))
      if (.not. ok) 
     > call errquit('paw_nonlocal:error popping stack ',0,2)

      return
      end

!      subroutine paw_Gop_print()
!        
!        implicit none
!
!#include "paw_geom.fh"
!#include "paw_basis.fh"
!#include "paw_proj.fh"
!#include "paw_nonlocal_data.fh"
!#include "bafdecls.fh"
!
!      double precision gen_gaunt_coeff
!      external gen_gaunt_coeff
!      integer i,j
!      integer in,ia
!      integer i_t,i_p
!      integer nion 
!      integer l,m
!      integer ng,ic,ig
!      integer i_rho0,i_rho
!      integer n_s
!      integer i_ylm,i_ylm1
!      integer nlmi
!      integer i_core_ps,i_core_ae
!      integer i_vxc0,i_vxc
!      integer i_pot0,i_pot
!      integer i_r
!      integer nb,nb2
!      integer i_phi_ae0,i_phi_ps0
!      integer i_mtr0,i_mtr
!      integer li,lj,mi,mj
!      integer nilm,njlm
!      double precision exc_tmp
!      double precision log_amesh
!      double precision tmp_ae,tmp_ps
!      complex*16 tmp
!      complex*16 ylm
!
!      nion = ion_nion() 
!
!      write(45,*) paw_Gop_u(3)
!
!      do in=1,nion
!        ia = ion_katm(in)
!        ng = paw_basis_ngrid(ia)
!        ic = paw_basis_i_cut(ia)
!        i_r = paw_basis_i_rgrid(ia)
!        log_amesh = paw_basis_log_amesh(ia)
!        nb = paw_basis_nbasis(ia)
!  
!*       *** print out exc potenetial ***
!        i_pot0 = int_mb(i_paw_Gop(1) + in - 1)
!        nilm = 0
!        do i=1,nb
!        li = paw_basis_orb_l(i,ia)
!        njlm = 0
!        do j=1,nb
!        lj = paw_basis_orb_l(j,ia)
!          do mi=-li,li
!          do mj=-lj,lj
!            i_pot = i_pot0-1+(njlm+lj+mj+1)+
!     >                (nilm+li+mi)*paw_proj_nbasis(ia)
!
!
!
!            write(45,*) i,mi,j,mj,in,
!     >                dcpl_mb(paw_Gop_u(1)+i_pot),       
!     >                dcpl_mb(paw_Gop_d(1)+i_pot)       
!             
!            
!
!            end do !mi
!          end do !mj
!        njlm = njlm + 2*lj+1
!        end do !j
!        nilm = nilm + 2*li+1
!        end do !i  
! 
!      end do !in
!
!      end
!

**************************************************
!
!       Name: 
!
!       Purpose: 
!
!       Created:        2/10/2002
!**************************************************
      subroutine paw_nonloc_Gop_u_ptr(ptr)
      implicit none
      integer ptr(4)

#include "paw_nonlocal_data.fh"

      integer i

      do i=1,4
      ptr(i) = paw_Gop_spin(i,1)
      end do
      return
      end 

**************************************************
!
!       Name: 
!
!       Purpose: 
!
!       Created:        2/10/2002
!**************************************************
      subroutine paw_nonloc_Gop_d_ptr(ptr)
      implicit none
      integer ptr(4)

#include "paw_nonlocal_data.fh"

      integer i

      do i=1,4
      ptr(i) = paw_Gop_spin(i,2)
      end do
      return
      end 


