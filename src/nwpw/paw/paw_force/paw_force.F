*
* $Id: paw_force.F,v 1.2 2003-03-11 00:52:50 marat Exp $
*

!**************************************************
!
!	Name: paw_ovlp_init
!
!	Purpose: initializes the paw overlap coefficients
!
!	Created:	8/07/2002
!**************************************************
      subroutine paw_force_init(ispin0,ne)
      implicit none
      integer ispin0
      integer ne(2)

#include "mafdecls.fh"
#include "paw_proj.fh"
#include "paw_force_data.fh" 
#include "paw_geom.fh" 
#include "paw_ma.fh" 
#include "paw_basis.fh" 
#include "paw_ovlp.fh" 
      
  
*     **** local variables ****
      logical value
      integer ii,ia,nlm,nlmii,n,m,mabs
      integer nbasis,i2,i1,isgn,nion
      integer paw_force_w_size
      integer paw_force_c_size
      integer ik
      integer in,li,lj,mi,mj,i,j


      nion = ion_nion()

      paw_force_c_size = paw_ovlp_coeff_size()

      value = value.and.my_alloc(mt_dcpl,paw_force_c_size,
     >                               'paw_force_cx_u',paw_force_cx_u)
      value = value.and.my_alloc(mt_dcpl,paw_force_c_size,
     >                               'paw_force_cx_s',paw_force_cx_s)
      value = value.and.my_alloc(mt_dcpl,paw_force_c_size,
     >                               'paw_force_cy_u',paw_force_cy_u)
      value = value.and.my_alloc(mt_dcpl,paw_force_c_size,
     >                               'paw_force_cy_s',paw_force_cy_s)
      value = value.and.my_alloc(mt_dcpl,paw_force_c_size,
     >                               'paw_force_cz_u',paw_force_cz_u)
      value = value.and.my_alloc(mt_dcpl,paw_force_c_size,
     >                               'paw_force_cz_s',paw_force_cz_s)



      paw_force_w_size = paw_ovlp_weight_size()
      value = value.and.my_alloc(mt_dcpl,paw_force_w_size,
     >                           "paw_force_wx",paw_force_wx)
      value = value.and.my_alloc(mt_dcpl,paw_force_w_size,
     >                           "paw_force_wx_u",paw_force_wx_u)
      value = value.and.my_alloc(mt_dcpl,paw_force_w_size,
     >                           "paw_force_wx_d",paw_force_wx_d)
      value = value.and.my_alloc(mt_dcpl,paw_force_w_size,
     >                           "paw_force_wy",paw_force_wy)
      value = value.and.my_alloc(mt_dcpl,paw_force_w_size,
     >                           "paw_force_wy_u",paw_force_wy_u)
      value = value.and.my_alloc(mt_dcpl,paw_force_w_size,
     >                           "paw_force_wy_d",paw_force_wy_d)
      value = value.and.my_alloc(mt_dcpl,paw_force_w_size,
     >                           "paw_force_wz",paw_force_wz)
      value = value.and.my_alloc(mt_dcpl,paw_force_w_size,
     >                           "paw_force_wz_u",paw_force_wz_u)
      value = value.and.my_alloc(mt_dcpl,paw_force_w_size,
     >                           "paw_force_wz_d",paw_force_wz_d)

      if (.not.value) call errquit('paw_force_init: alloc heap',1)

      return
      end

**************************************************
!
!	Name: paw_ovlp_coeff_set
!
!	Purpose: sets the paw overlap coefficients
!
!	Created:	8/07/2002
!**************************************************
      subroutine paw_force_solve(psi)
      implicit none
      complex*16 psi(*)
      
#include "mafdecls.fh"
#include "paw_force_data.fh"
#include "paw_proj.fh"
#include "paw_ovlp.fh"
#include "paw_geom.fh"
#include "paw_basis.fh"

      !*** local variables ***
      logical value
      integer ii,ia,nbasis,nlm,nlmii,npack1
      integer nion,prj,i1,i2,n
      integer indx_psi,indx_prj,indx_coeff
      integer exi(2),tmp(2)
      integer tmpx(2),tmpy(2),tmpz(2)
      integer ni,i,li,mi
      integer Gx(2),Gy(2),Gz(2)
      integer nfft3d
      integer indx2(4),sgn(4),coeff(4)
      integer ne_u,ne_d,neall
      double precision occup
      complex*16 wx_tmp,wy_tmp,wz_tmp
      integer i_nbasis,in,nlm1,indx_coeff1

 
      !**** external functions ****
      integer  control_ncut,G_indx
      real*8   control_rcut,double_factorial,lattice_omega
      real*8   lattice_unita
      external control_ncut,G_indx
      external control_rcut,double_factorial,lattice_omega
      external lattice_unita


      ne_u  = paw_ne_up()
      ne_d  = paw_ne_down()
      neall = paw_ne_total()
      call paw_ovlp_indx2_ptr(indx2)
      call paw_ovlp_sgn_ptr(sgn)
      call paw_ovlp_c_ptr(coeff)

      !*** allocate temporary space ***
      call D3dB_nfft3d(1,nfft3d)
      call Pack_npack(1,npack1)

      value =           MA_push_get(mt_dcpl,npack1,'exi',exi(2),exi(1))
      value = value.and.MA_push_get(mt_dcpl,npack1,'tmp',tmp(2),tmp(1))
      value = value.and.MA_push_get(mt_dcpl,npack1,'tmpx',
     >                   tmpx(2),tmpx(1))
      value = value.and.MA_push_get(mt_dcpl,npack1,'tmpy',
     >                   tmpy(2),tmpy(1))
      value = value.and.MA_push_get(mt_dcpl,npack1,'tmpz',
     >                   tmpz(2),tmpz(1))

      value = value.and. MA_push_get(mt_dbl,nfft3d,'Gx',Gx(2),Gx(1))
      value = value.and. MA_push_get(mt_dbl,nfft3d,'Gy',Gy(2),Gy(1))
      value = value.and. MA_push_get(mt_dbl,nfft3d,'Gz',Gz(2),Gz(1))
      if (.not.value) call errquit('paw_ovlp_coeff_set: push stack',0)

c    *** set reciprocal lattice arrays ***
      call D3dB_t_Copy(1,dbl_mb(G_indx(1)),dbl_mb(Gx(1)))
      call D3dB_t_Copy(1,dbl_mb(G_indx(2)),dbl_mb(Gy(1)))
      call D3dB_t_Copy(1,dbl_mb(G_indx(3)),dbl_mb(Gz(1)))
      call Pack_t_pack(1,dbl_mb(Gx(1)))
      call Pack_t_pack(1,dbl_mb(Gy(1)))
      call Pack_t_pack(1,dbl_mb(Gz(1)))


c     !*** compute overlap coefficients c(n=1:neall,nlm=1:nbasis(ii),ii=1:nion) ***
      nlmii = 0
      nion  = ion_nion()
      do ii=1,nion
        ia = ion_katm(ii)
        call strfac_pack(1,ii,dcpl_mb(exi(1)))
        nbasis = paw_proj_nbasis(ia)
        prj = paw_proj_i_prj_atom(ia)
        do nlm=1,nbasis
           
           nlmii = nlmii + 1
           indx_prj = (nlm-1)*npack1
           call Pack_cc_Mul(1,
     >                      dcpl_mb(exi(1)),
     >                      dcpl_mb(prj+indx_prj),
     >                      dcpl_mb(tmp(1)))

           call Pack_tc_Mul(1,
     >                      dbl_mb(Gx(1)),
     >                      dcpl_mb(tmp(1)),
     >                      dcpl_mb(tmpx(1)))

           call Pack_tc_Mul(1,
     >                      dbl_mb(Gy(1)),
     >                      dcpl_mb(tmp(1)),
     >                      dcpl_mb(tmpy(1)))

           call Pack_tc_Mul(1,
     >                      dbl_mb(Gz(1)),
     >                      dcpl_mb(tmp(1)),
     >                      dcpl_mb(tmpz(1)))

            do n=1,neall
             indx_psi   = (n-1)*npack1    + 1
             indx_coeff = (nlmii-1)*neall + (n-1)
              call Pack_cc_izdot(1,
     >                         dcpl_mb(tmpx(1)),
     >                         psi(indx_psi),
     >                         dcpl_mb(paw_force_cx_s(1)+indx_coeff))
              call Pack_cc_izdot(1,
     >                         dcpl_mb(tmpy(1)),
     >                         psi(indx_psi),
     >                         dcpl_mb(paw_force_cy_s(1)+indx_coeff))
              call Pack_cc_izdot(1,
     >                         dcpl_mb(tmpz(1)),
     >                         psi(indx_psi),
     >                         dcpl_mb(paw_force_cz_s(1)+indx_coeff))
           end do
        end do
      end do 

      call D3dB_Vector_SumAll(2*paw_force_cx_s(3),
     >                        dcpl_mb(paw_force_cy_s(1)))
      call D3dB_Vector_SumAll(2*paw_force_cy_s(3),
     >                        dcpl_mb(paw_force_cx_s(1)))
      call D3dB_Vector_SumAll(2*paw_force_cz_s(3),
     >                        dcpl_mb(paw_force_cz_s(1)))


      !*** unsrcramble to get true overlap coefficients ***
      do i1=1,paw_force_cx_s(3)
         i2 = int_mb(indx2(1)+i1-1)
         dcpl_mb(paw_force_cx_u(1)+i1-1) 
     >     = 0.5d0*( dcpl_mb(paw_force_cx_s(1)+i1-1)
     >             + int_mb(sgn(1)+i1-1)
     >              *dconjg(dcpl_mb(paw_force_cx_s(1)+i2)))
          dcpl_mb(paw_force_cy_u(1)+i1-1) 
     >     = 0.5d0*( dcpl_mb(paw_force_cy_s(1)+i1-1)
     >             + int_mb(sgn(1)+i1-1)
     >              *dconjg(dcpl_mb(paw_force_cy_s(1)+i2)))
          dcpl_mb(paw_force_cz_u(1)+i1-1) 
     >     = 0.5d0*( dcpl_mb(paw_force_cz_s(1)+i1-1)
     >             + int_mb(sgn(1)+i1-1)
     >              *dconjg(dcpl_mb(paw_force_cz_s(1)+i2)))
      end do     




      value = value.and.MA_pop_stack(Gz(2))
      value = value.and.MA_pop_stack(Gy(2))
      value = value.and.MA_pop_stack(Gx(2))
      value = value.and.MA_pop_stack(tmpz(2))
      value = value.and.MA_pop_stack(tmpy(2))
      value = value.and.MA_pop_stack(tmpx(2))
      value = value.and.MA_pop_stack(exi(2))
      if (.not.value) call errquit('paw_ovlp_coeff_set: pop stack',1)
      
c     *** calculate weights ***      

      if (paw_ispin().eq.1) then
        occup = 2.0d0
      else
        occup = 1.0d0
      end if

      ii = 0
      i_nbasis = 0
      do in=1,ion_nion()
        ia = ion_katm(in)
        nbasis = paw_proj_nbasis(ia)
        do nlm =1,nbasis
        do nlm1=1,nbasis

           wx_tmp = dcmplx(0.0d0,0.0d0)
           wy_tmp = dcmplx(0.0d0,0.0d0)
           wz_tmp = dcmplx(0.0d0,0.0d0)
           do n=1,ne_u

             indx_coeff  = (i_nbasis+nlm-1)*neall + (n-1)
             indx_coeff1 = (i_nbasis+nlm1-1)*neall + (n-1)

             wx_tmp = wx_tmp + conjg(dcpl_mb(coeff(1)+indx_coeff))*
     >        dcpl_mb(paw_force_cx_u(1)+indx_coeff1)+
     >        + conjg(dcpl_mb(coeff(1)+indx_coeff))*
     >        dcpl_mb(paw_force_cx_u(1)+indx_coeff1)
             
              wy_tmp = wy_tmp + conjg(dcpl_mb(coeff(1)+indx_coeff))*
     >        dcpl_mb(paw_force_cy_u(1)+indx_coeff1)+
     >        + conjg(dcpl_mb(coeff(1)+indx_coeff))*
     >        dcpl_mb(paw_force_cy_u(1)+indx_coeff1)
             
              wz_tmp = wz_tmp + conjg(dcpl_mb(coeff(1)+indx_coeff))*
     >        dcpl_mb(paw_force_cz_u(1)+indx_coeff1)+
     >        + conjg(dcpl_mb(coeff(1)+indx_coeff))*
     >        dcpl_mb(paw_force_cz_u(1)+indx_coeff1)
             
           end do
           dcpl_mb(paw_force_wx_u(1)+ii) = wx_tmp*occup
           dcpl_mb(paw_force_wy_u(1)+ii) = wy_tmp*occup
           dcpl_mb(paw_force_wz_u(1)+ii) = wz_tmp*occup


           !*** this loop is not computed when ispin=1 ***
           wx_tmp = dcmplx(0.0d0,0.0d0)
           wy_tmp = dcmplx(0.0d0,0.0d0)
           wz_tmp = dcmplx(0.0d0,0.0d0)
           do n=ne_u+1,ne_d+ne_u
             indx_coeff  = (i_nbasis+nlm-1)*neall + (n-1)
             indx_coeff1 = (i_nbasis+nlm1-1)*neall + (n-1)

             wx_tmp = wx_tmp + conjg(dcpl_mb(coeff(1)+indx_coeff))*
     >        dcpl_mb(paw_force_cx_u(1)+indx_coeff1)+
     >        + conjg(dcpl_mb(coeff(1)+indx_coeff))*
     >        dcpl_mb(paw_force_cx_u(1)+indx_coeff1)
             
              wy_tmp = wy_tmp + conjg(dcpl_mb(coeff(1)+indx_coeff))*
     >        dcpl_mb(paw_force_cy_u(1)+indx_coeff1)+
     >        + conjg(dcpl_mb(coeff(1)+indx_coeff))*
     >        dcpl_mb(paw_force_cy_u(1)+indx_coeff1)
             
              wz_tmp = wz_tmp + conjg(dcpl_mb(coeff(1)+indx_coeff))*
     >        dcpl_mb(paw_force_cz_u(1)+indx_coeff1)+
     >        + conjg(dcpl_mb(coeff(1)+indx_coeff))*
     >        dcpl_mb(paw_force_cz_u(1)+indx_coeff1)
 
           end do
           dcpl_mb(paw_force_wx_d(1)+ii) = wx_tmp*occup
           dcpl_mb(paw_force_wy_d(1)+ii) = wy_tmp*occup
           dcpl_mb(paw_force_wz_d(1)+ii) = wz_tmp*occup

           dcpl_mb(paw_force_wx(1)+ii) = dcpl_mb(paw_force_wx_u(1)+ii)+
     >                                 dcpl_mb(paw_force_wx_d(1)+ii)    
           dcpl_mb(paw_force_wy(1)+ii) = dcpl_mb(paw_force_wy_u(1)+ii)+
     >                                 dcpl_mb(paw_force_wy_d(1)+ii)    
           dcpl_mb(paw_force_wz(1)+ii) = dcpl_mb(paw_force_wz_u(1)+ii)+
     >                                 dcpl_mb(paw_force_wz_d(1)+ii)    
           
           ii=ii +1
        end do
        end do
        i_nbasis = i_nbasis + nbasis
      end do

      return
      end

c!**************************************************
c!
c!	Name: paw_ovlp_end
c!
c!	Purpose: removes space used by the paw overlap coefficients
c!
c!	Created:	8/07/2002
c!**************************************************
c      subroutine paw_ovlp_end()
c      implicit none      
c
c#include "mafdecls.fh"
c#include "paw_ma.fh" 
c#include "paw_ovlp_data.fh"
c
c      !*** local variables ***
c      logical value
c      
c      value =           MA_free_heap(coeff(2))
c      value = value.and.MA_free_heap(toeff(2))
c      value = value.and.MA_free_heap(indx2(2))
c      value = value.and.MA_free_heap(sgn(2))
c      value = value.and.MA_free_heap(paw_ovlp_w(2))
c      value = value.and.MA_free_heap(paw_ovlp_w_u(2))
c      value = value.and.MA_free_heap(paw_ovlp_w_d(2))
c      value = value.and.my_dealloc(i_paw_ovlp_w)
c
c      if (.not.value) call errquit('paw_ovlp_end: dealloc heap',0)
c      return
c      end
c
c
 
c**************************************************
c!
c!	Name: paw_ovlp_weigths_set
c!
c!	Purpose: sets the overlap weight coefficients
c!
c!	Created:	8/07/2002
c!**************************************************
c      subroutine paw_ovlp_weights_set()
c      
c      implicit none
c
c#include "mafdecls.fh"
c#include "paw_geom.fh"
c#include "paw_ovlp_data.fh"
c#include "paw_proj.fh"
c#include "paw_basis.fh"
c
c      integer i_nbasis,ii,ii1
c      integer in,nlm,nlm1,ia,n
c      integer nbasis
c      integer indx_coeff,indx_coeff1
c      integer i,j,mi,mj,li,lj
c      integer i_w0,i_w
c      integer nilm,njlm,niplm,njplm
c      real*8  occup
c      double complex w_tmp
c      
c      
c      !*** check for restricted calculation ****
c      occup = 1.0d0
c      if (ispin.eq.1) then
c        occup = 2.0d0
c      end if
c
c      ii = 0
c      i_nbasis = 0
c      do in=1,ion_nion()
c        ia = ion_katm(in)
c        nbasis = paw_proj_nbasis(ia)
c        do nlm =1,nbasis
c        do nlm1=1,nbasis
c
c           w_tmp = dcmplx(0.0d0,0.0d0)
c           do n=1,ne_u
c             indx_coeff  = (i_nbasis+nlm-1)*neall + (n-1)
c             indx_coeff1 = (i_nbasis+nlm1-1)*neall + (n-1)
c             w_tmp = w_tmp + conjg(dcpl_mb(coeff(1)+indx_coeff))*
c     >        dcpl_mb(coeff(1)+indx_coeff1)
c           end do
c           dcpl_mb(paw_ovlp_w_u(1)+ii) = w_tmp*occup
c
c
c           !*** this loop is not computed when ispin=1 ***
c           w_tmp = dcmplx(0.0d0,0.0d0)
c           do n=ne_u+1,ne_d+ne_u
c             indx_coeff  = (i_nbasis+nlm-1)*neall + (n-1)
c             indx_coeff1 = (i_nbasis+nlm1-1)*neall + (n-1)
c             w_tmp = w_tmp + conjg(dcpl_mb(coeff(1)+indx_coeff))*
c     >        dcpl_mb(coeff(1)+indx_coeff1)
c           end do
c           dcpl_mb(paw_ovlp_w_d(1)+ii) = w_tmp
c
c           dcpl_mb(paw_ovlp_w(1)+ii) = dcpl_mb(paw_ovlp_w_u(1)+ii)+
c     >                                 dcpl_mb(paw_ovlp_w_d(1)+ii)      
c           
c           ii =ii +1
c        end do
c        end do
c        i_nbasis = i_nbasis + nbasis
c      end do
c      
c
c!      write(*,*) " paw weights"
c!      ii = 0
c!      nilm = 0
c!      njlm = 0
c!      do in=1,ion_nion()
c!        ia = ion_katm(in)
c!        i_w0 = int_mb(i_paw_ovlp_w(1) + in - 1)
c!        nbasis = paw_basis_nbasis(ia)
c!        nilm = 0
c!        do i=1,nbasis
c!          li = paw_basis_orb_l(i,ia)
c!          do mi=-li,li
c!            njlm = 0
c!            do j=1,nbasis
c!              lj = paw_basis_orb_l(j,ia)
c!                do mj=-lj,lj
c!                  i_w   = i_w0-1+(njlm+lj+mj+1)+
c!     >                   (nilm+li+mi)*paw_proj_nbasis(ia)
c!                  write(*,*) i,mi,j,mj,in
c!                  write(*,*)"w=",dcpl_mb(paw_ovlp_w(1)+ii)
c!                  write(*,*)"w_up=",dcpl_mb(paw_ovlp_w_u(1)+ii)
c!                  write(*,*)"w_down=",dcpl_mb(paw_ovlp_w_d(1)+ii)
c!                
c!                  ii =ii +1
c!                end do
c!                njlm=njlm + 2*lj+1
c!             end do
c!           end do   
c!           nilm = nilm + 2*li+1
c!         end do
c!      end do
c    
c      return
c      end
c

