*
* $Id: paw_force.F,v 1.5 2003-03-15 00:05:20 bylaska Exp $
*

!**************************************************
!
!	Name: paw_force_init
!
!	Purpose: 
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
      integer nion
      integer paw_force_w_size
      integer paw_force_c_size


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

      value = value.and.my_alloc(mt_dbl,nion,
     >                           "paw_force_fx_atom",paw_force_fx_atom)
      value = value.and.my_alloc(mt_dbl,nion,
     >                           "paw_force_fy_atom",paw_force_fy_atom)
      value = value.and.my_alloc(mt_dbl,nion,
     >                           "paw_force_fz_atom",paw_force_fz_atom)


      if (.not.value) call errquit('paw_force_init: alloc heap',1)

      return
      end

**************************************************
!
!	Name: paw_force_solve
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
      integer exi(2),tmp(2),tmp0(2)
      integer tmpx(2),tmpy(2),tmpz(2)
      integer Gx(2),Gy(2),Gz(2)
      integer nfft3d
      integer indx2(4),sgn(4),coeff(4)
      integer ne_u,ne_d,neall
      double precision occup
      complex*16 wx_tmp,wy_tmp,wz_tmp
      integer i_nbasis,in,nlm1,indx_coeff1

      integer ni,i,li,mi,ik
 
      !**** external functions ****
      integer  G_indx
      external G_indx


c     *** number of electrons ***      
      ne_u  = paw_ne_up()
      ne_d  = paw_ne_down()
      neall = paw_ne_total()

c     *** pointers to paw_ovlp arrays ***     
      call paw_ovlp_indx2_ptr(indx2)
      call paw_ovlp_sgn_ptr(sgn)
      call paw_ovlp_c_ptr(coeff)

c     *** size of unpacked FFT grid ***      
      call D3dB_nfft3d(1,nfft3d)
c     *** size of packed grid ***      
      call Pack_npack(1,npack1)

c     *** allocate temporary space ***
      value =           MA_push_get(mt_dcpl,npack1,'exi',exi(2),exi(1))
      value = value.and.MA_push_get(mt_dcpl,npack1,'tmp0',
     >                     tmp0(2),tmp0(1))
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

c    *** get reciprocal lattice arrays ***
      call D3dB_t_Copy(1,dbl_mb(G_indx(1)),dbl_mb(Gx(1)))
      call D3dB_t_Copy(1,dbl_mb(G_indx(2)),dbl_mb(Gy(1)))
      call D3dB_t_Copy(1,dbl_mb(G_indx(3)),dbl_mb(Gz(1)))
      call Pack_t_pack(1,dbl_mb(Gx(1)))
      call Pack_t_pack(1,dbl_mb(Gy(1)))
      call Pack_t_pack(1,dbl_mb(Gz(1)))


c     *** compute force overlap coefficients ***
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
     >                      dcpl_mb(tmp0(1)))

           do ik=1,npack1
              dcpl_mb(tmpx(1)+ik-1)=dcpl_mb(tmp0(1)+ik-1)*
     >           dcmplx(0.0,-1.0)*dbl_mb(Gx(1)+ik-1)

              dcpl_mb(tmpy(1)+ik-1)=dcpl_mb(tmp0(1)+ik-1)*
     >           dcmplx(0.0,-1.0)*dbl_mb(Gy(1)+ik-1)

              dcpl_mb(tmpz(1)+ik-1)=dcpl_mb(tmp0(1)+ik-1)*
     >           dcmplx(0.0,-1.0)*dbl_mb(Gz(1)+ik-1)
           
           end do

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
     >                        dcpl_mb(paw_force_cx_s(1)))
      call D3dB_Vector_SumAll(2*paw_force_cy_s(3),
     >                        dcpl_mb(paw_force_cy_s(1)))
      call D3dB_Vector_SumAll(2*paw_force_cz_s(3),
     >                        dcpl_mb(paw_force_cz_s(1)))


      !*** unsrcramble to get true force overlap coefficients ***
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

c      write(*,*)
c      write(*,*) "paw force coefficients1:"
c      nlmii = 0
c      do ni=1,nion
c        ia = ion_katm(ni)
c        nbasis = paw_basis_nbasis(ia)
c        do i=1,nbasis
c          li = paw_basis_orb_l(i,ia)
c          do mi = -li,li
c           nlmii = nlmii + 1
c           do n=1,neall
c             indx_coeff = (nlmii-1)*neall + (n-1)
c             write(*,*) i,mi,ni,n
c             write(*,*) dcpl_mb(paw_force_cx_u(1)+indx_coeff),
c     >                  dcpl_mb(paw_force_cy_u(1)+indx_coeff),        
c     >                  dcpl_mb(paw_force_cz_u(1)+indx_coeff)        
c           end do
c           end do
c        end do
c      end do 




      value = value.and.MA_pop_stack(Gz(2))
      value = value.and.MA_pop_stack(Gy(2))
      value = value.and.MA_pop_stack(Gx(2))
      value = value.and.MA_pop_stack(tmpz(2))
      value = value.and.MA_pop_stack(tmpy(2))
      value = value.and.MA_pop_stack(tmpx(2))
      value = value.and.MA_pop_stack(tmp(2))
      value = value.and.MA_pop_stack(tmp0(2))
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
             
              wx_tmp = wx_tmp 
     >               + dconjg(dcpl_mb(coeff(1)+indx_coeff))
     >                *dcpl_mb(paw_force_cx_u(1)+indx_coeff1)
     >               + dconjg(dcpl_mb(coeff(1)+indx_coeff))
     >                *dcpl_mb(paw_force_cx_u(1)+indx_coeff1)
             
              wy_tmp = wy_tmp 
     >               + dconjg(dcpl_mb(coeff(1)+indx_coeff))
     >                *dcpl_mb(paw_force_cy_u(1)+indx_coeff1)
     >               + dconjg(dcpl_mb(coeff(1)+indx_coeff))
     >                *dcpl_mb(paw_force_cy_u(1)+indx_coeff1)
             
              wz_tmp = wz_tmp 
     >               + dconjg(dcpl_mb(coeff(1)+indx_coeff))
     >                *dcpl_mb(paw_force_cz_u(1)+indx_coeff1)
     >               + dconjg(dcpl_mb(coeff(1)+indx_coeff))
     >                *dcpl_mb(paw_force_cz_u(1)+indx_coeff1)
             
           end do
           dcpl_mb(paw_force_wx_u(1)+ii) = wx_tmp*occup
           dcpl_mb(paw_force_wy_u(1)+ii) = wy_tmp*occup
           dcpl_mb(paw_force_wz_u(1)+ii) = wz_tmp*occup


c           !*** this loop is not computed when ispin=1 ***
           wx_tmp = dcmplx(0.0d0,0.0d0)
           wy_tmp = dcmplx(0.0d0,0.0d0)
           wz_tmp = dcmplx(0.0d0,0.0d0)
           do n=ne_u+1,ne_d+ne_u
             indx_coeff  = (i_nbasis+nlm-1)*neall + (n-1)
             indx_coeff1 = (i_nbasis+nlm1-1)*neall + (n-1)

             wx_tmp = wx_tmp 
     >              + dconjg(dcpl_mb(coeff(1)+indx_coeff))
     >               *dcpl_mb(paw_force_cx_u(1)+indx_coeff1)
     >              + dconjg(dcpl_mb(coeff(1)+indx_coeff))
     >               *dcpl_mb(paw_force_cx_u(1)+indx_coeff1)
             
              wy_tmp = wy_tmp 
     >               + dconjg(dcpl_mb(coeff(1)+indx_coeff))
     >                *dcpl_mb(paw_force_cy_u(1)+indx_coeff1)
     >               + dconjg(dcpl_mb(coeff(1)+indx_coeff))
     >                *dcpl_mb(paw_force_cy_u(1)+indx_coeff1)
             
              wz_tmp = wz_tmp 
     >               + dconjg(dcpl_mb(coeff(1)+indx_coeff))
     >                *dcpl_mb(paw_force_cz_u(1)+indx_coeff1)
     >               + dconjg(dcpl_mb(coeff(1)+indx_coeff))
     >                *dcpl_mb(paw_force_cz_u(1)+indx_coeff1)
 
           end do
           dcpl_mb(paw_force_wx_d(1)+ii) = wx_tmp*occup
           dcpl_mb(paw_force_wy_d(1)+ii) = wy_tmp*occup
           dcpl_mb(paw_force_wz_d(1)+ii) = wz_tmp*occup

           dcpl_mb(paw_force_wx(1)+ii) = dcpl_mb(paw_force_wx_u(1)+ii)
     >                                 + dcpl_mb(paw_force_wx_d(1)+ii)    
           dcpl_mb(paw_force_wy(1)+ii) = dcpl_mb(paw_force_wy_u(1)+ii)
     >                                 + dcpl_mb(paw_force_wy_d(1)+ii)    
           dcpl_mb(paw_force_wz(1)+ii) = dcpl_mb(paw_force_wz_u(1)+ii)
     >                                 + dcpl_mb(paw_force_wz_d(1)+ii)    
           
           ii=ii +1
        end do
        end do
        i_nbasis = i_nbasis + nbasis
      end do

      call paw_force_atom()
      return
      end

**************************************************
!
!	Name: paw_ovlp_weigths_set
!
!	Purpose: sets the overlap weight coefficients
!
!	Created:	8/07/2002
!**************************************************
      subroutine paw_force_weights_write(nu)
      
      implicit none
      integer nu

#include "mafdecls.fh"
#include "paw_geom.fh"
#include "paw_force_data.fh"
#include "paw_proj.fh"
#include "paw_basis.fh"

      integer ii
      integer in,ia
      integer nbasis
      integer i,j,mi,mj,li,lj
      integer i_w0,i_w
      integer nilm,njlm
      integer i_paw_ovlp_w(4)
      
      call paw_ovlp_i_w_ptr(i_paw_ovlp_w)

      ii = 0
      nilm = 0
      njlm = 0
      do in=1,ion_nion()
        ia = ion_katm(in)
        i_w0 = int_mb(i_paw_ovlp_w(1) + in - 1)
        nbasis = paw_basis_nbasis(ia)
        nilm = 0
        do i=1,nbasis
          li = paw_basis_orb_l(i,ia)
          do mi=-li,li
            njlm = 0
            do j=1,nbasis
              lj = paw_basis_orb_l(j,ia)
                do mj=-lj,lj
                  i_w   = i_w0-1+(njlm+lj+mj+1)
     >                  +(nilm+li+mi)*paw_proj_nbasis(ia)

                  write(nu,*) i,mi,j,mj,in
                  write(nu,*) dcpl_mb(paw_force_wx_u(1)+ii)
                  write(nu,*) dcpl_mb(paw_force_wy_u(1)+ii)
                  write(nu,*) dcpl_mb(paw_force_wz_u(1)+ii)
                
                  ii =ii +1
                end do
                njlm=njlm + 2*lj+1
             end do
           end do   
           nilm = nilm + 2*li+1
         end do
      end do
    
      return
      end

!**************************************************
!
!	Name: paw_force_atom
!
!	Purpose: 
!
!	Created:	8/07/2002
!**************************************************
      subroutine paw_force_atom()
      implicit none

#include "mafdecls.fh"
#include "paw_force_data.fh"
#include "paw_geom.fh"
#include "paw_proj.fh"


      integer i
      integer gop_u(4),gop_d(4)
      complex*16 tmp_fx,tmp_fy,tmp_fz

      integer in,ii,ia
      integer i_w,nb2
      integer nion
      integer i_paw_ovlp_w(4)
      
      call paw_ovlp_i_w_ptr(i_paw_ovlp_w)

      call paw_nonloc_gop_u_ptr(gop_u)
      call paw_nonloc_gop_d_ptr(gop_d)

      write(*,*) " paw force, gop"
      do in=1,gop_u(3)
      write(*,*) 
     > dcpl_mb(gop_u(1)+in-1),
     > dcpl_mb(gop_d(1)+in-1),
     > dcpl_mb(paw_force_wx_u(1)+in-1)
      end do

      nion = ion_nion()
      do in=1,nion
        ia = ion_katm(in)
        nb2=paw_proj_nbasis(ia)**2

        i_w = int_mb(i_paw_ovlp_w(1) + in - 1)
 

        tmp_fx = dcmplx(0.0d0,0.0d0)
        tmp_fy = dcmplx(0.0d0,0.0d0)
        tmp_fz = dcmplx(0.0d0,0.0d0)
        do ii=1,nb2
           tmp_fx = tmp_fx + 
     >             dcpl_mb(gop_u(1)+i_w+ii-1)*
     >             dcpl_mb(paw_force_wx_u(1)+i_w+ii-1)     
     >             +
     >             dcpl_mb(gop_d(1)+i_w+ii-1)*
     >             dcpl_mb(paw_force_wx_d(1)+i_w+ii-1)     
           tmp_fy = tmp_fy + 
     >             dcpl_mb(gop_u(1)+i_w+ii-1)*
     >             dcpl_mb(paw_force_wy_u(1)+i_w+ii-1)     
     >             +
     >             dcpl_mb(gop_d(1)+i_w+ii-1)*
     >             dcpl_mb(paw_force_wy_d(1)+i_w+ii-1)     
           tmp_fz = tmp_fz + 
     >             dcpl_mb(gop_u(1)+i_w+ii-1)*
     >             dcpl_mb(paw_force_wz_u(1)+i_w+ii-1)     
     >             +
     >             dcpl_mb(gop_d(1)+i_w+ii-1)*
     >             dcpl_mb(paw_force_wz_d(1)+i_w+ii-1)     
        end do

        dbl_mb(paw_force_fx_atom(1)+in-1)=-dble(tmp_fx)
        dbl_mb(paw_force_fy_atom(1)+in-1)=-dble(tmp_fy)
        dbl_mb(paw_force_fz_atom(1)+in-1)=-dble(tmp_fz)

      end do


      write(*,*) "atomic force"
      do in=1,nion

        write(*,*) in, 
     >             dbl_mb(paw_force_fx_atom(1)+in-1),
     >             dbl_mb(paw_force_fy_atom(1)+in-1),
     >             dbl_mb(paw_force_fz_atom(1)+in-1)

      end do
      return
      end

!**************************************************
!
!	Name: paw_ovlp_end
!
!	Purpose: removes space used by the paw overlap coefficients
!
!	Created:	8/07/2002
!**************************************************
      subroutine paw_force_end()
      implicit none      

#include "mafdecls.fh"
#include "paw_ma.fh" 
#include "paw_ovlp_data.fh"

      !*** local variables ***
      logical value

      value = .true.
      
c      value =           MA_free_heap(coeff(2))
c      value = value.and.MA_free_heap(toeff(2))
c      value = value.and.MA_free_heap(indx2(2))
c      value = value.and.MA_free_heap(sgn(2))
c      value = value.and.MA_free_heap(paw_ovlp_w(2))
c      value = value.and.MA_free_heap(paw_ovlp_w_u(2))
c      value = value.and.MA_free_heap(paw_ovlp_w_d(2))
c      value = value.and.my_dealloc(i_paw_ovlp_w)
c
      if (.not.value) call errquit('paw_force_end: dealloc heap',0)
      return
      end


 


