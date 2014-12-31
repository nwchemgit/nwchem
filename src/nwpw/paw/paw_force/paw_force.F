*
* $Id$
*

!**************************************************
!
!       Name: paw_force_init
!
!       Purpose: 
!
!       Created:        8/07/2002
!**************************************************
      subroutine paw_force_init
cedo      subroutine paw_force_init(ispin0,ne)
      implicit none
cedo      integer ispin0
cedo      integer ne(2)

#include "bafdecls.fh"
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
      value =           my_alloc(mt_dcpl,paw_force_c_size,
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

      if (.not.value) call errquit('paw_force_init: alloc heap',0,1)

      return
      end

**************************************************
!
!       Name: paw_force_solve
!
!       Purpose: sets the paw overlap coefficients
!
!       Created:        8/07/2002
!**************************************************
      subroutine paw_force_solve(psi,fion)
      implicit none
      complex*16 psi(*)
      real*8     fion(3,*)
      
#include "bafdecls.fh"
#include "paw_force_data.fh"
#include "paw_proj.fh"
#include "paw_ovlp.fh"
#include "paw_geom.fh"
#include "paw_basis.fh"

      !*** local variables ***
      logical value
      integer taskid,np,tcount
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

      call Parallel_np(np)
      call Parallel_taskid(taskid)

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
      value =           BA_push_get(mt_dcpl,npack1,'exi',exi(2),exi(1))
      value = value.and.BA_push_get(mt_dcpl,npack1,'tmp0',
     >                     tmp0(2),tmp0(1))
      value = value.and.BA_push_get(mt_dcpl,npack1,'tmp',tmp(2),tmp(1))
      value = value.and.BA_push_get(mt_dcpl,npack1,'tmpx',
     >                   tmpx(2),tmpx(1))
      value = value.and.BA_push_get(mt_dcpl,npack1,'tmpy',
     >                   tmpy(2),tmpy(1))
      value = value.and.BA_push_get(mt_dcpl,npack1,'tmpz',
     >                   tmpz(2),tmpz(1))

      value = value.and. BA_push_get(mt_dbl,nfft3d,'Gx',Gx(2),Gx(1))
      value = value.and. BA_push_get(mt_dbl,nfft3d,'Gy',Gy(2),Gy(1))
      value = value.and. BA_push_get(mt_dbl,nfft3d,'Gz',Gz(2),Gz(1))
      if (.not.value) call errquit('paw_ovlp_coeff_set: push stack',0,0)

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


      value = value.and.BA_pop_stack(Gz(2))
      value = value.and.BA_pop_stack(Gy(2))
      value = value.and.BA_pop_stack(Gx(2))
      value = value.and.BA_pop_stack(tmpz(2))
      value = value.and.BA_pop_stack(tmpy(2))
      value = value.and.BA_pop_stack(tmpx(2))
      value = value.and.BA_pop_stack(tmp(2))
      value = value.and.BA_pop_stack(tmp0(2))
      value = value.and.BA_pop_stack(exi(2))
      if (.not.value) call errquit('paw_ovlp_coeff_set: pop stack',0,1)

      
c     *** calculate weights ***      
      if (paw_ispin().eq.1) then
        occup = 2.0d0
      else
        occup = 1.0d0
      end if

      if (np.gt.1) then
        tcount=2*paw_force_wx_u(3)
        call dcopy(tcount,0.0d0,0,dcpl_mb(paw_force_wx_u(1)),1)
        call dcopy(tcount,0.0d0,0,dcpl_mb(paw_force_wy_u(1)),1)
        call dcopy(tcount,0.0d0,0,dcpl_mb(paw_force_wz_u(1)),1)
        call dcopy(tcount,0.0d0,0,dcpl_mb(paw_force_wx_d(1)),1)
        call dcopy(tcount,0.0d0,0,dcpl_mb(paw_force_wy_d(1)),1)
        call dcopy(tcount,0.0d0,0,dcpl_mb(paw_force_wz_d(1)),1)
        call dcopy(tcount,0.0d0,0,dcpl_mb(paw_force_wx(1)),1)
        call dcopy(tcount,0.0d0,0,dcpl_mb(paw_force_wy(1)),1)
        call dcopy(tcount,0.0d0,0,dcpl_mb(paw_force_wz(1)),1)
      end if

      tcount = 0
      ii = 0
      i_nbasis = 0
      do in=1,ion_nion()

        ia = ion_katm(in)
        nbasis = paw_proj_nbasis(ia)
        do nlm =1,nbasis
        do nlm1=1,nbasis
        if (tcount.eq.taskid) then

           wx_tmp = dcmplx(0.0d0,0.0d0)
           wy_tmp = dcmplx(0.0d0,0.0d0)
           wz_tmp = dcmplx(0.0d0,0.0d0)
           do n=1,ne_u

             indx_coeff  = (i_nbasis+nlm-1)*neall + (n-1)
             indx_coeff1 = (i_nbasis+nlm1-1)*neall + (n-1)
             
              wx_tmp = wx_tmp 
     >               + dconjg(dcpl_mb(coeff(1)+indx_coeff))
     >                *dcpl_mb(paw_force_cx_u(1)+indx_coeff1)
     >               + dconjg(dcpl_mb(paw_force_cx_u(1)+indx_coeff))
     >                *dcpl_mb(coeff(1)+indx_coeff1)
             
              wy_tmp = wy_tmp 
     >               + dconjg(dcpl_mb(coeff(1)+indx_coeff))
     >                *dcpl_mb(paw_force_cy_u(1)+indx_coeff1)
     >               + dconjg(dcpl_mb(paw_force_cy_u(1)+indx_coeff))
     >                *dcpl_mb(coeff(1)+indx_coeff1)
             
              wz_tmp = wz_tmp 
     >               + dconjg(dcpl_mb(coeff(1)+indx_coeff))
     >                *dcpl_mb(paw_force_cz_u(1)+indx_coeff1)
     >               + dconjg(dcpl_mb(paw_force_cz_u(1)+indx_coeff))
     >                *dcpl_mb(coeff(1)+indx_coeff1)
             
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
     >              + dconjg(dcpl_mb(paw_force_cx_u(1)+indx_coeff))
     >               *dcpl_mb(coeff(1)+indx_coeff1)
             
              wy_tmp = wy_tmp 
     >               + dconjg(dcpl_mb(coeff(1)+indx_coeff))
     >                *dcpl_mb(paw_force_cy_u(1)+indx_coeff1)
     >               + dconjg(dcpl_mb(paw_force_cy_u(1)+indx_coeff))
     >                *dcpl_mb(coeff(1)+indx_coeff1)
             
              wz_tmp = wz_tmp 
     >               + dconjg(dcpl_mb(coeff(1)+indx_coeff))
     >                *dcpl_mb(paw_force_cz_u(1)+indx_coeff1)
     >               + dconjg(dcpl_mb(paw_force_cz_u(1)+indx_coeff))
     >                *dcpl_mb(coeff(1)+indx_coeff1)
 
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
           
        end if
        tcount=mod(tcount+1,np)
        ii=ii+1
        end do
        end do
        i_nbasis = i_nbasis + nbasis
      end do
      if (np.gt.1) then
        tcount=2*ii
        call D3dB_Vector_SumAll(tcount,dcpl_mb(paw_force_wx_u(1)))
        call D3dB_Vector_SumAll(tcount,dcpl_mb(paw_force_wy_u(1)))
        call D3dB_Vector_SumAll(tcount,dcpl_mb(paw_force_wz_u(1)))
        call D3dB_Vector_SumAll(tcount,dcpl_mb(paw_force_wx_d(1)))
        call D3dB_Vector_SumAll(tcount,dcpl_mb(paw_force_wy_d(1)))
        call D3dB_Vector_SumAll(tcount,dcpl_mb(paw_force_wz_d(1)))
        call D3dB_Vector_SumAll(tcount,dcpl_mb(paw_force_wx(1)))
        call D3dB_Vector_SumAll(tcount,dcpl_mb(paw_force_wy(1)))
        call D3dB_Vector_SumAll(tcount,dcpl_mb(paw_force_wz(1)))
      end if 

      call paw_force_atom(fion)
      call paw_force_mult(fion)


      return
      end

**************************************************
!
!       Name: paw_ovlp_weigths_set
!
!       Purpose: sets the overlap weight coefficients
!
!       Created:        8/07/2002
!**************************************************
      subroutine paw_force_weights_write(nu)
      
      implicit none
      integer nu

#include "bafdecls.fh"
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
!       Name: paw_force_atom
!
!       Purpose: 
!
!       Created:        8/07/2002
!**************************************************
      subroutine paw_force_atom(fion)
      implicit none
      real*8 fion(3,*)

#include "bafdecls.fh"
#include "paw_force_data.fh"
#include "paw_geom.fh"
#include "paw_proj.fh"


      logical value
      integer taskid,np,tcount,ft(2)
      integer i
      integer gop_u(4),gop_d(4)
      complex*16 tmp_fx,tmp_fy,tmp_fz

      integer in,ii,ia
      integer i_w,nb2
      integer nion
      integer i_paw_ovlp_w(4)
      
      call Parallel_np(np)
      call Parallel_taskid(taskid)

      nion = ion_nion()
      value = BA_push_get(mt_dbl,(3*nion),'ft',ft(2),ft(1))
      if (.not.value) call errquit('paw_force_atom: push stack',0,0)
      call dcopy(3*nion,0.0d0,0,dbl_mb(ft(1)),1)


      call paw_ovlp_i_w_ptr(i_paw_ovlp_w)

      call paw_nonloc_gop_u_ptr(gop_u)
      call paw_nonloc_gop_d_ptr(gop_d)


      tcount = 0
      do in=1,nion
        ia = ion_katm(in)
        nb2=paw_proj_nbasis(ia)**2

        i_w = int_mb(i_paw_ovlp_w(1) + in - 1)
 

        tmp_fx = dcmplx(0.0d0,0.0d0)
        tmp_fy = dcmplx(0.0d0,0.0d0)
        tmp_fz = dcmplx(0.0d0,0.0d0)
        do ii=1,nb2
        if (tcount.eq.taskid) then
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
        end if
        tcount = mod(tcount+1,np)
        end do

        !dbl_mb(paw_force_fx_atom(1)+in-1)=-dble(tmp_fx)
        !dbl_mb(paw_force_fy_atom(1)+in-1)=-dble(tmp_fy)
        !dbl_mb(paw_force_fz_atom(1)+in-1)=-dble(tmp_fz)

c        fion(1,in) = fion(1,in) - dble(tmp_fx)
c        fion(2,in) = fion(2,in) - dble(tmp_fy)
c        fion(3,in) = fion(3,in) - dble(tmp_fz)
        dbl_mb(ft(1)+3*(in-1))   = dbl_mb(ft(1)+3*(in-1))  -dble(tmp_fx)
        dbl_mb(ft(1)+3*(in-1)+1) = dbl_mb(ft(1)+3*(in-1)+1)-dble(tmp_fy)
        dbl_mb(ft(1)+3*(in-1)+2) = dbl_mb(ft(1)+3*(in-1)+2)-dble(tmp_fz)


      end do

      if (np.gt.1) call D3dB_Vector_SumAll(3*nion,dbl_mb(ft(1)))
      call daxpy(3*nion,1.0d0,dbl_mb(ft(1)),1,fion,1)

      value = BA_pop_stack(ft(2))
      if (.not.value) call errquit('paw_force_atom: pop stack',0,1)

      return
      end

!**************************************************
!
!       Name: paw_force_constraint
!
!       Purpose: 
!
!       Created:        8/07/2002
!**************************************************
      subroutine paw_force_constraint(lmbda,fion)
      implicit none
      real*8 lmbda(*)
      real*8 fion(3,*)

#include "bafdecls.fh"
#include "paw_force_data.fh"
#include "paw_proj.fh"
#include "paw_ovlp.fh"
#include "paw_geom.fh"
#include "paw_basis.fh"

      !**** local variables ****
      logical value
      integer taskid,np,tcount,ft(2)
      integer ii,ia,ms,n,m,a,b,ispin,n1(2),n2(2),nn,nion
      integer nbasis,basis_nbasis,i_nbasis
      integer la,ma,lb,mb,nlm1,nlm2,mlmbda,nlmbda
      integer matrix_ptr,overlap_ptr,ne_u,ne_d,neall
      integer indx_coeff1,indx_coeff,indx_lmbda
      integer coeff(4)
      real*8  occup,wx_tmp,wy_tmp,wz_tmp

      !**** local variables ****
      integer  paw_overlap_matrix_ptr
      external paw_overlap_matrix_ptr

      call Parallel_np(np)
      call Parallel_taskid(taskid)


      nion     = ion_nion()
      value = BA_push_get(mt_dbl,(3*nion),'ft',ft(2),ft(1))
      if (.not.value) 
     > call errquit('paw_force_constraint: push stack',0,0)
      call dcopy(3*nion,0.0d0,0,dbl_mb(ft(1)),1)


      call paw_ovlp_c_ptr(coeff)


c     *** number of electrons ***
      ne_u  = paw_ne_up()
      ne_d  = paw_ne_down()
      neall = paw_ne_total()

      ispin = paw_ispin()
      n1(1) = 1
      n2(1) = ne_u
      n1(2) = ne_u + 1
      n2(2) = ne_u + ne_d
      nn    = ne_u*ne_u
      occup = 1.0d0
      if (ispin.eq.1) occup = 2.0d0

      tcount = 0
      i_nbasis = 0
      do ii=1,nion
        ia           = ion_katm(ii)
        nbasis       = paw_proj_nbasis(ia)
        basis_nbasis = paw_basis_nbasis(ia)
        matrix_ptr   = paw_overlap_matrix_ptr(ia)
        
        nlm1 = 0
        do a=1,basis_nbasis
         la=paw_basis_orb_l(a,ia)
         do ma=-la,la
          nlm1 = nlm1+1

          nlm2 = 0
          do b=1,basis_nbasis
           lb=paw_basis_orb_l(b,ia)
           do mb=-lb,lb
            nlm2 = nlm2+1
         
            if ((mb.eq.ma).and.(la.eq.lb)) then
            if (tcount.eq.taskid) then

            overlap_ptr = matrix_ptr + (b-1)*basis_nbasis + (a-1)
            wx_tmp = 0.0d0
            wy_tmp = 0.0d0
            wz_tmp = 0.0d0
            do ms=1,ispin
            do n=n1(ms),n2(ms)
             indx_coeff  = (i_nbasis+nlm1- 1)*neall + (n-1)
             nlmbda      = n-n1(ms)
             do m=n1(ms),n2(ms)
              indx_coeff1 = (i_nbasis+nlm2-1)*neall + (m-1)
              mlmbda      = m-n1(ms)
              indx_lmbda = 1 + (mlmbda) + (nlmbda)*ne_u +(ms-1)*nn
  
              wx_tmp = wx_tmp
     >         + (  dble( dcpl_mb(coeff(1)         +indx_coeff))
     >             *dble( dcpl_mb(paw_force_cx_u(1)+indx_coeff1))
     >            + dimag(dcpl_mb(coeff(1)         +indx_coeff))
     >             *dimag(dcpl_mb(paw_force_cx_u(1)+indx_coeff1))
     >            + dble( dcpl_mb(paw_force_cx_u(1)+indx_coeff))
     >             *dble( dcpl_mb(coeff(1)         +indx_coeff1))
     >            + dimag(dcpl_mb(paw_force_cx_u(1)+indx_coeff))
     >             *dimag(dcpl_mb(coeff(1)         +indx_coeff1))
     >           )*lmbda(indx_lmbda)

             wy_tmp = wy_tmp
     >         + (  dble( dcpl_mb(coeff(1)         +indx_coeff))
     >             *dble( dcpl_mb(paw_force_cy_u(1)+indx_coeff1))
     >            + dimag(dcpl_mb(coeff(1)         +indx_coeff))
     >             *dimag(dcpl_mb(paw_force_cy_u(1)+indx_coeff1))
     >            + dble( dcpl_mb(paw_force_cy_u(1)+indx_coeff))
     >             *dble( dcpl_mb(coeff(1)         +indx_coeff1))
     >            + dimag(dcpl_mb(paw_force_cy_u(1)+indx_coeff))
     >             *dimag(dcpl_mb(coeff(1)         +indx_coeff1))
     >           )*lmbda(indx_lmbda)

             wz_tmp = wz_tmp
     >         + (  dble( dcpl_mb(coeff(1)         +indx_coeff))
     >             *dble( dcpl_mb(paw_force_cz_u(1)+indx_coeff1))
     >            + dimag(dcpl_mb(coeff(1)         +indx_coeff))
     >             *dimag(dcpl_mb(paw_force_cz_u(1)+indx_coeff1))
     >            + dble( dcpl_mb(paw_force_cz_u(1)+indx_coeff))
     >             *dble( dcpl_mb(coeff(1)         +indx_coeff1))
     >            + dimag(dcpl_mb(paw_force_cz_u(1)+indx_coeff))
     >             *dimag(dcpl_mb(coeff(1)         +indx_coeff1))
     >           )*lmbda(indx_lmbda)


             end do !*m*
            end do !*n*
            end do !*ms*
 
c            fion(1,ii) = fion(1,ii) + occup*wx_tmp*dbl_mb(overlap_ptr)
c            fion(2,ii) = fion(2,ii) + occup*wy_tmp*dbl_mb(overlap_ptr)
c            fion(3,ii) = fion(3,ii) + occup*wz_tmp*dbl_mb(overlap_ptr)

            dbl_mb(ft(1)+3*(ii-1))   = dbl_mb(ft(1)+3*(ii-1)) 
     >                               + occup*wx_tmp*dbl_mb(overlap_ptr)
            dbl_mb(ft(1)+3*(ii-1)+1) = dbl_mb(ft(1)+3*(ii-1)+1) 
     >                               + occup*wy_tmp*dbl_mb(overlap_ptr)
            dbl_mb(ft(1)+3*(ii-1)+2) = dbl_mb(ft(1)+3*(ii-1)+2) 
     >                               + occup*wz_tmp*dbl_mb(overlap_ptr)

            end if
            tcount = mod(tcount+1,np)
            end if !*ma.eq.mb.and.la.eq.lb*

           end do !*mb*
          end do !*lb*

         end do !*ma*
        end do !*la*

        i_nbasis = i_nbasis + nbasis
      end do !*ii*

      if (np.gt.1) call D3dB_Vector_SumAll(3*nion,dbl_mb(ft(1)))
      call daxpy(3*nion,1.0d0,dbl_mb(ft(1)),1,fion,1)

      value = BA_pop_stack(ft(2))
      if (.not.value) call errquit('paw_force_contraint: pop stack',0,1)

      return
      end



!**************************************************
!
!       Name: paw_force_mult
!
!       Purpose: 
!
!       Created:        8/07/2002
!**************************************************
      subroutine paw_force_mult(f_mult)

      implicit none

      double precision f_mult(3,*)

      


#include "bafdecls.fh"
#include "paw_force_data.fh"
#include "paw_mult_data.fh"
#include "paw_geom.fh"
#include "paw_proj.fh"
#include "paw_basis.fh"
#include "paw_comp_charge_data.fh"      

      !*** local variables ***
      logical value
      integer ft(2),taskid,np,tcount
      integer nion
      integer il,ilp
      integer i_lshell,lshell
      integer in,ia,jn,ja
      integer l1,l2,m1,m2
      integer mult_l,mult_l_prime
      integer i_qlm,i_qlmp
      integer k
      
      double precision sigma_a, sigma_a_prime
      double precision ra(3), ra_prime(3)
      double precision R
      
      complex*16 comp_charge_coeff
      complex*16 f_tmp(3)

      call Parallel_np(np)
      call Parallel_taskid(taskid)

      nion = ion_nion()
      value = BA_push_get(mt_dbl,(3*nion),'ft',ft(2),ft(1))
      if (.not.value) call errquit('paw_force_mult: push stack',0,0)
      call dcopy(3*nion,0.0d0,0,dbl_mb(ft(1)),1)


      !do in=1,nion
      !do k=1,3
      !  f_mult(k,in) = 0.0d0
      !end do
      !end do

      tcount=0
      do in=1,nion
      ra(1)   = ion_rion(1,in)
      ra(2)   = ion_rion(2,in)
      ra(3)   = ion_rion(3,in)
      ia = ion_katm(in)
      sigma_a = paw_basis_sigma(ia)
      mult_l  = paw_basis_mult_l(ia)
      i_qlm = int_mb(i_paw_qlm(1) + in - 1)

      do jn=1,nion
        i_qlmp = int_mb(i_paw_qlm(1) + jn - 1)
        ja = ion_katm(jn)
        sigma_a_prime = paw_basis_sigma(ja)
        mult_l_prime  = paw_basis_mult_l(ja)

        i_lshell = 0
        do lshell=1,nshl3d
         ra_prime(1)   = ion_rion(1,jn) + dbl_mb(rcell(1)+i_lshell)
         i_lshell = i_lshell+1
         ra_prime(2)   = ion_rion(2,jn) + dbl_mb(rcell(1)+i_lshell)
         i_lshell = i_lshell+1
         ra_prime(3)   = ion_rion(3,jn) + dbl_mb(rcell(1)+i_lshell)
         i_lshell = i_lshell+1

        
         R = dsqrt( (ra(1)-ra_prime(1))**2
     >              +(ra(2)-ra_prime(2))**2
     >              +(ra(3)-ra_prime(3))**2)

         if (R.lt.(4*sigma_smooth)) then
         if (tcount.eq.taskid) then

            il = 1
            do l1=0,mult_l
            do m1=-l1,l1

              ilp = 1
              do l2=0,mult_l_prime
              do m2=-l2,l2
                      
               comp_charge_coeff = dcpl_mb(paw_qlm(1)+i_qlm+il-1)*
     >                             dcpl_mb(paw_qlm(1)+i_qlmp+ilp-1)

               call paw_mult_force(l1,m1,sigma_a,
     >                             ra,
     >                             l2,m2,sigma_a_prime,
     >                             ra_prime,
     >                             f_tmp(1),
     >                             f_tmp(2),
     >                             f_tmp(3))

                do k=1,3
c                 f_mult(k,in) = f_mult(k,in)-
c     >                          dble(comp_charge_coeff*f_tmp(k))
                 dbl_mb(ft(1)+k-1+(in-1)*3) 
     >           = dbl_mb(ft(1)+k-1+(in-1)*3) 
     >           - dble(comp_charge_coeff*f_tmp(k))
                end do

               call paw_mult_force(l1,m1,sigma_smooth,
     >                             ra,
     >                             l2,m2,sigma_smooth,
     >                             ra_prime,
     >                             f_tmp(1),
     >                             f_tmp(2),
     >                             f_tmp(3))

                do k=1,3
c                 f_mult(k,in) = f_mult(k,in)-
c     >                          dble(comp_charge_coeff*f_tmp(k))
                 dbl_mb(ft(1)+k-1+(in-1)*3) 
     >           = dbl_mb(ft(1)+k-1+(in-1)*3) 
     >           - dble(comp_charge_coeff*f_tmp(k))
                end do


               call paw_mult_force(l1,m1,sigma_smooth,
     >                             ra,
     >                             l2,m2,sigma_a_prime,
     >                             ra_prime,
     >                             f_tmp(1),
     >                             f_tmp(2),
     >                             f_tmp(3))

                do k=1,3
c                 f_mult(k,in) = f_mult(k,in)+
c     >                          dble(comp_charge_coeff*f_tmp(k))
                 dbl_mb(ft(1)+k-1+(in-1)*3) 
     >           = dbl_mb(ft(1)+k-1+(in-1)*3) 
     >           + dble(comp_charge_coeff*f_tmp(k))
                end do


               call paw_mult_force(l1,m1,sigma_a,
     >                             ra,
     >                             l2,m2,sigma_smooth,
     >                             ra_prime,
     >                             f_tmp(1),
     >                             f_tmp(2),
     >                             f_tmp(3))

                do k=1,3
c                 f_mult(k,in) = f_mult(k,in)+
c     >                          dble(comp_charge_coeff*f_tmp(k))
                 dbl_mb(ft(1)+k-1+(in-1)*3) 
     >           = dbl_mb(ft(1)+k-1+(in-1)*3) 
     >           + dble(comp_charge_coeff*f_tmp(k))
                end do




                ilp = ilp + 1
              end do !*m2*
              end do !*l2*

              il = il + 1
            end do !*m1*
            end do !*l1*

          end if
          tcount = mod(tcount+1,np)
          end if

        end do !*lshell*
        end do !*jn*
        end do !*in*

      if (np.gt.1) call D3dB_Vector_SumAll(3*nion,dbl_mb(ft(1)))
      call daxpy(3*nion,1.0d0,dbl_mb(ft(1)),1,f_mult,1)

      value = BA_pop_stack(ft(2))
      if (.not.value) call errquit('paw_force_mult: pop stack',0,1)

      return
      end

!**************************************************
!
!       Name: paw_ovlp_end
!
!       Purpose: removes space used by the paw overlap coefficients
!
!       Created:        8/07/2002
!**************************************************
      subroutine paw_force_end()
      implicit none      

#include "bafdecls.fh"
#include "paw_ma.fh" 
#include "paw_force_data.fh"

      !*** local variables ***
      logical value

      value = .true.
      value = value.and.my_dealloc(paw_force_cx_u)
      value = value.and.my_dealloc(paw_force_cx_s)
      value = value.and.my_dealloc(paw_force_cy_u)
      value = value.and.my_dealloc(paw_force_cy_s)
      value = value.and.my_dealloc(paw_force_cz_u)
      value = value.and.my_dealloc(paw_force_cz_s)

      value = value.and.my_dealloc(paw_force_wx)
      value = value.and.my_dealloc(paw_force_wx_u)
      value = value.and.my_dealloc(paw_force_wx_d)
      value = value.and.my_dealloc(paw_force_wy)
      value = value.and.my_dealloc(paw_force_wy_u)
      value = value.and.my_dealloc(paw_force_wy_d)
      value = value.and.my_dealloc(paw_force_wz)
      value = value.and.my_dealloc(paw_force_wz_u)
      value = value.and.my_dealloc(paw_force_wz_d)

c
      if (.not.value) call errquit('paw_force_end: dealloc heap',0,0)
      return
      end


 


