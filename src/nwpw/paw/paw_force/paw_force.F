*
* $Id: paw_force.F,v 1.1 2003-03-10 05:45:24 marat Exp $
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

      ispin = ispin0
      neall = ne(1) + ne(2)
      ne_u=ne(1)
      ne_d=ne(2)

      nion = ion_nion()

c     !*** calculate size for force coeff arrays ***
      paw_force_c_size = paw_ovlp_coeff_size()
c      paw_force_c_size = 0
c      do ia=1,ion_nkatm()
c         paw_force_c_size = paw_force_c_size 
c     >                + paw_proj_nbasis(ia)*ion_natm(ia)
c      end do
c      paw_force_c_size = paw_force_c_size*neall

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


      value = value.and.my_alloc(mt_int,paw_force_c_size,
     >                               'indx2',indx2)
      value = value.and.my_alloc(mt_int,paw_force_c_size,
     >                               'sgn',sgn)


c     !*** calculate the size and indexing for weight coefficients  ***
      value = my_alloc(mt_int,ion_nion(),"i_paw_force_w",i_paw_force_w)
      if (.not.value) call errquit('paw_force_init: alloc heap',1)

      paw_force_w_size = paw_ovlp_weight_size()
c      paw_force_w_size = 0
c      do in=1,nion
c         int_mb(i_paw_force_w(1) + in - 1) = paw_force_w_size
c         ia = ion_katm(in)
c         paw_force_w_size = paw_force_w_size
c     >                + paw_proj_nbasis(ia)**2
c      end do


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
c!	Name: paw_ovlp_coeff_set
c!
c!	Purpose: sets the paw overlap coefficients
c!
c!	Created:	8/07/2002
c!**************************************************
c      subroutine paw_ovlp_coeff_set(psi)
c      implicit none
c      complex*16 psi(*)
c      
c#include "mafdecls.fh"
c#include "paw_ovlp_data.fh"
c#include "paw_proj.fh"
c#include "paw_geom.fh"
c#include "paw_basis.fh"
c
c      !*** local variables ***
c      logical value
c      integer ii,ia,nbasis,nlm,nlmii,npack1
c      integer nion,prj,i1,i2,n
c      integer indx_psi,indx_prj,indx_coeff
c      integer exi(2),tmp(2)
c      integer ni,i,li,mi
c
c
c
c      !*** allocate temporary space ***
c      call Pack_npack(1,npack1)
c      value =           MA_push_get(mt_dcpl,npack1,'exi',exi(2),exi(1))
c      value = value.and.MA_push_get(mt_dcpl,npack1,'tmp',tmp(2),tmp(1))
c      if (.not.value) call errquit('paw_ovlp_coeff_set: push stack',0)
c
c
c
cc      !*** compute overlap coefficients c(n=1:neall,nlm=1:nbasis(ii),ii=1:nion) ***
c      nlmii = 0
c      nion  = ion_nion()
c      do ii=1,nion
c        ia = ion_katm(ii)
c        call strfac_pack(1,ii,dcpl_mb(exi(1)))
c        nbasis = paw_proj_nbasis(ia)
c        prj = paw_proj_i_prj_atom(ia)
c        do nlm=1,nbasis
c           
c           nlmii = nlmii + 1
c           indx_prj = (nlm-1)*npack1
c           call Pack_cc_Mul(1,
c     >                      dcpl_mb(exi(1)),
c     >                      dcpl_mb(prj+indx_prj),
c     >                      dcpl_mb(tmp(1)))
c           do n=1,neall
c             indx_psi   = (n-1)*npack1    + 1
c             indx_coeff = (nlmii-1)*neall + (n-1)
c             call Pack_cc_izdot(1,
c     >                         dcpl_mb(tmp(1)),
c     >                         psi(indx_psi),
c     >                         dcpl_mb(toeff(1)+indx_coeff))
c           end do
c        end do
c      end do 
c      call D3dB_Vector_SumAll(2*neall*total_nbasis,
c     >                        dcpl_mb(toeff(1)))
c
c
c      !*** unsrcramble to get true overlap coefficients ***
c      do i1=1,total_nbasis*neall
c         i2 = int_mb(indx2(1)+i1-1)
c         dcpl_mb(coeff(1)+i1-1) 
c     >     = 0.5d0*( dcpl_mb(toeff(1)+i1-1)
c     >             + int_mb(sgn(1)+i1-1)
c     >              *dconjg(dcpl_mb(toeff(1)+i2)))
c      end do     
c
c
c
c
c
c      value =           MA_pop_stack(tmp(2))
c      value = value.and.MA_pop_stack(exi(2))
c      if (.not.value) call errquit('paw_ovlp_coeff_set: pop stack',1)
c      
c
c      return
c      end
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
c**************************************************
c!
c!	Name: paw_ovlp_calc_coef
c!
c!	Purpose: Calculate the projectors coefficients
c!
c!	Created:	2/10/2002
c!**************************************************
c      subroutine paw_ovlp_calc_coef(ne,psi,c1)
c      implicit none
c      integer ne
c      complex*16 psi(*)
c      complex*16 c1(ne,*)
c
c#include "mafdecls.fh"
c#include "paw_ovlp_data.fh"
c
c      !**** local variables ****
c      logical ok
c      integer ct(2)
c      !*** allocate space from stack ****
c      ok = MA_push_get(mt_dcpl,ne*total_nbasis,'ct',ct(2),ct(1))
c      if (.not.ok)
c     > call errquit('paw_ovlp_calc_coeffn: out of stack memory',0)
c
c
c      call paw_ovlp_calc_scramble_coef(ne,psi,dcpl_mb(ct(1)))
c      call paw_ovlp_unscramble_coef(ne,dcpl_mb(ct(1)),c1)
c
c
c      !*** deallocate temporary space ***
c      ok = MA_pop_stack(ct(2))
c      if (.not.ok)
c     > call errquit(
c     >  'paw_ovlp_calc_coeff: error freeing stack memory',1)
c
c      return
c      end
c
c
c**************************************************
c!
c!	Name: paw_ovlp_calc_scramble_coef
c!
c!	Purpose: Calculate the scrambled coefficients
c!
c!	Created:	2/10/2002
c!**************************************************
c      subroutine paw_ovlp_calc_scramble_coef(ne,psi,ct)
c      implicit none
c      integer ne
c      complex*16 psi(*)
c      complex*16 ct(ne,*)
c
c#include "mafdecls.fh"
c#include "paw_ovlp_data.fh"
c#include "paw_proj.fh"
c#include "paw_geom.fh"
c#include "paw_basis.fh"
c
c      !*** local variables ***
c      logical value
c      integer ii,ia,nbasis,nlm,nlmii,npack1
c      integer nion,prj,n
c      integer indx_psi,indx_prj
c      integer exi(2),tmp(2)
c
c
c      !*** allocate temporary space ***
c      call Pack_npack(1,npack1)
c      value =           MA_push_get(mt_dcpl,npack1,'exi',exi(2),exi(1))
c      value = value.and.MA_push_get(mt_dcpl,npack1,'tmp',tmp(2),tmp(1))
c      if (.not.value) 
c     > call errquit('paw_ovlp_coeff_set: out of stack memory',0)
c
c
cc     !*** compute overlap coefficients c(n=1:ne,nlm=1:nbasis(ii),ii=1:nion) ***
c      nlmii = 0
c      nion  = ion_nion()
c      do ii=1,nion
c        ia = ion_katm(ii)
c        call strfac_pack(1,ii,dcpl_mb(exi(1)))
c        nbasis = paw_proj_nbasis(ia)
c        prj = paw_proj_i_prj_atom(ia)
c        do nlm=1,nbasis
c           nlmii = nlmii + 1
c           indx_prj = (nlm-1)*npack1
c           call Pack_cc_Mul(1,
c     >                      dcpl_mb(exi(1)),
c     >                      dcpl_mb(prj+indx_prj),
c     >                      dcpl_mb(tmp(1)))
c           do n=1,ne
c             indx_psi   = (n-1)*npack1    + 1
c             call Pack_cc_izdot(1,
c     >                         dcpl_mb(tmp(1)),
c     >                         psi(indx_psi),
c     >                         ct(n,nlmii))
c           end do
c        end do
c      end do
c      call D3dB_Vector_SumAll(2*ne*total_nbasis,ct)
c
c
c      !*** deallocate temporary space ***
c      value =           MA_pop_stack(tmp(2))
c      value = value.and.MA_pop_stack(exi(2))
c      if (.not.value) 
c     > call errquit('paw_ovlp_coeff_set: error freeing stack memory',1)
c
c      return
c      end
c
c
c
c**************************************************
c!
c!	Name: paw_ovlp_unscramble_coef
c!
c!	Purpose: unscrambe the coefficients
c!
c!	Created:	2/10/2002
c!**************************************************
c      subroutine paw_ovlp_unscramble_coef(ne,ct,c1)
c      implicit none
c      integer ne
c      complex*16 ct(ne,*),c1(ne,*)
c
c
c#include "mafdecls.fh"
c#include "paw_ma.fh"
c#include "paw_proj.fh"
c#include "paw_ovlp_data.fh"
c#include "paw_geom.fh"
c
c*     **** local variables ****
c      integer ii,ia,nlm,nlmii,n,m,mabs
c      integer nbasis,nion,i2
c      real*8 isgn
c
c      nion = ion_nion()
c      nlmii = 0
c      do ii=1,nion
c         ia = ion_katm(ii)
c         nbasis = paw_proj_nbasis(ia)
c
c         do nlm=1,nbasis
c            nlmii = nlmii + 1
c            m = paw_proj_m(nlm,ia)
c            mabs = abs(m)
c            if (mod(mabs,2).eq.0) then
c              isgn = 1.0d0
c            else
c              isgn = -1.0d0
c            end if
c            i2 = nlmii-2*m
c            do n=1,ne
c              c1(n,nlmii) = 0.5d0*(ct(n,nlmii) 
c     >                   + isgn*dconjg(ct(n,i2)))
c            end do
c         end do
c      end do
c
c      return
c      end
c
c**************************************************
c!
c!	Name: paw_ovlp_S
c!
c!	Purpose: Apply the overlap operator
c!
c!	Created:	2/10/2002
c!**************************************************
c      subroutine paw_ovlp_S(ne,psi,Spsi)
c      implicit none
c      integer ne
c      complex*16 psi(*)
c      complex*16 Spsi(*)
c
c#include "mafdecls.fh"
c#include "paw_ovlp_data.fh"
c#include "paw_proj.fh"
c#include "paw_geom.fh"
c#include "paw_basis.fh"
c
c      !*** local variables ***
c      logical ok
c      integer ii,ia,nbasis,nlm,nlmii
c      integer nion,prj,n,npack1
c      integer indx_psi,indx_prj,indx
c      integer c1(2),pv(2),exi(2),tmp(2)
c
c
c      !*** allocate temporary space ***
c      call Pack_npack(1,npack1)
c      ok =        MA_push_get(mt_dcpl,npack1,'exi',exi(2),exi(1))
c      ok = ok.and.MA_push_get(mt_dcpl,npack1,'tmp',tmp(2),tmp(1))
c      ok = ok.and.
c     >     MA_push_get(mt_dcpl,ne*total_nbasis,'c1',c1(2),c1(1))
c      ok = ok.and.
c     >     MA_push_get(mt_dcpl,ne*total_nbasis,'pv',pv(2),pv(1))
c      if (.not.ok)
c     > call errquit('paw_ovlp_S: out of stack memory',0)
c
c
c      !**** calculate pv ****
c      call paw_ovlp_calc_coef(ne,psi,dcpl_mb(c1(1)))
c      call  paw_overlap_pv_gen(ne,total_nbasis,
c     >          dcpl_mb(c1(1)),
c     >          dcpl_mb(pv(1)))
c
c
c      call dcopy(2*ne*npack1,psi,1,Spsi,1)
c      nlmii = 0
c      nion  = ion_nion()
c      do ii=1,nion
c        ia = ion_katm(ii)
c        call strfac_pack(1,ii,dcpl_mb(exi(1)))
c        nbasis = paw_proj_nbasis(ia)
c        prj = paw_proj_i_prj_atom(ia)
c        do nlm=1,nbasis
c           nlmii = nlmii + 1
c           indx_prj = (nlm-1)*npack1
c           call Pack_cc_Mul(1,
c     >                      dcpl_mb(exi(1)),
c     >                      dcpl_mb(prj+indx_prj),
c     >                      dcpl_mb(tmp(1)))
c
c           do n=1,ne
c             indx_psi = (n-1)*npack1    + 1
c             indx     = (nlmii-1)*ne + (n-1)
c             call Pack_cc_zaxpy(1,
c     >                         dcpl_mb(pv(1)+indx),
c     >                         dcpl_mb(tmp(1)),
c     >                         Spsi(indx_psi))
c           end do
c        end do
c      end do
c
c      !*** deallocate temporary space ***
c      ok =        MA_pop_stack(pv(2))
c      ok = ok.and.MA_pop_stack(c1(2))
c      ok = ok.and.MA_pop_stack(tmp(2))
c      ok = ok.and.MA_pop_stack(exi(2))
c      if (.not.ok) 
c     > call errquit('paw_ovlp_S: error freeing stack memory',1)
c
c      return
c      end
c
c**************************************************
c!
c!	Name: paw_ovlp_coeff_ptr
c!
c!	Purpose: 
c!
c!	Created:	2/10/2002
c!**************************************************
c      subroutine paw_ovlp_coeff_ptr(ptr)
c      implicit none
c      integer ptr
c
c#include "paw_ovlp_data.fh"
c
c      ptr = coeff(1)
c      return
c      end 
c
