*
* $Id: paw_ovlp.F,v 1.33 2003-02-16 22:13:30 marat Exp $
*

!**************************************************
!
!	Name: paw_ovlp_init
!
!	Purpose: initializes the paw overlap coefficients
!
!	Created:	8/07/2002
!**************************************************
      subroutine paw_ovlp_init(ispin,ne)
      implicit none
      integer ispin
      integer ne(2)

#include "mafdecls.fh"
#include "paw_proj.fh"
#include "paw_ovlp_data.fh" 
#include "paw_geom.fh" 
#include "paw_ma.fh" 
#include "paw_basis.fh" 
  
*     **** local variables ****
      logical value
      integer ii,ia,nlm,nlmii,n,m,mabs
      integer nbasis,i2,i1,isgn,nion
      integer paw_ovlp_w_size
      integer in,li,lj,mi,mj,i,j


      nion = ion_nion()

c     !*** calculate total number of (n,l,m) projectors  ***
      total_nbasis = 0
      do ia=1,ion_nkatm()
c        !-- ion_natm(ia) is number of atoms of kind ia --
         total_nbasis = total_nbasis 
     >                + paw_proj_nbasis(ia)*ion_natm(ia)
      end do

c     !*** calculate the size and indexing for weight coefficients  ***
      value = my_alloc(mt_int,ion_nion(),"i_paw_ovlp_w",i_paw_ovlp_w)
      if (.not.value) call errquit('paw_ovlp_init: alloc heap',1)

      paw_ovlp_w_size = 0
      do in=1,nion
         int_mb(i_paw_ovlp_w(1) + in - 1) = paw_ovlp_w_size
         ia = ion_katm(in)
         paw_ovlp_w_size = paw_ovlp_w_size
     >                + paw_proj_nbasis(ia)**2
      end do



c     !*** allocate storage for paw_ovlp arrays*** 
c     !-- total number of electron states
      neall = ne(1) + ne(2)
      ne_u=ne(1)
      ne_d=ne(2)

      value =           MA_alloc_get(mt_dcpl,(total_nbasis*neall),
     >                               'coeff',coeff(2),coeff(1))
      value = value.and.MA_alloc_get(mt_dcpl,
     >                               (total_nbasis*neall),
     >                               'toeff',toeff(2),toeff(1))
      value = value.and.MA_alloc_get(mt_int,(total_nbasis*neall),
     >                               'indx2',indx2(2),indx2(1))
      value = value.and.MA_alloc_get(mt_int,(total_nbasis*neall),
     >                               'sgn',sgn(2),sgn(1))

      value = value.and.my_alloc(mt_dcpl,paw_ovlp_w_size,
     >                           "paw_ovlp_w",paw_ovlp_w)

      value = value.and.my_alloc(mt_dcpl,paw_ovlp_w_size,
     >                           "paw_ovlp_w",paw_ovlp_w_u)

      value = value.and.my_alloc(mt_dcpl,paw_ovlp_w_size,
     >                           "paw_ovlp_w",paw_ovlp_w_d)

      if (.not.value) call errquit('paw_ovlp_init: alloc heap',1)


c    !    atom=1                        atom=2
c    !l : 0  1  1  1  2  2  2  2  2 ... 0  1  1  1 ... 
c    !m : 0 -1  0  1 -2 -1  0  1  1 ... 0 -1  0  1 ...
c    !
c    !

c     !*** set sign array (isgn) and location of -m elements(indx2)
      nlmii = 0
      do ii=1,nion
c        !-- the atom kind is ia --
         ia = ion_katm(ii)
         nbasis = paw_proj_nbasis(ia)
c        !-- loop over projectors for atom of kind ia --
         do nlm=1,nbasis
            nlmii = nlmii + 1
            m = paw_proj_m(nlm,ia)
            mabs = abs(m)
c           !-- set sign coeff as (-1)**m --
            if (mod(mabs,2).eq.0) then
              isgn = 1
            else
              isgn = -1
            end if
            do n=1,neall
              i1 = (nlmii-1)*neall + (n-1)
              i2 = i1 - 2*m*neall
c              write(*,*) "i1,i2,total_nbasis,m:",i1,i2,total_nbasis,m
              int_mb(indx2(1)+i1) = i2     ! location of -m elements
              int_mb(sgn(1)  +i1) = isgn   ! sign array=(-1)**m 
            end do
         end do 
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
      subroutine paw_ovlp_end()
      implicit none      

#include "mafdecls.fh"
#include "paw_ma.fh" 
#include "paw_ovlp_data.fh"

      !*** local variables ***
      logical value
      
      value =           MA_free_heap(coeff(2))
      value = value.and.MA_free_heap(toeff(2))
      value = value.and.MA_free_heap(indx2(2))
      value = value.and.MA_free_heap(sgn(2))
      value = value.and.MA_free_heap(paw_ovlp_w(2))
      value = value.and.my_dealloc(i_paw_ovlp_w)
      if (.not.value) call errquit('paw_ovlp_end: dealloc heap',0)
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
      subroutine paw_ovlp_coeff_set(psi)
      implicit none
      complex*16 psi(*)
      
#include "mafdecls.fh"
#include "paw_ovlp_data.fh"
#include "paw_proj.fh"
#include "paw_geom.fh"
#include "paw_basis.fh"

      !*** local variables ***
      logical value
      integer ii,ia,nbasis,nlm,nlmii,npack1
      integer nion,prj,i1,i2,n
      integer indx_psi,indx_prj,indx_coeff
      integer exi(2),tmp(2)
      integer ni,i,li,mi





      !*** allocate temporary space ***
      call Pack_npack(1,npack1)
      value =           MA_push_get(mt_dcpl,npack1,'exi',exi(2),exi(1))
      value = value.and.MA_push_get(mt_dcpl,npack1,'tmp',tmp(2),tmp(1))
      if (.not.value) call errquit('paw_ovlp_coeff_set: push stack',0)



c      !*** compute overlap coefficients c(n=1:neall,nlm=1:nbasis(ii),ii=1:nion) ***
      nlmii = 0
      nion  = ion_nion()
      do ii=1,nion
        ia = ion_katm(ii)
        call strfac_pack(1,ii,dcpl_mb(exi(1)))
        nbasis = paw_proj_nbasis(ia)
        prj = paw_proj_i_prj_atom(ia)
        do nlm=1,nbasis
           write(*,*) "lm:",paw_proj_l(nlm,ia),paw_proj_m(nlm,ia)
           nlmii = nlmii + 1
           indx_prj = (nlm-1)*npack1
           call Pack_cc_Mul(1,
     >                      dcpl_mb(exi(1)),
     >                      dcpl_mb(prj+indx_prj),
     >                      dcpl_mb(tmp(1)))
           do n=1,neall
             indx_psi   = (n-1)*npack1    + 1
             indx_coeff = (nlmii-1)*neall + (n-1)
             call Pack_cc_izdot(1,
     >                         dcpl_mb(tmp(1)),
     >                         psi(indx_psi),
     >                         dcpl_mb(toeff(1)+indx_coeff))
           end do
        end do
      end do 
      call D3dB_Vector_SumAll(2*neall*total_nbasis,
     >                        dcpl_mb(toeff(1)))


      !*** unsrcramble to get true overlap coefficients ***
      do i1=1,total_nbasis*neall
         i2 = int_mb(indx2(1)+i1-1)
         dcpl_mb(coeff(1)+i1-1) 
     >     = 0.5d0*( dcpl_mb(toeff(1)+i1-1)
     >             + int_mb(sgn(1)+i1-1)
     >              *dconjg(dcpl_mb(toeff(1)+i2)))
      end do     

      write(*,*)
      write(*,*) "Overlap coefficients:"
      nlmii = 0
      do ni=1,nion
        ia = ion_katm(ni)
        nbasis = paw_basis_nbasis(ia)
        do i=1,nbasis
          li = paw_basis_orb_l(i,ia)
          do mi = -li,li
           nlmii = nlmii + 1
           do n=1,neall
             indx_coeff = (nlmii-1)*neall + (n-1)
             write(*,*) i,mi,ni,n
             write(*,*) dcpl_mb(coeff(1)+indx_coeff)
           end do
           end do
        end do
      end do 


c      !*** debug - write out the overlap coefficients ***
c      write(*,*)
c      write(*,*) "Overlap coefficients:",neall,total_nbasis
c      do i2=1,neall
c      do i1=1,total_nbasis
c         indx_coeff = (i2-1) + (i1-1)*neall
c         write(*,*) "nlm n c(n,nlm):",i1,i2,dcpl_mb(coeff(1)+indx_coeff)
c      end do
c      end do
c      write(*,*)

      value =           MA_pop_stack(tmp(2))
      value = value.and.MA_pop_stack(exi(2))
      if (.not.value) call errquit('paw_ovlp_coeff_set: pop stack',1)
      

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
      subroutine paw_ovlp_weights_set()
      
      implicit none

#include "mafdecls.fh"
#include "paw_geom.fh"
#include "paw_ovlp_data.fh"
#include "paw_proj.fh"
#include "paw_basis.fh"

      integer i_nbasis,ii,ii1
      integer in,nlm,nlm1,ia,n
      integer nbasis
      integer indx_coeff,indx_coeff1
      integer i,j,mi,mj,li,lj
      integer i_w0,i_w
      integer nilm,njlm,niplm,njplm
      double complex w_tmp
      
      

      ii = 0
      i_nbasis = 0
      do in=1,ion_nion()
        ia = ion_katm(in)
        nbasis = paw_proj_nbasis(ia)
        do nlm =1,nbasis
        do nlm1=1,nbasis

           w_tmp = cmplx(0.0d0,0.0d0)
           do n=1,ne_u
             indx_coeff  = (i_nbasis+nlm-1)*neall + (n-1)
             indx_coeff1 = (i_nbasis+nlm1-1)*neall + (n-1)
             w_tmp = w_tmp + conjg(dcpl_mb(coeff(1)+indx_coeff))*
     >        dcpl_mb(coeff(1)+indx_coeff1)
           end do
           dcpl_mb(paw_ovlp_w_u(1)+ii) = w_tmp

           w_tmp = cmplx(0.0d0,0.0d0)
           do n=ne_u+1,ne_d
             indx_coeff  = (i_nbasis+nlm-1)*neall + (n-1)
             indx_coeff1 = (i_nbasis+nlm1-1)*neall + (n-1)
             w_tmp = w_tmp + conjg(dcpl_mb(coeff(1)+indx_coeff))*
     >        dcpl_mb(coeff(1)+indx_coeff1)
           end do
           dcpl_mb(paw_ovlp_w_d(1)+ii) = w_tmp

           dcpl_mb(paw_ovlp_w(1)+ii) = dcpl_mb(paw_ovlp_w_u(1)+ii)+
     >                                 dcpl_mb(paw_ovlp_w_d(1)+ii)      
           
           ii =ii +1
        end do
        end do
        i_nbasis = i_nbasis + nbasis
      end do
      

      write(*,*) " paw weights"
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
                  i_w   = i_w0-1+(njlm+lj+mj+1)+
     >                   (nilm+li+mi)*paw_proj_nbasis(ia)
                  write(*,*) i,mi,j,mj,in
                  write(*,*)ii,i_w,i_w0,dcpl_mb(paw_ovlp_w(1)+ii),
     >                      dcpl_mb(paw_ovlp_w(1)+i_w)
                
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

**************************************************
!
!	Name: paw_ovlp_calc_coef
!
!	Purpose: Calculate the projectors coefficients
!
!	Created:	2/10/2002
!**************************************************
      subroutine paw_ovlp_calc_coef(ne,psi,c1)
      implicit none
      integer ne
      complex*16 psi(*)
      complex*16 c1(ne,*)

#include "mafdecls.fh"
#include "paw_ovlp_data.fh"

      !**** local variables ****
      logical ok
      integer ct(2)
      !*** allocate space from stack ****
      ok = MA_push_get(mt_dcpl,ne*total_nbasis,'ct',ct(2),ct(1))
      if (.not.ok)
     > call errquit('paw_ovlp_calc_coeffn: out of stack memory',0)


      call paw_ovlp_calc_scramble_coef(ne,psi,dcpl_mb(ct(1)))
      call paw_ovlp_unscramble_coef(ne,dcpl_mb(ct(1)),c1)


      !*** deallocate temporary space ***
      ok = MA_pop_stack(ct(2))
      if (.not.ok)
     > call errquit(
     >  'paw_ovlp_calc_coeff: error freeing stack memory',1)

      return
      end


**************************************************
!
!	Name: paw_ovlp_calc_scramble_coef
!
!	Purpose: Calculate the scrambled coefficients
!
!	Created:	2/10/2002
!**************************************************
      subroutine paw_ovlp_calc_scramble_coef(ne,psi,ct)
      implicit none
      integer ne
      complex*16 psi(*)
      complex*16 ct(ne,*)

#include "mafdecls.fh"
#include "paw_ovlp_data.fh"
#include "paw_proj.fh"
#include "paw_geom.fh"
#include "paw_basis.fh"

      !*** local variables ***
      logical value
      integer ii,ia,nbasis,nlm,nlmii,npack1
      integer nion,prj,n
      integer indx_psi,indx_prj
      integer exi(2),tmp(2)


      !*** allocate temporary space ***
      call Pack_npack(1,npack1)
      value =           MA_push_get(mt_dcpl,npack1,'exi',exi(2),exi(1))
      value = value.and.MA_push_get(mt_dcpl,npack1,'tmp',tmp(2),tmp(1))
      if (.not.value) 
     > call errquit('paw_ovlp_coeff_set: out of stack memory',0)


c     !*** compute overlap coefficients c(n=1:ne,nlm=1:nbasis(ii),ii=1:nion) ***
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
           do n=1,ne
             indx_psi   = (n-1)*npack1    + 1
             call Pack_cc_izdot(1,
     >                         dcpl_mb(tmp(1)),
     >                         psi(indx_psi),
     >                         ct(n,nlmii))
           end do
        end do
      end do
      call D3dB_Vector_SumAll(2*ne*total_nbasis,ct)


      !*** deallocate temporary space ***
      value =           MA_pop_stack(tmp(2))
      value = value.and.MA_pop_stack(exi(2))
      if (.not.value) 
     > call errquit('paw_ovlp_coeff_set: error freeing stack memory',1)

      return
      end



**************************************************
!
!	Name: paw_ovlp_unscramble_coef
!
!	Purpose: unscrambe the coefficients
!
!	Created:	2/10/2002
!**************************************************
      subroutine paw_ovlp_unscramble_coef(ne,ct,c1)
      implicit none
      integer ne
      complex*16 ct(ne,*),c1(ne,*)


#include "mafdecls.fh"
#include "paw_ma.fh"
#include "paw_proj.fh"
#include "paw_ovlp_data.fh"
#include "paw_geom.fh"

*     **** local variables ****
      integer ii,ia,nlm,nlmii,n,m,mabs
      integer nbasis,nion,i2
      real*8 isgn

      nion = ion_nion()
      nlmii = 0
      do ii=1,nion
         ia = ion_katm(ii)
         nbasis = paw_proj_nbasis(ia)

         do nlm=1,nbasis
            nlmii = nlmii + 1
            m = paw_proj_m(nlm,ia)
            mabs = abs(m)
            if (mod(mabs,2).eq.0) then
              isgn = 1.0d0
            else
              isgn = -1.0d0
            end if
            i2 = nlmii-2*m
            do n=1,ne
              c1(n,nlmii) = 0.5d0*(ct(n,nlmii) 
     >                   + isgn*dconjg(ct(n,i2)))
            end do
         end do
      end do

      return
      end

**************************************************
!
!	Name: paw_ovlp_S
!
!	Purpose: Apply the overlap operator
!
!	Created:	2/10/2002
!**************************************************
      subroutine paw_ovlp_S(ne,psi,Spsi)
      implicit none
      integer ne
      complex*16 psi(*)
      complex*16 Spsi(*)

#include "mafdecls.fh"
#include "paw_ovlp_data.fh"
#include "paw_proj.fh"
#include "paw_geom.fh"
#include "paw_basis.fh"

      !*** local variables ***
      logical ok
      integer ii,ia,nbasis,nlm,nlmii
      integer nion,prj,n,npack1
      integer indx_psi,indx_prj,indx
      integer c1(2),pv(2),exi(2),tmp(2)


      !*** allocate temporary space ***
      call Pack_npack(1,npack1)
      ok =        MA_push_get(mt_dcpl,npack1,'exi',exi(2),exi(1))
      ok = ok.and.MA_push_get(mt_dcpl,npack1,'tmp',tmp(2),tmp(1))
      ok = ok.and.
     >     MA_push_get(mt_dcpl,ne*total_nbasis,'c1',c1(2),c1(1))
      ok = ok.and.
     >     MA_push_get(mt_dcpl,ne*total_nbasis,'pv',pv(2),pv(1))
      if (.not.ok)
     > call errquit('paw_ovlp_S: out of stack memory',0)


      !**** calculate pv ****
      call paw_ovlp_calc_coef(ne,psi,dcpl_mb(c1(1)))
      call  paw_overlap_pv_gen(ne,total_nbasis,
     >          dcpl_mb(c1(1)),
     >          dcpl_mb(pv(1)))


      call dcopy(2*ne*npack1,psi,1,Spsi,1)
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

           do n=1,ne
             indx_psi = (n-1)*npack1    + 1
             indx     = (nlmii-1)*ne + (n-1)
             call Pack_cc_zaxpy(1,
     >                         dcpl_mb(pv(1)+indx),
     >                         dcpl_mb(tmp(1)),
     >                         Spsi(indx_psi))
           end do
        end do
      end do

      !*** deallocate temporary space ***
      ok =        MA_pop_stack(pv(2))
      ok = ok.and.MA_pop_stack(c1(2))
      ok = ok.and.MA_pop_stack(tmp(2))
      ok = ok.and.MA_pop_stack(exi(2))
      if (.not.ok) 
     > call errquit('paw_ovlp_S: error freeing stack memory',1)

      return
      end

