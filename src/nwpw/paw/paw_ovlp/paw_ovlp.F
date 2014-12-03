*
* $Id$
*

!**************************************************
!
!       Name: paw_ovlp_init
!
!       Purpose: initializes the paw overlap coefficients
!
!       Created:        8/07/2002
!**************************************************
      subroutine paw_ovlp_init(ispin0,ne)
      implicit none
      integer ispin0
      integer ne(2)

#include "bafdecls.fh"
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
      if (.not.value) call errquit('paw_ovlp_init: alloc heap',0,1)

      paw_ovlp_w_size = 0
      do in=1,nion
         int_mb(i_paw_ovlp_w(1) + in - 1) = paw_ovlp_w_size
         ia = ion_katm(in)
         paw_ovlp_w_size = paw_ovlp_w_size
     >                + paw_proj_nbasis(ia)**2
      end do



c     !*** allocate storage for paw_ovlp arrays*** 
c     !-- total number of electron states
      ispin = ispin0
      neall = ne(1) + ne(2)
      ne_u=ne(1)
      ne_d=ne(2)

      value =           my_alloc(mt_dcpl,(total_nbasis*neall),
     >                               'coeff',coeff)
      value = value.and.my_alloc(mt_dcpl,
     >                               (total_nbasis*neall),
     >                               'toeff',toeff)
      value = value.and.my_alloc(mt_int,(total_nbasis*neall),
     >                               'indx2',indx2)
      value = value.and.my_alloc(mt_int,(total_nbasis*neall),
     >                               'sgn',sgn)

      value = value.and.my_alloc(mt_dcpl,paw_ovlp_w_size,
     >                           "paw_ovlp_w",paw_ovlp_w)

      value = value.and.my_alloc(mt_dcpl,paw_ovlp_w_size,
     >                           "paw_ovlp_w_u",paw_ovlp_w_spin(1,1))

      if (ispin.eq.2)
     > value = value.and.my_alloc(mt_dcpl,paw_ovlp_w_size,
     >                           "paw_ovlp_w_d",paw_ovlp_w_spin(1,2))

      if (.not.value) call errquit('paw_ovlp_init: alloc heap',0,1)


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

              int_mb(indx2(1)+i1) = i2     ! location of -m elements
              int_mb(sgn(1)  +i1) = isgn   ! sign array=(-1)**m 
            end do
         end do 
      end do    


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
      subroutine paw_ovlp_end()
      implicit none      

#include "bafdecls.fh"
#include "paw_ma.fh" 
#include "paw_ovlp_data.fh"

      !*** local variables ***
      logical value
      
      value =           BA_free_heap(coeff(2))
      value = value.and.BA_free_heap(toeff(2))
      value = value.and.BA_free_heap(indx2(2))
      value = value.and.BA_free_heap(sgn(2))
      value = value.and.BA_free_heap(paw_ovlp_w(2))
      value = value.and.BA_free_heap(paw_ovlp_w_spin(2,1))
      if (ispin.eq.2)
     > value = value.and.BA_free_heap(paw_ovlp_w_spin(2,2))
      value = value.and.my_dealloc(i_paw_ovlp_w)

      if (.not.value) call errquit('paw_ovlp_end: dealloc heap',0,0)
      return
      end


**************************************************
!
!       Name: paw_ovlp_coeff_set
!
!       Purpose: sets the paw overlap coefficients
!
!       Created:        8/07/2002
!**************************************************
      subroutine paw_ovlp_coeff_set(psi)
      implicit none
      complex*16 psi(*)
      
#include "bafdecls.fh"
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


      call nwpw_timing_start(6)

      !*** allocate temporary space ***
      call Pack_npack(1,npack1)
      value =           BA_push_get(mt_dcpl,npack1,'exi',exi(2),exi(1))
      value = value.and.BA_push_get(mt_dcpl,npack1,'tmp',tmp(2),tmp(1))
      if (.not.value) call errquit('paw_ovlp_coeff_set: push stack',0,0)



c      !*** compute overlap coefficients c(n=1:neall,nlm=1:nbasis(ii),ii=1:nion) ***
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





      value =           BA_pop_stack(tmp(2))
      value = value.and.BA_pop_stack(exi(2))
      if (.not.value) call errquit('paw_ovlp_coeff_set: pop stack',0,1)
      
      call nwpw_timing_end(6)

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
      subroutine paw_ovlp_weights_set()
      
      implicit none

#include "bafdecls.fh"
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
      integer n1(2),n2(2),ms
      real*8  occup
      double complex w_tmp
      
      
      call nwpw_timing_start(6)
      !*** check for restricted calculation ****
c      occup = 1.0d0
      n1(1) = 1
      n2(1) = ne_u
      n1(2) = ne_u+1
      n2(2) = ne_u+ne_d

      ii = 0
      i_nbasis = 0
      do in=1,ion_nion()
        ia = ion_katm(in)
        nbasis = paw_proj_nbasis(ia)
        do nlm =1,nbasis
        do nlm1=1,nbasis

           do ms=1,ispin
             w_tmp = dcmplx(0.0d0,0.0d0)
             do n=n1(ms),n2(ms)
               indx_coeff  = (i_nbasis+nlm-1)*neall + (n-1)
               indx_coeff1 = (i_nbasis+nlm1-1)*neall + (n-1)
               w_tmp = w_tmp + conjg(dcpl_mb(coeff(1)+indx_coeff))
     >                         *dcpl_mb(coeff(1)+indx_coeff1)
             end do
             dcpl_mb(paw_ovlp_w_spin(1,ms)+ii) = w_tmp
           end do

           dcpl_mb(paw_ovlp_w(1)+ii) 
     >            = dcpl_mb(paw_ovlp_w_spin(1,1)+ii)
     >            + dcpl_mb(paw_ovlp_w_spin(1,ispin)+ii)      
           
           ii =ii +1
        end do
        end do
        i_nbasis = i_nbasis + nbasis
      end do
      

!      write(*,*) " paw weights"
!      ii = 0
!      nilm = 0
!      njlm = 0
!      do in=1,ion_nion()
!        ia = ion_katm(in)
!        i_w0 = int_mb(i_paw_ovlp_w(1) + in - 1)
!        nbasis = paw_basis_nbasis(ia)
!        nilm = 0
!        do i=1,nbasis
!          li = paw_basis_orb_l(i,ia)
!          do mi=-li,li
!            njlm = 0
!            do j=1,nbasis
!              lj = paw_basis_orb_l(j,ia)
!                do mj=-lj,lj
!                  i_w   = i_w0-1+(njlm+lj+mj+1)+
!     >                   (nilm+li+mi)*paw_proj_nbasis(ia)
!                  write(*,*) i,mi,j,mj,in
!                  write(*,*)"w=",dcpl_mb(paw_ovlp_w(1)+ii)
!                  write(*,*)"w_up=",dcpl_mb(paw_ovlp_w_u(1)+ii)
!                  write(*,*)"w_down=",dcpl_mb(paw_ovlp_w_d(1)+ii)
!                
!                  ii =ii +1
!                end do
!                njlm=njlm + 2*lj+1
!             end do
!           end do   
!           nilm = nilm + 2*li+1
!         end do
!      end do
      call nwpw_timing_end(6)
    
      return
      end

**************************************************
!
!       Name: paw_ovlp_calc_coef
!
!       Purpose: Calculate the projectors coefficients
!
!       Created:        2/10/2002
!**************************************************
      subroutine paw_ovlp_calc_coef(ne,psi,c1)
      implicit none
      integer ne
      complex*16 psi(*)
      complex*16 c1(ne,*)

#include "bafdecls.fh"
#include "paw_ovlp_data.fh"

      !**** local variables ****
      logical ok
      integer ct(2)
      !*** allocate space from stack ****
      ok = BA_push_get(mt_dcpl,ne*total_nbasis,'ct',ct(2),ct(1))
      if (.not.ok)
     > call errquit('paw_ovlp_calc_coeffn: out of stack memory',0,0)


      call paw_ovlp_calc_scramble_coef(ne,psi,dcpl_mb(ct(1)))
      call paw_ovlp_unscramble_coef(ne,dcpl_mb(ct(1)),c1)


      !*** deallocate temporary space ***
      ok = BA_pop_stack(ct(2))
      if (.not.ok)
     > call errquit(
     >  'paw_ovlp_calc_coeff: error freeing stack memory',1,0)

      return
      end


**************************************************
!
!       Name: paw_ovlp_calc_scramble_coef
!
!       Purpose: Calculate the scrambled coefficients
!
!       Created:        2/10/2002
!**************************************************
      subroutine paw_ovlp_calc_scramble_coef(ne,psi,ct)
      implicit none
      integer ne
      complex*16 psi(*)
      complex*16 ct(ne,*)

#include "bafdecls.fh"
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


      call nwpw_timing_start(6)
      !*** allocate temporary space ***
      call Pack_npack(1,npack1)
      value =           BA_push_get(mt_dcpl,npack1,'exi',exi(2),exi(1))
      value = value.and.BA_push_get(mt_dcpl,npack1,'tmp',tmp(2),tmp(1))
      if (.not.value) 
     > call errquit('paw_ovlp_coeff_set: out of stack memory',0,0)


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
      value =           BA_pop_stack(tmp(2))
      value = value.and.BA_pop_stack(exi(2))
      if (.not.value)  call errquit(
     > 'paw_ovlp_coeff_set: error freeing stack memory',0,1)

      call nwpw_timing_end(6)
      return
      end



**************************************************
!
!       Name: paw_ovlp_unscramble_coef
!
!       Purpose: unscrambe the coefficients
!
!       Created:        2/10/2002
!**************************************************
      subroutine paw_ovlp_unscramble_coef(ne,ct,c1)
      implicit none
      integer ne
      complex*16 ct(ne,*),c1(ne,*)


#include "bafdecls.fh"
#include "paw_ma.fh"
#include "paw_proj.fh"
#include "paw_ovlp_data.fh"
#include "paw_geom.fh"

*     **** local variables ****
      integer ii,ia,nlm,nlmii,n,m,mabs
      integer nbasis,nion,i2
      real*8 isgn

      call nwpw_timing_start(6)
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

      call nwpw_timing_end(6)
      return
      end

**************************************************
!
!       Name: paw_ovlp_S
!
!       Purpose: Apply the overlap operator
!
!       Created:        2/10/2002
!**************************************************
      subroutine paw_ovlp_S(ne,psi,Spsi)
      implicit none
      integer ne
      complex*16 psi(*)
      complex*16 Spsi(*)

#include "bafdecls.fh"
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
      ok =        BA_push_get(mt_dcpl,npack1,'exi',exi(2),exi(1))
      ok = ok.and.BA_push_get(mt_dcpl,npack1,'tmp',tmp(2),tmp(1))
      ok = ok.and.
     >     BA_push_get(mt_dcpl,ne*total_nbasis,'c1',c1(2),c1(1))
      ok = ok.and.
     >     BA_push_get(mt_dcpl,ne*total_nbasis,'pv',pv(2),pv(1))
      if (.not.ok)
     > call errquit('paw_ovlp_S: out of stack memory',0,0)


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
      ok =        BA_pop_stack(pv(2))
      ok = ok.and.BA_pop_stack(c1(2))
      ok = ok.and.BA_pop_stack(tmp(2))
      ok = ok.and.BA_pop_stack(exi(2))
      if (.not.ok) 
     > call errquit('paw_ovlp_S: error freeing stack memory',0,1)

      return
      end

**************************************************
!
!       Name: paw_ovlp_coeff_ptr
!
!       Purpose: 
!
!       Created:        2/10/2002
!**************************************************
      subroutine paw_ovlp_coeff_ptr(ptr)
      implicit none
      integer ptr

#include "paw_ovlp_data.fh"

      ptr = coeff(1)
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
      subroutine paw_ovlp_indx2_ptr(ptr)
      implicit none
      integer ptr(4)

#include "paw_ovlp_data.fh"
      integer i

      do i=1,4
      ptr(i) = indx2(i)
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
      subroutine paw_ovlp_sgn_ptr(ptr)
      implicit none
      integer ptr(4)

#include "paw_ovlp_data.fh"

      integer i

      do i=1,4
      ptr(i) = sgn(i)
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
      subroutine paw_ovlp_c_ptr(ptr)
      implicit none
      integer ptr(4)

#include "paw_ovlp_data.fh"

      integer i

      do i=1,4
      ptr(i) = coeff(i)
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
      subroutine paw_ovlp_i_w_ptr(ptr)
      implicit none
      integer ptr(4)

#include "paw_ovlp_data.fh"

      integer i

      do i=1,4
      ptr(i) = i_paw_ovlp_w(i)
      end do
      return
      end 



**************************************************
!
!       Name: paw_ovlp_coeff_size
!
!       Purpose:  returns the size of the overlap array
!
!       Created:        2/10/2002
!**************************************************
      function paw_ovlp_coeff_size()
      implicit none
       
      integer paw_ovlp_coeff_size 

#include "paw_ovlp_data.fh"

      paw_ovlp_coeff_size=coeff(3)
      return
      end 
**************************************************
!
!       Name: paw_ovlp_coeff_ptr
!
!       Purpose: 
!
!       Created:        2/10/2002
!**************************************************
      function paw_ovlp_weight_size()
      implicit none
       
      integer paw_ovlp_weight_size 

#include "paw_ovlp_data.fh"

      paw_ovlp_weight_size=paw_ovlp_w(3)
      return
      end 
**************************************************
!
!       Name: 
!
!       Purpose: 
!
!       Created:
!**************************************************
      function paw_ne_total()
      implicit none
       
      integer paw_ne_total

#include "paw_ovlp_data.fh"

      paw_ne_total=neall
      return
      end 
**************************************************
!
!       Name: 
!
!       Purpose: 
!
!       Created:
!**************************************************
      function paw_ne_up()
      implicit none
       
      integer paw_ne_up

#include "paw_ovlp_data.fh"

      paw_ne_up=ne_u
      return
      end 
**************************************************
!
!       Name: 
!
!       Purpose: 
!
!       Created:
!**************************************************
      function paw_ne_down()
      implicit none
       
      integer paw_ne_down

#include "paw_ovlp_data.fh"

      paw_ne_down=ne_d
      return
      end 
**************************************************
!
!       Name: 
!
!       Purpose: 
!
!       Created:
!**************************************************
      function paw_ispin()
      implicit none
       
      integer paw_ispin

#include "paw_ovlp_data.fh"

      paw_ispin = ispin
      return
      end 

**************************************************
!
!       Name: paw_ovlp_weigths_write
!
!       Purpose: dumps the overlap weight coefficients
!
!       Created:        8/07/2002
!**************************************************
      subroutine paw_ovlp_weights_write(nu)
      
      implicit none
      integer nu

#include "bafdecls.fh"
#include "paw_geom.fh"
#include "paw_ovlp_data.fh"
#include "paw_proj.fh"
#include "paw_basis.fh"

      integer ii
      integer in,ia
      integer nbasis
      integer i,j,mi,mj,li,lj
      integer i_w0,i_w
      integer nilm,njlm
      

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

                  write(nu,*) i,mi,j,mj,in
                  write(nu,*)dcpl_mb(paw_ovlp_w(1)+ii)
                
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
c
