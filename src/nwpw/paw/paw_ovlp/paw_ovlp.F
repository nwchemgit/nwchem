*
* $Id: paw_ovlp.F,v 1.18 2003-02-08 23:05:14 bylaska Exp $
*

!**************************************************
!
!	Name: paw_ovlp_init
!
!	Purpose: initializes the paw overlap coefficients
!
!	Created:	8/07/2002
!**************************************************
      subroutine paw_ovlp_init()
      implicit none

#include "mafdecls.fh"
#include "paw_proj.fh"
#include "paw_ovlp_data.fh" 
#include "paw_geom.fh" 
#include "paw_ma.fh" 
  
*     **** local variables ****
      logical value
      integer ii,ia,nlm,nlmii,n,m,mabs
      integer nbasis,i2,i1,isgn,nion
      integer w_size
      integer neall

*     **** external functions ****
      integer  psi_ne
      external psi_ne

c     !*** calculate total number of (n,l,m) projectors  ***
      total_nbasis = 0
      do ia=1,ion_nkatm()
c        !-- ion_natm(ia) is number of atoms of kind ia --
         total_nbasis = total_nbasis 
     >                + paw_proj_nbasis(ia)*ion_natm(ia)
      end do


c     !*** allocate storage for paw_ovlp arrays*** 
c     !-- total number of electron states
      neall = (psi_ne(1) + psi_ne(2))

      value =           MA_alloc_get(mt_dcpl,(total_nbasis*neall),
     >                               'coeff',coeff(2),coeff(1))
      value = value.and.MA_alloc_get(mt_dcpl,
     >                               (total_nbasis*neall),
     >                               'toeff',toeff(2),toeff(1))
      value = value.and.MA_alloc_get(mt_int,(total_nbasis*neall),
     >                               'indx2',indx2(2),indx2(1))
      value = value.and.MA_alloc_get(mt_int,(total_nbasis*neall),
     >                               'sgn',sgn(2),sgn(1))

      value = value.and.my_alloc(mt_dcpl,(total_nbasis**2),"w",w)
      if (.not.value) call errquit('paw_ovlp_init: alloc heap',1)


c    !    atom=1                        atom=2
c    !l : 0  1  1  1  2  2  2  2  2 ... 0  1  1  1 ... 
c    !m : 0 -1  0  1 -2 -1  0  1  1 ... 0 -1  0  1 ...
c    !
c    !

c     !*** set sign array (isgn) and location of -m elements(indx2)
      nion = ion_nion()
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
#include "paw_ovlp_data.fh"

      !*** local variables ***
      logical value
      
      value =           MA_free_heap(coeff(2))
      value = value.and.MA_free_heap(toeff(2))
      value = value.and.MA_free_heap(indx2(2))
      value = value.and.MA_free_heap(sgn(2))
      value = value.and.MA_free_heap(w(2))
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

      !*** local variables ***
      logical value
      integer ii,ia,nbasis,nlm,nlmii,neall,npack1
      integer nion,prj,i1,i2,n
      integer indx_psi,indx_prj,indx_coeff
      integer exi(2),tmp(2)
      integer nx,ny,nz
      real*8 scal_factor

      !*** nwpw external functions ****
      integer  psi_ne
      external psi_ne 
      real*8 lattice_omega
      external lattice_omega


c     !*** total number of electron states ***
      neall = (psi_ne(1) + psi_ne(2))

      call Pack_npack(1,npack1)

      !*** allocate temporary space ***
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
c           write(*,*) "toeff:",n,nlmii,'(',nlm,ii,')',
c     >                 dcpl_mb(toeff(1)+indx_coeff)
           end do
        end do
      end do 

      nlmii = 0
      do ii=1,nion
        ia = ion_katm(ii)
        nbasis = paw_proj_nbasis(ia)
        do nlm=1,nbasis
           nlmii = nlmii + 1
           do n=1,neall
             indx_coeff = (nlmii-1)*neall + (n-1)
             write(*,*) "toeff:",n,nlmii,'(',nlm,ii,')',
     >                 dcpl_mb(toeff(1)+indx_coeff)
           end do
        end do
      end do 



      !*** unsrcramble to get true overlap coefficients ***
      do i1=1,total_nbasis*neall
         i2 = int_mb(indx2(1)+i1-1)
         dcpl_mb(coeff(1)+i1-1) 
     >     = 0.5d0*( dcpl_mb(toeff(1)+i1-1)
     >             + int_mb(sgn(1)+i1-1)
     >              *dconjg(dcpl_mb(toeff(1)+i2)))
      end do     


      !*** debug - write out the overlap coefficients ***
      write(*,*)
      write(*,*) "Overlap coefficients:",neall,total_nbasis
      do i2=1,neall
      do i1=1,total_nbasis
         indx_coeff = (i2-1) + (i1-1)*neall
         write(*,*) "nlm n c(n,nlm):",i1,i2,dcpl_mb(coeff(1)+indx_coeff)
      end do
      end do
      write(*,*)
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

      integer nlmii
      integer i,nlm,ia,n
      integer neall,nbasis
      integer indx_coeff
      double complex w_tmp
      
c     **** external functions ****
      integer  psi_ne
      external psi_ne
      
c     !*** total number of electron states ***
      neall = (psi_ne(1) + psi_ne(2))

      nlmii = 0
      do i=1,ion_nion()
        nbasis = paw_proj_nbasis(ia)
        ia = ion_katm(i)
        do nlm=1,nbasis
           nlmii = nlmii + 1
           w_tmp = cmplx(0.0,0.0)
           do n=1,neall
             indx_coeff = (nlmii-1)*neall + (n-1)
             w_tmp = w_tmp + dcpl_mb(coeff(1)+indx_coeff)*
     >        conjg(dcpl_mb(coeff(1)+indx_coeff)) 
              write(*,*) "paw_ovlp:w",w_tmp
           end do
        end do
      end do 

      end
     
