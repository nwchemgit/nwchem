*
* $Id$
*

*     *************************************************
*     *                                               *
*     *              nwpw_compcharge_init             *
*     *                                               *
*     *************************************************
      subroutine nwpw_compcharge_init(nion0,nkatm0,
     >                                nprj,nbasis,psp_type,
     >                                lmax0,sigma,
     >                                nprj_max,l_prj,m_prj,b_prj,
     >                                comp_charge_matrix)
      implicit none
      integer nion0,nkatm0
      integer nprj(*),nbasis(*),psp_type(*),lmax0(*)
      real*8 sigma(*)
      integer nprj_max
      integer l_prj(nprj_max,*),m_prj(nprj_max,*),b_prj(nprj_max,*)
      integer comp_charge_matrix(*)

#include "mafdecls.fh"
#include "errquit.fh"
#include "nwpw_compcharge.fh"

*     ***** local variables ****
      logical value
      integer k,l,m,lm,ia,Gindx,Gall(2),G(3),tsize
      integer iasize,i,j,li,mi,lj,mj,bi,bj
      integer indx,comp_charge_tag,comp_indx
      real*8  scal,gg,fourpioveromega,sumall
      real*8  taunt

*     **** external functions ****
      integer  control_code,G_indx,c_G_indx
      real*8   double_factorial,lattice_omega,ewald_rcut
      external control_code,G_indx,c_G_indx
      external double_factorial,lattice_omega,ewald_rcut
      integer  psi_data_get_ptr
      external psi_data_get_ptr
      real*8   nwpw_gaunt
      external nwpw_gaunt

      nion         = nion0
      nkatm        = nkatm0
      sigma_smooth = ewald_rcut()
      value = MA_alloc_get(mt_int,nkatm,"mult_l",mult_l(2),mult_l(1))
      value = value.and.
     >        MA_alloc_get(mt_int,nkatm,"lm_size",lm_size(2),lm_size(1))
      if (.not.value)
     >   call errquit("nwpw_compcharge_end:allocate memory",0,MA_ERR)

      !**** band structure code ****
      if ((control_code().eq.5).or.
     >    (control_code().eq.13).or.
     >    (control_code().eq.14)) then
         isgamma = .false.
      else
         isgamma = .true.
      end if

      if (isgamma) then
         call D3dB_nfft3d(1,nfft3d)
         call Pack_npack(0,npack0)
         Gindx = G_indx(1)
      else
         call C3dB_nfft3d(1,nfft3d)
         call Cram_npack(0,npack0)
         Gindx = c_G_indx(1)
      end if

      mult_l_max = 0
      do ia=1,nkatm
         if (psp_type(ia).eq.4) then
            int_mb(mult_l(1)+ia-1)  = 2*lmax0(ia)
            int_mb(lm_size(1)+ia-1) = (2*lmax0(ia)+1)**2
            if (mult_l_max.lt.(2*lmax0(ia))) mult_l_max = 2*lmax0(ia)
         end if
      end do

      !*** allocate gk_smooth, gk,and glm ***
      lm_size_max = (mult_l_max+1)**2
      value = MA_alloc_get(mt_dbl,npack0,
     >                     "gk_smooth",gk_smooth(2),gk_smooth(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,npack0*nkatm,"gk",gk(2),gk(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,npack0*lm_size_max,
     >                     "glm",glm(2),glm(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,2*nion*lm_size_max,
     >                     "Qlm",Qlm(2),Qlm(1))
      value = value.and.
     >        MA_push_get(mt_dbl,3*nfft3d,"Gall",Gall(2),Gall(1))
      if (.not.value)
     >   call errquit("nwpw_compcharge_end:allocate memory",1,MA_ERR)

      call dcopy(2*nion*lm_size_max,0.0d0,0,dbl_mb(Qlm(1)),1)
      call dcopy(3*nfft3d,dbl_mb(Gindx),1,dbl_mb(Gall(1)),1)
      G(1) = Gall(1)
      G(2) = Gall(1)+nfft3d
      G(3) = Gall(1)+2*nfft3d
      if (isgamma) then
         call Pack_t_pack(0,dbl_mb(G(1)))
         call Pack_t_pack(0,dbl_mb(G(2)))
         call Pack_t_pack(0,dbl_mb(G(3)))
      else
         call Cram_r_pack(0,dbl_mb(G(1)))
         call Cram_r_pack(0,dbl_mb(G(2)))
         call Cram_r_pack(0,dbl_mb(G(3)))
      end if

      lm = 0
      do l =0,mult_l_max

*        **** define  |k|**l / (2*l+1)!! ****
         scal = 1.0d0/double_factorial(2*l+1)
         if (l.gt.0) then
            do k=1,npack0
               gg = dbl_mb(G(1)+k-1)**2 
     >            + dbl_mb(G(2)+k-1)**2 
     >            + dbl_mb(G(3)+k-1)**2
               dbl_mb(gk(1)+k-1) = scal*dsqrt(gg)**l
            end do
         else
            call dcopy(npack0,scal,0,dbl_mb(gk(1)),1)
         end if

*        **** define glm =  |k|**l * Tlm(k)/ (2*l+1)!! ****
         do m =-l,l
            call Tesseral3_vector_lm(l,m,npack0,
     >                               dbl_mb(G(1)),
     >                               dbl_mb(G(2)),
     >                               dbl_mb(G(3)),
     >                               dbl_mb(glm(1)+lm*npack0))
            do k=1,npack0
               dbl_mb(glm(1)+lm*npack0+k-1) 
     >            = dbl_mb(glm(1)+lm*npack0+k-1)*dbl_mb(gk(1)+k-1)
            end do

            lm = lm + 1
         end do

      end do

*     **** define gk_smooth(k)  = 4*pi * Exp[-k*k*sigma_smooth**2 / 4] ****
      scal            = 0.25d0*sigma_smooth**2
      fourpioveromega = 16.0d0*datan(1.0d0)/lattice_omega()
      do k=1,npack0
         gg = dbl_mb(G(1)+k-1)**2 
     >      + dbl_mb(G(2)+k-1)**2 
     >      + dbl_mb(G(3)+k-1)**2
         dbl_mb(gk_smooth(1)+k-1) = fourpioveromega*dexp(-gg*scal)
      end do

*     **** define gk(k,ia)  = 4*pi * Exp[-k*k*sigma(ia**2 / 4] ****
      do ia=1,nkatm
         if (psp_type(ia).eq.4) then
            do k=1,npack0
               scal = 0.25d0*sigma(ia)**2
               gg   = dbl_mb(G(1)+k-1)**2 
     >              + dbl_mb(G(2)+k-1)**2 
     >              + dbl_mb(G(3)+k-1)**2
               dbl_mb(gk(1)+(ia-1)*npack0+k-1) 
     >              = fourpioveromega*dexp(-gg*scal)
            end do
         else
            call dcopy(npack0,0.0d0,0,dbl_mb(gk(1)+(ia-1)*npack0),1)
         end if
      end do
      if (.not.MA_pop_stack(Gall(2)))
     >   call errquit("nwpw_compcharge_init: pop stack",0,MA_ERR)



      if (.not.MA_alloc_get(mt_int,nkatm,"nindx_Tndiff",
     >                      nindx_Tndiff(2),nindx_Tndiff(1)))
     >  call errquit("nwpw_compcharge_init:allocate fail",3,MA_ERR)

      if (.not.MA_alloc_get(mt_int,nkatm,"shift_Tndiff",
     >                      shift_Tndiff(2),shift_Tndiff(1)))
     >  call errquit("nwpw_compcharge_init:allocate fail",3,MA_ERR)

      indx = 0
      do ia=1,nkatm
         int_mb(shift_Tndiff(1)+ia-1) = indx
         if (psp_type(ia).eq.4) then
            comp_charge_tag = comp_charge_matrix(ia)

            do l=0,int_mb(mult_l(1)+ia-1)
               comp_indx = psi_data_get_ptr(comp_charge_tag,l+1)
               do m=-l,l

                  do j=1,nprj(ia)
                     lj = l_prj(j,ia)
                     mj = l_prj(j,ia)
                     bj = b_prj(j,ia)

                     do i=1,nprj(ia)
                        li = l_prj(i,ia)
                        mi = l_prj(i,ia)
                        bi = b_prj(i,ia)
                        taunt = nwpw_gaunt(.false.,l,m,li,mi,lj,mj)
     >                       *dbl_mb(comp_indx+(bi-1)+(bj-1)*nbasis(ia))
                        if (dabs(taunt).gt.1.0d-11) then
                           indx = indx + 1
                        end if
                     end do
                  end do

               end do
            end do
         end if
         int_mb(nindx_Tndiff(1)+ia-1)=indx-int_mb(shift_Tndiff(1)+ia-1)
      end do

      value = MA_alloc_get(mt_int,indx,"lm_Tndiff",
     >                     lm_Tndiff(2),lm_Tndiff(1))
      value = value.and.
     >        MA_alloc_get(mt_int,indx,"iprj_Tndiff",
     >                     iprj_Tndiff(2),iprj_Tndiff(1))
      value = value.and.
     >        MA_alloc_get(mt_int,indx,"jprj_Tndiff",
     >                     jprj_Tndiff(2),jprj_Tndiff(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,indx,"coeff_Tndiff",
     >                     coeff_Tndiff(2),coeff_Tndiff(1))
      if (.not.value)
     > call errquit(
     >  "nwpw_compcharge_init: error allocating work arrays",0,MA_ERR)

      indx = 0
      do ia=1,nkatm
         if (psp_type(ia).eq.4) then
            comp_charge_tag = comp_charge_matrix(ia)

            lm = 0
            do l=0,int_mb(mult_l(1)+ia-1)
               comp_indx = psi_data_get_ptr(comp_charge_tag,l+1)
               do m=-l,l

                  do j=1,nprj(ia)
                     lj = l_prj(j,ia)
                     mj = l_prj(j,ia)
                     bj = b_prj(j,ia)

                     do i=1,nprj(ia)
                        li = l_prj(i,ia)
                        mi = l_prj(i,ia)
                        bi = b_prj(i,ia)
                        taunt = nwpw_gaunt(.false.,l,m,li,mi,lj,mj)
     >                       *dbl_mb(comp_indx+(bi-1)+(bj-1)*nbasis(ia))
                        if (dabs(taunt).gt.1.0d-11) then
                           int_mb(lm_Tndiff(1)+indx)   = lm
                           int_mb(iprj_Tndiff(1)+indx) = i
                           int_mb(jprj_Tndiff(1)+indx) = j
                           dbl_mb(coeff_Tndiff(1)+indx) = taunt
                           indx = indx + 1
                        end if
                     end do
                  end do
                  lm = lm + 1
               end do
            end do
         end if
      end do


      return
      end


*     *************************************************
*     *                                               *
*     *              nwpw_compcharge_end              *
*     *                                               *
*     *************************************************
      subroutine nwpw_compcharge_end()
      implicit none

#include "mafdecls.fh"
#include "errquit.fh"
#include "nwpw_compcharge.fh"

      logical value
      value =           MA_free_heap(mult_l(2))
      value = value.and.MA_free_heap(lm_size(2))
      value = value.and.MA_free_heap(gk_smooth(2))
      value = value.and.MA_free_heap(gk(2))
      value = value.and.MA_free_heap(glm(2))
      value = value.and.MA_free_heap(Qlm(2))
      value = value.and.MA_free_heap(shift_Tndiff(2))
      value = value.and.MA_free_heap(nindx_Tndiff(2))
      value = value.and.MA_free_heap(lm_Tndiff(2))
      value = value.and.MA_free_heap(iprj_Tndiff(2))
      value = value.and.MA_free_heap(jprj_Tndiff(2))
      value = value.and.MA_free_heap(coeff_Tndiff(2))
      if (.not.value)
     >   call errquit("nwpw_compcharge_end: freeing heap",0,MA_ERR)
      return
      end 

*     *************************************************
*     *                                               *
*     *              nwpw_compcharge_gen_Qlm          *
*     *                                               *
*     *************************************************

*     This routine computes Qlm for atom ii.

      subroutine nwpw_compcharge_gen_Qlm(ii,ia,ispin,ne,nprj,sw1)
      implicit none
      integer ii,ia
      integer ispin,ne(2),nprj
      real*8  sw1(ne(1)+ne(2),nprj)

#include "mafdecls.fh"
#include "errquit.fh"
#include "nwpw_compcharge.fh"

*     **** local variables ****
      integer n,k,ms,lm,iprj,jprj,shift,indx,n1(2),n2(2)
      real*8  scal,coeff,w

*     **** external functions ****
      real*8   lattice_omega
      external lattice_omega

      call dcopy(ispin*nion*lm_size_max,0.0d0,0,dbl_mb(Qlm(1)+indx),1)
      n1(1) = 1
      n1(2) = ne(1)+1
      n2(1) = ne(1)
      n2(2) = ne(1)+ne(2)
      scal = 1.0d0/lattice_omega()

      shift = int_mb(shift_Tndiff(1)+ia-1)
      do k=1,int_mb(nindx_Tndiff(1)+ia-1)
         lm    = int_mb(lm_Tndiff(1)+shift+k-1)
         iprj  = int_mb(iprj_Tndiff(1)+shift+k-1)
         jprj  = int_mb(jprj_Tndiff(1)+shift+k-1)
         coeff = scal*dbl_mb(coeff_Tndiff(1)+shift+k-1)
         do ms=1,ispin
            w = 0.0d0
            do n=n1(ms),n2(ms)
               w = w + sw1(n,iprj)*sw1(n,jprj)
            end do
            indx = (ms-1)*nion*lm_size_max+(ii-1)*lm_size_max + lm
            dbl_mb(Qlm(1)+indx) = dbl_mb(Qlm(1)+indx) + coeff*w
         end do
      end do


      return
      end 

c*     *************************************************
c*     *                                               *
c*     *              nwpw_compcharge_Qlm_sub          *
c*     *                                               *
c*     *************************************************
c      subroutine nwpw_compcharge_Qlm_sub(ispin,ne,nprj,Tndiff,sw1,qlm)
c      implicit none
c      integer ispin,ne(2),nprj
c      real*8  Tndiff(nprj,nprj)
c      real*8  sw1(ne(1)+ne(2),nprj)
c      real*8  qlm(*)
c
c*     **** local variables ****
c      integer ms,i,j,n,n1(2),n2(2)
c      real*8 psum
c
c      n1(1) = 1
c      n1(2) = ne(1)+1
c      n2(1) = ne(1)
c      n2(2) = ne(1)+ne(2)
c      do ms=1,ispin
c         qlm(ms) = 0.0d0
c      end do
c      do j=1,nprj
c         do i=1,nprj
c            do ms=1,ispin
c               psum = 0.0d0
c               do n=n1(ms),n2(ms)
c                  psum = psum + sw1(n,i)*sw1(n,j)
c               end do
c               qlm(ms) = qlm(ms) + Tndiff(i,j)*psum
c            end do
c         end do
c      end do
c      return
c      end


*     *************************************************
*     *                                               *
*     *         nwpw_compcharge_gen_rgaussian         *
*     *                                               *
*     *************************************************
      subroutine nwpw_compcharge_gen_rgaussian(l,sigma,nr,r,gl)
      implicit none
      integer  l
      integer  nr
      real*8 sigma
      double precision r(*)
      double precision gl(*)

      !*** local variables ***
      integer i
      real*8  c,pi

*     ***** external functions ****
      real*8   double_factorial
      external double_factorial

      pi = 4.0d0*datan(1.0d0)
      c  = 2**(l+2)/(dsqrt(pi)*double_factorial(2*l+1)*sigma**(2*l+3))

      !*** this fixes possible underflow error ***
      do i=1,nr
         gl(i) = 0.0d0
      end do

      do i=1,nr
        if (dabs(r(i)) .lt. (8.0d0*sigma))
     >     gl(i) = c*r(i)**l*dexp(-(r(i)/sigma)**2)
      end do

      return
      end

*     *************************************************
*     *                                               *
*     *     nwpw_compcharge_gen_dn_cmp_smooth         *
*     *                                               *
*     *************************************************
      subroutine nwpw_compcharge_gen_dn_cmp_smooth(ms,dng_smooth)
      implicit none
      integer ms
      complex*16 dng_smooth(*)

#include "mafdecls.fh"
#include "errquit.fh"
#include "nwpw_compcharge.fh"

*     **** local variables ****
      logical value
      integer ii,ia,lm,l,m,indx,QTlm(2),tmp(2),exi(2)
      real*8  q,sumall

*     **** external functions ****
      integer  ion_katm
      external ion_katm

      value = MA_push_get(mt_dcpl,npack0,"QTlm",QTlm(2),QTlm(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,npack0,"tmp",tmp(2),tmp(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,npack0,"exi",exi(2),exi(1))
      if (.not.value)
     >  call errquit("nwpw_compcharge_gen_dn_cmp_smooth:stack",0,MA_ERR)

      call dcopy(2*npack0,0.0d0,0,dng_smooth,1)
      do ii=1,nion
         ia = ion_katm(ii)
         call strfac_pack(0,ii,dcpl_mb(exi(1)))

         call dcopy(2*npack0,0.0d0,0,dcpl_mb(QTlm(1)),1)
         lm = 0
         do l=0,int_mb(mult_l(1)+ia-1)
            do m=-l,l
               indx = (ms-1)*nion*lm_size_max + (ii-1)*lm_size_max + lm 
               q = dbl_mb(Qlm(1)+indx)
               if (mod(l,2).eq.0) then
                  call Pack_tc_Mul(0,dbl_mb(glm(1)+lm*npack0),
     >                               dcpl_mb(exi(1)),
     >                               dcpl_mb(tmp(1)))
               else
                  call Pack_tc_iMul(0,dbl_mb(glm(1)+lm*npack0),
     >                                dcpl_mb(exi(1)),
     >                                dcpl_mb(tmp(1)))
               end if

               call Pack_cc_daxpy(0,q,dcpl_mb(tmp(1)),
     >                                dcpl_mb(QTlm(1)))

               lm = lm + 1
            end do
         end do
         call Pack_cc_Sum2(0,dcpl_mb(QTlm(1)),dng_smooth)

      end do
      call Pack_tc_Mul2(0,dbl_mb(gk_smooth(1)),dng_smooth)


      value =           MA_pop_stack(exi(2))
      value = value.and.MA_pop_stack(tmp(2))
      value = value.and.MA_pop_stack(QTlm(2))
      if (.not.value)
     >  call errquit("nwpw_compcharge_gen_dn_cmp_smooth:stack",1,MA_ERR)

      return
      end


