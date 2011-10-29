*
* $Id: nwpw_compcharge.F 19707 2010-10-29 17:59:36Z d3y133 $
*

*     *************************************************
*     *                                               *
*     *              nwpw_compcharge_init             *
*     *                                               *
*     *************************************************
      subroutine nwpw_compcharge_init(nion0,nkatm0,nprj,nbasis,psp_type,
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
      real*8  scal,gg,fourpioveromega
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
     >        MA_alloc_get(mt_dbl,nion*lm_size_max,"Qlm",Qlm(2),Qlm(1))
      value = value.and.
     >        MA_push_get(mt_dbl,3*nfft3d,"Gall",Gall(2),Gall(1))
      if (.not.value)
     >   call errquit("nwpw_compcharge_end:allocate memory",1,MA_ERR)

      call dcopy(nion*lm_size_max,0.0d0,0,dbl_mb(Qlm(1)),1)
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




      if (.not.MA_alloc_get(mt_int,nkatm,"Tndiff_shift",
     >                     Tndiff_shift(2),Tndiff_shift(1)))
     >  call errquit("nwpw_compcharge_init:allocate fail",3,MA_ERR)

      tsize = 0
      do ia=1,nkatm
         int_mb(Tndiff_shift(1)+ia-1) = tsize
         if (psp_type(ia).eq.4) then
            iasize = int_mb(lm_size(1)+ia-1)*nprj(ia)*nprj(ia)
            tsize  = tsize + iasize
         end if
      end do
      if (.not.MA_alloc_get(mt_dbl,tsize,"Tndiff",
     >                     Tndiff(2),Tndiff(1)))
     >  call errquit("nwpw_compcharge_init:allocate fail",4,MA_ERR)

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

                        taunt = nwpw_gaunt(.true.,l,m,li,mi,lj,mj)
                        dbl_mb(Tndiff(1)+indx)
     >                    = taunt*dbl_mb(comp_indx+bi+(bj-1)*nbasis(ia))

c                  write(*,'(A,7I4,3E13.3)') "i,j,l,m,indx,bi,bj",
c     >                         i,j,l,m,indx,
c     >                         bi,bj,
c     >                         taunt,dbl_mb(Tndiff(1)+indx),
c     >                     dbl_mb(comp_indx+(bi-1)+(bj-1)*nbasis(ia))

                        indx = indx + 1
                     end do
                  end do

               end do
            end do
         end if
      end do

c      lm = 0
c      ia = 1
c      do l=0,int_mb(mult_l(1)+ia-1)
c      do m=-l,l
c          write(*,*) "Tndiff, l,m=",l,m
c          do i=1,nprj(ia)
c             indx = (i-1) + lm*nprj(ia)*nprj(ia)
c             write(*,'(12E13.3)') 
c     >       (dbl_mb(Tndiff(1)+indx+(j-1)*nprj(ia)),j=1,nprj(ia))
c          end do
c          write(*,*)
c          lm = lm + 1
c      end do
c      end do

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
      value = value.and.MA_free_heap(Tndiff_shift(2))
      value = value.and.MA_free_heap(Tndiff(2))
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

      subroutine nwpw_compcharge_gen_Qlm(ii,ia,nprj,neall,sw1)
      implicit none
      integer ii,ia
      integer nprj,neall
      real*8  sw1(neall,nprj)

#include "mafdecls.fh"
#include "errquit.fh"
#include "nwpw_compcharge.fh"

*     **** local variables ****
      integer l,m,lm,shiftia
      real*8  scal

*     **** external functions ****
      real*8   nwpw_compcharge_Qlm_sub,lattice_omega
      external nwpw_compcharge_Qlm_sub,lattice_omega

      shiftia = int_mb(Tndiff_shift(1)+ia-1)

c      write(*,*) "sw1="
c      do m=1,neall
c         write(*,'(10E13.3)') (sw1(m,l),l=1,nprj)
c      end do
c      write(*,*)
      scal = 1.0d0/lattice_omega()

      lm = 0
      do l=0,int_mb(mult_l(1)+ia-1)
         do m=-l,l
            dbl_mb(Qlm(1)+(ii-1)*lm_size_max + lm) 
     >         = scal*nwpw_compcharge_Qlm_sub(nprj,neall,
     >                    dbl_mb(Tndiff(1)+shiftia+lm*nprj*nprj),
     >                    sw1)

            write(*,*) "l,m,qlm,Tndiff=",l,m,
     >                 dbl_mb(Qlm(1)+(ii-1)*lm_size_max+lm),
     >                 dbl_mb(Tndiff(1)+shiftia+(lm-1)*nprj*nprj)
     
            lm = lm + 1
         end do
      end do

      return
      end 

*     *************************************************
*     *                                               *
*     *              nwpw_compcharge_Qlm_sub          *
*     *                                               *
*     *************************************************
      real*8 function nwpw_compcharge_Qlm_sub(nprj,neall,Tndiff,sw1)
      implicit none
      integer nprj,neall
      real*8  Tndiff(nprj,nprj)
      real*8  sw1(neall,nprj)

*     **** local variables ****
      integer i,j,n
      real*8 qlm,psum

      qlm = 0.0d0
      do j=1,nprj
         do i=1,nprj
            psum = 0.0d0
            do n=1,neall
               psum = psum + sw1(n,i)*sw1(n,j)
            end do
            qlm = qlm + Tndiff(i,j)*psum
         end do
      end do

      nwpw_compcharge_Qlm_sub = qlm
      return
      end


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




