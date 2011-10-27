*
* $Id: nwpw_compcharge.F 19707 2010-10-29 17:59:36Z d3y133 $
*

*     *************************************************
*     *                                               *
*     *              nwpw_compcharge_init             *
*     *                                               *
*     *************************************************
      subroutine nwpw_compcharge_init()
      implicit none

#include "mafdecls.fh"
#include "errquit.fh"
#include "nwpw_compcharge.fh"

*     ***** local variables ****
      logical value
      integer k,l,m,lm,ia,Gindx,Gall(2),G(3),lmax0,psp_type
      real*8  scal,gg,fourpioveromega

*     **** external functions ****
      integer  control_code,G_indx,c_G_indx,ion_nkatm_qm
      integer  psp_psp_type,cpsp_psp_type,psp_lmax,cpsp_lmax
      real*8   double_factorial,lattice_omega,ewald_rcut,psp_sigma
      external control_code,G_indx,c_G_indx,ion_nkatm_qm
      external psp_psp_type,cpsp_psp_type,psp_lmax,cpsp_lmax
      external double_factorial,lattice_omega,ewald_rcut,psp_sigma

      nkatm        = ion_nkatm_qm()
      sigma_smooth = ewald_rcut()

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
         lmax = 0
         do ia=1,nkatm
            if (psp_psp_type(ia).eq.4) then
               lmax0 = psp_lmax(ia)
               if (lmax0.gt.lmax) lmax = lmax0
            end if
         end do
      else
         call C3dB_nfft3d(1,nfft3d)
         call Cram_npack(0,npack0)
         Gindx = c_G_indx(1)
         lmax = 0
         do ia=1,nkatm
            if (cpsp_psp_type(ia).eq.4) then
               lmax0 = cpsp_lmax(ia)
               if (lmax0.gt.lmax) lmax = lmax0
            end if
         end do
      end if

      !*** allocate gk_smooth, gk,and glm ***
      lm_size = (lmax+1)**2
      value = MA_alloc_get(mt_dbl,npack0,
     >                     "gk_smooth",gk_smooth(2),gk_smooth(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,npack0*nkatm,"gk",gk(2),gk(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,npack0*lm_size,"glm",glm(2),glm(1))
      value = value.and.
     >        MA_push_get(mt_dbl,3*nfft3d,"Gall",Gall(2),Gall(1))
      if (.not.value)
     >   call errquit("nwpw_compcharge_end:allocate memory",1,MA_ERR)

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
      do l =0,lmax

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
         psp_type = psp_psp_type(ia)
         if (psp_type.eq.4) then
            do k=1,npack0
               scal = 0.25d0*psp_sigma(ia)**2
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
      value =           MA_free_heap(gk_smooth(2))
      value = value.and.MA_free_heap(gk(2))
      value = value.and.MA_free_heap(glm(2))
      if (.not.value)
     >   call errquit("nwpw_compcharge_end: freeing heap",0,MA_ERR)

      return
      end 
