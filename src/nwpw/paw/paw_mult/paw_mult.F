*
* $Id: paw_mult.F,v 1.2 2003-02-17 01:06:23 bylaska Exp $
*


!**************************************************
!
!       Name: paw_mult_init
!
!       Purpose: initializes  paw_mult
!
!       Created:        2/16/2003
!**************************************************
      subroutine paw_mult_init()
      implicit none

#include "mafdecls.fh"
#include "paw_mult_data.fh"
#include "paw_geom.fh"
#include "paw_ma.fh"
#include "paw_basis.fh"

      !**** local variables ***
      logical ok
      integer k,nfft3d,npack0,nion,nkatm
      integer lm_size
      integer ia,ii
      integer l,m
      integer Gx(2),Gy(2),Gz(2),Ylm(2)
      real*8  scal,gg,fourpi
      complex*16 itol,cscal
    
      !**** external functions ****
      integer  control_ncut,G_indx
      real*8   control_rcut,double_factorial
      external control_ncut,G_indx
      external control_rcut,double_factorial


      fourpi = 16.0d0*datan(1.0d0)

      !*** allocate paw mult memory from heap ***
      call D3dB_nfft3d(1,nfft3d)
      nion = ion_nion()
      nkatm = ion_nkatm()
      call Pack_npack(0,npack0)
      lm_size = (paw_basis_max_mult_l()+1)**2

      !*** allocate gk_smooth,gk,and glm ***
      ok =        my_alloc(mt_dbl,npack0,"gk_smooth",gk_smooth)
      ok = ok.and.my_alloc(mt_dbl,npack0*nkatm,"gk",gk)
      ok = ok.and.my_alloc(mt_dcpl,npack0*lm_size,"g_lm",g_lm)

      !*** allocate paw mult arrays ***
      ok = ok.and.my_alloc(mt_int,nion,"i_v_mult",i_v_mult) !** same as i_paw_qlm **
      if (.not.ok) 
     > call errquit("paw_mult_init: out of heap memory",0)

      lm_size = 0
      do ii=1,nion
        ia = ion_katm(ii)

        int_mb(i_v_mult(1) + ii - 1) = lm_size
        lm_size = lm_size + (paw_basis_mult_l(ia) + 1)**2
      end do
      ok = my_alloc(mt_dcpl,lm_size,"v_mult",v_mult)
      if (.not.ok) 
     > call errquit("paw_mult_init: out of heap memory",1)



      !**** initialize sigma_smooth and ncut***
      sigma_smooth = control_rcut()
      ncut         = control_ncut()


      !**** initialize gk_smooth, gk, and g_lm ****

      !**** allocate stack memory ****
      ok = MA_push_get(mt_dbl,nfft3d,'Gx',Gx(2),Gx(1))
      ok = ok.and.
     >     MA_push_get(mt_dbl,nfft3d,'Gy',Gy(2),Gy(1))
      ok = ok.and.
     >     MA_push_get(mt_dbl,nfft3d,'Gz',Gz(2),Gz(1))
      ok = ok.and.
     >     MA_push_get(mt_dcpl,npack0,'Ylm',Ylm(2),Ylm(1))
      if (.not.ok) 
     > call errquit("paw_mult_init: out of stack memory",2)

      call dcopy(nfft3d,dbl_mb(G_indx(1)),1,dbl_mb(Gx(1)))
      call dcopy(nfft3d,dbl_mb(G_indx(2)),1,dbl_mb(Gy(1)))
      call dcopy(nfft3d,dbl_mb(G_indx(3)),1,dbl_mb(Gz(1)))
      call Pack_t_pack(0,dbl_mb(Gx(1)))
      call Pack_t_pack(0,dbl_mb(Gy(1)))
      call Pack_t_pack(0,dbl_mb(Gz(1)))

      !**** g_lm ****
      itol = dcmplx(1.0d0,0.0d0)
      do l=0,paw_basis_max_mult_l()
        cscal = itol/double_factorial(l)
        itol = itol*dcmplx(0.0d0,-1.0d0)

        do m = -l,l

          !*** generate Ylm ***
          call spher_harmonics_generate(l,m,npack0,
     >                      dbl_mb(Gx(1)),
     >                      dbl_mb(Gy(1)),
     >                      dbl_mb(Gz(1)),
     >                      dcpl_mb(Ylm(1)))

          !**** (-i)**l*Ylm(k)*|k|**l/(2l+1)!! ****
          do k=1,npack0
            gg = dbl_mb(Gx(1)+k-1)**2
     >         + dbl_mb(Gy(1)+k-1)**2
     >         + dbl_mb(Gz(1)+k-1)**2
            dcpl_mb(g_lm(1)+k-1)=cscal
     >                          *dsqrt(gg)**l
     >                          *dcpl_mb(Ylm(1)+k-1)
     >
          end do !*k*

        end do !*m*
      end do !*l*

      !**** gk_smooth and gk ****
      do k=1,npack0
         gg = dbl_mb(Gx(1)+k-1)**2
     >      + dbl_mb(Gy(1)+k-1)**2
     >      + dbl_mb(Gz(1)+k-1)**2
         scal = 0.25d0 * sigma_smooth**2
         dbl_mb(gk_smooth(1)+k-1) = fourpi*dexp(-scal*gg)
         do ia=1,nkatm
           scal = 0.25d0 * paw_basis_sigma(ia)**2
           dbl_mb(gk(1)+k-1+(ia-1)*npack0) = fourpi*dexp(-scal*gg)
         end do !*ia*
      end do !*k*

      !**** deallocate stack memory ****
      ok =        MA_pop_stack(Ylm(2))
      ok = ok.and.MA_pop_stack(Gz(2))
      ok = ok.and.MA_pop_stack(Gy(2))
      ok = ok.and.MA_pop_stack(Gx(2))
      if (.not.ok) 
     > call errquit("paw_mult_init: out of stack memory",2)


      return
      end


*************************************************
!
!       Name: paw_mult_end
!
!       Purpose: deallocates heap memory
!
!       Created:        2/16/2003
!**************************************************
      subroutine paw_mult_end()
      implicit none

#include "paw_mult_data.fh"
#include "paw_ma.fh"

      !**** local varables ****
      logical ok

      ok =        my_dealloc(i_v_mult)
      ok = ok.and.my_dealloc(v_mult)
      ok = ok.and.my_dealloc(g_lm)
      ok = ok.and.my_dealloc(gk)
      ok = ok.and.my_dealloc(gk_smooth)
      if (.not.ok)
     > call errquit("paw_mult_end: error freeing heap",0)

      return
      end 




