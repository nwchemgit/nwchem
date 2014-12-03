#define TCGMSG


*     ***************************
*     *                         *
*     *    c_semicore_xc_F      *
*     *                         *
*     ***************************

      subroutine c_semicore_xc_F(ispin,xcp,fion)
      implicit none
#include "errquit.fh"
      integer ispin
      real*8  xcp(*)

      real*8 fion(3,*)

#include "bafdecls.fh"
#include "c_semicore_common.fh"

*     **** local variables ****
      logical value
      integer npack0,nfft3d
      integer ii,ia,nx,ny,nz
      real*8  sumx,sumy,sumz
      real*8  scal1,scal2
      integer exi(2),vxcG(2)
      integer Gx(2),Gy(2),Gz(2)
      integer dng(2)
      integer dngx(2),dngy(2),dngz(2)
      integer cngx(2),cngy(2),cngz(2)

*     **** external functions ****
      integer  ion_nion,ion_katm,c_G_indx
      real*8   lattice_omega
      external ion_nion,ion_katm,c_G_indx
      external lattice_omega


      call C3dB_nx(1,nx)
      call C3dB_ny(1,ny)
      call C3dB_nz(1,nz)
      scal1 = 1.0d0/dble(nx*ny*nz)
      scal2 = 1.0d0/lattice_omega()

      call C3dB_nfft3d(1,nfft3d)
      call Cram_npack(0,npack0)

      value = BA_push_get(mt_dcpl,nfft3d,'exi', exi(2), exi(1))
      value = value.and.
     >        BA_push_get(mt_dbl,nfft3d,'vxcG',vxcG(2),vxcG(1))

      value = value.and.
     >        BA_push_get(mt_dbl, nfft3d,'Gx',Gx(2),Gx(1))
      value = value.and.
     >        BA_push_get(mt_dbl, nfft3d,'Gy',Gy(2),Gy(1))
      value = value.and.
     >        BA_push_get(mt_dbl, nfft3d,'Gz',Gz(2),Gz(1))

      value = value.and.
     >        BA_push_get(mt_dcpl, nfft3d,'dng',dng(2),dng(1))
      value = value.and.
     >        BA_push_get(mt_dcpl, nfft3d,'cngx',cngx(2),cngx(1))
      value = value.and.
     >        BA_push_get(mt_dcpl, nfft3d,'cngy',cngy(2),cngy(1))
      value = value.and.
     >        BA_push_get(mt_dcpl, nfft3d,'cngz',cngz(2),cngz(1))

      value = value.and.
     >        BA_push_get(mt_dbl, nfft3d,'dngx',dngx(2),dngx(1))
      value = value.and.
     >        BA_push_get(mt_dbl, nfft3d,'dngy',dngy(2),dngy(1))
      value = value.and.
     >        BA_push_get(mt_dbl, nfft3d,'dngz',dngz(2),dngz(1))

      if (.not. value) 
     > call errquit(' c_semicore_xc_F:out of stack memory',0,MA_ERR)


      call C3dB_rr_Sum(1,xcp(1),
     >                   xcp(1+(ispin-1)*nfft3d),
     >                   dbl_mb(vxcG(1)))
c     call C3dB_r_SMul(1,0.5d0,dcpl_mb(vxcG(1)),dcpl_mb(vxcG(1)))


*     **** define Gx, Gy, and Gz in packed space ****
      call C3dB_t_Copy(1,dbl_mb(c_G_indx(1)),dbl_mb(Gx(1)))
      call C3dB_t_Copy(1,dbl_mb(c_G_indx(2)),dbl_mb(Gy(1)))
      call C3dB_t_Copy(1,dbl_mb(c_G_indx(3)),dbl_mb(Gz(1)))
      call Cram_r_pack(0,dbl_mb(Gx(1)))
      call Cram_r_pack(0,dbl_mb(Gy(1)))
      call Cram_r_pack(0,dbl_mb(Gz(1)))


      do ii=1,ion_nion()
         ia = ion_katm(ii)

         if (log_mb(semicore(1)+ia)) then

*          **** structure factor and local pseudopotential ****
           call cstrfac(ii,dcpl_mb(exi(1)))
           call Cram_c_pack(0,dcpl_mb(exi(1)))

*          **** put sqrt(core-density) at atom position ****
           call Cram_rc_Mul(0,
     >               dbl_mb(ncore(1)+(ia-1)*5*npack0),
     >              dcpl_mb(exi(1)),
     >              dcpl_mb(dng(1)))

           call Cram_rc_iMul(0,dbl_mb(Gx(1)),dcpl_mb(dng(1)),
     >                                       dcpl_mb(cngx(1)))
           call Cram_rc_iMul(0,dbl_mb(Gy(1)),dcpl_mb(dng(1)),
     >                                       dcpl_mb(cngy(1)))
           call Cram_rc_iMul(0,dbl_mb(Gz(1)),dcpl_mb(dng(1)),
     >                                       dcpl_mb(cngz(1)))

*          **** put dng,dngx,dngy,dngz in real space ****
           call Cram_c_unpack(0,dcpl_mb(dng(1)))
           call Cram_c_unpack(0,dcpl_mb(cngx(1)))
           call Cram_c_unpack(0,dcpl_mb(cngy(1)))
           call Cram_c_unpack(0,dcpl_mb(cngz(1)))

           !call C3dB_cr_fft3b(1,dcpl_mb(dng(1)))
           !call C3dB_cr_fft3b(1,dcpl_mb(cngx(1)))
           !call C3dB_cr_fft3b(1,dcpl_mb(cngy(1)))
           !call C3dB_cr_fft3b(1,dcpl_mb(cngz(1)))
           call C3dB_cr_pfft3b(1,0,dcpl_mb(dng(1)))
           call C3dB_cr_pfft3b(1,0,dcpl_mb(cngx(1)))
           call C3dB_cr_pfft3b(1,0,dcpl_mb(cngy(1)))
           call C3dB_cr_pfft3b(1,0,dcpl_mb(cngz(1)))

           call C3dB_ccr_Mul(1,dcpl_mb(dng(1)),
     >                         dcpl_mb(cngx(1)),
     >                         dbl_mb(dngx(1)))
           call C3dB_ccr_Mul(1,dcpl_mb(dng(1)),
     >                         dcpl_mb(cngy(1)),
     >                         dbl_mb(dngy(1)))
           call C3dB_ccr_Mul(1,dcpl_mb(dng(1)),
     >                         dcpl_mb(cngz(1)),
     >                         dbl_mb(dngz(1)))

           call C3dB_rr_dot(1,dbl_mb(dngx(1)),dbl_mb(vxcG(1)),sumx)
           call C3dB_rr_dot(1,dbl_mb(dngy(1)),dbl_mb(vxcG(1)),sumy)
           call C3dB_rr_dot(1,dbl_mb(dngz(1)),dbl_mb(vxcG(1)),sumz)

           fion(1,ii) = fion(1,ii) + sumx*scal2*scal1
           fion(2,ii) = fion(2,ii) + sumy*scal2*scal1
           fion(3,ii) = fion(3,ii) + sumz*scal2*scal1

         end if

      end do

      value = BA_pop_stack(dngz(2))
      value = value.and.BA_pop_stack(dngy(2))
      value = value.and.BA_pop_stack(dngx(2))
      value = value.and.BA_pop_stack(cngz(2))
      value = value.and.BA_pop_stack(cngy(2))
      value = value.and.BA_pop_stack(cngx(2))
      value = value.and.BA_pop_stack(dng(2))
      value = value.and.BA_pop_stack(Gz(2))
      value = value.and.BA_pop_stack(Gy(2))
      value = value.and.BA_pop_stack(Gx(2))
      value = value.and.BA_pop_stack(vxcG(2))
      value = value.and.BA_pop_stack(exi(2))
      if (.not. value) 
     > call errquit('c_semicore_force:error popping stack',0, MA_ERR)


      return
      end


c $Id$
