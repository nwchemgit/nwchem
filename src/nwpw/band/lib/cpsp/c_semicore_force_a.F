#define TCGMSG


*     ***************************
*     *                         *
*     *    semicore_xc_F        *
*     *                         *
*     ***************************

      subroutine semicore_xc_F(ispin,xcp,fion)
      implicit none
#include "errquit.fh"
      integer ispin
      real*8  xcp(*)

      real*8 fion(3,*)

#include "bafdecls.fh"

*     **** semicore common block ****
c     real*8  ncore(nfft3d,nkatmx),rcore(nkatmx)
c     logocal semicore(0:nkatmx)
      integer ncore(2),rcore(2)
      integer semicore(2)
      common / ccore / ncore,rcore,semicore

*     **** local variables ****
      logical value
      integer npack0,nfft3d,n2ft3d
      integer ii,ia,nx,ny,nz
      real*8  sumx,sumy,sumz
      real*8  scal1,scal2
      integer exi(2),vxcG(2)
      integer tmpx(2),tmpy(2),tmpz(2)

*     **** external functions ****
      integer  ion_nion,ion_katm
      real*8   lattice_omega
      external ion_nion,ion_katm
      external lattice_omega


      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)
      scal1 = 1.0d0/dble(nx*ny*nz)
      scal2 = 1.0d0/lattice_omega()

      call D3dB_nfft3d(1,nfft3d)
      call D3dB_n2ft3d(1,n2ft3d)
      call Pack_npack(0,npack0)

      value = BA_push_get(mt_dcpl,nfft3d,'exi', exi(2), exi(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,nfft3d,'vxcG',vxcG(2),vxcG(1))
      value = value.and.
     >        BA_push_get(mt_dcpl, npack0,'tmpx',tmpx(2),tmpx(1))
      value = value.and.
     >        BA_push_get(mt_dcpl, npack0,'tmpy',tmpy(2),tmpy(1))
      value = value.and.
     >        BA_push_get(mt_dcpl, npack0,'tmpz',tmpz(2),tmpz(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)


      call dcopy(n2ft3d,0.0d0,0,dcpl_mb(vxcG(1)),1)
      call D3dB_rr_Sum(1,xcp(1),
     >                   xcp(1+(ispin-1)*n2ft3d),
     >                   dcpl_mb(vxcG(1)))
      call D3dB_rc_fft3f(1,dcpl_mb(vxcG(1)))
      call Pack_c_pack(0,dcpl_mb(vxcG(1)))


      do ii=1,ion_nion()
         ia = ion_katm(ii)

         if (log_mb(semicore(1)+ia)) then

*          **** structure factor and local pseudopotential ****
           call strfac(ii,dcpl_mb(exi(1)))
           call Pack_c_pack(0,dcpl_mb(exi(1)))

*          **** put (core-density) at atom position ****
           call Pack_tc_iMul(0,
     >               dbl_mb(ncore(1)+(ia-1)*5*npack0+2*npack0),
     >              dcpl_mb(exi(1)),
     >              dcpl_mb(tmpx(1)))
           call Pack_tc_iMul(0,
     >               dbl_mb(ncore(1)+(ia-1)*5*npack0+3*npack0),
     >              dcpl_mb(exi(1)),
     >              dcpl_mb(tmpy(1)))
           call Pack_tc_iMul(0,
     >               dbl_mb(ncore(1)+(ia-1)*5*npack0+4*npack0),
     >              dcpl_mb(exi(1)),
     >              dcpl_mb(tmpz(1)))

           call Pack_cc_dot(0,dcpl_mb(tmpx(1)),dcpl_mb(vxcG(1)),sumx)
           call Pack_cc_dot(0,dcpl_mb(tmpy(1)),dcpl_mb(vxcG(1)),sumy)
           call Pack_cc_dot(0,dcpl_mb(tmpz(1)),dcpl_mb(vxcG(1)),sumz)

           fion(1,ii) = fion(1,ii) + sumx*dsqrt(scal1*scal2)
           fion(2,ii) = fion(2,ii) + sumy*dsqrt(scal1*scal2)
           fion(3,ii) = fion(3,ii) + sumz*dsqrt(scal1*scal2)
         end if

      end do

      value = BA_pop_stack(tmpz(2))
      value = BA_pop_stack(tmpy(2))
      value = BA_pop_stack(tmpx(2))
      value = BA_pop_stack(vxcG(2))
      value = BA_pop_stack(exi(2))


      return
      end


c $Id$
