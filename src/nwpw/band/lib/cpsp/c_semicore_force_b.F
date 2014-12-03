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
      integer Gx(2),Gy(2),Gz(2)
      integer xtmp(2),dng(2)

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
     >        BA_push_get(mt_dbl, nfft3d,'Gx',Gx(2),Gx(1))
      value = value.and.
     >        BA_push_get(mt_dbl, nfft3d,'Gy',Gy(2),Gy(1))
      value = value.and.
     >        BA_push_get(mt_dbl, nfft3d,'Gz',Gz(2),Gz(1))

      value = value.and.
     >        BA_push_get(mt_dbl, npack0,'xtmp',xtmp(2),xtmp(1))

      value = value.and.
     >        BA_push_get(mt_dcpl, nfft3d,'dng',dng(2),dng(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)


      call dcopy(n2ft3d,0.0d0,0,dcpl_mb(vxcG(1)),1)
      call D3dB_rr_Sum(1,xcp(1),
     >                   xcp(1+(ispin-1)*n2ft3d),
     >                   dcpl_mb(vxcG(1)))
      call D3dB_rc_fft3f(1,dcpl_mb(vxcG(1)))
      call Pack_c_pack(0,dcpl_mb(vxcG(1)))

*     **** define Gx, Gy, and Gz in packed space ****
      call D3dB_t_Copy(1,dbl_mb(G_indx(1)),dbl_mb(Gx(1)))
      call D3dB_t_Copy(1,dbl_mb(G_indx(2)),dbl_mb(Gy(1)))
      call D3dB_t_Copy(1,dbl_mb(G_indx(3)),dbl_mb(Gz(1)))
      call Pack_t_pack(0,dbl_mb(Gx(1)))
      call Pack_t_pack(0,dbl_mb(Gy(1)))
      call Pack_t_pack(0,dbl_mb(Gz(1)))


      do ii=1,ion_nion()
         ia = ion_katm(ii)

         if (log_mb(semicore(1)+ia)) then

*          **** structure factor and local pseudopotential ****
           call strfac(ii,dcpl_mb(exi(1)))
           call Pack_c_pack(0,dcpl_mb(exi(1)))

*          **** put sqrt(core-density) at atom position ****
           call Pack_tc_iMul(0,
     >               dbl_mb(ncore(1)+(ia-1)*5*npack0),
     >              dcpl_mb(exi(1)),
     >              dcpl_mb(dng(1)))

*          **** put dng in real space and square it ****
           call Pack_c_unpack(0,dcpl_mb(dng(1)))
           call D3dB_cr_fft3b(1,dcpl_mb(dng(1)))
           call D3dB_rr_Sqr(1,dcpl_mb(dng(1)),dcpl_mb(dng(1)))

*          **** put dng back in complex space ****
           call D3dB_r_Zero_Ends(1,dcpl_mb(dng(1)))
           call D3dB_rc_fft3f(1,dcpl_mb(dng(1)))
           call Pack_c_pack(0,dcpl_mb(dng(1)))

           do k=1,npack0
              dbl_mb(xtmp(1)+k-1) =
     >         dimag(dcpl_mb(VxcG(1)+k-1))* dble(dcpl_mb(dng(1)+k-1))
     >        - dble(dcpl_mb(VxcG(1)+k-1))*dimag(dcpl_mb(dng(1)+k-1))
           end do

           call Pack_tt_dot(0,dbl_mb(Gx(1)),dbl_mb(xtmp(1)),sumx)
           call Pack_tt_dot(0,dbl_mb(Gy(1)),dbl_mb(xtmp(1)),sumy)
           call Pack_tt_dot(0,dbl_mb(Gz(1)),dbl_mb(xtmp(1)),sumz)

           fion(1,ii) = fion(1,ii) + sumx*scal2*scal2*dsqrt(scal1)
           fion(2,ii) = fion(2,ii) + sumy*scal2*scal2*dsqrt(scal1)
           fion(3,ii) = fion(3,ii) + sumz*scal2*scal2*dsqrt(scal1)
         end if

      end do

      value = BA_pop_stack(dng(2))
      value = BA_pop_stack(xtmp(2))
      value = BA_pop_stack(Gz(2))
      value = BA_pop_stack(Gy(2))
      value = BA_pop_stack(Gx(2))
      value = BA_pop_stack(vxcG(2))
      value = BA_pop_stack(exi(2))


      return
      end


c $Id$
