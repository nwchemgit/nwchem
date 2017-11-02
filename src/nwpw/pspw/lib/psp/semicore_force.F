*
* $Id$
*

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
      real*8  scal1,scal2
      integer exi(2),vxcG(2)
      integer Gx(2),Gy(2),Gz(2)
      integer dng(2)
      integer dngx(2),dngy(2),dngz(2)
      real*8  sumx,sumy,sumz
      common /semicore_xc_F_sums/ sumx,sumy,sumz

*     **** external functions ****
      integer  ion_nion,ion_katm,G_indx
      real*8   lattice_omega
      external ion_nion,ion_katm,G_indx
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
     >        BA_push_get(mt_dcpl, nfft3d,'dng',dng(2),dng(1))
      value = value.and.
     >        BA_push_get(mt_dcpl, nfft3d,'dngx',dngx(2),dngx(1))
      value = value.and.
     >        BA_push_get(mt_dcpl, nfft3d,'dngy',dngy(2),dngy(1))
      value = value.and.
     >        BA_push_get(mt_dcpl, nfft3d,'dngz',dngz(2),dngz(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)


      call D3dB_rr_Sum(1,xcp(1),
     >                   xcp(1+(ispin-1)*n2ft3d),
     >                   dcpl_mb(vxcG(1)))
c     call D3dB_r_SMul(1,0.5d0,dcpl_mb(vxcG(1)),dcpl_mb(vxcG(1)))
      call D3dB_r_Zero_Ends(1,dcpl_mb(vxcG(1)))


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
           call Pack_tc_Mul(0,
     >               dbl_mb(ncore(1)+(ia-1)*5*npack0),
     >              dcpl_mb(exi(1)),
     >              dcpl_mb(dng(1)))

           call Pack_tc_iMul(0,dbl_mb(Gx(1)),dcpl_mb(dng(1)),
     >                                       dcpl_mb(dngx(1)))
           call Pack_tc_iMul(0,dbl_mb(Gy(1)),dcpl_mb(dng(1)),
     >                                       dcpl_mb(dngy(1)))
           call Pack_tc_iMul(0,dbl_mb(Gz(1)),dcpl_mb(dng(1)),
     >                                       dcpl_mb(dngz(1)))
c          call Pack_c_SMul(0,(-1.0d0),dcpl_mb(dngx(1)),
c    >                                 dcpl_mb(dngx(1)))
c          call Pack_c_SMul(0,(-1.0d0),dcpl_mb(dngy(1)),
c    >                                 dcpl_mb(dngy(1)))
c          call Pack_c_SMul(0,(-1.0d0),dcpl_mb(dngz(1)),
c    >                                 dcpl_mb(dngz(1)))

*          **** put dng,dngx,dngy,dngz in real space ****
           call Pack_c_unpack(0,dcpl_mb(dng(1)))
           call Pack_c_unpack(0,dcpl_mb(dngx(1)))
           call Pack_c_unpack(0,dcpl_mb(dngy(1)))
           call Pack_c_unpack(0,dcpl_mb(dngz(1)))

           call D3dB_cr_fft3b(1,dcpl_mb(dng(1)))
           call D3dB_cr_fft3b(1,dcpl_mb(dngx(1)))
           call D3dB_cr_fft3b(1,dcpl_mb(dngy(1)))
           call D3dB_cr_fft3b(1,dcpl_mb(dngz(1)))
c           call D3dB_rr_Mul(1,dcpl_mb(dng(1)),
c     >                        dcpl_mb(dngx(1)),
c     >                        dcpl_mb(dngx(1)))
c           call D3dB_rr_Mul(1,dcpl_mb(dng(1)),
c     >                        dcpl_mb(dngy(1)),
c     >                        dcpl_mb(dngy(1)))
c           call D3dB_rr_Mul(1,dcpl_mb(dng(1)),
c     >                        dcpl_mb(dngz(1)),
c     >                        dcpl_mb(dngz(1)))
           call D3dB_rr_Mul2(1,dcpl_mb(dng(1)),dcpl_mb(dngx(1)))
           call D3dB_rr_Mul2(1,dcpl_mb(dng(1)),dcpl_mb(dngy(1)))
           call D3dB_rr_Mul2(1,dcpl_mb(dng(1)),dcpl_mb(dngz(1)))

           call D3dB_rr_dot(1,dcpl_mb(dngx(1)),dcpl_mb(vxcG(1)),sumx)
           call D3dB_rr_dot(1,dcpl_mb(dngy(1)),dcpl_mb(vxcG(1)),sumy)
           call D3dB_rr_dot(1,dcpl_mb(dngz(1)),dcpl_mb(vxcG(1)),sumz)

!$OMP MASTER
           fion(1,ii) = fion(1,ii) + sumx*scal2*scal1
           fion(2,ii) = fion(2,ii) + sumy*scal2*scal1
           fion(3,ii) = fion(3,ii) + sumz*scal2*scal1
!$OMP END MASTER

c          write(*,*) "semicore-force:",sumx*scal2*scal1,
c    >                                  sumy*scal2*scal1,
c    >                                  sumz*scal2*scal1
         end if

      end do
!$OMP BARRIER

      value = BA_pop_stack(dngz(2))
      value = value.and.BA_pop_stack(dngy(2))
      value = value.and.BA_pop_stack(dngx(2))
      value = value.and.BA_pop_stack(dng(2))
      value = value.and.BA_pop_stack(Gz(2))
      value = value.and.BA_pop_stack(Gy(2))
      value = value.and.BA_pop_stack(Gx(2))
      value = value.and.BA_pop_stack(vxcG(2))
      value = value.and.BA_pop_stack(exi(2))
      if (.not. value) 
     > call errquit('semicore_force:error popping stack',0, MA_ERR)

      return
      end


