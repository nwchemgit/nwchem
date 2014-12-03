*
* $Id$
*
      subroutine lcao_generate_vlmatrix(vlmatrix,phi1,phi2)
      implicit none
#include "errquit.fh"
      integer vlmatrix
      complex*16 phi1(*)
      complex*16 phi2(*)

#include "bafdecls.fh"
#include "global.fh"

*     **** local variables ****
      logical move,value
      integer i,j,nbasis,npack1,n2ft3d,nfft3d
      integer nx,ny,nz
      real*8  sum,fion(3),scal1,scal2
      complex*16 junk(2)
      integer vl(2),vlr_l(2),rgrid(2),vpsi

*     ***** external functions ****
      integer  aorbs_nbasis,control_version
      real*8   lattice_omega
      external aorbs_nbasis,control_version
      external lattice_omega

      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)
      scal1 = 1.0d0/dble(nx*ny*nz)
      scal2 = 1.0d0/lattice_omega()

      call Pack_npack(1,npack1)
      call D3dB_nfft3d(1,nfft3d)
      n2ft3d = 2*nfft3d
      value = BA_push_get(mt_dcpl,nfft3d,'vall',vl(2),vl(1))
      value = value.and.
     >        BA_push_get(mt_dbl,n2ft3d,'vlr_l',vlr_l(2),vlr_l(1))
      vpsi = vlr_l(1)
      if (control_version().eq.4)
     >  value = value.and.
     >        BA_push_get(mt_dbl,(3*n2ft3d),'rgrid',rgrid(2),rgrid(1))
      if (.not. value) 
     > call errquit('lcao_generate_vlmatrix: out of stack memory',0,
     &       MA_ERR)


*     *****************************************
*     **** generate local pseudopotential  ****
*     *****************************************
      move =.false.
      call v_local(dcpl_mb(vl(1)),
     >             move,
     >             junk,
     >             fion)
*     *** long-range psp for charge systems ***
      if (control_version().eq.4) then
        call lattice_r_grid(dbl_mb(rgrid(1)))
        call v_lr_local(dbl_mb(rgrid(1)),
     >                  dbl_mb(vlr_l(1)))
        value = BA_pop_stack(rgrid(2))
      end if

*     **********************************************************
*     **** add up k-space potentials, vall = scal2*vsr_l    ****
*     **********************************************************
c      call Pack_c_SMul(0,scal2,dcpl_mb(vl(1)),dcpl_mb(vl(1)))
      call Pack_c_SMul1(0,scal2,dcpl_mb(vl(1)))


*     **********************************************
*     **** fourier transform k-space potentials ****
*     **********************************************
      call Pack_c_unpack(0,dcpl_mb(vl(1)))
      call D3dB_cr_pfft3b(1,0,dcpl_mb(vl(1)))

*     **********************************
*     **** add vlr_l                ****
*     **********************************
      if (control_version().eq.4) then
c        call D3dB_rr_Sum(1,dcpl_mb(vl(1)),
c     >                   dbl_mb(vlr_l(1)),
c     >                   dcpl_mb(vl(1)))
        call D3dB_rr_Sum2(1,dbl_mb(vlr_l(1)),dcpl_mb(vl(1)))
      end if



      nbasis = aorbs_nbasis()
      do j=1,nbasis

*        **** get phi2 ****
         call aorbs_aorb(j,phi2)
   
*       **** get psi_r ****
         call Pack_c_Copy(1,phi2,dbl_mb(vpsi))
         call Pack_c_unpack(1,dbl_mb(vpsi))
         call D3dB_cr_pfft3b(1,1,dbl_mb(vpsi))
         call D3dB_r_Zero_Ends(1,dbl_mb(vpsi))

c         call D3dB_rr_Mul(1,dcpl_mb(vl(1)),
c     >                      dbl_mb(vpsi),
c     >                      dbl_mb(vpsi))
         call D3dB_rr_Mul2(1,dcpl_mb(vl(1)),dbl_mb(vpsi))

*        *******************************
*        **** get vpsi_r in k-space ****
*        *******************************
         call D3dB_rc_pfft3f(1,1,dbl_mb(vpsi))
         call Pack_c_pack(1,dbl_mb(vpsi))
         call Pack_c_Copy(1,dbl_mb(vpsi),phi2)

         do i=1,j

*           **** get phi1 ****
            call aorbs_aorb(i,phi1)
            call Pack_cc_dot(1,phi1,phi2,sum)
            sum = sum*scal1

            call ga_fill_patch(vlmatrix,i,i,j,j,sum)
            call ga_fill_patch(vlmatrix,j,j,i,i,sum)
         end do
      end do

      value = BA_pop_stack(vlr_l(2))
      value = BA_pop_stack(vl(2))
     
      return
      end


