*
* $Id: lcao_generate_vhxcmatrix.F,v 1.1 2001-12-01 23:13:26 bylaska Exp $
*
      subroutine lcao_generate_vhxcmatrix(ispin,n2ft3d,dn,
     >                                    vhxcmatrix,phi1,phi2)
      implicit none
      integer ispin,n2ft3d
      real*8  dn(n2ft3d,2)
      integer vhxcmatrix(2)
      complex*16 phi1(*)
      complex*16 phi2(*)

#include "mafdecls.fh"
#include "global.fh"

*     **** local variables ****
      logical value
      integer ms,i,j,nbasis,npack1,nfft3d
      integer nx,ny,nz
      real*8  sum,scal1,scal2,dv
      integer xcp(2),vc(2),dng(2),vpsi,vall

*     ***** external functions ****
      integer  aorbs_nbasis,control_version,control_gga
      real*8   lattice_omega
      external aorbs_nbasis,control_version,control_gga
      external lattice_omega

      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)
      scal1 = 1.0d0/dble(nx*ny*nz)
      scal2 = 1.0d0/lattice_omega()
      dv    = lattice_omega()*scal1

      call Pack_npack(1,npack1)
      call D3dB_nfft3d(1,nfft3d)
      value = MA_push_get(mt_dbl,(2*n2ft3d),'xcp',xcp(2),xcp(1))
      value = value.and.
     >        MA_push_get(mt_dbl,(2*n2ft3d),'vc',vc(2),vc(1))
      value = value.and.
     >        MA_push_get(mt_dbl,(2*n2ft3d),'dng',dng(2),dng(1))
      vpsi = dng(1)
      vall = dng(1)+n2ft3d
      if (.not. value) 
     > call errquit('lcao_geneate_vhxcmatrix: out of stack memory',0)


*     *************************************************
*     **** generate exchange-correlation potential ****
*     *************************************************
      if (control_gga().eq.0) then
        call vxc(n2ft3d,ispin,dn,
     >                    dbl_mb(xcp(1)),
     >                    dbl_mb(vc(1)),
     >                    dbl_mb(dng(1)))
      else
        call v_bwexc(n2ft3d,ispin,dn,
     >                    dbl_mb(xcp(1)),
     >                    dbl_mb(vc(1)))
      end if


*     ************************************
*     **** generate coulomb potential ****
*     ************************************
      call D3dB_rr_Sum(1,dn(1,1),dn(1,ispin),dbl_mb(vall))

      if (control_version().eq.3) then
         call D3dB_r_SMul(1,scal1,dbl_mb(vall),dbl_mb(dng(1)))
         call D3dB_rc_fft3f(1,dbl_mb(dng(1)))
         call Pack_c_pack(0,dbl_mb(dng(1)))
         call coulomb_v(dbl_mb(dng(1)),dbl_mb(vc(1)))
         call Pack_c_unpack(0,dbl_mb(vc(1)))
         call D3dB_cr_fft3b(1,dbl_mb(vc(1)))
      else
         call coulomb2_v(dbl_mb(vall),dbl_mb(vc(1)))
      end if


*     **********************************************
*     **** fourier transform k-space potentials ****
*     **********************************************

      nbasis = aorbs_nbasis()
      do j=1,nbasis

*        **** get phi2 ****
         call aorbs_aorb(j,phi2)
   
*       **** get psi_r ****
         call Pack_c_Copy(1,phi2,dbl_mb(vpsi))
         call Pack_c_unpack(1,dbl_mb(vpsi))
         call D3dB_cr_fft3b(1,dbl_mb(vpsi))
         call D3dB_r_Zero_Ends(1,dbl_mb(vpsi))

         do ms=1,ispin
            call D3dB_rr_Sum(1,dbl_mb(vc(1)),
     >                         dbl_mb(xcp(1)+(ms-1)*n2ft3d),
     >                         dbl_mb(vall))

            call D3dB_rr_Mul(1,dbl_mb(vall),
     >                         dbl_mb(vpsi),
     >                         dbl_mb(vpsi))

*           *******************************
*           **** get vpsi_r in k-space ****
*           *******************************
            call D3dB_rc_fft3f(1,dbl_mb(vpsi))
            call Pack_c_pack(1,dbl_mb(vpsi))
            call Pack_c_Copy(1,dbl_mb(vpsi),phi2)

            do i=1,j

*             **** get phi1 ****
              call aorbs_aorb(i,phi1)
              call Pack_cc_dot(1,phi1,phi2,sum)
              sum = sum*scal1

              call ga_fill_patch(vhxcmatrix(ms),i,i,j,j,sum)
              call ga_fill_patch(vhxcmatrix(ms),j,j,i,i,sum)
            end do
         end do
      end do

      value = MA_pop_stack(dng(2))
      value = MA_pop_stack(vc(2))
      value = MA_pop_stack(xcp(2))
     
      return
      end


