*
* $Id$
*
      subroutine lcao_generate_vhxcmatrix(ispin,n2ft3d,dn,
     >                                    vhxcmatrix,phi1,phi2,
     >                                    ehartr,exc,pxc)
      implicit none
      integer ispin,n2ft3d
      real*8  dn(n2ft3d,2)
      integer vhxcmatrix(2)
      complex*16 phi1(*)
      complex*16 phi2(*)
      real*8     ehartr,exc,pxc

#include "bafdecls.fh"
#include "global.fh"
#include "nwpwxc.fh"
#include "util.fh"
#include "errquit.fh"

*     **** local variables ****
      logical value, use_lda
      integer ms,i,j,nbasis,npack1,nfft3d,gga
      integer nx,ny,nz
      real*8  sum,scal1,scal2,dv,exc2,pxc2
      integer xcp(2),tmp1(2),tmp2(2),dnall(2),vpsi,vall,xvc,rho

*     ***** external functions ****
      logical  psp_semicore
      integer  aorbs_nbasis,control_version,control_gga
      real*8   lattice_omega,coulomb_e
      external psp_semicore
      external aorbs_nbasis,control_version,control_gga
      external lattice_omega,coulomb_e

      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)
      scal1 = 1.0d0/dble(nx*ny*nz)
      scal2 = 1.0d0/lattice_omega()
      dv    = lattice_omega()*scal1

      call Pack_npack(1,npack1)
      call D3dB_nfft3d(1,nfft3d)
      value = BA_push_get(mt_dbl,(2*n2ft3d),'xcp',xcp(2),xcp(1))
      value = value.and.
     >        BA_push_get(mt_dbl,(2*n2ft3d),'tmp1',tmp1(2),tmp1(1))
      value = value.and.
     >        BA_push_get(mt_dbl,(2*n2ft3d),'tmp2',tmp2(2),tmp2(1))
      value = value.and.
     >        BA_push_get(mt_dbl,(2*n2ft3d),'dnall',dnall(2),dnall(1))
      call dcopy(2*n2ft3d,0.0d0,0,dbl_mb(xcp(1)),1)
      call dcopy(2*n2ft3d,0.0d0,0,dbl_mb(tmp1(1)),1)
      call dcopy(2*n2ft3d,0.0d0,0,dbl_mb(tmp2(1)),1)
      call dcopy(2*n2ft3d,0.0d0,0,dbl_mb(dnall(1)),1)
      xvc  = tmp1(1)
      rho  = tmp1(1)+n2ft3d
      vpsi = tmp2(1)
      vall = tmp2(1)+n2ft3d
      if (.not. value) 
     > call errquit('lcao_geneate_vhxcmatrix: out of stack memory',0,
     &       MA_ERR)


*     *************************************************
*     **** generate exchange-correlation potential ****
*     *************************************************
      if (psp_semicore(0)) then
           call semicore_density(dbl_mb(tmp1(1)))
c           call D3dB_r_SMul(1,0.5d0,dbl_mb(tmp1(1)),dbl_mb(tmp1(1)))
           call D3dB_r_SMul1(1,0.5d0,dbl_mb(tmp1(1)))
      else
           call dcopy(n2ft3d,0.0d0,0,dbl_mb(tmp1(1)),1)
      end if
      do ms=1,ispin
          call D3dB_rr_Sum(1,dn(1,ms),
     >                     dbl_mb(tmp1(1)),
     >                     dbl_mb(dnall(1) +(ms-1)*n2ft3d))
      end do


      gga = control_gga()
      use_lda = (.not.nwpwxc_is_on().and.gga.eq.0).or.
     &          (nwpwxc_is_on().and.nwpwxc_is_lda())
      
      if (use_lda) then
        call vxc(n2ft3d,ispin,dbl_mb(dnall(1)),
     >                    dbl_mb(xcp(1)),
     >                    dbl_mb(tmp1(1)),
     >                    dbl_mb(tmp2(1)))
      else 
        call v_bwexc(gga,n2ft3d,ispin,dbl_mb(dnall(1)),
     >                    1.0d0,1.0d0,
     >                    dbl_mb(xcp(1)),
     >                    dbl_mb(tmp1(1)))
      end if

*     **** get exchange-correlation energy ****
      call D3dB_rr_dot(1,dbl_mb(dnall(1)),dbl_mb(tmp1(1)),exc)
      call D3dB_rr_dot(1,dn(1,1),dbl_mb(xcp(1)),pxc)

      if (ispin.eq.1) then
         exc= exc + exc
         pxc= pxc + pxc
      else
         call D3dB_rr_dot(1,dbl_mb(dnall(1)+n2ft3d),
     >                      dbl_mb(tmp1(1)),exc2)
         call D3dB_rr_dot(1,dn(1,2),dbl_mb(xcp(1)+n2ft3d),pxc2)
         exc= exc + exc2
         pxc= pxc + pxc2
      end if
      exc = exc*dv
      pxc = pxc*dv

      value = BA_pop_stack(dnall(2))
      if (.not. value) 
     > call errquit('lcao_geneate_vhxcmatrix: popping stack memory',1,
     &       MA_ERR)

*     ************************************
*     **** generate coulomb potential ****
*     ************************************
      call D3dB_rr_Sum(1,dn(1,1),dn(1,ispin),dbl_mb(rho))

      if (control_version().eq.3) then
c         call D3dB_r_SMul(1,scal1,dbl_mb(rho),dbl_mb(rho))
         call D3dB_r_SMul1(1,scal1,dbl_mb(rho))
         call D3dB_rc_pfft3f(1,0,dbl_mb(rho))
         call Pack_c_pack(0,dbl_mb(rho))
         ehartr = coulomb_e(dbl_mb(rho))
         call coulomb_v(dbl_mb(rho),dbl_mb(xvc))
         call Pack_c_unpack(0,dbl_mb(xvc))
         call D3dB_cr_pfft3b(1,0,dbl_mb(xvc))
      else
         call coulomb2_v(dbl_mb(rho),dbl_mb(xvc))
         call D3dB_rr_dot(1,dbl_mb(rho),dbl_mb(xvc),ehartr)
         ehartr = 0.5d0*ehartr*dv
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
         call D3dB_cr_pfft3b(1,1,dbl_mb(vpsi))
         call D3dB_r_Zero_Ends(1,dbl_mb(vpsi))

         do ms=1,ispin
            call D3dB_rr_Sum(1,dbl_mb(xvc),
     >                         dbl_mb(xcp(1)+(ms-1)*n2ft3d),
     >                         dbl_mb(vall))

c            call D3dB_rr_Mul(1,dbl_mb(vall),
c     >                         dbl_mb(vpsi),
c     >                         dbl_mb(vall))
            call D3dB_rr_Mul2(1,dbl_mb(vpsi),dbl_mb(vall))

*           *******************************
*           **** get vpsi_r in k-space ****
*           *******************************
            call D3dB_rc_pfft3f(1,1,dbl_mb(vall))
            call Pack_c_pack(1,dbl_mb(vall))
            call Pack_c_Copy(1,dbl_mb(vall),phi2)

            do i=1,j

*             **** get phi1 ****
              call aorbs_aorb(i,phi1)
              call Pack_cc_dot(1,phi1,phi2,sum)
              sum = sum*scal1
c              write(*,*) "vhxc:",i,j,sum

              call ga_fill_patch(vhxcmatrix(ms),i,i,j,j,sum)
              call ga_fill_patch(vhxcmatrix(ms),j,j,i,i,sum)
            end do
         end do
      end do

      value =           BA_pop_stack(tmp2(2))
      value = value.and.BA_pop_stack(tmp1(2))
      value = value.and.BA_pop_stack(xcp(2))
      if (.not. value) 
     > call errquit('lcao_geneate_vhxcmatrix: popping stack memory',2,
     &       MA_ERR)
     
      return
      end


