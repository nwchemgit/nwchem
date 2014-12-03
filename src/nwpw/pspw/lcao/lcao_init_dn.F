*
* $Id$
*
      subroutine lcao_init_dn(ispin,ne,n2ft3d,dn,phi)
      implicit none
#include "errquit.fh"
      integer ispin,ne(2)
      integer n2ft3d
      real*8 dn(n2ft3d,2)
      complex*16 phi(*)

#include "bafdecls.fh"
#include "global.fh"

*     **** local variables ****
      logical value
      integer i,k,nbasis,ms
      integer nx,ny,nz
      real*8  sum,scal,scal1,scal2,dv
      real*8 dnscal(2)
      integer tmp(2)


*     ***** external functions ****
      integer  aorbs_nbasis
      real*8   aorbs_weight,lattice_omega,util_random
      external aorbs_nbasis
      external aorbs_weight,lattice_omega,util_random

      dnscal(1) = dble(ne(1))/dble(ne(1)+ne(2))
      dnscal(2) = dble(ne(2))/dble(ne(1)+ne(2))

      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)
      scal1 = 1.0d0/dble(nx*ny*nz)
      scal2 = 1.0d0/lattice_omega()
      dv    = scal1*lattice_omega()

      value = BA_push_get(mt_dbl,n2ft3d,'tmp',tmp(2),tmp(1))
      if (.not. value) 
     >  call errquit('lcao_init_dn:out of stack memory',0, MA_ERR)


      call dcopy(ispin*n2ft3d,0.0d0,0,dn,1)
      nbasis = aorbs_nbasis()
      do i=1,nbasis

*        **** get phi1 ****
         call aorbs_aorb(i,phi)
         call Pack_c_Copy(1,phi,dbl_mb(tmp(1)))
         call Pack_c_unpack(1,dbl_mb(tmp(1)))
         call D3dB_cr_pfft3b(1,1,dbl_mb(tmp(1)))

        do ms=1,ispin
           do k=1,n2ft3d
              scal = aorbs_weight(i)*scal2*dnscal(ms)
              dn(i,ms) = dn(i,ms) + scal*(dbl_mb(tmp(1)+k-1)**2)
           end do
        end do

      end do


*     **** randomize dn(*,2) if ispin=2 and ne(1)==ne(2) ****
      if ((ispin.eq.2).and.(ne(1).eq.ne(2))) then
           scal = util_random(9)
           do k=1,n2ft3d
              scal = (0.5d0-util_random(0))/dsqrt(dble(n2ft3d))
              dn(i,2) = dn(i,2) + scal
           end do
      end if


*     **** normalize densities ***
      do ms=1,ispin
        call D3dB_r_Zero_Ends(1,dn(1,ms))
        call D3dB_r_dsum(1,dn(1,ms),sum)
        sum = sum*dv
        sum = dble(ne(ms))/sum
c        call D3dB_r_SMul(1,sum,dn(1,ms),dn(1,ms))
        call D3dB_r_SMul1(1,sum,dn(1,ms))
      end do

      value = BA_pop_stack(tmp(2))
      if (.not. value) 
     >  call errquit('lcao_init_dn:popping stack memory',0, MA_ERR)
      return
      end


