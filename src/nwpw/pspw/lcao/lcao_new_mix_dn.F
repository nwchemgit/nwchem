*
* $Id: lcao_new_mix_dn.F,v 1.3 2001-12-04 02:39:31 bylaska Exp $
*
      subroutine lcao_new_mix_dn(alpha,ispin,ne,n2ft3d,dn,
     >                           psimatrix,phi,rho_error)
      implicit none
      real*8  alpha
      integer ispin,ne(2)
      integer n2ft3d
      real*8 dn(n2ft3d,2)
      integer psimatrix(*)
      complex*16 phi(*)
      real*8     rho_error

#include "mafdecls.fh"
#include "global.fh"

*     **** local variables ****
      logical value
      integer n,k,nbasis,ms
      integer nx,ny,nz
      real*8  scal1,scal2,dv
      integer tmp(2),newdn(2),coef(2)


*     ***** external functions ****
      integer  aorbs_nbasis
      real*8   lattice_omega
      external aorbs_nbasis
      external lattice_omega


      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)
      scal1 = 1.0d0/dble(nx*ny*nz)
      scal2 = 1.0d0/lattice_omega()
      dv    = scal1*lattice_omega()

      nbasis = aorbs_nbasis()

      value = MA_push_get(mt_dbl,n2ft3d,'tmp',tmp(2),tmp(1))
      value = value.and.
     >        MA_push_get(mt_dbl,(2*n2ft3d),'newdn',newdn(2),newdn(1))
      value = value.and.
     >        MA_push_get(mt_dbl,(nbasis),'coef',coef(2),coef(1))
      if (.not. value) 
     >  call errquit('lcao_new_mix_dn:out of stack memory',0)


      call dcopy(2*n2ft3d,0.0d0,0,dbl_mb(newdn(1)),1)
      do ms=1,ispin
         do n=1,ne(ms)

*           **** calculate electronic wavefunction ****
            call ga_get(psimatrix(ms),1,nbasis,n,n,
     >                  dbl_mb(coef(1)),nbasis) 
            call lcao_make_psi(dbl_mb(coef(1)),dbl_mb(tmp(1)),phi)

*           **** put in real-space ***
            call Pack_c_unpack(1,dbl_mb(tmp(1)))
            call D3dB_cr_fft3b(1,dbl_mb(tmp(1)))

*           **** add to density ****
            do k=1,n2ft3d
              dbl_mb(newdn(1)+(ms-1)*n2ft3d+k-1)  
     >        = dbl_mb(newdn(1)+(ms-1)*n2ft3d+k-1) 
     >        + scal2*(dbl_mb(tmp(1)+k-1)**2)
            end do

         end do
      end do

*     **** calculate rho_error ****
      rho_error = 0.0d0
      do ms=1,ispin
         call D3dB_r_Zero_Ends(1,dbl_mb(newdn(1)+(ms-1)*n2ft3d))
         do k=1,n2ft3d
            rho_error = rho_error
     >                + dabs( dbl_mb(newdn(1)+(ms-1)*n2ft3d+k-1) 
     >                      - dn(k,ms))
         end do
      end do
      rho_error = rho_error*dv
      call D3dB_SumAll(rho_error)


*     **** Calculate New Density = using DIIS *****
      call lcao_diis_dn_generate(ispin,n2ft3d,
     >                           dbl_mb(newdn(1)),
     >                           rho_error,
     >                           dn)

c*     **** Calculate New Density = combination of MO DEN and old DN *
c      do ms=1,ispin
c         do k=1,n2ft3d
c               dn(k,ms) = alpha*dbl_mb(newdn(1)+(ms-1)*n2ft3d+k-1) 
c     >                  + (1.0d0-alpha)*dn(k,ms)
c         end do
c      end do


c*     **** normalize densities ***
c      do ms=1,ispin
c        call D3dB_r_dsum(1,dn(1,ms),sum)
c        sum = sum*dv
c        sum = dble(ne(ms))/sum
c        call D3dB_r_SMul(1,sum,dn(1,ms),dn(1,ms))
c      end do

      value = MA_pop_stack(coef(2))
      value = value.and.
     >        MA_pop_stack(newdn(2))
      value = value.and.
     >        MA_pop_stack(tmp(2))
      if (.not. value) 
     >  call errquit('lcao_init_dn:popping stack memory',0)
      return
      end


