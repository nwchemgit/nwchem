*
* $Id$
*
      subroutine lcao_new_mix_dn(alpha,ispin,ne,n2ft3d,dn,
     >                           psimatrix,phi,rho_error,reset_diis)
      implicit none
#include "errquit.fh"
      real*8  alpha
      integer ispin,ne(2)
      integer n2ft3d
      real*8 dn(n2ft3d,2)
      integer psimatrix(*)
      complex*16 phi(*)
      real*8     rho_error
      logical    reset_diis

#include "bafdecls.fh"
#include "global.fh"

*     **** local variables ****
      logical value
      integer n,k,nbasis,ms
      integer nx,ny,nz
      real*8  scal1,scal2,dv,sum
      integer tmp(2),newdn(2),coef(2),dnerr(2)
      real*8  rho_error_old


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

      value = BA_push_get(mt_dbl,n2ft3d,'tmp',tmp(2),tmp(1))
      value = value.and.
     >        BA_push_get(mt_dbl,(2*n2ft3d),'newdn',newdn(2),newdn(1))
      value = value.and.
     >        BA_push_get(mt_dbl,(2*n2ft3d),'dnerr',dnerr(2),dnerr(1))
      value = value.and.
     >        BA_push_get(mt_dbl,(nbasis),'coef',coef(2),coef(1))
      if (.not. value) 
     >  call errquit('lcao_new_mix_dn:out of stack memory',0, MA_ERR)


      call dcopy(2*n2ft3d,0.0d0,0,dbl_mb(newdn(1)),1)
      do ms=1,ispin
         do n=1,ne(ms)

*           **** calculate electronic wavefunction ****
            call ga_get(psimatrix(ms),1,nbasis,n,n,
     >                  dbl_mb(coef(1)),nbasis) 
            call lcao_make_psi(dbl_mb(coef(1)),dbl_mb(tmp(1)),phi)

*           **** put in real-space ***
            call Pack_c_unpack(1,dbl_mb(tmp(1)))
            call D3dB_cr_pfft3b(1,1,dbl_mb(tmp(1)))

*           **** add to density ****
            do k=1,n2ft3d
              dbl_mb(newdn(1)+(ms-1)*n2ft3d+k-1)  
     >        = dbl_mb(newdn(1)+(ms-1)*n2ft3d+k-1) 
     >        + scal2*(dbl_mb(tmp(1)+k-1)**2)
            end do

         end do
      end do


*     **** calculate rho_error ****
      rho_error_old = rho_error
      rho_error     = 0.0d0
      do ms=1,ispin
         call D3dB_r_Zero_Ends(1,dbl_mb(newdn(1)+(ms-1)*n2ft3d))
         do k=1,n2ft3d
            call D3dB_rr_Minus(1,
     >                         dbl_mb(newdn(1)+(ms-1)*n2ft3d),
     >                         dn(1,ms),
     >                         dbl_mb(dnerr(1)+(ms-1)*n2ft3d))

         end do
         call D3dB_rr_dot(1,dbl_mb(dnerr(1)+(ms-1)*n2ft3d),
     >                      dbl_mb(dnerr(1)+(ms-1)*n2ft3d),
     >                      sum)
         rho_error = rho_error + sum
      end do
      rho_error = rho_error*dv
      if (rho_error.ge.rho_error_old) alpha=0.2d0*alpha

c     reset_diis =    (rho_error.ge.rho_error_old)
c    >           .or.(rho_error.gt.1.0d0)
      reset_diis =    (rho_error.ge.rho_error_old)
     >           .or.(rho_error.gt.0.10d0)

      if (.not. reset_diis) then
*       **** Calculate New Density = using DIIS *****
        call lcao_diis_dn_generate(alpha,ispin,n2ft3d,
     >                           dbl_mb(newdn(1)),
     >                           dn,reset_diis)
      end if


*     **** soften the new density by the old density ****
      if (reset_diis) then
      do ms=1,ispin
         do k=1,n2ft3d
            dn(k,ms) = alpha*dbl_mb(newdn(1)+(ms-1)*n2ft3d+k-1) 
     >                  + (1.0d0-alpha)*dn(k,ms)
         end do
      end do
      end if


c*     **** normalize densities ***
c      do ms=1,ispin
c        call D3dB_r_dsum(1,dn(1,ms),sum)
c        sum = sum*dv
c        sum = dble(ne(ms))/sum
c        call D3dB_r_SMul(1,sum,dn(1,ms),dn(1,ms))
c      end do

      value = BA_pop_stack(coef(2))
      value = value.and.
     >        BA_pop_stack(dnerr(2))
      value = value.and.
     >        BA_pop_stack(newdn(2))
      value = value.and.
     >        BA_pop_stack(tmp(2))
      if (.not. value) 
     >  call errquit('lcao_init_dn:popping stack memory',0, MA_ERR)
      return
      end


