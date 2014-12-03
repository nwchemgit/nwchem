*
* $Id$
*
      subroutine lcao_new_dn(ispin,ne,n2ft3d,dn,
     >                           psimatrix,phi)
      implicit none
#include "errquit.fh"
      integer ispin,ne(2)
      integer n2ft3d
      real*8 dn(n2ft3d,2)
      integer psimatrix(*)
      complex*16 phi(*)

#include "bafdecls.fh"
#include "global.fh"

*     **** local variables ****
      logical value
      integer n,k,nbasis,ms
      real*8  scal2
      integer tmp(2),coef(2)


*     ***** external functions ****
      integer  aorbs_nbasis
      real*8   lattice_omega
      external aorbs_nbasis
      external lattice_omega


      scal2 = 1.0d0/lattice_omega()

      nbasis = aorbs_nbasis()

      value = BA_push_get(mt_dbl,n2ft3d,'tmp',tmp(2),tmp(1))
      value = value.and.
     >        BA_push_get(mt_dbl,(nbasis),'coef',coef(2),coef(1))
      if (.not. value) 
     >  call errquit('lcao_new_dn:out of stack memory',0, MA_ERR)


      call dcopy(ispin*n2ft3d,0.0d0,0,dn,1)
      do ms=1,ispin
         do n=1,ne(ms)

*           **** calculate electronic wavefunction ****
            call ga_get(psimatrix(ms),1,nbasis,n,n,
     >                  dbl_mb(coef(1)),nbasis) 
            call ga_sync()
            call lcao_make_psi(dbl_mb(coef(1)),dbl_mb(tmp(1)),phi)

*           **** put in real-space ***
            call Pack_c_unpack(1,dbl_mb(tmp(1)))
            call D3dB_cr_pfft3b(1,1,dbl_mb(tmp(1)))

*           **** add to density ****
            do k=1,n2ft3d
              dn(k,ms)= dn(k,ms) + scal2*(dbl_mb(tmp(1)+k-1)**2)
            end do

         end do
      end do

      value = BA_pop_stack(coef(2))
	  value = value.and.BA_pop_stack(tmp(2))
      if (.not. value) 
     >  call errquit('lcao_new_dn:popping stack memory',0, MA_ERR)
      return
      end


