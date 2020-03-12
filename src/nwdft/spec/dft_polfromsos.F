      subroutine dft_polfromsos(alpha,n,dpx,dpy,dpz,deltae,ipol)
      implicit none
c     alpha=2*\sum_{i,f} dipole_{if}^2/(E_f-E_i)
      double precision alpha(3,3)           ! polarizability [out]
      integer n                             ! no. of states [in]
      double precision dpx(*),dpy(*),dpz(*) ! dipoles [in]
      double precision deltae(*) ! DeltaEigenv [in]
      integer ipol ! [in]
c
      double precision alpha_ev(3)          ! polarizability evec  [out]
      double precision alpha_iso            ! polarizability isotr [out]
      double precision alpha_ani            ! polarizability aniso [out]
      
#include "util_params.fh"      
#include "mafdecls.fh"      
c
      integer i_f
c      
      call dfill(9,0d0,alpha,1)
      do i_f=1,n
         alpha(1,1)=alpha(1,1)+dpx(i_f)*dpx(i_f)/deltae(i_f)
         alpha(2,1)=alpha(2,1)+dpy(i_f)*dpx(i_f)/deltae(i_f)
         alpha(3,1)=alpha(3,1)+dpz(i_f)*dpx(i_f)/deltae(i_f)
         alpha(1,2)=alpha(1,2)+dpx(i_f)*dpy(i_f)/deltae(i_f)
         alpha(2,2)=alpha(2,2)+dpy(i_f)*dpy(i_f)/deltae(i_f)
         alpha(3,2)=alpha(3,2)+dpz(i_f)*dpy(i_f)/deltae(i_f)
         alpha(1,3)=alpha(1,3)+dpx(i_f)*dpz(i_f)/deltae(i_f)
         alpha(2,3)=alpha(2,3)+dpy(i_f)*dpz(i_f)/deltae(i_f)
         alpha(3,3)=alpha(3,3)+dpz(i_f)*dpz(i_f)/deltae(i_f)
      enddo
      call dscal(9,4d0/ipol,alpha,1)
      call util_poltensor(alpha,alpha_ev,alpha_iso,alpha_ani,' SOS')
        
      return
      end
