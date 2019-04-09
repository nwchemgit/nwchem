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
      integer if
c      
      call dfill(9,0d0,alpha,1)
      do if=1,n
         alpha(1,1)=alpha(1,1)+dpx(if)*dpx(if)/deltae(if)
         alpha(2,1)=alpha(2,1)+dpy(if)*dpx(if)/deltae(if)
         alpha(3,1)=alpha(3,1)+dpz(if)*dpx(if)/deltae(if)
         alpha(1,2)=alpha(1,2)+dpx(if)*dpy(if)/deltae(if)
         alpha(2,2)=alpha(2,2)+dpy(if)*dpy(if)/deltae(if)
         alpha(3,2)=alpha(3,2)+dpz(if)*dpy(if)/deltae(if)
         alpha(1,3)=alpha(1,3)+dpx(if)*dpz(if)/deltae(if)
         alpha(2,3)=alpha(2,3)+dpy(if)*dpz(if)/deltae(if)
         alpha(3,3)=alpha(3,3)+dpz(if)*dpz(if)/deltae(if)
      enddo
      if(ipol.eq.1) 
     c     call dscal(9,2d0,alpha,1)
      call util_poltensor(alpha,alpha_ev,alpha_iso,alpha_ani,' SOS')
        
      return
      end
