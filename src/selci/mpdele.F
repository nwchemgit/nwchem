C> \ingroup selci
C> @{
      subroutine selci_mpdele(numf, nroot, deltae_mp, vc,
     $     roots_mp, nfmax, dnorm_mp,iocc)
*
* $Id$
*
#include "implicit.fh"
      integer iocc(*)
c
c     deltae(root) =  - <root|H|*><*|H|root> / (D(Eroot))
c     D(Eroot) = (E0(root) - E0(occupation))
c     dnorm = ( <*|H|root>/(D(Eroot)) ) **2
c
      dimension deltae_mp(nroot), vc(nfmax, nroot), roots_mp(nroot),
     $     dnorm_mp(nroot)
c
      call dfill(nroot, 0.0d0, deltae_mp, 1)
      call dfill(nroot, 0.0d0, dnorm_mp, 1)
c
      e0j = selci_e0mp2(iocc)
      do 10 ir = 1,nroot
         denom = roots_mp(ir) - e0j
         do 20 i = 1,numf
            coeff = vc(i,ir)/denom
c:old            deltae_mp(ir) = deltae_mp(ir) - vc(i,ir)*coeff
            deltae_mp(ir) = deltae_mp(ir) + vc(i,ir)*coeff
            dnorm_mp(ir)  = dnorm_mp(ir) + coeff*coeff
 20      continue
 10   continue
c
      end
C> @}
