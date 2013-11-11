C> \ingroup selci
C> @{
      subroutine selci_mknsneed(ioconf,noconf,nintpo,nbitpi,
     $     norbs,nsneed)
*
* $Id$
*
#include "implicit.fh"
c
      dimension ioconf(nintpo,noconf), nsneed(3)
c
c     Determine the number of unpaired electrons that are necessary
c     in the coupling coefficients.  This is equal to (by mode)
c
c     1    for initial configuration generation (conf) we just
c              want the dimension of the spin spaces up to max(ns)
c     2    for a CI max(ns) + an additional 2 if there are also
c              doubly occupied orbitals in the config. with this ns
c     3    for a selection max(ns) over doubly substituted
c              configurations plus an additional 2 if there are also
c              doubly occupied orbitals in this config.
c
      dimension iocc(255)
c
      nsneed(1) = 0
      nsneed(2) = 0
      nsneed(3) = 0
      do iconf = 1,noconf
         call selci_upkcon(norbs, iocc, ioconf(1,iconf), nintpo, nbitpi)
         ns = 0
	 nd = 0
         do j = 1,norbs
            if (iocc(j).eq.1) ns = ns + 1
	    if (iocc(j).eq.3) nd = nd + 1
         enddo
         nsneed(1) = max(nsneed(1),ns) 
c
         if (nd .gt. 0 ) then
            nd = nd - 1
            ns = ns + 2
         endif
         nsneed(2) = max(nsneed(2),ns) 
c
         if (nd .gt. 0 ) then
            nd = nd - 1
            ns = ns + 2
         endif
         if (nd .gt. 0 ) then
            nd = nd - 1
            ns = ns + 2
         endif
	 nsneed(3) = max(nsneed(3),ns)
      enddo
c
      end
C> @}
