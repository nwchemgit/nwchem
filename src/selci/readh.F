C> \ingroup selci
C> @{
      subroutine selci_readh(iflham,hh,nnci)
*
* $Id$
*
#include "implicit.fh"
      dimension hh(nnci)
c
      read (iflham) hh
c
      end
C> @}
