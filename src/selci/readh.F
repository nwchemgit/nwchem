      subroutine selci_readh(iflham,hh,nnci)
*
* $Id: readh.F,v 1.2 1997-10-31 23:42:24 d3e129 Exp $
*
#include "implicit.fh"
      dimension hh(nnci)
c
      read (iflham) hh
c
      end
