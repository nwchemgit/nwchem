      subroutine selci_readh(iflham,hh,nnci)
#include "implicit.fh"
      dimension hh(nnci)
c
      read (iflham) hh
c
      end
