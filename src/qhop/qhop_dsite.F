      subroutine qhop_dsite(ga,rsa,rind)
c
c $Id: qhop_dsite.F,v 1.1 2003-07-17 19:36:51 d3j191 Exp $
c
      implicit none
c
#include "qhop_common.fh"
c
      integer ga,rind,rsa
c
      rind=arat(ptarat(ptseq(rsa,1))+ga-
     + ptseq(rsa,2),1)
c      
      return
      end
