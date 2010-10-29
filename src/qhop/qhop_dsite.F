      subroutine qhop_dsite(ga,rsa,rind)
c
c $Id$
c
      implicit none
c
#include "qhop_common.fh"
#include "mafdecls.fh"
c
      integer ga,rind,rsa
c
      call qhop_dsite2(ga,rsa,rind,
     + int_mb(i_ptseq),int_mb(i_ptarat),int_mb(i_iarat),
     + dbl_mb(i_racs),dbl_mb(i_racs+mxar),int_mb(i_ptpar),
     + dbl_mb(i_par),dbl_mb(i_par+7*maxpar),dbl_mb(i_par+14*maxpar),
     + dbl_mb(i_par+17*maxpar),dbl_mb(i_par+23*maxpar),
     + dbl_mb(i_par+28*maxpar))
c
      return
      end
      subroutine qhop_dsite2(ga,rsa,rind,
     + ptseq,ptarat,arat,deq,aneq,ptpar,tunnel,defe12,zpef,
     + tdsgl,tstval,e12fxy)
c
c $Id$
c
      implicit none
c
#include "qhop_common.fh"
c
      integer ga,rind,rsa
c
      real*8 tunnel(maxpar,7),defe12(maxpar,7),zpef(maxpar,3),
     + tdsgl(maxpar,6),tstval(maxpar,5),e12fxy(maxpar,3)
      real*8 deq(mxar),aneq(mxar)
      integer ptseq(maxseq,2),ptpar(maxpar),arat(mxar,4),ptarat(mxseq)
c
      rind=arat(ptarat(ptseq(rsa,1))+ga-
     + ptseq(rsa,2),1)
c      
      return
      end
