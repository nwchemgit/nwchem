      subroutine riscf_dovecs ( doit, g_vec )
C$Id: riscf_dovecs.F,v 1.1 1995-12-29 02:04:12 vg038 Exp $
C     switches on use of MO vectors (Ken Dyalls proposal) for RI-SCF
      implicit none
#include "cfock.fh"
      logical doit
      integer g_vec

      dovecs = doit
      g_vecs = g_vec

      return
      end
