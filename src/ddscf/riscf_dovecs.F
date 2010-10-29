      subroutine riscf_dovecs ( doit, g_vec )
C$Id$
C     switches on use of MO vectors (Ken Dyalls proposal) for RI-SCF
      implicit none
#include "cfock.fh"
      logical doit
      integer g_vec

      dovecs = doit
      g_vecs = g_vec

      return
      end
