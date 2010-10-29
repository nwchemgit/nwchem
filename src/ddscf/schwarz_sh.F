      double precision function schwarz_shell(ish, jsh)
C$Id$
      implicit none
#include "cschwarz.fh"
#include "mafdecls.fh"
      integer ish, jsh
      double precision value
#include "itri.fh"
c
      call comp4_extract(int_mb(k_csh), itri(ish,jsh), value, nb_per_i)
      schwarz_shell = value
c
      end
