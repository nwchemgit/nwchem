      double precision function schwarz_atom(iat, jat)
C$Id: schwarz_atom.F,v 1.4 1995-02-02 23:12:10 d3g681 Exp $
      implicit none
#include "cschwarz.fh"
#include "mafdecls.fh"
      integer iat, jat
      double precision value
#include "itri.fh"
c
      call comp4_extract(int_mb(k_cat), itri(iat,jat), value, nb_per_i)
      schwarz_atom = value
c
      end
