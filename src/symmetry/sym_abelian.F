      logical function sym_abelian_group(geom)
*
* $Id: sym_abelian.F,v 1.2 1999-07-27 21:06:58 d3e129 Exp $
*
      implicit none
      integer geom              ! [input]
#include "inp.fh"
c
c     Return true if the point group of the geometry is abelian
c
      character*8 group
      integer nab, ind
      parameter (nab = 8)
      character*4 ab(nab)
      data ab/ 'C1','Cs','Ci','C2','D2','C2v','C2h','D2h'/
c
      call sym_group_name(geom, group)
c
      sym_abelian_group = inp_match(nab,.false.,group,ab,ind)
c
      end
