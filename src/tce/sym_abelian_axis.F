      integer function sym_abelian_axis(geom,axis)
c
c $Id: sym_abelian_axis.F,v 1.1 2003-07-25 17:56:55 sohirata Exp $
c
c Returns irrep of a given axis (x=1, y=2, or z=3) if real Abelian,
c otherwise returns -1.
c
c Move this to symmetry directory!
c
      implicit none
#include "inp.fh"
      integer geom              ! [input]
      integer axis              ! [input] 1 = x, 2 = y, 3 = z
      character*8 group
      integer nab, ind
      parameter (nab = 8)
      character*4 ab(nab)
      integer air(3,8)
      data ab/'C1','Cs','Ci','C2','D2','C2v','C2h','D2h'/
      data air/0,0,0,
     1         0,0,1,
     2         1,1,1,
     3         1,1,0,
     4         3,2,1,
     5         2,3,0,
     6         3,3,1,
     7         7,5,3/
c
      if ((axis.lt.1).or.(axis.gt.3)) then
         sym_abelian_axis = -1
         return
      endif
      call sym_group_name(geom,group)
c
      if (inp_match(nab,.false.,group,ab,ind)) then
c
c        real Abelian
c
         sym_abelian_axis = air(axis,ind)  
      else
c
c        non-Abelian or complex Abelian
c
         sym_abelian_axis = -1
         return
      endif
      end
