      integer function sym_abelian_axis(geom,axis)
!
! $Id$
!
! Returns irrep of a given axis if real Abelian,
! otherwise returns -1.
!
! Dipole:
! x=1 y=2 z=3
!
! Quadrupole:
! xx=4 yy=5 zz=6 xy=7 xz=8 yz=9
!
! Octupole:
! xxx=10 yyy=11 zzz=12 xxy=13 xxz=14 yyx=15 yyz=16 zzx=17 zzy=18 xyz=19
!
! Move this to symmetry directory!
!
      implicit none
#include "inp.fh"
      integer geom              ! [input]
      integer axis              ! [input]
      character*8 group
      integer nab, ind
      parameter (nab = 8)
      character*4 ab(nab)
      integer air(19,8)
      data ab/'C1','Cs','Ci','C2','D2','C2v','C2h','D2h'/
!              x y z x y z x x y x y z x x y y z z x
!                    x y z y z z x y z x x y y z z y
!                                x y z y z x z x y z
      data air/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     1         0,0,1,0,0,0,0,1,1,0,0,0,0,1,0,1,0,0,1,
     2         1,1,1,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,
     3         1,1,0,0,0,0,0,1,1,1,1,0,1,0,1,0,1,1,0,
     4         3,2,1,0,0,0,1,2,3,3,2,1,2,1,3,1,3,2,0,
     5         2,3,0,0,0,0,1,2,3,2,3,0,3,0,2,0,2,3,1,
     6         3,3,1,0,0,0,0,2,2,3,3,1,3,1,3,1,3,3,1,
     7         7,5,3,0,0,0,2,4,6,7,5,3,5,3,7,3,7,5,1/
!
      if ((axis.lt.1).or.(axis.gt.19)) then
         sym_abelian_axis = -1
         return
      endif
      call sym_group_name(geom,group)
!
      if (inp_match(nab,.false.,group,ab,ind)) then
!
!        real Abelian
!
         sym_abelian_axis = air(axis,ind)
      else
!
!        non-Abelian or complex Abelian
!
         sym_abelian_axis = -1
         return
      endif
      end
