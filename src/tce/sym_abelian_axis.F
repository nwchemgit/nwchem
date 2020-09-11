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
! xx=4 xy=5 xz=6 yy=7 yz=8 zz=9
!
! Octupole:
! xxx=10 xxy=11 xxz=12 xyy=13 xyz=14 xzz=15 yyy=16 yyz=17 yzz=18 zzz=19
!
! Move this to symmetry directory!
!
      implicit none
#include "inp.fh"
#include "global.fh"            /* ga_nodeid() */
#include "util.fh"              /* util_print(..) */
#include "stdio.fh"             /* LuOut=6 */
      integer geom              ! [input]
      integer axis              ! [input]
      character*8 group
      integer nab, nac, ind
      parameter (nab = 8)
      parameter (nac = 19)
      character*4 ab(nab)
      character*4 ac(nac)
      integer air(nac,nab)
      data ab/'C1','Cs','Ci','C2','D2','C2v','C2h','D2h'/
!              x y z x x x y y z x x x x x x y y y z
!                    x y z y z z x x x y y z y y z z
!                                x y z y z z y z z z
      data ac/'x  ','y  ','z  ',
     1        'xx ','xy ','xz ','yy ','yz ','zz ',
     2        'xxx','xxy','xxz','xyy','xyz',
     3        'xzz','yyy','yyz','yzz','zzz'/
      data air/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     1         0,0,1,0,0,1,0,1,0,0,0,1,0,1,0,0,1,0,0,
     2         1,1,1,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,
     3         1,1,0,0,0,1,0,1,0,1,1,0,1,0,1,1,0,1,0,
     4         3,2,1,0,1,2,0,3,0,3,2,1,3,0,3,2,1,2,1,
     5         2,3,0,0,1,2,0,3,0,2,3,0,2,1,2,3,0,3,0,
     6         3,3,1,0,0,2,0,2,0,3,3,1,3,1,3,3,1,3,1,
     7         7,5,3,0,2,4,0,6,0,7,5,3,7,1,7,5,3,5,3/
      if ((axis.lt.1).or.(axis.gt.19)) then
         sym_abelian_axis = -1
         return
      endif
      call sym_group_name(geom,group)
      if (inp_match(nab,.false.,group,ab,ind)) then
!        real Abelian
         sym_abelian_axis = air(axis,ind)
         if ((ga_nodeid().eq.0).and.
     1       util_print('tcesym',print_debug)) then
             write(LuOut,100) axis,ab(ind),ac(axis),air(axis,:)
             call util_flush(LuOut)
         endif
      else
!        non-Abelian or complex Abelian
         sym_abelian_axis = -1
         return
      endif
  100 format(1x,'sym_abelian_axis: axis=',i2,' group=',a3,
     1          ' char=',a3,' irrep=',8i1)
      end
