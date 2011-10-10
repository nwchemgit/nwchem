      subroutine argos_prop_rdrest(lfn,fil)
c
      implicit none
c
#include "argos_prop.fh"
c
      integer lfn
      character*255 fil
c
      character*18 string
      integer i
c
      if(me.ne.0) return
c
      open(unit=lfn,file=fil(1:index(fil,' ')-1),
     + status='old',form='formatted',err=9999)
      rewind(lfn)
c
    1 continue
      read(lfn,1000,end=9998) string
 1000 format(a18)
      if(string.ne.'restart properties') goto 1
c
      read(lfn,1001) maxp,nsum,nsumt
 1001 format(3i7)
      read(lfn,1002) tsum,t2sum,tsumt,t2sumt
 1002 format(4e20.12)
      read(lfn,1003) (psum(i),i=1,maxp)
      read(lfn,1003) (p2sum(i),i=1,maxp)
      read(lfn,1003) (pslop(i),i=1,maxp)
      read(lfn,1003) (psumt(i),i=1,maxp)
      read(lfn,1003) (p2sumt(i),i=1,maxp)
      read(lfn,1003) (pslopt(i),i=1,maxp)
 1003 format(4e20.12)
c
      return
 9998 continue
      call md_abort
     + ('Unable to read restart file in argos_prop_rdrest ',me)
      return
 9999 continue
      call md_abort
     + ('Unable to open restart file in argos_prop_rdrest ',me)
      return
      end
c $Id$
