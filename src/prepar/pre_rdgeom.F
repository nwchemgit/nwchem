      logical function pre_rdgeom(irtdb,csa,isar,isgm,isfnd,xs,vs,
     + msa,nsa,nmis,msm,nsm)
c
c $Id$
c
      implicit none
c
#include "geom.fh"
c
      real*8 cau2nm
      parameter(cau2nm=5.29177249d-02)
c
      integer irtdb,msa,nsa,nmis,msm,nsm
      character*16 csa(msa)
      integer isar(msa),isgm(msa),isfnd(msa)
      real*8 xs(3,msa),vs(3,msa)
c
      integer igeom,i
      character*16 tag
      real*8 chrg,coord(3)
c
      if(.not.geom_create(igeom,'geometry'))
     + call md_abort('pre_rdgeom: geom_create failed',9999)
      if(.not.geom_rtdb_load(irtdb,igeom,'geometry'))
     + call md_abort('pre_rdgeom: geom_rtdb_load failed',9999)
c
      do 1 i=1,msa
      if(.not.geom_cent_get(igeom,i,tag,coord,chrg))
     + call md_abort('pre_rdgeom: geom_cent_get failed',9999)
      xs(1,i)=cau2nm*coord(1)
      xs(2,i)=cau2nm*coord(2)
      xs(3,i)=cau2nm*coord(3)
      vs(1,i)=0.0d0
      vs(2,i)=0.0d0
      vs(3,i)=0.0d0
      isfnd(i)=1
      isar(i)=0
    1 continue
      nsa=msa
      nmis=0
c
      if(.not.geom_destroy(igeom))
     + call md_abort('pre_rdgeom: geom_destroy failed',9999)
c
      pre_rdgeom=.true.
      return
c
 9999 continue
      pre_rdgeom=.false.
      return
      end
