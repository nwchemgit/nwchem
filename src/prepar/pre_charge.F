      logical function pre_charge(irtdb,source,latm,qatm,matm,natm,
     + lbnd,mbnd,nbnd)
c
      implicit none
c
#include "msgids.fh"
#include "mafdecls.fh"
#include "geom.fh"
#include "rtdb.fh"
c
      logical pre_task1
      externall pre_task1
c
c     in  : latm(2,matm)      = atomic number
c           matm              = dimension atom list
c           lbnd(1:2,mbnd)    = bond indices
c           mbnd              = dimension bond list
c           nbnd              = length bond list
c
c     out : qatm(matm)        = guestimated atomic partial charges
c
      character*80 source
      integer irtdb
      integer matm,natm
      integer latm(5,matm)
      real*8 qatm(matm)
      integer mbnd,nbnd
      integer lbnd(2,mbnd)
c
      integer i,j,ia,ja,it,jt,itmp,iab(20),nb,nox,nhx,itask
      integer igeom2
      character*16 tag
      real*8 coord(3),charge
c
      do 1 i=1,matm
      qatm(i)=0.0d0
    1 continue
c
      if(source(1:8).eq.'geometry') then
      itask=1
      call ga_brdcst(pre_i00,itask,ma_sizeof(mt_int,1,mt_byte),0)
      if(itask.eq.1) then
      if(.not.pre_task1(irtdb)) call errquit('pre_task1 failed',9999)
      if(.not.geom_create(igeom2,'resp:geometry'))
     + call errquit('pre_charge: geom_create failed',0)
      if(.not.geom_rtdb_load(irtdb,igeom2,'resp:geometry'))
     + call errquit('pre_charge: geom_rtdb_load failed',0)
      do 9 i=1,natm
      if(.not.geom_cent_get(igeom2,i,tag,coord,charge))
     + call errquit('pre_charge: geom_cent_get failed',i)
      qatm(i)=charge
    9 continue
      if(.not.geom_destroy(igeom2))
     + call errquit('pre_charge: geom_destroy failed',9999)
      endif
      return
      endif
c
      do 2 i=1,nbnd
      ia=lbnd(1,i)
      ja=lbnd(2,i)
      it=latm(2,ia)
      jt=latm(2,ja)
      if(it.ne.jt) then
c
c     smaller atom number in it/ia
c
      if(it.gt.jt) then
      itmp=ia
      ia=ja
      ja=itmp
      itmp=it
      it=jt
      jt=itmp
      endif
c
c     H - C : 0.05  : -0.05
c
      if(it.eq.1.and.jt.eq.6) then
      qatm(ia)=qatm(ia)+0.05d0
      qatm(ja)=qatm(ja)-0.05d0
      endif
c
c     H - N : 0.27  : -0.27
c
      if(it.eq.1.and.jt.eq.7) then
      qatm(ia)=qatm(ia)+0.27d0
      qatm(ja)=qatm(ja)-0.27d0
      endif
c
c     H - S : 0.19  : -0.19
c
      if(it.eq.1.and.jt.eq.16) then
      qatm(ia)=qatm(ia)+0.19d0
      qatm(ja)=qatm(ja)-0.19d0
      endif
c
c     C - N : 0.07  : -0.07
c
      if(it.eq.6.and.jt.eq.7) then
      qatm(ia)=qatm(ia)+0.07d0
      qatm(ja)=qatm(ja)-0.07d0
      endif
c
c     C - O : 0.30  : -0.30
c
      if(it.eq.6.and.jt.eq.8.and.latm(3,ja).eq.2) then
      qatm(ia)=qatm(ia)+0.30d0
      qatm(ja)=qatm(ja)-0.30d0
      endif
c
c     C = O : 0.57  : -0.57
c
      if(it.eq.6.and.jt.eq.8.and.latm(3,ja).eq.1) then
      qatm(ia)=qatm(ia)+0.57d0
      qatm(ja)=qatm(ja)-0.57d0
      endif
c
c     C - S : 0.11  : -0.11
c
      if(it.eq.6.and.jt.eq.16) then
      qatm(ia)=qatm(ia)+0.11d0
      qatm(ja)=qatm(ja)-0.11d0
      endif
c
      endif
    2 continue
c
      do 3 i=1,natm
      nb=0
      do 4 j=1,nbnd
      ia=lbnd(1,j)
      ja=lbnd(2,j)
      if(ia.eq.i) then
      nb=nb+1
      iab(nb)=ja
      else
      if(ja.eq.i) then
      nb=nb+1
      iab(nb)=ia
      endif
      endif
    4 continue
c
c     -COO (-)
c
      if(latm(2,i).eq.6.and.nb.eq.3) then
      nox=0
      do 5 j=1,nb
      if(latm(2,iab(j)).eq.8.and.latm(3,iab(j)).eq.1) nox=nox+1
    5 continue
      if(nox.eq.2) then
      do 6 j=1,nb
      if(latm(2,iab(j)).eq.8.and.latm(3,iab(j)).eq.1)
     + qatm(iab(j))=qatm(iab(j))-0.23
    6 continue
      qatm(i)=qatm(i)-0.54
      endif
      endif
c
c     -NH3 (+)
c
      if(latm(2,i).eq.7.and.nb.eq.4) then
      nhx=0
      do 7 j=1,nb
      if(latm(2,iab(j)).eq.1.and.latm(3,iab(j)).eq.1) nhx=nhx+1
    7 continue
      if(nhx.eq.3) then
      do 8 j=1,nb
      if(latm(2,iab(j)).eq.1.and.latm(3,iab(j)).eq.1)
     + qatm(iab(j))=qatm(iab(j))+0.07
    8 continue
      qatm(i)=qatm(i)+0.79
      endif
      endif
c
    3 continue
c
      pre_charge=.true.
      return
c
 9999 pre_charge=.false.
      return
      end
