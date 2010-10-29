      subroutine dia_evaluate()
c
c $Id$
c
c     evaluate the defined properties
c     -------------------------------
c
      implicit none
c
#include "dia_common.fh"
#include "dia_params.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "msgids.fh"
#include "util.fh"
c
      integer itval
      integer ilt,iht,jlt,jht,ndxt,ldt
c
      filord=card(8:80)
      fildat=filord(1:index(filord,' ')-1)//'.dat'
      filord=filord(1:index(filord,' ')-1)//'.ord'
c
      call ga_distribution(ga_trj,me,il,ih,jl,jh)
      call ga_access(ga_trj,il,ih,jl,jh,ndx,ld)
c
c     allocate global array for the properties
c     ----------------------------------------
c
      if(ntval.gt.0) then
      if(.not.ga_destroy(ga_tval))
     + call md_abort('Failed to destroy ga_tval',0)
      endif
      ntval=ndist+nangl+ntors+nimps
c
      if(ntval.gt.0) then
c
      if(.not.ga_create(mt_dbl,ntval,nfrtot,'tval',ntval,nfrdim,
     + ga_tval)) call md_abort('Failed to create ga_tval',0)
c
      itval=1
c
      call ga_distribution(ga_tval,me,ilt,iht,jlt,jht)
      call ga_access(ga_tval,ilt,iht,jlt,jht,ndxt,ldt)
c
c     evaluate distances
c     ------------------
c
      if(ndist.gt.0)
     + call dia_evaldist(dbl_mb(ndx),dbl_mb(ndxt),itval)
      itval=itval+ndist
c
c     evaluate angles
c     ---------------
c
      if(nangl.gt.0)
     + call dia_evalangl(dbl_mb(ndx),dbl_mb(ndxt),itval)
      itval=itval+nangl
c
c     evaluate torsions
c     -----------------
c
      if(ntors.gt.0)
     + call dia_evaltors(dbl_mb(ndx),dbl_mb(ndxt),itval)
      itval=itval+ntors
c
      endif
c
      if(nord.gt.0)
     + call dia_evalord(dbl_mb(ndx))
c
      call dia_export()
c
c     reset
c
      ndist=0
      nangl=0
      ntors=0
      nord=0
c
      return
      end
      subroutine dia_evaldist(xtrj,tval,itval)
c
      implicit none
c
#include "msgids.fh"
#include "global.fh"
#include "dia_common.fh"
c
      real*8 dia_distance
      external dia_distance
c
      integer itval
      real*8 xtrj(ntrj,3,nfrdim),tval(ntval,nfrdim)
c
      real*8 wtrj(3)
c
      integer i,k
c
      do 1 i=1,ndist
      do 2 k=1,nfrme
      tval(itval+i-1,k)=dia_distance(xtrj(1,1,k),ntrj,wtrj,1,3,
     + idist(i,1,1),idist(i,1,2),idist(i,2,1),idist(i,2,2),box)
    2 continue
    1 continue
c
      return
      end
      subroutine dia_evalangl(xtrj,tval,itval)
c
      implicit none
c
#include "msgids.fh"
#include "global.fh"
#include "dia_common.fh"
c
      real*8 dia_angle
      external dia_angle
c
      integer itval
      real*8 xtrj(ntrj,3,nfrdim),tval(ntval,nfrdim)
c
      real*8 wtrj(3)
c
      integer i,k
c
      do 1 i=1,nangl
      do 2 k=1,nfrme
      tval(itval+i-1,k)=dia_angle(xtrj(1,1,k),ntrj,wtrj,1,3,
     + iangl(i,1,1),iangl(i,1,2),iangl(i,2,1),iangl(i,2,2),
     + iangl(i,3,1),iangl(i,3,2))
    2 continue
    1 continue
c
      return
      end
      subroutine dia_evaltors(xtrj,tval,itval)
c
      implicit none
c
#include "msgids.fh"
#include "global.fh"
#include "dia_common.fh"
c
      real*8 dia_torsion
      external dia_torsion
c
      integer itval
      real*8 xtrj(ntrj,3,nfrdim),tval(ntval,nfrdim)
c
      real*8 wtrj(3)
c
      integer i,k
c
      do 1 i=1,ntors
      do 2 k=1,nfrme
      tval(itval+i-1,k)=dia_torsion(xtrj(1,1,k),ntrj,wtrj,1,3,
     + itors(i,1,1),itors(i,1,2),itors(i,2,1),itors(i,2,2),
     + itors(i,3,1),itors(i,3,2),itors(i,4,1),itors(i,4,2))
    2 continue
    1 continue
c
      return
      end
      subroutine dia_export()
c
      implicit none
c
#include "msgids.fh"
#include "global.fh"
#include "dia_common.fh"
#include "mafdecls.fh"
c
      integer i,j,ilt,iht,jlt,jht,nfr
      integer id1,id2,idexp,jtors
      logical lopen
      integer ntdim
c
      ntdim=ntval+1
      if(ntval.le.0) return
c
      call ga_sync()
c
      if(me.eq.0) then
c
      inquire(unit=lfndat,opened=lopen)
      if(.not.lopen)
     + open(unit=lfndat, file=fildat(1:index(fildat,' ')-1),
     + form='formatted',status='unknown')
c
      write(lfndat,1000) ntdim,nfrtot
 1000 format(i5,i10,/,'Time')
      if(ndist.gt.0) then
      do 1 i=1,ndist
      write(lfndat,1001) 'Distance ',cdist(i)
 1001 format(a,a)
    1 continue
      endif
      if(nangl.gt.0) then
      do 2 i=1,nangl
      write(lfndat,1002) 'Angle    ',cangl(i)
 1002 format(a,a)
    2 continue
      endif
      if(ntors.gt.0) then
      do 3 i=1,ntors
      write(lfndat,1003) 'Torsion  ',ctors(i)
 1003 format(a,a)
    3 continue
      endif
      if(nimps.gt.0) then
      endif
c
      if(.not.ma_push_get(mt_dbl,nfrdim*ntdim,'tv',l_tval,i_tval))
     + call md_abort('Failed to allocate tv',0)
c
      do 6 i=0,np-1
      call ga_distribution(ga_tval,i,ilt,iht,jlt,jht)
      call ga_get(ga_tval,ilt,iht,jlt,jht,dbl_mb(i_tval+1),ntdim)
      call ga_distribution(ga_trt,i,ilt,iht,jlt,jht)
      call ga_get(ga_trt,1,1,jlt,jht,dbl_mb(i_tval),ntdim)
      nfr=min(nfrdim,nfrtot-i*nfrdim)
      call dia_wrtdat(dbl_mb(i_tval),ntdim,nfr)
    6 continue
c
      if(.not.ma_pop_stack(l_tval))
     + call md_abort('Failed to deallocate tv',0)
c
      endif
c
      return
      end
      subroutine dia_wrtdat(tval,ntdim,nfr)
c
      implicit none
c
#include "msgids.fh"
#include "global.fh"
#include "dia_common.fh"
#include "mafdecls.fh"
c
      integer ntdim,nfr
      real*8 tval(ntdim,nfrdim)
c
      integer i,j
c
      do 1 j=1,nfr
      write(lfndat,1000)
 1000 format('frame')
      write(lfndat,1001) (tval(i,j),i=1,ntdim)
 1001 format(6e12.5)
    1 continue
c
      return
      end
      subroutine dia_evalord(xtrj)
c
      implicit none
c
#include "msgids.fh"
#include "global.fh"
#include "dia_common.fh"
c
      real*8 xtrj(ntrj,3,nfrdim)
c
      integer i,j,k,ndordr,ndxo,inum
      real*8 s,d
      logical lopen
c
      do 4 i=1,nord
      rord(i)=0.0d0
    4 continue
c
      do 1 i=1,nord
      ndxo=iord(i,1,3)
      ndordr=iord(i,3,3)
      do 2 k=1,nfrme
      do 3 j=2,ndordr+1
      d=(xtrj(iord(i,j,1),1,k)-xtrj(iord(i,1,1),1,k))**2+
     +(xtrj(iord(i,j,1),2,k)-xtrj(iord(i,1,1),2,k))**2+
     +(xtrj(iord(i,j,1),3,k)-xtrj(iord(i,1,1),3,k))**2
      s=3.0d0*(xtrj(iord(i,j,1),3,k)-xtrj(iord(i,1,1),3,k))**2/d-1.0d0
      rord(ndxo)=rord(ndxo)+s
    3 continue
    2 continue
    1 continue
c
      call ga_dgop(mag_d09,rord,mord,'+')
c
      if(me.eq.0) then
      inquire(unit=lfnord,opened=lopen)
      if(.not.lopen)
     + open(unit=lfnord, file=filord(1:index(filord,' ')-1),
     + form='formatted',status='unknown')
      k=0
      j=0
      do 5 i=1,nord
      if(iord(i,1,3).ne.k) then
      k=iord(i,1,3)
      j=j+1
      endif
    5 continue
c
      write(lfnord,1000) j,nord,nfrtot
 1000 format(2i5,i10)
c
      k=0
      j=0
      do 6 i=1,nord
      if(iord(i,1,3).ne.k) then
      k=iord(i,1,3)
      j=0
      endif
      j=j+1
      write(lfnord,1001) k,j,cord(i)
 1001 format(2i5,a)
    6 continue
c
      k=0
      do 7 i=1,nord
      if(iord(i,1,3).ne.k) then
      inum=2*iord(i,3,3)*iord(i,2,3)*nfrtot
      k=iord(i,1,3)
      write(lfnord,1002) k,iord(i,3,3),iord(i,2,3),rord(k)/dble(inum)
 1002 format(3i5,f12.5)
      endif
    7 continue
c
      close(unit=lfnord,status='keep')
c
      endif
c
      return
      end      
