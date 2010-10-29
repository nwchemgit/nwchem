      subroutine dia_rama()
c
c $Id$
c
      implicit none
c
#include "dia_params.fh"
#include "dia_common.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c
      character*255 string
      integer ilt,iht,jlt,jht,ndxt,ldt
c
      string=card(8:80)
      strram=string
c
c     get access to local trajectory and time
c
      call ga_distribution(ga_trj,me,il,ih,jl,jh)
      call ga_access(ga_trj,il,ih,jl,jh,ndx,ld)
      call ga_distribution(ga_trt,me,ilt,iht,jlt,jht)
      call ga_access(ga_trt,ilt,iht,jlt,jht,ndxt,ldt)
c
      call dia_rama1(int_mb(i_iram),dbl_mb(ndx),dbl_mb(ndxt))
c
      return
      end
      subroutine dia_rama1(iram,xtrj,trt)
c
      implicit none
c
#include "dia_params.fh"
#include "dia_common.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c
      integer iram(msgm,7)
      real*8 xtrj(ntrj,3,nfrdim),trt(6,nfrdim)
c
      real*8 dia_torsion
      external dia_torsion
c
      real*8 phi,psi,rhb1,rhb2,dist
      integer i,j,k,ihb1,ihb2,ndxram
      logical lphi,lpsi
c
      write(filram,1000) strram(1:index(strram,' ')-1),
     + '.ram '
 1000 format(a,a)
      open(unit=lfnram,file=filram(1:index(filram,' ')-1),
     + form='formatted',status='unknown')
      rewind(lfnram)
      ndxram=me*nfrdim*(73*nsgm+20)
      call fseek(lfnram,ndxram,0)
      do 1 k=1,nfrme
c      write(filram,1000) strram(1:index(strram,' ')-1),
c     + me*nfrdim+k,'.ram '
c 1000 format(a,i5.5,a)
      call dia_rama2(iram,xtrj(1,1,k),trt(1,k))
    1 continue
c
      close(unit=lfnram)
c
      return
      end
      subroutine dia_rama2(iram,x,t)
c
      implicit none
c
#include "dia_params.fh"
#include "dia_common.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c
      integer iram(msgm,7)
      real*8 x(ntrj,3),t(6)
c
      real*8 dia_torsion
      external dia_torsion
c
      real*8 phi,psi,rhb1,rhb2,dist
      integer i,j,k,ihb1,ihb2
      logical lphi,lpsi
c
      write(lfnram,1000) t(1),nsgm
 1000 format(f12.6,i7)
c
      rhb1=0.0d0
      rhb2=0.0d0
c
      do 101 i=1,nsgm
      lphi=iram(i,1).gt.0
      lpsi=iram(i,5).gt.0
      phi=0.0d0
      psi=0.0d0
      if(lphi) phi=dia_torsion(dbl_mb(i_xdat),msa,
     + dbl_mb(i_wdat),mwm,mwa,
     + iram(i,1),0,iram(i,2),0,iram(i,3),0,iram(i,4),0)
      if(lpsi) psi=dia_torsion(dbl_mb(i_xdat),msa,
     + dbl_mb(i_wdat),mwm,mwa,
     + iram(i,2),0,iram(i,3),0,iram(i,4),0,iram(i,5),0)
      ihb1=0
      ihb2=0
c
      if(iram(i,6).gt.0) then
      do 2 j=1,nsgm
      if(i.ne.j.and.iram(j,7).gt.0) then
      dist=(x(iram(i,6),1)-x(iram(j,7),1))**2+
     + (x(iram(i,6),2)-x(iram(j,7),2))**2+
     + (x(iram(i,6),3)-x(iram(j,7),3))**2
      if(ihb1.eq.0.or.dist.lt.rhb1) then
      ihb1=j
      rhb1=dist
      endif
      endif
    2 continue
      endif
c
      if(iram(i,7).gt.0) then
      do 3 j=1,nsgm
      if(i.ne.j.and.iram(j,6).gt.0) then
      dist=(x(iram(i,7),1)-x(iram(j,6),1))**2+
     + (x(iram(i,7),2)-x(iram(j,6),2))**2+
     + (x(iram(i,7),3)-x(iram(j,6),3))**2
      if(ihb2.eq.0.or.dist.lt.rhb2) then
      ihb2=j
      rhb2=dist
      endif
      endif
    3 continue
      endif
c
      if(ihb1.gt.0) then
      rhb1=sqrt(rhb1)
      if(rhb1.lt.rhbmin.or.rhb1.gt.rhbmax) then
      ihb1=0
      rhb1=0.0d0
      endif
      endif
c
      if(ihb2.gt.0) then
      rhb2=sqrt(rhb2)
      if(rhb2.lt.rhbmin.or.rhb2.gt.rhbmax) then
      ihb2=0
      rhb2=0.0d0
      endif
      endif
c
      if(me.eq.0) then
      if(lphi.or.lpsi) then
      write(lfnram,1001) i,lphi,lpsi,phi,psi,ihb1,ihb2,rhb1,rhb2
 1001 format(i7,1x,2l1,2f12.6,2i7,2f12.6)
      endif
      endif
  101 continue
c
      return
      end
      subroutine dia_ramhdr(sgmnam,iram)
c
      implicit none
c
#include "dia_common.fh"
#include "global.fh"
#include "msgids.fh"
c
      character*16 sgmnam(nsa)
      integer iram(msgm,7)
c
      integer i,j,number
c
      number=0
      do 1 i=1,nsa
      if(i.eq.nsa) goto 2
      if(sgmnam(i)(11:16).ne.sgmnam(i+1)(11:16)) goto 2
      goto 1
    2 continue
      read(sgmnam(i)(11:16),'(i6)') j
      if(iram(j,1).gt.0.or.iram(j,5).gt.0) number=number+1
    1 continue
c
      write(lfnram,1000) number
 1000 format(i5)
c
      do 3 i=1,nsa
      if(i.eq.nsa) goto 4
      if(sgmnam(i)(11:16).ne.sgmnam(i+1)(11:16)) goto 4
      goto 3
    4 continue
      read(sgmnam(i)(11:16),'(i6)') j
      if(iram(j,1).gt.0.or.iram(j,5).gt.0)
     + write(lfnram,1002) sgmnam(i)(1:5),sgmnam(i)(11:16)
 1002 format(a5,a6)
    3 continue
c
      return
      end
