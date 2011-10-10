      subroutine argos_diana_rdchg(x,q,t)
c
      implicit none
c
#include "argos_diana_common.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c
      real*8 x(nsa,3),q(nsa,nq)
      character*2 t(nsa)
c
      integer i,j
      real*8 xd(3)
c
      if(me.eq.0) then
c
      do 1 i=1,nsa
      read(lfnchg,1000) t(i),(x(i,j),j=1,3),(q(i,j),j=1,nq)
 1000 format(a2,3f12.6,4f12.6)
      x(i,1)=0.1d0*x(i,1)
      x(i,2)=0.1d0*x(i,2)
      x(i,3)=0.1d0*x(i,3)
      if(i.eq.1) then
      do 2 j=1,3
      xmin(j)=x(i,j)
      xmax(j)=x(i,j)
    2 continue
      else
      do 3 j=1,3
      if(xmin(j).gt.x(i,j)) xmin(j)=x(i,j)
      if(xmax(j).lt.x(i,j)) xmax(j)=x(i,j)
    3 continue
      endif
    1 continue
c
      do 4 j=1,3
      if(iesppb.eq.0.or.rcut.lt.0.0d0.or.
     + (iesppb.eq.1.and.j.ne.3).or.
     + (iesppb.eq.2.and.j.eq.3)) then
      xd(j)=xmax(j)-xmin(j)+2.0d0*abs(rcut)
      xmin(j)=xmin(j)-abs(rcut)
      else
      xd(j)=box(j)
      xmin(j)=0.5d0*(xmin(j)+xmax(j)-box(j))
      endif
    4 continue
c
      do 5 j=1,3
      spac(j)=1.0d0/(ngrid-1)
    5 continue
c
      do 6 j=1,3
      ngrd(j)=xd(j)/spac(j)+1
    6 continue
c
      do 7 j=1,3
      spac(j)=xd(j)/dble(ngrd(j))
    7 continue
c
      endif
c
c     broadcast to all nodes
c
      if(np.gt.1) then
      call ga_brdcst(mag_d01,x,3*nsa*ma_sizeof(mt_dbl,1,mt_byte),0)
      call ga_brdcst(mag_d02,xmin,3*ma_sizeof(mt_dbl,1,mt_byte),0)
      call ga_brdcst(mag_d03,ngrd,3*ma_sizeof(mt_int,1,mt_byte),0)
      call ga_brdcst(mag_d04,spac,3*ma_sizeof(mt_dbl,1,mt_byte),0)
      call ga_brdcst(mag_d05,q,nsa*nq*ma_sizeof(mt_dbl,1,mt_byte),0)
      endif
c
      return
      end
c $Id$
