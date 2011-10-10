      subroutine argos_cafe_dlist(nums,xs,fs,isgan,idrs,rdrs,xdrs)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c
      integer nums
      real*8 xs(msa,3),fs(msa,3)
      integer isgan(msa)
      integer idrs(ndrs,2)
      real*8 rdrs(ndrs,6),xdrs(ndrs,6)
c
      integer i,j,k
      real*8 d,di,dd,dn,dx(3),dfor,df(3),en
c
      call timer_start(17)
      edrs=zero
c
      do 1 i=1,ndrs
      do 2 j=1,6
      xdrs(i,j)=zero
    2 continue
      do 3 j=1,nums
      if(isgan(j).eq.idrs(i,1)) then
      xdrs(i,1)=xs(j,1)
      xdrs(i,2)=xs(j,2)
      xdrs(i,3)=xs(j,3)
      endif
      if(isgan(j).eq.idrs(i,2)) then
      xdrs(i,4)=xs(j,1)
      xdrs(i,5)=xs(j,2)
      xdrs(i,6)=xs(j,3)
      endif
    3 continue
    1 continue
      call timer_stop(17)
c
      call timer_start(18)
      call ga_dgop(mcf_59,xdrs,6*ndrs,'+')
      call timer_stop(18)
c
      call timer_start(19)
      do 4 i=1,ndrs
      dx(1)=xdrs(i,1)-xdrs(i,4)
      dx(2)=xdrs(i,2)-xdrs(i,5)
      dx(3)=xdrs(i,3)-xdrs(i,6)
      d=sqrt(dx(1)*dx(1)+dx(2)*dx(2)+dx(3)*dx(3))
      if(d.lt.tiny) then
      di=one
      else
      di=one/d
      endif
      dn=d-rdrs(i,2)
      k=0
      if(d.lt.rdrs(i,1)) then
      k=1
      dd=rdrs(i,1)-rdrs(i,2)
      en=rdrs(i,4)*(dn-half*dd)*dd
      dfor=rdrs(i,4)*di*dd
      elseif(d.gt.rdrs(i,3)) then
      k=4
      dd=rdrs(i,3)-rdrs(i,2)
      en=rdrs(i,5)*(dn-half*dd)*dd
      dfor=rdrs(i,5)*di*dd
      elseif(d.lt.rdrs(i,2)) then
      k=2
      en=half*rdrs(i,4)*dn*dn
      dfor=rdrs(i,4)*di*dn
      else
      k=3
      en=half*rdrs(i,5)*dn*dn
      dfor=rdrs(i,5)*di*dn
      endif
      edrs=edrs+en
      df(1)=dfor*dx(1)
      df(2)=dfor*dx(2)
      df(3)=dfor*dx(3)
c      if(me.eq.0) then
c      write(*,2222) k,d,(rdrs(i,j),j=1,3),en,dfor,df
c 2222 format(i3,4f12.6,e12.5,4f15.3)
c      endif
      do 5 j=1,nums
      if(isgan(j).eq.idrs(i,1)) then
      fs(j,1)=fs(j,1)-df(1)
      fs(j,2)=fs(j,2)-df(2)
      fs(j,3)=fs(j,3)-df(3)
      endif
      if(isgan(j).eq.idrs(i,2)) then
      fs(j,1)=fs(j,1)+df(1)
      fs(j,2)=fs(j,2)+df(2)
      fs(j,3)=fs(j,3)+df(3)
      endif
    5 continue
    4 continue
      call timer_stop(19)
c
      return
      end
c $Id$
