      subroutine argos_cafe_shkp(nums,xs,ys,isgan,isat,ishp,
     + ipmf,jpmf,rpmf,xpmf,ypmf,wgt)
c
      implicit none
c
#include "argos_cafe_common.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c
      integer nums
      real*8 xs(msa,3),ys(msa,3)
      integer isgan(msa),isat(msa),ishp(msa)
      integer ipmf(numpmf,8),jpmf(numpmf,4,npmfa)
      real*8 rpmf(3,6,numpmf),xpmf(4,4,numpmf),ypmf(4,3,numpmf)
      real*8 wgt(mat,mset),ft(4,3),rv(2),w(4)
      integer ip(3),it(2),ign(4),iml(4),idt(4),ich(4),ihp(4),isg(4)
c
      integer i,j,k,l,m
      logical skip
c
c      epmf=zero
      skip=.true.
c
      do 1 i=1,numpmf
      if(ipmf(i,6).eq.2.and.
     + (ipmf(i,1).eq.1.or.ipmf(i,1).eq.5.or.ipmf(i,1).eq.6)) then
      skip=.false.
      do 2 j=1,4
      xpmf(j,1,i)=zero
      xpmf(j,2,i)=zero
      xpmf(j,3,i)=zero
      xpmf(j,4,i)=zero
      ypmf(j,1,i)=zero
      ypmf(j,2,i)=zero
      ypmf(j,3,i)=zero
    2 continue
      m=0
      if(ipmf(i,1).eq.1.or.ipmf(i,1).eq.5.or.ipmf(i,1).eq.6) m=2
      if(ipmf(i,1).eq.2) m=3
      if(ipmf(i,1).eq.3) m=4
      if(ipmf(i,1).eq.4) m=4
      if(ipmf(i,1).ne.1) call md_abort('SHAKE on distance only',me)
      do 3 j=1,m
      do 4 k=1,ipmf(i,j+1)
      do 5 l=1,nums
      if(isgan(l).eq.jpmf(i,j,k)) then
      xpmf(j,1,i)=xpmf(j,1,i)+xs(l,1)
      xpmf(j,2,i)=xpmf(j,2,i)+xs(l,2)
      xpmf(j,3,i)=xpmf(j,3,i)+xs(l,3)
      xpmf(j,4,i)=xpmf(j,4,i)+wgt(isat(l),iset)
      ypmf(j,1,i)=ypmf(j,1,i)+ys(l,1)
      ypmf(j,2,i)=ypmf(j,2,i)+ys(l,2)
      ypmf(j,3,i)=ypmf(j,3,i)+ys(l,3)
      endif
    5 continue
    4 continue
      xpmf(j,1,i)=xpmf(j,1,i)/dble(ipmf(i,j+1))
      xpmf(j,2,i)=xpmf(j,2,i)/dble(ipmf(i,j+1))
      xpmf(j,3,i)=xpmf(j,3,i)/dble(ipmf(i,j+1))
      ypmf(j,1,i)=ypmf(j,1,i)/dble(ipmf(i,j+1))
      ypmf(j,2,i)=ypmf(j,2,i)/dble(ipmf(i,j+1))
      ypmf(j,3,i)=ypmf(j,3,i)/dble(ipmf(i,j+1))
      do 14 k=1,ipmf(i,j+1)
      do 15 l=1,nums
      if(isgan(l).eq.jpmf(i,j,k)) then
      xs(l,1)=xs(l,1)-xpmf(j,1,i)
      xs(l,2)=xs(l,2)-xpmf(j,2,i)
      xs(l,3)=xs(l,3)-xpmf(j,3,i)
      endif
   15 continue
   14 continue
    3 continue
      endif
    1 continue
c
      if(skip) return
c
      call ga_dgop(mcf_59,xpmf,16*numpmf,'+')
      call ga_dgop(mcf_63,ypmf,12*numpmf,'+')
c
      do 6 i=1,numpmf
c
      if(ipmf(i,6).eq.2.and.ipmf(i,1).eq.1) then
      do 7 j=1,2
      w(j)=xpmf(j,4,i)
      ign(j)=j
      iml(j)=j
      idt(j)=0
      isg(j)=1
      ihp(i)=0
      ich(j)=1
      do 8 k=1,3
      ft(j,k)=zero
    8 continue
    7 continue
      it(1)=1
      ip(1)=1
      ip(2)=2
      ip(3)=1
      call argos_cafe_shks(1,it,1,3,ip,rpmf,rv,2,4,ign,iml,isg,idt,ihp,
     + xpmf(1,1,i),ypmf(1,1,i),w,int_mb(i_lseq))
      m=0
      if(ipmf(i,1).eq.1.or.ipmf(i,1).eq.5.or.ipmf(i,1).eq.6) m=2
      if(ipmf(i,1).eq.2) m=3
      if(ipmf(i,1).eq.3) m=4
      if(ipmf(i,1).eq.4) m=4
      if(ipmf(i,1).ne.1) call md_abort('SHAKE on distance only',me)
      do 9 j=1,m
      do 10 k=1,ipmf(i,j+1)
      do 11 l=1,nums
      if(isgan(l).eq.jpmf(i,j,k)) then
      xs(l,1)=xs(l,1)+xpmf(j,1,i)
      xs(l,2)=xs(l,2)+xpmf(j,2,i)
      xs(l,3)=xs(l,3)+xpmf(j,3,i)
      endif
   11 continue
   10 continue
    9 continue
c
      endif
c      if(ipmf(i,1).eq.3) then
c      call md_abort('Angle pmf not implemented',0)
c      endif
c      if(ipmf(i,1).eq.4) then
c      call md_abort('Torsion pmf not implemented',0)
c      endif
    6 continue
c
      return
      end
c $Id$
