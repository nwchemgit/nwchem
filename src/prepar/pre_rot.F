      logical function pre_rot(irot,jrot,arot,xrot,is,xs,msa,nsa)
c
c $Id$
c
      implicit none
c
      integer msa,nsa,irot,jrot
      integer is(msa)
      real*8 xs(3,msa),xrot(3),arot
c
      real*8 cog(3),xr(3),x(3),y(3)
c
      integer i,j,number
c
      if(nsa.le.0.or.irot.le.0.or.jrot.gt.nsa) then
      pre_rot=.false.
      return
      endif
c
c     find center of geometry
c
c      number=0
c      do 1 i=1,nsa
c      if(is(i).ge.irot.and.is(i).le.jrot) then
c      number=number+1
c      do 2 j=1,3
c      cog(j)=cog(j)+xs(j,i)
c    2 continue
c      endif
c    1 continue
c      cog(1)=cog(1)/dble(number)
c      cog(2)=cog(2)/dble(number)
c      cog(3)=cog(3)/dble(number)
c      xr(1)=xrot(1)+cog(1)
c      xr(2)=xrot(2)+cog(2)
c      xr(3)=xrot(3)+cog(3)
c
c     rotate
c
      cog(1)=0.0d0
      cog(2)=0.0d0
      cog(3)=0.0d0
      do 3 i=1,nsa
      if(is(i).ge.irot.and.is(i).le.jrot) then
      x(1)=xs(1,i)
      x(2)=xs(2,i)
      x(3)=xs(3,i)
      call rotate(cog,xrot,arot,x,y)
      xs(1,i)=y(1)
      xs(2,i)=y(2)
      xs(3,i)=y(3)
      endif
    3 continue
c
c      print*,'rotation done ',arot
      pre_rot=.true.
c
      return
      end
