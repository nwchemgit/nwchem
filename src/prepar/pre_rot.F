      logical function pre_rot(irot,jrot,arot,xrot,is,xs,msa,nsa)
c
c $Id: pre_rot.F,v 1.2 2001-08-29 23:04:05 d3j191 Exp $
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
      number=0
      do 1 i=1,nsa
      if(is(i).ge.irot.and.is(i).le.jrot) then
      number=number+1
      do 2 j=1,3
      cog(j)=cog(j)+xs(j,i)
    2 continue
      endif
    1 continue
      cog(1)=cog(1)/dble(number)
      cog(2)=cog(2)/dble(number)
      cog(3)=cog(3)/dble(number)
      xr(1)=xrot(1)+cog(1)
      xr(2)=xrot(2)+cog(2)
      xr(3)=xrot(3)+cog(3)
c
c     rotate
c
      do 3 i=1,nsa
      if(is(i).ge.irot.and.is(i).le.jrot) then
      x(1)=xs(1,i)
      x(2)=xs(2,i)
      x(3)=xs(3,i)
      call rotate(cog,xr,arot,x,y)
      xs(1,i)=y(1)
      xs(2,i)=y(2)
      xs(3,i)=y(3)
      endif
    3 continue
c
      print*,'rotation done ',arot
      pre_rot=.true.
c
      return
      end
