      SUBROUTINE tool_randm(iseed,x)

      implicit none

      integer i,iseed,itoz,itozz,mz,mult

      real*8 x,add,dimax,ddimax
c     real*8 rand

      logical newjob

      dimension mz(250)

      save newjob,itoz,dimax,ddimax

      data newjob/.true./

      if(newjob)then
       if(mod(iseed,2).eq.0)iseed=iseed+1
       mult=65539
       add=2147483648.0d00
       dimax=1.0d00/add
       ddimax=0.50d00*dimax
       do i=1,250
        x=rand(iseed)
        mz(i)=x*iseed
       enddo
       itoz=1
       newjob=.false.
      else
       itoz=itoz+1
       if(itoz.gt.250)itoz=itoz-250
       itozz=itoz+103
       if(itozz.gt.250)itozz=itozz-250
       mz(itoz)=ieor(mz(itoz),mz(itozz))
       x=mz(itoz)*dimax+ddimax
       x=2.0d00*x
      endif

      return

      END
c $Id$
