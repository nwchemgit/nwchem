      subroutine pre_collaps(xs,id,msa,nsa,touch,mode)
c
      implicit none
c
      integer msa,nsa,mode
      integer id(msa)
      real*8 xs(3,msa)
      real*8 touch
c
      integer mid,i,j,k,l,m,n
      real*8 dist,f,fmax,fmin,dx(3)
c
      mid=0
      do 1 i=1,nsa
      mid=max(mid,id(i))
    1 continue
c
      do 8 m=1,10
      do 2 i=1,mid
      dist=0.0d0
      dx(1)=0.0d0
      dx(2)=0.0d0
      dx(3)=0.0d0
      n=1
      do 3 j=1,nsa
      if(id(j).eq.i) then
      n=n+1
      dx(1)=dx(1)-xs(1,j)
      dx(2)=dx(2)-xs(2,j)
      dx(3)=dx(3)-xs(3,j)
      endif
    3 continue
      dx(1)=dx(1)/dble(n)
      dx(2)=dx(2)/dble(n)
      dx(3)=dx(3)/dble(n)
      if(mode.eq.3) dx(3)=0.0d0
      fmax=1.0d0
      fmin=0.0d0
      do 6 l=1,10
      f=0.5d0*(fmax-fmin)
      do 4 j=1,nsa
      if(id(j).eq.i) then
      dist=0.0d0
      do 5 k=1,nsa
      if(id(k).ne.i) then
      dist=min(dist,(xs(1,j)+f*dx(1)-xs(1,k))**2
     + +(xs(2,j)+f*dx(2)-xs(2,k))**2
     + +(xs(3,j)+f*dx(3)-xs(3,k))**2)
      endif
    5 continue
      if(sqrt(dist).gt.touch) then
      fmin=f
      else
      fmax=f
      endif
      endif
    4 continue
      write(*,'(10i5)') m,i,l
    6 continue
      if(fmin.gt.0.0d0) then
      do 7 j=1,nsa
      if(id(j).eq.i) then
      xs(1,j)=xs(1,j)+fmin*dx(1)
      xs(2,j)=xs(2,j)+fmin*dx(2)
      xs(3,j)=xs(3,j)+fmin*dx(3)
      endif
    7 continue
      endif
    2 continue
    8 continue
c
      return
      end
