      subroutine pre_collaps(xs,id,msa,nsa,touch,mode)
c
      implicit none
c
      integer msa,nsa,mode
      integer id(msa)
      real*8 xs(3,msa)
      real*8 touch
c
      logical first,found
      integer i,j,k,m
      integer nmol,large
      real*8 rt
      real*8 a,b,c,d,r,r1,r2,root,dt,dtran,tr(3),tran(3),xmax(3),xmin(3)
c
      write(*,'(/,a,i4,f12.6,/)') ' Collapsing using mode ',mode,touch
c
c     find total number of molecules
c
      nmol=0
      do 1 i=1,nsa
      nmol=max(nmol,id(i))
    1 continue
c
    2 continue
c
c     find molecule with largest possible tranlation
c
      large=0
      dtran=0.0d0
      root=0.0d0
c
      do 3 m=1,nmol
c
c     determine translation vector
c
      rt=-1.0d0
      found=.false.
      first=.true.
      do 4 i=1,nsa
      if(id(i).eq.m) then
      found=.true.
      if(first) then
      do 5 k=1,3
      xmax(k)=xs(k,i)
      xmin(k)=xs(k,i)
    5 continue
      first=.false.
      else
      do 6 k=1,3
      xmax(k)=max(xmax(k),xs(k,i))
      xmin(k)=min(xmin(k),xs(k,i))
    6 continue
      endif
      endif
    4 continue
      if(found) then
      do 7 k=1,3
      tr(k)=-0.5d0*(xmin(k)+xmax(k))
    7 continue
      if(mode.eq.3) tr(3)=0.0d0
c
      dt=sqrt(tr(1)*tr(1)+tr(2)*tr(2)+tr(3)*tr(3))
c
      do 8 i=1,nsa
      if(id(i).eq.m) then
      do 9 j=1,nsa
      if(id(j).ne.m) then
      a=0.0d0
      b=0.0d0
      c=0.0d0
      do 10 k=1,3
      a=a+tr(k)*tr(k)
      b=b+(xs(k,i)-xs(k,j))*tr(k)
      c=c+xs(k,i)*xs(k,i)+xs(k,j)*xs(k,j)-2.0*xs(k,i)*xs(k,j)
   10 continue
      b=2.0d0*b
      c=c-touch*touch
      if(abs(a).gt.1.0d-5) then
c
c     find smallest positive root
c
      d=(b*b-4.0d0*a*c)/(4.0d0*a*a)
      if(d.gt.0.0d0) then
      r1=-b/(2.0d0*a)
      if(r1.ge.0) then
      r1=r1+sqrt(d)
      else
      r1=r1-sqrt(d)
      endif
      r2=c/(r1*a)
      r=-1.0d0
      if(r1.gt.0.0d0) r=r1
      if(r2.gt.0.0d0.and.r2.lt.r1) r=r2
      if(r.gt.1.0d0) r=-1.0d0
c
      if(r.gt.0.0d0) then
      d=r*sqrt(tr(1)*tr(1)+tr(2)*tr(2)+tr(3)*tr(3))
      if(d.lt.dt) then
      dt=d
      rt=r
      endif
      endif
c
      endif
      endif
      endif
    9 continue
      endif
    8 continue
c
      endif
c
      if(rt.gt.0.and.dt.gt.dtran) then
      large=m
      root=rt
      tran(1)=tr(1)
      tran(2)=tr(2)
      tran(3)=tr(3)
      dtran=dt
cx      write(*,'(a,i4,a,f12.6,a,f12.6)') 
cx     + '             Molecule ',large,' by ',dtran,' nm',root
      endif
c
    3 continue
c
      write(*,'(a,i4,a,f12.6,a,f12.6)') 
     + ' Translating molecule ',large,' by ',dtran,' nm',root
c
c     process largest possible translation
c
      if(large.gt.0) then
      do 11 i=1,nsa
      if(id(i).eq.large) then
      do 12 k=1,3
      xs(k,i)=xs(k,i)+root*tran(k)
   12 continue
      endif
   11 continue
      if(dtran.gt.touch) goto 2
      endif
c
      return
      end

