      subroutine pre_collaps(xs,id,ndx,msa,nsa,touch,mode,nmoves)
c
      implicit none
c
      integer msa,nsa,mode,nmoves
      integer id(msa),ndx(msa,2)
      real*8 xs(3,msa)
      real*8 touch
c
      logical first,found
      integer i,j,k,m
      integer nmol,large,nummov
      real*8 rt,rk,dist
      real*8 a,b,c,d,r,r1,r2,root,dt,dtran,tr(3),tran(3),xmax(3),xmin(3)
c
      write(*,'(/,a,i4,f12.6,/)') ' Collapsing using mode ',mode,touch
c
c     find total number of molecules
c
c     on entry, id(i) contains the molecule number of atom i
c
      do 13 i=1,nsa
      ndx(i,1)=0
      ndx(i,2)=0
   13 continue
c
c     store in ndx(i,1) the number of atoms in molecule i
c
      nmol=0
      do 1 i=1,nsa
      nmol=max(nmol,id(i))
      ndx(id(i),1)=ndx(id(i),1)+1
    1 continue
c
      do 14 i=1,nsa
      ndx(i,2)=id(i)
   14 continue
c
c     group single atom molecules with nearest multi-atom molecules
c
      do 15 i=1,nsa
      if(ndx(id(i),1).eq.1) then
      k=0
      rk=1.0d10
      do 16 j=1,nsa
      if(ndx(id(j),1).gt.1) then
      dist=(xs(1,i)-xs(1,j))**2+(xs(2,i)-xs(2,j))**2+
     + (xs(3,i)-xs(3,j))**2
      if(dist.lt.rk) then
      k=j
      rk=dist
      endif
      endif
   16 continue
      if(k.gt.0) then
      ndx(id(i),1)=0
      ndx(i,2)=id(k)
      write(*,'(a,i8,a,i4,a,f12.6)') 
     + ' Single ion ',i,' grouped with molecule ',
     + id(k),', distance is ',sqrt(rk)
      endif
      endif
   15 continue
c
      nummov=0
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
      if(ndx(m,1).gt.0) then
c
c     determine translation vector
c
      rt=-1.0d0
      found=.false.
      first=.true.
      do 4 i=1,nsa
      if(ndx(i,2).eq.m) then
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
      if(ndx(i,2).eq.m) then
      do 9 j=1,nsa
      if(ndx(j,2).ne.m) then
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
      if(ndx(i,2).eq.large) then
      do 12 k=1,3
      xs(k,i)=xs(k,i)+root*tran(k)
   12 continue
      endif
   11 continue
      nummov=nummov+1
      if(dtran.gt.touch.and.nummov.lt.nmoves) goto 2
      endif
c
      return
      end

