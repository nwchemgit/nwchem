      logical function argos_prep_hybrid(xatm,latm,matm,natm,
     + lbnd,mbnd,nbnd)
c
c $Id$
c
c     in  : xatm(3,matm) = atomic coordinates
c         : latm(2,matm) = atomic numbers
c         : matm         = dimension of coordinate list
c         : natm         = number of atoms
c
      implicit none
c
      real*8 angl
      external angl
c
      integer matm,natm
      real*8 xatm(3,matm)
      integer latm(5,matm)
      integer mbnd,nbnd
      integer lbnd(2,mbnd)
c
      integer i,j,k,l,m,nn,nl(100),ia,ib,ic,id,ie,ig
      real*8 angle1,angle2,angle3,angle4,angle5,angle6
c
      do 1 i=1,natm
      if(latm(2,i).gt.0) then
      nn=0
      do 2 j=1,natm
      if(latm(2,j).gt.0.and.i.ne.j) then
      if(i.eq.j) goto 2
      do 6 k=1,nbnd
      if((lbnd(1,k).eq.i.and.lbnd(2,k).eq.j).or.
     + (lbnd(1,k).eq.j.and.lbnd(2,k).eq.i)) then
      nn=nn+1
      nl(nn)=j
      goto 2
      endif
    6 continue
      endif
    2 continue
      latm(3,i)=0
      if(nn.eq.1) then
      if(latm(2,i).eq.7) latm(3,i)=4
      elseif(nn.eq.2) then
      angle1=angl(xatm(1,nl(1)),xatm(1,i),xatm(1,nl(2)))
      ia=0
      if(angle1.ge.1.48d0.and.angle1.lt.1.66d0) ia=6
      if(angle1.ge.1.83d0.and.angle1.lt.2.01d0) ia=4
      if(angle1.ge.2.01d0.and.angle1.lt.2.18d0) ia=3
      if(angle1.ge.2.96d0.and.angle1.lt.3.15d0) ia=2
      latm(3,i)=ia
      if(latm(2,i).eq.7) latm(3,i)=3
      if(latm(2,i).eq.8) latm(3,i)=2
      elseif(nn.eq.3) then
      angle1=angl(xatm(1,nl(1)),xatm(1,i),xatm(1,nl(2)))
      ia=0
      if(angle1.ge.1.48d0.and.angle1.lt.1.66d0) ia=6
      if(angle1.ge.1.83d0.and.angle1.lt.2.01d0) ia=4
      if(angle1.ge.2.01d0.and.angle1.lt.2.18d0) ia=3
      if(angle1.ge.2.96d0.and.angle1.lt.3.15d0) ia=2
      angle2=angl(xatm(1,nl(1)),xatm(1,i),xatm(1,nl(3)))
      ib=0
      if(angle2.ge.1.48d0.and.angle2.lt.1.66d0) ib=6
      if(angle2.ge.1.83d0.and.angle2.lt.2.01d0) ib=4
      if(angle2.ge.2.01d0.and.angle2.lt.2.18d0) ib=3
      if(angle2.ge.2.96d0.and.angle2.lt.3.15d0) ib=2
      angle3=angl(xatm(1,nl(2)),xatm(1,i),xatm(1,nl(3)))
      ic=0
      if(angle3.ge.1.48d0.and.angle3.lt.1.66d0) ic=6
      if(angle3.ge.1.83d0.and.angle3.lt.2.01d0) ic=4
      if(angle3.ge.2.01d0.and.angle3.lt.2.18d0) ic=3
      if(angle3.ge.2.96d0.and.angle3.lt.3.15d0) ic=2
      if(ia.eq.ib.and.ib.eq.ic) latm(3,i)=ia
      latm(3,i)=4
      if(angle1+angle2+angle3.gt.6.0) latm(3,i)=3
      elseif(nn.eq.4) then
      angle1=angl(xatm(1,nl(1)),xatm(1,i),xatm(1,nl(2)))
      ia=0
      if(angle1.ge.1.48d0.and.angle1.lt.1.66d0) ia=6
      if(angle1.ge.1.83d0.and.angle1.lt.2.01d0) ia=4
      if(angle1.ge.2.01d0.and.angle1.lt.2.18d0) ia=3
      if(angle1.ge.2.96d0.and.angle1.lt.3.15d0) ia=2
      angle2=angl(xatm(1,nl(1)),xatm(1,i),xatm(1,nl(3)))
      ib=0
      if(angle2.ge.1.48d0.and.angle2.lt.1.66d0) ib=6
      if(angle2.ge.1.83d0.and.angle2.lt.2.01d0) ib=4
      if(angle2.ge.2.01d0.and.angle2.lt.2.18d0) ib=3
      if(angle2.ge.2.96d0.and.angle2.lt.3.15d0) ib=2
      angle3=angl(xatm(1,nl(1)),xatm(1,i),xatm(1,nl(4)))
      ic=0
      if(angle3.ge.1.48d0.and.angle3.lt.1.66d0) ic=6
      if(angle3.ge.1.83d0.and.angle3.lt.2.01d0) ic=4
      if(angle3.ge.2.01d0.and.angle3.lt.2.18d0) ic=3
      if(angle3.ge.2.96d0.and.angle3.lt.3.15d0) ic=2
      angle4=angl(xatm(1,nl(2)),xatm(1,i),xatm(1,nl(3)))
      id=0
      if(angle4.ge.1.48d0.and.angle4.lt.1.66d0) id=6
      if(angle4.ge.1.83d0.and.angle4.lt.2.01d0) id=4
      if(angle4.ge.2.01d0.and.angle4.lt.2.18d0) id=3
      if(angle4.ge.2.96d0.and.angle4.lt.3.15d0) id=2
      angle5=angl(xatm(1,nl(2)),xatm(1,i),xatm(1,nl(4)))
      ie=0
      if(angle5.ge.1.48d0.and.angle5.lt.1.66d0) ie=6
      if(angle5.ge.1.83d0.and.angle5.lt.2.01d0) ie=4
      if(angle5.ge.2.01d0.and.angle5.lt.2.18d0) ie=3
      if(angle5.ge.2.96d0.and.angle5.lt.3.15d0) ie=2
      angle6=angl(xatm(1,nl(3)),xatm(1,i),xatm(1,nl(4)))
      ig=0
      if(angle6.ge.1.48d0.and.angle6.lt.1.66d0) ig=6
      if(angle6.ge.1.83d0.and.angle6.lt.2.01d0) ig=4
      if(angle6.ge.2.01d0.and.angle6.lt.2.18d0) ig=3
      if(angle6.ge.2.96d0.and.angle6.lt.3.15d0) ig=2
      if(ia.eq.ib.and.ia.eq.ic.and.ia.eq.id.and.ia.eq.ie.and.ia.eq.ig)
     + latm(3,i)=ia
      if(latm(2,i).eq.6) latm(3,i)=4
      endif
      endif
    1 continue
c
c     sanity checks
c     -------------
c
      do 3 i=1,natm
      nn=0
      if(latm(2,i).gt.0) then
      do 4 j=1,natm
      if(latm(2,j).gt.0.and.i.ne.j) then
      if(i.eq.j) goto 4
      do 7 k=1,nbnd
      if((lbnd(1,k).eq.i.and.lbnd(2,k).eq.j).or.
     + (lbnd(1,k).eq.j.and.lbnd(2,k).eq.i)) then
      nn=nn+1
      nl(nn)=j
      goto 4
      endif
    7 continue
      endif
    4 continue
      endif
c
c     for following checks:
c
c     latm(2,i) : atomic number of atom i
c     latm(3,i) : hybridization of atom i: number of neighboring atoms
c     nn        : number of non-hydrogen neighbors of atom i
c     nl(1:nn)  : id's of atoms nn non-hydrogen neighbors of atom i
c
c     c attached to peptide c should be methyl group
c
      if(latm(2,i).eq.6.and.nn.eq.3) then
      l=0
      m=0
      do 5 k=1,3
      if(latm(2,nl(k)).eq.8.and.latm(3,nl(k)).eq.1) l=l+1
      if(latm(2,nl(k)).eq.7.and.latm(3,nl(k)).eq.3) l=l+1
      if(latm(2,nl(k)).eq.6.and.latm(3,nl(k)).eq.0) then
      l=l+1
      m=k
      endif
    5 continue
      if(l.eq.3.and.m.ne.0) latm(3,nl(m))=4
      endif
c
c     c between saturated atoms must be saturated
c
      if(latm(2,i).eq.6.and.nn.eq.2) then
      if(latm(2,nl(1)).eq.6.and.latm(3,nl(1)).eq.4.and.
     + latm(2,nl(2)).eq.6.and.latm(3,nl(2)).eq.4) latm(3,i)=4
      endif
c
c     NH2-C(sp2)
c
      if(latm(2,i).eq.7.and.nn.eq.1) then
      if(latm(2,nl(1)).eq.6.and.latm(3,nl(1)).eq.3) latm(3,i)=3
      endif
c
c     N(sp2)-CH-N(sp2)
c
      if(latm(2,i).eq.6.and.nn.eq.2) then
      if(latm(2,nl(1)).eq.7.and.latm(3,nl(1)).eq.3.and.
     + latm(2,nl(2)).eq.7.and.latm(3,nl(2)).eq.3) latm(3,i)=3
      endif
c
    3 continue
c
      argos_prep_hybrid=.true.
      return
      end
