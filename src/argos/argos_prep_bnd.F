      logical function argos_prep_bnd(xatm,latm,catm,matm,natm,ilo,ihi,
     + xato,lato,cato,mato,nato,lbnd,mbnd,nbnd,maxscf,
     + llnk,clnk,mlnk,nlnk,iconst,mang,lang,mdih,ldih,itmp)
c
c $Id$
c
c     function to append the list of bonds from atomic coordinates
c              for all atoms i in [ilo,ihi]
c
c     in  : xatm(3,matm) = atomic coordinates
c         : latm(2,matm) = atomic numbers
c         : matm         = dimension of coordinate list
c         : natm         = number of atoms
c         : mbnd         = dimension of the bond index list
c     out : lbnd(2,mbnd) = list of bond indices
c         : nbnd         = number of bonds
c         : mang         = dimension of the angle index list
c
c     function returns .false. if mbnd < nbnd
c
      implicit none
c
      real*8 argos_prep_atsiz
      logical argos_prep_ang,argos_prep_dih,argos_prep_ring
      external argos_prep_atsiz,argos_prep_ang,argos_prep_dih,
     + argos_prep_ring
c
      integer matm,mbnd,mato
      integer natm,nbnd,nato,maxscf,iconst
      real*8 xatm(3,matm),xato(3,mato)
      integer latm(5,matm),lato(5,mato),lbnd(2,mbnd),itmp(12,mato)
      character*6 catm(3,matm),cato(3,mato)
      integer ilo,ihi
      integer i,j,k,l,nn,nl(10),nadd,ndx,neigh
      integer mlnk,nlnk
      integer llnk(4,mlnk)
      character*4 clnk(2,mlnk)
      integer mang,mdih
      integer lang(3,mang),ldih(4,mdih)
      real*8 dist
      logical linker
      integer nang,ndih
c
      argos_prep_bnd=.false.
c
      nato=0
      do 1 i=ilo,ihi
      neigh=0
      if(latm(2,i).gt.0) then
      nn=0
      do 2 j=1,natm
      if(latm(2,j).gt.0.and.i.ne.j) then
      dist=sqrt((xatm(1,j)-xatm(1,i))*(xatm(1,j)-xatm(1,i))+
     + (xatm(2,j)-xatm(2,i))*(xatm(2,j)-xatm(2,i))+
     + (xatm(3,j)-xatm(3,i))*(xatm(3,j)-xatm(3,i)))
      linker=.false.
      do 23 k=1,nlnk
      if(i.eq.llnk(3,k).and.j.eq.llnk(4,k)) linker=.true.
      if(j.eq.llnk(3,k).and.i.eq.llnk(4,k)) linker=.true.
   23 continue
      if(linker.or.(dist.gt.0.05d0.and.
     + dist.le.argos_prep_atsiz(latm(2,j))+
     + argos_prep_atsiz(latm(2,i)))) then
      nn=nn+1
      nl(nn)=j
      if(j.gt.i.and.j.le.ihi) then
      nbnd=nbnd+1
      if(nbnd.gt.mbnd) call md_abort('Dimension problem',0)
      lbnd(1,nbnd)=(i-ilo+1)
      lbnd(2,nbnd)=(j-ilo+1)
      endif
      if(j.lt.ilo.or.j.gt.ihi) then
      neigh=j
      endif
      endif
      endif
    2 continue
      nato=nato+1
      do 3 j=1,5
      lato(j,nato)=latm(j,i)
    3 continue
      itmp(11,nato)=nn
      itmp(12,nato)=neigh
      do 40 j=1,nn
      itmp(j,nato)=nl(j)
   40 continue
      do 4 j=1,3
      cato(j,nato)=catm(j,i)
      xato(j,nato)=xatm(j,i)
    4 continue
      lato(1,nato)=0
      endif
    1 continue
c
      if(iconst.gt.0) then
c
      nang=0
      if(.not.argos_prep_ang(lbnd,mbnd,nbnd,lang,mang,nang))
     + call md_abort('argos_prep_ang failed',9999)
      ndih=0
      if(.not.argos_prep_dih(lang,mang,nang,ldih,mdih,ndih))
     + call md_abort('argos_prep_dih failed',9999)
      if(.not.argos_prep_ring(lato,mato,nato,lang,mang,nang,
     + ldih,mdih,ndih))
     + call md_abort('argos_prep_ring failed',9999)
c
      do 20 i=1,nato
      nadd=0
      nn=itmp(11,i)
      do 41 j=1,nn
      nl(j)=itmp(j,i)
   41 continue
      if(lato(2,i).eq.6) then
      if(nn.eq.0) then
      nadd=4
      elseif(nn.eq.1) then
      nadd=3
      if(lato(2,nl(1)).eq.6.and.lato(3,nl(1)).eq.0) nadd=3
      if(lato(3,i).eq.1) nadd=0
      if(lato(3,i).eq.2) nadd=1
      if(lato(3,i).eq.3) nadd=2
      if(lato(3,i).eq.4) nadd=3
      elseif(nn.eq.2) then
      nadd=2
      if(lato(3,i).eq.3) nadd=1
      elseif(nn.eq.3) then
      if(lato(3,i).eq.4) nadd=1
      endif
      endif
c
      if(lato(2,i).eq.7) then
      if(nn.eq.0) then
      nadd=4
      elseif(nn.eq.1) then
      nadd=3
      if(lato(3,i).eq.3) nadd=2
      elseif(nn.eq.2) then
      nadd=1
      if(lato(4,i).eq.6) nadd=0
      endif
      endif
c
      if(lato(2,i).eq.8) then
      if(nn.eq.0) then
      nadd=2
      elseif(nn.eq.1) then
      if(lato(3,nl(1)).eq.3) nadd=0
      if(lato(3,nl(1)).eq.4) nadd=1
      endif
      endif
c
      lato(1,i)=nadd
c
c     if atoms are added ther are no coordinates to do
c     an scf/esp calculation
c
      if(nadd.gt.0) maxscf=0
   20 continue
c
c     make corrections based on the aromitic ring information
c     -------------------------------------------------------
c
      do 21 i=1,nato
      nn=itmp(11,i)
      do 42 j=1,nn
      nl(j)=itmp(j,i)
   42 continue
      if(lato(2,i).eq.6) then
      if(nn.eq.2.and.lato(4,i).eq.5) then
      if(lato(2,nl(1)).eq.7.and.lato(4,nl(1)).eq.5.and.
     + lato(2,nl(2)).eq.7.and.lato(4,nl(2)).eq.5) then
      if(itmp(11,nl(1)).eq.3.and.itmp(11,nl(2)).eq.3)
     + call md_abort('N(3)-C(3)-N(3) in 5R',0)
      lato(1,nl(2))=0
      if(itmp(11,nl(1)).eq.2.and.itmp(11,nl(2)).eq.3) then
      lato(1,nl(1))=0
      elseif(itmp(11,nl(1)).eq.3.and.itmp(11,nl(2)).eq.2) then
      lato(1,nl(1))=0
      elseif(itmp(11,nl(1)).eq.2.and.itmp(11,nl(2)).eq.2) then
      lato(1,nl(1))=1
      endif
      endif
      endif
      endif
   21 continue
      endif
c
      do 30 i=1,nato
      lato(3,i)=itmp(12,i)
   30 continue
c
c     add hydrogens
c     -------------
c
      ndx=1
      do 6 i=1,nato
      do 7 j=1,nbnd
      if(lbnd(1,j).eq.i) lbnd(1,j)=-ndx
      if(lbnd(2,j).eq.i) lbnd(2,j)=-ndx
    7 continue
      nadd=lato(1,i)
      lato(1,i)=ndx
      ndx=ndx+nadd+1
    6 continue
      do 77 i=1,nbnd
      lbnd(1,i)=abs(lbnd(1,i))
      lbnd(2,i)=abs(lbnd(2,i))
   77 continue
c
      do 8 i=nato,1,-1
      if(lato(1,i).ne.i) then
      j=lato(1,i)
      lato(1,j)=lato(1,i)
      lato(2,j)=lato(2,i)
      lato(3,j)=lato(3,i)
      lato(4,j)=lato(4,i)
      lato(5,j)=lato(5,i)
      xato(1,j)=xato(1,i)
      xato(2,j)=xato(2,i)
      xato(3,j)=xato(3,i)
      cato(1,j)=cato(1,i)
      cato(2,j)=cato(2,i)
      cato(3,j)=cato(3,i)
      lato(1,i)=0
      endif
    8 continue
      nato=ndx-1
c
      j=1
      k=1
      do 9 i=1,nato
      if(lato(1,i).ne.0) then
      j=i
      k=1
      else
      lato(1,i)=0
      lato(2,i)=1
      lato(3,i)=0
      lato(4,i)=0
      lato(5,i)=0
      xato(1,i)=0.0d0
      xato(2,i)=0.0d0
      xato(3,i)=0.0d0
      cato(1,i)=' H'//cato(1,j)(3:5)
      cato(2,i)=' H'//cato(2,j)(3:5)
      cato(3,i)='     '
   11 continue
      k=k+1
      if(k.ge.10) call md_abort('Too many hydrogens',0)
      write(cato(1,i)(1:1),'(i1)') k
      write(cato(2,i)(1:1),'(i1)') k
      do 10 l=1,i-1
      if(cato(1,i).eq.cato(1,l)) goto 11
   10 continue
      nbnd=nbnd+1
      lbnd(1,nbnd)=j
      lbnd(2,nbnd)=i
      endif
    9 continue
c
      do 12 i=1,nbnd
      if(lbnd(1,i).gt.lbnd(2,i)) then
      j=lbnd(1,i)
      lbnd(1,i)=lbnd(2,i)
      lbnd(2,i)=j
      endif
   12 continue
c
      do 13 i=1,nbnd-1
      do 14 j=i+1,nbnd
      if(lbnd(1,i).gt.lbnd(1,j).or.
     + (lbnd(1,i).eq.lbnd(1,j).and.lbnd(2,i).gt.lbnd(2,j))) then
      k=lbnd(1,i)
      lbnd(1,i)=lbnd(1,j)
      lbnd(1,j)=k
      k=lbnd(2,i)
      lbnd(2,i)=lbnd(2,j)
      lbnd(2,j)=k
      endif
   14 continue
   13 continue
c
      argos_prep_bnd=.true.
      return
      end
