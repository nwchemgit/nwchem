      subroutine pre_collaps(xs,id,ndx,vec,msa,nsa,touch,mode,nmoves,
     + ncolgr,icolgr,lfnout,lfnpdb,lrgpdb,sysnam,filpdb,iopt,box,
     + num,amass,mat,
     + csa,isat,isgm,imol,ifra,vs,
     + cwa,iwat,xw,vw,mwm,mwa,nwm,nwa,
     + xwc,vwc,mwmc,nwmc,slvnam,iropt,irrand,nxrep,nyrep,nzrep,drep,
     + msb,nsb,idsb,rdist,nskip,iskip,lang,lfnmrg,nmerge,xmerge,filmrg,
     + irenum,invert,ihop,ips)
c
      implicit none
c
#include "mafdecls.fh"
      logical pre_wrtpdb
      external pre_wrtpdb
c
      integer msa,nsa,mode,nmoves
      integer id(msa),ndx(msa,3)
      real*8 xs(3,msa),vec(6,msa)
      real*8 touch
      integer ncolgr
      integer icolgr(2,ncolgr)
c
      integer lfnout,lfnpdb,iopt,mwm,nwm,mwa,nwa,mat,mwmc,nwmc
      integer lfnmrg,lrgpdb
      character*255 filpdb,filmrg(100)
      integer num(mat),isat(msa),isgm(msa),imol(msa),ifra(msa),iwat(mwa)
      integer ihop(msa),ips(msa)
      character*16 cwa(mwa),csa(msa)
      character*3 slvnam
      real*8 amass(mat),vs(3,msa),xw(3,mwa,mwm),vw(3,mwa,mwm)
      real*8 xwc(3,mwa,mwmc),vwc(3,mwa,mwmc),box(3)
      integer iropt,irrand,nxrep,nyrep,nzrep,nskip
      real*8 rdist,drep(3),xmerge(3,100)
      integer msb,nsb,nmols,irenum,invert,nmerge
      integer idsb(2,msb),iskip(3,nskip),lang(*)
      character*80 sysnam
c
      character*255 filnam
c
      logical first,found
      integer i,j,k,m
      integer nmol,large,nummov
      real*8 rt,rk,dist
      real*8 a,b,c,d,e,r,r1,r2,root,dt,dtran
      real*8 tr(3),tran(3),xmax(3),xmin(3)
      real*8 a2,b2,c2,d2,p(3),factor
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
      ndx(i,3)=0
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
      ndx(id(k),1)=ndx(id(k),1)+1
      ndx(i,2)=id(k)
      write(*,'(a,i8,a,i4,a,f12.6)') 
     + ' Single ion ',i,' grouped with molecule ',
     + id(k),', distance is ',sqrt(rk)
      endif
      endif
   15 continue
c
      if(ncolgr.gt.0) then
      do 26 j=1,ncolgr
      do 27 i=1,nsa
      if(ndx(i,2).eq.icolgr(2,j)) ndx(i,2)=icolgr(1,j)
   27 continue
      write(*,'(a,i8,a,i4)') 
     + ' Molecule   ',icolgr(2,j),
     + ' grouped with molecule ',icolgr(1,j)
   26 continue
      endif
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
c     determine centers of geometry
c
      do 31 m=1,nmol
      vec(1,m)=0.0d0
      vec(2,m)=0.0d0
      vec(3,m)=0.0d0
      vec(4,m)=0.0d0
      vec(5,m)=0.0d0
      vec(6,m)=0.0d0
   31 continue
      do 32 i=1,nsa
      m=ndx(i,2)
      vec(1,m)=xs(1,i)
      vec(2,m)=xs(2,i)
      vec(3,m)=xs(3,i)
      vec(4,m)=xs(1,i)
      vec(5,m)=xs(2,i)
      vec(6,m)=xs(3,i)
   32 continue
      do 33 i=1,nsa
      m=ndx(i,2)
      vec(1,m)=min(vec(1,m),xs(1,i))
      vec(2,m)=min(vec(2,m),xs(2,i))
      vec(3,m)=min(vec(3,m),xs(3,i))
      vec(4,m)=max(vec(4,m),xs(1,i))
      vec(5,m)=max(vec(5,m),xs(2,i))
      vec(6,m)=max(vec(6,m),xs(3,i))
   33 continue
      do 34 m=1,nmol
      if(ndx(m,1).gt.0) then
      vec(1,m)=0.5d0*(vec(1,m)+vec(4,m))
      vec(2,m)=0.5d0*(vec(2,m)+vec(5,m))
      vec(3,m)=0.5d0*(vec(3,m)+vec(6,m))
      vec(4,m)=vec(4,m)-vec(1,m)
      vec(5,m)=vec(5,m)-vec(2,m)
      vec(6,m)=vec(6,m)-vec(3,m)
      endif
   34 continue
c
      do 3 m=1,nmol
      if(ndx(m,1).gt.0) then
c
      tr(1)=-vec(1,m)
      tr(2)=-vec(2,m)
      tr(3)=-vec(3,m)
      if(mode.eq.3) tr(3)=0.0d0
c
      do 41 i=1,nmol
      ndx(i,3)=1
   41 continue
      ndx(m,3)=0
      b2=vec(1,m)*vec(1,m)+vec(2,m)*vec(2,m)+vec(3,m)*vec(3,m)
      d=sqrt(vec(4,m)*vec(4,m)+vec(5,m)*vec(5,m)+vec(6,m)*vec(6,m))
      e=sqrt(vec(4,m)*vec(4,m)+vec(5,m)*vec(5,m))
      do 42 i=1,nmol
      k=1
      if(ndx(i,1).eq.0.or.i.eq.m) then
      k=0
      else
      if(abs(vec(1,i)-(vec(1,m)+0.5d0*tr(1))).gt.
     +   abs(vec(1,m)+0.5d0*tr(1))+vec(4,i)+vec(4,m)) k=0
      if(abs(vec(2,i)-(vec(2,m)+0.5d0*tr(2))).gt.
     +   abs(vec(2,m)+0.5d0*tr(2))+vec(5,i)+vec(5,m)) k=0
      if(abs(vec(3,i)-(vec(3,m)+0.5d0*tr(3))).gt.
     +   abs(vec(3,m)+0.5d0*tr(3))+vec(6,i)+vec(6,m)) k=0
      if(k.eq.1) then
      a2=vec(1,i)*vec(1,i)+vec(2,i)*vec(2,i)+vec(3,i)*vec(3,i)
      c2=(vec(1,m)-vec(1,i))**2+(vec(2,m)-vec(2,i))**2+
     + (vec(3,m)-vec(3,i))**2
      factor=(1.0d0-0.5d0*(b2+c2-a2)/b2)
      p(1)=factor*vec(1,m)
      p(2)=factor*vec(2,m)
      p(3)=factor*vec(3,m)
      a=sqrt((p(1)-vec(1,i))**2+(p(2)-vec(2,i))**2+(p(3)-vec(3,i))**2)
      c=sqrt(vec(4,i)**2+vec(5,i)**2+vec(6,i)**2)
      if(a.gt.c+d) k=0
      if(mode.eq.3.and. k.eq.1) then
      if(abs(p(3)-vec(3,i)).gt.abs(vec(6,i))+abs(vec(6,m))) k=0
      endif
      endif
      endif
      if(k.eq.0) ndx(i,3)=0
   42 continue
c
c     determine translation vector
c
      rt=-1.0d0
c
      dt=sqrt(tr(1)*tr(1)+tr(2)*tr(2)+tr(3)*tr(3))
c
      do 8 i=1,nsa
      if(ndx(i,2).eq.m) then
      do 9 j=1,nsa
      if(ndx(ndx(j,2),3).eq.1) then
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
      if(rt.gt.0.and.dt.gt.dtran) then
      large=m
      root=rt
      tran(1)=tr(1)
      tran(2)=tr(2)
      tran(3)=tr(3)
      dtran=dt
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
      write(filnam,'(a,i5.5,a)')
     + sysnam(1:index(sysnam,' ')-1)//'-c',nummov,'.pdb'
c
      if(.not.pre_wrtpdb(lfnout,lfnpdb,lrgpdb,filnam,iopt,box,num,amass,
     + mat,csa,isat,isgm,imol,ifra,xs,vs,msa,nsa,cwa,iwat,xw,vw,mwm,mwa,
     + nwm,nwa,xwc,vwc,mwmc,nwmc,slvnam,iropt,irrand,nxrep,nyrep,nzrep,
     + drep,msb,nsb,idsb,rdist,nskip,iskip,lang,lfnmrg,nmerge,xmerge,
     + filmrg,irenum,invert,ihop,ips))
     +  call md_abort('pre_wrtpdb failed',9999)
c
      if(dtran.gt.touch.and.(nummov.lt.nmoves.or.nmoves.lt.0)) goto 2
      endif
c
      return
      end

c $Id$
