      logical function pre_misfit(lfnout,isfnd,xs,csa,qsa,msa,nsa,
     + idsb,cdsb,msb,nsb,grid,mgrid,ngrid,gdist,iwater,iwfnd,xw,
     + qwa,mwm,mwa,nwm,nwa,npbtyp,box)
c
c $Id: pre_misfit.F,v 1.19 2000-03-18 00:08:53 d3j191 Exp $
c
      implicit none
c
#include "util.fh"
c
#if defined(CRAY_T3D) || defined(CRAY_T3E)
      integer arg_nint
      external arg_nint
#endif
c
      integer lfnout,msa,nsa,msb,nsb,mgrid,iwater,mwm,mwa,nwm,nwa
      integer isfnd(msa),idsb(2,msb),iwfnd(mwm),npbtyp
      real*8 cdsb(msb),xs(3,msa),qsa(msa)
      character*16 csa(msa)
      real*8 xw(3,mwa,mwm),qwa(mwa),box(3)
c
      real*8 grid(4,mgrid),distx,disty,distz
      integer i,ic,j,nn,nndx(10),nnn,nnndx(10),nm,miss
      real*8 d(3),c(3),p(3),t(3),r,dr,dd,px(3),pn(3),pz(3),pw(3)
      real*8 xmax(3),xmin(3),dx(3),dist2,gdist,gdist2,angle,boxh(3)
      integer ix,iy,iz,k,imax,imin,ngrid,iwm,isa,iwm2,iwa,num
c
      gdist2=gdist*gdist
      boxh(1)=0.5d0*box(1)
      boxh(2)=0.5d0*box(2)
      boxh(3)=0.5d0*box(3)
c
c     loop over all solute atoms
c     --------------------------
c
      ngrid=0      
      miss=0
c
      do 1 i=1,nsa
c
c     for each atom without coordinates
c     ---------------------------------
c
      if(isfnd(i).ne.1) then
c
      if(util_print('coordinates',print_debug)) then
      write(lfnout,'(a,i7,5x,a16)') 'Not found: atom ',i,csa(i)
      endif
      nn=0
c
c     find all bonded atoms
c     ---------------------
c
      do 2 j=1,nsb
      if(idsb(1,j).eq.i) then
      nn=nn+1
      nndx(nn)=idsb(2,j)
      endif
      if(idsb(2,j).eq.i) then
      nn=nn+1
      nndx(nn)=idsb(1,j)
      endif
    2 continue
c
      if(nn.eq.0) goto 1
c
c     if bound to single other existing atom
c
      if(nn.eq.1.and.isfnd(nndx(1)).eq.1) then
c
c     find atoms bound to that atom
c
      nnn=0
      do 3 j=1,nsb
      if(idsb(1,j).eq.nndx(1)) then
      nnn=nnn+1
      nnndx(nnn)=idsb(2,j)
      endif
      if(idsb(2,j).eq.nndx(1)) then
      nnn=nnn+1
      nnndx(nnn)=idsb(1,j)
      endif
      if(idsb(1,j).eq.i.or.idsb(2,j).eq.i) dr=cdsb(j)
    3 continue
c
c     count the number of missing atoms
c
      nm=0
      do 4 j=1,nnn
      if(isfnd(nnndx(j)).ne.1) nm=nm+1
    4 continue
c
c     if j is only missing atom
c
      if(nm.eq.1) then
      c(1)=xs(1,nndx(1))
      c(2)=xs(2,nndx(1))
      c(3)=xs(3,nndx(1))
      d(1)=0.0d0
      d(2)=0.0d0
      d(3)=0.0d0
      num=0
      do 5 j=1,nnn
      if(nnndx(j).ne.i) then
      d(1)=d(1)+xs(1,nnndx(j))
      d(2)=d(2)+xs(2,nnndx(j))
      d(3)=d(3)+xs(3,nnndx(j))
      num=num+1
      endif
    5 continue
      d(1)=d(1)/dble(num)-c(1)
      d(2)=d(2)/dble(num)-c(2)
      d(3)=d(3)/dble(num)-c(3)
      r=sqrt(d(1)*d(1)+d(2)*d(2)+d(3)*d(3))
      if(r.lt.0.0001) call errquit('pre_misfit: vector too small',9001)
      xs(1,i)=c(1)-dr*d(1)/r
      xs(2,i)=c(2)-dr*d(2)/r
      xs(3,i)=c(3)-dr*d(3)/r
      isfnd(i)=1
      endif
c
c     if j is one of two missing atoms
c
      if(nm.eq.2) then
      c(1)=xs(1,nndx(1))
      c(2)=xs(2,nndx(1))
      c(3)=xs(3,nndx(1))
      d(1)=0.0d0
      d(2)=0.0d0
      d(3)=0.0d0
      ic=0
      num=0
      do 6 j=1,nnn
      if(nnndx(j).ne.i.and.isfnd(nnndx(j)).eq.1) then
      d(1)=d(1)+xs(1,nnndx(j))
      d(2)=d(2)+xs(2,nnndx(j))
      d(3)=d(3)+xs(3,nnndx(j))
      ic=nnndx(j)
      num=num+1
      endif
    6 continue
      d(1)=d(1)/dble(num)-c(1)
      d(2)=d(2)/dble(num)-c(2)
      d(3)=d(3)/dble(num)-c(3)
      r=sqrt(d(1)*d(1)+d(2)*d(2)+d(3)*d(3))
      if(r.lt.0.0001) call errquit('pre_misfit: vector too small',9002)
      d(1)=dr*d(1)/r
      d(2)=dr*d(2)/r
      d(3)=dr*d(3)/r
      if(nnn.eq.3) then
      if(abs(d(1)).le.abs(d(2)).and.abs(d(1)).le.abs(d(3))) then
      t(1)=c(1)
      t(2)=d(3)+c(2)
      t(3)=-d(2)+c(3)
      elseif(abs(d(2)).le.abs(d(1)).and.abs(d(2)).le.abs(d(3))) then
      t(1)=d(3)+c(1)
      t(2)=c(2)
      t(3)=-d(1)+c(3)
      else
      t(1)=d(2)+c(1)
      t(2)=-d(1)+c(2)
      t(3)=c(3)
      endif
      d(1)=d(1)+c(1)
      d(2)=d(2)+c(2)
      d(3)=d(3)+c(3)
      angle=2.094395102d0
      call rotate(c,t,angle,d,p)
      xs(1,i)=p(1)
      xs(2,i)=p(2)
      xs(3,i)=p(3)
      isfnd(i)=1
      else
      d(1)=d(1)+c(1)
      d(2)=d(2)+c(2)
      d(3)=d(3)+c(3)
      r=sqrt((xs(1,ic)-c(1))**2+(xs(2,ic)-c(2))**2+(xs(3,ic)-c(3))**2)
      if(r.lt.0.0001) call errquit('pre_misfit: vector too small',9004)
      t(1)=dr*(xs(1,ic)-c(1))/r+c(1)
      t(2)=dr*(xs(2,ic)-c(2))/r+c(2)
      t(3)=dr*(xs(3,ic)-c(3))/r+c(3)
      angle=1.570796327d0
      call rotate(c,d,angle,t,p)
      d(1)=d(1)-c(1)
      d(2)=d(2)-c(2)
      d(3)=d(3)-c(3)
      t(1)=t(1)-c(1)
      t(2)=t(2)-c(2)
      t(3)=t(3)-c(3)
      c(1)=d(2)*t(3)-t(2)*d(3)+xs(1,nndx(1))
      c(2)=d(3)*t(1)-t(3)*d(1)+xs(2,nndx(1))
      c(3)=d(1)*t(2)-t(1)*d(2)+xs(3,nndx(1))
      d(1)=xs(1,nndx(1))
      d(2)=xs(2,nndx(1))
      d(3)=xs(3,nndx(1))
      angle=3.141592654d0
      call rotate(d,c,angle,p,t)
      xs(1,i)=t(1)
      xs(2,i)=t(2)
      xs(3,i)=t(3)
      isfnd(i)=1
      endif
      endif
c
c     if j is one of three missing atoms
c
      if(nm.eq.3) then
      c(1)=xs(1,nndx(1))
      c(2)=xs(2,nndx(1))
      c(3)=xs(3,nndx(1))
      d(1)=0.0d0
      d(2)=0.0d0
      d(3)=0.0d0
      ic=0
      num=0
      do 7 j=1,nnn
      if(nnndx(j).ne.i.and.isfnd(nnndx(j)).eq.1) then
      d(1)=d(1)+xs(1,nnndx(j))
      d(2)=d(2)+xs(2,nnndx(j))
      d(3)=d(3)+xs(3,nnndx(j))
      ic=nnndx(j)
      num=num+1
      endif
    7 continue
      d(1)=d(1)/dble(num)-c(1)
      d(2)=d(2)/dble(num)-c(2)
      d(3)=d(3)/dble(num)-c(3)
      r=sqrt(d(1)*d(1)+d(2)*d(2)+d(3)*d(3))
      if(r.lt.0.0001) call errquit('pre_misfit: vector too small',9005)
      d(1)=dr*d(1)/r
      d(2)=dr*d(2)/r
      d(3)=dr*d(3)/r
      if(nnn.eq.4) then
      if(abs(d(1)).le.abs(d(2)).and.abs(d(1)).le.abs(d(3))) then
      t(1)=c(1)
      t(2)=d(3)+c(2)
      t(3)=-d(2)+c(3)
      elseif(abs(d(2)).le.abs(d(1)).and.abs(d(2)).le.abs(d(3))) then
      t(1)=d(3)+c(1)
      t(2)=c(2)
      t(3)=-d(1)+c(3)
      else
      t(1)=d(2)+c(1)
      t(2)=-d(1)+c(2)
      t(3)=c(3)
      endif
      d(1)=d(1)+c(1)
      d(2)=d(2)+c(2)
      d(3)=d(3)+c(3)
      angle=1.910633236d0
      call rotate(c,t,angle,d,p)
      xs(1,i)=p(1)
      xs(2,i)=p(2)
      xs(3,i)=p(3)
      isfnd(i)=1
      else
c
c     three missing hydrogens from atom with other than 4 neighbors
c
      call errquit('pre_misfit: case not implemented involving atom',i)
      endif
      endif
c
      endif
      endif
c
    1 continue
c
c
      do 8 i=1,nsa
c
c     for each atom without coordinates
c     ---------------------------------
c
      if(isfnd(i).ne.1) then
      nn=0
c
c     find all bonded atoms
c     ---------------------
c
      do 9 j=1,nsb
      if(idsb(1,j).eq.i) then
      nn=nn+1
      nndx(nn)=idsb(2,j)
      endif
      if(idsb(2,j).eq.i) then
      nn=nn+1
      nndx(nn)=idsb(1,j)
      endif
    9 continue
c
c     if not bound to any atom
c     ------------------------
c
      if(nn.eq.0) then
      if(util_print('coordinates',print_debug)) then
      write(lfnout,'(a,i7)') 'Setting up grid for atom ',i
      endif
c
c     setup grid
c
      if(ngrid.eq.0) then
      do 10 j=1,nsa
      if(isfnd(j).eq.1) then
      xmax(1)=xs(1,j)
      xmax(2)=xs(2,j)
      xmax(3)=xs(3,j)
      xmin(1)=xs(1,j)
      xmin(2)=xs(2,j)
      xmin(3)=xs(3,j)
      goto 11
      endif
   10 continue
   11 continue
      do 12 k=1,3
      do 13 j=1,nsa
      if(isfnd(j).eq.1) then
      if(xmax(k).lt.xs(k,j)) xmax(k)=xs(k,j)
      if(xmin(k).gt.xs(k,j)) xmin(k)=xs(k,j)
      endif
   13 continue
      dx(k)=(xmax(k)-xmin(k)+4.0d0*gdist)/(mgrid-1)
   12 continue
c
      ngrid=1
      do 14 ix=1,mgrid
      do 15 iy=1,mgrid
      do 16 iz=1,mgrid
      grid(1,ngrid)=xmin(1)+dble(ix-1)*dx(1)
      grid(2,ngrid)=xmin(2)+dble(iy-1)*dx(2)
      grid(3,ngrid)=xmin(3)+dble(iz-1)*dx(3)
      do 17 k=1,nsa
      if(isfnd(k).eq.1) then
      if(npbtyp.eq.0) then
      dist2=(grid(1,ngrid)-xs(1,k))**2+(grid(2,ngrid)-xs(2,k))**2+
     + (grid(3,ngrid)-xs(3,k))**2
      else
      distx=grid(1,ngrid)-xs(1,k)
      disty=grid(2,ngrid)-xs(2,k)
      distz=grid(3,ngrid)-xs(3,k)
#if defined(CRAY_T3D) || defined(CRAY_T3E)
      if(abs(distx).gt.boxh(1))
     + distx=distx-arg_nint(distx/box(1))*box(1)
      if(abs(disty).gt.boxh(2))
     + disty=disty-arg_nint(distx/box(2))*box(2)
      if(abs(distz).gt.boxh(3))
     + distz=distz-arg_nint(distx/box(3))*box(3)
#else
      if(abs(distx).gt.boxh(1)) distx=distx-nint(distx/box(1))*box(1)
      if(abs(disty).gt.boxh(2)) disty=disty-nint(distx/box(2))*box(2)
      if(abs(distz).gt.boxh(3)) distz=distz-nint(distx/box(3))*box(3)
#endif
      dist2=distx**2+disty**2+distz**2
      endif
      if(dist2.lt.gdist2) goto 16
      endif
   17 continue
      ngrid=ngrid+1
   16 continue
   15 continue
   14 continue
      ngrid=ngrid-1
c
      if(ngrid.le.0) call errquit('No grid points',9999)
c
c     calculate grid potential
c
      do 18 k=1,ngrid
      grid(4,k)=0.0d0
      do 19 j=1,nsa
      if(isfnd(j).eq.1) then
      if(npbtyp.eq.0) then
      dist2=dsqrt((grid(1,k)-xs(1,j))**2+
     + (grid(2,k)-xs(2,j))**2+(grid(3,k)-xs(3,j))**2)
      else
      distx=grid(1,ngrid)-xs(1,k)
      disty=grid(2,ngrid)-xs(2,k)
      distz=grid(3,ngrid)-xs(3,k)
#if defined(CRAY_T3D) || defined(CRAY_T3E)
      if(abs(distx).gt.boxh(1))
     + distx=distx-arg_nint(distx/box(1))*box(1)
      if(abs(disty).gt.boxh(2))
     + disty=disty-arg_nint(distx/box(2))*box(2)
      if(abs(distz).gt.boxh(3))
     + distz=distz-arg_nint(distx/box(3))*box(3)
#else
      if(abs(distx).gt.boxh(1)) distx=distx-nint(distx/box(1))*box(1)
      if(abs(disty).gt.boxh(2)) disty=disty-nint(distx/box(2))*box(2)
      if(abs(distz).gt.boxh(3)) distz=distz-nint(distx/box(3))*box(3)
#endif
      dist2=distx**2+disty**2+distz**2
      endif
      grid(4,k)=grid(4,k)+qsa(j)/dist2
      endif
   19 continue
   18 continue
c
      endif
c
      if(ngrid.le.0) call errquit('No grid points',9999)
c
c     find extrema
c
      imax=1
      imin=1
      do 20 j=1,ngrid
      if(grid(4,j).gt.grid(4,imax)) imax=j
      if(grid(4,j).lt.grid(4,imin)) imin=j
   20 continue
c
c     choose best grid point
c
      if(qsa(i).gt.0.0d0) then
      xs(1,i)=grid(1,imin)
      xs(2,i)=grid(2,imin)
      xs(3,i)=grid(3,imin)
      else
      xs(1,i)=grid(1,imax)
      xs(2,i)=grid(2,imax)
      xs(3,i)=grid(3,imax)
      imin=imax
      endif
      isfnd(i)=1
      write(6,6000) qsa(i)*grid(4,imin)*138.9354d0
 6000 format(' Added counterion with energy ',f12.3,' kJ/mol')
c
c     remove grid point from list
c
      do 21 j=imin,ngrid
      do 22 k=1,4
      grid(k,j)=grid(k,j+1)
   22 continue
   21 continue
      ngrid=ngrid-1
c
      if(ngrid.le.0) call errquit('No grid points',9999)
c
c     update grid potential
c
      do 27 k=1,ngrid
      dist2=dsqrt((grid(1,k)-xs(1,i))**2+
     + (grid(2,k)-xs(2,i))**2+(grid(3,k)-xs(3,i))**2)
      grid(4,k)=grid(4,k)+qsa(i)/dist2
   27 continue
c
      else
      miss=miss+1
      if(util_print('coordinates',print_none)) then
      write(lfnout,9998) i,csa(i)
 9998 format(' Coordinates missing for atom ',i5,':',a)
      endif
      endif
      endif
    8 continue
c
      if(miss.ne.0) then
      call errquit('pre_misfit: missing non hydrogen coordinates',9999)
      endif
c
c     crystal waters
c
      if(iwater.gt.0) then
      do 23 iwm=1,nwm
      if(iwfnd(iwm).eq.1) then
      px(1)=0.0d0
      px(2)=0.0d0
      px(3)=0.0d0
      do 24 isa=1,nsa
      dx(1)=(xs(1,isa)-xw(1,1,iwm))
      dx(2)=(xs(2,isa)-xw(2,1,iwm))
      dx(3)=(xs(3,isa)-xw(3,1,iwm))
      dd=sqrt(dx(1)*dx(1)+dx(2)*dx(2)+dx(3)*dx(3))
      px(1)=px(1)+qsa(isa)*dx(1)/dd
      px(2)=px(2)+qsa(isa)*dx(2)/dd
      px(3)=px(3)+qsa(isa)*dx(3)/dd
   24 continue
      do 25 iwm2=1,nwm
      if(iwfnd(iwm).eq.3) then
      do 26 iwa=1,3
      dx(1)=(xw(1,iwa,iwm2)-xw(1,1,iwm))
      dx(2)=(xw(2,iwa,iwm2)-xw(2,1,iwm))
      dx(3)=(xw(3,iwa,iwm2)-xw(3,1,iwm))
      dd=sqrt(dx(1)*dx(1)+dx(2)*dx(2)+dx(3)*dx(3))
      px(1)=px(1)+qwa(iwa)*dx(1)/dd
      px(2)=px(2)+qwa(iwa)*dx(2)/dd
      px(3)=px(3)+qwa(iwa)*dx(3)/dd
   26 continue
      endif
   25 continue
      pw(1)=0.0d0
      pw(2)=0.0d0
      pw(3)=0.055d0
      angle=acos((pw(1)*px(1)+pw(2)*px(2)+pw(3)*px(3))/
     + sqrt((pw(1)*pw(1)+pw(2)*pw(2)+pw(3)*pw(3))*
     + (px(1)*px(1)+px(2)*px(2)+px(3)*px(3))))
      pn(1)=pw(2)*px(3)-pw(3)*px(2)
      pn(2)=pw(3)*px(1)-pw(1)*px(3)
      pn(3)=pw(1)*px(2)-pw(2)*px(1)
      pz(1)=0.0d0
      pz(2)=0.0d0
      pz(3)=0.0d0
      pw(1)=0.0d0
      pw(2)=-0.083d0
      pw(3)=0.055d0
      call rotate(pz,pn,angle,pw,dx)
      xw(1,2,iwm)=xw(1,1,iwm)+dx(1)
      xw(2,2,iwm)=xw(2,1,iwm)+dx(2)
      xw(3,2,iwm)=xw(3,1,iwm)+dx(3)
      pw(1)=0.0d0
      pw(2)=0.083d0
      pw(3)=0.055d0
      call rotate(pz,pn,angle,pw,dx)
      xw(1,3,iwm)=xw(1,1,iwm)+dx(1)
      xw(2,3,iwm)=xw(2,1,iwm)+dx(2)
      xw(3,3,iwm)=xw(3,1,iwm)+dx(3)
      endif
   23 continue
      endif
c
      pre_misfit=.true.
      return
 9999 continue
      pre_misfit=.false.
      return
      end
