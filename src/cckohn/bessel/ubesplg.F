      subroutine ubesplg(rtdb,npd,lmax,nper,nint,alpha,
     &     xr,wwl,derl,csc,csl)
c
c     bessel fuction spline code with new prescription for
c     the outgoing wave continuum function.
c
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
      integer rtdb
      integer nbig, lbig, nsmall,nsmall2,nchmx,nemx,np,npd
      parameter (nbig=15000,lbig=6,nsmall=4015)
      parameter (nchmx=20,nemx=100)
      parameter (nsmall2=2*nsmall)
      integer ipow(0:lbig)
      double precision echan(nchmx),energy(nemx)
      double precision xstart,ystart,ypstart,al,mu,h,xmax,
     &     e,eta,ak,rd,r,w,cs,rd26
      integer l
      common /integ/ xstart,ystart,ypstart,al,mu,h,xmax
      common /energy/ e,eta,ak
      common/spl/rd,r(nsmall),w(nsmall,0:lbig),cs(nsmall),rd26,l
      complex*16 csc(nsmall,0:lbig),scc(nsmall)
      double precision aj(0:(lbig+1)),ay(0:(lbig+1))
      double precision aj1(0:(lbig+1)),ay1(0:(lbig+1))
      double precision csl(nsmall,0:lbig),derl(nsmall,0:lbig)
      complex*16 x1,x2,f1,f2,f,fp,gp,sig(0:lbig),tv
      complex*16 yc1,yc2,zl
c      double precision css(nsmall)
      complex*16 c1,c2,ai,ww(nsmall),wwl(nsmall,0:lbig)
      double precision xr(*)
      double precision zero
      integer gg
      external zero,gg
      double precision der(nsmall),xx(nbig),yy(nbig),zz(nbig),
     &     scr(nsmall),wwlr(nsmall2,0:lbig),cscr(nsmall2,0:lbig)
c      equivalence (wwl(1,0),wwlr(1,0)),(csc(1,0),cscr(1,0))
c
      integer i,lmax,nper,nint,nl,ncomp,j,index,i1,i2,np2,k
      double precision znuc,alpha,aint,x,yp1,yp2,det,c2r
c
      ai=cmplx(0.,1.)
c
c..unicos
c      call link ("unit5=infreeg,unit6=(outfree,hc,create)
c..unicos     1,unit8=bessplr//")
c
c      open(5,file='infree')
c      open(6,file='outfree')
c      open(8,file='bessplr',form='unformatted')
c
      ipow(0)=0
      do 66 i=1,lbig
         ipow(i)=1
 66   continue
c     call keep80(4hwave)
c     call frid(10hxerox+film,1,4,1)
c     call dders(-1)
      eta=0.
      znuc=0.
      mu=1.
c
c lmax=max. l value desired
c nper=no. of points per unit a.u. for splined functions
c xstart = first r-point (.00001 is o.k.)
c xmax = last r-point(r=k*x!)
c nint = multiple of nper for integration(3 or 4 is reasonable)
c alpha = exponential parameter in test function
c
c      read(5,*)lmax,nper,xstart,xmax,nint,alpha
      if(.not.rtdb_get(rtdb,'cckohn:lmax',mt_int,1,lmax))
     &     call errquit("can not get lmax from rtdb", 100, RTDB_ERR)
      if(.not.rtdb_get(rtdb,'cckohn:nper',mt_int,1,nper))
     &     call errquit("can not get nper from rtdb", 100, RTDB_ERR)
      if(.not.rtdb_get(rtdb,'cckohn:xstart',mt_dbl,1,xstart))
     &     call errquit("can not get xstart from rtdb", 100, RTDB_ERR)
      if(.not.rtdb_get(rtdb,'cckohn:xmax',mt_dbl,1,xmax))
     &     call errquit("can not get xmax from rtdb", 100, RTDB_ERR)
      if(.not.rtdb_get(rtdb,'cckohn:nint',mt_int,1,nint))
     &     call errquit("can not get nint from rtdb", 100, RTDB_ERR)
      if(.not.rtdb_get(rtdb,'cckohn:alpha',mt_dbl,1,alpha))
     &     call errquit("can not get alpha from rtdb", 100, RTDB_ERR)
      write(6,*)' Greens Function Code '
      write(6,101)lmax,nper,xstart,xmax,nint,alpha
 101  format(/,' L-Max  = ',i5,/,
     2     ' Nper   = ',i5,/,
     3     ' XStart = ',f20.12,/,
     4     ' XMax   = ',f20.12,/,
     5     ' Nint   = ',i5,/,
     6     ' Alpha  = ',f20.12)
      np=nper*xmax+1
      if(np.gt.nsmall.or.npd.ne.nsmall)then
         write(6,100)
 100     format(" np.gt.nsmall")
         stop
      endif
      rd=(xmax-xstart)/(np-1)
      rd26=rd*rd/6.
      do i=1,np
         r(i)=rd*(i-1.)+xstart
      enddo
      aint=nint
      h=rd/aint
      nl=lmax+1
      do 667 i=1,np
         x=r(i)
         call sjymec(x,aj,ay,ncomp,lmax)
         do 668 j=0,lmax
            w(i,j)=exp(-x*alpha)*aj(j)*x
 668     continue
 667  continue
      do 33 l=0,lmax
         al=l
         ystart=xstart**(l+1)
         ypstart=(al+1.)*xstart**l
         yp1=(w(2,l)-w(1,l))/rd
         yp2=(w(np,l)-w(np-1,l))/rd
         call cckohn_spline(r,w(1,l),np,yp1,yp2,scr,cs)
         call outward(xx,yy,index,zero)
         call outward(xx,zz,index,gg)
         i1=index
         i2=index-1./h
         x1=xx(i1)
         x2=xx(i2)
         det=yy(i1)*zz(i2)-yy(i2)*zz(i1)
c         write(*,*)"l,x1,x2,xx(i1),xx(i2),yy(i1),yy(i2),zz(i1),
c     &        zz(i2),det",
c     &        l,x1,x2,xx(i1),xx(i2),yy(i1),yy(i2),zz(i1),zz(i2),det   
c      write(*,*)"l,i1,i2,x1,x2,det",l,i1,i2,x1,x2,det
         call sjymec(x1,aj,ay,ncomp,l)
         call sjymec(x2,aj1,ay1,ncomp,l)
         f1=x1*(-ay(l)+ai*aj(l))
         f2=x2*(-ay1(l)+ai*aj1(l))
c      write(*,*)"l,x1,x2,ay(l),ai,aj(l),ay1(l),aj1(l)",f1,f2,
c     &     l,x1,x2,ay(l),ai,aj(l),ay1(l),aj1(l),f1,f2
         c1=f1*zz(i2)-f2*zz(i1)
         c1=c1/det
         c2=f2*yy(i1)-f1*yy(i2)
         c2=c2/det
         j=0
         do 3 i=1,index,nint
            j=j+1
            ww(j)=c1*yy(i)+c2*zz(i)
c            if(l.eq.1.and.j.eq.2)then
c               write(*,*)ww(j),c1,c2,xx(i),yy(i),zz(i)
c            endif
            ww(j)=ww(j)/xx(i)
c            if(l.eq.1.and.j.eq.2)then
c               write(*,*)ww(j),c1,c2,xx(i),yy(i),zz(i)
c            endif
            wwl(j,l)=ww(j)
c            if(l.eq.1.and.j.eq.2)then
c               write(*,*)ww(j),c1,c2,xx(i),yy(i),zz(i)
c            endif
 3       continue
         yc1=(ww(2)-ww(1))/rd
         yc2=(ww(np)-ww(np-1))/rd
         call splinec(r,ww,np,yc1,yc2,scc,csc(1,l))
         c2r=c2
         do 5 i=1,np
            der(i)=w(i,l)*c2r
            der(i)=der(i)/r(i)**ipow(l)
            derl(i,l)=der(i)
 5       continue
         yp1=(der(2)-der(1))/rd
         yp2=(der(np)-der(np-1))/rd
         call cckohn_spline(r,der,np,yp1,yp2,scr,csl(1,l))
 33   continue
      call dcopy(np,r,1,xr,1)
c      write(*,*)"r"
c      write(*,"(6f12.5)")(r(i),i=1,np)
c      write(*,*)"xr"
c      write(*,"(6f12.5)")(xr(i),i=1,np)
c      write(*,*)"wwl"
c      write(*,"(6f12.5)")((wwl(i,k),i=1,np),k=0,lmax)
c      write(*,*)"derl"
c      write(*,"(6f12.5)")((derl(i,k),i=1,np),k=0,lmax)
c      write(*,*)"csc"
c      write(*,"(6f12.5)")((csc(i,k),i=1,np),k=0,lmax)
c      write(*,*)"csl"
c      write(*,"(6f12.5)")((csl(i,k),i=1,np),k=0,lmax)
c      write(8)lmax,np,xstart,rd,alpha
c      write(8)(r(i),i=1,np)
c      np2=2*np
c      write(8)((wwlr(i,k),i=1,np2),k=0,lmax)
c      write(8)((derl(i,k),i=1,np),k=0,lmax)
c      write(8)((cscr(i,k),i=1,np2),k=0,lmax)
c      write(8)((csl(i,k),i=1,np),k=0,lmax)
      write(6,*)' BessplG Finished '
      return
      end
c $Id$
