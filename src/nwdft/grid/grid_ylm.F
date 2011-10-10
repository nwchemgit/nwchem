      subroutine grid_ylm(qxyz,narg,lmax,mmax,ylm)
      implicit none
c
c     computes ylms on grid
c     
c     real valued ylms normed properly to 1 - 10/6/88
c     
c     
      double precision qxyz(*)
      integer narg,lmax,mmax
      double precision ylm(1:narg,0:lmax,0:2*mmax)
c
      double precision fourpi,tol
      integer llmax
      parameter (llmax=30)
      parameter (fourpi=4.*3.14159265358979,tol=1.e-8)
c
      complex*16 cmp(narg)
c      double precision scr
c      common/stuff/scr(7*margs)
c      double precision buff(narg*2*llmax),biff(4*narg)
      double precision p(narg,0:llmax),x(narg)
      double precision xx(narg),yy(narg),zz(narg),r(narg),
     &     sphi(narg),cmphi(narg),smphi(narg),cphi(narg)
      double precision pp(narg,0:llmax),pm(narg,0:llmax),fac(0:60)
      double precision xsq, xqq, const
c
      integer lmmax,i,j,mu,mumax
c
c      open(6,file='outass',status='unknown')
c      open(5,file='inass',status='unknown')
c      open(10,file='ylms',form='unformatted',status='unknown')
c      read(5,*,end=999)npt,mumax,lmax
c      irecl=4*npt*8
c      call openabs(9,'grid',irecl)
c      call rdabs(9,it,3,0)
c      write(6,*)it
c      nbig=it(3)
c      if(it(1).ne.npt)then
c         write(6,*)" stopping because grid buffer length is wrong"
c         write(6,*)it(1),it(2),it(3)
c         write(6,*)npt,irecl
c         stop
c      endif
c
      lmmax=lmax+mmax
      fac(0)=1.
      do  i=1,lmmax
         fac(i)=i*fac(i-1)
      enddo
c
c read in a block of grid points and transfer to a temporary location
c
c      iset=1
c      iread=0
c      iquit=0
c
c      write(10)nbig,npt,lmax,mumax
c      marg=min0(nbig,npt)
c      nread=4*marg
c      call rdabs(9,biff(1),nread,iset)
c      iset=iset+1
c      iread=iread+marg
c 32    continue
c      narg=marg
c
      call dcopy(narg,qxyz(1),3,xx(1),1)
      call dcopy(narg,qxyz(2),3,yy(1),1)
      call dcopy(narg,qxyz(3),3,zz(1),1)
c
c      iremn=nbig-iread
c      if(iremn.eq.0)then
c        iquit=1
c        go to 34
c      endif
c      marg=min0(iremn,npt)
c      nread=4*marg
c      call rdabs(9,biff(1),nread,iset)
c      iset=iset+1
c      iread=iread+marg
c34    continue
c
      do 3 i=1,narg
         r(i)=sqrt(xx(i)**2+yy(i)**2+zz(i)**2)
         x(i)=zz(i)/r(i)
         xsq=x(i)**2
         xqq=abs(1.-xsq)
         if(xqq.lt.tol)then
            cphi(i)=1.   
            sphi(i)=0.
            if(x(i))66,67,68
 66         x(i)=-1.
            go to 3
 67         write(6,*)" screwup in ass1"
            stop
 68         x(i)=1.
         else
            cphi(i)=xx(i)/r(i)/sqrt(1.-xsq)
            sphi(i)=yy(i)/r(i)/sqrt(1.-xsq)
         endif
 3    continue
      mumax = mmax
      do 1 mu=0,mumax
         call plm(x,narg,mu,lmax,p)
         if(mu.eq.0)then
            do 101 i=0,lmax
               const=sqrt((2*i+1)/fourpi)
               do 10 j=1,narg
                  p(j,i)=p(j,i)*const
 10            continue
 101        continue
            do 11 i=0,lmax
               call dcopy(narg,p(1,i),1,ylm(1,i,0),1)
 11         continue
c            jbuf=narg*(lmax+1)
c     write(10)(buff(i),i=1,jbuf)
         else
            do 4 i=1,narg
               cmp(i)=(cmplx(cphi(i),sphi(i)))**mu
               cmphi(i)=dble(cmp(i))
               smphi(i)=dimag(cmp(i))
 4          continue
            do 121 i=mu,lmax
               const=sqrt((2*i+1)/fourpi*fac(i-mu)/fac(i+mu))
c     sqrt(2) factor added to normalize "real valued" ylms
               const=const*sqrt(2.0)
               do 12 j=1,narg
                  pp(j,i)=p(j,i)*const*cmphi(j)
                  pm(j,i)=p(j,i)*const*smphi(j)
 12            continue
 121        continue
            do 13 i=mu,lmax
               call dcopy(narg,pp(1,i),1,ylm(1,i,2*mu-1),1)
               call dcopy(narg,pm(1,i),1,ylm(1,i,2*mu),  1)
 13         continue
c            jbuf=2*narg*(1+lmax-mu)
c     write(10)(buff(i),i=1,jbuf)
         endif
 1    continue
      do j = 0,2*mmax
         do i = 0,lmax
c            write(*,*)"ylm(narg,llmax,mmax),l=,m=", i, j
c            write(*,"(2i5,4f10.4)")
c     &           (i,j,xx(mu),yy(mu),zz(mu),ylm(mu,i,j), mu=1,narg)
         enddo
      enddo
c
c return to start a fetch another block of points
c
c      if(iquit.eq.0)go to 32
c      do 2 i=1,narg
c2     write(6,100)x(i),mu,(p(i,j),j=mu,lmax)
c100   format(" x=",f10.5,"  mu = ",i3/(10e12.4))
c999   continue
c      call closeabs(9)
      return
      end
c
c     subroutine for calculating associated Legendre function
c
      subroutine plm(x,ngrid,mu,lmax,p)
      implicit none
      integer ngrid,n,mu,lmax,i,l,mu1,l0,l1,m,mu2
      double precision p(ngrid,0:2),x(ngrid),p0(ngrid),p1(ngrid),
     &     p2(ngrid), arg(ngrid),q0(ngrid),q1(ngrid),q2(ngrid)
c     
c     this routine calculates a vector of associated legengre 
c     polynomials P(l,mu) for fixed mu and l=mu,,,lmax for an array 
c     of n points, x(i),i=1,,n.
c
c     mu = 0 case
c     
      n=ngrid
      if(mu.eq.0)then
         do 2 i=1,n
            p(i,0)=1.
            p(i,1)=x(i)
 2       continue
         if(lmax.lt.2)return
         do  l=2,lmax
            do  i=1,n
               p(i,l)=((2*l-1)*x(i)*p(i,l-1)-(l-1)*p(i,l-2))/l
            enddo
         enddo
         return
      endif
c     
c     mu = 1 case
c     
      if(mu.eq.1)then
         do  i=1,n
            arg(i)=sqrt(1.-x(i)*x(i))
            p(i,1)=-arg(i)
         enddo
         if(lmax.lt.2)return
         do  i=1,n
            p(i,2)=-3.*x(i)*arg(i)
         enddo
         if(lmax.lt.3)return
         do  l=3,lmax
            do  i=1,n
               p(i,l)=((2*l-1)*x(i)*p(i,l-1)-l*p(i,l-2))/(l-1)
            enddo
         enddo
         return
      endif
c
c mu must be larger than 1
c
c
c     recurr across to l=mu+1, with m=0
c     
      mu1=mu+1
      do i=1,n
         p0(i)=1.
         p1(i)=x(i)
      enddo
      do l=2,mu1
         do i=1,n
            p2(i)=((2*l-1)*x(i)*p1(i)-(l-1)*p0(i))/l
            p0(i)=p1(i)
            p1(i)=p2(i)
         enddo
      enddo
c     
c     recurr across to l=mu+1, with m=1
c     
      do i=1,n
         arg(i)=sqrt(1.-x(i)*x(i))
         q0(i)=-arg(i)
         q1(i)=-3.*arg(i)*x(i)
      enddo
      do l=3,mu1
         do i=1,n
            q2(i)=((2*l-1)*x(i)*q1(i)-l*q0(i))/(l-1)
            q0(i)=q1(i)
            q1(i)=q2(i)
         enddo
      enddo
c     
c     with l fixed at mu and mu+1, recurr down to m=mu
c     
      l0=mu
      l1=mu1
      do m=2,mu
         do i=1,n
            if(arg(i).eq.0.)then   
               q0(i)=0.
               q1(i)=0.
            else
               p2(i)=-2.*(m-1)*x(i)*q0(i)/arg(i)-
     &              (l0-m+2)*(l0+m-1)*p0(i)
               q2(i)=-2.*(m-1)*x(i)*q1(i)/arg(i)-
     &              (l1-m+2)*(l1+m-1)*p1(i)
               p0(i)=q0(i)
               p1(i)=q1(i)
               q0(i)=p2(i)
               q1(i)=q2(i)
            endif
         enddo
      enddo
c     
c     compute the desired vector and quit
c     
      do 12 i=1,n
         p(i,mu)=q0(i)
 12   continue
      if(lmax.eq.mu)return
      do 13 i=1,n
         p(i,mu1)=q1(i)
 13   continue
      if(lmax.eq.mu1)return
      mu2=mu1+1
      do l=mu2,lmax
         do i=1,n
            p(i,l)=((2*l-1)*x(i)*p(i,l-1)-(l+mu-1)*p(i,l-2))/(l-mu)
         enddo
      enddo
      return
      end
      
c $Id$
